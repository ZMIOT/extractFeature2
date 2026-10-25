#include "SeamSmoother.h"
#include "SeamDetector.h"
#include "UniformSeamWidth.h"
#include "AlgorithmParameters.h"
#include "HullSurfaces.h"
#include "PointCloudIO.h"
#include "ProfileCurvature.h"
#include <fstream>
#include <cstdio>

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

#include "wall_bead_profiles.inc"

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

void testIndependentAlgorithmParameters() {
  const auto curve = weld::ProfileLocalizationMethod::CurveLineIntersection;
  const auto curvature = weld::ProfileLocalizationMethod::CurvatureCandidates;
  weld::AlgorithmParameters banks;
  auto first = banks.get(curve);
  first.profileCurvatureRadius = 5;
  first.finalSmoothingLength = 27;
  first.curveLineFallback.splineRadius = 3;
  first.minimumFeatureDepth = 2;
  first.uniformWidthFilterLength=21;first.uniformWidthFilterEnabled=false;
  banks.set(curve, first);
  auto second = banks.get(curvature);
  require(second.profileCurvatureRadius == 8 && second.finalSmoothingLength == 12 &&
          second.minimumFeatureDepth == 1 && second.minSeamSeparation == 6 &&
          second.uniformWidthFilterLength==12 && second.uniformWidthFilterEnabled &&
          banks.get(curve).minSeamSeparation == 8,
          "Algorithm parameter banks share values");
  second.profileCurvatureRadius = 16;
  second.curvatureSearch.continuityWeight = 9;
  banks.set(curvature, second);
  banks.reset(curvature);
  require(banks.get(curve).profileCurvatureRadius == 5 &&
          banks.get(curve).finalSmoothingLength == 27 &&
          banks.get(curve).curveLineFallback.splineRadius == 3 &&
          banks.get(curve).minimumFeatureDepth == 2,
          "Resetting curvature modified curve-line parameters");
}

void testUniformWidthFinalFilter() {
  std::vector<weld::UniformWidthSection> sections;
  auto makeResult=[] {weld::SeamResult r;r.frameBased=true;return r;};
  auto raw=makeResult(),filtered=makeResult();
  for(int i=0;i<100;++i) {
    const float station=.5f*i;
    const float shift=.003f*station*station+(i>=50?8:0)+.35f*std::sin(2.f*i);
    const float noise=.3f*std::cos(1.8f*i);
    auto point=[&](float t){const float y=shift+t;return weld::PointT(station+.15f*y,y,.02f*t*t+.1f*station+noise);};
    weld::UniformWidthSection s;s.station=station;s.bottomCenterU=shift;
    for(int k=-30;k<=30;++k)s.contour.push_back(point(float(k)));
    s.sourceLeft=point(-5);s.sourceRight=point(5);sections.push_back(s);
    for(auto *r:{&raw,&filtered}) {
      r->leftSeam->push_back(s.sourceLeft);r->rightSeam->push_back(s.sourceRight);
      r->leftFrameIndices.push_back(i);r->rightFrameIndices.push_back(i);
    }
  }
  raw.leftSegmentStarts=raw.rightSegmentStarts=filtered.leftSegmentStarts=filtered.rightSegmentStarts={0,50};
  weld::DetectorConfig c;c.uniformWidthEnabled=true;c.uniformWidthFilterEnabled=false;
  weld::applyUniformSeamWidth(raw,c,sections);
  c.uniformWidthFilterEnabled=true;weld::applyUniformSeamWidth(filtered,c,sections);
  require(filtered.uniformWidthFilterApplied && !raw.uniformWidthFilterApplied,"Post-width filter toggle was ignored");
  require(raw.leftFrameIndices==filtered.leftFrameIndices && filtered.leftSegmentStarts==raw.leftSegmentStarts,
          "Post-width filter changed frame identity or segment boundaries");
  auto roughness=[](const weld::SeamResult &r) {
    double sum=0;for(int i=1;i<99;++i)if(i!=49 && i!=50) {
      const auto a=weld::toVector((*r.leftSeam)[i-1]),b=weld::toVector((*r.leftSeam)[i]),d=weld::toVector((*r.leftSeam)[i+1]);
      sum+=(a-2*b+d).squaredNorm();
    }return sum;
  };
  require(roughness(filtered)<.05*roughness(raw),"Post-width filtering did not remove frame oscillation");
  for(int i=0;i<100;++i) {
    const auto l=weld::toVector((*filtered.leftSeam)[i]),r=weld::toVector((*filtered.rightSeam)[i]);
    require(std::abs((r-l).norm()-raw.uniformWidthTarget)<1.e-4f,"Filtering changed the fixed chord width");
    require(std::abs(l.x()-.15f*l.y()-.5f*i)<1.e-4f && std::abs(r.x()-.15f*r.y()-.5f*i)<1.e-4f,
            "Filtering moved points off their tilted scan plane");
    require(l.y()<sections[i].bottomCenterU && r.y()>sections[i].bottomCenterU,"Filtering crossed the bottom midpoint");
  }
  require((*filtered.leftSeam)[50].y-(*filtered.leftSeam)[49].y>7,"Filtering bridged a real segment gap");
}

void testBottomMidpointGuard() {
  for (int opening : {1,-1}) {
    std::vector<Eigen::Vector2f> points;
    for (int j=-600;j<=260;++j) {
      const float u=.1f*j;
      float height=std::min(18.f,std::max(0.f,2*(std::fabs(u-7.f)-9.f)));
      if (j==-400) height=-100; // isolated outlier must not move the floor
      points.emplace_back(u,opening*height);
    }
    const float middle=weld::profileBottomMidpoint(points,opening,.6f);
    require(std::fabs(middle-7.f)<.15f,"Bottom midpoint followed scan margins or an isolated spike");
    require(weld::profilePairStraddlesBottom(-2,16,middle),"Valid pair rejected by bottom split");
    require(!weld::profilePairStraddlesBottom(8,16,middle) &&
            !weld::profilePairStraddlesBottom(-2,6,middle),"Same-side pair crossed the bottom midpoint");
    require(!weld::profilePairStraddlesBottom(middle,16,middle) &&
            !weld::profilePairStraddlesBottom(-2,middle,middle),"Midpoint itself became a side candidate");
    for (auto axis : {weld::Axis::X,weld::Axis::Y,weld::Axis::Z}) {
      const int uAxis=axis==weld::Axis::X?1:0;
      const int vAxis=axis==weld::Axis::Z?1:2;
      weld::CloudT::Ptr cloud(new weld::CloudT);
      std::vector<weld::ScanFrameRange> frames;
      for (int frame=0;frame<10;++frame) {
        const auto begin=cloud->size();
        for (const auto &uv:points) {
          weld::PointT p;
          p.data[int(axis)]=.5f*frame; p.data[uAxis]=uv.x(); p.data[vAxis]=uv.y();
          // Use a clean measured fixture; isolated noise was covered above.
          if (std::fabs(uv.x()+40.f)<.01f) p.data[vAxis]=opening*18.f;
          cloud->push_back(p);
        }
        frames.push_back({begin,cloud->size()-begin,frame});
      }
      for (auto method : {weld::ProfileLocalizationMethod::CurveLineIntersection,
                          weld::ProfileLocalizationMethod::CurvatureCandidates}) {
        weld::DetectorConfig c;
        c.extractMainBody=false; c.voxelSize=0; c.outlierMeanK=0;
        c.trimUnstableEndpoints=false; c.longitudinalAxis=axis;
        c.interiorOpeningDirection=opening; c.profileLocalizationMethod=method;
        weld::DebugTrace trace;
        const auto result=weld::SeamDetector().detect(cloud,c,nullptr,nullptr,&trace,&frames);
        require(result.ok(),"Midpoint guard lost a shifted groove with asymmetric scan margins");
        std::size_t checked=0;
        for (const auto &slice:trace.slices) for (const auto &pair:slice.pairs) {
          require(pair.hasBottomMidpoint && std::fabs(pair.bottomMidpointU-7)<1,
                  "Candidate pair lost its independently measured bottom reference");
          require(pair.points[0].data[uAxis]<pair.bottomMidpointU &&
                  pair.points[1].data[uAxis]>pair.bottomMidpointU,
                  "Accepted candidates crossed the bottom midpoint");
          ++checked;
        }
        require(checked>0,"No candidates exercised the midpoint guard");
      }
    }
  }
  std::vector<Eigen::Vector2f> ramp;
  for (int i=0;i<=200;++i) ramp.emplace_back(.1f*i,.2f*i);
  const float middle=weld::profileBottomMidpoint(ramp,1,.6f);
  require(middle<.5f && !weld::profilePairStraddlesBottom(8,9,middle),
          "Two points on one rising wall passed as opposite floor junctions");
}

void testShallowFeatureCutoff() {
  for (int opening : {1, -1}) {
    weld::CloudT::Ptr cloud(new weld::CloudT);
    std::vector<weld::ScanFrameRange> frames;
    for (int f=0; f<54; ++f) {
      const auto begin=cloud->size();
      const bool deep=(f>=6 && f<=35) || (f>=37 && f<=46);
      for (int j=-200; j<=200; ++j) {
        const float u=.1f*j;
        float v=std::min(deep ? 6.f : .25f, std::max(0.f, 2*(std::fabs(u)-8)));
        v += .05f*u + .02f*std::sin(3*u);
        if (!deep && j==100) v+=8; // isolated noise cannot restore a lost groove
        cloud->push_back({f*.5f,u,opening*v});
      }
      frames.push_back({begin,cloud->size()-begin,f});
    }
    for (auto method : {weld::ProfileLocalizationMethod::CurveLineIntersection,
                        weld::ProfileLocalizationMethod::CurvatureCandidates}) {
      weld::DetectorConfig c;
      c.extractMainBody=false; c.voxelSize=0; c.outlierMeanK=0;
      c.profileLocalizationMethod=method; c.interiorOpeningDirection=opening;
      c.trimUnstableEndpoints=false; c.mergeProfileFrames=true; c.profileMergeFrameCount=5;
      c.maxInterpolatedGapStations=50; // generic gap recovery must not bypass a depth stop
      weld::DebugTrace trace;
      auto result=weld::SeamDetector().detect(cloud,c,nullptr,nullptr,&trace,&frames);
      require(result.ok(),"Depth cutoff lost the established deep groove");
      require(result.leftFrameIndices.front()>=6 && result.leftFrameIndices.back()<=35 &&
              result.leftFrameIndices.back()>=30,
              "Output crossed a shallow section or discarded the established groove");
      require(result.leftFrameIndices==result.rightFrameIndices,"Depth cutoff split left/right output ranges");
      bool checkedBoundary=false;
      for (const auto &slice:trace.slices) {
        if (slice.frameIndex==36 || slice.frameIndex<6 || slice.frameIndex>=47) {
          require(slice.pairs.empty() && slice.featureDepth>=0 && slice.featureDepth<c.minimumFeatureDepth &&
                  slice.gates.count("feature_depth_cutoff"),
                  "Merged neighbors, noise or fallback recovered a planar section");
          require(!slice.noCandidateReason.empty(),"Depth cutoff has no debug explanation");
          checkedBoundary=true;
        }
      }
      require(checkedBoundary,"Depth cutoff test did not inspect its boundary");
      c.minimumFeatureDepth=100;
      require(!weld::SeamDetector().detect(cloud,c,nullptr,nullptr,nullptr,&frames).ok(),
              "All-shallow sections reached path optimization");
      c.minimumFeatureDepth=-1;
      require(!weld::SeamDetector().detect(cloud,c,nullptr,nullptr,nullptr,&frames).ok(),
              "Invalid minimum feature depth was accepted");
    }
  }
}

void testSupplementalCandidates() {
  std::vector<weld::ProfileCurvatureSample> samples;
  for (int i=-100;i<=100;++i) {
    const float x=.05f*i;
    weld::ProfileCurvatureSample s;
    s.point={x,.25f*x*x}; s.valid=true;
    s.curvature=.5f/std::pow(1.f+.25f*x*x,1.5f);
    samples.push_back(s);
  }
  const auto kinds=weld::supplementalProfileCandidates(samples,2,.08f,1);
  require((kinds[100]&2)!=0,"Physical-window valley candidate missing");
  require(std::any_of(kinds.begin(),kinds.end(),[](unsigned k){return (k&4)!=0;}),
          "Curvature-change candidates missing");
  for (auto &s:samples) { s.point.y()=-s.point.y(); s.curvature=-s.curvature; }
  const auto reverse=weld::supplementalProfileCandidates(samples,2,.08f,-1);
  require(reverse==kinds,"Candidate search does not mirror with groove opening");
  const auto wrong=weld::supplementalProfileCandidates(samples,2,.08f,1);
  require((wrong[100]&2)==0,"Valley search accepted opposite opening");
  samples[100].valid=false;
  const auto gap=weld::supplementalProfileCandidates(samples,2,.08f,-1);
  require(gap[100]==0 && gap[99]==0 && gap[101]==0,"Candidate window crossed invalid samples");
  for (auto &s:samples) { s.valid=true; s.point.y()=0; s.curvature=0; }
  const auto flat=weld::supplementalProfileCandidates(samples,2,.08f,1);
  require(std::all_of(flat.begin(),flat.end(),[](unsigned k){return k==0;}),
          "Flat profile produced supplemental features");
}

void testWallTurnCentroid() {
  std::vector<weld::ProfileCurvatureSample> samples(641);
  auto locate=[&](float balance,int opening) {
    for (int i=0;i<int(samples.size());++i) {
      const float x=(i-320)*.025f;
      auto &p=samples[i]; p.valid=true; p.point={x,0};
      p.curvature=opening*((1-balance)*std::exp(-4*(x+1)*(x+1))+
                           (1+balance)*std::exp(-4*(x-1)*(x-1)));
    }
    return weld::wallSupportedTurnCentroid(samples,{0,0},{0,1},0,.6f,4,opening,-5,1);
  };
  const int before=locate(-.02f,1), after=locate(.02f,1);
  require(before>=0 && after>=0 && std::abs(after-before)*.025f<.15f,
          "Nearly equal neighboring lobes made the turning centroid jump");
  require(locate(.02f,-1)==after,"Reversed opening changed turning centroid identity");
  locate(.02f,1);
  samples[600].curvature=1000; // Outside the physical wall neighborhood.
  require(weld::wallSupportedTurnCentroid(samples,{0,0},{0,1},0,.6f,4,1,-5,1)==after,
          "Remote sharp feature captured the wall turning centroid");
  for (auto &p:samples) p.curvature=0;
  require(weld::wallSupportedTurnCentroid(samples,{0,0},{0,1},0,.6f,4,1,-5,1)<0,
          "Flat plane generated a turning centroid");
}

void testMeasuredCurvatureStep() {
  std::string error;
  std::vector<weld::ScanFrameRange> frames;
  auto cloud=weld::PointCloudIO::load(
      std::string(SEAM_TEST_DATA_DIR)+"/curvature_2_frames_2040_2110.txt",&error,&frames);
  require(cloud && !cloud->empty(),"Two-stage wall regression fixture is missing");
  auto c=weld::AlgorithmParameters().get(weld::ProfileLocalizationMethod::CurvatureCandidates);
  c.extractMainBody=false; c.voxelSize=0; c.outlierMeanK=0;
  c.trimUnstableEndpoints=false; c.finalSmoothingEnabled=false;
  weld::DebugTrace trace;
  require(weld::SeamDetector().detect(cloud,c,nullptr,nullptr,&trace,&frames).ok(),
          "Two-stage wall lost the measured curvature path");
  std::vector<float> widths;
  const weld::DebugPair *previous=nullptr;
  for (const auto &slice:trace.slices) if (slice.frameIndex>=2045 && slice.frameIndex<=2105) {
    require(slice.pairs.size()==1 && slice.selected==0,"Stable mode retained competing wall turns");
    const auto &pair=slice.pairs[0];
    require(pair.hasBottomMidpoint && pair.points[0].y<pair.bottomMidpointU &&
            pair.points[1].y>pair.bottomMidpointU,"Turning centroid crossed bottom midpoint");
    for (int side=0;side<2;++side) {
      const auto &p=pair.points[side];
      const auto &wall=pair.candidateWalls[side];
      require(std::fabs(-wall.dv*(p.y-wall.u)+wall.du*(p.z-wall.v))<=c.profileWallTolerance+1.e-4f,
              "Turning centroid lost independent wall support");
      if (previous) {
        const auto &q=previous->points[side];
        require(std::hypot(p.y-q.y,p.z-q.z)<1.f,"Two-stage wall jumped between curvature lobes");
      }
    }
    widths.push_back(pair.points[1].y-pair.points[0].y); previous=&pair;
  }
  require(widths.size()==61,"Step regression did not cover the complete neighborhood");
  for (std::size_t i=12;i<widths.size();++i)
    require(std::fabs(widths[i]-widths[i-12])<1.5f,"One-sided width step survived stable localization");
}

void testUniformBottomWidth() {
  for(auto metric:{weld::UniformWidthMetric::BottomArc,weld::UniformWidthMetric::Chord})
  for(auto axis:{weld::Axis::X,weld::Axis::Y,weld::Axis::Z}) {
    weld::DetectorConfig c;c.uniformWidthEnabled=true;c.uniformWidthMetric=metric;c.longitudinalAxis=axis;
    weld::SeamResult result;result.frameBased=true;
    std::vector<weld::UniformWidthSection> sections;
    const int a=int(axis),u=a==0?1:0,v=a==2?1:2;
    const float scale=std::sqrt(metric==weld::UniformWidthMetric::BottomArc?1.29f:1.04f);
    for(int frame=0;frame<20;++frame) {
      auto point=[&](float y){Eigen::Vector3f p;p[a]=frame+.2f*y;p[u]=y;p[v]=.5f*std::abs(y);return weld::toPoint(p);};
      const float width=frame==0?1.f:frame==1?2.f:frame==18?30.f:frame==19?40.f:10.f;
      weld::UniformWidthSection s;
      s.contour={point(-30),point(0),point(30)};s.bottomCenterU=0;
      s.sourceLeft=point(-width/(2*scale));s.sourceRight=point(width/(2*scale));
      sections.push_back(s);result.leftSeam->push_back(s.sourceLeft);result.rightSeam->push_back(s.sourceRight);
      result.leftFrameIndices.push_back(frame+100);result.rightFrameIndices.push_back(frame+100);
      result.leftPointQuality.emplace_back();result.rightPointQuality.emplace_back();
      result.leftPointConfidence.push_back(1);result.rightPointConfidence.push_back(1);
    }
    result.leftSegmentStarts=result.rightSegmentStarts={0,15};
    weld::applyUniformSeamWidth(result,c,sections);
    require(result.uniformWidthSamples==20 && result.uniformWidthTrimEach==2 &&
            std::abs(result.uniformWidthTarget-10)<1.e-4f,"Trimmed mean did not remove each ten-percent tail");
    require(result.uniformWidthAdjusted==20 && result.uniformWidthOmitted==0,"Supported contour lost points");
    for(int i=0;i<20;++i) {
      const auto l=weld::toVector((*result.leftSeam)[i]),r=weld::toVector((*result.rightSeam)[i]);
      require(std::abs(l[u]+5/scale)<1.e-4f && std::abs(r[u]-5/scale)<1.e-4f,
              "Uniform pair did not extend symmetrically along the measured bottom");
      require(std::abs(l[a]-i-.2f*l[u])<1.e-4f && std::abs(r[v]-.5f*r[u])<1.e-4f,
              "Uniform pair left the original tilted scan plane or contour");
      if(metric==weld::UniformWidthMetric::Chord) require(std::abs((r-l).norm()-10)<1.e-4f,"Chord width is not constant");
    }
    sections[9].contour.clear();
    weld::applyUniformSeamWidth(result,c,sections);
    require(result.uniformWidthOmitted==1 && result.leftSeam->size()==19 &&
            result.leftFrameIndices[9]==110 && result.rightPointQuality.size()==19,
            "Missing contour extrapolated or broke metadata alignment");
    require(result.leftSegmentStarts==std::vector<std::size_t>({0,9,14}),"Uniform-width omission bridged a missing section");
    require(result.leftPointQuality[0].flags&weld::QualityUniformWidth,"Uniform adjustment was not recorded");
  }
}

void testMeasuredCurveLineTail() {
  std::string error;
  std::vector<weld::ScanFrameRange> frames;
  auto cloud=weld::PointCloudIO::load(
      std::string(SEAM_TEST_DATA_DIR)+"/curve_line_2_frames_2450_2540.txt",&error,&frames);
  require(cloud && !cloud->empty(),"Measured curve-line tail fixture is missing");
  weld::DetectorConfig c;
  c.extractMainBody=false; c.voxelSize=0; c.outlierMeanK=0;
  c.finalSmoothingEnabled=false;
  weld::DebugTrace trace;
  const auto result=weld::SeamDetector().detect(cloud,c,nullptr,nullptr,&trace,&frames);
  require(result.ok(),"Measured wavy floor lost its curve-line path");
  const weld::DebugPair *previous=nullptr;
  int checked=0; bool selectedOuterAlternative=false;
  for (const auto &slice:trace.slices) if (slice.frameIndex>=2459 && slice.frameIndex<=2529) {
    require(slice.selected>=0,"Reported curve-line frame lost its selected pair");
    const auto &pair=slice.pairs[slice.selected];
    require(pair.source!="central_bead_valleys","Marginal bead stole the outer junction path");
    selectedOuterAlternative|=pair.source=="curve_line_outer_alternative";
    require(pair.hasBottomMidpoint && pair.points[0].y<pair.bottomMidpointU &&
            pair.points[1].y>pair.bottomMidpointU,"Curve-line candidate crossed the bottom midpoint");
    if (previous) for (int side=0;side<2;++side) {
      const auto &p=pair.points[side], &q=previous->points[side];
      require(std::hypot(p.y-q.y,p.z-q.z)<1.2f,"Wavy floor line switched levels between frames");
    }
    previous=&pair; ++checked;
  }
  require(checked==71 && selectedOuterAlternative,"Reported tail/ambiguous bead was not exercised");
  require(result.leftFrameIndices.front()<=2459 && result.leftFrameIndices.back()>=2529,
          "Reported curve-line interval was truncated");
}

void testMeasuredCurveLineBranchContinuity() {
  std::string error;
  std::vector<weld::ScanFrameRange> frames;
  auto cloud=weld::PointCloudIO::load(
      std::string(SEAM_TEST_DATA_DIR)+"/curve_line_3_frames_1430_1560.txt",&error,&frames);
  require(cloud && !cloud->empty(),"Measured inner/outer branch fixture is missing");
  weld::DetectorConfig config;
  config.extractMainBody=false; config.voxelSize=0; config.outlierMeanK=0;
  config.finalSmoothingEnabled=false; config.trimUnstableEndpoints=false;
  weld::DebugTrace trace;
  const auto result=weld::SeamDetector().detect(cloud,config,nullptr,nullptr,&trace,&frames);
  require(result.ok(),"Measured branch fixture lost its path");
  int ambiguous=0, checked=0;
  const weld::DebugPair *previous=nullptr;
  for (const auto &slice:trace.slices) if (slice.frameIndex>=1440 && slice.frameIndex<=1555) {
    require(slice.selected>=0,"Measured branch interval lost a frame");
    const auto &pair=slice.pairs[slice.selected];
    require(pair.source!="central_bead_valleys","Inner shape bonus displaced a continuous outer track");
    for (const auto &alternative:slice.pairs)
      ambiguous+=alternative.source=="central_bead_valleys";
    require(pair.hasBottomMidpoint && pair.points[0].y<pair.bottomMidpointU &&
            pair.points[1].y>pair.bottomMidpointU,"Branch crossed bottom midpoint");
    if (previous) for (int side=0;side<2;++side) {
      const auto &p=pair.points[side], &q=previous->points[side];
      require(std::hypot(p.y-q.y,p.z-q.z)<1.2f,"Outer path jumped at bead onset/exit");
    }
    previous=&pair; ++checked;
  }
  require(checked==116 && ambiguous>=80,"Regression did not exercise sustained branch ambiguity");
  for(auto method:{weld::ProfileLocalizationMethod::CurveLineIntersection,weld::ProfileLocalizationMethod::CurvatureCandidates}) {
    auto uniform=weld::AlgorithmParameters().get(method);
    uniform.extractMainBody=false;uniform.voxelSize=0;uniform.outlierMeanK=0;
    uniform.trimUnstableEndpoints=false;uniform.uniformWidthEnabled=true;
    uniform.uniformWidthMetric=weld::UniformWidthMetric::Chord;
    // Production no-debug path, including final smoothing and restoration.
    const auto fixed=weld::SeamDetector().detect(cloud,uniform,nullptr,nullptr,nullptr,&frames);
    require(fixed.ok() && fixed.uniformWidthAdjusted>100,"Measured uniform-width mode lost the groove");
    require(fixed.uniformWidthSamples>=fixed.uniformWidthAdjusted && fixed.uniformWidthTrimEach==fixed.uniformWidthSamples/10,
            "Measured trimmed mean used the wrong sample count");
    for(std::size_t i=0;i<fixed.leftSeam->size();++i)
      require(std::abs((weld::toVector((*fixed.rightSeam)[i])-weld::toVector((*fixed.leftSeam)[i])).norm()-fixed.uniformWidthTarget)<.002f,
              "Original scan-plane restoration changed the uniform chord width");
  }
}

void testMeasuredCurvatureStability() {
  // Measured body points from 2.txt: debug ordinals 790-792 are source frame
  // IDs 2049-2051. Include neighboring frames for the normal five-frame fit.
  std::string error;
  std::vector<weld::ScanFrameRange> frames;
  auto cloud=weld::PointCloudIO::load(
      std::string(SEAM_TEST_DATA_DIR)+"/curvature_2_frames_2042_2058.txt", &error, &frames);
  require(cloud && !cloud->empty(),"Measured curvature regression fixture is missing");
  auto c=weld::AlgorithmParameters().get(weld::ProfileLocalizationMethod::CurvatureCandidates);
  c.extractMainBody=false; c.voxelSize=0; c.outlierMeanK=0;
  c.trimUnstableEndpoints=false; c.finalSmoothingEnabled=false;
  weld::DebugTrace trace;
  const auto result=weld::SeamDetector().detect(cloud,c,nullptr,nullptr,&trace,&frames);
  require(result.ok(),"Measured narrow groove lost its curvature path");
  const weld::DebugPair *previous=nullptr;
  int checked=0;
  for (const auto &slice:trace.slices) if (slice.frameIndex>=2049 && slice.frameIndex<=2051) {
    require(!slice.pairs.empty() && slice.pairs.size()<=2,
            "Similar real sections produced a broad redundant candidate set");
    require(slice.selected>=0,"Measured curvature frame has no selected candidate");
    const auto &pair=slice.pairs[slice.selected];
    if (previous) for (int side=0;side<2;++side) {
      const auto &p=pair.points[side], &q=previous->points[side];
      require(std::hypot(p.y-q.y,p.z-q.z)<1.f,
              "Selected real curvature candidate jumped to another wall turn");
    }
    previous=&pair; ++checked;
  }
  require(checked==3,"Measured regression did not cover all three reported frames");
}

void testShallowLowerTurns() {
  // Test a class of short-wall grooves, including tilt, inverted opening,
  // unequal wall slopes and phase-varying measurement noise. Assert geometry
  // before endpoint trimming/final smoothing can hide a wrong shoulder pick.
  for (float depth:{.8f,1.5f,3.f,10.f})
  for (int opening:{-1,1})
  for (float tilt:{-.18f,.23f}) {
    weld::CloudT::Ptr cloud(new weld::CloudT);
    std::vector<weld::ScanFrameRange> frames;
    const float co=std::cos(tilt),si=std::sin(tilt);
    for (int frame=0;frame<12;++frame) {
      const auto begin=cloud->size();
      for (int j=-250;j<=250;++j) {
        const float u=.1f*j;
        const float height=std::min(depth,std::max(0.f,(std::fabs(u)-8.f)*(u<0 ? 1.2f : .85f)));
        const float v=opening*(height+.02f*std::sin(11*u+.27f*frame));
        cloud->push_back({.5f*frame,co*u-si*v,si*u+co*v});
      }
      frames.push_back({begin,cloud->size()-begin,frame});
    }
    for (auto method:{weld::ProfileLocalizationMethod::CurvatureCandidates,
                     weld::ProfileLocalizationMethod::CurveLineIntersection}) {
      auto config=weld::AlgorithmParameters().get(method);
      config.extractMainBody=false;config.voxelSize=0;config.outlierMeanK=0;
      config.minimumFeatureDepth=0;config.trimUnstableEndpoints=false;
      config.finalSmoothingEnabled=false;config.interiorOpeningDirection=opening;
      weld::DebugTrace trace;
      const auto result=weld::SeamDetector().detect(cloud,config,nullptr,nullptr,&trace,&frames);
      if (!result.ok()) std::cerr<<"shallow depth="<<depth<<" opening="<<opening<<" tilt="<<tilt<<" method="<<int(method)<<'\n';
      require(result.ok(),"Measured shallow lower turns were not recovered");
      int selected=0;
      std::array<weld::PointT,2> previous;
      for (const auto &slice:trace.slices) if (slice.selected>=0) {
        const auto &pair=slice.pairs[slice.selected];
        for (int side=0;side<2;++side) {
          const auto &p=pair.points[side];
          const float u=co*p.y+si*p.z;
          const float height=opening*(-si*p.y+co*p.z);
          if (!(height<.55f*depth+.06f && std::fabs(u-(side ? 8.f : -8.f))<2.f))
            std::cerr<<"bad foot depth="<<depth<<" opening="<<opening<<" tilt="<<tilt<<" u="<<u<<" h="<<height<<" source="<<pair.source<<'\n';
          require(height<.55f*depth+.06f && std::fabs(u-(side ? 8.f : -8.f))<2.f,
                  "Shallow groove selected its upper shoulder instead of the lower foot");
          if (pair.source.find("emerging_profile")!=std::string::npos)
            require(opening*pair.profileCurvatures[side]>0,
                    "Fallback emitted a feature endpoint without measured lower curvature");
          if (selected && depth<=3.f) {
            const float jump=std::hypot(p.y-previous[side].y,p.z-previous[side].z);
            if (jump>=.4f) std::cerr<<"jump depth="<<depth<<" opening="<<opening<<" tilt="<<tilt
                <<" method="<<int(method)<<" frame="<<slice.frameIndex<<" side="<<side
                <<" jump="<<jump<<" source="<<pair.source<<'\n';
            require(jump<.4f,"Noisy shallow groove switched between upper and lower branches");
          }
          previous[side]=p;
        }
        ++selected;
      }
      require(selected>=10,"Shallow lower-turn path lost most of its measured frames");
    }
  }
}

void testUnsupportedShallowTurns() {
  for (int shape=0;shape<5;++shape) {
    weld::CloudT::Ptr cloud(new weld::CloudT);
    std::vector<weld::ScanFrameRange> frames;
    for (int frame=0;frame<9;++frame) {
      const auto begin=cloud->size();
      for (int j=-250;j<=250;++j) {
        const float u=.1f*j;
        if (shape==2 && std::fabs(u)<3.5f) continue;
        float v=.015f*std::sin(9*u+.2f*frame);
        if (shape==1) v+=std::min(1.5f,std::max(0.f,u)); // only one wall
        if (shape==2) v+=std::min(1.5f,std::max(0.f,(std::fabs(u)-8.f)*1.2f));
        if (shape==3 && j==0) v+=1.f; // an isolated point is not a pair of feet
        if (shape==4) v+=std::min(1.f,std::max(0.f,(std::fabs(u)-.8f)*1.2f));
        cloud->push_back({.5f*frame,u,v});
      }
      frames.push_back({begin,cloud->size()-begin,frame});
    }
    auto config=weld::AlgorithmParameters().get(weld::ProfileLocalizationMethod::CurvatureCandidates);
    config.extractMainBody=false;config.voxelSize=0;config.outlierMeanK=0;
    config.minimumFeatureDepth=0;config.trimUnstableEndpoints=false;config.finalSmoothingEnabled=false;
    weld::DebugTrace trace;
    const auto result=weld::SeamDetector().detect(cloud,config,nullptr,nullptr,&trace,&frames);
    if (result.ok()) std::cerr<<"unsupported shallow shape="<<shape<<'\n';
    require(!result.ok(),"Unsupported or too-narrow contour acquired a shallow seam");
  }
}

void testMeasuredShallowLowerTurns() {
  for (int specimen:{2,3}) {
  std::string error;
  std::vector<weld::ScanFrameRange> frames;
  const auto cloud=weld::PointCloudIO::load(
      std::string(SEAM_TEST_DATA_DIR)+(specimen==3 ? "/shallow_3_frames_375_445.txt" :
                                                    "/shallow_2_frames_1335_1445.txt"),&error,&frames);
  require(cloud && !cloud->empty(),"Measured shallow-groove fixture is missing");
  auto config=weld::AlgorithmParameters().get(weld::ProfileLocalizationMethod::CurvatureCandidates);
  config.extractMainBody=false;config.voxelSize=0;config.outlierMeanK=0;
  config.trimUnstableEndpoints=false;config.finalSmoothingEnabled=false;
  for (float minimumDepth:{0.f,1.f}) {
    config.minimumFeatureDepth=minimumDepth;
    weld::DebugTrace trace;
    const auto result=weld::SeamDetector().detect(cloud,config,nullptr,nullptr,&trace,&frames);
    require(result.ok(),"Measured short-wall groove was discarded");
    int checked=0;
    const weld::DebugPair *previous=nullptr;
    const int first=specimen==3 ? 395 : 1350, last=specimen==3 ? 425 : 1385;
    for (const auto &slice:trace.slices) if (slice.frameIndex>=first && slice.frameIndex<=last) {
      require(slice.selected>=0,"Reported shallow neighborhood lost a frame");
      const auto &pair=slice.pairs[slice.selected];
      for (int side=0;side<2;++side) {
        require(pair.profileCurvatures[side]>0,"Measured shallow fallback still uses a shoulder endpoint");
        if (previous) require(std::hypot(pair.points[side].y-previous->points[side].y,
                                        pair.points[side].z-previous->points[side].z)<1.f,
                              "Measured shallow lower turns jump between frames");
      }
      if (slice.frameIndex==401) {
        require(pair.points[0].y>387 && pair.points[0].y<391 && pair.points[0].z<-606.6f &&
                pair.points[1].y>398 && pair.points[1].y<401 && pair.points[1].z<-606.1f,
                "Frame 401 did not select the two observed lower junctions");
      }
      if (slice.frameIndex==1360) {
        require(pair.points[0].y>388 && pair.points[0].y<390 && pair.points[0].z<-610.1f &&
                pair.points[1].y>396 && pair.points[1].y<398 && pair.points[1].z<-610.1f,
                "Specimen 2 shallow groove did not select its measured lower junctions");
      }
      previous=&pair;++checked;
    }
    require(checked==last-first+1,"Measured shallow regression did not cover its neighborhood");
  }
  }
}

void testShortWallFallback() {
  weld::CloudT::Ptr cloud(new weld::CloudT);
  std::vector<weld::ScanFrameRange> frames;
  for (int frame=0;frame<24;++frame) {
    const auto begin=cloud->size();
    for (int j=-200;j<=200;++j) {
      const float y=.1f*j;
      const float z=std::min(1.2f,std::max(0.f,2.f*(std::fabs(y)-8.f)));
      cloud->push_back({frame*.5f,y,z+.01f*frame});
    }
    frames.push_back({begin,cloud->size()-begin,frame});
  }
  weld::DetectorConfig config;
  config.extractMainBody=false; config.outlierMeanK=0; config.voxelSize=0;
  config.trimUnstableEndpoints=false;
  config.curveLineFallback.enabled=false;
  const auto direct=weld::SeamDetector().detect(cloud,config,nullptr,nullptr,nullptr,&frames);
  config.curveLineFallback.enabled=true;
  weld::DebugTrace trace;
  const auto recovered=weld::SeamDetector().detect(cloud,config,nullptr,nullptr,&trace,&frames);
  if (!recovered.ok() || recovered.leftSeam->size()<=direct.leftSeam->size()) {
    std::cerr << "short wall direct=" << direct.leftSeam->size() << " fallback=" << recovered.leftSeam->size()
              << " summary=" << recovered.summary << '\n';
    if (!trace.slices.empty()) for (const auto &gate:trace.slices.front().gates)
      std::cerr << gate.first << '=' << gate.second << '\n';
  }
  require(recovered.ok() && recovered.leftSeam->size()>direct.leftSeam->size(),
          "Short-wall curvature fallback did not recover missing sections");
  bool used=false;
  for (const auto &slice:trace.slices) for (const auto &pair:slice.pairs) {
    used=used || pair.source.find("curve_line_fallback/")==0;
    require(pair.hasBottomMidpoint && pair.points[0].y<pair.bottomMidpointU &&
            pair.points[1].y>pair.bottomMidpointU,"Curvature fallback bypassed the bottom split");
    require(pair.points[0].y<-4 && pair.points[1].y>4,
            "Short-wall recovery selected a one-sided residual instead of the true groove");
  }
  require(used,"Recovered candidates have no fallback provenance");
  const auto silent=weld::SeamDetector().detect(cloud,config,nullptr,nullptr,nullptr,&frames);
  require(silent.leftSeam->size()==recovered.leftSeam->size(),"Debug mode changed fallback coverage");
  for (std::size_t i=0;i<silent.leftSeam->size();++i)
    require((weld::toVector((*silent.leftSeam)[i])-weld::toVector((*recovered.leftSeam)[i])).norm()<1.e-6f,
            "Debug mode changed fallback coordinates");
  config.curvatureSearch.windowLength=-1; // unrelated standalone settings
  require(weld::SeamDetector().detect(cloud,config,nullptr,nullptr,nullptr,&frames).ok(),
          "Curve-line fallback read standalone curvature settings");
  config.curveLineFallback.splineRadius=0;
  require(!weld::SeamDetector().detect(cloud,config,nullptr,nullptr,nullptr,&frames).ok(),
          "Invalid active fallback settings were accepted");
  config.curveLineFallback.splineRadius=8;
  config.interiorOpeningDirection=-1;
  require(!weld::SeamDetector().detect(cloud,config,nullptr,nullptr,nullptr,&frames).ok(),
          "Fallback accepted a groove with the wrong opening direction");
}

void testFloorSideTurnBoundary() {
  {
    std::vector<weld::ProfileCurvatureSample> lobe(5);
    for (int i=0;i<5;++i) {
      lobe[i].point={.1f*i,.9f-.1f*i}; lobe[i].curvature=1; lobe[i].valid=true;
    }
    weld::ProfileTurnRegion region;
    region.begin=0; region.end=4; region.representative=2;
    const Eigen::Vector2f normal(0,1), foot(.4f,.5f);
    const int chosen=weld::supportedProfileTurnRepresentative(lobe,region,normal,0,foot,.61f,2);
    require(chosen==3,"A barely unsupported median discarded the entire supported turn");
    require(weld::supportedProfileTurnRepresentative(lobe,region,normal,0,foot,.71f,2)==2,
            "A supported nominal representative was moved");
    region.end=2;
    require(weld::supportedProfileTurnRepresentative(lobe,region,normal,0,foot,.61f,2)==-1,
            "Turn representative escaped its own curvature region");
    region.end=4; lobe[3].valid=false; lobe[4].curvature=-1;
    require(weld::supportedProfileTurnRepresentative(lobe,region,normal,0,foot,.61f,2)==-1,
            "Turn representative crossed invalid samples or changed curvature sign");
  }
  std::vector<weld::ProfileCurvatureSample> samples(5);
  const float v[] = {6, 4, 1, 4, 6};
  for (int i = 0; i < 5; ++i) {
    samples[i].point = Eigen::Vector2f(float(i - 2), v[i]);
    samples[i].valid = true;
  }
  weld::ProfileTurnRegion left;
  left.begin = 0; left.end = 2; left.representative = 1;
  weld::ProfileTurnRegion right;
  right.begin = 2; right.end = 4; right.representative = 3;
  require(weld::profileTurnFloorBoundary(samples, left, 1) == 2 &&
          weld::profileTurnFloorBoundary(samples, right, 1) == 2,
          "Positive-V trough did not choose the lower floor-side endpoints");

  for (auto &sample : samples) sample.point.y() = -sample.point.y();
  require(weld::profileTurnFloorBoundary(samples, left, -1) == 2 &&
          weld::profileTurnFloorBoundary(samples, right, -1) == 2,
          "Negative-V trough did not choose the upper floor-side endpoints");
}

float roughness(const weld::CloudT &cloud) {
  float total = 0.0f;
  for (std::size_t i = 1; i + 1 < cloud.size(); ++i) {
    const Eigen::Vector3f second = weld::toVector(cloud[i + 1]) -
                                   2.0f * weld::toVector(cloud[i]) +
                                   weld::toVector(cloud[i - 1]);
    total += second.tail<2>().norm();
  }
  return total;
}

weld::SeamResult fixture() {
  weld::SeamResult result;
  for (int i = 0; i < 31; ++i) {
    const float alternating = i % 2 ? 0.35f : -0.35f;
    result.leftSeam->push_back(
        weld::PointT(float(i), -10.0f + alternating,
                     0.18f * float(i) + alternating));
    result.rightSeam->push_back(
        weld::PointT(float(i), 10.0f - alternating,
                     0.18f * float(i) - alternating));
    result.leftPointConfidence.push_back(0.9f);
    result.rightPointConfidence.push_back(0.9f);
    weld::SeamPointQuality left, right;
    left.sourceConfidence = right.sourceConfidence = 0.9f;
    result.leftPointQuality.push_back(left);
    result.rightPointQuality.push_back(right);
  }
  (*result.leftSeam)[15].z += 4.0f;
  result.leftSegmentStarts = {0};
  result.rightSegmentStarts = {0};
  return result;
}

struct ProfileFixture {
  weld::CloudT::Ptr cloud{new weld::CloudT};
  std::vector<weld::ScanFrameRange> frames;
};

ProfileFixture profileFixture() {
  ProfileFixture fixture;
  for (int frame = 0; frame < 24; ++frame) {
    const std::size_t begin = fixture.cloud->size();
    for (int sample = 0; sample <= 466; ++sample) {
      const float y = -35.0f + 0.15f * static_cast<float>(sample);
      const float absoluteY = std::fabs(y);
      const float z = absoluteY < 12.0f
                          ? 20.0f + 0.2f * std::cos(y / 5.0f)
                          : absoluteY < 22.0f
                                ? 2.0f * (22.0f - absoluteY)
                                : 0.0f;
      fixture.cloud->push_back(
          weld::PointT(2.0f * frame, y, z + 0.04f * frame));
    }
    fixture.frames.push_back(
        {begin, fixture.cloud->size() - begin, frame});
  }
  return fixture;
}

void testMultiFrameProjection() {
  const weld::DetectorConfig defaults;
  require(defaults.profileMergeFrameCount == 5 && defaults.mergeProfileFrames,
          "Default projection is not five frames");
  require(defaults.profileTurnLocation == weld::ProfileTurnLocation::RegionMedian,
          "Candidate localization does not default to the region median");
  require(defaults.profileLocalizationMethod ==
              weld::ProfileLocalizationMethod::CurveLineIntersection,
          "Profile localization does not default to curve-line intersections");
  require(std::fabs(defaults.centralBeadMinHeight - .8f) < 1.e-6f &&
              std::fabs(defaults.centralBeadMinWidth - 2.f) < 1.e-6f,
          "Central-bead classification defaults changed");
  require(!defaults.useSurfaceRefinement,
          "Surface refinement is enabled by default");
  require(!defaults.adaptiveSeamCount,
          "Seam count does not default to two tracks");
  require(defaults.debugStartSlice == 1,
          "Debug stepping does not default to the first frame");
  require(defaults.extractMainBody,
          "Main-body extraction is not enabled by default");
  // Tilted scan planes distinguish true orthogonal projection from simply
  // replacing the longitudinal coordinate. Exercise all longitudinal axes.
  for (int axis = 0; axis < 3; ++axis) {
    const int u = axis == 0 ? 1 : 0, v = axis == 2 ? 1 : 2;
    ProfileFixture fixture;
    const std::vector<int> ids{0,1,2,3,4,5,6,10,11,12,13,14};
    for (int id : ids) {
      const auto begin = fixture.cloud->size();
      for (int j = -250; j <= 250; ++j) {
        weld::PointT p;
        p.data[u] = j * .1f;
        p.data[v] = std::min(15.f, std::max(0.f, 2 * (std::fabs(p.data[u]) - 8.f)));
        p.data[axis] = id * .8f + .15f * p.data[u] + .07f * p.data[v];
        fixture.cloud->push_back(p);
      }
      fixture.frames.push_back({begin,fixture.cloud->size()-begin,id});
    }
    const weld::CloudT original(*fixture.cloud);
    weld::DetectorConfig config;
    config.longitudinalAxis = static_cast<weld::Axis>(axis);
    config.extractMainBody = false; config.voxelSize = 0; config.outlierMeanK = 0;
    config.useSurfaceRefinement = false; config.finalSmoothingEnabled = false;
    config.trimUnstableEndpoints = false; config.adaptiveSeamCount = false;
    weld::DebugTrace merged;
    bool sawOverviewCandidates = false;
    if (axis == 0) {
      config.debugFlag = true;
      config.debugFrameBatchSize = 1;
      merged.onStep = [&](const weld::DebugStep &step) {
        if (step.title.compare(0, 1, "F") == 0 &&
            step.overviewCandidates[0] && step.overviewCandidates[1] &&
            !step.overviewCandidates[0]->empty() &&
            !step.overviewCandidates[1]->empty())
          sawOverviewCandidates = true;
        return true;
      };
    }
    const auto result = weld::SeamDetector().detect(fixture.cloud,config,nullptr,nullptr,&merged,&fixture.frames);
    require(result.ok(), "Projected tilted groove detection failed");
    require(axis != 0 || sawOverviewCandidates,
            "Frame debug steps did not publish accumulated overview candidates");
    require(merged.slices.size() == ids.size(), "Merge dropped output sections");
    for (std::size_t f = 0; f < ids.size(); ++f) {
      const auto &slice = merged.slices[f];
      std::vector<std::size_t> neighbors;
      for (std::size_t k = 0; k < ids.size(); ++k)
        if (std::abs(ids[k]-ids[f]) <= 2 && ((ids[k] < 7) == (ids[f] < 7))) neighbors.push_back(k);
      require(slice.raw.size() == neighbors.size() * 501, "Wrong centered window or missing-frame boundary");
      require(slice.sourceRaw.size() == 501,
              "Overview highlight does not contain exactly the selected source frame");
      for (std::size_t j = 0; j < 501; ++j)
        require((weld::toVector(slice.sourceRaw[j]) -
                 weld::toVector(original[fixture.frames[f].begin + j])).norm() == 0,
                "Overview highlight was not kept in original point-cloud coordinates");
      require(slice.gates.at("profile_merged_frames") == neighbors.size(), "Wrong merged-frame diagnostic");
      require(!slice.fittedProfile.empty(), "Merged profile was not fitted with B-spline");
      for (const auto &side : slice.sourceCandidates)
        for (const auto &candidate : side)
          require(std::fabs(candidate.data[axis]-.15f*candidate.data[u]-
                                .07f*candidate.data[v]-.8f*ids[f]) < 2e-4f,
                  "Overview candidate was not restored to its original scan plane");
      std::size_t index = 0;
      // Plane equation: axis - .15 U - .07 V = .8 frameId.
      for (auto k : neighbors) for (std::size_t j = 0; j < 501; ++j) {
        const auto &p = original[fixture.frames[k].begin+j];
        const float distance = (p.data[axis]-.15f*p.data[u]-.07f*p.data[v]-.8f*ids[f]) / (1+.15f*.15f+.07f*.07f);
        const auto &q = slice.raw[index++];
        require(std::fabs(q.data[u]-(p.data[u]+.15f*distance)) < 2e-4f &&
                std::fabs(q.data[v]-(p.data[v]+.07f*distance)) < 2e-4f,
                "Multi-frame points were not orthogonally projected onto target plane");
        require(q.data[axis] == slice.station, "Projected input has mixed section stations");
      }
    }
    for (std::size_t i=0;i<result.leftSeam->size();++i) {
      const auto &p=(*result.leftSeam)[i];
      require(std::fabs(p.data[axis]-.15f*p.data[u]-.07f*p.data[v]-.8f*result.leftFrameIndices[i]) < 2e-4f,
              "Output was not restored to its original frame plane");
    }
    config.mergeProfileFrames = false;
    weld::DebugTrace single;
    const auto singleResult = weld::SeamDetector().detect(fixture.cloud,config,nullptr,nullptr,&single,&fixture.frames);
    config.mergeProfileFrames = true; config.profileMergeFrameCount = 1;
    const auto oneResult = weld::SeamDetector().detect(fixture.cloud,config,nullptr,nullptr,nullptr,&fixture.frames);
    require(singleResult.ok() && oneResult.ok() && singleResult.leftFrameIndices == oneResult.leftFrameIndices,
            "Disabled merge and one-frame mode differ");
    for (std::size_t i=0;i<singleResult.leftSeam->size();++i)
      require((weld::toVector((*singleResult.leftSeam)[i])-weld::toVector((*oneResult.leftSeam)[i])).norm()==0,
              "One-frame option changed the single-frame result");
    for (std::size_t f=0;f<ids.size();++f) {
      require(single.slices[f].raw.size()==501, "Disabled merge still borrowed neighbors");
      for (std::size_t j=0;j<501;++j) {
        const auto &p=original[fixture.frames[f].begin+j], &q=single.slices[f].raw[j];
        require(p.data[u]==q.data[u] && p.data[v]==q.data[v], "Disabled merge altered measured U/V");
      }
    }
    for (std::size_t i=0;i<original.size();++i)
      require((weld::toVector(original[i])-weld::toVector((*fixture.cloud)[i])).norm()==0, "Projection modified input cloud");
    config.profileMergeFrameCount=4;
    require(!weld::SeamDetector().detect(fixture.cloud,config,nullptr,nullptr,nullptr,&fixture.frames).ok(),
            "Even merge window was accepted");
  }
}

void testTurnRegions() {
  // A broad lobe with two nearly equal maxima: reversing a small perturbation
  // should not move its representative from one maximum to the other.
  std::vector<float> peakPositions, medianPositions, turnPositions;
  for (int sign : {-1,1}) {
    std::vector<weld::ProfileCurvatureSample> samples;
    for (int i=-40;i<=40;++i) {
      const float s=.05f*i;
      const float left=std::exp(-std::pow((s+.45f)/.18f,2.f));
      const float right=std::exp(-std::pow((s-.45f)/.18f,2.f));
      weld::ProfileCurvatureSample p; p.point={s,0}; p.valid=true;
      p.curvature=std::fabs(s)<=1.2f ? .8f+.18f*(left+right)+sign*.025f*(right-left) : 0;
      samples.push_back(p);
    }
    for (bool weighted : {false,true}) {
      const auto regions=weld::locateProfileTurnRegions(samples,weighted,.6f,4.f);
      require(regions.size()==1,"One high-curvature lobe produced multiple representatives");
      const auto &r=regions.front();
      require(r.arcLength<=4.f && r.turningAngle>1,"Turn-region size/integral is invalid");
      (weighted ? turnPositions : medianPositions).push_back(samples[r.representative].point.x());
      if (weighted) peakPositions.push_back(samples[r.peak].point.x());
    }
    samples[40].valid=false;
    const auto broken=weld::locateProfileTurnRegions(samples,true,.6f,4);
    require(broken.size()==2,"Invalid samples did not split a turn region");
    for (const auto &r:broken) require(r.end<40 || r.begin>40,"Region bridged an invalid sample");
    samples[40].valid=true;
    auto valley=samples;
    for (int i=36;i<=44;++i) valley[i].curvature=.1f;
    const auto separate=weld::locateProfileTurnRegions(valley,true,.6f,4);
    require(separate.size()==2,"Low-curvature valley did not separate two turns");
    for (const auto &r:separate) require(r.end<36 || r.begin>44,"Region crossed the curvature valley");
    for (std::size_t i=40;i<samples.size();++i) samples[i].curvature=-samples[i].curvature;
    const auto reversed=weld::locateProfileTurnRegions(samples,true,.6f,4);
    require(reversed.size()==2,"Opposite turns were merged");
    for (const auto &r:reversed)
      require(samples[r.begin].curvature*samples[r.end].curvature>0,"Region crosses a curvature sign change");
    const auto bounded=weld::locateProfileTurnRegions(samples,true,.6f,.5f);
    for (const auto &r:bounded) require(r.arcLength<=.50001f,"Turn region exceeded its physical length bound");
  }
  require(std::fabs(peakPositions[1]-peakPositions[0])>.8f,"Peak-switch fixture did not switch maxima");
  require(std::fabs(medianPositions[1]-medianPositions[0])<.1f &&
          std::fabs(turnPositions[1]-turnPositions[0])<.1f,"Region localization retained a peak-height switch");
  // Known integral: for k(s)=1+s on [0,1], half the turn occurs at
  // sqrt(2.5)-1, while the arc median is .5. Unequal sampling must not bias it.
  for (bool dense : {false,true}) {
    std::vector<weld::ProfileCurvatureSample> samples;
    auto add=[&](float s,float k) { weld::ProfileCurvatureSample p; p.point={s,0};p.curvature=k;p.valid=true;samples.push_back(p); };
    add(-.01f,0);
    for (int i=0;i<=100;++i) {
      const float s=i*.01f;
      add(s,1+s);
      if (dense && i<25) for (int j=1;j<5;++j) add(s+.002f*j,1+s+.002f*j);
    }
    add(1.01f,0);
    for (bool weighted : {false,true}) {
      const auto regions=weld::locateProfileTurnRegions(samples,weighted,.5f,4);
      require(regions.size()==1,"Known monotone curvature lobe was lost");
      const float expected=weighted ? std::sqrt(2.5f)-1 : .5f;
      require(std::fabs(samples[regions[0].representative].point.x()-expected)<.011f,
              "Region median depends on sampling density or integrates curvature incorrectly");
    }
    require(weld::locateProfileTurnRegions(samples,true,0,4).empty(),"Invalid region threshold was accepted");
  }
}

void testProfileCurvature() {
  // Unequal sampling on the two faces must not drag a sharp knee onto the
  // dense floor. Check location, not just smoothness or detection count.
  std::vector<Eigen::Vector2f> knee;
  for (int i=-100;i<=0;++i) knee.emplace_back(i*.1f,0.f);
  for (int i=1;i<=25;++i) knee.emplace_back(i*.1f,i*.5f);
  const auto turns=weld::measureProfileTurns(knee,3.6f,.3f);
  std::vector<bool> eligible(knee.size(),false);
  for (std::size_t i=0;i<eligible.size();++i) eligible[i]=i%3==0;
  const auto sparseTurns=weld::measureProfileTurns(knee,3.6f,.3f,&eligible);
  std::size_t retained=0;
  for (const auto &turn:turns) if (eligible[turn.index]) {
    require(retained<sparseTurns.size(),"Peak-only tangent fitting lost an eligible probe");
    const auto &sparse=sparseTurns[retained++];
    require(turn.index==sparse.index && turn.angle==sparse.angle && turn.curvature==sparse.curvature &&
            turn.rms==sparse.rms && (turn.point-sparse.point).squaredNorm()==0 &&
            (turn.before-sparse.before).squaredNorm()==0 && (turn.after-sparse.after).squaredNorm()==0,
            "Peak-only tangent fitting changed a retained probe");
  }
  require(retained==sparseTurns.size() && retained<turns.size(),"Peak-only tangent fitting ignored its mask");
  require(!turns.empty(), "Measured turn search produced no probes");
  const auto strongest=std::max_element(turns.begin(),turns.end(),
      [](const weld::ProfileTurnSample &a,const weld::ProfileTurnSample &b) {
        return std::fabs(a.curvature)<std::fabs(b.curvature);
      });
  require(strongest->point.norm()<.1f && strongest->angle>1.f,
          "Unequal wall/floor sampling displaced the observed sharp turn");
  std::vector<Eigen::Vector2f> vertical, circle, reversed;
  for (int i=0; i<101; ++i) {
    const float t = .004f*i*i; // Nonuniform sampling, including a vertical tangent.
    vertical.emplace_back(30.f, t);
    const float a = .01f*i;
    circle.emplace_back(20*std::cos(a), 20*std::sin(a));
    reversed.emplace_back(20*std::cos(a), -20*std::sin(a));
  }
  for (int shape=0;shape<3;++shape) {
    const auto fitted=weld::fitParametricProfileSpline(shape==0?vertical:shape==1?circle:reversed,2,.1f);
    require(fitted.size()>30,"Full profile B-spline produced no curve");
    const float expected=shape==0?0.f:shape==1?.05f:-.05f;
    for (std::size_t i=fitted.size()/4;i<3*fitted.size()/4;++i)
      require(fitted[i].valid && std::fabs(fitted[i].curvature-expected)<.005f,
              "Global B-spline analytic curvature disagrees with known geometry");
  }
  const auto line = weld::fitProfileCurvature(vertical, 8);
  require(line[50].valid && std::fabs(line[50].curvature)<1.e-5f,
          "Vertical straight wall acquired curvature");
  const auto arc = weld::fitProfileCurvature(circle, 3);
  const auto reverse = weld::fitProfileCurvature(reversed, 3);
  require(arc[50].valid && std::fabs(arc[50].curvature-.05f)<.002f,
          "B-spline curvature missed the analytic circle");
  require(reverse[50].valid && std::fabs(reverse[50].curvature+.05f)<.002f,
          "B-spline curvature lost its turn direction");
  require(!arc.front().valid && !arc.back().valid,
          "Unsupported curve endpoints became curvature candidates");
}

void testFinalCurve() {
  weld::DetectorConfig config;
  weld::SeamResult result;
  weld::CloudT truth;
  for (int i=0; i<601; ++i) {
    const float x=.3f*i;
    const float y=-10+.01f*x;
    const float z=.001f*x*x+3*std::sin(x/25);
    truth.push_back({x,y,z});
    const float spike = i%53==26 ? 8.f : (i>=250 && i<256 ? 6.f : 0.f);
    const float noise = (i%2 ? .4f : -.4f) + spike;
    result.leftSeam->push_back({x,y+noise,z+noise});
    result.rightSeam->push_back({x,y+20-noise,z-noise});
    result.leftPointConfidence.push_back(.9f);
    result.rightPointConfidence.push_back(.9f);
    result.leftFrameIndices.push_back(i+100);
    result.rightFrameIndices.push_back(i+100);
  }
  const weld::CloudT before(*result.leftSeam);
  const auto ids = result.leftFrameIndices;
  weld::applyFinalCurveSmoothing(result, config);
  require(result.leftFinalCurve.applied, "Final continuous curve was not fitted");
  double oldError=0, newError=0;
  for (std::size_t i=0; i<truth.size(); ++i) {
    oldError += (weld::toVector(before[i])-weld::toVector(truth[i])).squaredNorm();
    newError += (weld::toVector((*result.leftSeam)[i])-weld::toVector(truth[i])).squaredNorm();
    require((*result.leftSeam)[i].x == truth[i].x, "Final fit changed scan station");
    require((*result.rightSeam)[i].y-(*result.leftSeam)[i].y >= config.minSeamSeparation-1.e-4f,
            "Final fit violated seam separation");
  }
  require(newError < .1*oldError, "Final fit failed to recover the known curved seam");
  require(roughness(*result.leftSeam) < .05f*roughness(before), "Final fit retained alternating spikes");
  require(result.leftFrameIndices == ids, "Final fit changed frame identity");
  config.finalSmoothingEnabled=false;
  const weld::CloudT disabled(*result.leftSeam);
  weld::applyFinalCurveSmoothing(result,config);
  for (std::size_t i=0; i<disabled.size(); ++i)
    require((weld::toVector(disabled[i])-weld::toVector((*result.leftSeam)[i])).norm()==0,
            "Disabled final fit changed coordinates");

  // Explicit segment breaks and long missing runs must not pull separate
  // sloped seams toward each other. Nonuniform sampling preserves a line.
  config.finalSmoothingEnabled=true;
  for (bool explicitBreak : {false,true}) {
    weld::SeamResult split;
    for (int i=0; i<40; ++i) {
      const float x=.03f*i*i + (!explicitBreak && i>=20 ? 100.f : 0.f);
      const float offset=i>=20 ? 35.f : 0.f;
      split.leftSeam->push_back({x,-10+offset+.1f*x,2*x+offset});
      split.rightSeam->push_back({x,10+offset+.1f*x,2*x+offset});
    }
    if (explicitBreak) split.leftSegmentStarts=split.rightSegmentStarts={0,20};
    const weld::CloudT original(*split.leftSeam);
    weld::applyFinalCurveSmoothing(split,config);
    for (std::size_t i=0;i<original.size();++i)
      require((weld::toVector(original[i])-weld::toVector((*split.leftSeam)[i])).norm()<.002f,
              "Final fit bridged a gap or flattened an actual slope");
  }
}

void requireSameTracks(const weld::SeamResult &a,const weld::SeamResult &b) {
  require(a.ok()==b.ok() && a.leftFrameIndices==b.leftFrameIndices &&
          a.rightFrameIndices==b.rightFrameIndices && a.centerFrameIndices==b.centerFrameIndices,
          "Wall-region pruning changed frame identity");
  require(a.endpointTrimmedStart==b.endpointTrimmedStart && a.endpointTrimmedEnd==b.endpointTrimmedEnd &&
          a.leftPointConfidence==b.leftPointConfidence && a.rightPointConfidence==b.rightPointConfidence &&
          a.centerPointConfidence==b.centerPointConfidence,"Wall-region pruning changed output quality or endpoints");
  const std::array<weld::CloudT::Ptr,3> first{{a.leftSeam,a.rightSeam,a.centerSeam}},second{{b.leftSeam,b.rightSeam,b.centerSeam}};
  for (int side=0;side<3;++side) {
    require(first[side]->size()==second[side]->size(),"Wall-region pruning changed track size");
    for (std::size_t i=0;i<first[side]->size();++i)
      require((weld::toVector((*first[side])[i])-weld::toVector((*second[side])[i])).squaredNorm()==0,
              "Wall-region pruning changed a trajectory coordinate");
  }
}

void testCurvatureTrajectory() {
  ProfileFixture fixture;
  for (int frame=0; frame<60; ++frame) {
    const auto begin = fixture.cloud->size();
    const float jitter = frame%2 ? .22f : -.22f;
    const float shift = .3f*std::sin(frame*.06f);
    for (int j=-350; j<=350; ++j) {
      const float y = j*.1f;
      const float across = y-shift-jitter;
      const float distance=std::fabs(across)-5.5f;
      // Some walls bend gradually before the stronger lower junction. The
      // remote straight wall must not force a candidate onto that upper bend.
      const float height=frame%3==0 ? distance+std::max(0.f,distance-3.f) : 2*distance;
      const float z = std::min(16.f, std::max(0.f,height));
      // Corrugated floor and alternating sampling offsets must not make the
      // final trajectory jump.
      const float floorNoise = std::max(0.f, 1-std::fabs(across)/5.5f) *
                               .6f*std::sin(across*2+frame*.7f);
      fixture.cloud->push_back({frame*.5f, y, z+floorNoise+.02f*frame});
    }
    fixture.frames.push_back({begin,fixture.cloud->size()-begin,frame});
  }
  weld::DetectorConfig config;
  config.extractMainBody = false;
  config.voxelSize = 0;
  config.outlierMeanK = 0;
  config.useSurfaceRefinement = false;
  config.adaptiveSeamCount = false;
  config.profileLocalizationMethod =
      weld::ProfileLocalizationMethod::CurvatureCandidates;
  weld::DebugTrace trace;
  const auto result = weld::SeamDetector().detect(fixture.cloud,config,nullptr,nullptr,&trace,&fixture.frames);
  // Compare raw path candidates, before final smoothing, against the known
  // slowly curving junction, including the stable wall-turn centroid.
  for (auto mode : {weld::ProfileTurnLocation::Peak,weld::ProfileTurnLocation::RegionMedian,weld::ProfileTurnLocation::TurnMedian,weld::ProfileTurnLocation::FloorSideBoundary,weld::ProfileTurnLocation::WallTurnCentroid}) {
    auto compare=config; compare.profileTurnLocation=mode; compare.finalSmoothingEnabled=false;
    weld::DebugTrace candidates;
    const auto trial=weld::SeamDetector().detect(fixture.cloud,compare,nullptr,nullptr,&candidates,&fixture.frames);
    require(trial.ok(),"A turn-location mode lost the synthetic groove");
    double error2=0,jitter2=0; std::size_t points=0,triples=0;
    for (int side=0;side<2;++side) {
      std::vector<Eigen::Vector2f> residuals;
      for (const auto &p:candidates.selectedPath[side]) {
        const float frame=p.x/.5f;
        const Eigen::Vector2f residual(p.y-((side?5.5f:-5.5f)+.3f*std::sin(frame*.06f)),p.z-.02f*frame);
        error2+=residual.squaredNorm();++points;residuals.push_back(residual);
      }
      for (std::size_t i=1;i+1<residuals.size();++i) {
        jitter2+=(residuals[i+1]-2*residuals[i]+residuals[i-1]).squaredNorm();++triples;
      }
    }
    std::cout<<"turn_mode="<<int(mode)<<" points="<<points<<" candidate_rms_error="<<std::sqrt(error2/points)
             <<" candidate_second_difference_rms="<<std::sqrt(jitter2/triples)<<'\n';
    require(points==120 && std::sqrt(error2/points)<1.1,"Turn-region candidate localization lost frames or acquired excessive bias");
  }
  auto exhaustive=config; exhaustive.useWallCandidatePruning=false;
  weld::DebugTrace exhaustiveTrace;
  const auto reference=weld::SeamDetector().detect(fixture.cloud,exhaustive,nullptr,nullptr,&exhaustiveTrace,&fixture.frames);
  requireSameTracks(result,reference);
  std::size_t pruned=0;
  require(trace.slices.size()==exhaustiveTrace.slices.size(),"Wall-region pruning changed sections");
  for (std::size_t i=0;i<trace.slices.size();++i) {
    const auto &a=trace.slices[i]; const auto &b=exhaustiveTrace.slices[i];
    auto gate=a.gates.find("wall_roi_pruned"); if (gate!=a.gates.end()) pruned+=gate->second;
    require(a.pairs.size()==b.pairs.size(),"Wall-region pruning removed a valid candidate pair");
    for (std::size_t j=0;j<a.pairs.size();++j) {
      require(a.pairs[j].source==b.pairs[j].source && a.pairs[j].proposalCost==b.pairs[j].proposalCost,
              "Wall-region pruning changed candidate scoring/order");
      for (int side=0;side<2;++side)
        require((weld::toVector(a.pairs[j].points[side])-weld::toVector(b.pairs[j].points[side])).squaredNorm()==0,
                "Wall-region pruning changed a curvature peak");
    }
  }
  require(pruned>0,"Wall-region optimization did not skip any surface fits");
  // A different bin grid changes the search radius in indices. It must still
  // preserve the exhaustive result, rather than pruning with a fixed bin count.
  auto stable=config; stable.profileSamplingMode=weld::ProfileSamplingMode::StableGrid;
  stable.profileBinSize=.6f;
  const auto stableResult=weld::SeamDetector().detect(fixture.cloud,stable,nullptr,nullptr,nullptr,&fixture.frames);
  stable.useWallCandidatePruning=false;
  const auto stableReference=weld::SeamDetector().detect(fixture.cloud,stable,nullptr,nullptr,nullptr,&fixture.frames);
  requireSameTracks(stableResult,stableReference);
  require(result.ok() && result.leftSeam->size()==60, "Curvature lost corrugated-floor frames");
  require(result.leftFinalCurve.applied && result.rightFinalCurve.applied,
          "Curvature detector did not fit the complete trajectory");
  int curvatureSides = 0;
  for (const auto &slice : trace.slices)
    if (slice.selected>=0) for (float value : slice.pairs[slice.selected].profileCurvatures)
      curvatureSides += std::fabs(value)>1.e-5f;
  for (const auto &slice : trace.slices) for (const auto &pair : slice.pairs) {
    require(pair.source=="spline_curvature", "Unverified fallback was displayed as a curvature candidate");
    for (int side=0;side<2;++side) {
      const auto &wall=pair.candidateWalls[side];
      const auto &p=pair.points[side];
      const float wallDistance=std::fabs(-wall.dv*(p.y-wall.u)+wall.du*(p.z-wall.v));
      require(wallDistance<=config.profileWallTolerance+1.e-4f,
              "A displayed candidate escaped the hard wall band");
      const float center=.3f*std::sin(slice.frameIndex*.06f)+(slice.frameIndex%2 ? .22f : -.22f);
      const float localSlope=slice.frameIndex%3==0 ? 1.f : 2.f;
      const float truthDistance=std::fabs(p.z-.02f*slice.frameIndex-
           localSlope*((side ? 1.f : -1.f)*(p.y-center)-5.5f))/std::sqrt(1+localSlope*localSlope);
      if (truthDistance>=.7f) std::cerr << "wall error frame=" << slice.frameIndex << " side=" << side
          << " point=" << p.y << ',' << p.z << " center=" << center << " distance=" << truthDistance << '\n';
      require(truthDistance<.7f, "A candidate left the analytic wall near its junction");
      require(std::fabs(p.y-center-(side ? 5.5f : -5.5f))<1.1f,
              "A displayed candidate passed beyond the known local turn");
    }
  }
  require(curvatureSides>60, "Curvature mode silently fell back on most sides");
  for (int side=0; side<2; ++side) {
    const auto &cloud = side ? *result.rightSeam : *result.leftSeam;
    for (std::size_t i=0; i<cloud.size(); ++i) {
      const float expected = (side ? 5.5f : -5.5f)+.3f*std::sin(i*.06f);
      require(std::fabs(cloud[i].y-expected)<1.5f, "Curvature path left the known wall junction");
      if (i) require(std::hypot(cloud[i].y-cloud[i-1].y,cloud[i].z-cloud[i-1].z)<.3f,
                     "Final curvature trajectory retained a frame-to-frame jump");
    }
  }
}

void testNarrowWallSpikes() {
  std::vector<Eigen::Vector2f> clean,spiked;
  for (int i=0;i<=200;++i) {
    const float x=-10+.1f*i;
    clean.emplace_back(x,std::max(0.f,2*(std::fabs(x)-3)));
  }
  spiked=clean; spiked[30].y()-=3;
  const auto cleanMask=weld::wallProfileSpikeMask(clean,4.8f,2.4f,.6f);
  const auto spikeMask=weld::wallProfileSpikeMask(spiked,4.8f,2.4f,.6f);
  require(std::count(cleanMask.begin(),cleanMask.end(),true)==0,
          "Wall filtering erased a real wall-floor corner");
  require(spikeMask[30],"An isolated wall outlier was not removed before spline fitting");
  require(!spikeMask[70] && !spikeMask[130],"Spike removal erased the true groove feet");
  ProfileFixture fixture;
  for (int frame=0;frame<50;++frame) {
    const auto begin=fixture.cloud->size();
    const bool spiked=frame>=20 && frame<30;
    const auto *profile=spiked ? kSpikedWallProfile : kCleanWallProfile;
    const auto count=spiked ? sizeof(kSpikedWallProfile)/sizeof(kSpikedWallProfile[0]) :
                             sizeof(kCleanWallProfile)/sizeof(kCleanWallProfile[0]);
    for (std::size_t j=0;j<count;++j)
      fixture.cloud->push_back({frame*.3f,profile[j][0],profile[j][1]});
    fixture.frames.push_back({begin,fixture.cloud->size()-begin,frame});
  }
  weld::DetectorConfig config;
  config.extractMainBody=false; config.voxelSize=0; config.outlierMeanK=0;
  config.useSurfaceRefinement=false; config.trimUnstableEndpoints=false;
  config.adaptiveSeamCount=false;
  const auto result=weld::SeamDetector().detect(fixture.cloud,config,nullptr,nullptr,nullptr,&fixture.frames);
  require(result.ok() && result.leftSeam->size()==50 && result.rightSeam->size()==50,
          "A narrow wall spike split a physically continuous groove into separate runs");
  for (const auto &p:*result.leftSeam)
    require(p.y>415 && p.y<422,"Wall spike replaced the actual left wall foot");
  require(result.leftFinalCurve.maximumWallResidual<=config.profileWallTolerance+2.e-4f &&
          result.rightFinalCurve.maximumWallResidual<=config.profileWallTolerance+2.e-4f,
          "Recovered wall-spike frames escaped the observed wall constraints");
  auto omitFrames=[&](int first,int last) {
    ProfileFixture subset;
    for (const auto &frame:fixture.frames) {
      if (frame.frameIndex>=first && frame.frameIndex<last) continue;
      const auto begin=subset.cloud->size();
      for (std::size_t j=frame.begin;j<frame.begin+frame.count;++j) subset.cloud->push_back((*fixture.cloud)[j]);
      subset.frames.push_back({begin,subset.cloud->size()-begin,frame.frameIndex});
    }
    return subset;
  };
  auto shortGap=omitFrames(5,9);
  const auto joined=weld::SeamDetector().detect(shortGap.cloud,config,nullptr,nullptr,nullptr,&shortGap.frames);
  if (joined.leftSeam->size()!=46 || joined.leftFrameIndices.front()!=0 || joined.leftFrameIndices.back()!=49)
    std::cerr<<"Short gap: "<<joined.summary<<" points="<<joined.leftSeam->size()<<'\n';
  require(joined.leftSeam->size()==46 && joined.leftFrameIndices.front()==0 && joined.leftFrameIndices.back()==49,
          "Compatible measured walls were discarded after a short physical scan gap");
  require(std::find(joined.leftFrameIndices.begin(),joined.leftFrameIndices.end(),6)==joined.leftFrameIndices.end(),
          "Short-gap handling invented an absent input frame");
  ProfileFixture rotatingWalls;
  for (int frame=0;frame<30;++frame) {
    if (frame>=5 && frame<9) continue;
    const auto begin=rotatingWalls.cloud->size();
    const float slope=frame<5 ? .7f : 2.f;
    for (int sample=-200;sample<=200;++sample) {
      const float y=.1f*sample;
      rotatingWalls.cloud->push_back(
          {frame*.3f,y,slope*std::max(0.f,std::fabs(y)-5.5f)});
    }
    rotatingWalls.frames.push_back(
        {begin,rotatingWalls.cloud->size()-begin,frame});
  }
  const auto rotatedJoin=weld::SeamDetector().detect(
      rotatingWalls.cloud,config,nullptr,nullptr,nullptr,&rotatingWalls.frames);
  require(rotatedJoin.ok() && rotatedJoin.leftSeam->size()==26 &&
              rotatedJoin.leftFrameIndices.front()==0 &&
              rotatedJoin.leftFrameIndices.back()==29,
          "A short gap split stable junctions when local wall tangents rotated");
  auto longGap=omitFrames(15,35);
  const auto separated=weld::SeamDetector().detect(longGap.cloud,config,nullptr,nullptr,nullptr,&longGap.frames);
  require(separated.leftSeam->size()<30,"A long unsupported gap was joined as one continuous groove");
}

void testIndependentWalls() {
  ProfileFixture fixture;
  for (int frame=0;frame<60;++frame) {
    const auto begin=fixture.cloud->size();
    for (int j=-450;j<=450;++j) {
      const float y=j*.1f, d=std::max(0.f,std::fabs(y)-15.f);
      const float wall=std::min(30.f,d<6 ? d*d/6 : 2*d-6);
      const float depth=frame>=15 && frame<45 ? 4.f : .2f;
      const float pit=-depth*std::max(0.f,1-std::fabs(y+5)/3);
      const float bead=frame>=25 && frame<=28 ? .9f*std::exp(-std::pow((y-18)/.25f,2)) : 0;
      fixture.cloud->push_back({frame*.5f,y,wall+pit+bead});
    }
    fixture.frames.push_back({begin,fixture.cloud->size()-begin,frame});
  }
  weld::DetectorConfig config;
  config.extractMainBody=false; config.voxelSize=0; config.outlierMeanK=0;
  config.useSurfaceRefinement=false; config.trimUnstableEndpoints=false;
  weld::DebugTrace directTrace;
  const auto result=weld::SeamDetector().detect(
      fixture.cloud,config,nullptr,nullptr,&directTrace,&fixture.frames);
  auto exhaustive=config; exhaustive.useWallCandidatePruning=false;
  requireSameTracks(result,weld::SeamDetector().detect(fixture.cloud,exhaustive,nullptr,nullptr,nullptr,&fixture.frames));
  require(result.ok() && result.leftSeam->size()==60,"Independent walls lost the groove beside a floor pit");
  bool usedCurveLine=false;
  for (const auto &slice:directTrace.slices)
    for (const auto &pair:slice.pairs)
      usedCurveLine = usedCurveLine ||
                      pair.source == "curve_line_flat_floor" ||
                      pair.source == "curve_line_wavy_floor" ||
                      pair.source == "central_bead_valleys";
  require(usedCurveLine,"Default localization did not use curve-line intersections");
  require(result.leftFinalCurve.wallConstrainedPoints==60 && result.rightFinalCurve.wallConstrainedPoints==60,
          "The final fit lost the independent wall constraints");
  require(result.leftFinalCurve.maximumWallResidual<=config.profileWallTolerance+2.e-4f &&
          result.rightFinalCurve.maximumWallResidual<=config.profileWallTolerance+2.e-4f,
          "Reported final wall residual exceeds the allowed strip");
  for (const auto &p:*result.leftSeam) require(p.y < -12,"A floor pit replaced the independently observed left wall");
  for (const auto &p:*result.rightSeam) require(p.y > 12,"A floor pit replaced the independently observed right wall");

  // Observations deliberately favor the interior. The final fit must retain
  // both oblique wall strips and the width constraint, in every coordinate axis.
  for (const auto axis:{weld::Axis::X,weld::Axis::Y,weld::Axis::Z}) {
    weld::SeamResult r; config.longitudinalAxis=axis;
    const std::array<int,2> axes=axis==weld::Axis::X ? std::array<int,2>{{1,2}} :
        axis==weld::Axis::Y ? std::array<int,2>{{0,2}} : std::array<int,2>{{0,1}};
    for (int i=0;i<101;++i) for (int side=0;side<2;++side) {
      const float station=.5f*i, u=(side ? 10.f:-10.f)+.02f*station;
      weld::PointT p; p.data[int(axis)]=station;
      p.data[axes[0]]=u+(side ? -4.f:4.f); p.data[axes[1]]=.1f*station+(i%2 ? .3f:-.3f);
      (side?r.rightSeam:r.leftSeam)->push_back(p);
      weld::SeamPointQuality q; q.trackedWall=true;
      q.trackedNormal=Eigen::Vector2f(1,.2f).normalized();
      q.trackedOffset=q.trackedNormal.dot(Eigen::Vector2f(u,.1f*station));
      (side?r.rightPointQuality:r.leftPointQuality).push_back(q);
    }
    weld::applyFinalCurveSmoothing(r,config);
    for (int side=0;side<2;++side) for (int i=0;i<101;++i) {
      const auto &p=(*(side?r.rightSeam:r.leftSeam))[i];
      const auto &q=(side?r.rightPointQuality:r.leftPointQuality)[i];
      require(std::fabs(q.trackedNormal.dot(Eigen::Vector2f(p.data[axes[0]],p.data[axes[1]]))-q.trackedOffset)<=config.profileWallTolerance+2.e-4f,
              "Final curve escaped the independently observed wall");
      require((*r.rightSeam)[i].data[axes[0]]-(*r.leftSeam)[i].data[axes[0]]>=config.minSeamSeparation-2.e-4f,
              "Wall-constrained final curve violated pair width");
    }
  }
}

void testBottomShapeClassification() {
  auto makeFixture = [](int shape) {
    ProfileFixture fixture;
    for (int frame = 0; frame < 14; ++frame) {
      const std::size_t begin = fixture.cloud->size();
      for (int sample = -300; sample <= 300; ++sample) {
        const float y = .1f * sample;
        const float absoluteY = std::fabs(y);
        const float wall = absoluteY > 12.f
                               ? std::min(24.f, 2.f * (absoluteY - 12.f))
                               : 0.f;
        float bottom = 0.f;
        if (shape == 1)
          bottom = 2.6f * std::exp(-.5f * std::pow(y / 2.f, 2.f));
        else if (shape == 2)
          bottom = .45f * std::cos(3.14159265f * y / 3.f);
        fixture.cloud->push_back(
            {frame * .6f, y, wall + bottom + .02f * frame});
      }
      fixture.frames.push_back(
          {begin, fixture.cloud->size() - begin, frame});
    }
    return fixture;
  };
  auto classify = [&](int shape, weld::ProfileLocalizationMethod method,
                      const char *expectedSource,
                      const char *expectedShape) {
    const auto fixture = makeFixture(shape);
    weld::DetectorConfig config;
    config.extractMainBody = false;
    config.voxelSize = 0;
    config.outlierMeanK = 0;
    config.trimUnstableEndpoints = false;
    config.finalSmoothingEnabled = false;
    config.useSurfaceRefinement = false;
    config.profileLocalizationMethod = method;
    weld::DebugTrace trace;
    const auto result = weld::SeamDetector().detect(
        fixture.cloud, config, nullptr, nullptr, &trace, &fixture.frames);
    require(result.ok(), "Bottom-shape fixture was rejected");
    int selectedMatches = 0;
    std::map<std::string, int> selectedSources;
    for (const auto &slice : trace.slices)
      if (slice.selected >= 0 &&
          slice.selected < static_cast<int>(slice.pairs.size())) {
        const auto &source = slice.pairs[slice.selected].source;
        ++selectedSources[source];
        if (source == expectedSource &&
            slice.pairs[slice.selected].profileShape == expectedShape)
          ++selectedMatches;
      }
    if (selectedMatches < 10) {
      std::cerr << "shape=" << shape << " expected=" << expectedSource
                << '/' << expectedShape;
      for (const auto &entry : selectedSources)
        std::cerr << ' ' << entry.first << '=' << entry.second;
      std::cerr << '\n';
    }
    require(selectedMatches >= 10,
            "Curve extraction did not classify the expected bottom shape");
    return std::make_pair(result, trace);
  };

  const auto direct = weld::ProfileLocalizationMethod::CurveLineIntersection;
  const auto curvature = weld::ProfileLocalizationMethod::CurvatureCandidates;
  const auto flat = classify(0, direct, "curve_line_flat_floor", "flat_floor");
  const auto bead = classify(1, direct, "central_bead_valleys", "central_bead");
  const auto wavy = classify(2, direct, "curve_line_wavy_floor", "wavy_floor");
  classify(0, curvature, "spline_curvature", "flat_floor");
  classify(1, curvature, "central_bead_valleys", "central_bead");
  classify(2, curvature, "spline_curvature", "wavy_floor");
  require(!flat.first.leftSeam->empty() && !wavy.first.leftSeam->empty(),
          "Flat or wavy floor lost its outer junctions");
  for (std::size_t i = 0; i < bead.first.leftSeam->size(); ++i) {
    const float left = (*bead.first.leftSeam)[i].y;
    const float right = (*bead.first.rightSeam)[i].y;
    require(left < -1.f && left > -9.f && right > 1.f && right < 9.f,
            "Central bead did not select the two adjacent valley regions");
  }
}

void testEmergingEndpointTrim() {
  // A slowly moving weak feature can still accumulate into a large terminal
  // hook. Exercise candidate generation AND final trimming, not just the
  // standalone trimmer: the old emerging-prefix cap hid this geometry.
  auto make=[](bool bent) {
    ProfileFixture input;
    for (int frame=0;frame<240;++frame) {
      const float x=.5f*frame;
      const float offset=bent ? .05f*std::pow(std::max(0.f,30.f-x),2) : 0.f;
      const auto begin=input.cloud->size();
      for (int j=-140;j<=140;++j) {
        const float y=.25f*j;
        float z=.01f*y;
        if (frame<12)
          z-=(.5f+.08f*frame)*std::max(0.f,std::min(1.f,(8.f-std::fabs(y))/3.f));
        else if (frame<24)
          z+=(.6f+.06f*(frame-12))*std::max(0.f,1.f-std::fabs(y)/5.f);
        else
          z+=std::min(std::min(10.f,.45f*(frame-23)),
                      std::max(0.f,(std::fabs(y)-5.f)*1.5f));
        input.cloud->push_back({x,y+offset,z});
      }
      input.frames.push_back({begin,input.cloud->size()-begin,frame});
    }
    return input;
  };
  weld::DetectorConfig config;
  config.profileLocalizationMethod=weld::ProfileLocalizationMethod::CurvatureCandidates;
  config.interiorOpeningDirection=0; // measured notch -> bead -> groove
  config.extractMainBody=false; config.voxelSize=0; config.outlierMeanK=0;
  config.minimumFeatureDepth=0; // explicitly exercise weak onset recovery
  config.profileMergeFrameCount=5;
  for (bool bent:{false,true}) {
    const auto input=make(bent);
    config.trimUnstableEndpoints=false;
    weld::DebugTrace trace;
    const auto uncut=weld::SeamDetector().detect(input.cloud,config,nullptr,nullptr,&trace,&input.frames);
    require(uncut.ok() && !uncut.leftFrameIndices.empty(),"Emerging endpoint fixture was not detected");
    const auto first=std::find_if(trace.slices.begin(),trace.slices.end(),
        [&](const weld::DebugSlice &s){return s.frameIndex==uncut.leftFrameIndices.front();});
    std::cout << "emerging bent=" << bent << " first=" << uncut.leftFrameIndices.front()
              << " source=" << (first!=trace.slices.end() && first->selected>=0 ?
                  first->pairs[first->selected].source : "missing") << '\n';
    require(first!=trace.slices.end() && first->selected>=0 &&
            first->pairs[first->selected].source=="emerging_profile",
            "Emerging endpoint fixture did not start with a fallback observation");
    config.trimUnstableEndpoints=true;
    const auto trimmed=weld::SeamDetector().detect(input.cloud,config,nullptr,nullptr,nullptr,&input.frames);
    std::cout << "emerging trim=" << trimmed.endpointTrimmedStart << '\n';
    require(trimmed.ok(),"Endpoint trim removed the established groove");
    if (bent) {
      require(trimmed.endpointTrimmedStart>0 &&
              trimmed.leftFrameIndices.front()>uncut.leftFrameIndices.front(),
              "Continuous emerging-prefix hook bypassed endpoint trimming");
    } else {
      require(trimmed.endpointTrimmedStart==0 &&
              trimmed.leftFrameIndices.front()==uncut.leftFrameIndices.front(),
              "Straight shallow onset was removed merely for using fallback candidates");
    }
    require(trimmed.leftFrameIndices==trimmed.rightFrameIndices,
            "Emerging endpoint trim desynchronized left and right frames");
    for (std::size_t i=0;i<trimmed.leftSeam->size();++i) {
      const auto source=i+trimmed.endpointTrimmedStart;
      require((weld::toVector((*trimmed.leftSeam)[i])-weld::toVector((*uncut.leftSeam)[source])).norm()==0 &&
              (weld::toVector((*trimmed.rightSeam)[i])-weld::toVector((*uncut.rightSeam)[source])).norm()==0,
              "Emerging endpoint trim moved retained measured coordinates");
    }
  }
}

void testEndpointTrim() {
  weld::DetectorConfig config;
  auto make=[](int shape) {
    weld::SeamResult r; r.frameBased=true;
    for (int i=0;i<=400;++i) {
      const float x=.5f*i;
      float bend=0;
      if (shape==1) bend=.03f*std::pow(std::max(0.f,25-x),2)+.03f*std::pow(std::max(0.f,x-175),2);
      if (shape==2) bend=.0005f*x*x;
      if (shape==3 && i==20) bend=8;
      if (shape==4 && x>180) bend=8;
      if (shape==5) bend=.03f*std::pow(std::max(0.f,x-140),2);
      if (shape==6) bend=8*std::exp(-std::pow((x-55)/5,2));
      if (shape==7) bend=8/(1+std::exp(-(x-100)/3));
      if (shape==8 && x<20) bend=8;
      r.leftSeam->push_back({x,-10+.1f*x+bend,.02f*x});
      r.rightSeam->push_back({x,10+.1f*x,.02f*x});
      r.leftFrameIndices.push_back(1000+i); r.rightFrameIndices.push_back(1000+i);
      r.leftPointConfidence.push_back(.9f); r.rightPointConfidence.push_back(.8f);
      r.leftPointQuality.emplace_back(); r.rightPointQuality.emplace_back();
    }
    return r;
  };
  for (int shape:{0,2,3,6,7}) {
    auto r=make(shape); weld::trimSeamEndpoints(r,config);
    require(r.leftSeam->size()==401 && r.endpointTrimmedStart==0 && r.endpointTrimmedEnd==0,
            "Endpoint trimming removed a straight/gently curved/noisy stable seam");
  }
  auto bent=make(1); const weld::CloudT original(*bent.leftSeam);
  bent.leftSegmentStarts={0,200}; bent.rightSegmentStarts={0,200};
  weld::trimSeamEndpoints(bent,config);
  require(bent.endpointTrimmedStart>0 && bent.endpointTrimmedEnd>0,"Large hooks survived endpoint trimming");
  require(bent.leftSeam->front().x>=15 && bent.leftSeam->front().x<=28 &&
          bent.leftSeam->back().x>=172 && bent.leftSeam->back().x<=185,"Endpoint cut missed the known bend onset");
  require(bent.leftFrameIndices==bent.rightFrameIndices && bent.leftSeam->size()==bent.rightSeam->size() &&
          bent.leftPointQuality.size()==bent.leftSeam->size() && bent.rightPointConfidence.size()==bent.leftSeam->size(),
          "Paired endpoint trimming lost aligned IDs or metadata");
  require(bent.leftSegmentStarts==std::vector<std::size_t>({0,200-bent.endpointTrimmedStart}),"Trimming lost an interior segment break");
  for (std::size_t i=0;i<bent.leftSeam->size();++i)
    require((weld::toVector((*bent.leftSeam)[i])-weld::toVector(original[i+bent.endpointTrimmedStart])).norm()==0,
            "Endpoint trimming moved retained coordinates");
  auto jumped=make(4); weld::trimSeamEndpoints(jumped,config);
  require(jumped.endpointTrimmedStart==0 && jumped.leftSeam->back().x<=180 && jumped.leftSeam->back().x>=179.5f,
          "Step trimming removed the good prefix or retained the displaced tail");
  auto headJump=make(8); weld::trimSeamEndpoints(headJump,config);
  require(headJump.endpointTrimmedEnd==0 && headJump.leftSeam->front().x>=20 &&
          headJump.leftSeam->front().x<=20.5f,
          "Head step trimming removed the good suffix or retained displaced points");
  auto longTail=make(5); weld::trimSeamEndpoints(longTail,config);
  require(longTail.endpointTrimmedStart==0 && longTail.leftSeam->back().x>=132 && longTail.leftSeam->back().x<=150,
          "Long bad tail was clipped at a fixed terminal percentage rather than its onset");
  auto rotated=make(1); std::vector<float> stations;
  const auto rotatedDecision=weld::classifySeamWidth(rotated,config);
  for (const auto &p:*rotated.leftSeam) stations.push_back(p.x);
  for (auto *cloud:{rotated.leftSeam.get(),rotated.rightSeam.get()}) for (auto &p:*cloud) {
    const float oldX=p.x; p.x=p.y; p.y=-oldX;
  }
  rotated.leftSegmentStarts={0,200}; rotated.rightSegmentStarts={0,200};
  weld::trimSeamEndpoints(rotated,config,&stations,&rotatedDecision);
  require(rotated.endpointTrimmedStart==bent.endpointTrimmedStart && rotated.endpointTrimmedEnd==bent.endpointTrimmedEnd,
          "Trimming changed with world orientation or nonmonotonic world X");
  auto shortCurve=make(1);
  shortCurve.leftSeam->resize(20); shortCurve.rightSeam->resize(20);
  weld::trimSeamEndpoints(shortCurve,config);
  require(shortCurve.leftSeam->size()==20,"Short unsupported trajectory was trimmed");
  auto disabled=make(1); config.trimUnstableEndpoints=false;
  weld::trimSeamEndpoints(disabled,config);
  require(disabled.leftSeam->size()==401,"Disabled endpoint trimming removed points");
}

void testSingleEndpointTrim() {
  weld::DetectorConfig config;
  config.adaptiveSeamCount = true;
  auto make=[](bool centerBend) {
    weld::SeamResult r; r.frameBased=true;
    for (int i=0;i<=400;++i) {
      const float station=.5f*i;
      // A minority of wide frames: each wall turns, but their midpoint stays
      // continuous. This is still one whole-seam single-track decision.
      const float width=8+12/(1+std::exp(-(station-155)/1.5f));
      const float center=.02f*station+(centerBend ? .03f*std::pow(std::max(0.f,station-175),2) : 0);
      r.leftSeam->push_back({station,center-width/2,.01f*station});
      r.rightSeam->push_back({station,center+width/2,.01f*station});
      r.leftFrameIndices.push_back(1000+i); r.rightFrameIndices.push_back(1000+i);
      r.leftPointConfidence.push_back(.9f); r.rightPointConfidence.push_back(.8f);
      r.leftPointQuality.emplace_back(); r.rightPointQuality.emplace_back();
    }
    return r;
  };
  auto stable=make(false);
  const auto decision=weld::classifySeamWidth(stable,config);
  require(decision.mode==1 && decision.singleVotes>decision.doubleVotes,"Single endpoint fixture mode wrong");
  weld::trimSeamEndpoints(stable,config,nullptr,&decision);
  require(stable.endpointTrimmedStart==0 && stable.endpointTrimmedEnd==0 && stable.leftSeam->size()==401,
          "Changing groove width truncated a stable single centerline");
  weld::applySeamWidthRule(stable,config,&decision);
  require(stable.centerSeam->size()==401 && stable.centerFrameIndices.back()==1400,
          "Stable single endpoint or IDs were lost");
  auto dual=make(false); auto dualConfig=config; dualConfig.adaptiveSeamCount=false;
  weld::trimSeamEndpoints(dual,dualConfig);
  require(dual.endpointTrimmedEnd>0,"Double-track wall bend no longer triggers endpoint trimming");
  auto bent=make(true);
  const weld::CloudT left(*bent.leftSeam),right(*bent.rightSeam);
  const auto bentDecision=weld::classifySeamWidth(bent,config);
  weld::trimSeamEndpoints(bent,config,nullptr,&bentDecision);
  require(bent.endpointTrimmedEnd>0 && bent.leftSeam->back().x>=170 && bent.leftSeam->back().x<=185,
          "A genuine single centerline tail bend survived or was cut at the wall flare");
  require(bent.leftFrameIndices==bent.rightFrameIndices && bent.leftPointQuality.size()==bent.leftSeam->size(),
          "Single endpoint trim lost paired metadata");
  weld::applySeamWidthRule(bent,config,&bentDecision);
  for (std::size_t i=0;i<bent.centerSeam->size();++i) {
    const auto index=std::size_t(bent.centerFrameIndices[i]-1000);
    require((weld::toVector((*bent.centerSeam)[i])-.5f*(weld::toVector(left[index])+weld::toVector(right[index]))).norm()==0,
            "Single endpoint trimming extrapolated or moved an observed center point");
  }
  // Removing stations must not change the already selected whole-seam mode.
  weld::SeamResult changedVotes;
  for (float width:{8.f,8.f,8.f,20.f,20.f}) {
    const float station=float(changedVotes.leftSeam->size());
    changedVotes.leftSeam->push_back({station,-width/2,0});
    changedVotes.rightSeam->push_back({station,width/2,0});
  }
  const auto fixed=weld::classifySeamWidth(changedVotes,config);
  changedVotes.leftSeam->erase(changedVotes.leftSeam->begin(),changedVotes.leftSeam->begin()+3);
  changedVotes.rightSeam->erase(changedVotes.rightSeam->begin(),changedVotes.rightSeam->begin()+3);
  weld::applySeamWidthRule(changedVotes,config,&fixed);
  require(changedVotes.centerSeam->size()==2 && changedVotes.leftSeam->empty(),
          "Endpoint cuts changed the whole-seam mode after its decision");
}

void testSeamWidthRule() {
  weld::DetectorConfig config;
  config.adaptiveSeamCount = true;
  auto makeWidths=[](const std::vector<float> &widths) {
    weld::SeamResult r; r.frameBased=true;
    for (std::size_t i=0;i<widths.size();++i) {
      r.leftSeam->push_back({float(i),-widths[i]/2,0});
      r.rightSeam->push_back({float(i),widths[i]/2,0});
      r.leftFrameIndices.push_back(100+int(i)); r.rightFrameIndices.push_back(100+int(i));
    }
    return r;
  };
  auto mixed=makeWidths({5,10,12.99f,13,14,9,8});
  weld::applySeamWidthRule(mixed,config);
  require(mixed.ok() && mixed.centerFrameIndices==std::vector<int>({100,101,102,103,104,105,106}) &&
          mixed.leftSeam->empty() && mixed.rightSeam->empty(),"Single majority did not apply to entire seam");
  require(mixed.centerSegmentStarts==std::vector<std::size_t>({0}),"Minority wide frames split the center track");
  const std::string filename="width_rule_test.csv";
  std::string error;
  require(weld::PointCloudIO::saveFrameSeamCsv(filename,mixed,&error),"Whole-seam output CSV failed");
  std::ifstream input(filename); std::string line; int centers=0;
  std::getline(input,line);
  while (std::getline(input,line)) {
    require(line.find("C,")==0,"Single majority CSV retained a side track"); ++centers;
  }
  input.close(); std::remove(filename.c_str());
  require(centers==7,"CSV omitted minority-width frames from center output");
  auto wide=makeWidths({13,14,10,16,12});
  const weld::CloudT originalLeft(*wide.leftSeam);
  const auto originalIds=wide.leftFrameIndices;
  wide.leftSegmentStarts={0,3}; wide.rightSegmentStarts={0,3};
  weld::applySeamWidthRule(wide,config);
  require(wide.centerSeam->empty() && wide.leftFrameIndices==originalIds && wide.rightFrameIndices==originalIds,
          "Double majority did not preserve both tracks for all frames");
  require(wide.leftSegmentStarts==std::vector<std::size_t>({0,3}) && wide.rightSegmentStarts==wide.leftSegmentStarts,
          "Whole-seam mode lost actual trajectory gaps");
  for (std::size_t i=0;i<originalLeft.size();++i)
    require((weld::toVector(originalLeft[i])-weld::toVector((*wide.leftSeam)[i])).norm()==0,
            "Whole-seam double mode changed a narrow frame position");
  auto split=makeWidths({5,6,18,7,8}); split.leftSegmentStarts={0,2}; split.rightSegmentStarts={0,2};
  weld::applySeamWidthRule(split,config);
  require(split.centerSegmentStarts==std::vector<std::size_t>({0,2}),"Single majority bridged an actual gap");
  auto tied=makeWidths({12.99f,13});
  weld::applySeamWidthRule(tied,config);
  require(tied.centerSeam->size()==2 && tied.leftSeam->empty(),"Equal votes did not select single mode");
  auto disabled=makeWidths({5,6,15}); config.adaptiveSeamCount=false;
  weld::applySeamWidthRule(disabled,config);
  require(disabled.leftSeam->size()==3 && disabled.centerSeam->empty(),"Disabled rule changed the output mode");
  config.adaptiveSeamCount=true;
  for (float width:{5.f,10.f,12.f,16.f}) {
    ProfileFixture narrow;
    for (int frame=0;frame<8;++frame) {
      const auto begin=narrow.cloud->size();
      for (int j=-250;j<=250;++j) {
        const float y=j*.1f;
        narrow.cloud->push_back({frame*.5f,y,std::min(15.f,std::max(0.f,2*(std::fabs(y)-width/2)))});
      }
      narrow.frames.push_back({begin,narrow.cloud->size()-begin,frame});
    }
    config.extractMainBody=false; config.voxelSize=0; config.outlierMeanK=0;
    config.useSurfaceRefinement=false;
    const auto result=weld::SeamDetector().detect(narrow.cloud,config,nullptr,nullptr,nullptr,&narrow.frames);
    if (!result.ok()) std::cerr<<"width="<<width<<" "<<result.summary<<'\n';
    require(result.ok(),"Adaptive width detector lost a known groove");
    require(width<13 ? result.centerSeam->size()==8 && result.leftSeam->empty() :
                      result.leftSeam->size()==8 && result.centerSeam->empty(),"End-to-end groove mode wrong");
    for (const auto &p:*result.centerSeam)
      require(std::fabs(p.y)<std::max(.75f,.75f*config.profileBinSize),
              "Narrow center displaced beyond one profile-bin allowance");
  }
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc==2 && std::string(argv[1])=="--shallow-turns-only") {
      testShallowLowerTurns();testUnsupportedShallowTurns();testMeasuredShallowLowerTurns();testShortWallFallback();
      return 0;
    }
    if (argc==2 && std::string(argv[1])=="--endpoint-trim-only") {
      testEndpointTrim();
      testSingleEndpointTrim();
      testEmergingEndpointTrim();
      return 0;
    }
    testIndependentAlgorithmParameters();
    testBottomMidpointGuard();
    testWallTurnCentroid();
    testMeasuredCurvatureStability();
    testMeasuredCurvatureStep();
    testMeasuredCurveLineTail();
    testUniformBottomWidth();
    testUniformWidthFinalFilter();
    testMeasuredCurveLineBranchContinuity();
    testShallowFeatureCutoff();
    testShallowLowerTurns();
    testUnsupportedShallowTurns();
    testMeasuredShallowLowerTurns();
    testSupplementalCandidates();
    testShortWallFallback();
    testFloorSideTurnBoundary();
    testTurnRegions();
    testMultiFrameProjection();
    testProfileCurvature();
    testFinalCurve();
    testCurvatureTrajectory();
    testIndependentWalls();
    testBottomShapeClassification();
    testNarrowWallSpikes();
    testSeamWidthRule();
    testEndpointTrim();
    testSingleEndpointTrim();
    testEmergingEndpointTrim();
    weld::DetectorConfig legacy;
    auto unchanged = fixture();
    const weld::CloudT originalLeft(*unchanged.leftSeam);
    weld::applySeamPostprocess(unchanged, legacy);
    require(!unchanged.leftSmoothing.applied, "Legacy mode unexpectedly changed");
    for (std::size_t i = 0; i < originalLeft.size(); ++i)
      require((weld::toVector(originalLeft[i]) -
               weld::toVector((*unchanged.leftSeam)[i]))
                      .norm() == 0.0f,
              "Legacy coordinates changed");

    auto smoothed = fixture();
    const weld::CloudT noisyLeft(*smoothed.leftSeam);
    weld::DetectorConfig config;
    config.postprocessMode = weld::SeamPostprocessMode::RobustRts;
    config.rtsProcessNoise = 0.015f;
    config.rtsMeasurementStddev = 0.45f;
    config.rtsMaxCorrection = 2.0f;
    weld::applySeamPostprocess(smoothed, config);
    require(smoothed.leftSmoothing.applied &&
                smoothed.rightSmoothing.applied,
            "RTS mode was not applied");
    require(roughness(*smoothed.leftSeam) < roughness(noisyLeft),
            "RTS did not reduce synthetic roughness");
    require(smoothed.leftSmoothing.innovationOutliers > 0,
            "Synthetic outlier was not reported");
    for (std::size_t i = 0; i < smoothed.leftSeam->size(); ++i) {
      require((*smoothed.leftSeam)[i].x == float(i),
              "Longitudinal station changed");
      require((*smoothed.rightSeam)[i].y - (*smoothed.leftSeam)[i].y >=
                  config.minSeamSeparation - 1.0e-4f,
              "Pair separation was violated");
      require(smoothed.leftPointQuality[i].rtsCorrection <=
                  config.rtsMaxCorrection + 1.0e-3f,
              "Correction limit was violated");
    }

    // The stable profile branch must remain usable across sampling-bin
    // changes while its physical detector scale stays fixed.
    const auto profiles = profileFixture();
    for (float bin : {0.8f, 1.2f, 1.6f}) {
      weld::DetectorConfig stable;
      stable.interiorOpeningDirection = -1; // This fixture is deliberately inverted.
      stable.extractMainBody = false;
      stable.profileSamplingMode = weld::ProfileSamplingMode::StableGrid;
      stable.profileBinSize = bin;
      stable.profileReferenceScale = 1.2f;
      const auto detected = weld::SeamDetector().detect(
          profiles.cloud, stable, nullptr, nullptr, nullptr, &profiles.frames);
      require(detected.ok(), "Stable-grid profile detection failed");
      require(detected.leftSeam->size() == detected.rightSeam->size(),
              "Stable-grid pair sizes differ");
    }
    // A dense but smaller disconnected plane must not defeat voxel support.
    weld::CloudT::Ptr hull(new weld::CloudT);
    for (float x = 0; x <= 60; x += .5f) {
      for (float y = -25; y <= 25; y += .5f) hull->push_back({x, y, 0});
      for (int side : {-1, 1}) {
        for (float angle = 0; angle <= 1.570797f; angle += .05f)
          hull->push_back({x, side * (25 + 5 * std::sin(angle)), 5 - 5 * std::cos(angle)});
        for (float z = 5; z <= 30; z += .5f) hull->push_back({x, side * 30.f, z});
      }
    }
    const auto bodyCount = hull->size();
    for (int repeat = 0; repeat < 100; ++repeat)
      for (int x = 0; x <= 15; ++x)
        for (int y = 0; y <= 15; ++y) hull->push_back({float(x), float(y), -20});
    const auto body = weld::extractHullSurfaces(hull, weld::DetectorConfig());
    require(body.ok(), "Hull extraction failed");
    require(body.cloud->size() > .95 * bodyCount, "Lost more than 5% of synthetic hull");
    for (std::size_t i = 0; i < body.cloud->size(); ++i) {
      require(body.sourceIndices[i] < bodyCount, "Retained disconnected interference");
      require((weld::toVector((*body.cloud)[i]) - weld::toVector((*hull)[body.sourceIndices[i]])).norm() == 0,
              "Hull extraction modified raw coordinates");
    }
    weld::DetectorConfig framed;
    framed.interiorOpeningDirection = -1;
    const auto baseline = weld::SeamDetector().detect(profiles.cloud, framed, nullptr, nullptr, nullptr, &profiles.frames);
    require(baseline.ok(), "Protected hull frame detection failed");
    require(!baseline.mainBodyCloud->empty(), "No first-stage body exposed");
    for (int batch : {1, 5, 0}) {
      framed.debugFlag = true;
      framed.debugFrameBatchSize = batch;
      framed.debugStartSlice = batch == 5 ? 7 : 1;
      framed.debugSelectedFrameIndex = batch == 5
          ? profiles.frames[6].frameIndex : -1;
      weld::DebugTrace trace;
      int groups = 0, total = 0, last = -1, firstFrame = -1;
      int firstOrdinal = 0;
      double finalProgress = 0;
      std::string first;
      trace.onStep = [&](const weld::DebugStep &step) {
        if (first.empty()) first = step.title;
        if (step.title.compare(0, 1, "F") == 0) {
          require(step.hasFrameView, "Frame view is not facing the cross-section");
          require(step.overviewCloud && !step.overviewCloud->empty(),
                  "Frame debug view does not include the complete-cloud overview");
          require(step.positionCloud && !step.positionCloud->empty(),
                  "Frame debug view does not highlight the selected slice in the cloud");
          const int count = int(step.metrics.at("batch_count"));
          require(count > 0 && count <= batch, "Invalid batch size");
          require(step.metrics.at("first_frame") > last, "Repeated/out-of-order frame batch");
          if (!groups) firstFrame = int(step.metrics.at("first_frame"));
          last = int(step.metrics.at("last_frame"));
          if (!groups) firstOrdinal = int(step.metrics.at("first_ordinal"));
          require(int(step.metrics.at("total_slices")) >= int(step.metrics.at("last_ordinal")),
                  "Debug progress exceeded the complete point cloud");
          finalProgress = step.metrics.at("progress_percent");
          total += count; ++groups;
        }
        return true;
      };
      const auto observed = weld::SeamDetector().detect(profiles.cloud, framed, nullptr, nullptr, &trace, &profiles.frames);
      require(observed.ok(), "Debug changed detection success");
      require(observed.leftFrameIndices == baseline.leftFrameIndices, "Debug changed original frame IDs");
      require(observed.leftSeam->size() == baseline.leftSeam->size(), "Debug changed output size");
      for (std::size_t i = 0; i < baseline.leftSeam->size(); ++i) {
        require((weld::toVector((*observed.leftSeam)[i]) - weld::toVector((*baseline.leftSeam)[i])).norm() == 0,
                "Debug changed left seam");
        require((weld::toVector((*observed.rightSeam)[i]) - weld::toVector((*baseline.rightSeam)[i])).norm() == 0,
                "Debug changed right seam");
      }
      if (batch > 0) {
        require(first.compare(0, 1, "H") == 0, "First pause is not hull extraction");
        const int expectedStart = std::min(framed.debugStartSlice, int(trace.slices.size()));
        const int expected = int(trace.slices.size()) - expectedStart + 1;
        require(firstOrdinal == expectedStart, "Selected debug start slice was ignored");
        if (framed.debugSelectedFrameIndex >= 0)
          require(firstFrame == framed.debugSelectedFrameIndex,
                  "Selected original frame was not executed by slice debugging");
        require(total == expected, "Last partial batch or selected start position was lost");
        require(groups == (expected + batch - 1) / batch, "Incorrect batch count");
        require(std::fabs(finalProgress - 100.0) < 1e-6,
                "Final debug group did not report complete-cloud progress");
      }
    }
    framed.debugFrameBatchSize = 1;
    framed.debugStartSlice = 1;
    framed.debugSelectedFrameIndex = -1;
    weld::DebugTrace skip;
    int calls = 0;
    skip.onStep = [&](const weld::DebugStep &) { ++calls; return false; };
    require(weld::SeamDetector().detect(profiles.cloud, framed, nullptr, nullptr, &skip, &profiles.frames).ok(),
            "Skipping debug stopped detection");
    require(calls == 1, "Skip did not suppress remaining callbacks");
    // Every frame that reaches the debug viewer without a candidate pair must
    // retain a concrete explanation, including when that frame is replayed.
    ProfileFixture flat;
    for (int frame = 0; frame < 3; ++frame) {
      const auto begin = flat.cloud->size();
      for (int j = 0; j < 40; ++j)
        flat.cloud->push_back({float(frame), .25f * j, 0.f});
      flat.frames.push_back({begin, flat.cloud->size() - begin, frame});
    }
    auto noCandidate = framed;
    noCandidate.extractMainBody = false;
    noCandidate.voxelSize = 0;
    noCandidate.outlierMeanK = 0;
    noCandidate.profileMergeFrameCount = 1;
    noCandidate.debugFrameBatchSize = 1;
    weld::DebugTrace noCandidateTrace;
    int reasonMessages = 0;
    noCandidateTrace.onStep = [&](const weld::DebugStep &step) {
      if (step.title.compare(0, 1, "F") == 0 &&
          step.message.find("【无候选原因】") != std::string::npos)
        ++reasonMessages;
      return true;
    };
    weld::SeamDetector().detect(flat.cloud, noCandidate, nullptr, nullptr,
                                &noCandidateTrace, &flat.frames);
    require(noCandidateTrace.slices.size() == flat.frames.size(),
            "No-candidate debug fixture lost frames");
    for (const auto &slice : noCandidateTrace.slices) {
      require(slice.pairs.empty(), "Flat debug fixture unexpectedly found a pair");
      require(!slice.noCandidateReason.empty(),
              "No-candidate frame did not retain its rejection reason");
    }
    require(reasonMessages == int(flat.frames.size()),
            "Live frame debug message did not display each no-candidate reason");
    // A seam can enter the scan as a shallow notch or bead before it develops
    // enough vertical wall height for the normal wall-constrained detector.
    ProfileFixture emerging;
    for (int frame = 0; frame < 42; ++frame) {
      const auto begin = emerging.cloud->size();
      for (int j = -140; j <= 140; ++j) {
        const float y = .25f * j;
        float z = .01f * y;
        if (frame < 12) {
          const float amplitude = .5f + .08f * frame;
          // A real shallow flat-bottom groove has two distinct lower feet.
          // Its width must satisfy the configured pair separation; a single
          // V-bottom cannot be turned into two shoulder points to pass it.
          z -= amplitude * std::max(0.f, std::min(1.f,(8.f-std::fabs(y))/3.f));
        } else if (frame < 24) {
          const float amplitude = .6f + .06f * (frame - 12);
          z += amplitude * std::max(0.f, 1.f - std::fabs(y) / 5.f);
        } else {
          const float depth = .45f * (frame - 23);
          z += std::min(depth, std::max(0.f, (std::fabs(y) - 5.f) * 1.5f));
        }
        emerging.cloud->push_back({.5f * frame, y, z});
      }
      emerging.frames.push_back(
          {begin, emerging.cloud->size() - begin, frame});
    }
    auto emergingConfig = noCandidate;
    // Explicit compatibility test for shallow recovery; production defaults
    // now stop below minimumFeatureDepth, covered by testShallowFeatureCutoff.
    emergingConfig.minimumFeatureDepth = 0;
    // This sequence contains both a notch and a bead; do not inherit the
    // inverted opening of the unrelated framed fixture above.
    emergingConfig.interiorOpeningDirection = 0;
    emergingConfig.profileLocalizationMethod =
        weld::ProfileLocalizationMethod::CurvatureCandidates;
    emergingConfig.profileMergeFrameCount = 5;
    emergingConfig.trimUnstableEndpoints = true;
    emergingConfig.debugFrameBatchSize = 0;
    weld::DebugTrace emergingTrace;
    const auto emergingResult = weld::SeamDetector().detect(
        emerging.cloud, emergingConfig, nullptr, nullptr, &emergingTrace,
        &emerging.frames);
    require(emergingResult.ok(), "Shallow seam onset was not detected");
    if (emergingResult.leftFrameIndices.empty() || emergingResult.leftFrameIndices.front() >= 12) {
      std::cerr << "emerging first=" << (emergingResult.leftFrameIndices.empty() ? -1 : emergingResult.leftFrameIndices.front()) << '\n';
      for (const auto &slice:emergingTrace.slices) {
        std::cerr << "frame=" << slice.frameIndex << " pairs=" << slice.pairs.size();
        for (const auto &p:slice.pairs) std::cerr << " [" << p.points[0].y << ',' << p.bottomMidpointU << ',' << p.points[1].y << "]";
        for (const auto &g:slice.gates) if (g.first=="candidate_crosses_bottom_midpoint" || g.first=="emerging_profile_fallback") std::cerr << ' ' << g.first << '=' << g.second;
        std::cerr << '\n';
      }
    }
    require(!emergingResult.leftFrameIndices.empty() &&
                emergingResult.leftFrameIndices.front() < 12,
            "Shallow onset was discarded until stable walls appeared");
    bool usedEmergingProfile = false;
    for (const auto &slice : emergingTrace.slices) {
      usedEmergingProfile |=
          slice.gates.find("emerging_profile_fallback") != slice.gates.end();
      for (const auto &pair : slice.pairs) {
        require(pair.hasBottomMidpoint && std::fabs(pair.bottomMidpointU) < 1.f,
                "Tilted shallow feature moved the bottom reference toward a scan edge");
        require(pair.points[0].y < pair.bottomMidpointU &&
                pair.points[1].y > pair.bottomMidpointU,
                "Emerging profile bypassed the bottom midpoint guard");
      }
    }
    require(usedEmergingProfile,
            "Shallow onset did not use the emerging-profile fallback");
    // A known 11-unit trough is narrower than the former 16-unit gate.
    // Assert localization against its analytic lower corners, not point count.
    ProfileFixture narrow;
    for (int frame = 0; frame < 40; ++frame) {
      const auto begin = narrow.cloud->size();
      for (int j = -350; j <= 350; ++j) {
        const float y = j * .1f;
        const float z = std::min(16.f, std::max(0.f, (std::abs(y) - 5.5f) * 2));
        narrow.cloud->push_back({2.f * frame, y, z + .03f * frame});
      }
      narrow.frames.push_back({begin, narrow.cloud->size() - begin, frame});
    }
    weld::DetectorConfig trough;
    trough.adaptiveSeamCount = false;
    // Explicit opt-in: the product default now skips S12-S16 refinement.
    trough.useSurfaceRefinement = true;
    // This block exercises the legacy candidate trust-region guard itself.
    trough.finalSmoothingEnabled = false;
    trough.extractMainBody = false;
    trough.voxelSize = 0;
    trough.outlierMeanK = 0;
    weld::DebugTrace troughTrace;
    const auto narrowResult = weld::SeamDetector().detect(narrow.cloud, trough, nullptr, nullptr, &troughTrace, &narrow.frames);
    require(narrowResult.ok(), "Narrow positive-V trough was rejected");
    for (std::size_t i = 0; i < narrowResult.leftSeam->size(); ++i) {
      const float z = .03f * narrowResult.leftFrameIndices[i];
      require(std::hypot((*narrowResult.leftSeam)[i].y + 5.5f, (*narrowResult.leftSeam)[i].z - z) < 1.5f,
              "Left result missed the known lower corner");
      require(std::hypot((*narrowResult.rightSeam)[i].y - 5.5f, (*narrowResult.rightSeam)[i].z - z) < 1.5f,
              "Right result missed the known lower corner");
      require((*narrowResult.rightSeam)[i].y - (*narrowResult.leftSeam)[i].y >= trough.minSeamSeparation,
              "Correction guard broke pair separation");
      require((weld::toVector((*narrowResult.leftSeam)[i]) - weld::toVector(troughTrace.interpolatedPath[0][i])).norm() <= trough.maxCandidateCorrection + 1e-3f,
              "Final correction escaped the candidate trust region");
    }
    require(!troughTrace.wallPoints.empty(), "Enabled surface refinement did not run");
    // Switching off geometry refinement must preserve frame/candidate identity,
    // bypass both fitting stages and still localize the analytic trough.
    for (bool surfaceEnabled : {true, false}) {
      auto optional = trough;
      optional.useSurfaceRefinement = surfaceEnabled;
      optional.useLocalFloorConstraint = false;
      weld::DebugTrace optionalTrace;
      const auto optionalResult = weld::SeamDetector().detect(
          narrow.cloud, optional, nullptr, nullptr, &optionalTrace, &narrow.frames);
      require(optionalResult.ok(), "Optional refinement lost the narrow trough");
      require(optionalResult.leftFrameIndices == narrowResult.leftFrameIndices &&
              optionalResult.rightFrameIndices == narrowResult.rightFrameIndices,
              "Refinement switch changed output frame identity");
      require(surfaceEnabled || optionalTrace.wallPoints.empty(),
              "Disabled refinement still fitted surfaces");
      for (int side = 0; side < 2; ++side) {
        const auto &points = side ? *optionalResult.rightSeam : *optionalResult.leftSeam;
        const auto &quality = side ? optionalResult.rightPointQuality : optionalResult.leftPointQuality;
        for (std::size_t i = 0; i < points.size(); ++i) {
          require(!quality[i].floorAccepted, "Disabled floor constraint still applied");
          require(surfaceEnabled || (!quality[i].wallAccepted && quality[i].geometryCorrection == 0),
                  "Disabled refinement still moved points onto a wall");
          const float expectedY = side ? 5.5f : -5.5f;
          const float expectedZ = .03f * optionalResult.leftFrameIndices[i];
          require(std::hypot(points[i].y - expectedY, points[i].z - expectedZ) < 1.5f,
                  "Optional refinement missed the analytic corner");
          require((weld::toVector(optionalTrace.interpolatedPath[side][i]) -
                   weld::toVector(troughTrace.interpolatedPath[side][i])).norm() == 0,
                  "Refinement switch changed upstream candidates");
        }
      }
    }
    // The master switch also overrides a retained, checked floor option.
    auto masterOff = trough;
    masterOff.useSurfaceRefinement = false;
    weld::DebugTrace masterOffTrace;
    const auto masterOffResult = weld::SeamDetector().detect(
        narrow.cloud, masterOff, nullptr, nullptr, &masterOffTrace, &narrow.frames);
    require(masterOffResult.ok() && masterOffTrace.wallPoints.empty(),
            "Floor option bypassed the disabled master switch");
    for (const auto &quality : masterOffResult.leftPointQuality)
      require(!quality.wallAccepted && !quality.floorAccepted,
              "Master switch left an active geometry constraint");
    trough.interiorOpeningDirection = -1;
    require(!weld::SeamDetector().detect(narrow.cloud, trough, nullptr, nullptr, nullptr, &narrow.frames).ok(),
            "Wrong opening direction accepted the opposite trough");
    auto clippedInput = weld::CloudT::Ptr(new weld::CloudT(*profiles.cloud));
    auto clippedFrames = profiles.frames;
    clippedFrames.push_back({clippedInput->size(), 1, 100});
    clippedInput->push_back({70, 80, -30});
    framed.debugFlag = false;
    const auto clipped = weld::SeamDetector().detect(clippedInput, framed, nullptr, nullptr, nullptr, &clippedFrames);
    require(clipped.ok() && clipped.inputFrameCount == 25, "Completely clipped frame broke extraction");
    require(clipped.leftFrameIndices == baseline.leftFrameIndices, "Clipped frame renumbered original IDs");
    clippedFrames.back().begin += 1;
    require(!weld::SeamDetector().detect(clippedInput, framed, nullptr, nullptr, nullptr, &clippedFrames).ok(),
            "Malformed frame map was accepted before body masking");
    const char *textPath = "target_header_test.txt";
    { std::ofstream file(textPath); file << "index x y z target\n1 1 2 3 7\n2 4 5 6 9\n"; }
    std::string error;
    std::vector<weld::ScanFrameRange> loadedFrames;
    const auto loaded = weld::PointCloudIO::load(textPath, &error, &loadedFrames);
    std::remove(textPath);
    require(loaded->size() == 2 && loadedFrames.size() == 2 && loadedFrames.back().frameIndex == 9,
            "target header lost original frames");
    std::cout << "RTS, profiles, hull body and frame stepping tests passed.\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
