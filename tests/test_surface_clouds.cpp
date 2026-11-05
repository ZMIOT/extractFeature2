#include "SurfaceClouds.h"
#include "SeamDetector.h"
#include "AlgorithmParameters.h"
#include "SeamSmoother.h"
#include "PointCloudIO.h"
#include <cmath>
#include <iostream>
#include <set>
#include <map>
#include <tuple>
#include <stdexcept>
#include <fstream>
#include <iomanip>

namespace {
void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}
using Key = std::tuple<float,float,float>;
Key key(const weld::PointT &p) { return {p.x,p.y,p.z}; }
void checkMeasured(const weld::CloudT &input, const weld::SeamResult &result) {
  std::set<Key> originals, used;
  for (const auto &p : input) originals.insert(key(p));
  for (const auto &cloud : {result.leftSurface,result.rightSurface,result.weldSurface}) {
    require(!cloud->empty(), "A surface is empty");
    std::set<Key> current;
    for (const auto &p : *cloud) {
      require(originals.count(key(p))!=0, "Surface coordinate was generated or changed");
      require(!used.count(key(p)), "Surface clouds overlap");
      current.insert(key(p));
    }
    used.insert(current.begin(),current.end());
  }
}

void geometryTest(weld::Axis axis, int opening, int wallShape, bool framed, bool testEnds=false) {
  weld::DetectorConfig config;
  config.longitudinalAxis=axis;
  config.interiorOpeningDirection=opening;
  config.sliceStep=2;
  const int a=int(axis),u=a==0?1:0,v=a==2?1:2;
  auto wallWidth=[&](float height) {
    return 8+(wallShape==1 ? 0 : wallShape==2 ? 3*(1-std::exp(-height/3)) : .5f*height);
  };
  auto point=[&](float station,float width,float height) {
    weld::PointT p;
    p.data[a]=framed ? -station+.15f*width : station;
    p.data[u]=width;p.data[v]=opening*height;
    return p;
  };
  weld::CloudT::Ptr input(new weld::CloudT);
  std::vector<weld::ScanFrameRange> frames;
  std::vector<weld::SurfaceSection> sections;
  for(int f=0;f<12;++f) {
    const float station=2.f*f;
    const auto begin=input->size();
    // A bowed weld face proves this is measured surface extraction, not a
    // fitted bottom plane. The extra sheet must not join the weld output.
    for(int j=-80;j<=80;++j) {
      const float y=.1f*j;
      input->push_back(point(station,y,.02f*(64-y*y)));
      if(j%4==0) input->push_back(point(station,y,-12));
    }
    for(int side : {-1,1}) {
      for(int j=1;j<=150;++j) {
        const float z=.1f*j,y=side*wallWidth(z);
        input->push_back(point(station,y,z));
        input->push_back(point(station,y+side*20,z)); // disconnected parallel wall
      }
      for(int j=1;j<=80;++j)
        input->push_back(point(station,side*(wallWidth(15)+.1f*j),15));
    }
    frames.push_back({begin,input->size()-begin,100+f});
    if(f==5 || f==6) continue; // No support may leak across absent sections.
    weld::SurfaceSection section;
    section.station=station;section.frameIndex=100+f;
    section.left=point(station,-8,0);section.right=point(station,8,0);
    const float slope=wallShape==1?0:wallShape==2?1:.5f;
    section.leftNormal=Eigen::Vector2f(1,slope*opening).normalized();
    section.rightNormal=Eigen::Vector2f(1,-slope*opening).normalized();
    sections.push_back(section);
  }
  weld::SeamResult result;
  weld::extractSurfaceClouds(input,sections,config,result,framed?&frames:nullptr);
  std::cout<<"geometry axis="<<a<<" opening="<<opening<<" shape="<<wallShape<<" framed="<<framed
           <<": "<<result.leftSurface->size()<<"/"<<result.rightSurface->size()<<"/"<<result.weldSurface->size()<<'\n';
  checkMeasured(*input,result);
  if(testEnds) {
    std::array<weld::SurfaceEnd,2> ends;
    ends[0].valid=ends[1].valid=true;
    ends[0].station=8;ends[0].point=point(8,-8,0);
    ends[1].station=18;ends[1].point=point(18,8,0);
    weld::SeamResult limited;
    weld::extractSurfaceClouds(input,sections,config,limited,framed?&frames:nullptr,true,&ends);
    checkMeasured(*input,limited);
    for(int side=0;side<3;++side) {
      const auto before=side==0?result.leftSurface:side==1?result.rightSurface:result.weldSurface;
      const auto after=side==0?limited.leftSurface:side==1?limited.rightSurface:limited.weldSurface;
      std::set<Key> expected;
      for(const auto &p:*before) {
        const float station=framed?-(p.data[a]-.15f*p.data[u]):p.data[a];
        const float stop=side==1?18:8;
        if(side==2 || station<=stop+1.e-4f)expected.insert(key(p));
      }
      require(after->size()==expected.size(),"Independent surface endpoints lost valid coverage");
      for(const auto &p:*after)require(expected.count(key(p))!=0,"Surface exceeded its independent endpoint");
      if(side<2)require(after->size()<before->size(),"Endpoint limit did not remove the wall tail");
      else require(after->size()==before->size(),"Wall endpoint limit removed measured floor");
    }
    require(limited.weldSurfaceModel.valid(),"Independent floor model was not fitted");
    require(limited.weldSurfaceModel.patches.back().lastFrame==result.weldSurfaceModel.patches.back().lastFrame,
            "Wall endpoint limit shortened floor model coverage");
    config.limitSurfaceEndsToTracks=false;
    weld::SeamResult restored;
    weld::extractSurfaceClouds(input,sections,config,restored,framed?&frames:nullptr,false,&ends);
    require(restored.weldSurface->size()==result.weldSurface->size(),"Surface end opt-out did not restore coverage");
  }
  require(result.leftSurface->size()>1200 && result.rightSurface->size()>1200,
          "Too much side-wall coverage lost");
  require(result.weldSurface->size()>1300, "Too much bowed weld-face coverage lost");
  for(int side : {0,1,2}) {
    const auto cloud=side==0?result.leftSurface:side==1?result.rightSurface:result.weldSurface;
    for(const auto &p : *cloud) {
      const float y=p.data[u], z=opening*p.data[v];
      const float station=framed?-(p.data[a]-.15f*y):p.data[a];
      require(std::abs(station-10)>1.e-3f && std::abs(station-12)>1.e-3f,
              "Filled a missing frame/section");
      if(side==2) {
        // For rounded feet, the transition inside the two fitted big faces
        // belongs to the weld face; upper walls and shoulders do not.
        const float limit=wallShape==2?wallWidth(4):8.001f;
        require(std::abs(y)<=limit && z>=-.001f && z<(wallShape==2?4:1.3f),
                "Weld face includes a wall, shoulder or interference sheet");
      } else {
        require((side==0?y<0:y>0) && z>0 && z<=15.001f,
                "Side identity or opening direction incorrect");
        require(std::abs(std::abs(y)-wallWidth(z))<1.e-3f,
                "Outer shoulder included in a side face");
      }
    }
  }
}

void detectorTest() {
  std::string error;
  std::vector<weld::ScanFrameRange> frames;
  const auto input=weld::PointCloudIO::load(
      std::string(SEAM_TEST_DATA_DIR)+"/curvature_2_frames_2040_2110.txt",&error,&frames);
  require(input && !input->empty(), "Measured fixture did not load");
  for(bool framed : {true}) {
    auto config=weld::AlgorithmParameters().get(weld::ProfileLocalizationMethod::CurvatureCandidates);
    // The fixture is already body-masked and cropped from a complete scan.
    config.extractMainBody=false;config.voxelSize=0;config.outlierMeanK=0;
    config.trimUnstableEndpoints=false;
    const auto result=weld::SeamDetector().detect(input,config,nullptr,nullptr,nullptr,framed?&frames:nullptr);
    if(!result.ok()) std::cerr<<result.summary<<'\n';
    require(result.ok(), "Detector failed measured fixture");
    checkMeasured(*input,result);
    require(result.leftSurface->size()>100 && result.rightSurface->size()>100 && result.weldSurface->size()>100,
            "Measured surface support unexpectedly small");
    std::cout<<(framed?"framed":"spatial")<<" measured: "
      <<result.leftSurface->size()<<" / "<<result.rightSurface->size()<<" / "
      <<result.weldSurface->size()<<" surface points; tracks "<<result.leftSeam->size()<<'\n';
  }
  weld::CloudT::Ptr synthetic(new weld::CloudT);
  for(int x=0;x<80;++x) for(int j=-250;j<=250;++j) {
    const float y=.1f*j;
    synthetic->push_back(weld::PointT(.8f*x,y,std::min(15.f,std::max(0.f,2*(std::abs(y)-8)))));
  }
  const auto config=weld::AlgorithmParameters().get(weld::ProfileLocalizationMethod::CurvatureCandidates);
  const auto result=weld::SeamDetector().detect(synthetic,config);
  if(!result.ok()) std::cerr<<result.summary<<'\n';
  require(result.ok(), "Unmodified defaults failed complete synthetic groove");
  checkMeasured(*synthetic,result);
  std::cout<<"unmodified defaults: "<<result.leftSurface->size()<<" / "<<result.rightSurface->size()
           <<" / "<<result.weldSurface->size()<<" surface points\n";
}

void roughWeldTest() {
  weld::CloudT::Ptr input(new weld::CloudT);
  for(int j=-150;j<=150;++j) for(int k=0;k<4;++k) {
    const float y=.05f*j;
    input->push_back(weld::PointT(0,y,.35f*std::sin(2.31f*j+1.37f*k)));
  }
  for(int side:{-1,1})for(int j=0;j<=300;++j)
    input->push_back(weld::PointT(0,side*8.f,.05f*j));
  weld::SurfaceSection section;
  section.left=weld::PointT(0,-8,0);section.right=weld::PointT(0,8,0);
  section.leftNormal=section.rightNormal=Eigen::Vector2f(1,0);
  weld::SeamResult result;
  weld::extractSurfaceClouds(input,{section},weld::DetectorConfig(),result);
  checkMeasured(*input,result);
  require(result.weldSurface->size()>1100,"Rough weld face was rejected for lacking planar normals");
  for(const auto &p:*result.weldSurface)require(std::abs(p.z)<=.351f,"Side leaked into rough weld face");
}

void shallowTipTest() {
  std::string error;
  std::vector<weld::ScanFrameRange> frames;
  const auto input=weld::PointCloudIO::load(
      std::string(SEAM_TEST_DATA_DIR)+"/shallow_3_frames_375_445.txt",&error,&frames);
  require(input && !input->empty(), "Shallow measured fixture did not load");
  auto config=weld::AlgorithmParameters().get(weld::ProfileLocalizationMethod::CurvatureCandidates);
  config.extractMainBody=false;config.voxelSize=0;config.outlierMeanK=0;
  for(float bin:{1.2f,.8f}) {
  config.profileBinSize=bin;
  const auto result=weld::SeamDetector().detect(input,config,nullptr,nullptr,nullptr,&frames);
  require(result.ok(), "Shallow measured detection failed");
  checkMeasured(*input,result);
  std::map<Key,int> sourceFrames;
  for(const auto &frame:frames)for(std::size_t j=frame.begin;j<frame.begin+frame.count;++j)
    sourceFrames.emplace(key((*input)[j]),frame.frameIndex);
  std::map<int,std::array<int,3>> coverage;
  for(int side=0;side<3;++side) {
    const auto cloud=side==0?result.leftSurface:side==1?result.rightSurface:result.weldSurface;
    for(const auto &p:*cloud)++coverage[sourceFrames.at(key(p))][side];
  }
  for(std::size_t i=0;i<result.leftFrameIndices.size();++i) {
    const int frame=result.leftFrameIndices[i];
    const auto counts=coverage[frame];
    if(counts[0]<5 || counts[1]<5 || counts[2]<50) {
      std::cerr<<"Shallow tip frame "<<frame<<" coverage "<<counts[0]<<'/'<<counts[1]<<'/'<<counts[2]
               <<" normals "<<result.leftPointQuality[i].trackedNormal.transpose()<<" / "
               <<result.rightPointQuality[i].trackedNormal.transpose()<<'\n';
      require(false,"Measured tip surfaces do not cover the detected trajectory");
    }
  }
  std::cout<<"Shallow tip: all three surfaces cover "<<result.leftFrameIndices.size()<<" trajectory frames\n";
  }
}

void shortTransitionTest(weld::Axis axis,int opening) {
  weld::DetectorConfig config;
  config.longitudinalAxis=axis;config.interiorOpeningDirection=opening;
  const int a=int(axis),u=a==0?1:0,v=a==2?1:2;
  auto point=[&](float station,float width,float height) {
    weld::PointT p;p.data[a]=station;p.data[u]=width;p.data[v]=opening*height;return p;
  };
  weld::CloudT::Ptr input(new weld::CloudT);
  std::vector<weld::SurfaceSection> sections;
  for(int f=0;f<5;++f) {
    weld::SurfaceSection section;section.station=2.f*f;
    section.left=point(section.station,-8,0);section.right=point(section.station,8,0);
    // End sections have a shallow measured transition but no reliable tracked
    // plane. The absent middle section must still remain absent in all outputs.
    if(f!=2)sections.push_back(section);
    for(int j=-200;j<=200;++j) {
      const float width=.05f*j, height=.25f*std::max(0.f,std::abs(width)-8);
      input->push_back(point(section.station,width,height));
      input->push_back(point(section.station,width,height+8)); // detached sheet
    }
  }
  weld::SeamResult result;
  weld::extractSurfaceClouds(input,sections,config,result);
  checkMeasured(*input,result);
  for(int side=0;side<3;++side) {
    std::map<float,int> counts;
    const auto cloud=side==0?result.leftSurface:side==1?result.rightSurface:result.weldSurface;
    for(const auto &p:*cloud) {
      ++counts[p.data[a]];
      require(opening*p.data[v]<=.501f,"Detached sheet entered short-wall recovery");
      require(side==2?std::abs(p.data[u])<=8.001f:
              (side==0?p.data[u]<=-8:p.data[u]>=8),"Short transition face identity incorrect");
    }
    require(counts.size()==4 && !counts.count(4),"Short-wall recovery crossed a missing section");
    for(const auto &section:sections)require(counts[section.station]>20,"Short surfaces missed a trajectory section");
  }
  // A flat sheet through both supplied toes provides no side-wall evidence.
  for(auto &p:*input)p.data[v]=0;
  weld::extractSurfaceClouds(input,sections,config,result);
  require(result.leftSurface->empty() && result.rightSurface->empty() && result.weldSurface->empty(),
          "Flat floor was invented into three surfaces");
}

void sparseSideKeepsFloorTest(weld::Axis axis,int opening,int sparseSide) {
  weld::DetectorConfig config;
  config.longitudinalAxis=axis;config.interiorOpeningDirection=opening;
  const int a=int(axis),u=a==0?1:0,v=a==2?1:2;
  auto point=[&](float station,float width,float height) {
    weld::PointT p;p.data[a]=station;p.data[u]=width;p.data[v]=opening*height;return p;
  };
  weld::CloudT::Ptr input(new weld::CloudT);
  std::vector<weld::SurfaceSection> sections;
  for(int f=0;f<8;++f) {
    const float station=2.f*f;
    weld::SurfaceSection section;section.station=station;
    section.left=point(station,-8,0);section.right=point(station,8,0);
    section.leftNormal=section.rightNormal=Eigen::Vector2f(1,0);
    sections.push_back(section);
    for(int j=-80;j<=80;++j)input->push_back(point(station,.1f*j,0));
    for(int side=0;side<2;++side) {
      // The first three sections model the shallow end in the screenshot: a
      // tracked foot exists, but one side has too few measurements to fit.
      // Later sections make both walls valid, exercising the common-span trim.
      const int samples=side==sparseSide && f<3?3:150;
      for(int j=1;j<=samples;++j)
        input->push_back(point(station,side?8:-8,.1f*j));
    }
  }
  weld::SeamResult result;
  weld::extractSurfaceClouds(input,sections,config,result);
  const auto sparse=sparseSide?result.rightSurface:result.leftSurface;
  const auto supported=sparseSide?result.leftSurface:result.rightSurface;
  for(const auto &p:*sparse)
    require(p.data[a]>=6,"Sparse endpoint side was fabricated into a fitted wall");
  require(supported->size()>800,"Supported side wall disappeared");
  require(result.weldSurface->size()>1200,"Sparse side erased the measured bottom");
  int sparseEndFloor=0;
  for(const auto &p:*result.weldSurface)
    if(p.data[a]<6)++sparseEndFloor;
  require(sparseEndFloor>400,"Common wall span trimmed the sparse-end bottom");
  for(const auto &p:*result.weldSurface)
    require(std::abs(p.data[u])<8.001f && std::abs(p.data[v])<1.e-4f,
            "Sparse-side bottom escaped its measured feet or included a wall");
}

void singleWallFloorTest(weld::Axis axis,int opening,int side) {
  weld::DetectorConfig config;
  config.longitudinalAxis=axis;config.interiorOpeningDirection=opening;
  const int a=int(axis),u=a==0?1:0,v=a==2?1:2,inward=side==0?1:-1;
  auto point=[&](float station,float distance,float height) {
    weld::PointT p;p.data[a]=station+.1f*distance;
    p.data[u]=inward*distance;p.data[v]=opening*height;return p;
  };
  weld::CloudT::Ptr input(new weld::CloudT);
  std::vector<weld::ScanFrameRange> frames;
  std::vector<weld::SurfaceSection> sections;
  for(int f=0;f<12;++f) {
    const auto begin=input->size();
    for(int j=0;j<=300;++j) {
      const float distance=.1f*j;
      // A curved measured floor, a detached parallel sheet, and a disconnected
      // continuation beyond a real gap. None requires an opposite trajectory.
      const float height=.001f*distance*distance+(f>=7?2*std::max(0.f,distance-20):0);
      input->push_back(point(float(f),distance,height));
      input->push_back(point(float(f),distance,8+.001f*distance*distance));
      input->push_back(point(float(f),distance+40,0));
    }
    // At the first section the accepted tracked foot remains, but the wall
    // itself has no fit support. The measured bottom must still be retained.
    for(int j=1;f!=0 && j<=150;++j) {
      const float height=.1f*j;
      input->push_back(point(float(f),-.5f*height,height));
      input->push_back(point(float(f),-8-.1f*j,15)); // outer shoulder
    }
    frames.push_back({begin,input->size()-begin,100+f});
    if(f==6)continue;
    weld::SurfaceSection section;section.station=float(f);section.frameIndex=100+f;
    section.leftValid=side==0;section.rightValid=side==1;
    section.left=point(float(f),side==0?0:8,0);
    section.right=point(float(f),side==0?8:0,0);
    (side==0?section.leftNormal:section.rightNormal)=Eigen::Vector2f(1,inward*opening*.5f).normalized();
    sections.push_back(section);
  }
  weld::SeamResult result;
  weld::extractSurfaceClouds(input,sections,config,result,&frames);
  require((side==0?result.rightSurface:result.leftSurface)->empty(),"Single wall invented the opposite surface");
  require(!(side==0?result.leftSurface:result.rightSurface)->empty(),"Single wall surface disappeared");
  require(result.weldSurface->size()>2600,"Single wall lost its connected measured floor");
  require(result.weldSurfaceModel.valid(),"Single wall floor model was not fitted");
  require(result.weldSurfaceModel.patches.front().firstFrame==100,
          "Sparse tracked wall foot lost the first measured floor section");
  std::set<Key> measured,used;
  for(const auto &p:*input)measured.insert(key(p));
  for(const auto &cloud:{result.leftSurface,result.rightSurface,result.weldSurface})
    for(const auto &p:*cloud) {
      require(measured.count(key(p)) && used.insert(key(p)).second,"Single wall surfaces changed or reused a measurement");
      require(std::abs(p.data[a]-.1f*inward*p.data[u]-6)>1.e-4f,"Single wall filled an absent section");
    }
  for(const auto &p:*result.weldSurface) {
    const float distance=inward*p.data[u],height=opening*p.data[v];
    require(distance>0 && distance<=30.001f && height<1,"Single wall floor crossed a gap, wall or detached sheet");
  }
  std::array<weld::SurfaceEnd,2> ends;
  ends[side].valid=true;ends[side].station=4;
  weld::SeamResult limited;
  weld::extractSurfaceClouds(input,sections,config,limited,&frames,true,&ends);
  require(limited.weldSurface->size()==result.weldSurface->size(),"Single wall endpoint clipped measured floor");
  for(std::size_t i=0;i<result.weldSurface->size();++i)
    require(key((*limited.weldSurface)[i])==key((*result.weldSurface)[i]),"Single wall endpoint changed floor support");
  for(const auto &p:*(side==0?limited.leftSurface:limited.rightSurface))
    require(p.data[a]-.1f*inward*p.data[u]<=4+1.e-4f,"Single wall ignored its trajectory endpoint");
  require(limited.weldSurfaceModel.valid() && limited.weldSurfaceModel.patches.back().lastFrame==111,
          "Single wall endpoint clipped floor model");
  for(auto &section:sections)section.leftValid=section.rightValid=false;
  weld::extractSurfaceClouds(input,sections,config,result,&frames,false);
  require(result.weldSurface->empty(),"Unsupported section generated a floor");
}

void boundedTwoWallFloorTest(weld::Axis axis,bool framed) {
  weld::DetectorConfig config;config.longitudinalAxis=axis;config.sliceStep=2;
  const int a=int(axis),u=a==0?1:0,v=a==2?1:2;
  auto point=[&](float station,float width,float height) {
    weld::PointT p;p.data[a]=framed?-station+.15f*width:station;
    p.data[u]=width;p.data[v]=height;return p;
  };
  weld::CloudT::Ptr input(new weld::CloudT);
  std::vector<weld::SurfaceSection> sections;
  std::vector<weld::ScanFrameRange> frames;
  for(int f=0;f<20;++f) {
    const float station=2.f*f;const auto begin=input->size();
    weld::SurfaceSection s;s.station=station;s.frameIndex=framed?100+f:-1;
    s.left=point(station,-8,0);s.right=point(station,8,0);
    s.leftValid=f>=3 && f<=17;s.rightValid=f<=14;
    s.leftNormal=Eigen::Vector2f(1,.5f).normalized();
    s.rightNormal=Eigen::Vector2f(1,-.5f).normalized();
    // A transient lost left wall must not turn an exterior floor sheet into
    // weld material just because the right wall remains visible (6.txt).
    if(f==7)s.leftValid=false;
    for(int j=(f==7?-220:-80);j<=80;++j)input->push_back(point(station,.1f*j,0));
    for(int side:{0,1})if(side?s.rightValid:s.leftValid)
      for(int j=1;j<=150;++j)input->push_back(point(station,(side?1:-1)*(8+.05f*j),.1f*j));
    sections.push_back(s);frames.push_back({begin,input->size()-begin,100+f});
  }
  weld::SeamResult result;
  weld::extractSurfaceClouds(input,sections,config,result,framed?&frames:nullptr);
  checkMeasured(*input,result);
  require(result.weldSurface->size()>1500,"Two-wall bound discarded supported interior floor");
  for(const auto &p:*result.weldSurface) {
    const float station=framed?-p.data[a]+.15f*p.data[u]:p.data[a];
    require(station>=-1.e-4f && station<=34+1.e-4f,"Floor protrudes past all measured wall support");
    require(p.data[u]>=-8 && p.data[u]<=8,"Floor crossed the measured floor edge");
    if(station>=6-1.e-4f && station<=28+1.e-4f)
      require(p.data[u]>-8 && p.data[u]<8,"Interior floor crossed a measured wall foot");
  }
  require(result.weldSurfaceModel.valid(),"Bounded floor model is missing");
  require(result.weldSurfaceModel.patches.front().firstStation==0 &&
          result.weldSurfaceModel.patches.back().lastStation==34,
          "One missing endpoint trajectory clipped the connected floor model");
  for(const auto &patch:result.weldSurfaceModel.patches) {
    require(patch.firstStation>=0 && patch.lastStation<=34,"Floor model extrapolated its measured support");
    for(int i=0;i<=40;++i)for(int j=0;j<=20;++j) {
      Eigen::Vector3d p,n;require(patch.evaluate(i/40.,j/20.,p,&n),"Bounded floor has a singular fit");
      const double station=framed?-p[a]+.15*p[u]:p[a];
      // The open measured edge and plane-clipped interior have slightly
      // different widths. A spline can overshoot by its fitting tolerance;
      // exact corridor membership is checked on the measured points above.
      const float fitTolerance=2*std::max(config.hullSurfaceDistance,config.profileWallTolerance);
      require(station>=-1.e-3 && station<=34+1.e-3 && p[u]>=-8-fitTolerance && p[u]<=8+fitTolerance,
              "Fitted bottom protrudes beyond the two-wall corridor");
    }
  }
}

void curvedToeStripTest(weld::Axis axis,int opening,int side,bool gap) {
  weld::DetectorConfig config;config.longitudinalAxis=axis;config.interiorOpeningDirection=opening;
  config.sliceStep=2;
  const int a=int(axis),u=a==0?1:0,v=a==2?1:2;
  const float sign=side?1.f:-1.f;
  auto point=[&](float station,float across,float height) {
    weld::PointT p;p.data[a]=station;p.data[u]=across;p.data[v]=opening*height;return p;
  };
  weld::CloudT::Ptr input(new weld::CloudT);
  std::vector<weld::SurfaceSection> sections;
  for(int f=0;f<8;++f) {
    const float station=2.f*f;
    weld::SurfaceSection s;s.station=station;s.left=point(station,-8,0);s.right=point(station,8,0);
    s.leftNormal=s.rightNormal=Eigen::Vector2f(1,0);sections.push_back(s);
    for(int j=-80;j<=80;++j)input->push_back(point(station,.1f*j,0));
    for(int j=1;j<=150;++j) {
      const float h=.1f*j,x=8+2.5f*(1-std::exp(-h/3));
      if(!gap || h<1 || h>3)input->push_back(point(station,sign*x,h));
      input->push_back(point(station,-sign*8,h));
      input->push_back(point(station,sign*(x+15),h)); // disconnected parallel wall
    }
  }
  weld::SeamResult result;
  weld::extractSurfaceClouds(input,sections,config,result);
  checkMeasured(*input,result);
  int nearToe=0;
  for(const auto &p:*(side?result.rightSurface:result.leftSurface)) {
    const float h=opening*p.data[v];
    require(std::abs(p.data[u])<12,"Curved toe recovery reached a detached wall");
    require(h>0,"Curved toe recovery stole flat floor");
    if(h<1)++nearToe;
    if(gap)require(h<1 || h>3,"Curved toe recovery invented missing measurements");
  }
  if(gap)require(nearToe==0,"Curved toe recovery crossed an unsupported gap to the wall");
  else require(nearToe>20,"Curved measured wall foot was rejected by the broad plane");
}

void measuredFootContinuationTest() {
  weld::DetectorConfig config;
  weld::CloudT candidate,physical;
  std::vector<float> stations;
  std::vector<weld::SurfaceEnd> feet;
  for(int i=0;i<=400;++i) {
    const float s=.5f*i,hook=.03f*std::pow(std::max(0.f,s-170),2);
    weld::PointT p{s,8,hook};physical.push_back(p);
    weld::SurfaceEnd foot;foot.valid=true;foot.station=s;foot.point=p;feet.push_back(foot);
    p.y-=std::min(18.f,std::max(0.f,(s-120)*1.2f));
    candidate.push_back(p);stations.push_back(s);
  }
  // Missing single frames and a detached contact must not terminate the run,
  // nor may a replacement point be generated at those indices.
  feet[250].valid=feet[251].valid=false;feet[270].point.y+=30;
  const auto recovered=weld::selectMeasuredFootContinuation(candidate,stations,feet,240,config);
  require(!recovered.recovered.empty() && recovered.kept.back()==400,"Supported wall-foot tail was not recovered");
  for(auto index:{250u,251u,270u})
    require(std::find(recovered.kept.begin(),recovered.kept.end(),index)==recovered.kept.end(),"Unsupported foot was retained or invented");
  require(weld::selectMeasuredFootContinuation(physical,stations,feet,340,config).kept.empty(),
          "A genuine physical hook was misidentified as candidate drift");
  auto remote=feet;for(auto &foot:remote)foot.point.y+=40;
  require(weld::selectMeasuredFootContinuation(candidate,stations,remote,240,config).kept.empty(),
          "Recovery jumped to a disconnected parallel wall");
  auto gap=feet;for(int i=300;i<=320;++i)gap[i].valid=false;
  const auto broken=weld::selectMeasuredFootContinuation(candidate,stations,gap,240,config);
  require(!broken.kept.empty() && broken.kept.back()==299,"Recovery crossed a sustained measurement gap");
  // Once repaired, the real right hook still trims independently. The left
  // track is straight and must remain complete.
  weld::SeamResult result;std::vector<float> rightStations;
  for(std::size_t i=0;i<stations.size();++i)result.leftSeam->push_back({stations[i],-8,0});
  for(auto i:recovered.kept) {
    result.rightSeam->push_back(physical[i]);rightStations.push_back(stations[i]);
  }
  weld::trimSeamEndpoints(result,config,&stations,nullptr,&rightStations);
  require(result.leftSeam->size()==401 && result.rightSeam->back().x>165 && result.rightSeam->back().x<195,
          "Recovered track did not stop independently at the real physical hook");
}

void endpointSurfaceTest() {
  weld::DetectorConfig config;
  config.limitSurfaceEndsToTracks=false;
  config.extractMainBody=false;config.inputVoxelSize=0;config.voxelSize=0;config.outlierMeanK=0;
  config.useSurfaceRefinement=false;config.trimUnstableEndpoints=false;
  weld::CloudT::Ptr input(new weld::CloudT);
  std::vector<weld::ScanFrameRange> frames;
  auto bend=[](float x){const float t=std::max(0.f,x-170);return .015f*t*t;};
  for(int f=0;f<=240;++f) {
    if(f==210)continue; // Missing measurements must remain absent in surfaces.
    const auto begin=input->size();
    for(int j=-250;j<=250;++j) {
      const float y=.1f*j,z=std::min(16.f,std::max(0.f,(std::abs(y)-8)*2));
      input->push_back({float(f),y+bend(float(f)),z});
    }
    frames.push_back({begin,input->size()-begin,1000+f});
  }
  const auto untrimmed=weld::SeamDetector().detect(input,config,nullptr,nullptr,nullptr,&frames);
  config.trimUnstableEndpoints=true;
  const auto trimmed=weld::SeamDetector().detect(input,config,nullptr,nullptr,nullptr,&frames);
  require(untrimmed.ok() && trimmed.ok(),"Bent endpoint fixture did not detect");
  require(trimmed.sideEndpointTrimmedEnd[0]>0 && trimmed.sideEndpointTrimmedEnd[1]>0,
          "Bent endpoint fixture did not exercise trajectory cuts");
  checkMeasured(*input,trimmed);
  for(int side=0;side<2;++side) {
    const auto &cut=side?*trimmed.rightSeam:*trimmed.leftSeam;
    const auto &all=side?*untrimmed.rightSeam:*untrimmed.leftSeam;
    const auto begin=trimmed.sideEndpointTrimmedStart[side];
    require(cut.size()+begin+trimmed.sideEndpointTrimmedEnd[side]==all.size(),"Surface recovery changed trajectory cuts");
    for(std::size_t i=0;i<cut.size();++i)
      require(key(cut[i])==key(all[begin+i]),"Surface recovery moved a trimmed trajectory point");
  }
  for(int face=0;face<3;++face) {
    const auto &before=face==0?*untrimmed.leftSurface:face==1?*untrimmed.rightSurface:*untrimmed.weldSurface;
    const auto &after=face==0?*trimmed.leftSurface:face==1?*trimmed.rightSurface:*trimmed.weldSurface;
    require(after.size()>=.95*before.size(),"Trajectory endpoint cuts discarded measured surface sections");
    std::set<Key> originalSurface;
    for(const auto &p:before)originalSurface.insert(key(p));
    const float retainedEnd=std::min(trimmed.leftSeam->back().x,trimmed.rightSeam->back().x);
    float last=-1;
    for(std::size_t i=0;i<after.size();++i) {
      if(after[i].x<=retainedEnd)require(originalSurface.count(key(after[i]))!=0,"Endpoint recovery altered retained surface labels");
      require(after[i].x!=210,"Surface recovery filled an unmeasured frame");
      last=std::max(last,after[i].x);
      if(face==2)require(std::abs(after[i].z)<1.e-4f,"Recovered floor contains upper walls");
    }
    require(last>=230,"Supported bent surface endpoint was lost");
  }
  require(trimmed.weldSurfaceModel.valid(),"Recovered endpoint floor model was not fitted");
  config.limitSurfaceEndsToTracks=true;
  const auto limited=weld::SeamDetector().detect(input,config,nullptr,nullptr,nullptr,&frames);
  require(limited.ok(),"Synchronized endpoint fixture did not detect");
  checkMeasured(*input,limited);
  for(int side=0;side<2;++side) {
    const auto &before=side?*trimmed.rightSeam:*trimmed.leftSeam;
    const auto &after=side?*limited.rightSeam:*limited.leftSeam;
    require(before.size()==after.size(),"Surface limits changed trajectory length");
    for(std::size_t i=0;i<before.size();++i)require(key(before[i])==key(after[i]),"Surface limits moved trajectory");
    for(const auto &p:*(side?limited.rightSurface:limited.leftSurface))
      require(p.x<=after.back().x+1.e-4f,"Side surface exceeded its trajectory end");
  }
  require(limited.weldSurface->size()==trimmed.weldSurface->size(),"Wall endpoint limits discarded measured floor");
  for(std::size_t i=0;i<trimmed.weldSurface->size();++i)
    require(key((*limited.weldSurface)[i])==key((*trimmed.weldSurface)[i]),"Wall endpoint limits changed floor points");
  require(limited.weldSurfaceModel.valid() &&
          limited.weldSurfaceModel.patches.back().lastFrame==trimmed.weldSurfaceModel.patches.back().lastFrame,
          "Wall endpoint limits shortened recovered floor model");
  std::cout<<"Endpoint trim retains measured surfaces and models; trajectory cuts and frame gaps preserved.\n";
}

void rejectedBoundaryTest() {
  for(int opening:{-1,1})for(bool single:{false,true}) {
    weld::DetectorConfig config;config.interiorOpeningDirection=opening;
    weld::CloudT::Ptr input(new weld::CloudT);
    std::vector<weld::SurfaceSection> sections;
    for(int f=0;f<8;++f) {
      const float x=2.f*f;
      for(int j=-300;j<=300;++j) {
        input->push_back({x,.1f*j,0});
        input->push_back({x,.1f*j,-12.f*opening}); // detached floor
      }
      for(int side:{-1,1})if(!single || side==1)for(int j=1;j<=150;++j) {
        const float height=.1f*j;
        input->push_back({x,side*(30+.5f*height),opening*height});
        input->push_back({x,side*(50+.5f*height),opening*height}); // detached wall
      }
      weld::SurfaceSection section;section.station=x;
      // Rejected track points lie on a weak interior feature, far from the
      // true walls. They cannot be trusted as toe patches or as floor limits.
      section.left={x,-4,-4.f*opening};section.right={x,4,-4.f*opening};
      section.leftBoundaryUncertain=section.rightBoundaryUncertain=true;
      sections.push_back(section);
    }
    weld::SeamResult result;
    std::vector<std::array<weld::SurfaceEnd,2>> feet;
    weld::extractSurfaceClouds(input,sections,config,result,nullptr,false,nullptr,&feet);
    require(feet.size()==sections.size(),"Measured feet lost section alignment");
    std::set<Key> original;
    for(const auto &p:*input)original.insert(key(p));
    for(const auto &pair:feet)for(int side=0;side<2;++side) {
      require(pair[side].valid==(!single || side==1),"Measured foot invented a missing wall or lost a supported one");
      if(pair[side].valid) {
        require(original.count(key(pair[side].point))!=0,"Recovered wall foot was not measured");
        require(std::abs(std::abs(pair[side].point.y)-30)<3 && pair[side].point.z==0,
                "Recovered wall foot follows an interior ripple or detached sheet");
      }
    }
    std::cout<<"Rejected boundary opening="<<opening<<" single="<<single<<" surfaces="<<result.leftSurface->size()<<'/'<<result.rightSurface->size()<<'/'<<result.weldSurface->size()<<'\n';
    require(result.weldSurface->size()>4400,"Rejected interior boundary hid the measured broad floor");
    require(result.rightSurface->size()>800,"Rejected boundary hid the measured right wall");
    if(single)require(result.leftSurface->empty(),"Rejected boundary created a missing wall");
    else checkMeasured(*input,result);
    for(const auto &p:*result.weldSurface)
      require(std::abs(p.y)<=30.001f && p.z==0,"Boundary recovery included a detached sheet or wall");
    for(const auto &cloud:{result.leftSurface,result.rightSurface})for(const auto &p:*cloud)
      require(std::abs(p.y)<40,"Boundary recovery selected a detached parallel wall");
  }
}
}
int main(int argc,char **argv) {
  try {
    if(argc==2 && std::string(argv[1])=="--endpoint-surfaces") {
      measuredFootContinuationTest();
      for(auto axis:{weld::Axis::X,weld::Axis::Y,weld::Axis::Z})
        for(int opening:{-1,1})for(int side:{0,1})for(bool gap:{false,true})curvedToeStripTest(axis,opening,side,gap);
      for(auto axis:{weld::Axis::X,weld::Axis::Y,weld::Axis::Z})
        for(int opening:{-1,1})for(bool framed:{false,true})geometryTest(axis,opening,0,framed,true);
      rejectedBoundaryTest();
      endpointSurfaceTest();
      return 0;
    }
    if(argc==2 && std::string(argv[1])=="--single-wall-floor") {
      for(auto axis:{weld::Axis::X,weld::Axis::Y,weld::Axis::Z})
        for(int opening:{-1,1})for(int side:{0,1})singleWallFloorTest(axis,opening,side);
      std::cout<<"Single-wall measured floor, model, axes, gaps and detached sheets passed.\n";
      return 0;
    }
    if(argc==3) {
      std::string error;
      const auto input=weld::PointCloudIO::load(argv[1],&error);
      require(input && !input->empty(), "Full scan could not be loaded");
      const auto config=weld::AlgorithmParameters().get(weld::ProfileLocalizationMethod::CurvatureCandidates);
      std::cout<<"Default cloud-only full scan: "<<input->size()<<" input points"<<std::endl;
      const auto result=weld::SeamDetector().detect(input,config);
      std::cout<<result.summary<<std::endl;
      require(result.ok(), "Full scan detection failed");
      for(const auto &stage:result.stageMilliseconds) std::cout<<stage.first<<": "<<stage.second<<" ms\n";
      const std::pair<const char *,weld::CloudT::Ptr> outputs[]={
        {"left_surface",result.leftSurface},{"right_surface",result.rightSurface},{"weld_surface",result.weldSurface}};
      std::vector<bool> used(input->size(),false);
      for(const auto &item:outputs) {
        require(!item.second->empty(), "Full scan surface empty");
        std::size_t sourceIndex=0;
        for(const auto &p:*item.second) {
          while(sourceIndex<input->size() && key((*input)[sourceIndex])!=key(p)) ++sourceIndex;
          require(sourceIndex<input->size(), "Full scan surface is not an original ordered subset");
          require(!used[sourceIndex], "Full scan surface outputs overlap");
          used[sourceIndex++]=true;
        }
        std::ofstream file(std::string(argv[2])+"_"+item.first+".pcd");
        require(bool(file), "Cannot write surface validation output");
        file<<"# .PCD v0.7\nVERSION 0.7\nFIELDS x y z\nSIZE 4 4 4\nTYPE F F F\nCOUNT 1 1 1\nWIDTH "
            <<item.second->size()<<"\nHEIGHT 1\nVIEWPOINT 0 0 0 1 0 0 0\nPOINTS "<<item.second->size()<<"\nDATA ascii\n";
        file<<std::setprecision(9);
        for(const auto &p:*item.second)file<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
        require(bool(file), "Surface validation write failed");
      }
      return 0;
    }
    for(auto axis : {weld::Axis::X,weld::Axis::Y,weld::Axis::Z})
      for(int opening : {-1,1}) for(int shape : {0,1,2}) for(bool framed : {false,true})
        geometryTest(axis,opening,shape,framed);
    detectorTest();
    for(auto axis:{weld::Axis::X,weld::Axis::Y,weld::Axis::Z})
      for(bool framed:{false,true})boundedTwoWallFloorTest(axis,framed);
    roughWeldTest();
    shallowTipTest();
    for(auto axis:{weld::Axis::X,weld::Axis::Y,weld::Axis::Z})
      for(int opening:{-1,1}) {
        shortTransitionTest(axis,opening);
        for(int side:{0,1})sparseSideKeepsFloorTest(axis,opening,side);
        for(int side:{0,1})singleWallFloorTest(axis,opening,side);
      }
    std::cout<<"Three measured surfaces: geometry, coverage, axes, reversed tilted frames, gaps and defaults passed.\n";
    return 0;
  } catch(const std::exception &error) { std::cerr<<error.what()<<'\n';return 1; }
}
