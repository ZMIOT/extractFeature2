#include "SeamDetector.h"
#include "PointCloudIO.h"
#include <Eigen/Eigenvalues>
#include <fstream>
#include <iostream>
#include <chrono>
#include <limits>
#include <algorithm>
#include <iomanip>

int main(int argc, char **argv) {
  if (argc < 3) return 2;
  std::string error;
  std::vector<weld::ScanFrameRange> frames;
  auto input = weld::PointCloudIO::load(argv[1], &error, &frames);
  if (!input || input->empty()) { std::cerr << error; return 3; }
  std::cout << "Loaded " << input->size() << " points, " << frames.size() << " frames" << std::endl;
  weld::DetectorConfig config;
  if (argc > 3) config.minSeamSeparation = std::stof(argv[3]);
  if (argc > 4) config.localAngleWindow = std::stoi(argv[4]);
  if (argc > 5) config.interiorOpeningDirection = std::stoi(argv[5]);
  const bool prepared = argc > 6 && std::string(argv[6]) == "--prepared";
  if (prepared) { config.extractMainBody=false; config.voxelSize=0; config.outlierMeanK=0; }
  for (int i=7;i<argc;++i) if (std::string(argv[i])=="--no-endpoint-trim") config.trimUnstableEndpoints=false;
  for (int i=7;i<argc;++i) if (std::string(argv[i])=="--exhaustive-candidates") config.useWallCandidatePruning=false;
  for (int i=7;i<argc;++i) {
    const std::string arg=argv[i];
    if (arg=="--no-final-smoothing") config.finalSmoothingEnabled=false;
    if (arg=="--no-uniform-width-filter") config.uniformWidthFilterEnabled=false;
    if (arg=="--uniform-width-filter-length" && i+1<argc)config.uniformWidthFilterLength=std::stof(argv[++i]);
    if (arg=="--uniform-width" && i+1<argc) {
      const std::string value=argv[++i];
      if(value!="arc" && value!="chord")return 2;
      config.uniformWidthEnabled=true;
      config.uniformWidthMetric=value=="arc"?weld::UniformWidthMetric::BottomArc:weld::UniformWidthMetric::Chord;
    }
    if (arg=="--central-bead-height" && i+1<argc)
      config.centralBeadMinHeight=std::stof(argv[++i]);
    if (arg=="--minimum-feature-depth" && i+1<argc)
      config.minimumFeatureDepth=std::stof(argv[++i]);
    if (arg=="--central-bead-width" && i+1<argc)
      config.centralBeadMinWidth=std::stof(argv[++i]);
    if (arg=="--localization" && i+1<argc) {
      const std::string mode=argv[++i];
      if (mode=="curve-line") config.profileLocalizationMethod=weld::ProfileLocalizationMethod::CurveLineIntersection;
      else if (mode=="curvature") config.profileLocalizationMethod=weld::ProfileLocalizationMethod::CurvatureCandidates;
      else { std::cerr<<"Invalid localization method\n"; return 2; }
    }
    if (arg=="--profile-turn-location" && i+1<argc) {
      const std::string mode=argv[++i];
      if (mode=="peak") config.profileTurnLocation=weld::ProfileTurnLocation::Peak;
      else if (mode=="region-median") config.profileTurnLocation=weld::ProfileTurnLocation::RegionMedian;
      else if (mode=="turn-median") config.profileTurnLocation=weld::ProfileTurnLocation::TurnMedian;
      else if (mode=="floor-side-boundary") config.profileTurnLocation=weld::ProfileTurnLocation::FloorSideBoundary;
      else if (mode=="wall-turn-centroid") config.profileTurnLocation=weld::ProfileTurnLocation::WallTurnCentroid;
      else { std::cerr<<"Invalid turn location mode\n"; return 2; }
    }
    if (arg=="--profile-merge-frames" && i+1<argc) {
      config.mergeProfileFrames=true;
      config.profileMergeFrameCount=std::stoi(argv[++i]);
    }
    if (arg=="--profile-curvature-radius" && i+1<argc)
      config.profileCurvatureRadius=std::stof(argv[++i]);
    if (arg=="--profile-turn-region-fraction" && i+1<argc)
      config.profileTurnRegionFraction=std::stof(argv[++i]);
    if (arg=="--profile-turn-region-length" && i+1<argc)
      config.profileTurnRegionMaxLength=std::stof(argv[++i]);
    if (arg=="--profile-wall-tolerance" && i+1<argc)
      config.profileWallTolerance=std::stof(argv[++i]);
  }
  weld::DebugTrace trace;
  const auto start = std::chrono::steady_clock::now();
  const auto result = weld::SeamDetector().detect(input, config, nullptr, nullptr, &trace, &frames);
  for (int i=7;i+3<argc;++i) if (std::string(argv[i])=="--export-frames") {
    const int first=std::stoi(argv[i+1]), last=std::stoi(argv[i+2]);
    std::ofstream fixture(argv[i+3]);
    fixture << "index x y z target\n" << std::setprecision(9);
    std::size_t index=0;
    for (const auto &slice:trace.slices) if (slice.frameIndex>=first && slice.frameIndex<=last)
      for (const auto &p:slice.sourceRaw)
        fixture << index++ << ' ' << p.x << ' ' << p.y << ' ' << p.z << ' ' << slice.frameIndex << '\n';
  }
  const std::string base = argv[2];
  std::ofstream report(base + "_report.txt"), profiles(base + "_profiles.csv"), pairs(base + "_pairs.csv"), stages(base + "_stages.csv"), support(base + "_support.csv");
  report << result.summary << '\n';
  report << "profile_turn_location," << int(config.profileTurnLocation) << '\n';
  std::ofstream regions(base+"_turn_regions.csv");
  std::ofstream candidateWalls(base+"_candidate_walls.csv");
  candidateWalls << "frame,pair,side,u,v,du,dv,rms\n" << std::setprecision(9);
  for (const auto &slice:trace.slices) for (std::size_t i=0;i<slice.pairs.size();++i)
    for (int side=0;side<2;++side) {
      const auto &wall=slice.pairs[i].candidateWalls[side];
      candidateWalls << slice.frameIndex << ',' << i << ',' << side << ',' << wall.u << ',' << wall.v
                     << ',' << wall.du << ',' << wall.dv << ',' << wall.rms << '\n';
    }
  regions << "frame,begin_u,begin_v,end_u,end_v,peak_u,peak_v,representative_u,representative_v,arc_length,turning_angle\n";
  for (const auto &slice:trace.slices) for (const auto &r:slice.turnRegions)
    regions<<slice.frameIndex<<','<<r.begin.y<<','<<r.begin.z<<','<<r.end.y<<','<<r.end.z<<','
           <<r.peak.y<<','<<r.peak.z<<','<<r.representative.y<<','<<r.representative.z<<','<<r.arcLength<<','<<r.turningAngle<<'\n';
  std::size_t enumerated=0,pruned=0,spikes=0;
  for (const auto &slice:trace.slices) {
    auto e=slice.gates.find("enumerated"),p=slice.gates.find("wall_roi_pruned");
    if (e!=slice.gates.end()) enumerated+=e->second;
    if (p!=slice.gates.end()) pruned+=p->second;
    auto s=slice.gates.find("profile_spikes_removed");
    if (s!=slice.gates.end()) spikes+=s->second;
  }
  report << "enumerated_pairs," << enumerated << "\npruned_pairs," << pruned << '\n';
  report << "profile_spikes_removed," << spikes << '\n';
  for (const auto &stage : result.stageMilliseconds) { report << stage.first << ',' << stage.second << '\n'; std::cout << stage.first << ',' << stage.second << std::endl; }
  report << "total_ms," << std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count() << '\n';
  profiles << "frame,station,u,v\n";
  std::ofstream localProfiles(base+"_local_profiles.csv");
  localProfiles << "frame,station,u,v\n";
  pairs << std::setprecision(9);
  pairs << "frame,station,pair,selected,lu,lv,ru,rv,cost,source,left_curvature,right_curvature,left_wall_distance,right_wall_distance,bottom_midpoint_u,has_bottom_midpoint\n";
  if (argc > 6 && std::string(argv[6]) == "--cache-prepared") {
    std::ofstream cache(base+"_prepared.txt");
    cache << "index x y z target\n" << std::setprecision(9);
    std::size_t index=0;
    for (const auto &slice:trace.slices) for (const auto &p:slice.raw)
      cache << index++ << ' ' << slice.station << ' ' << p.y << ' ' << p.z << ' ' << slice.frameIndex << '\n';
  }
  std::ofstream splines(base+"_splines.csv");
  splines << "frame,u,v,curvature\n";
  std::ofstream sliceStatus(base+"_slice_status.csv");
  std::ofstream gateRows(base+"_gates.csv");
  gateRows << "frame,gate,count\n";
  sliceStatus << "frame,status,pairs,feature_depth";
  const std::vector<std::string> gateNames{
      "invalid_fit","slope_or_broad_angle","invalid_local_fit",
      "local_angle","turn_direction","wrong_opening_direction",
      "parallel_lines","intersection_or_separation","clearance",
      "curve_line_no_intersection","curve_line_geometry","candidate_crosses_bottom_midpoint","feature_depth_cutoff","accepted"};
  for (const auto &name:gateNames) sliceStatus << ',' << name;
  sliceStatus << '\n';
  for (const auto &slice : trace.slices) {
    for (const auto &gate : slice.gates)
      gateRows << slice.frameIndex << ',' << gate.first << ',' << gate.second << '\n';
    sliceStatus << slice.frameIndex << ',' << slice.status << ',' << slice.pairs.size() << ',' << slice.featureDepth;
    for (const auto &name:gateNames) {
      const auto found=slice.gates.find(name);
      sliceStatus << ',' << (found==slice.gates.end()?0:found->second);
    }
    sliceStatus << '\n';
    for (std::size_t i=0;i<slice.fittedProfile.size();++i) { const auto &p=slice.fittedProfile[i]; splines<<slice.frameIndex<<','<<p.y<<','<<p.z<<','<<slice.fittedCurvature[i]<<'\n'; }
    for (const auto &p:slice.localizationProfile) localProfiles << slice.frameIndex << ',' << slice.station << ',' << p.y << ',' << p.z << '\n';
    for (const auto &p : slice.profile) profiles << slice.frameIndex << ',' << slice.station << ',' << p.u << ',' << p.v << '\n';
    for (std::size_t i=0;i<slice.pairs.size();++i) {
      const auto &p=slice.pairs[i];
      auto wallDistance=[&](int side) { const auto &l=p.source.find("spline_curvature")!=std::string::npos ? p.candidateWalls[side] : p.lines[side?2:0]; const auto &q=p.points[side]; return std::fabs(-l.dv*(q.y-l.u)+l.du*(q.z-l.v)); };
      pairs << slice.frameIndex << ',' << slice.station << ',' << i << ',' << (slice.selected==int(i)) << ',' << p.points[0].y << ',' << p.points[0].z << ',' << p.points[1].y << ',' << p.points[1].z << ',' << p.proposalCost << ',' << p.source << ',' << p.profileCurvatures[0] << ',' << p.profileCurvatures[1] << ',' << wallDistance(0) << ',' << wallDistance(1) << ',' << p.bottomMidpointU << ',' << p.hasBottomMidpoint << '\n';
    }
  }
  stages << "stage,side,station,u,v\n";
  auto dump=[&](const char *name,int side,const weld::PointT &p){stages<<name<<','<<side<<','<<p.x<<','<<p.y<<','<<p.z<<'\n';};
  for(int side=0;side<2;++side) {
    for(const auto&p:trace.selectedPath[side])dump("candidate",side,p);
    for(const auto&p:trace.interpolatedPath[side])dump("interpolated",side,p);
    for(const auto&p:trace.beforeWallPath[side])dump("regularized",side,p);
  }
  std::ofstream walls(base+"_walls.csv");
  walls<<"side,station,normal_u,normal_v,offset,rms,status,method\n";
  for(const auto &p:trace.wallPoints) {
    dump("wall",p.side,p.beforeFloor);dump("floor",p.side,p.after);
    walls<<p.side<<','<<p.before.x<<','<<p.normal.x()<<','<<p.normal.y()<<','<<p.offset<<','<<p.rms<<','<<p.status<<','<<p.wallFitMethod<<'\n';
  }
  support << "side,frame,x,y,z,nearest_raw,normal_x,tracked_wall,tracked_residual\n";
  for(int side=0;side<3;++side) {
    const auto &cloud=side==2?result.centerSeam:side?result.rightSeam:result.leftSeam;
    const auto &ids=side==2?result.centerFrameIndices:side?result.rightFrameIndices:result.leftFrameIndices;
    for(std::size_t i=0;i<ids.size();++i) {
      auto frame=std::lower_bound(frames.begin(),frames.end(),ids[i],[](const weld::ScanFrameRange &f,int id){return f.frameIndex<id;});
      if(frame==frames.end()||frame->frameIndex!=ids[i])return 4;
      Eigen::Vector3d mean=Eigen::Vector3d::Zero();
      double nearest=std::numeric_limits<double>::infinity();
      for(std::size_t j=frame->begin;j<frame->begin+frame->count;++j){mean+=weld::toVector((*input)[j]).cast<double>();nearest=std::min(nearest,double((weld::toVector((*input)[j])-weld::toVector((*cloud)[i])).norm()));}
      mean/=double(frame->count);Eigen::Matrix3d cov=Eigen::Matrix3d::Zero();
      for(std::size_t j=frame->begin;j<frame->begin+frame->count;++j){const Eigen::Vector3d d=weld::toVector((*input)[j]).cast<double>()-mean;cov+=d*d.transpose();}
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(cov);
      const auto &quality=side==2?result.centerPointQuality:side?result.rightPointQuality:result.leftPointQuality;
      const auto &q=quality[i];
      const auto &p=(*cloud)[i];support<<side<<','<<ids[i]<<','<<p.x<<','<<p.y<<','<<p.z<<','<<nearest<<','<<std::abs(eig.eigenvectors()(0,0))<<','<<q.trackedWall<<','<<
          (q.trackedWall ? std::fabs(q.trackedNormal.dot(Eigen::Vector2f(p.y,p.z))-q.trackedOffset) : -1)<<'\n';
    }
  }
  if(result.ok())weld::PointCloudIO::saveFrameSeamCsv(base+"_seams.csv",result,&error);
  std::cout<<result.summary<<std::endl;
  return result.ok()?0:5;
}
