#include "SeamDirection.h"
#include "SeamDetector.h"
#include "PointCloudIO.h"
#include <Eigen/Geometry>
#include <iostream>
#include <limits>
#include <stdexcept>

void require(bool condition,const char *message) {if (!condition) throw std::runtime_error(message);}
weld::CloudT::Ptr groove(weld::Axis axis,int opening,std::vector<weld::ScanFrameRange> &frames) {
  weld::CloudT::Ptr cloud(new weld::CloudT);
  const int a=int(axis),u=a==0?1:0,v=a==2?1:2;
  for (int f=0;f<41;++f) {
    const auto begin=cloud->size();
    for (int i=-160;i<=160;++i) {
      const float x=.2f*i;
      weld::PointT p;
      p.data[a]=5.f*f;
      p.data[u]=x;
      p.data[v]=opening*std::min(16.f,std::max(0.f,2.f*(std::abs(x)-10.f)));
      cloud->push_back(p);
    }
    frames.push_back({begin,cloud->size()-begin,f});
  }
  return cloud;
}
int main(int argc,char **argv) {
  try {
    if (argc>1) {
      for (int i=1;i<argc;++i) {
        std::string error;std::vector<weld::ScanFrameRange> frames;
        const auto cloud=weld::PointCloudIO::load(argv[i],&error,&frames);
        require(cloud && !cloud->empty(),"Fixture failed to load");
        const auto d=weld::estimateSeamDirection(*cloud,weld::DetectorConfig{},&frames);
        std::cout<<argv[i]<<" axis="<<int(d.axis)<<" pca="<<d.principalDirection.transpose()<<" gap="<<d.principalGap
                 <<" opening="<<d.opening<<" votes="<<d.positiveVotes<<"/"<<d.negativeVotes<<'\n';
      }
      return 0;
    }
    for (const auto axis:{weld::Axis::X,weld::Axis::Y,weld::Axis::Z}) for (int sign:{1,-1}) {
      std::vector<weld::ScanFrameRange> frames;
      const auto cloud=groove(axis,sign,frames);
      weld::DetectorConfig config;
      config.longitudinalAxis=axis==weld::Axis::X?weld::Axis::Z:weld::Axis::X;
      config.interiorOpeningDirection=-sign;
      for (bool framed:{false,true}) {
        const auto d=weld::estimateSeamDirection(*cloud,config,framed?&frames:nullptr);
        require(d.axisReliable && d.axis==axis,"PCA did not correct the configured axis");
        require(d.openingReliable && d.opening==sign,"Opening vote did not correct the configured sign");
      }
      config.extractMainBody=false;config.inputVoxelSize=0;config.voxelSize=0;config.outlierMeanK=0;
      config.trimUnstableEndpoints=false;
      weld::DebugTrace trace;
      const auto result=weld::SeamDetector().detect(cloud,config,nullptr,nullptr,&trace,&frames,false);
      require(result.ok(),"Automatic directions failed in the complete detection pipeline");
      require(result.resolvedLongitudinalAxis==axis && result.resolvedOpeningDirection==sign,"Resolved directions not propagated");
      require(trace.config.longitudinalAxis==axis && trace.config.interiorOpeningDirection==sign,"Debug config retained stale directions");
      require((*cloud)[321].data[int(axis)]==5,"Direction detection changed input coordinates");
    }
    weld::DetectorConfig config;
    for (const auto *name:{"curvature_2_frames_2040_2110.txt","curvature_2_frames_2042_2058.txt"}) {
      std::string error;std::vector<weld::ScanFrameRange> frames;
      const auto cloud=weld::PointCloudIO::load(std::string(SEAM_TEST_DATA_DIR)+"/"+name,&error,&frames);
      require(cloud && !cloud->empty(),"Measured direction fixture missing");
      const auto d=weld::estimateSeamDirection(*cloud,config,&frames);
      require(d.axisReliable && d.axis==weld::Axis::X && d.openingReliable && d.opening==1,
              "Short oblique scan confused profile tilt with longitudinal direction");
    }
    weld::CloudT constant;
    for (int i=0;i<100;++i) constant.push_back(weld::PointT(1,2,3));
    constant.push_back(weld::PointT(std::numeric_limits<float>::quiet_NaN(),0,0));
    auto d=weld::estimateSeamDirection(constant,config);
    require(!d.axisReliable && !d.openingReliable,"Degenerate cloud reported a reliable direction");
    std::vector<weld::ScanFrameRange> oneFrame{{0,constant.size(),0}};
    d=weld::estimateSeamDirection(constant,config,&oneFrame);
    require(!d.axisReliable && !d.openingReliable,"Single frame invented a longitudinal progression");
    weld::CloudT plane;
    for (int a=0;a<100;++a) for (int u=-30;u<=30;++u) plane.push_back(weld::PointT(4.f*a,float(u),.3f*u));
    d=weld::estimateSeamDirection(plane,config);
    require(d.axisReliable && !d.openingReliable,"Tilted plane was mistaken for a groove");
    std::vector<weld::ScanFrameRange> frames;
    auto cloud=groove(weld::Axis::X,1,frames);
    const Eigen::Matrix3f rotation=Eigen::AngleAxisf(.35f,Eigen::Vector3f::UnitZ()).toRotationMatrix();
    for (auto &p:*cloud) p=weld::toPoint(rotation*weld::toVector(p)+Eigen::Vector3f(10000,20000,-30000));
    d=weld::estimateSeamDirection(*cloud,config,&frames);
    require(d.axisReliable && d.axis==weld::Axis::X && d.principalDirection.dot(rotation.col(0))>.999f,"Tilted PCA direction is inaccurate");
    config.autoLongitudinalAxis=false;config.autoOpeningDirection=false;config.longitudinalAxis=weld::Axis::Z;
    d=weld::estimateSeamDirection(*cloud,config);
    require(d.axis==weld::Axis::Z && !d.openingReliable,"Manual mode was overridden");
    std::vector<weld::ScanFrameRange> bodyFrames;
    auto noisyBody=groove(weld::Axis::Y,-1,bodyFrames);
    const weld::CloudT sparseBody=*noisyBody;
    for (const auto &p:sparseBody) if (p.y<200) for (int offset=1;offset<5;++offset)
      noisyBody->push_back(weld::PointT(p.x,p.y+offset,p.z));
    for (int i=0;i<50;++i) noisyBody->push_back(weld::PointT(1000+100.f*i,2000,0));
    weld::DetectorConfig bodyConfig;
    require(weld::estimateSeamDirection(*noisyBody,bodyConfig).axis==weld::Axis::X,"Noise fixture did not bias the raw PCA");
    const auto bodyResult=weld::SeamDetector().detect(noisyBody,bodyConfig,nullptr,nullptr,nullptr,nullptr,false);
    std::cout<<"Filtered body: "<<bodyResult.directionSummary<<" removed="<<bodyResult.removedBodyPoints
             <<" points="<<bodyResult.mainBodyCloud->size()<<"; "<<bodyResult.summary<<'\n';
    require(bodyResult.ok() && bodyResult.removedBodyPoints>=50 && bodyResult.resolvedLongitudinalAxis==weld::Axis::Y &&
            bodyResult.resolvedOpeningDirection==-1,"Direction was not estimated from the filtered main body");
    std::cout << "Automatic seam direction tests passed\n";
    return 0;
  } catch (const std::exception &e) {std::cerr << e.what() << '\n';return 1;}
}
