#include "SurfaceClouds.h"
#include "SeamDetector.h"
#include "AlgorithmParameters.h"
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

void geometryTest(weld::Axis axis, int opening, int wallShape, bool framed) {
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
}
int main(int argc,char **argv) {
  try {
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
    roughWeldTest();
    shallowTipTest();
    for(auto axis:{weld::Axis::X,weld::Axis::Y,weld::Axis::Z})
      for(int opening:{-1,1})shortTransitionTest(axis,opening);
    std::cout<<"Three measured surfaces: geometry, coverage, axes, reversed tilted frames, gaps and defaults passed.\n";
    return 0;
  } catch(const std::exception &error) { std::cerr<<error.what()<<'\n';return 1; }
}
