#include "SurfaceFitting.h"
#include <Eigen/Geometry>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool condition,const char *message){if(!condition)throw std::runtime_error(message);}
void curvedTest(weld::Axis axis,bool vertical,bool framed,bool gaps) {
  weld::CloudT input;std::vector<weld::SurfaceSection> sections;
  std::vector<std::vector<std::size_t>> groups;std::vector<unsigned char> labels;
  weld::DetectorConfig config;config.longitudinalAxis=axis;config.sliceStep=1;
  const int a=int(axis),b=a==0?1:0,c=a==2?1:2;
  // Nonplanar surface, with reversed/tilted scans and a large world translation.
  auto point=[&](double x,double y) {
    const double h=.003*x*x+.02*y*y;
    weld::PointT p;p.data[a]=float(1000+(framed?-x+.12*y:x));
    p.data[b]=float(-200+(vertical?h:y));p.data[c]=float(500+(vertical?y:h));return p;
  };
  for(int r=0;r<40;++r) {
    if(gaps && r>=15 && r<20)continue;
    weld::SurfaceSection s;s.station=float(r);s.frameIndex=framed?100+r:-1;sections.push_back(s);
    groups.emplace_back();
    for(int j=0;j<=80;++j) {
      const double y=-5+.125*j;
      groups.back().push_back(input.size());input.push_back(point(r,y));labels.push_back(1);
    }
  }
  auto model=weld::fitSurfaceModel(input,groups,sections,labels,1,config);
  require(model.valid(),"Curved surface fitting failed");
  require(model.patches.size()==(gaps?2:1),"Fit bridged a missing scan range");
  for(const auto &patch:model.patches) {
    require(patch.rmsError<.06,"Curved fit residual too large");
    for(double u:{0.,.17,.5,.91,1.})for(double v:{0.,.23,.5,.87,1.}) {
      Eigen::Vector3d p,n;
      require(patch.evaluate(u,v,p,&n),"Evaluation or analytic normal failed");
      const double y=vertical?p[c]-500:p[b]+200;
      const double x=framed?-(p[a]-1000-.12*y):p[a]-1000;
      const double h=vertical?p[b]+200:p[c]-500;
      require(std::abs(h-.003*x*x-.02*y*y)<.08,"Spline does not reproduce analytic surface");
      require(std::abs(n.norm()-1)<1e-10,"Normal is not normalized");
      if(u>.01 && u<.99 && v>.01 && v<.99) {
        Eigen::Vector3d pu0,pu1,pv0,pv1;
        patch.evaluate(u-1e-5,v,pu0);patch.evaluate(u+1e-5,v,pu1);
        patch.evaluate(u,v-1e-5,pv0);patch.evaluate(u,v+1e-5,pv1);
        const Eigen::Vector3d numerical=(pu1-pu0).cross(pv1-pv0).normalized();
        require(n.dot(numerical)>.999999,"Derivative normal differs from finite differences");
      }
    }
  }
  Eigen::Vector3d sentinel(1,2,3);
  require(!model.evaluate(99,.5,.5,sentinel) && sentinel==Eigen::Vector3d(1,2,3),"Bad patch index changed output");
  require(!model.evaluate(0,-.01,.5,sentinel),"Out-of-domain parameter was extrapolated");
  require(!model.evaluate(0,.5,std::numeric_limits<double>::quiet_NaN(),sentinel),"NaN accepted");
  auto broken=model.patches.front();broken.knotsU.clear();
  require(!broken.evaluate(.5,.5,sentinel),"Malformed spline accepted");
  // One profile only contains a curve, regardless of how many points it has.
  groups.resize(1);sections.resize(1);
  require(!weld::fitSurfaceModel(input,groups,sections,labels,1,config).valid(),"Single scan invented a surface");
  std::cout<<"curved axis="<<a<<" vertical="<<vertical<<" framed="<<framed<<" gaps="<<gaps<<" passed\n";
}
void noisyBoundaryTest() {
  weld::CloudT cloud;std::vector<weld::SurfaceSection> sections;
  std::vector<std::vector<std::size_t>> groups;std::vector<unsigned char> labels;
  weld::DetectorConfig config;config.profileWallTolerance=config.hullSurfaceDistance=.05f;
  for(int r=0;r<50;++r) {
    weld::SurfaceSection section;section.frameIndex=r;section.station=float(r);
    sections.push_back(section);groups.emplace_back();
    for(int j=0;j<=100;++j) {
      const float y=(r<25?5.f:12.f)*j/100;
      const float noise=.01f*std::sin(2.31f*j+1.4f*r)+(j==17 && r%7==0?.3f:0.f);
      groups.back().push_back(cloud.size());labels.push_back(1);cloud.push_back({float(r),y,noise});
    }
  }
  const auto model=weld::fitSurfaceModel(cloud,groups,sections,labels,1,config);
  require(model.valid() && model.patches.size()>1,"Abrupt measured boundary did not subdivide");
  for(const auto &patch:model.patches) {
    require(patch.rmsError<=.101,"Low-quality patch was returned as valid");
    Eigen::Vector3d p;
    require(patch.evaluate(.5,.5,p) && std::abs(p.z())<.05,"Noise/outliers displaced fitted plane");
  }
  // A missing side in a present frame must split this face independently.
  for(auto index:groups[10])labels[index]=2;
  const auto missing=weld::fitSurfaceModel(cloud,groups,sections,labels,1,config);
  for(const auto &patch:missing.patches)
    require(!(patch.firstFrame<10 && patch.lastFrame>10),"Missing side was bridged");
}
}
int main() {
  try {
    for(auto axis:{weld::Axis::X,weld::Axis::Y,weld::Axis::Z})
      for(bool vertical:{false,true})for(bool framed:{false,true})curvedTest(axis,vertical,framed,true);
    curvedTest(weld::Axis::X,false,false,false);
    noisyBoundaryTest();
    std::cout<<"B-spline analytic geometry, normals, gaps, reversed tilted scans, vertical walls and invalid input passed.\n";
    return 0;
  }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
