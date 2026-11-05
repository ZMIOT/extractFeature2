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
  config.profileCurvatureRadius=4;
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
  for(std::size_t i=1;i<model.patches.size();++i) {
    const auto &a=model.patches[i-1],&b=model.patches[i];
    if(a.lastFrame!=b.firstFrame)continue;
    for(int j=0;j<=100;++j) {
      Eigen::Vector3d p,q;
      require(a.evaluate(1,j/100.,p) && b.evaluate(0,j/100.,q) && (p-q).norm()<1.e-8,
              "Adjacent fitted subdivisions have a crack");
    }
  }
  // An explicitly invalid face must still split independently.
  for(auto index:groups[10])labels[index]=2;
  sections[10].leftValid=false;
  const auto missing=weld::fitSurfaceModel(cloud,groups,sections,labels,1,config);
  for(const auto &patch:missing.patches)
    require(!(patch.firstFrame<10 && patch.lastFrame>10),"Missing side was bridged");
}
void holesTest(bool framed,bool vertical) {
  weld::CloudT cloud;std::vector<weld::SurfaceSection> sections;
  std::vector<std::vector<std::size_t>> groups;std::vector<unsigned char> labels;
  weld::DetectorConfig config;config.sliceStep=1;config.profileCurvatureRadius=4;
  for(int r=0;r<32;++r) {
    if(r==18 || r==19)continue; // absent scan identities
    weld::SurfaceSection s;s.station=float(r);s.frameIndex=framed?r:-1;
    sections.push_back(s);groups.emplace_back();
    if(r==8 || r==9)continue; // whole missing face in present scans
    for(int j=0;j<=80;++j) {
      if(r>=12 && r<=16 && j>=30 && j<=50)continue;
      const float y=.125f*j,h=.001f*r*r+.01f*y*y;
      groups.back().push_back(cloud.size());labels.push_back(1);
      cloud.push_back({framed?-float(r)+.12f*y:float(r),vertical?h:y,vertical?y:h});
    }
  }
  const auto model=weld::fitSurfaceModel(cloud,groups,sections,labels,1,config);
  require(model.patches.size()==1,"Small supported holes split the fitted face");
  require(model.patches.front().sampleCount==cloud.size(),"Synthetic constraints polluted measured sample count");
  for(int i=0;i<=80;++i)for(int j=0;j<=20;++j) {
    Eigen::Vector3d p,n;
    require(model.evaluate(0,i/80.,j/20.,p,&n),"Hole repair produced a singular surface");
    const double y=vertical?p.z():p.y(),h=vertical?p.y():p.z();
    const double x=framed?-p.x()+.12*y:p.x();
    require(std::abs(h-.001*x*x-.01*y*y)<.08,"Hole interpolation left the measured surface");
  }
  for(std::size_t i=0;i<sections.size();++i)if(sections[i].station>=23 && sections[i].station<=28)
    for(auto index:groups[i])labels[index]=0;
  const auto large=weld::fitSurfaceModel(cloud,groups,sections,labels,1,config);
  for(const auto &patch:large.patches)
    require(!(patch.firstStation<23 && patch.lastStation>28),"Large unsupported gap was filled");
}
void denseFrameHoleTest() {
  weld::CloudT cloud;std::vector<weld::SurfaceSection> sections;
  std::vector<std::vector<std::size_t>> groups;std::vector<unsigned char> labels;
  weld::DetectorConfig config;config.sliceStep=.25f;config.profileCurvatureRadius=8;
  for(int r=0;r<96;++r) {
    weld::SurfaceSection section;section.station=.25f*r;section.frameIndex=r;
    sections.push_back(section);groups.emplace_back();
    if(r>=35 && r<=52)continue;
    for(int j=0;j<=40;++j) {
      groups.back().push_back(cloud.size());labels.push_back(1);
      cloud.push_back({section.station,.25f*j,.001f*section.station*section.station});
    }
  }
  const auto model=weld::fitSurfaceModel(cloud,groups,sections,labels,1,config);
  require(model.patches.size()==1,"Dense scan hole was split by frame count instead of physical distance");
  for(int i=0;i<=100;++i) {
    Eigen::Vector3d p,n;
    require(model.evaluate(0,i/100.,.5,p,&n) && std::abs(p.z()-.001*p.x()*p.x())<.03,
            "Physical-neighborhood repair changed surface shape");
  }
}
void wallFootGapTest(weld::Axis axis,int opening,bool framed,bool curved=false) {
  weld::CloudT cloud;std::vector<weld::SurfaceSection> sections;
  std::vector<std::vector<std::size_t>> groups;std::vector<unsigned char> labels;
  weld::DetectorConfig config;config.longitudinalAxis=axis;config.sliceStep=1;
  config.interiorOpeningDirection=opening;
  const int a=int(axis),b=a==0?1:0,c=a==2?1:2;
  auto point=[&](float station,float h) {
    weld::PointT p;p.data[b]=8+(curved?.2f*h*h:.5f*h);p.data[c]=opening*h;
    p.data[a]=framed?-station+.15f*p.data[b]:station;return p;
  };
  for(int r=0;r<40;++r) {
    weld::SurfaceSection s;s.station=float(r);s.frameIndex=framed?r:-1;
    s.right=point(float(r),0);s.rightNormal=Eigen::Vector2f(1,curved?0:-.5f*opening).normalized();
    sections.push_back(s);groups.emplace_back();
    const int begin=r>=12 && r<=18?30:2;
    for(int j=begin;j<=100;++j) {
      groups.back().push_back(cloud.size());labels.push_back(2);cloud.push_back(point(float(r),.1f*j));
    }
  }
  const auto model=weld::fitSurfaceModel(cloud,groups,sections,labels,2,config);
  require(model.patches.size()==1 && model.rejectedPatchCount==0,"Supported wall-foot notch did not close continuously");
  const auto &patch=model.patches.front();
  require(patch.sampleCount==cloud.size(),"Inferred wall strip counted as measured points");
  for(int i=0;i<=100;++i) {
    Eigen::Vector3d p,q,n;
    require(patch.evaluate(i/100.,0,p,&n) && patch.evaluate(i/100.,1,q),"Repaired wall has a singular surface");
    const auto &foot=std::abs(p[c])<std::abs(q[c])?p:q;
    require(std::abs(foot[c])<.2 && std::abs(foot[b]-8)<.2,"Wall repair left a gap at the supported foot");
    for(int j=0;j<=20;++j) {
      require(patch.evaluate(i/100.,j/20.,p),"Repaired wall cannot be evaluated");
      require(opening*p[c]>=-.2 && opening*p[c]<=10.2,"Repaired wall extends past its foot or outer rim");
      const double expected=curved?.2*p[c]*p[c]:.5*opening*p[c];
      require(std::abs(p[b]-8-expected)<(curved?.65:.05),"Gap fill left the neighboring wall surface");
    }
  }
  // A remote boundary is not a license to bridge a large unsupported strip.
  for(auto &s:sections)s.right=point(s.station,-30);
  const auto remote=weld::fitSurfaceModel(cloud,groups,sections,labels,2,config);
  require(remote.valid(),"Rejected remote foot discarded the measured wall");
  for(const auto &p:remote.patches)for(double u:{0.,.5,1.})for(double v:{0.,1.}) {
    Eigen::Vector3d q;require(p.evaluate(u,v,q) && opening*q[c]>-1,"Wall repair bridged a remote boundary");
  }
}

void outerRimRepairTest(weld::Axis axis,int opening,bool framed,int missingKind) {
  weld::CloudT cloud;std::vector<weld::SurfaceSection> sections;
  std::vector<std::vector<std::size_t>> groups;std::vector<unsigned char> labels;
  weld::DetectorConfig config;config.longitudinalAxis=axis;config.interiorOpeningDirection=opening;
  config.sliceStep=1;config.profileCurvatureRadius=8;
  const int a=int(axis),b=a==0?1:0,c=a==2?1:2;
  auto point=[&](float station,float h) {
    weld::PointT p;p.data[b]=8+.5f*h;p.data[c]=opening*h;
    p.data[a]=framed?-station+.15f*p.data[b]:station;return p;
  };
  for(int r=0;r<40;++r) {
    weld::SurfaceSection s;s.station=float(r);s.frameIndex=r;
    s.right=point(float(r),0);s.rightNormal=Eigen::Vector2f(1,-.5f*opening).normalized();
    const bool missing=missingKind==1?(r>=10 && r<=29):missingKind==2?r<=5:(r>=12 && r<=18);
    if(missingKind==3 && missing)s.rightBoundaryUncertain=true;
    sections.push_back(s);groups.emplace_back();
    for(int j=2;j<=(missing?25:100);++j) {
      groups.back().push_back(cloud.size());labels.push_back(2);cloud.push_back(point(float(r),.1f*j));
    }
  }
  config.surfaceRimRepairLength=0;
  const auto original=weld::fitSurfaceModel(cloud,groups,sections,labels,2,config);
  require(original.valid(),"Rim fixture did not fit");
  config.surfaceRimRepairLength=20;
  const auto repaired=weld::fitSurfaceModel(cloud,groups,sections,labels,2,config);
  require(repaired.valid() && repaired.rejectedPatchCount==0,"Rim repair lost supported patches");
  if(missingKind!=0) {
    require(repaired.patches.size()==original.patches.size(),"Protected rim boundary changed patch coverage");
    for(std::size_t i=0;i<original.patches.size();++i)
      require(repaired.patches[i].controlPoints==original.patches[i].controlPoints,
              "Rim repair altered an endpoint, long gap or uncertain boundary");
  }
  bool lowMiddle=false;
  for(const auto &patch:repaired.patches)for(int i=0;i<=80;++i) {
    Eigen::Vector3d p,q;
    require(patch.evaluate(i/80.,0,p) && patch.evaluate(i/80.,1,q),"Rim edge is singular");
    const auto &rim=opening*p[c]>opening*q[c]?p:q;
    const auto &foot=opening*p[c]>opening*q[c]?q:p;
    const double station=framed?-rim[a]+.15*rim[b]:rim[a];
    require(std::abs(foot[c])<.3,"Rim repair displaced the wall foot");
    if(missingKind==0)require(opening*rim[c]<=10.4,"Rim repair overshot supported outer boundary");
    require(std::abs(rim[b]-8-.5*opening*rim[c])<.1,"Rim repair left neighboring wall geometry");
    if(missingKind==0)require(opening*rim[c]>9.5,"Bounded outer rim notch was not repaired");
    else if((missingKind==2?station<3:missingKind==1?(station>15 && station<24):(station>13 && station<17)) && opening*rim[c]<5)
      lowMiddle=true;
  }
  if(missingKind!=0)require(lowMiddle,"Rim repair bridged an endpoint, long gap or uncertain boundary");
  if(missingKind==0) {
    std::size_t count=0;for(const auto &p:repaired.patches)count+=p.sampleCount;
    require(repaired.patches.size()==1 && count==cloud.size(),"Inferred rim samples were counted as measurements");
    Eigen::Vector3d p,q;
    require(original.patches.front().evaluate(.38,0,p) && original.patches.front().evaluate(.38,1,q),"Original rim evaluation failed");
    require(std::max(opening*p[c],opening*q[c])<8,"Disabled rim option altered the original notch");
  }
}

void shortSupportedStripTest() {
  weld::CloudT cloud;std::vector<weld::SurfaceSection> sections;
  std::vector<std::vector<std::size_t>> groups;std::vector<unsigned char> labels;
  weld::DetectorConfig config;config.sliceStep=2;
  config.profileWallTolerance=config.hullSurfaceDistance=.05f;
  for(int r=0;r<5;++r) {
    weld::SurfaceSection s;s.station=2.f*r;sections.push_back(s);groups.emplace_back();
    const float width=r==2?13.f:10.f;
    for(int j=0;j<=100;++j) {
      groups.back().push_back(cloud.size());labels.push_back(1);
      cloud.push_back({s.station,width*j/100.f,0});
    }
  }
  const auto model=weld::fitSurfaceModel(cloud,groups,sections,labels,1,config);
  require(model.patches.size()==1 && model.rejectedPatchCount==0,
          "Five supported sections were discarded because the coarse cubic basis underfit a local boundary");
  const auto &patch=model.patches.front();
  require(patch.countU==5 && patch.rmsError<=.1 && patch.sampleCount==cloud.size(),
          "Short-patch retry relaxed quality or invented measured samples");
  require(patch.firstStation==0 && patch.lastStation==8,"Short-patch retry lost supported coverage");
  for(double u:{0.,.25,.5,.75,1.})for(double v:{0.,.5,1.}) {
    Eigen::Vector3d p,n;
    require(patch.evaluate(u,v,p,&n) && std::abs(p.z())<1.e-7,"Short strip retry left the measured plane");
  }
}
}
int main() {
  try {
    for(auto axis:{weld::Axis::X,weld::Axis::Y,weld::Axis::Z})
      for(bool vertical:{false,true})for(bool framed:{false,true})curvedTest(axis,vertical,framed,true);
    curvedTest(weld::Axis::X,false,false,false);
    noisyBoundaryTest();
    for(bool framed:{false,true})for(bool vertical:{false,true})holesTest(framed,vertical);
    denseFrameHoleTest();
    shortSupportedStripTest();
    for(auto axis:{weld::Axis::X,weld::Axis::Y,weld::Axis::Z})
      for(int opening:{-1,1})for(bool framed:{false,true})outerRimRepairTest(axis,opening,framed,0);
    for(int kind:{1,2,3})outerRimRepairTest(weld::Axis::X,1,false,kind);
    for(auto axis:{weld::Axis::X,weld::Axis::Y,weld::Axis::Z})
      for(int opening:{-1,1})for(bool framed:{false,true})wallFootGapTest(axis,opening,framed);
    wallFootGapTest(weld::Axis::X,1,true,true);
    std::cout<<"B-spline analytic geometry, normals, gaps, reversed tilted scans, vertical walls and invalid input passed.\n";
    return 0;
  }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
