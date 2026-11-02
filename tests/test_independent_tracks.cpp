#include "SeamDetector.h"
#include "SeamSmoother.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

static void require(bool value,const char *message) {if(!value)throw std::runtime_error(message);}

int main() {
  try {
    weld::DetectorConfig config;
    config.extractMainBody=false;config.inputVoxelSize=0;config.voxelSize=0;config.outlierMeanK=0;
    config.trimUnstableEndpoints=false;config.useSurfaceRefinement=false;
    auto detect=[&](bool onlyLeft,bool mirrored) {
      weld::CloudT::Ptr cloud(new weld::CloudT);
      std::vector<weld::ScanFrameRange> frames;
      for(int frame=0;frame<80;++frame) {
        const auto begin=cloud->size();
        for(int j=-350;j<=350;++j) {
          float y=.1f*j;
          const bool wall=y<0 ? frame>=10 : (!onlyLeft && frame<50);
          const float z=wall?std::min(16.f,std::max(0.f,(std::abs(y)-8)*2)):0;
          if(mirrored)y=-y;
          cloud->push_back({float(frame),y,z});
        }
        frames.push_back({begin,cloud->size()-begin,100+frame});
      }
      return weld::SeamDetector().detect(cloud,config,nullptr,nullptr,nullptr,&frames,false);
    };
    const auto result=detect(false,false);
    std::cout<<result.summary<<'\n';
    require(result.ok() && !result.leftFrameIndices.empty() && !result.rightFrameIndices.empty(),"Independent asymmetric fixture was not detected");
    std::cout<<"ranges "<<result.leftFrameIndices.front()<<' '<<result.leftFrameIndices.back()<<" / "<<result.rightFrameIndices.front()<<' '<<result.rightFrameIndices.back()<<'\n';
    require(result.leftFrameIndices.back()>=175 && result.rightFrameIndices.back()<=155,"Missing right wall shortened the left track or fabricated a right tail");
    require(result.rightFrameIndices.front()<=105 && result.leftFrameIndices.front()>=105,"Missing left wall shortened the right prefix");
    for(int side=0;side<2;++side) {
      const auto &cloud=side?*result.rightSeam:*result.leftSeam;
      const auto &ids=side?result.rightFrameIndices:result.leftFrameIndices;
      const auto &quality=side?result.rightPointQuality:result.leftPointQuality;
      require(ids.size()==cloud.size() && quality.size()==cloud.size(),"Unequal tracks lost metadata alignment");
      for(std::size_t i=0;i<cloud.size();++i) {
        require(std::fabs(cloud[i].x-(ids[i]-100))<1.e-4f,"Independent frame restoration used the other side's frame IDs");
        require(std::fabs(cloud[i].y-(side?8:-8))<2,"Independent candidate left the measured wall foot");
      }
    }
    const auto single=detect(true,false),mirrored=detect(true,true);
    require(single.ok() && single.rightSeam->empty() && single.leftSeam->size()>50,"Single left wall required an opposite trajectory");
    require(mirrored.ok() && mirrored.leftSeam->empty() && mirrored.rightSeam->size()>50,"Single right wall required an opposite trajectory");
    for(const auto *oneSide:{&single,&mirrored}) {
      require(oneSide->weldSurface->size()>10000,"Single-side detection discarded the measured floor");
      require((oneSide==&single?oneSide->rightSurface:oneSide->leftSurface)->empty(),"Floor recovery invented a missing wall surface");
      bool reachedOpenEdge=false;
      for(const auto &p:*oneSide->weldSurface) {
        require(std::abs(p.z)<1.e-4f,"Single-side floor contains upper wall points");
        reachedOpenEdge|=(oneSide==&single?p.y:-p.y)>30;
      }
      require(reachedOpenEdge,"Single-side floor was clipped to a fabricated opposite boundary");
    }
    // Endpoint tests include deliberately unequal lengths and original frame IDs.
    weld::SeamResult trimmed;
    for(int i=0;i<=400;++i) {
      const float x=.5f*i,bend=.03f*std::pow(std::max(0.f,x-175),2);
      trimmed.leftSeam->push_back({x,-10+bend,0});trimmed.leftFrameIndices.push_back(1000+i);
      trimmed.leftPointConfidence.push_back(.9f);trimmed.leftPointQuality.emplace_back();
      if(i<350) {
        trimmed.rightSeam->push_back({x,10,0});trimmed.rightFrameIndices.push_back(1000+i);
        trimmed.rightPointConfidence.push_back(.8f);trimmed.rightPointQuality.emplace_back();
      }
    }
    const auto originalRight=*trimmed.rightSeam;
    config.trimUnstableEndpoints=true;
    weld::trimSeamEndpoints(trimmed,config);
    require(trimmed.sideEndpointTrimmedEnd[0]>0 && trimmed.sideEndpointTrimmedEnd[1]==0,"A left hook trimmed the straight right side");
    require(trimmed.rightSeam->size()==originalRight.size() && trimmed.rightFrameIndices.back()==1349,"Independent trim changed right IDs/count");
    for(std::size_t i=0;i<originalRight.size();++i)require((weld::toVector(originalRight[i])-weld::toVector((*trimmed.rightSeam)[i])).norm()==0,"Independent trim moved right points");
    std::cout<<"Independent trajectory, single wall, metadata and endpoint checks passed\n";
    return 0;
  } catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
}
