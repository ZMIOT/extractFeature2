#include "SeamDetector.h"
#include "AlgorithmParameters.h"
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>

namespace {
void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}
void synthetic(bool bead) {
  weld::CloudT::Ptr cloud(new weld::CloudT);
  std::vector<weld::ScanFrameRange> frames;
  for (int f=0; f<32; ++f) {
    const auto begin=cloud->size();
    for (int j=-300; j<=300; ++j) {
      const float y=.1f*j;
      const float z=std::min(24.f,std::max(0.f,2*(std::fabs(y)-12))) +
          (bead ? 2.6f*std::exp(-.5f*std::pow(y/2,2)) : 0);
      // A tilted acquisition plane: density changes can otherwise move the
      // median longitudinal coordinate and corrupt the station pitch.
      cloud->push_back({.6f*f+.25f*y+.1f*z,y,z});
      if (j < -50 && (j+f)%3==0) cloud->push_back(cloud->back());
    }
    frames.push_back({begin,cloud->size()-begin,100+f});
  }
  auto config=weld::AlgorithmParameters().get(weld::ProfileLocalizationMethod::CurvatureCandidates);
  config.extractMainBody=false;config.voxelSize=0;config.outlierMeanK=0;
  config.trimUnstableEndpoints=false;config.finalSmoothingEnabled=false;
  std::map<int,float> stations;
  weld::SeamResult reference;
  for (float size : {0.f,.2f,.5f,.8f}) {
    config.inputVoxelSize=size;
    weld::DebugTrace trace;
    const auto result=weld::SeamDetector().detect(cloud,config,nullptr,nullptr,&trace,&frames,false);
    require(result.ok(),"Sampling rejected supported synthetic geometry");
    if(size==0) reference=result;
    else {
      require(result.leftFrameIndices==reference.leftFrameIndices && result.rightFrameIndices==reference.rightFrameIndices,
              "Sampling changed supported frame coverage");
      for(int side=0;side<2;++side) {
        const auto &a=side?result.rightSeam:result.leftSeam;
        const auto &b=side?reference.rightSeam:reference.leftSeam;
        for(std::size_t i=0;i<a->size();++i)
          require(weld::toVector((*a)[i])==weld::toVector((*b)[i]),
                  "Sampling changed full-resolution framed feature localization");
      }
    }
    int selected=0,beadSelections=0;
    for (const auto &slice:trace.slices) {
      if (size==0) stations[slice.frameIndex]=slice.station;
      else require(stations.at(slice.frameIndex)==slice.station,
                   "Sampling changed the original frame station clock");
      if (slice.selected>=0) {
        ++selected;
        beadSelections+=slice.pairs[slice.selected].source=="central_bead_valleys";
      }
    }
    std::cout<<"bead="<<bead<<" voxel="<<size<<" selected="<<selected
             <<" valleys="<<beadSelections<<'\n';
    require(selected>=28,"Sampling lost supported synthetic frames");
    require(bead ? beadSelections>=28 : beadSelections==0,
            "Sampling stability fix suppressed real beads or invented a bead");
    for (int side=0; side<2; ++side) {
      const auto &points=side?result.rightSeam:result.leftSeam;
      const auto &ids=side?result.rightFrameIndices:result.leftFrameIndices;
      for (std::size_t i=0;i<points->size();++i) {
        const auto &p=(*points)[i];
        require(std::fabs(p.x-.25f*p.y-.1f*p.z-.6f*(ids[i]-100))<.001f,
                "Sampling rotated or moved the acquisition plane");
        if (bead) require(std::fabs(p.y)>1 && std::fabs(p.y)<9,
                          "True central bead was replaced by outer wall feet");
      }
    }
  }
}
}
int main() try {
  synthetic(false);synthetic(true);
  return 0;
} catch (const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
