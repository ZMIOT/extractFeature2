#include "InputSampling.h"
#include "SeamDetector.h"
#include "AlgorithmParameters.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool condition,const char* message) {
  if(!condition)throw std::runtime_error(message);
}
bool same(const weld::CloudT &a,const weld::CloudT &b) {
  if(a.size()!=b.size())return false;
  for(std::size_t i=0;i<a.size();++i)
    if(a[i].x!=b[i].x || a[i].y!=b[i].y || a[i].z!=b[i].z)return false;
  return true;
}
void samplingRules() {
  weld::CloudT::Ptr input(new weld::CloudT);
  for(int f=0;f<2;++f) {
    input->push_back({.01f,.01f,.01f});
    input->push_back({.04f,.04f,.04f});
    input->push_back({.05f,.05f,.05f});
    input->push_back({.09f,.09f,.09f});
  }
  const weld::CloudT original=*input;
  const std::vector<weld::ScanFrameRange> frames{{0,4,10},{4,4,20}};
  const auto result=weld::sampleMeasuredInput(input,.1f,&frames);
  require(result.error.empty() && result.cloud->size()==6,"Unexpected representative count");
  require(result.frames.size()==2 && result.frames[0].frameIndex==10 &&
          result.frames[1].frameIndex==20 && result.frames[1].begin==3 &&
          result.frames[0].count==3 && result.frames[1].count==3,"Frames were mixed or renumbered");
  for(int f=0;f<2;++f) {
    require((*result.cloud)[3*f].x==.01f && (*result.cloud)[3*f+1].x==.05f &&
            (*result.cloud)[3*f+2].x==.09f,"Endpoints, nearest measured point or input order changed");
  }
  require(same(*input,original),"Sampling modified the input");
  require(weld::sampleMeasuredInput(input,0,&frames).cloud==input,"Disabled sampling copied the cloud");
  auto invalid=frames;invalid[1].begin=3;
  require(!weld::sampleMeasuredInput(input,.1f,&invalid).error.empty(),"Malformed frame map accepted");
  require(!weld::sampleMeasuredInput(input,-1).error.empty(),"Negative sampling size accepted");
  require(!weld::sampleMeasuredInput(input,std::numeric_limits<float>::infinity()).error.empty(),
          "Nonfinite sampling size accepted");
  weld::CloudT::Ptr extreme(new weld::CloudT);
  extreme->push_back({-1.e30f,0,0});extreme->push_back({1.e30f,0,0});
  const auto large=weld::sampleMeasuredInput(extreme,1.e-6f);
  require(large.error.empty() && same(*large.cloud,*extreme),"Extreme coordinates overflowed voxel keys");
}
void commonPipeline(bool framed,bool mainBody) {
  weld::CloudT::Ptr input(new weld::CloudT);
  std::vector<weld::ScanFrameRange> frames;
  for(int f=0;f<80;++f) {
    const auto begin=input->size();
    for(int j=-500;j<=500;++j) {
      const float y=.05f*j;
      input->push_back({.8f*f,y,std::min(15.f,std::max(0.f,2*(std::abs(y)-8)))});
    }
    frames.push_back({begin,input->size()-begin,100+f});
  }
  const auto sampled=weld::sampleMeasuredInput(input,.2f,framed?&frames:nullptr);
  require(sampled.error.empty() && sampled.cloud->size()<input->size(),"Fixture was not sampled");
  auto config=weld::AlgorithmParameters().get(weld::ProfileLocalizationMethod::CurvatureCandidates);
  config.extractMainBody=mainBody;config.outlierMeanK=0;
  config.voxelSize=.8f;config.inputVoxelSize=.2f;
  const auto actual=weld::SeamDetector().detect(input,config,nullptr,nullptr,nullptr,framed?&frames:nullptr);
  config.inputVoxelSize=0;config.voxelSize=0;
  const auto expected=weld::SeamDetector().detect(sampled.cloud,config,nullptr,nullptr,nullptr,
                                                framed?&sampled.frames:nullptr);
  require(actual.ok() && expected.ok(),"Sampled groove detection failed");
  require(actual.inputSamplingPoints==sampled.cloud->size() &&
          actual.inputSamplingRemovedPoints==input->size()-sampled.cloud->size(),"Sampling counts incorrect");
  require(same(*actual.leftSeam,*expected.leftSeam) && same(*actual.rightSeam,*expected.rightSeam) &&
          same(*actual.leftSurface,*expected.leftSurface) && same(*actual.rightSurface,*expected.rightSurface) &&
          same(*actual.weldSurface,*expected.weldSurface) && same(*actual.mainBodyCloud,*expected.mainBodyCloud),
          "A downstream stage used the original cloud or applied another voxel pass");
  require(actual.leftFrameIndices==expected.leftFrameIndices && actual.rightFrameIndices==expected.rightFrameIndices,
          "Sampling changed frame identity inside the pipeline");
}
}
int main() {
  try {
    samplingRules();
    for(bool framed:{false,true})for(bool mainBody:{false,true})commonPipeline(framed,mainBody);
    std::cout<<"Input sampling: measured subsets, frame identity and shared downstream cloud passed\n";
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
