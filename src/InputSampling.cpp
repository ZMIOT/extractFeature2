#include "InputSampling.h"
#include <pcl/common/point_tests.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>

namespace weld {
namespace {
using Cell=std::array<std::int64_t,3>;
struct CellHash {
  std::size_t operator()(const Cell &cell) const {
    std::size_t hash=0;
    for(auto value:cell)
      hash^=std::hash<std::int64_t>{}(value)+std::size_t(0x9e3779b9u)+(hash<<6)+(hash>>2);
    return hash;
  }
};
struct Choice {std::size_t index;double distance;};
}

InputSamplingResult sampleMeasuredInput(const CloudT::ConstPtr &input, float size,
                                       const std::vector<ScanFrameRange> *frames) {
  InputSamplingResult result;
  if(!input || !std::isfinite(size) || size<0) {
    result.error="Invalid input sampling cloud or voxel size.";return result;
  }
  const bool framed=frames && !frames->empty();
  if(framed) {
    std::size_t expected=0;int previous=-1;
    for(const auto &frame:*frames) {
      if(frame.begin!=expected || !frame.count || expected>input->size() ||
         frame.count>input->size()-expected || frame.frameIndex<=previous) {
        result.error="Invalid frame ranges for input sampling.";return result;
      }
      expected+=frame.count;previous=frame.frameIndex;
    }
    if(expected!=input->size()) {
      result.error="Input sampling frames do not cover the cloud.";return result;
    }
  }
  if(size==0) {
    result.cloud=input;if(framed)result.frames=*frames;return result;
  }
  CloudT::Ptr sampled(new CloudT);
  sampled->reserve(input->size());
  std::unordered_map<Cell,Choice,CellHash> cells;
  std::vector<unsigned char> keep;
  auto sampleRange=[&](std::size_t begin,std::size_t count) {
    cells.clear();cells.reserve(count);keep.assign(count,0);
    if(count)keep.front()=keep.back()=1;
    for(std::size_t offset=0;offset<count;++offset) {
      const auto &point=(*input)[begin+offset];
      if(!pcl::isFinite(point)) {keep[offset]=1;continue;}
      Cell cell;double squaredDistance=0;bool representable=true;
      for(int axis=0;axis<3;++axis) {
        const double coordinate=point.data[axis];
        const double index=std::floor(coordinate/double(size));
        // Pass through extreme coordinates instead of overflowing integer keys.
        if(index < -9223372036854775808.0 || index >= 9223372036854775808.0) {
          representable=false;break;
        }
        cell[axis]=static_cast<std::int64_t>(index);
        const double delta=coordinate-(index+.5)*double(size);
        squaredDistance+=delta*delta;
      }
      if(!representable) {keep[offset]=1;continue;}
      const auto inserted=cells.emplace(cell,Choice{offset,squaredDistance});
      if(!inserted.second && squaredDistance<inserted.first->second.distance)
        inserted.first->second={offset,squaredDistance};
    }
    for(const auto &entry:cells)keep[entry.second.index]=1;
    for(std::size_t offset=0;offset<count;++offset)
      if(keep[offset])sampled->push_back((*input)[begin+offset]);
  };
  if(framed) {
    result.frames.reserve(frames->size());
    for(const auto &frame:*frames) {
      const auto begin=sampled->size();
      sampleRange(frame.begin,frame.count);
      result.frames.push_back({begin,sampled->size()-begin,frame.frameIndex});
    }
  } else sampleRange(0,input->size());
  sampled->header=input->header;
  sampled->sensor_origin_=input->sensor_origin_;
  sampled->sensor_orientation_=input->sensor_orientation_;
  sampled->is_dense=input->is_dense;
  result.cloud=sampled;
  return result;
}
}
