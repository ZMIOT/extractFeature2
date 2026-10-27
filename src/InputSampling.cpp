#include "InputSampling.h"
#include <pcl/common/point_tests.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <unordered_set>

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
bool cellOf(const PointT &p, float size, Cell &cell) {
  for (int axis=0; axis<3; ++axis) {
    const double index=std::floor(double(p.data[axis])/double(size));
    if (!std::isfinite(index) || index < -9223372036854775808.0 || index >= 9223372036854775808.0)
      return false;
    cell[axis]=static_cast<std::int64_t>(index);
  }
  return true;
}
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

InputSamplingResult recoverMeasuredFrameSupport(
    const CloudT::ConstPtr &input, const std::vector<ScanFrameRange> &frames,
    const CloudT::ConstPtr &support, const std::vector<ScanFrameRange> &supportFrames,
    float size) {
  InputSamplingResult result;
  if (!input || !support || frames.empty() || supportFrames.empty() || !std::isfinite(size) || size<=0) {
    result.error="Invalid measured frame support input."; return result;
  }
  // Reuse the zero-size path's complete range validation without copying points.
  if (!sampleMeasuredInput(input,0,&frames).error.empty() ||
      !sampleMeasuredInput(support,0,&supportFrames).error.empty()) {
    result.error="Invalid measured frame support ranges."; return result;
  }
  CloudT::Ptr recovered(new CloudT);
  recovered->reserve(input->size());
  std::unordered_set<Cell,CellHash> cells;
  for (const auto &kept:supportFrames) {
    const auto original=std::lower_bound(frames.begin(),frames.end(),kept.frameIndex,
        [](const ScanFrameRange &r,int id){return r.frameIndex<id;});
    if (original==frames.end() || original->frameIndex!=kept.frameIndex) {
      result.error="Body frame has no original support."; return result;
    }
    cells.clear();cells.reserve(kept.count);
    std::vector<std::array<float,3>> extreme;
    for (std::size_t i=kept.begin;i<kept.begin+kept.count;++i) {
      Cell cell;
      const auto &p=(*support)[i];
      if (cellOf(p,size,cell)) cells.insert(cell);
      else extreme.push_back({p.x,p.y,p.z});
    }
    const auto begin=recovered->size();
    for (std::size_t i=original->begin;i<original->begin+original->count;++i) {
      const auto &p=(*input)[i];Cell cell;
      if (!pcl::isFinite(p)) continue;
      const bool retained=cellOf(p,size,cell) ? cells.count(cell)!=0 :
          std::find(extreme.begin(),extreme.end(),std::array<float,3>{{p.x,p.y,p.z}})!=extreme.end();
      if (retained) recovered->push_back(p);
    }
    if (recovered->size()>begin)
      result.frames.push_back({begin,recovered->size()-begin,kept.frameIndex});
  }
  recovered->header=input->header;
  recovered->sensor_origin_=input->sensor_origin_;
  recovered->sensor_orientation_=input->sensor_orientation_;
  recovered->is_dense=true;
  result.cloud=recovered;
  return result;
}
}
