#pragma once
#include "SeamTypes.h"

namespace weld {
struct InputSamplingResult {
  CloudT::ConstPtr cloud;
  std::vector<ScanFrameRange> frames;
  std::string error;
};

// Keep a measured point nearest each voxel center, in input order. Frames are
// sampled independently; their identifiers and both endpoint samples survive.
// A zero size returns the original cloud. No averaged coordinates are generated.
InputSamplingResult sampleMeasuredInput(const CloudT::ConstPtr &input, float size,
                                       const std::vector<ScanFrameRange> *frames = nullptr);

// Restore measured detail only inside cells retained by the sampled body mask.
// This is a separate trajectory-support cloud, not the sampled output cloud.
InputSamplingResult recoverMeasuredFrameSupport(
    const CloudT::ConstPtr &input, const std::vector<ScanFrameRange> &frames,
    const CloudT::ConstPtr &support, const std::vector<ScanFrameRange> &supportFrames,
    float size);
}
