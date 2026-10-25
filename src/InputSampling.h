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
}
