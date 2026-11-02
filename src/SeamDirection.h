#pragma once
#include "SeamTypes.h"

namespace weld {
// Analyze retained, measured body points without changing their coordinates.
struct SeamDirectionEstimate {
  Axis axis = Axis::X;
  Eigen::Vector3f principalDirection = Eigen::Vector3f::Zero();
  float principalGap = 0;
  bool axisReliable = false;
  int opening = 0, positiveVotes = 0, negativeVotes = 0;
  bool openingReliable = false;
};
SeamDirectionEstimate estimateSeamDirection(const CloudT &cloud,
    const DetectorConfig &config, const std::vector<ScanFrameRange> *frames = nullptr);
}
