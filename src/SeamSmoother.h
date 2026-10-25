#pragma once

#include "SeamTypes.h"

namespace weld {

// Applies the optional final robust constant-slope Kalman filter and
// fixed-interval RTS smoother.  Legacy mode is an intentional no-op.
void applySeamPostprocess(SeamResult &result, const DetectorConfig &config);

// Call after per-point geometry/candidate guards, before restoring scan planes.
// Fits a robust continuous trend within each segment and preserves station IDs.
void applyFinalCurveSmoothing(SeamResult &result, const DetectorConfig &config);
struct SeamWidthDecision {
  int mode = 2;
  std::size_t singleVotes = 0, doubleVotes = 0;
};
// Decide once from the complete supported trajectory, before endpoint cuts.
SeamWidthDecision classifySeamWidth(const SeamResult &result, const DetectorConfig &config);

// Trim the center in single mode, or either boundary in double mode.
// Both observed boundaries retain the same range for later output generation.
// Use monotonic scan stations as the
// parameter (world X may turn back). Without stations, use the configured axis.
void trimSeamEndpoints(SeamResult &result, const DetectorConfig &config,
                       const std::vector<float> *stations = nullptr,
                       const SeamWidthDecision *decision = nullptr);
// Choose one mode for the complete fitted trajectory by station-width majority.
// Equal votes select the center track; preserve existing physical segment gaps.
void applySeamWidthRule(SeamResult &result, const DetectorConfig &config,
                       const SeamWidthDecision *decision = nullptr);

} // namespace weld
