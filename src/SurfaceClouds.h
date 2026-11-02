#pragma once
#include "SeamTypes.h"
#include <array>

namespace weld {
// Internal geometric boundaries, captured before single/constant-width track
// conversion. XYZ coordinates are original measurements' coordinates; station
// is the detector's ordered working station (not necessarily original X).
struct SurfaceSection {
  float station = 0;
  int frameIndex = -1;
  PointT left, right;
  bool leftValid = true, rightValid = true;
  // Trajectory endpoints may be rejected while the measured face continues.
  // Re-establish those boundaries from this section's broad wall and floor.
  bool leftBoundaryUncertain = false, rightBoundaryUncertain = false;
  Eigen::Vector2f leftNormal = Eigen::Vector2f::Zero();
  Eigen::Vector2f rightNormal = Eigen::Vector2f::Zero();
};

struct SurfaceEnd {
  bool valid = false;
  float station = 0;
  PointT point;
};

struct SurfaceFootRecovery {
  std::vector<std::size_t> kept, recovered;
};
// Validate a measured continuation of an apparently rejected trajectory tail.
// Empty output means there is no evidence of candidate drift; keep normal trim.
SurfaceFootRecovery selectMeasuredFootContinuation(const CloudT &track,
    const std::vector<float> &stations, const std::vector<SurfaceEnd> &feet,
    std::size_t firstRejected, const DetectorConfig &config);

// Select measured points only. A point belongs to at most one output. Missing
// frames/sections are not extrapolated; no synthetic points are generated.
void extractSurfaceClouds(const CloudT::ConstPtr &source,
                          const std::vector<SurfaceSection> &sections,
                          const DetectorConfig &config, SeamResult &result,
                          const std::vector<ScanFrameRange> *frames = nullptr,
                          bool generateSurfaceModels = true,
                          const std::array<SurfaceEnd,2> *ends = nullptr,
                          std::vector<std::array<SurfaceEnd,2>> *measuredToes = nullptr);
}
