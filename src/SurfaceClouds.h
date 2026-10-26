#pragma once
#include "SeamTypes.h"

namespace weld {
// Internal geometric boundaries, captured before single/constant-width track
// conversion. XYZ coordinates are original measurements' coordinates; station
// is the detector's ordered working station (not necessarily original X).
struct SurfaceSection {
  float station = 0;
  int frameIndex = -1;
  PointT left, right;
  Eigen::Vector2f leftNormal = Eigen::Vector2f::Zero();
  Eigen::Vector2f rightNormal = Eigen::Vector2f::Zero();
};

// Select measured points only. A point belongs to at most one output. Missing
// frames/sections are not extrapolated; no synthetic points are generated.
void extractSurfaceClouds(const CloudT::ConstPtr &source,
                          const std::vector<SurfaceSection> &sections,
                          const DetectorConfig &config, SeamResult &result,
                          const std::vector<ScanFrameRange> *frames = nullptr,
                          bool generateSurfaceModels = true);
}
