#pragma once
#include "SeamDebug.h"

namespace weld {
struct HullSurfacesResult {
  CloudT::Ptr cloud{new CloudT};
  std::vector<std::size_t> sourceIndices;
  std::vector<std::uint32_t> surfaceIds;
  std::vector<std::size_t> regionSizes;
  std::vector<Eigen::Vector3f> regionCenters, regionNormals;
  std::size_t analysisPoints = 0;
  std::size_t referencePatches = 0;
  // Input indices recovered by bounded multi-patch support, not generated points.
  std::vector<std::size_t> recoveredIndices;
  // Existing raw points retained by the wider bounded cleanup ROI, without
  // claiming that they agree with a fitted surface or that they are noise-free.
  std::vector<std::size_t> protectedUncertainIndices;
  std::size_t rejectedNearSurfacePoints = 0;
  std::string error;
  bool ok() const { return error.empty() && !sourceIndices.empty(); }
};
HullSurfacesResult extractHullSurfaces(const CloudT::ConstPtr &input,
                                     const DetectorConfig &config,
                                     const std::function<bool()> &cancelled = {});
}
