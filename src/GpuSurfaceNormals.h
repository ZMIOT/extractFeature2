#pragma once

#include "SeamTypes.h"

#include <Eigen/Core>
#include <string>
#include <vector>

namespace weld {
namespace detail {

struct GpuSurfaceNormals {
  std::vector<std::vector<Eigen::Vector2f>> normals;
  std::vector<std::vector<unsigned char>> supported;
  std::string device;
  double milliseconds = 0.0;
};

// Optional OpenCL acceleration for the independent radius-neighbour covariance
// calculations used by surface-cloud extraction. Failure is a normal CPU
// fallback, so the public SDK has no mandatory GPU runtime dependency.
bool computeSurfaceNormalsGpu(
    const CloudT &source, const std::vector<std::vector<std::size_t>> &groups,
    int uAxis, int vAxis, float radius, GpuSurfaceNormals &result);

} // namespace detail
} // namespace weld
