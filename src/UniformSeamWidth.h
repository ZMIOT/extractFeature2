#pragma once
#include "SeamTypes.h"

namespace weld {
// Contour and candidate points are in the original scan plane, including its
// longitudinal tilt. Missing/unsupported sections have an empty contour.
struct UniformWidthSection {
  double station = 0;
  std::vector<PointT> contour;
  float bottomCenterU = 0;
  PointT sourceLeft, sourceRight;
};
void applyUniformSeamWidth(SeamResult &result, const DetectorConfig &config,
                           const std::vector<UniformWidthSection> &sections);
}
