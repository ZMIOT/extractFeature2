#pragma once
#include "SurfaceClouds.h"
namespace weld {
// Uses measured section identities and labels, before those are discarded.
BSplineSurfaceModel fitSurfaceModel(const CloudT &source,
    const std::vector<std::vector<std::size_t>> &groups,
    const std::vector<SurfaceSection> &sections,
    const std::vector<unsigned char> &labels, unsigned char face,
    const DetectorConfig &config);
}
