#include "HullSurfaces.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <thread>

namespace {
void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}
void same(const weld::HullSurfacesResult &a, const weld::HullSurfacesResult &b) {
  require(a.error == b.error && a.sourceIndices == b.sourceIndices &&
          a.surfaceIds == b.surfaceIds && a.regionSizes == b.regionSizes &&
          a.recoveredIndices == b.recoveredIndices &&
          a.protectedUncertainIndices == b.protectedUncertainIndices &&
          a.analysisPoints == b.analysisPoints &&
          a.referencePatches == b.referencePatches &&
          a.rejectedNearSurfacePoints == b.rejectedNearSurfacePoints,
          "Parallel hull classification differs from serial output");
  require(a.regionCenters.size() == b.regionCenters.size() &&
          a.regionNormals.size() == b.regionNormals.size() &&
          a.cloud->size() == b.cloud->size(), "Hull output dimensions differ");
  for (std::size_t i = 0; i < a.regionCenters.size(); ++i) {
    require(a.regionCenters[i] == b.regionCenters[i] &&
            a.regionNormals[i] == b.regionNormals[i], "Region geometry changed");
  }
  for (std::size_t i = 0; i < a.cloud->size(); ++i) {
    require(weld::toVector((*a.cloud)[i]) == weld::toVector((*b.cloud)[i]),
            "Parallel hull extraction changed point coordinates/order");
  }
}
}

int main() try {
  weld::CloudT::Ptr hull(new weld::CloudT);
  for (int xi = 0; xi <= 200; ++xi) {
    const float x = .5f * xi;
    for (int yi = -50; yi <= 50; ++yi) hull->push_back({x, .5f * yi, 0});
    for (int side : {-1, 1}) {
      for (int ai = 0; ai <= 31; ++ai) {
        const float angle = .05f * ai;
        hull->push_back({x, side * (25 + 5 * std::sin(angle)), 5 - 5 * std::cos(angle)});
      }
      for (int zi = 10; zi <= 60; ++zi) hull->push_back({x, side * 30.f, .5f * zi});
    }
  }
  const auto bodyCount = hull->size();
  for (int xi = 0; xi < 20; ++xi)
    for (int yi = 0; yi < 20; ++yi) hull->push_back({float(xi), float(yi), -20});
  hull->push_back({200, 200, 200});
  hull->push_back({std::numeric_limits<float>::quiet_NaN(), 0, 0});
  weld::DetectorConfig config;
  const auto caller = std::this_thread::get_id();
  std::size_t callbacks = 0;
  const auto serial = weld::extractHullSurfaces(hull, config, [&] {
    require(std::this_thread::get_id() == caller, "Cancellation callback moved threads");
    ++callbacks;
    return false;
  });
  require(serial.ok() && serial.analysisPoints >= 8192 && callbacks > 0,
          "Fixture does not exercise parallel-sized hull extraction");
  require(serial.cloud->size() > .95 * bodyCount, "Lost synthetic hull support");
  for (auto index : serial.sourceIndices)
    require(index < bodyCount, "Kept disconnected/nonfinite interference");
  for (int repeat = 0; repeat < 2; ++repeat)
    same(serial, weld::extractHullSurfaces(hull, config));
  bool cancelled = false;
  try { weld::extractHullSurfaces(hull, config, [] { return true; }); }
  catch (const std::runtime_error &e) {
    cancelled = std::string(e.what()) == "Hull extraction cancelled";
  }
  require(cancelled, "Cancellation no longer interrupts hull extraction");
  std::cout << "Hull parallel/serial equality and cancellation passed\n";
  return 0;
} catch (const std::exception &e) {
  std::cerr << e.what() << '\n';
  return 1;
}
