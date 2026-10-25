#pragma once
#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace weld {
/** Non-rational tensor-product B-spline, in original XYZ coordinates/units.
 * S(u,v) = sum N_i,degreeU(u) N_j,degreeV(v) controlPoints[i*countV+j].
 * Clamped full knot vectors (including all endpoint repetitions), u,v in [0,1].
 * u follows scan order; v follows the measured cross-section contour.
 * Each patch has independent parameters; never connect separate patches.
 */
struct BSplineSurfacePatch {
  int degreeU = 3, degreeV = 3;
  int countU = 0, countV = 0;
  std::vector<double> knotsU, knotsV;
  std::vector<std::array<double,3>> controlPoints;
  int firstFrame = -1, lastFrame = -1;
  double firstStation = 0, lastStation = 0;
  std::size_t sampleCount = 0;
  // Residual at assigned scan/contour parameters, NOT closest-point distance.
  double rmsError = 0, maxError = 0;
  bool valid() const;
  // False on invalid model, non-finite/out-of-range parameters or singular normal.
  // Output arguments are left unchanged on failure. No extrapolation.
  bool evaluate(double u, double v, Eigen::Vector3d &point,
                Eigen::Vector3d *normal = nullptr) const;
};

/** One named face; gaps or missing face support split it into separate patches.
 * A valid model may cover only part of the measured face. rejectedPatchCount
 * counts runs lacking sufficient 2-D support or failing fitting/quality checks.
 */
struct BSplineSurfaceModel {
  std::vector<BSplineSurfacePatch> patches;
  std::size_t rejectedPatchCount = 0;
  std::string message = "not_fitted";
  bool valid() const { return !patches.empty(); }
  bool evaluate(std::size_t patch, double u, double v, Eigen::Vector3d &point,
                Eigen::Vector3d *normal = nullptr) const;
};
}
