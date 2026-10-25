#pragma once

#include <Eigen/Core>
#include <vector>

namespace weld {

// Locate the midpoint of the measured low-height support band (+V opening;
// mirrored for -V). A short median suppresses isolated spikes. The reference
// is independent of the candidate pair and of asymmetric scan margins.
float profileBottomMidpoint(const std::vector<Eigen::Vector2f> &points,
                            int openingDirection, float heightTolerance);

// Strict side identity: a point on the dividing line is not a left/right
// junction. Never clamp an unsupported candidate onto that line.
bool profilePairStraddlesBottom(float leftU, float rightU, float midpointU);

// Mark short spikes returning to a locally supported steep wall. True corners,
// long features and unsupported ends remain unchanged; no replacement points.
std::vector<bool> wallProfileSpikeMask(const std::vector<Eigen::Vector2f> &points,
                                      float radius, float maximumWidth, float tolerance);

struct ProfileCurvatureSample {
  Eigen::Vector2f point = Eigen::Vector2f::Zero();
  Eigen::Vector2f tangent = Eigen::Vector2f::Zero();
  float curvature = 0; // Signed d(theta)/ds, in inverse point-cloud units.
  float rms = 0;
  bool valid = false;
};

// Overlapping cubic B-spline spans (clamped single-span/Bernstein basis),
// fitted in chord-length coordinates. Derivatives are evaluated analytically.
// A physical radius smooths noise without dividing by transverse dx near a
// vertical wall. Discontinuities must be supplied as separate segments.
std::vector<ProfileCurvatureSample> fitProfileCurvature(
    const std::vector<Eigen::Vector2f> &points, float radius);

// One continuous, regularized parametric cubic B-spline for the complete
// connected profile. Samples include analytic first/second-derivative curvature.
std::vector<ProfileCurvatureSample> fitParametricProfileSpline(
    const std::vector<Eigen::Vector2f> &points, float knotSpacing, float sampleStep);

struct ProfileTurnSample {
  unsigned candidateKinds = 0; // 1: turn, 2: local valley, 4: curvature change
  Eigen::Vector2f point = Eigen::Vector2f::Zero();
  Eigen::Vector2f before = Eigen::Vector2f::Zero(), after = Eigen::Vector2f::Zero();
  Eigen::Vector2f beforeCenter = Eigen::Vector2f::Zero(), afterCenter = Eigen::Vector2f::Zero();
  float beforeRms = 0, afterRms = 0;
  float angle = 0, curvature = 0, rms = 0;
  int index = -1;
};

// Physical arc-length windows; invalid samples delimit components. Return a
// mask so coincident feature types merge before one-sided support validation.
std::vector<unsigned> supplementalProfileCandidates(
    const std::vector<ProfileCurvatureSample> &samples, float windowLength,
    float prominence, int openingDirection);

struct ProfileTurnRegion {
  int begin = 0, end = 0, peak = 0, representative = 0; // inclusive sample indices
  float arcLength = 0, turningAngle = 0;
};

// Integrate positive turning mass in a measured wall's local support band.
// Return an actual supported spline sample nearest its arc-length centroid;
// two nearby curvature lobes contribute continuously instead of competing.
int wallSupportedTurnCentroid(const std::vector<ProfileCurvatureSample> &samples,
    const Eigen::Vector2f &center, const Eigen::Vector2f &normal, float offset,
    float tolerance, float radius, int opening, float midpoint, int side);

// Keep a turn's identity when its nominal median just misses wall support:
// return the nearest admissible sample in THIS lobe, or -1. Never project a
// point off the spline or borrow from another curvature lobe.
int supportedProfileTurnRepresentative(
    const std::vector<ProfileCurvatureSample> &samples, const ProfileTurnRegion &region,
    const Eigen::Vector2f &normal, float offset, const Eigen::Vector2f &foot,
    float wallTolerance, float footRadius);

// Find connected same-sign high-curvature lobes and return an actual spline
// sample at the arc-length median or the median of integrated |curvature| ds.
// Regions cannot cross invalid samples, low-curvature valleys or maxArcLength.
std::vector<ProfileTurnRegion> locateProfileTurnRegions(
    const std::vector<ProfileCurvatureSample> &samples, bool angleWeighted,
    float relativeThreshold, float maxArcLength);

// Return the inclusive region endpoint nearest the groove floor. For a trough
// opening toward +V the floor is at the lower V endpoint; the rule reverses
// for -V. With either direction, infer the opening from the complete profile.
int profileTurnFloorBoundary(
    const std::vector<ProfileCurvatureSample> &samples,
    const ProfileTurnRegion &region, int openingDirection);

// Estimate one-sided tangent directions on a sampled curve at equal arc length.
// The spline detector uses these only to validate wall orientation and support;
// it replaces the finite-window turn rate with analytic spline curvature.
// An optional per-vertex mask skips tangent fits outside eligible indices,
// preserving probe coordinates, index rounding and all retained values.
std::vector<ProfileTurnSample> measureProfileTurns(
    const std::vector<Eigen::Vector2f> &points, float span, float step,
    const std::vector<bool> *eligibleIndices = nullptr);

} // namespace weld
