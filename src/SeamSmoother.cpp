#include "SeamSmoother.h"

#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <limits>
#include <stdexcept>

namespace weld {
namespace {

using State = Eigen::Vector2f;
using Covariance = Eigen::Matrix2f;

float clamp01(float value) {
  return std::max(0.0f, std::min(1.0f, value));
}

std::array<int, 2> localTransverseAxes(Axis axis) {
  switch (axis) {
  case Axis::X:
    return {{1, 2}};
  case Axis::Y:
    return {{0, 2}};
  case Axis::Z:
    return {{0, 1}};
  }
  return {{1, 2}};
}

std::vector<std::size_t> normalizedStarts(const std::vector<std::size_t> &starts,
                                          std::size_t count) {
  std::vector<std::size_t> result;
  if (count == 0)
    return result;
  result.push_back(0);
  for (std::size_t start : starts)
    if (start > 0 && start < count && start > result.back())
      result.push_back(start);
  return result;
}

float observationVariance(const SeamPointQuality &quality,
                          const DetectorConfig &config) {
  const float confidence = std::max(0.20f, quality.sourceConfidence);
  float deviation = config.rtsMeasurementStddev / confidence;
  if (quality.flags & QualityInterpolated)
    deviation *= 3.0f;
  if (!quality.wallAccepted)
    deviation *= 1.35f;
  if (config.useSurfaceRefinement && config.useLocalFloorConstraint && !quality.floorAccepted)
    deviation *= 1.20f;
  return std::max(1.0e-6f, deviation * deviation);
}

void smoothCoordinate(const CloudT &points,
                      const std::vector<SeamPointQuality> &quality,
                      const std::vector<std::size_t> &segmentStarts, int axis,
                      int stationAxis, const DetectorConfig &config,
                      std::vector<float> &smoothed, std::vector<float> &nis,
                      std::vector<float> &sigma) {
  const std::size_t count = points.size();
  smoothed.resize(count);
  nis.assign(count, 0.0f);
  sigma.assign(count, 0.0f);
  for (std::size_t i = 0; i < count; ++i)
    smoothed[i] = points[i].data[axis];

  const auto starts = normalizedStarts(segmentStarts, count);
  for (std::size_t segment = 0; segment < starts.size(); ++segment) {
    const std::size_t begin = starts[segment];
    const std::size_t end = segment + 1 < starts.size() ? starts[segment + 1]
                                                        : count;
    if (end - begin < 3)
      continue;

    const std::size_t length = end - begin;
    std::vector<State> predicted(length), filtered(length), result(length);
    std::vector<Covariance> predictedCovariance(length),
        filteredCovariance(length), resultCovariance(length), transition(length);

    const float firstGap = std::max(
        1.0e-4f, points[begin + 1].data[stationAxis] -
                       points[begin].data[stationAxis]);
    filtered[0] << points[begin].data[axis],
        (points[begin + 1].data[axis] - points[begin].data[axis]) / firstGap;
    filteredCovariance[0] = Covariance::Zero();
    filteredCovariance[0](0, 0) = observationVariance(quality[begin], config);
    filteredCovariance[0](1, 1) =
        4.0f * filteredCovariance[0](0, 0) / (firstGap * firstGap);
    predicted[0] = filtered[0];
    predictedCovariance[0] = filteredCovariance[0];
    transition[0] = Covariance::Identity();

    for (std::size_t local = 1; local < length; ++local) {
      const std::size_t i = begin + local;
      const float ds = std::max(1.0e-4f, points[i].data[stationAxis] -
                                             points[i - 1].data[stationAxis]);
      Covariance f;
      f << 1.0f, ds, 0.0f, 1.0f;
      transition[local] = f;
      const float ds2 = ds * ds;
      Covariance q;
      q << 0.25f * ds2 * ds2, 0.5f * ds2 * ds,
          0.5f * ds2 * ds, ds2;
      q *= std::max(1.0e-8f, config.rtsProcessNoise);

      predicted[local] = f * filtered[local - 1];
      predictedCovariance[local] =
          f * filteredCovariance[local - 1] * f.transpose() + q;

      float variance = observationVariance(quality[i], config);
      const float innovation = points[i].data[axis] - predicted[local][0];
      float innovationVariance = predictedCovariance[local](0, 0) + variance;
      const float rawNis = innovation * innovation /
                           std::max(1.0e-8f, innovationVariance);
      nis[i] = rawNis;
      if (rawNis > config.rtsInnovationGate) {
        const float inflation = std::min(
            100.0f, rawNis / std::max(1.0e-4f, config.rtsInnovationGate));
        variance *= inflation;
        innovationVariance = predictedCovariance[local](0, 0) + variance;
      }

      const State gain = predictedCovariance[local].col(0) /
                         std::max(1.0e-8f, innovationVariance);
      filtered[local] = predicted[local] + gain * innovation;
      Covariance correction = Covariance::Identity();
      correction.col(0) -= gain;
      filteredCovariance[local] =
          correction * predictedCovariance[local] * correction.transpose() +
          variance * gain * gain.transpose();
    }

    result.back() = filtered.back();
    resultCovariance.back() = filteredCovariance.back();
    for (std::size_t local = length - 1; local-- > 0;) {
      const Covariance inverse =
          predictedCovariance[local + 1].ldlt().solve(Covariance::Identity());
      const Covariance gain = filteredCovariance[local] *
                              transition[local + 1].transpose() * inverse;
      result[local] = filtered[local] +
                      gain * (result[local + 1] - predicted[local + 1]);
      resultCovariance[local] =
          filteredCovariance[local] +
          gain * (resultCovariance[local + 1] -
                  predictedCovariance[local + 1]) *
              gain.transpose();
    }
    for (std::size_t local = 0; local < length; ++local) {
      const std::size_t i = begin + local;
      if (result[local].allFinite() && resultCovariance[local].allFinite()) {
        smoothed[i] = result[local][0];
        sigma[i] = std::sqrt(std::max(0.0f, resultCovariance[local](0, 0)));
      }
    }
  }
}

bool satisfiesGeometry(const Eigen::Vector3f &point,
                       const SeamPointQuality &quality,
                       const DetectorConfig &config) {
  const auto axes = localTransverseAxes(config.longitudinalAxis);
  const float tolerance = std::max(1.0e-4f, config.rtsGeometryTolerance);
  if (quality.wallAccepted) {
    const Eigen::Vector2f uv(point[axes[0]], point[axes[1]]);
    if (std::fabs(quality.wallNormal.dot(uv) - quality.wallOffset) > tolerance)
      return false;
  }
  if (quality.floorAccepted &&
      std::fabs(quality.floorNormal.dot(point) - quality.floorOffset) >
          tolerance)
    return false;
  return true;
}

Eigen::Vector3f projectSingleConstraint(Eigen::Vector3f point,
                                        const SeamPointQuality &quality,
                                        const DetectorConfig &config) {
  const auto axes = localTransverseAxes(config.longitudinalAxis);
  if (quality.wallAccepted) {
    Eigen::Vector2f uv(point[axes[0]], point[axes[1]]);
    const float squaredNorm = quality.wallNormal.squaredNorm();
    if (squaredNorm > 1.0e-8f)
      uv -= quality.wallNormal *
            ((quality.wallNormal.dot(uv) - quality.wallOffset) / squaredNorm);
    point[axes[0]] = uv.x();
    point[axes[1]] = uv.y();
  } else if (quality.floorAccepted) {
    Eigen::Vector2f normal(quality.floorNormal[axes[0]],
                           quality.floorNormal[axes[1]]);
    const float squaredNorm = normal.squaredNorm();
    if (squaredNorm > 1.0e-8f) {
      const int stationAxis = static_cast<int>(config.longitudinalAxis);
      Eigen::Vector2f uv(point[axes[0]], point[axes[1]]);
      const float offset = quality.floorOffset -
                           quality.floorNormal[stationAxis] * point[stationAxis];
      uv -= normal * ((normal.dot(uv) - offset) / squaredNorm);
      point[axes[0]] = uv.x();
      point[axes[1]] = uv.y();
    }
  }
  return point;
}

void smoothSide(CloudT &points, std::vector<float> &confidence,
                std::vector<SeamPointQuality> &quality,
                const std::vector<std::size_t> &segmentStarts,
                const DetectorConfig &config) {
  if (points.size() < 3 || quality.size() != points.size())
    return;
  const int stationAxis = static_cast<int>(config.longitudinalAxis);
  const auto axes = localTransverseAxes(config.longitudinalAxis);
  std::array<std::vector<float>, 2> values, nis, sigma;
  for (int component = 0; component < 2; ++component)
    smoothCoordinate(points, quality, segmentStarts, axes[component], stationAxis,
                     config, values[component], nis[component], sigma[component]);

  for (std::size_t i = 0; i < points.size(); ++i) {
    const Eigen::Vector3f original = toVector(points[i]);
    Eigen::Vector3f target = original;
    target[axes[0]] = values[0][i];
    target[axes[1]] = values[1][i];
    quality[i].innovationNis = std::max(nis[0][i], nis[1][i]);
    quality[i].positionSigma = std::hypot(sigma[0][i], sigma[1][i]);
    if (quality[i].innovationNis > config.rtsInnovationGate)
      quality[i].flags |= QualityInnovationOutlier;

    Eigen::Vector3f delta = target - original;
    delta[stationAxis] = 0.0f;
    const float requested = delta.norm();
    if (requested > config.rtsMaxCorrection && requested > 1.0e-8f) {
      delta *= config.rtsMaxCorrection / requested;
      quality[i].flags |= QualityCorrectionLimited;
    }
    target = projectSingleConstraint(original + delta, quality[i], config);

    if (!satisfiesGeometry(target, quality[i], config)) {
      const Eigen::Vector3f projected = target;
      float scale = 1.0f;
      for (int pass = 0; pass < 14; ++pass) {
        target = projectSingleConstraint(
            original + scale * (projected - original), quality[i], config);
        if (satisfiesGeometry(target, quality[i], config))
          break;
        scale *= 0.5f;
      }
      quality[i].flags |= QualityGeometryLimited;
      if (!satisfiesGeometry(target, quality[i], config))
        target = original;
    }
    points[i] = toPoint(target);
    quality[i].rtsCorrection = (target - original).norm();

    float evidence = 1.0f;
    if (quality[i].wallRms >= 0.0f)
      evidence = std::min(evidence, std::exp(-quality[i].wallRms /
                                             std::max(0.1f, config.surfaceRmsScale)));
    if (quality[i].floorRms >= 0.0f)
      evidence = std::min(evidence, std::exp(-quality[i].floorRms /
                                             std::max(0.1f, config.floorPatchMaxRms)));
    if (quality[i].innovationNis > config.rtsInnovationGate)
      evidence *= std::sqrt(config.rtsInnovationGate /
                            quality[i].innovationNis);
    if (i < confidence.size())
      confidence[i] = clamp01(confidence[i] * (0.5f + 0.5f * evidence));
  }
}

void enforcePairSeparation(SeamResult &result, const CloudT &leftBefore,
                           const CloudT &rightBefore,
                           const DetectorConfig &config) {
  const int sideAxis = localTransverseAxes(config.longitudinalAxis)[0];
  const std::size_t count = std::min(result.leftSeam->size(),
                                     result.rightSeam->size());
  for (std::size_t i = 0; i < count; ++i) {
    const float oldWidth = rightBefore[i].data[sideAxis] -
                           leftBefore[i].data[sideAxis];
    const float newWidth = (*result.rightSeam)[i].data[sideAxis] -
                           (*result.leftSeam)[i].data[sideAxis];
    if (newWidth >= config.minSeamSeparation)
      continue;
    float scale = 0.0f;
    if (oldWidth > config.minSeamSeparation && oldWidth > newWidth)
      scale = clamp01((oldWidth - config.minSeamSeparation) /
                      (oldWidth - newWidth));
    (*result.leftSeam)[i].data[sideAxis] =
        leftBefore[i].data[sideAxis] +
        scale * ((*result.leftSeam)[i].data[sideAxis] -
                 leftBefore[i].data[sideAxis]);
    (*result.rightSeam)[i].data[sideAxis] =
        rightBefore[i].data[sideAxis] +
        scale * ((*result.rightSeam)[i].data[sideAxis] -
                 rightBefore[i].data[sideAxis]);
    result.leftPointQuality[i].flags |= QualityPairConstraintLimited;
    result.rightPointQuality[i].flags |= QualityPairConstraintLimited;
  }
}

SeamSmoothingSummary summarize(const CloudT &before, const CloudT &after,
                               std::vector<SeamPointQuality> &quality) {
  SeamSmoothingSummary summary;
  summary.applied = true;
  float total = 0.0f;
  const std::size_t count = std::min(before.size(), after.size());
  for (std::size_t i = 0; i < count; ++i) {
    quality[i].rtsCorrection = (toVector(after[i]) - toVector(before[i])).norm();
    total += quality[i].rtsCorrection;
    summary.maxCorrection = std::max(summary.maxCorrection,
                                     quality[i].rtsCorrection);
    summary.innovationOutliers +=
        (quality[i].flags & QualityInnovationOutlier) != 0;
    summary.correctionLimited +=
        (quality[i].flags & QualityCorrectionLimited) != 0;
    summary.geometryLimited +=
        (quality[i].flags & QualityGeometryLimited) != 0;
    summary.pairConstraintLimited +=
        (quality[i].flags & QualityPairConstraintLimited) != 0;
  }
  summary.meanCorrection = count ? total / static_cast<float>(count) : 0.0f;
  return summary;
}

} // namespace

namespace {
// Solve both curves together: curvature regularization, robust observations,
// independently observed wall strips, and pair width in the same optimization.
bool fitTrackedCurves(SeamResult &result, const DetectorConfig &config) {
  const int n=int(result.leftSeam->size());
  if (!n || result.rightSeam->size()!=std::size_t(n)) return false;
  const auto axes=localTransverseAxes(config.longitudinalAxis);
  const int station=int(config.longitudinalAxis), size=4*n;
  const std::array<CloudT*,2> clouds{{result.leftSeam.get(),result.rightSeam.get()}};
  const std::array<std::vector<SeamPointQuality>*,2> qualities{{&result.leftPointQuality,&result.rightPointQuality}};
  const std::array<const std::vector<float>*,2> confidences{{&result.leftPointConfidence,&result.rightPointConfidence}};
  const std::array<const std::vector<std::size_t>*,2> segments{{&result.leftSegmentStarts,&result.rightSegmentStarts}};
  bool tracked=false;
  for (const auto *q:qualities) for (const auto &p:*q) tracked|=p.trackedWall;
  if (!tracked) return false;
  Eigen::VectorXd observed(size),base(size),weights(size);
  const Eigen::Vector2d origin((*clouds[0])[0].data[axes[0]],(*clouds[0])[0].data[axes[1]]);
  std::vector<Eigen::Triplet<double>> penalty, constraints;
  std::vector<double> lower,upper;
  const double lambda=std::pow(double(config.finalSmoothingLength)/(2*3.141592653589793),4);
  for (int side=0;side<2;++side) {
    const auto &cloud=*clouds[side];
    auto starts=normalizedStarts(*segments[side],n);
    std::vector<float> gaps;
    for (int i=1;i<n;++i) if (cloud[i].data[station]>cloud[i-1].data[station])
      gaps.push_back(cloud[i].data[station]-cloud[i-1].data[station]);
    std::sort(gaps.begin(),gaps.end());
    const double maxGap=gaps.empty()?0:gaps[gaps.size()/2]*std::max(2.,double(config.maxInterpolatedGapStations)+1.5);
    for (int i=0;i<n;++i) {
      const int index=4*i+2*side;
      const double score=i<int(confidences[side]->size()) ? (*confidences[side])[i] : 1;
      for (int dim=0;dim<2;++dim) { observed[index+dim]=double(cloud[i].data[axes[dim]])-origin[dim]; base[index+dim]=std::max(.1,std::min(1.,score)); }
      if (i<int(qualities[side]->size())) {
        const auto &q=(*qualities[side])[i];
        if (q.trackedWall && q.trackedNormal.norm()>.9f) {
          const int row=int(lower.size());
          for (int dim=0;dim<2;++dim) constraints.emplace_back(row,index+dim,q.trackedNormal[dim]);
          const double offset=q.trackedOffset-q.trackedNormal.cast<double>().dot(origin);
          lower.push_back(offset-config.profileWallTolerance); upper.push_back(offset+config.profileWallTolerance);
        }
      }
      if (i<1 || i+1>=n || std::binary_search(starts.begin(),starts.end(),i) ||
          std::binary_search(starts.begin(),starts.end(),i+1)) continue;
      const double h0=double(cloud[i].data[station])-cloud[i-1].data[station];
      const double h1=double(cloud[i+1].data[station])-cloud[i].data[station];
      if (h0<=1.e-5 || h1<=1.e-5 || h0>maxGap || h1>maxGap) continue;
      const double d[3]={2/(h0*(h0+h1)),-2/(h0*h1),2/(h1*(h0+h1))};
      for (int dim=0;dim<2;++dim) for (int r=0;r<3;++r) for (int c=0;c<3;++c)
        penalty.emplace_back(4*(i+r-1)+2*side+dim,4*(i+c-1)+2*side+dim,lambda*d[r]*d[c]);
    }
  }
  for (int i=0;i<n;++i) {
    const int row=int(lower.size());
    constraints.emplace_back(row,4*i,-1); constraints.emplace_back(row,4*i+2,1);
    lower.push_back(config.minSeamSeparation); upper.push_back(std::numeric_limits<double>::infinity());
  }
  Eigen::SparseMatrix<double> a(int(lower.size()),size),regularizer(size,size);
  a.setFromTriplets(constraints.begin(),constraints.end()); regularizer.setFromTriplets(penalty.begin(),penalty.end());
  const double rho=100.;
  const Eigen::SparseMatrix<double> normal=rho*(a.transpose()*a);
  Eigen::VectorXd x=observed,z=a*x,dual=Eigen::VectorXd::Zero(a.rows());
  auto clip=[&](Eigen::VectorXd &v){for (int j=0;j<v.size();++j) v[j]=std::max(lower[j],std::min(upper[j],v[j]));};
  clip(z); weights=base;
  for (int pass=0;pass<8;++pass) {
    Eigen::SparseMatrix<double> system=regularizer+normal;
    for (int j=0;j<size;++j) system.coeffRef(j,j)+=weights[j];
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(system);
    if (solver.info()!=Eigen::Success) throw std::runtime_error("Wall-constrained curve factorization failed");
    const Eigen::VectorXd data=weights.cwiseProduct(observed);
    for (int iteration=0;iteration<600;++iteration) {
      x=solver.solve(data+rho*a.transpose()*(z-dual));
      if (solver.info()!=Eigen::Success || !x.allFinite()) throw std::runtime_error("Wall-constrained curve solve failed");
      const Eigen::VectorXd value=a*x,old=z;
      z=value+dual; clip(z); dual+=value-z;
      if ((value-z).lpNorm<Eigen::Infinity>()<1.e-5 && (z-old).lpNorm<Eigen::Infinity>()<1.e-5) break;
    }
    for (int j=0;j<size;j+=2) {
      const double residual=(x.segment<2>(j)-observed.segment<2>(j)).norm();
      weights[j]=weights[j+1]=base[j]*std::min(1.,double(config.finalSmoothingRobustScale)/std::max(1.e-9,residual));
    }
  }
  Eigen::VectorXd check=a*x,bounded=check; clip(bounded);
  if ((check-bounded).lpNorm<Eigen::Infinity>()>1.e-3)
    throw std::runtime_error("Observed walls and requested seam width have no converged feasible curve");
  for (int side=0;side<2;++side) for (int i=0;i<n;++i) for (int dim=0;dim<2;++dim)
    (*clouds[side])[i].data[axes[dim]]=float(x[4*i+2*side+dim]+origin[dim]);
  return true;
}

void fitFinalCurve(CloudT &points, const std::vector<float> &confidence,
                   const std::vector<std::size_t> &segmentStarts,
                   const DetectorConfig &config) {
  if (points.size() < 3) return;
  const int stationAxis = int(config.longitudinalAxis);
  const auto axes = localTransverseAxes(config.longitudinalAxis);
  auto starts = normalizedStarts(segmentStarts, points.size());
  // Never smooth across missing runs or a change of scan direction.
  std::vector<float> gaps;
  for (std::size_t i = 1; i < points.size(); ++i) {
    const float gap = points[i].data[stationAxis] - points[i-1].data[stationAxis];
    if (std::isfinite(gap) && gap > 1.e-5f) gaps.push_back(gap);
  }
  if (gaps.empty()) return;
  std::sort(gaps.begin(), gaps.end());
  const double maxGap = gaps[gaps.size()/2] *
      std::max(2.0, double(config.maxInterpolatedGapStations) + 1.5);
  for (std::size_t i = 1; i < points.size(); ++i) {
    const double gap = double(points[i].data[stationAxis]) - points[i-1].data[stationAxis];
    if (!std::isfinite(gap) || gap <= 1.e-5 || gap > maxGap) starts.push_back(i);
  }
  std::sort(starts.begin(), starts.end());
  starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
  starts.push_back(points.size());
  const double lambda = std::pow(double(config.finalSmoothingLength) / (2.0 * 3.141592653589793), 4);
  for (std::size_t segment = 1; segment < starts.size(); ++segment) {
    const std::size_t begin = starts[segment-1];
    const int n = int(starts[segment] - begin);
    if (n < 3) continue;
    Eigen::MatrixXd observed(n, 2), fitted(n, 2);
    Eigen::VectorXd base(n), weights(n);
    // Subtract an origin for numerical stability with large world coordinates.
    const Eigen::RowVector2d origin(points[begin].data[axes[0]], points[begin].data[axes[1]]);
    for (int i = 0; i < n; ++i) {
      observed(i,0) = double(points[begin+i].data[axes[0]]) - origin[0];
      observed(i,1) = double(points[begin+i].data[axes[1]]) - origin[1];
      const float score = begin+i < confidence.size() ? confidence[begin+i] : 1.0f;
      base[i] = std::isfinite(score) ? std::max(.1, std::min(1.0, double(score))) : .1;
    }
    weights = base;
    std::vector<Eigen::Triplet<double>> penalty;
    penalty.reserve(9 * n);
    for (int i = 1; i + 1 < n; ++i) {
      const double h0 = double(points[begin+i].data[stationAxis]) - points[begin+i-1].data[stationAxis];
      const double h1 = double(points[begin+i+1].data[stationAxis]) - points[begin+i].data[stationAxis];
      // Physical second derivative: exact zero for a sloped straight seam,
      // including nonuniform station spacing.
      const double d[3] = {2/(h0*(h0+h1)), -2/(h0*h1), 2/(h1*(h0+h1))};
      for (int row=0; row<3; ++row)
        for (int col=0; col<3; ++col)
          penalty.emplace_back(i+row-1, i+col-1, lambda*d[row]*d[col]);
    }
    bool valid = false;
    for (int pass = 0; pass < 10; ++pass) {
      auto entries = penalty;
      for (int i=0; i<n; ++i) entries.emplace_back(i,i,weights[i]);
      Eigen::SparseMatrix<double> system(n,n);
      system.setFromTriplets(entries.begin(), entries.end());
      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(system);
      if (solver.info() != Eigen::Success) { valid = false; break; }
      fitted = solver.solve(weights.asDiagonal() * observed);
      valid = solver.info() == Eigen::Success && fitted.allFinite();
      if (!valid) break;
      for (int i=0; i<n; ++i) {
        const double residual = (fitted.row(i)-observed.row(i)).norm();
        // Vector Huber weights share the outlier decision across U and V.
        weights[i] = base[i] * std::min(1.0, double(config.finalSmoothingRobustScale) / std::max(1.e-9,residual));
      }
    }
    if (valid) for (int i=0; i<n; ++i) {
      points[begin+i].data[axes[0]] = float(fitted(i,0)+origin[0]);
      points[begin+i].data[axes[1]] = float(fitted(i,1)+origin[1]);
    }
  }
}

FinalCurveSummary finalCurveSummary(const CloudT &before, const CloudT &after,
                                    std::vector<SeamPointQuality> &quality,
                                    const DetectorConfig &config) {
  FinalCurveSummary summary;
  summary.applied = before.size() >= 3;
  double total = 0;
  for (std::size_t i=0; i<after.size(); ++i) {
    const float distance = (toVector(after[i])-toVector(before[i])).norm();
    total += distance;
    summary.maxCorrection = std::max(summary.maxCorrection,distance);
    summary.outlierPoints += distance > config.finalSmoothingRobustScale;
    if (i < quality.size()) {
      quality[i].finalCurveCorrection = distance;
      if (quality[i].trackedWall) {
        const auto axes=localTransverseAxes(config.longitudinalAxis);
        const Eigen::Vector2f point(after[i].data[axes[0]],after[i].data[axes[1]]);
        ++summary.wallConstrainedPoints;
        summary.maximumWallResidual=std::max(summary.maximumWallResidual,
            std::fabs(quality[i].trackedNormal.dot(point)-quality[i].trackedOffset));
      }
      if (distance > 1.e-5f) {
        quality[i].flags |= QualityFinalCurveSmoothed;
        // Optional local refits are distinct from the independent wall strips
        // retained above; moving the final point invalidates their acceptance.
        quality[i].wallAccepted = quality[i].floorAccepted = false;
        quality[i].supportDistance = -1;
      }
    }
  }
  summary.meanCorrection = after.empty() ? 0 : float(total/after.size());
  return summary;
}
} // namespace

void applyFinalCurveSmoothing(SeamResult &result, const DetectorConfig &config) {
  if (!config.finalSmoothingEnabled || !result.leftSeam || !result.rightSeam ||
      !std::isfinite(config.finalSmoothingLength) || config.finalSmoothingLength <= 0 ||
      !std::isfinite(config.finalSmoothingRobustScale) || config.finalSmoothingRobustScale <= 0)
    return;
  const CloudT leftBefore(*result.leftSeam), rightBefore(*result.rightSeam);
  const bool tracked=fitTrackedCurves(result,config);
  if (!tracked) {
    fitFinalCurve(*result.leftSeam, result.leftPointConfidence, result.leftSegmentStarts, config);
    fitFinalCurve(*result.rightSeam, result.rightPointConfidence, result.rightSegmentStarts, config);
  }
  const int across = localTransverseAxes(config.longitudinalAxis)[0];
  const int station = int(config.longitudinalAxis);
  // Enforce width around the smooth center, never restore individual noisy points.
  for (std::size_t i=0; !tracked && i<std::min(result.leftSeam->size(),result.rightSeam->size()); ++i) {
    auto &left = (*result.leftSeam)[i]; auto &right = (*result.rightSeam)[i];
    if (std::abs(left.data[station]-right.data[station]) > 1.e-4f) continue;
    if (right.data[across]-left.data[across] >= config.minSeamSeparation) continue;
    const float center = .5f*(left.data[across]+right.data[across]);
    left.data[across] = center-.5f*config.minSeamSeparation;
    right.data[across] = center+.5f*config.minSeamSeparation;
    if (i<result.leftPointQuality.size()) result.leftPointQuality[i].flags |= QualityPairConstraintLimited;
    if (i<result.rightPointQuality.size()) result.rightPointQuality[i].flags |= QualityPairConstraintLimited;
  }
  result.leftFinalCurve = finalCurveSummary(leftBefore,*result.leftSeam,result.leftPointQuality,config);
  result.rightFinalCurve = finalCurveSummary(rightBefore,*result.rightSeam,result.rightPointQuality,config);
}

namespace {
struct EndpointLine {
  Eigen::Vector3d point=Eigen::Vector3d::Zero(), direction=Eigen::Vector3d::Zero();
  double rms=0;
  bool valid=false;
};
EndpointLine endpointLine(const CloudT &cloud, std::size_t begin, std::size_t end,
                          const std::vector<float> &stations, double minimumSpan) {
  EndpointLine line;
  if (end-begin<3 || stations[end-1]-stations[begin]<minimumSpan) return line;
  const Eigen::Vector3d origin=toVector(cloud[begin]).cast<double>();
  double sum=0, squares=0; Eigen::Vector3d mean=Eigen::Vector3d::Zero(), cross=Eigen::Vector3d::Zero();
  for (auto i=begin;i<end;++i) {
    const Eigen::Vector3d delta=toVector(cloud[i]).cast<double>()-origin;
    const double t=stations[i]-stations[begin]; sum+=t; squares+=t*t; mean+=delta; cross+=t*delta;
  }
  const double count=double(end-begin), variance=squares-sum*sum/count;
  if (variance<=1.e-8) return line;
  line.point=origin+mean/count;
  const Eigen::Vector3d velocity=(cross-sum*mean/count)/variance;
  if (velocity.norm()<1.e-8) return line;
  line.direction=velocity.normalized();
  for (auto i=begin;i<end;++i) {
    const Eigen::Vector3d delta=toVector(cloud[i]).cast<double>()-line.point;
    line.rms+=(delta-line.direction*delta.dot(line.direction)).squaredNorm();
  }
  line.rms=std::sqrt(line.rms/count);
  line.valid=line.point.allFinite() && line.direction.allFinite();
  return line;
}

// Walk from the stable inner portion toward one endpoint. A run must persist
// over a physical distance, so isolated points or sampling density do not vote.
std::size_t endpointCut(const CloudT &left, const CloudT &right,
                        std::size_t begin, std::size_t end, bool head,
                        const DetectorConfig &config, const std::vector<float> &stations) {
  const auto uncut=head ? begin : end;
  if (end-begin<9) return uncut;
  const double window=config.endpointTrimWindow;
  const double first=stations[begin], last=stations[end-1];
  if (last-first<4*window) return uncut;
  auto lower=[&](double station) {
    return std::size_t(std::lower_bound(stations.begin()+begin,stations.begin()+end,station)-stations.begin());
  };
  // Start beside the central reference span, not an arbitrary terminal
  // percentage: a long defective tail may start well before the last quarter.
  const double referenceHalfLength=std::max(window,.1*(last-first));
  const double boundary=.5*(first+last)+(head ? -referenceHalfLength : referenceHalfLength);
  auto unstableAt=[&](std::size_t index, std::size_t *cutIndex) {
    const double station=stations[index];
    const auto innerBegin=head ? index : lower(station-2*window);
    const auto innerEnd=head ? std::min(end,lower(station+2*window)+1) : index+1;
    const auto outerBegin=head ? lower(station-window) : index;
    const auto outerEnd=head ? index+1 : std::min(end,lower(station+window)+1);
    bool unstable=false;
    for (const CloudT *cloud:{&left,&right}) {
      const auto inner=endpointLine(*cloud,innerBegin,innerEnd,stations,window);
      const auto outer=endpointLine(*cloud,outerBegin,outerEnd,stations,.25*window);
      if (!inner.valid || !outer.valid || inner.rms>.5*config.endpointTrimDeviation) continue;
      const double angle=std::acos(std::max(-1.0,std::min(1.0,inner.direction.dot(outer.direction))))*180/3.141592653589793;
      std::vector<double> deviations;
      for (auto i=outerBegin;i<outerEnd;++i) {
        const Eigen::Vector3d delta=toVector((*cloud)[i]).cast<double>()-inner.point;
        deviations.push_back((delta-inner.direction*delta.dot(inner.direction)).norm());
      }
      const auto middle=deviations.begin()+deviations.size()/2;
      std::nth_element(deviations.begin(),middle,deviations.end());
      // A parallel but displaced tail is a jump even when its tangent agrees.
      const bool changed=(angle>=config.endpointTrimAngleDegrees && *middle>=config.endpointTrimDeviation) ||
                         *middle>=3*config.endpointTrimDeviation;
      if (!changed) continue;
      // An interior bump that returns to the reference is not a bad terminal
      // run. Require a physical interval of recovery, not a single crossing.
      const auto recoveryBegin=head ? lower(station-3*window) : lower(station+window);
      const auto recoveryEnd=head ? lower(station-window) : std::min(end,lower(station+3*window)+1);
      std::size_t nearStart=end; bool recovered=false;
      for (auto j=recoveryBegin;j<recoveryEnd;++j) {
        const Eigen::Vector3d delta=toVector((*cloud)[j]).cast<double>()-inner.point;
        const double distance=(delta-inner.direction*delta.dot(inner.direction)).norm();
        if (distance<config.endpointTrimDeviation) {
          if (nearStart==end) nearStart=j;
          if (stations[j]-stations[nearStart]>=.25*window) { recovered=true; break; }
        } else nearStart=end;
      }
      if (recovered) continue;
      unstable=true;
      if (cutIndex) {
        // The outer window detects a change ahead of this index. Locate the
        // actual sustained departure against THIS stable inner reference;
        // moving the reference into the bad run would mask steps and hooks.
        std::size_t departure=end;
        std::size_t j=head ? outerEnd-1 : outerBegin;
        for (;;) {
          const Eigen::Vector3d delta=toVector((*cloud)[j]).cast<double>()-inner.point;
          const double distance=(delta-inner.direction*delta.dot(inner.direction)).norm();
          if (distance>=config.endpointTrimDeviation) {
            if (departure==end) departure=j;
            if (std::fabs(stations[j]-stations[departure])>=.25*window) {
              const auto cut=head ? departure+1 : departure;
              *cutIndex=head ? std::max(*cutIndex,cut) : std::min(*cutIndex,cut);
              break;
            }
          } else departure=end;
          if ((head && j==outerBegin) || (!head && j+1==outerEnd)) break;
          if (head) --j; else ++j;
        }
      }
    }
    return unstable;
  };
  std::size_t index=std::min(end-1,lower(boundary)), run=end;
  std::size_t cutIndex=uncut;
  bool armed=!unstableAt(index,nullptr);
  // If the reference-zone edge already lies on a bend, look inward by at
  // most one window for its onset. Never use that zone edge as an artificial cut.
  if (!armed) {
    const auto initial=index;
    while ((head ? index+1<end : index>begin) && std::fabs(stations[index]-stations[initial])<window) {
      if (head) ++index; else --index;
      if (!unstableAt(index,nullptr)) { armed=true; break; }
    }
    if (!armed) index=initial;
  }
  for (;;) {
    const double station=stations[index];
    std::size_t observedCut=uncut;
    const bool unstable=unstableAt(index,&observedCut);
    if (!unstable) armed=true;
    if (unstable && armed) {
      if (run==end) run=index;
      cutIndex=head ? std::max(cutIndex,observedCut) : std::min(cutIndex,observedCut);
      if (cutIndex!=uncut && std::fabs(station-stations[run])>=.25*window)
        return cutIndex;
    } else { run=end; cutIndex=uncut; }
    if ((head && index==begin) || (!head && index+1==end)) break;
    if (head) --index; else ++index;
  }
  return uncut;
}

template<class T>
void retainEndpointRange(std::vector<T> &values, std::size_t begin, std::size_t end) {
  if (values.empty()) return;
  const auto first=std::min(begin,values.size()), last=std::min(end,values.size());
  values=std::vector<T>(values.begin()+first,values.begin()+last);
}
void retainEndpointStarts(std::vector<std::size_t> &starts,std::size_t begin,std::size_t end) {
  std::vector<std::size_t> kept;
  if (end>begin) kept.push_back(0);
  for (auto start:starts) if (start>begin && start<end) kept.push_back(start-begin);
  starts=std::move(kept);
}
} // namespace

SeamWidthDecision classifySeamWidth(const SeamResult &result, const DetectorConfig &config) {
  SeamWidthDecision decision;
  if (!config.adaptiveSeamCount || !result.leftSeam || !result.rightSeam ||
      result.leftSeam->empty() || result.leftSeam->size()!=result.rightSeam->size()) return decision;
  const int axis=localTransverseAxes(config.longitudinalAxis)[0];
  for (std::size_t i=0;i<result.leftSeam->size();++i) {
    const float width=(*result.rightSeam)[i].data[axis]-(*result.leftSeam)[i].data[axis];
    if (width>=config.doubleSeamWidth) ++decision.doubleVotes;
    else ++decision.singleVotes;
  }
  decision.mode=decision.doubleVotes>decision.singleVotes ? 2 : 1;
  return decision;
}

void trimSeamEndpoints(SeamResult &result, const DetectorConfig &config, const std::vector<float> *scanStations,
                       const SeamWidthDecision *decision) {
  result.endpointTrimmedStart=result.endpointTrimmedEnd=0;
  if (!config.trimUnstableEndpoints || !result.leftSeam || !result.rightSeam ||
      result.leftSeam->size()!=result.rightSeam->size() || result.leftSeam->size()<9 ||
      !std::isfinite(config.endpointTrimWindow) || config.endpointTrimWindow<=0 ||
      !std::isfinite(config.endpointTrimDeviation) || config.endpointTrimDeviation<=0 ||
      !std::isfinite(config.endpointTrimAngleDegrees) || config.endpointTrimAngleDegrees<=0 ||
      config.endpointTrimAngleDegrees>=90) return;
  const auto count=result.leftSeam->size();
  const auto &left=*result.leftSeam; const auto &right=*result.rightSeam;
  const int axis=int(config.longitudinalAxis);
  std::vector<float> stations;
  if (scanStations) {
    if (scanStations->size()!=count) return;
    stations=*scanStations;
  } else for (const auto &p:left) stations.push_back(p.data[axis]);
  for (float value:stations) if (!std::isfinite(value)) return;
  auto starts=normalizedStarts(result.leftSegmentStarts,count);
  auto rightStarts=normalizedStarts(result.rightSegmentStarts,count);
  starts.insert(starts.end(),rightStarts.begin(),rightStarts.end());
  std::vector<float> gaps;
  for (std::size_t i=1;i<count;++i) {
    const float gap=stations[i]-stations[i-1];
    if (gap>1.e-5f && std::isfinite(gap)) gaps.push_back(gap);
  }
  if (gaps.empty()) return;
  std::sort(gaps.begin(),gaps.end());
  const double maxGap=gaps[gaps.size()/2]*std::max(2.0,double(config.maxInterpolatedGapStations)+1.5);
  for (std::size_t i=1;i<count;++i) {
    const double gap=stations[i]-stations[i-1];
    if (gap<=1.e-5 || gap>maxGap) starts.push_back(i);
  }
  starts.push_back(count);
  std::sort(starts.begin(),starts.end()); starts.erase(std::unique(starts.begin(),starts.end()),starts.end());
  const auto widthDecision=decision ? *decision : classifySeamWidth(result,config);
  CloudT center;
  if (widthDecision.mode==1) {
    center.reserve(count);
    for (std::size_t i=0;i<count;++i)
      center.push_back(toPoint(.5f*(toVector(left[i])+toVector(right[i]))));
  }
  // A single seam follows the center. Width changes alone must not stop it
  // when both independently supported boundaries still define a stable center.
  const auto &trimLeft=widthDecision.mode==1 ? center : left;
  const auto &trimRight=widthDecision.mode==1 ? center : right;
  const auto begin=endpointCut(trimLeft,trimRight,0,starts[1],true,config,stations);
  const auto end=endpointCut(trimLeft,trimRight,starts[starts.size()-2],count,false,config,stations);
  if (begin==0 && end==count) return;
  // Only remove observations; do not extrapolate replacement endpoints.
  CloudT::Ptr keptLeft(new CloudT), keptRight(new CloudT);
  for (auto i=begin;i<end;++i) { keptLeft->push_back(left[i]); keptRight->push_back(right[i]); }
  result.leftSeam=keptLeft; result.rightSeam=keptRight;
  retainEndpointRange(result.leftFrameIndices,begin,end); retainEndpointRange(result.rightFrameIndices,begin,end);
  retainEndpointRange(result.leftPointConfidence,begin,end); retainEndpointRange(result.rightPointConfidence,begin,end);
  retainEndpointRange(result.leftPointQuality,begin,end); retainEndpointRange(result.rightPointQuality,begin,end);
  retainEndpointStarts(result.leftSegmentStarts,begin,end); retainEndpointStarts(result.rightSegmentStarts,begin,end);
  if (!result.leftPointConfidence.empty()) result.leftConfidence=
      std::accumulate(result.leftPointConfidence.begin(),result.leftPointConfidence.end(),0.f)/result.leftPointConfidence.size();
  if (!result.rightPointConfidence.empty()) result.rightConfidence=
      std::accumulate(result.rightPointConfidence.begin(),result.rightPointConfidence.end(),0.f)/result.rightPointConfidence.size();
  result.endpointTrimmedStart=begin; result.endpointTrimmedEnd=count-end;
}

void applySeamWidthRule(SeamResult &result, const DetectorConfig &config, const SeamWidthDecision *decision) {
  if (!config.adaptiveSeamCount || result.leftSeam->size()!=result.rightSeam->size() || result.leftSeam->empty()) return;
  const auto widthDecision=decision ? *decision : classifySeamWidth(result,config);
  const int axis=localTransverseAxes(config.longitudinalAxis)[0];
  const CloudT left(*result.leftSeam), right(*result.rightSeam);
  const auto leftIds=result.leftFrameIndices, rightIds=result.rightFrameIndices;
  const auto leftScores=result.leftPointConfidence, rightScores=result.rightPointConfidence;
  const auto leftQuality=result.leftPointQuality, rightQuality=result.rightPointQuality;
  const auto leftStarts=result.leftSegmentStarts, rightStarts=result.rightSegmentStarts;
  result.leftSeam.reset(new CloudT); result.rightSeam.reset(new CloudT); result.centerSeam.reset(new CloudT);
  result.leftFrameIndices.clear(); result.rightFrameIndices.clear(); result.centerFrameIndices.clear();
  result.leftPointConfidence.clear(); result.rightPointConfidence.clear(); result.centerPointConfidence.clear();
  result.leftPointQuality.clear(); result.rightPointQuality.clear(); result.centerPointQuality.clear();
  result.leftSegmentStarts.clear(); result.rightSegmentStarts.clear(); result.centerSegmentStarts.clear();
  result.minimumGrooveWidth=std::numeric_limits<float>::infinity(); result.maximumGrooveWidth=0;
  for (std::size_t i=0;i<left.size();++i) {
    const float width=right[i].data[axis]-left[i].data[axis];
    result.minimumGrooveWidth=std::min(result.minimumGrooveWidth,width);
    result.maximumGrooveWidth=std::max(result.maximumGrooveWidth,width);
  }
  // One mode for the complete seam. Equal votes default to a single track.
  const int mode=widthDecision.mode;
  for (std::size_t i=0;i<left.size();++i) {
    const bool split=i==0 ||
        std::find(leftStarts.begin(),leftStarts.end(),i)!=leftStarts.end() ||
        std::find(rightStarts.begin(),rightStarts.end(),i)!=rightStarts.end();
    const float l=i<leftScores.size()?leftScores[i]:0, r=i<rightScores.size()?rightScores[i]:0;
    if (mode==1) {
      if (split) result.centerSegmentStarts.push_back(result.centerSeam->size());
      result.centerSeam->push_back(toPoint(.5f*(toVector(left[i])+toVector(right[i]))));
      result.centerPointConfidence.push_back(.5f*(l+r));
      if (i<leftIds.size()) result.centerFrameIndices.push_back(leftIds[i]);
      SeamPointQuality quality; quality.sourceConfidence=.5f*(l+r);
      if (i<leftQuality.size() && i<rightQuality.size()) quality.flags=leftQuality[i].flags|rightQuality[i].flags;
      result.centerPointQuality.push_back(quality);
    } else {
      if (split) {
        result.leftSegmentStarts.push_back(result.leftSeam->size());
        result.rightSegmentStarts.push_back(result.rightSeam->size());
      }
      result.leftSeam->push_back(left[i]); result.rightSeam->push_back(right[i]);
      result.leftPointConfidence.push_back(l); result.rightPointConfidence.push_back(r);
      if (i<leftIds.size()) result.leftFrameIndices.push_back(leftIds[i]);
      if (i<rightIds.size()) result.rightFrameIndices.push_back(rightIds[i]);
      result.leftPointQuality.push_back(i<leftQuality.size()?leftQuality[i]:SeamPointQuality());
      result.rightPointQuality.push_back(i<rightQuality.size()?rightQuality[i]:SeamPointQuality());
    }
  }
  result.summary += " Whole-seam majority width rule: threshold="+std::to_string(config.doubleSeamWidth)+
      ", single/double votes="+std::to_string(widthDecision.singleVotes)+"/"+std::to_string(widthDecision.doubleVotes)+
      ", mode="+(mode==1 ? std::string("single") : std::string("double"))+
      " (ties single); center="+std::to_string(result.centerSeam->size())+
      ", left="+std::to_string(result.leftSeam->size())+", right="+std::to_string(result.rightSeam->size())+".";
}

void applySeamPostprocess(SeamResult &result, const DetectorConfig &config) {
  if (config.postprocessMode != SeamPostprocessMode::RobustRts ||
      !result.leftSeam || !result.rightSeam)
    return;
  const CloudT leftBefore(*result.leftSeam), rightBefore(*result.rightSeam);
  smoothSide(*result.leftSeam, result.leftPointConfidence,
             result.leftPointQuality, result.leftSegmentStarts, config);
  smoothSide(*result.rightSeam, result.rightPointConfidence,
             result.rightPointQuality, result.rightSegmentStarts, config);
  if (leftBefore.size() == result.leftSeam->size() &&
      rightBefore.size() == result.rightSeam->size() &&
      result.leftPointQuality.size() == result.leftSeam->size() &&
      result.rightPointQuality.size() == result.rightSeam->size())
    enforcePairSeparation(result, leftBefore, rightBefore, config);
  result.leftSmoothing =
      summarize(leftBefore, *result.leftSeam, result.leftPointQuality);
  result.rightSmoothing =
      summarize(rightBefore, *result.rightSeam, result.rightPointQuality);
}

} // namespace weld
