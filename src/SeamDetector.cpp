/**
 * 方案一算法阅读入口
 *
 * 目标：从槽形点云中输出左右两条焊缝，并使可靠点同时位于局部底面和侧壁。
 * 坐标：A 为 longitudinalAxis；U/V 是 transverseAxes(A) 返回的两个截面轴。
 *       L 是 U 较小的一侧，R 是 U 较大的一侧，与屏幕左右、全局上下无关。
 * 数据：入口采样保留实测点；后续主体、候选、精修和面提取共用采样点云。
 *       filteredCloud 用于搜索结构；局部侧壁/底面精修读取采样后的实测源点。
 *       不能把降采样中值点当作最后的实测平面支持，否则会再次产生偏移。
 *
 * 主流程索引（在编辑器中搜索对应 S 编号；文档使用相同编号）：
 * S00 输入检查/逐帧入口采样 -> H 船体主体提取 -> S01 有限点 -> S02 离群滤波 -> S03 体素降采样
 * 主体保护模式跳过二次滤波；debugFrameBatchSize>0 时只在 H、帧组 F 和最终结果暂停。
 * -> S04 纵向切片 -> S05 截面分箱 -> S06 三面结构候选 -> S07 局部折角
 * -> S08 保留主体 -> S09 成对路径优化 -> S10 内部短缺口插补
 * -> S11 曲线正则 -> S12 侧壁局部拟合 -> S13 高残差侧壁重拟合
 * -> S14 贴壁投影 -> S15 正方形底面拟合 -> S16 三约束交点 -> S17 输出。
 *
 * 关键约束：纵向站位固定；左右成对选取；没有可靠支撑就保留上一阶段结果；
 *           最终连续轨迹在候选保护后另行拟合，不再强制贴回逐帧局部面。
 * 阅读入口：detect -> collectFrameSlices/collectSlices -> retainMainSegment -> optimizePairedPaths
 * -> interpolateInternalGaps -> fillSideResult。单侧处理在 fillSideResult
 * 中展开。 debugFlag 开启时，界面/CLI 经 detectWithDebugViewer 在 GUI
 * 线程逐步显示 VTK。 算法只发布 DebugStep；窗口逻辑位于 SeamDebugViewer.cpp。
 */

#include "SeamDetector.h"
#include "UniformSeamWidth.h"
#include "SeamSmoother.h"
#include "ProfileCurvature.h"
#include "HullSurfaces.h"
#include "SurfaceClouds.h"
#include "InputSampling.h"

#include <pcl/common/common.h>
#include <pcl/common/point_tests.h>
#ifndef WELD_NO_PCL_FILTERS
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/filters/voxel_grid.h>
#endif

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <utility>
#include <vector>

namespace weld {
namespace {

constexpr float kPi = 3.14159265358979323846f;

struct ProfilePoint {
  PointT point;
  float u = 0.0f;
  float v = 0.0f;
  // 当前横向分箱内包含的原始点数。大面的扫描点通常比边缘噪声密集，
  // 拟合时使用这个权重可以避免稀疏的外轮廓把侧面方向拉偏。
  float support = 1.0f;
  // 分箱只用于建立稳定的横向顺序；均值和箱内离差保留原始二维分布，
  // 因此同一个 Y 箱内沿 Z 展开的近竖直侧面不会被压扁成一个点。
  Eigen::Vector2f fitMean = Eigen::Vector2f::Zero();
  Eigen::Matrix2f fitScatter = Eigen::Matrix2f::Zero();
};

struct Candidate {
  bool curvatureFallback = false;
  PointT point;
  Eigen::Vector2f wallNormal = Eigen::Vector2f::Zero();
  float wallOffset = 0, wallRms = 0;
  bool hasWall = false;
  bool trackedWall = false;
  Eigen::Vector2f trackedNormal = Eigen::Vector2f::Zero();
  float trackedOffset = 0;
  Eigen::Vector2f floorDirection = Eigen::Vector2f::UnitX();
  Eigen::Vector3f localOffset = Eigen::Vector3f::Zero();
  float feature = 0.0f;
  float profileCurvature = 0.0f;
  bool fromRegionEdges = false;
  bool fromLocalCorner = false;
  bool emergingProfile = false;
  std::array<float, 2> historyDistance{
      {std::numeric_limits<float>::infinity(),
       std::numeric_limits<float>::infinity()}};
  std::array<bool, 2> historyAvailable{{false, false}};
  std::array<float, 2> nodeCost{{0.0f, 0.0f}};
};

struct SurfacePair {
  float bottomMidpointU = std::numeric_limits<float>::quiet_NaN();
  bool centralBead = false;
  Candidate left;
  Candidate right;
  std::shared_ptr<DebugPair> debug;
};

struct Slice {
  UniformWidthSection widthSection;
  float station = 0.0f;
  int frameIndex = -1;
  std::size_t debugIndex = 0;
  bool hasGeometry = false;
  bool featureDepthBoundary = false;
  std::vector<SurfacePair> pairs;
};

// 二维剖面中的直线使用“单位方向 + 法向距离”表示。与 v=k*u+b 相比，
// 这种正交拟合即使面对近乎竖直的侧面也不会出现无穷大斜率或巨大的伪残差。
struct Line2D {
  Eigen::Vector2f center = Eigen::Vector2f::Zero();
  Eigen::Vector2f direction = Eigen::Vector2f::UnitX();
  Eigen::Vector2f normal = Eigen::Vector2f::UnitY();
  float offset = 0.0f;
  float rms = std::numeric_limits<float>::infinity();
  float support = 0.0f;
  bool valid = false;
};

// profile 中的半开区间 [begin, end)。区域生长只沿已经按横向坐标排序的
// 剖面相邻点进行，因而不会把三维点云中另一个纵向站位的近邻误接进来。
struct ProfileRegion {
  int begin = 0;
  int end = 0;
  Line2D seedLine;
  float quality = 0.0f;
  bool valid = false;
};

bool stepDebugEnabled(const DebugTrace *debug) {
  return debug && debug->config.debugFlag && bool(debug->onStep);
}

void publishDebugStep(DebugTrace *debug, const std::string &title,
                      std::vector<DebugCloudStage> layers,
                      const std::string &message = {},
                      const std::map<std::string, double> &metrics = {},
                      bool frameView = false,
                      CloudT::ConstPtr overviewCloud = {},
                      CloudT::ConstPtr positionCloud = {},
                      std::array<CloudT::ConstPtr, 2> overviewCandidates = {}) {
  if (!stepDebugEnabled(debug))
    return;
  if (debug->config.debugFrameBatchSize > 0 &&
      title.compare(0, 1, "F") != 0 && title.compare(0, 1, "H") != 0 &&
      title.compare(0, 3, "S17") != 0)
    return;
  DebugStep step;
  step.title = title;
  step.status = "completed";
  step.message = message;
  step.metrics = metrics;
  step.layers = std::move(layers);
  step.overviewCloud = std::move(overviewCloud);
  step.positionCloud = std::move(positionCloud);
  step.overviewCandidates = std::move(overviewCandidates);
  if (frameView) {
    const int axis = int(debug->config.longitudinalAxis);
    step.hasFrameView = true;
    step.viewU = Eigen::Vector3f::Zero();
    step.viewV = Eigen::Vector3f::Zero();
    step.viewU[axis == 0 ? 1 : 0] = 1;
    step.viewV[axis == 2 ? 1 : 2] = 1;
    if (!step.layers.empty() && step.layers.front().cloud && !step.layers.front().cloud->empty()) {
      Eigen::Vector3d center = Eigen::Vector3d::Zero();
      for (const auto &point : *step.layers.front().cloud) center += toVector(point).cast<double>();
      step.viewOrigin = (center / double(step.layers.front().cloud->size())).cast<float>();
    }
  }
  if (!debug->onStep(step))
    debug->onStep = {};
}

CloudT::Ptr pointsCloud(const std::vector<PointT> &points) {
  CloudT::Ptr cloud(new CloudT);
  cloud->reserve(points.size());
  for (const auto &point : points)
    cloud->push_back(point);
  return cloud;
}

CloudT::Ptr candidateCloud(const std::vector<Candidate> &path) {
  CloudT::Ptr cloud(new CloudT);
  cloud->reserve(path.size());
  for (const auto &candidate : path)
    cloud->push_back(candidate.point);
  return cloud;
}

void publishPathStep(DebugTrace *debug, const std::string &title,
                     const std::vector<Candidate> &left,
                     const std::vector<Candidate> &right) {
  if (stepDebugEnabled(debug))
    publishDebugStep(
        debug, title,
        {{"左焊缝", candidateCloud(left)}, {"右焊缝", candidateCloud(right)}});
}

void publishProfileStep(DebugTrace *debug,
                        const std::vector<ProfilePoint> &profile,
                        float station, const DetectorConfig &config) {
  if (!stepDebugEnabled(debug))
    return;
  CloudT::Ptr cloud(new CloudT);
  double supportSum = 0.0;
  int lowSupportBins = 0;
  for (const auto &point : profile) {
    cloud->push_back(point.point);
    supportSum += point.support;
    lowSupportBins += point.support < 3.0f;
  }
  const double geometryScale =
      config.profileSamplingMode == ProfileSamplingMode::StableGrid
          ? config.profileReferenceScale
          : config.profileBinSize;
  publishDebugStep(debug, "S05 截面分箱 · " + std::to_string(station),
                   {{"分箱轮廓", cloud}}, {},
                   {{"profile_bins", static_cast<double>(profile.size())},
                    {"mean_bin_support",
                     profile.empty() ? 0.0
                                     : supportSum / double(profile.size())},
                    {"low_support_bins", double(lowSupportBins)},
                    {"profile_bin_size", config.profileBinSize},
                    {"profile_geometry_scale", geometryScale}});
}

void publishPairStep(DebugTrace *debug, const std::vector<SurfacePair> &pairs,
                     float station) {
  if (!stepDebugEnabled(debug))
    return;
  std::vector<Candidate> left, right, corners;
  for (const auto &pair : pairs) {
    left.push_back(pair.left);
    right.push_back(pair.right);
    if (pair.left.fromLocalCorner)
      corners.push_back(pair.left);
    if (pair.right.fromLocalCorner)
      corners.push_back(pair.right);
  }
  publishPathStep(debug, "S06 三面结构候选 · " + std::to_string(station), left,
                  right);
  if (stepDebugEnabled(debug))
    publishDebugStep(debug, "S07 曲率 / 局部折角定位 · " + std::to_string(station),
                     {{"轮廓特征候选", candidateCloud(corners)}});
}

void publishRetainedSlices(DebugTrace *debug,
                           const std::vector<Slice> &slices) {
  if (!stepDebugEnabled(debug))
    return;
  std::vector<Candidate> left, right;
  for (const auto &slice : slices)
    for (const auto &pair : slice.pairs) {
      left.push_back(pair.left);
      right.push_back(pair.right);
    }
  publishPathStep(debug, "S08 保留主体连续段", left, right);
}

std::string sideStepTitle(const char *stage, int side, const PointT &point,
                          Axis axis) {
  return std::string(stage) + (side == 0 ? " · L · " : " · R · ") +
         std::to_string(coordinate(point, axis));
}

float clamp01(float value) { return std::max(0.0f, std::min(1.0f, value)); }

float median(std::vector<float> values) {
  if (values.empty()) {
    return 0.0f;
  }
  const std::size_t middle = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + middle, values.end());
  float result = values[middle];
  if (values.size() % 2 == 0) {
    const auto lower =
        std::max_element(values.begin(), values.begin() + middle);
    result = 0.5f * (result + *lower);
  }
  return result;
}

std::array<int, 2> transverseAxes(Axis longitudinal) {
  switch (longitudinal) {
  case Axis::X:
    return {{1, 2}};
  case Axis::Y:
    return {{0, 2}};
  case Axis::Z:
    return {{0, 1}};
  }
  return {{1, 2}};
}

float component(const PointT &point, int axis) {
  if (axis == 0)
    return point.x;
  if (axis == 1)
    return point.y;
  return point.z;
}

float transverseDistance(const PointT &lhs, const PointT &rhs,
                         Axis longitudinal) {
  Eigen::Vector3f difference = toVector(lhs) - toVector(rhs);
  difference[static_cast<int>(longitudinal)] = 0.0f;
  return difference.norm();
}

float transverseDistance(const PointT &lhs, const Eigen::Vector3f &rhs,
                         Axis longitudinal) {
  Eigen::Vector3f difference = toVector(lhs) - rhs;
  difference[static_cast<int>(longitudinal)] = 0.0f;
  return difference.norm();
}

float transverseOffsetDistance(const Candidate &lhs, const Candidate &rhs,
                               Axis longitudinal) {
  Eigen::Vector3f difference = lhs.localOffset - rhs.localOffset;
  difference[static_cast<int>(longitudinal)] = 0.0f;
  return difference.norm();
}

bool stableProfileSampling(const DetectorConfig &config) {
  return config.profileSamplingMode == ProfileSamplingMode::StableGrid;
}

// In compatibility mode profileBinSize keeps its original dual role.  The
// stable mode separates sampling resolution from the physical scale used by
// the detector's downstream geometry gates.
float profileGeometryScale(const DetectorConfig &config) {
  return std::max(
      0.05f, stableProfileSampling(config) ? config.profileReferenceScale
                                           : config.profileBinSize);
}

int scaledProfileCount(const DetectorConfig &config, int configuredCount,
                       int legacyMinimum, int stableMinimum = 3) {
  const int legacyCount = std::max(legacyMinimum, configuredCount);
  if (!stableProfileSampling(config))
    return legacyCount;
  const float physicalSpan =
      static_cast<float>(legacyCount) * profileGeometryScale(config);
  return std::max(stableMinimum,
                  static_cast<int>(std::lround(
                      physicalSpan / std::max(0.05f, config.profileBinSize))));
}

float maximumProfileGap(const DetectorConfig &config) {
  return std::max(0.2f, 3.0f * profileGeometryScale(config));
}

// S06 辅助：按每个分箱的原始支持点数、均值和箱内散布重建协方差。
// 最小特征值对应直线法向，sqrt(最小特征值) 为正交 RMS。
// 不采用普通 y=kx+b 拟合，因为侧壁可能接近竖直，斜率会发散。
Line2D fitOrthogonalLine(const std::vector<ProfilePoint> &profile, int begin,
                         int end, bool contourOnly = false) {
  Line2D line;
  if (begin < 0 || end > static_cast<int>(profile.size()) || end - begin < 3) {
    return line;
  }

  // 第一步：用每个横向分箱的原始点数加权求中心和协方差。
  for (int i = begin; i < end; ++i) {
    const ProfilePoint &point = profile[static_cast<std::size_t>(i)];
    const float weight = contourOnly ? 1.0f : std::max(1.0f, point.support);
    const Eigen::Vector2f mean =
        contourOnly ? Eigen::Vector2f(point.u, point.v) : point.fitMean;
    line.center += weight * mean;
    line.support += weight;
  }
  if (line.support < 3.0f) {
    return line;
  }
  line.center /= line.support;

  Eigen::Matrix2f covariance = Eigen::Matrix2f::Zero();
  for (int i = begin; i < end; ++i) {
    const ProfilePoint &point = profile[static_cast<std::size_t>(i)];
    const float weight = contourOnly ? 1.0f : std::max(1.0f, point.support);
    const Eigen::Vector2f mean =
        contourOnly ? Eigen::Vector2f(point.u, point.v) : point.fitMean;
    const Eigen::Vector2f difference = mean - line.center;
    if (contourOnly) {
      covariance += weight * difference * difference.transpose();
    } else {
      covariance +=
          point.fitScatter + weight * difference * difference.transpose();
    }
  }
  covariance /= line.support;

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> solver(covariance);
  if (solver.info() != Eigen::Success || solver.eigenvalues().y() < 1.0e-5f) {
    return line;
  }
  line.direction = solver.eigenvectors().col(1).normalized();
  // 统一让方向朝横向坐标增大的方向，后面才能判断左右侧面是否向相反方向弯折。
  if (line.direction.x() < 0.0f) {
    line.direction = -line.direction;
  }
  line.normal = Eigen::Vector2f(-line.direction.y(), line.direction.x());
  line.offset = line.normal.dot(line.center);

  line.rms =
      std::sqrt(std::max(0.0f, line.normal.dot(covariance * line.normal)));
  line.valid = std::isfinite(line.rms);
  return line;
}

Line2D fitRobustBottomLine(const std::vector<ProfilePoint> &profile, int begin,
                           int end, const DetectorConfig &config) {
  if (begin < 0 || end > int(profile.size()) || end - begin < 3)
    return {};
  const int count = end - begin;
  const int minimumSeparation = std::max(2, count / 3);
  const float inlierDistance = std::max(
      config.profileWallTolerance, .5f * profileGeometryScale(config));
  int bestInliers = 0;
  float bestResidual = std::numeric_limits<float>::infinity();
  Line2D bestSeed;
  for (int first = begin; first < end; ++first) {
    for (int second = first + minimumSeparation; second < end; ++second) {
      Eigen::Vector2f direction(profile[second].u-profile[first].u,
                                profile[second].v-profile[first].v);
      if (direction.norm() < 1.e-4f)
        continue;
      direction.normalize();
      if (direction.x() < 0) direction = -direction;
      if (std::fabs(direction.y()) /
              std::max(1.e-3f,std::fabs(direction.x())) >
          config.maxBottomSlope)
        continue;
      const Eigen::Vector2f normal(-direction.y(),direction.x());
      const float offset = normal.dot(
          Eigen::Vector2f(profile[first].u,profile[first].v));
      int inliers = 0;
      float residual = 0;
      for (int i=begin;i<end;++i) {
        const float distance = std::fabs(
            normal.dot(Eigen::Vector2f(profile[i].u,profile[i].v))-offset);
        if (distance <= inlierDistance) {
          ++inliers;
          residual += distance*distance;
        }
      }
      if (inliers > bestInliers ||
          (inliers == bestInliers && residual < bestResidual)) {
        bestInliers = inliers;
        bestResidual = residual;
        bestSeed.center = Eigen::Vector2f(profile[first].u,profile[first].v);
        bestSeed.direction = direction;
        bestSeed.normal = normal;
        bestSeed.offset = offset;
        bestSeed.valid = true;
      }
    }
  }
  if (!bestSeed.valid || bestInliers < 3)
    return fitOrthogonalLine(profile,begin,end);
  std::vector<ProfilePoint> inliers;
  for (int i=begin;i<end;++i) {
    const float distance = std::fabs(
        bestSeed.normal.dot(Eigen::Vector2f(profile[i].u,profile[i].v))-
        bestSeed.offset);
    if (distance <= inlierDistance)
      inliers.push_back(profile[i]);
  }
  const auto fitted=fitOrthogonalLine(inliers,0,int(inliers.size()),true);
  return fitted.valid ? fitted : bestSeed;
}

// A wavy floor can contain two similarly supported levels. Hard inlier-count
// voting switches between them when a single bin crosses the distance gate.
// Smooth L1 reweighting retains the complete floor while limiting outliers.
Line2D fitContinuousBottomLine(const std::vector<ProfilePoint> &profile,
                               int begin,int end,const DetectorConfig &config) {
  auto line=fitOrthogonalLine(profile,begin,end,true);
  if (!line.valid) return line;
  const float scale=std::max(config.profileWallTolerance,.5f*profileGeometryScale(config));
  for (int iteration=0;iteration<8;++iteration) {
    Eigen::Vector2f center=Eigen::Vector2f::Zero();
    std::vector<float> weights;
    float total=0;
    for (int i=begin;i<end;++i) {
      const Eigen::Vector2f p(profile[i].u,profile[i].v);
      const float residual=(line.normal.dot(p)-line.offset)/scale;
      const float weight=1/std::sqrt(1+residual*residual);
      weights.push_back(weight); total+=weight; center+=weight*p;
    }
    center/=total;
    Eigen::Matrix2f covariance=Eigen::Matrix2f::Zero();
    for (int i=begin;i<end;++i) {
      const Eigen::Vector2f delta=Eigen::Vector2f(profile[i].u,profile[i].v)-center;
      covariance+=weights[i-begin]*delta*delta.transpose();
    }
    covariance/=total;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> solver(covariance);
    if (solver.info()!=Eigen::Success || solver.eigenvalues().y()<1.e-5f) return {};
    line.center=center; line.support=total;
    line.direction=solver.eigenvectors().col(1).normalized();
    if (line.direction.x()<0) line.direction=-line.direction;
    line.normal={-line.direction.y(),line.direction.x()};
    line.offset=line.normal.dot(center);
    line.rms=std::sqrt(std::max(0.f,solver.eigenvalues().x()));
  }
  return line;
}

float acuteAngle(const Line2D &first, const Line2D &second) {
  const float cosine = std::max(
      0.0f, std::min(1.0f, std::fabs(first.direction.dot(second.direction))));
  return std::acos(cosine);
}

float cross2D(const Eigen::Vector2f &first, const Eigen::Vector2f &second) {
  return first.x() * second.y() - first.y() * second.x();
}

bool intersectLines(const Line2D &first, const Line2D &second,
                    Eigen::Vector2f *intersection) {
  if (!intersection || !first.valid || !second.valid) {
    return false;
  }
  Eigen::Matrix2f normals;
  normals.row(0) = first.normal.transpose();
  normals.row(1) = second.normal.transpose();
  if (std::fabs(normals.determinant()) < 1.0e-3f) {
    return false;
  }
  *intersection =
      normals.inverse() * Eigen::Vector2f(first.offset, second.offset);
  return intersection->allFinite();
}

// S05 截面轮廓：输入当前厚截面的原始采样，输出按 U 递增排列的分箱。
// 代表点取各坐标中位数供定位；另外保留 support/fitMean/fitScatter 供面拟合。
// 当箱内纵向跨度足够时扣除局部纵向线性趋势，减轻厚截面在弯曲区的混叠。
// 断点观察：points.size()、result.size()、每箱 support 和 fitScatter。
std::vector<ProfilePoint> makeProfile(const std::vector<PointT> &points,
                                      Axis longitudinal,
                                      const DetectorConfig &config,
                                      float station = 0.0f,
                                      float detrendMinimumSpan = 0.0f) {
  std::vector<ProfilePoint> result;
  if (points.size() < 5) {
    return result;
  }

  const std::array<int, 2> axes = transverseAxes(longitudinal);

  // 坐标系固定时，整条轨迹必须使用同一剖面定义。之前逐片选择变化范围较大
  // 的坐标，船体高度范围超过宽度时会从 Y 突然切换到 Z，产生分叉和跨面。
  // 当前 X 向扫描因此始终按 Y 排序并拟合 Z；其它纵轴也使用固定的剩余轴序。
  const int uAxis = axes[0];
  const int vAxis = axes[1];

  std::vector<PointT> sorted = points;
  std::sort(sorted.begin(), sorted.end(),
            [uAxis](const PointT &lhs, const PointT &rhs) {
              return component(lhs, uAxis) < component(rhs, uAxis);
            });

  const float safeBinSize = std::max(0.05f, config.profileBinSize);
  std::size_t begin = 0;
  while (begin < sorted.size()) {
    const float binStart = component(sorted[begin], uAxis);
    std::size_t end = begin + 1;
    if (stableProfileSampling(config)) {
      // All slices share profileGridOrigin. A point crossing one fixed cell
      // boundary can no longer move the start of every following bin.
      const auto binIndex = static_cast<long long>(
          std::floor((binStart - config.profileGridOrigin) / safeBinSize));
      while (end < sorted.size()) {
        const float coordinateU = component(sorted[end], uAxis);
        const auto nextIndex = static_cast<long long>(std::floor(
            (coordinateU - config.profileGridOrigin) / safeBinSize));
        if (nextIndex != binIndex)
          break;
        ++end;
      }
    } else {
      // Original rolling binning is intentionally retained byte-for-byte.
      while (end < sorted.size() &&
             component(sorted[end], uAxis) - binStart < safeBinSize) {
        ++end;
      }
    }

    std::vector<PointT> corrected(sorted.begin() + begin, sorted.begin() + end);
    // 同一横向箱里的点来自截面厚度内的多个纵向位置。分别估计每个箱的
    // Theil-Sen 斜率，可去掉船体纵向陡坡，又不会受截面覆盖范围变化影响。
    if (detrendMinimumSpan > 0.0f && corrected.size() >= 3) {
      const int longitudinalAxis = static_cast<int>(longitudinal);
      std::vector<float> slopes;
      for (std::size_t i = 0; i < corrected.size(); ++i) {
        for (std::size_t j = i + 1; j < corrected.size(); ++j) {
          const float delta = corrected[j].data[longitudinalAxis] -
                              corrected[i].data[longitudinalAxis];
          if (std::fabs(delta) >= detrendMinimumSpan) {
            slopes.push_back(
                (corrected[j].data[vAxis] - corrected[i].data[vAxis]) / delta);
          }
        }
      }
      if (!slopes.empty()) {
        const float slope = median(std::move(slopes));
        if (std::isfinite(slope) && std::fabs(slope) < 100.0f) {
          for (PointT &point : corrected) {
            point.data[vAxis] -=
                slope * (point.data[longitudinalAxis] - station);
          }
        }
      }
    }

    std::vector<float> xs;
    std::vector<float> ys;
    std::vector<float> zs;
    xs.reserve(end - begin);
    ys.reserve(end - begin);
    zs.reserve(end - begin);
    for (const PointT &point : corrected) {
      xs.push_back(point.x);
      ys.push_back(point.y);
      zs.push_back(point.z);
    }
    ProfilePoint sample;
    sample.point = PointT(median(xs), median(ys), median(zs));
    sample.u = component(sample.point, uAxis);
    sample.v = component(sample.point, vAxis);
    sample.support = static_cast<float>(end - begin);
    for (const PointT &point : corrected) {
      sample.fitMean +=
          Eigen::Vector2f(component(point, uAxis), component(point, vAxis));
    }
    sample.fitMean /= sample.support;
    for (const PointT &point : corrected) {
      const Eigen::Vector2f difference(
          component(point, uAxis) - sample.fitMean.x(),
          component(point, vAxis) - sample.fitMean.y());
      sample.fitScatter += difference * difference.transpose();
    }
    result.push_back(sample);
    begin = end;
  }
  return result;
}

// S06 辅助：从交界向底面内侧逐步扩大拟合窗口。
// 法向漂移过大时停止扩大，防止把局部弯曲底面压成跨槽统一平面。
Line2D fitAdaptiveBoundaryLine(const std::vector<ProfilePoint> &profile,
                               int anchor, int direction, int minimumWindow,
                               int maximumWindow,
                               const DetectorConfig &config) {
  Line2D selected;
  Line2D previous;
  const float maximumNormalChange =
      config.maxLocalFitNormalChangeDegrees * kPi / 180.0f;
  const float maximumRms = std::max(0.25f, config.surfaceRmsScale);
  for (int window = std::max(3, minimumWindow); window <= maximumWindow;
       ++window) {
    const int begin = direction > 0 ? anchor : anchor - window;
    const int end = direction > 0 ? anchor + window : anchor;
    const Line2D current = fitOrthogonalLine(profile, begin, end);
    if (!current.valid || current.rms > maximumRms) {
      continue;
    }
    if (previous.valid && acuteAngle(previous, current) > maximumNormalChange) {
      break;
    }
    selected = current;
    previous = current;
  }
  return selected;
}

float pointToLineDistance(const ProfilePoint &point, const Line2D &line) {
  if (!line.valid) {
    return std::numeric_limits<float>::infinity();
  }
  return std::fabs(line.normal.dot(point.fitMean) - line.offset);
}

// S07 辅助：从可靠面种子沿轮廓生长，逐点检查切线变化、到面距离和总法向漂移。
// 底面与侧壁分别生长；面内连续起伏可以保留，发生真实折角时停止。
ProfileRegion growProfileRegion(const std::vector<ProfilePoint> &profile,
                                int seedBegin, int seedEnd, int minimumIndex,
                                int maximumIndex, const DetectorConfig &config,
                                bool contourOnly = false) {
  ProfileRegion region;
  const int pointCount = static_cast<int>(profile.size());
  minimumIndex = std::max(0, minimumIndex);
  maximumIndex = std::min(pointCount, maximumIndex);
  seedBegin = std::max(minimumIndex, seedBegin);
  seedEnd = std::min(maximumIndex, seedEnd);
  if (seedEnd - seedBegin < 3) {
    return region;
  }

  region.begin = seedBegin;
  region.end = seedEnd;
  region.seedLine = fitOrthogonalLine(profile, seedBegin, seedEnd, contourOnly);
  if (!region.seedLine.valid) {
    return region;
  }

  // 轮廓中值仍有离散噪声。底面用较宽的局部窗口估计趋势，避免四个点的
  // 微小起伏使生长反复停止；真实侧面的方向突变仍由步进角和总漂移限制。
  const int window =
      scaledProfileCount(config, config.regionGrowingWindow, 3) *
      (contourOnly ? 2 : 1);
  const float maximumStepAngle =
      config.regionGrowingSmoothnessDegrees * kPi / 180.0f;
  const float maximumDrift =
      config.regionGrowingMaxNormalDriftDegrees * kPi / 180.0f;
  const float maximumDistance = std::max(0.1f, config.regionGrowingDistance);
  const float maximumGap = maximumProfileGap(config);
  float acceptedQuality = 0.0f;
  int acceptedCount = 0;

  auto acceptAt = [&](int index, int direction) {
    int localBegin = 0;
    int localEnd = 0;
    if (direction < 0) {
      localBegin = region.begin;
      localEnd = std::min(region.end, region.begin + window);
    } else {
      localBegin = std::max(region.begin, region.end - window);
      localEnd = region.end;
    }
    const Line2D local =
        fitOrthogonalLine(profile, localBegin, localEnd, contourOnly);
    if (!local.valid) {
      return false;
    }

    const int expandedBegin = direction < 0 ? index : localBegin;
    const int expandedEnd = direction < 0 ? localEnd : index + 1;
    const Line2D expanded =
        fitOrthogonalLine(profile, expandedBegin, expandedEnd, contourOnly);
    if (!expanded.valid || acuteAngle(local, expanded) > maximumStepAngle ||
        acuteAngle(region.seedLine, expanded) > maximumDrift) {
      return false;
    }

    const int adjacent = direction < 0 ? index + 1 : index - 1;
    if (std::fabs(profile[static_cast<std::size_t>(index)].u -
                  profile[static_cast<std::size_t>(adjacent)].u) > maximumGap) {
      return false;
    }
    const auto &sample = profile[static_cast<std::size_t>(index)];
    const float distance =
        contourOnly
            ? std::fabs(local.normal.dot(Eigen::Vector2f(sample.u, sample.v)) -
                        local.offset)
            : pointToLineDistance(sample, local);
    const float distanceTolerance =
        contourOnly ? std::max(maximumDistance, std::min(2.0f * maximumDistance,
                                                         2.5f * local.rms))
                    : maximumDistance;
    if (distance > distanceTolerance) {
      return false;
    }
    acceptedQuality += std::exp(-distance / distanceTolerance) *
                       std::exp(-acuteAngle(local, expanded) /
                                std::max(1.0e-3f, maximumStepAngle));
    ++acceptedCount;
    return true;
  };

  while (region.begin > minimumIndex && acceptAt(region.begin - 1, -1)) {
    --region.begin;
  }
  while (region.end < maximumIndex && acceptAt(region.end, +1)) {
    ++region.end;
  }

  region.quality =
      acceptedCount > 0
          ? clamp01(acceptedQuality / static_cast<float>(acceptedCount))
          : 0.5f;
  region.valid = region.end - region.begin >= 3;
  return region;
}

struct RegionEdgePair {
  Eigen::Vector2f sidePoint = Eigen::Vector2f::Zero();
  Eigen::Vector2f groundPoint = Eigen::Vector2f::Zero();
  float distance = std::numeric_limits<float>::infinity();
  bool valid = false;
};

// 对区域已经相交的圆角，取两个生长前沿之间的中间位置；存在缺口时，
// 则在两个前沿的小邻域内寻找真正的二维最近点。搜索始终局限在同一截面。
RegionEdgePair matchRegionEdges(const std::vector<ProfilePoint> &profile,
                                const ProfileRegion &side,
                                const ProfileRegion &ground, bool left,
                                const DetectorConfig &config) {
  RegionEdgePair result;
  if (!side.valid || !ground.valid) {
    return result;
  }

  const int sideFront = left ? side.end - 1 : side.begin;
  const int groundFront = left ? ground.begin : ground.end - 1;
  const bool overlaps =
      left ? sideFront >= groundFront : groundFront >= sideFront;
  if (overlaps) {
    const int middle = (sideFront + groundFront) / 2;
    const ProfilePoint &point = profile[static_cast<std::size_t>(middle)];
    result.sidePoint = Eigen::Vector2f(point.u, point.v);
    result.groundPoint = result.sidePoint;
    result.distance = 0.0f;
    result.valid = true;
    return result;
  }

  const int edgeSamples =
      scaledProfileCount(config, config.regionGrowingWindow, 2, 2);
  const int sideBegin =
      left ? std::max(side.begin, side.end - edgeSamples) : side.begin;
  const int sideEnd =
      left ? side.end : std::min(side.end, side.begin + edgeSamples);
  const int groundBegin =
      left ? ground.begin : std::max(ground.begin, ground.end - edgeSamples);
  const int groundEnd =
      left ? std::min(ground.end, ground.begin + edgeSamples) : ground.end;
  for (int sideIndex = sideBegin; sideIndex < sideEnd; ++sideIndex) {
    const Eigen::Vector2f sidePoint(
        profile[static_cast<std::size_t>(sideIndex)].u,
        profile[static_cast<std::size_t>(sideIndex)].v);
    for (int groundIndex = groundBegin; groundIndex < groundEnd;
         ++groundIndex) {
      const Eigen::Vector2f groundPoint(
          profile[static_cast<std::size_t>(groundIndex)].u,
          profile[static_cast<std::size_t>(groundIndex)].v);
      const float distance = (sidePoint - groundPoint).norm();
      if (distance < result.distance) {
        result.sidePoint = sidePoint;
        result.groundPoint = groundPoint;
        result.distance = distance;
      }
    }
  }
  result.valid = std::isfinite(result.distance) &&
                 result.distance <= config.regionGrowingMaxEdgeGap;
  return result;
}

Candidate candidateFromRegionEdges(const std::vector<ProfilePoint> &profile,
                                   const RegionEdgePair &edges,
                                   Axis longitudinal, float feature) {
  Candidate candidate;
  if (!edges.valid || profile.empty()) {
    return candidate;
  }
  const std::array<int, 2> axes = transverseAxes(longitudinal);
  const Eigen::Vector2f midpoint = 0.5f * (edges.sidePoint + edges.groundPoint);
  const auto nearest = std::min_element(
      profile.begin(), profile.end(),
      [&midpoint](const ProfilePoint &lhs, const ProfilePoint &rhs) {
        return (Eigen::Vector2f(lhs.u, lhs.v) - midpoint).squaredNorm() <
               (Eigen::Vector2f(rhs.u, rhs.v) - midpoint).squaredNorm();
      });
  candidate.point = nearest->point;
  candidate.point.data[axes[0]] = midpoint.x();
  candidate.point.data[axes[1]] = midpoint.y();
  candidate.feature = clamp01(feature);
  candidate.fromRegionEdges = true;
  return candidate;
}

Candidate candidateFromIntersection(const std::vector<ProfilePoint> &profile,
                                    const Eigen::Vector2f &intersection,
                                    Axis longitudinal, float feature,
                                    float profileBinSize) {
  Candidate candidate;
  if (profile.empty()) {
    return candidate;
  }
  const std::array<int, 2> axes = transverseAxes(longitudinal);
  std::vector<std::size_t> nearestIndices(profile.size());
  std::iota(nearestIndices.begin(), nearestIndices.end(), 0);
  std::sort(
      nearestIndices.begin(), nearestIndices.end(),
      [&profile, &intersection](std::size_t lhs, std::size_t rhs) {
        return (Eigen::Vector2f(profile[lhs].u, profile[lhs].v) - intersection)
                   .squaredNorm() <
               (Eigen::Vector2f(profile[rhs].u, profile[rhs].v) - intersection)
                   .squaredNorm();
      });
  const ProfilePoint &nearest = profile[nearestIndices.front()];

  candidate.point = nearest.point;
  // 单个最近点会随采样密度和毛刺在相邻箱之间跳动。取交点附近最多五个
  // 样本的坐标中值作为弱观测，理论交点仍保持主导。
  std::vector<float> nearbyU;
  std::vector<float> nearbyV;
  const float refineRadius = std::max(1.0f, 2.5f * profileBinSize);
  for (std::size_t index : nearestIndices) {
    const ProfilePoint &sample = profile[index];
    if ((Eigen::Vector2f(sample.u, sample.v) - intersection).norm() >
        refineRadius) {
      break;
    }
    nearbyU.push_back(sample.u);
    nearbyV.push_back(sample.v);
    if (nearbyU.size() >= 5) {
      break;
    }
  }
  const Eigen::Vector2f measured(
      nearbyU.empty() ? nearest.u : median(std::move(nearbyU)),
      nearbyV.empty() ? nearest.v : median(std::move(nearbyV)));
  const Eigen::Vector2f refined = 0.90f * intersection + 0.10f * measured;
  candidate.point.data[axes[0]] = refined.x();
  candidate.point.data[axes[1]] = refined.y();
  candidate.feature = clamp01(feature);
  return candidate;
}

Candidate nearestHistoryCandidate(const std::vector<ProfilePoint> &profile,
                                  const HistoryPrediction &prediction,
                                  const DetectorConfig &config) {
  Candidate candidate;
  if (!prediction.valid || profile.empty()) {
    return candidate;
  }
  const auto nearest = std::min_element(
      profile.begin(), profile.end(),
      [&prediction, &config](const ProfilePoint &lhs, const ProfilePoint &rhs) {
        return transverseDistance(lhs.point, prediction.point,
                                  config.longitudinalAxis) <
               transverseDistance(rhs.point, prediction.point,
                                  config.longitudinalAxis);
      });
  candidate.point = nearest->point;
  candidate.feature = 0.05f;
  return candidate;
}

struct LocalContourCorner {
  int index = -1;
  float angle = 0.0f;
  float quality = 0.0f;
  float curvature = 0.0f;
  Eigen::Vector2f point = Eigen::Vector2f::Zero();
  Line2D wall;
};

struct ObservedWall {
  Line2D line;
  Eigen::Vector2f foot = Eigen::Vector2f::Zero();
  float height = 0;
  float station = 0;
};

// Establish the wall independently of any curvature peak. A small depression
// in the floor cannot replace the long, steep, physically observed side wall.
void observeWalls(const std::vector<ProfileCurvatureSample> &curve,
                  const DetectorConfig &config, std::array<ObservedWall,2> &walls,
                  float bottomMidpoint) {
  if (curve.empty()) return;
  const int opening=config.interiorOpeningDirection ? config.interiorOpeningDirection :
      (curve[curve.size()/2].point.y()<=.5f*(curve.front().point.y()+curve.back().point.y()) ? 1 : -1);
  std::vector<float> arc(curve.size(),0);
  for (std::size_t i=1;i<curve.size();++i) arc[i]=arc[i-1]+(curve[i].point-curve[i-1].point).norm();
  for (int side=0;side<2;++side) {
    const float sign=float(opening*(side ? 1 : -1));
    for (std::size_t begin=0;begin<curve.size();) {
      if (!curve[begin].valid || sign*curve[begin].tangent.y()<.45f) { ++begin; continue; }
      std::size_t end=begin+1;
      while (end<curve.size()) {
        if (curve[end].valid && sign*curve[end].tangent.y()>=.45f) { ++end; continue; }
        // A short measured bead/ripple does not change the wall's identity.
        // Bound by arc length, so a broad floor cannot bridge to an inner pit.
        std::size_t next=end;
        while (next<curve.size() && arc[next]-arc[end-1]<=2*profileGeometryScale(config) &&
               (!curve[next].valid || sign*curve[next].tangent.y()<.45f)) ++next;
        if (next>=curve.size() || arc[next]-arc[end-1]>2*profileGeometryScale(config)) break;
        std::size_t nextEnd=next+1;
        while (nextEnd<curve.size() && curve[nextEnd].valid && sign*curve[nextEnd].tangent.y()>=.45f) ++nextEnd;
        // Both pieces must independently have wall-height support. Otherwise
        // this would extend a real wall into a tiny neighboring floor ripple.
        const float support=2*profileGeometryScale(config);
        if (sign*(curve[end-1].point.y()-curve[begin].point.y())<support ||
            sign*(curve[nextEnd-1].point.y()-curve[next].point.y())<support) break;
        end=nextEnd;
      }
      const float height=sign*(curve[end-1].point.y()-curve[begin].point.y());
      const std::size_t foot=side ? begin : end-1;
      const bool correctSide = side ? curve[foot].point.x()>bottomMidpoint
                                    : curve[foot].point.x()<bottomMidpoint;
      if (correctSide && height>=2*profileGeometryScale(config) && height>walls[side].height) {
        std::vector<ProfilePoint> patch;
        for (std::size_t j=begin;j<end;++j) {
          if ((curve[j].point-curve[foot].point).norm()>2*profileGeometryScale(config)) continue;
          ProfilePoint p; p.u=curve[j].point.x(); p.v=curve[j].point.y(); p.support=1;
          patch.push_back(p);
        }
        if (patch.size()>=3) {
          const auto line=fitOrthogonalLine(patch,0,int(patch.size()),true);
          if (line.valid) { walls[side].line=line; walls[side].height=height; walls[side].foot=curve[foot].point; }
        }
      }
      begin=end;
    }
  }
}

bool nearObservedWall(const ProfileTurnSample &turn, const ObservedWall &wall,
                      const DetectorConfig &config) {
  return wall.line.valid &&
      (turn.point-wall.foot).norm()<=std::min(config.profileCurvatureRadius,5*profileGeometryScale(config)) &&
      std::fabs(wall.line.normal.dot(turn.point)-wall.line.offset)<=config.profileWallTolerance;
}

int contourBoundaryRadius(const DetectorConfig &config) {
  return scaledProfileCount(config,config.fitHalfWindow,5)+
         scaledProfileCount(config,config.fitGuardPoints,1,1);
}

// 区域前沿之间有缺口或重叠时，中点未必是板面交界。只在已经通过三面
// 结构确认的候选附近，寻找“稳定侧面 -> 局部底面轮廓”的实测转折。
// 外侧短窗口必须与原始分布拟合的侧面同向，内侧短窗口必须明显转向；
// 代表点还必须贴近侧面且接近原候选，避免选到槽内起伏或外侧肩部。
LocalContourCorner
findLocalContourCorner(const std::vector<ProfilePoint> &profile,
                       const Line2D &sideSurface, const Line2D &broadBottom,
                       const Eigen::Vector2f &reference, int boundary,
                       bool left, const DetectorConfig &config,
                       const std::vector<ProfileCurvatureSample> &curvatures,
                       const std::vector<ProfileTurnSample> &turns,
                       const std::array<ObservedWall,2> &walls) {
  LocalContourCorner best;
  if (!config.useProfileCurvature && acuteAngle(sideSurface, broadBottom) <
      config.minContourCornerAngleDegrees * kPi / 180.0f)
    return best;
  const int window = scaledProfileCount(
      config, std::min(5, config.localAngleWindow), 3);
  const int n = static_cast<int>(profile.size());
  const int radius = contourBoundaryRadius(config);
  const float maxDistance = std::max(0.5f, config.regionGrowingDistance);
  const float minTurn = std::max(config.minLocalSurfaceAngleDegrees,
                                 config.minSurfaceAngleDegrees) *
                        kPi / 180.0f;
  const float maxSideAngle =
      std::max(12.0f, config.regionGrowingSmoothnessDegrees) * kPi / 180.0f;
  if (config.useProfileCurvature) {
    const float expectedTurn = left ? cross2D(sideSurface.direction, broadBottom.direction)
                                    : cross2D(broadBottom.direction, sideSurface.direction);
    float bestMagnitude = 0, bestWallDistance = std::numeric_limits<float>::infinity();
    for (const auto &turn : turns) {
      if (!nearObservedWall(turn,walls[left ? 0 : 1],config)) continue;
      if (std::abs(turn.index-boundary)>radius || turn.angle*expectedTurn<=0 ||
          std::fabs(turn.angle)<config.minLocalSurfaceAngleDegrees*kPi/180.f ||
          walls[left ? 0 : 1].line.rms>config.profileWallTolerance) continue;
      const auto &observed=walls[left ? 0 : 1].line;
      const auto &outer=observed.direction;
      const float wallDistance=std::fabs(observed.normal.dot(turn.point)-observed.offset);
      if ((turn.point-reference).norm()>config.surfaceIntersectionTolerance) continue;
      const float sideAngle=std::acos(std::min(1.f,std::fabs(outer.dot(sideSurface.direction))));
      if (sideAngle>std::max(25.f*kPi/180.f,maxSideAngle)) continue;
      const float magnitude=std::fabs(turn.curvature);
      // Maximize analytic B-spline curvature among supported wall peaks.
      if (magnitude<bestMagnitude-1.e-4f ||
          (std::fabs(magnitude-bestMagnitude)<=1.e-4f && wallDistance>=bestWallDistance)) continue;
      bestMagnitude=magnitude;
      bestWallDistance=wallDistance;
      best.index=turn.index;
      best.point=turn.point;
      best.curvature=turn.curvature;
      best.angle=std::fabs(turn.angle)*180.f/kPi;
      best.wall=walls[left ? 0 : 1].line;
      best.quality=clamp01(std::fabs(turn.angle)/(0.5f*kPi))*
                   std::exp(-observed.rms/std::max(.05f,config.profileWallTolerance));
    }
    return best;
  }
  for (int k = std::max(window - 1, boundary - radius);
       k <= std::min(n - window, boundary + radius); ++k) {
    const auto &sample = profile[static_cast<std::size_t>(k)];
    const Eigen::Vector2f point(sample.u, sample.v);
    const float sideDistance =
        std::fabs(sideSurface.normal.dot(point) - sideSurface.offset);
    if (sideDistance > maxDistance ||
        (point - reference).norm() > config.surfaceIntersectionTolerance)
      continue;
    // 不把没有采样的孔洞当成观测到的折角。
    bool gap = false;
    for (int j = k - window + 2; j < k + window; ++j) {
      if (profile[j].u - profile[j - 1].u > maximumProfileGap(config)) {
        gap = true;
        break;
      }
    }
    if (gap)
      continue;
    const Line2D before =
        fitOrthogonalLine(profile, k - window + 1, k + 1, true);
    const Line2D after = fitOrthogonalLine(profile, k, k + window, true);
    if (!before.valid || !after.valid)
      continue;
    const Line2D &outer = left ? before : after;
    const Line2D &inner = left ? after : before;
    const float sideAngle = acuteAngle(outer, sideSurface);
    const float turn = acuteAngle(outer, inner);
    if (sideAngle > maxSideAngle || turn < minTurn || outer.rms > maxDistance ||
        inner.rms > maxDistance)
      continue;
    const float quality =
        clamp01(turn / (0.5f * kPi)) * std::exp(-sideAngle / maxSideAngle) *
        std::exp(-sideDistance / maxDistance) *
        std::exp(-(outer.rms + inner.rms) / (2.0f * maxDistance));
    if (quality > best.quality) {
      best.index = k;
      best.angle = turn * 180.0f / kPi;
      best.quality = quality;
      best.point = point;
    }
  }
  return best;
}

// 跨槽底拟合只作为稳定的结构基准；左右交界附近再分别做多尺度局部底面
// 拟合。局部结果与基准方向、交点都一致时才有限修正左右交点，从而兼顾
// 焊接翘曲和缺口/焊道附近的稳定性。
// S06 三面结构候选：枚举左/右边界，拟合左侧壁、跨槽底面和右侧壁。
// 依次检查有效拟合、底面斜率、两侧交角、局部法向变化、转向和左右间距。
// 再进入 S07 的区域边缘/轮廓折角定位；有限数量候选按质量排序去重。
// 调试重点：哪一个 continue 门槛丢弃了候选，而不是先放大所有阈值。
struct SurfaceProposal {
  bool centralBead = false;
  SurfacePair pair;
  float cost = std::numeric_limits<float>::infinity();
};

enum class ProfileBottomShape { Flat, Wavy, CentralBead };

// 一组边界候选的局部状态，只在本次评估中使用，不跨截面共享。
// 各阶段依次完成结构验证、交点定位、区域生长和局部折角修正。
class SurfacePairBuilder {
public:
  SurfacePairBuilder(const std::vector<ProfilePoint> &profile,
                     const std::vector<ProfilePoint> &localizationProfile,
                     const DetectorConfig &config, int leftBoundary,
                     int rightBoundary, DebugSlice *debug,
                     const std::vector<ProfileCurvatureSample> &curvatures,
                     const std::vector<ProfileTurnSample> &turns,
                     const std::array<ObservedWall,2> &walls)
      : profile(profile), localizationProfile(localizationProfile),
        config(config), leftBoundary(leftBoundary),
        rightBoundary(rightBoundary), debug(debug), curvatures(curvatures), turns(turns), walls(walls),
        sideWindow(scaledProfileCount(config, config.fitHalfWindow, 5)),
        guard(scaledProfileCount(config, config.fitGuardPoints, 1, 1)),
        localWindow(std::max(
            3, std::min(sideWindow,
                        scaledProfileCount(config, config.localAngleWindow,
                                           3)))),
        pointCount(static_cast<int>(profile.size())) {}

  bool build(SurfaceProposal &result, bool retainOuterJunctions = false) {
    if (!fitSurfaces() || !validateLocalAngles() || !locateIntersections())
      return false;
    if (config.profileLocalizationMethod ==
        ProfileLocalizationMethod::CurveLineIntersection) {
      if (!locateCurveLineIntersections())
        return false;
      if (!retainOuterJunctions) locateBottomShape();
      scoreGeometry();
      recordSurfaceFits();
      proposal.pair.left = candidateFromIntersection(
          profile, leftIntersection, config.longitudinalAxis, feature,
          profileGeometryScale(config));
      proposal.pair.right = candidateFromIntersection(
          profile, rightIntersection, config.longitudinalAxis, feature,
          profileGeometryScale(config));
      if (bottomShape != ProfileBottomShape::CentralBead)
        attachSurfaceConstraints();
      else
        proposal.pair.left.floorDirection =
            proposal.pair.right.floorDirection = broadBottom.direction;
      auto attachCurveTangent = [](Candidate &candidate,
                                   const Eigen::Vector2f &point,
                                   const Eigen::Vector2f &tangent,
                                   float rms) {
        if (!tangent.allFinite() || tangent.squaredNorm() < .5f)
          return;
        const Eigen::Vector2f normal(-tangent.y(), tangent.x());
        candidate.hasWall = true;
        candidate.wallNormal = normal;
        candidate.wallOffset = normal.dot(point);
        candidate.wallRms = rms;
        candidate.trackedWall = true;
        candidate.trackedNormal = normal;
        candidate.trackedOffset = candidate.wallOffset;
      };
      if (bottomShape != ProfileBottomShape::CentralBead) {
        attachCurveTangent(proposal.pair.left, leftIntersection,
                           leftCurveTangent, leftSurface.rms);
        attachCurveTangent(proposal.pair.right, rightIntersection,
                           rightCurveTangent, rightSurface.rms);
      }
      finishProposal();
      result = std::move(proposal);
      return true;
    }
    scoreGeometry();
    recordSurfaceFits();
    locateBottomShape();
    if (bottomShape == ProfileBottomShape::CentralBead) {
      proposal.pair.left = candidateFromIntersection(
          profile, leftIntersection, config.longitudinalAxis, feature,
          profileGeometryScale(config));
      proposal.pair.right = candidateFromIntersection(
          profile, rightIntersection, config.longitudinalAxis, feature,
          profileGeometryScale(config));
      proposal.pair.left.floorDirection =
          proposal.pair.right.floorDirection = broadBottom.direction;
      finishProposal();
      result = std::move(proposal);
      return true;
    }
    growRegionCandidates();
    refineContourCandidates();
    if (config.useProfileCurvature &&
        (!proposal.pair.left.fromLocalCorner || !proposal.pair.right.fromLocalCorner)) {
      countGate("no_spline_curvature");
      return false;
    }
    attachSurfaceConstraints();
    finishProposal();
    result = std::move(proposal);
    return true;
  }

private:
  void countGate(const char *gate) {
    if (debug)
      ++debug->gates[gate];
  }
  bool fitSurfaces() {
    // 侧面使用交接点外侧的一段；底面从每个交界向槽内分别做多尺度
    // 局部拟合。guard 区域故意不参加拟合，因为这里通常是焊道冠部。
    leftSurface = fitOrthogonalLine(profile, leftBoundary - guard - sideWindow,
                                    leftBoundary - guard);
    const bool direct = config.profileLocalizationMethod ==
                        ProfileLocalizationMethod::CurveLineIntersection;
    const int bottomBegin = direct ? leftBoundary : leftBoundary + guard;
    const int bottomEnd = direct ? rightBoundary + 1 : rightBoundary - guard;
    broadBottom = direct
                      ? fitContinuousBottomLine(profile, bottomBegin, bottomEnd,
                                            config)
                      : fitOrthogonalLine(profile, bottomBegin, bottomEnd);
    rightSurface = fitOrthogonalLine(profile, rightBoundary + guard,
                                     rightBoundary + guard + sideWindow);
    const int availableBottomPoints = rightBoundary - leftBoundary - 2 * guard;
    const int maximumBottomWindow = std::min(
        availableBottomPoints / 2, std::max(localWindow, 2 * sideWindow));
    if (!leftSurface.valid || !broadBottom.valid || !rightSurface.valid ||
        (!direct && maximumBottomWindow < localWindow)) {
      countGate("invalid_fit");
      return false;
    }
    broadBottomSlope = std::fabs(broadBottom.direction.y()) /
                       std::max(1.0e-3f, std::fabs(broadBottom.direction.x()));
    leftAngle = acuteAngle(leftSurface, broadBottom);
    rightAngle = acuteAngle(rightSurface, broadBottom);
    if (broadBottomSlope > config.maxBottomSlope ||
        leftAngle < config.minSurfaceAngleDegrees * kPi / 180.0f ||
        rightAngle < config.minSurfaceAngleDegrees * kPi / 180.0f) {
      countGate("slope_or_broad_angle");
      return false;
    }
    if (direct)
      return true;
    leftBottom =
        fitAdaptiveBoundaryLine(profile, leftBoundary + guard, +1, localWindow,
                                maximumBottomWindow, config);
    rightBottom =
        fitAdaptiveBoundaryLine(profile, rightBoundary - guard, -1, localWindow,
                                maximumBottomWindow, config);
    const float maximumBottomAngleDifference =
        config.maxLocalBottomAngleDifferenceDegrees * kPi / 180.0f;
    leftLocalUsable = leftBottom.valid && acuteAngle(leftBottom, broadBottom) <=
                                              maximumBottomAngleDifference;
    rightLocalUsable =
        rightBottom.valid &&
        acuteAngle(rightBottom, broadBottom) <= maximumBottomAngleDifference;
    return true;
  }

  bool validateLocalAngles() {
    if (config.profileLocalizationMethod ==
        ProfileLocalizationMethod::CurveLineIntersection) {
      leftLocalAngle = leftAngle;
      rightLocalAngle = rightAngle;
      const float leftTurn =
          cross2D(broadBottom.direction, leftSurface.direction);
      const float rightTurn =
          cross2D(broadBottom.direction, rightSurface.direction);
      if (leftTurn * rightTurn >= 0 ||
          (config.interiorOpeningDirection != 0 &&
           leftTurn * config.interiorOpeningDirection >= 0)) {
        countGate("turn_direction");
        return false;
      }
      return true;
    }
    // 仅比较远处大面的方向会把平滑弧面两端的累计转角误判为折角。
    // 因此在每个理论交界的两侧各取一个短窗口：真正的板面交接会在很短
    // 距离内改变法向，而同一弧面的切向只会缓慢连续变化。
    const Line2D leftLocalOuter = fitOrthogonalLine(
        profile, leftBoundary - guard - localWindow, leftBoundary - guard);
    const Line2D leftLocalInner = fitOrthogonalLine(
        profile, leftBoundary + guard, leftBoundary + guard + localWindow);
    const Line2D rightLocalInner = fitOrthogonalLine(
        profile, rightBoundary - guard - localWindow, rightBoundary - guard);
    const Line2D rightLocalOuter = fitOrthogonalLine(
        profile, rightBoundary + guard, rightBoundary + guard + localWindow);
    if (!leftLocalOuter.valid || !leftLocalInner.valid ||
        !rightLocalInner.valid || !rightLocalOuter.valid) {
      countGate("invalid_local_fit");
      return false;
    }
    leftLocalAngle = acuteAngle(leftLocalOuter, leftLocalInner);
    rightLocalAngle = acuteAngle(rightLocalInner, rightLocalOuter);
    // 平均交角负责保留一侧轻微遮挡的截面，同时每一侧还必须达到一半
    // 下限，防止一侧噪声替另一侧完全平滑的弧面“代为通过”。
    if (std::min(leftLocalAngle, rightLocalAngle) <
            0.5f * (config.minLocalSurfaceAngleDegrees * kPi / 180.0f) ||
        0.5f * (leftLocalAngle + rightLocalAngle) <
            (config.minLocalSurfaceAngleDegrees * kPi / 180.0f)) {
      countGate("local_angle");
      return false;
    }

    // U 形结构的两个侧面必须位于底面方向的两侧。这个符号约束会排除
    // 同一块侧面上的两个纹理折角，但不依赖 U 形朝上还是朝下。
    const float leftTurn =
        cross2D(broadBottom.direction, leftSurface.direction);
    const float rightTurn =
        cross2D(broadBottom.direction, rightSurface.direction);
    if (leftTurn * rightTurn >= 0.0f) {
      countGate("turn_direction");
      return false;
    }
    if (config.interiorOpeningDirection != 0 &&
        leftTurn * config.interiorOpeningDirection >= 0.0f) {
      countGate("wrong_opening_direction");
      return false;
    }

    return true;
  }

  bool locateIntersections() {
    if (!intersectLines(leftSurface, broadBottom, &leftIntersection) ||
        !intersectLines(rightSurface, broadBottom, &rightIntersection)) {
      countGate("parallel_lines");
      return false;
    }
    // 跨槽底拟合提供稳定基准；只有局部底面与基准方向一致、且算出的
    // 交点没有远离基准时，才允许它做有限修正。这样既能跟随缓慢翘曲，
    // 又不会因焊道或缺口让交点突然跳到另一处。
    auto refineWithLocalBottom = [&](const Line2D &sideSurface,
                                     const Line2D &localBottom,
                                     bool localUsable,
                                     Eigen::Vector2f *intersection) {
      if (!localUsable || !intersection) {
        return;
      }
      Eigen::Vector2f localIntersection;
      if (!intersectLines(sideSurface, localBottom, &localIntersection)) {
        return;
      }
      const float difference = (localIntersection - *intersection).norm();
      if (difference > config.surfaceIntersectionTolerance) {
        return;
      }
      const float weight =
          0.15f * std::exp(-difference /
                           std::max(1.0f, config.surfaceIntersectionTolerance));
      *intersection =
          (1.0f - weight) * *intersection + weight * localIntersection;
    };
    refineWithLocalBottom(leftSurface, leftBottom, leftLocalUsable,
                          &leftIntersection);
    refineWithLocalBottom(rightSurface, rightBottom, rightLocalUsable,
                          &rightIntersection);
    leftBoundaryError =
        std::fabs(leftIntersection.x() -
                  profile[static_cast<std::size_t>(leftBoundary)].u);
    rightBoundaryError =
        std::fabs(rightIntersection.x() -
                  profile[static_cast<std::size_t>(rightBoundary)].u);
    if (leftBoundaryError > config.surfaceIntersectionTolerance ||
        rightBoundaryError > config.surfaceIntersectionTolerance ||
        rightIntersection.x() - leftIntersection.x() <
            (config.adaptiveSeamCount
                 ? std::max(1.f, 2 * profileGeometryScale(config))
                 : config.minSeamSeparation)) {
      countGate("intersection_or_separation");
      return false;
    }

    // 不比较绝对 Z。基准底面是当前站位的局部切线，点到切线的正交
    // 距离会随弧面姿态一起变化。
    leftClearance = std::fabs(broadBottom.normal.dot(leftSurface.center) -
                              broadBottom.offset);
    rightClearance = std::fabs(broadBottom.normal.dot(rightSurface.center) -
                               broadBottom.offset);
    if (leftClearance < config.minSurfaceClearance ||
        rightClearance < config.minSurfaceClearance) {
      countGate("clearance");
      return false;
    }

    return true;
  }

  bool locateCurveLineIntersections() {
    auto fitSide = [&](bool left, const Eigen::Vector2f &reference,
                       Eigen::Vector2f *intersection,
                       Eigen::Vector2f *intersectionTangent) {
      const int begin = left ? leftBoundary - guard - sideWindow
                             : rightBoundary - guard - localWindow;
      const int end = left ? leftBoundary + guard + localWindow
                           : rightBoundary + guard + sideWindow;
      if (begin < 0 || end > pointCount || end - begin < 5)
        return false;
      std::vector<Eigen::Vector2f> points;
      const float minimumU = profile[begin].u;
      const float maximumU = profile[end-1].u;
      for (const auto &point : localizationProfile)
        if (point.u >= minimumU && point.u <= maximumU)
          points.emplace_back(point.u,point.v);
      if (points.size() < 5) {
        points.clear();
        points.reserve(end-begin);
        for (int i=begin;i<end;++i)
          points.emplace_back(profile[i].u,profile[i].v);
      }
      const float sampleStep =
          std::max(.05f, profileGeometryScale(config) / 8);
      const auto curve = fitParametricProfileSpline(
          points, std::max(.2f, config.profileCurvatureRadius / 4),
          sampleStep);
      if (curve.size() < 2)
        return false;

      float bestScore = std::numeric_limits<float>::infinity();
      Eigen::Vector2f best = Eigen::Vector2f::Zero();
      Eigen::Vector2f bestTangent = Eigen::Vector2f::Zero();
      bool found = false;
      const float minimumCross = std::sin(
          .5f * config.minSurfaceAngleDegrees * kPi / 180.f);
      for (std::size_t i = 1; i < curve.size(); ++i) {
        if (!curve[i-1].valid || !curve[i].valid)
          continue;
        const float firstDistance =
            broadBottom.normal.dot(curve[i-1].point) - broadBottom.offset;
        const float secondDistance =
            broadBottom.normal.dot(curve[i].point) - broadBottom.offset;
        const float denominator = firstDistance - secondDistance;
        if (firstDistance * secondDistance > 0 ||
            std::fabs(denominator) < 1.e-7f)
          continue;
        const float fraction =
            std::max(0.f, std::min(1.f, firstDistance / denominator));
        const Eigen::Vector2f tangent =
            ((1-fraction)*curve[i-1].tangent + fraction*curve[i].tangent)
                .normalized();
        if (!tangent.allFinite() ||
            std::fabs(cross2D(tangent, broadBottom.direction)) < minimumCross)
          continue;
        const Eigen::Vector2f point =
            (1-fraction)*curve[i-1].point + fraction*curve[i].point;
        const float score = (point-reference).norm();
        if (score < bestScore) {
          bestScore = score;
          best = point;
          bestTangent = tangent;
          found = true;
        }
      }
      // A smooth fillet can approach the fitted floor from one side without
      // crossing it inside the measured side-window. Extend the local spline
      // tangent by a tightly bounded distance and intersect that continuation
      // with the floor line; this preserves the curve direction instead of
      // falling back to the coarse straight-wall fit.
      if (!found) {
        const float maximumExtension = std::min(
            config.surfaceIntersectionTolerance,
            2 * profileGeometryScale(config));
        for (const auto &sample:curve) {
          if (!sample.valid || !sample.tangent.allFinite())
            continue;
          const float denominator = broadBottom.normal.dot(sample.tangent);
          if (std::fabs(denominator) < minimumCross)
            continue;
          const float extension =
              (broadBottom.offset-broadBottom.normal.dot(sample.point)) /
              denominator;
          if (std::fabs(extension) > maximumExtension)
            continue;
          const Eigen::Vector2f point = sample.point+extension*sample.tangent;
          const float score = (point-reference).norm()+.25f*std::fabs(extension);
          if (score < bestScore) {
            bestScore=score;
            best=point;
            bestTangent=sample.tangent;
            found=true;
          }
        }
      }
      if (!found || bestScore > config.surfaceIntersectionTolerance)
        return false;
      *intersection = best;
      *intersectionTangent = bestTangent;
      return true;
    };

    const Eigen::Vector2f lineLeft = leftIntersection;
    const Eigen::Vector2f lineRight = rightIntersection;
    if (!fitSide(true, lineLeft, &leftIntersection, &leftCurveTangent) ||
        !fitSide(false, lineRight, &rightIntersection, &rightCurveTangent)) {
      countGate("curve_line_no_intersection");
      return false;
    }
    leftBoundaryError = std::fabs(
        leftIntersection.x() - profile[leftBoundary].u);
    rightBoundaryError = std::fabs(
        rightIntersection.x() - profile[rightBoundary].u);
    if (leftBoundaryError > config.surfaceIntersectionTolerance ||
        rightBoundaryError > config.surfaceIntersectionTolerance ||
        rightIntersection.x() - leftIntersection.x() <
            (config.profileLocalizationMethod ==
                         ProfileLocalizationMethod::CurveLineIntersection &&
                     config.adaptiveSeamCount
                 ? std::max(1.f, 2 * profileGeometryScale(config))
                 : config.minSeamSeparation)) {
      countGate("curve_line_geometry");
      return false;
    }
    return true;
  }

  void locateBottomShape() {
    // The side-curve/floor intersections define the complete groove floor.
    // Fit that interval once more and decide whether it contains one supported
    // central bead. Only that case changes the output locations to the two
    // valleys beside the bead; ordinary ripples retain the outer junctions.
    const float outerLeft = leftIntersection.x();
    const float outerRight = rightIntersection.x();
    const float floorWidth = outerRight - outerLeft;
    if (floorWidth <= std::max(2.f, config.centralBeadMinWidth * 1.5f))
      return;

    std::vector<Eigen::Vector2f> points;
    for (const auto &point : localizationProfile)
      if (point.u >= outerLeft && point.u <= outerRight)
        points.emplace_back(point.u, point.v);
    if (points.size() < 7)
      return;

    const float sampleStep =
        std::max(.05f, profileGeometryScale(config) / 8.f);
    const auto curve = fitParametricProfileSpline(
        points, std::max(.2f, config.profileCurvatureRadius / 4.f),
        sampleStep);
    if (curve.size() < 9)
      return;

    const Line2D robustShapeBottom =
        config.profileLocalizationMethod == ProfileLocalizationMethod::CurveLineIntersection
            ? broadBottom
            : fitRobustBottomLine(profile, leftBoundary, rightBoundary + 1, config);
    const Line2D &shapeBottom = robustShapeBottom.valid
                                    ? robustShapeBottom
                                    : broadBottom;

    const int opening = config.interiorOpeningDirection
                            ? config.interiorOpeningDirection
                            : (curve[curve.size() / 2].point.y() <=
                                       .5f * (curve.front().point.y() +
                                              curve.back().point.y())
                                   ? 1
                                   : -1);
    std::vector<float> height(curve.size(), 0.f);
    std::vector<int> usable;
    usable.reserve(curve.size());
    const float centralLeft = outerLeft + .2f * floorWidth;
    const float centralRight = outerRight - .2f * floorWidth;
    float squaredResidual = 0.f;
    for (std::size_t i = 0; i < curve.size(); ++i) {
      if (!curve[i].valid)
        continue;
      height[i] = opening *
                  (shapeBottom.normal.dot(curve[i].point) -
                   shapeBottom.offset);
      if (curve[i].point.x() >= centralLeft &&
          curve[i].point.x() <= centralRight) {
        usable.push_back(static_cast<int>(i));
        squaredResidual += height[i] * height[i];
      }
    }
    if (usable.size() < 7)
      return;

    const float rms = std::sqrt(squaredResidual / usable.size());
    const float wavyThreshold =
        std::max(.1f, .15f * config.centralBeadMinHeight);
    bottomShape = rms >= wavyThreshold ? ProfileBottomShape::Wavy
                                      : ProfileBottomShape::Flat;

    int peak = usable.front();
    for (int index : usable)
      if (height[index] > height[peak])
        peak = index;

    // Search for the nearest supported valleys on both sides. Limiting the
    // range to one quarter of the floor prevents the outer wall transitions
    // from replacing the valleys adjacent to a bead.
    const float valleySpan =
        std::max(2.f * config.centralBeadMinWidth, .25f * floorWidth);
    int leftMinimum = -1, rightMinimum = -1;
    for (int index : usable) {
      const float delta = curve[index].point.x() - curve[peak].point.x();
      if (delta < 0.f && delta >= -valleySpan &&
          (leftMinimum < 0 || height[index] < height[leftMinimum]))
        leftMinimum = index;
      if (delta > 0.f && delta <= valleySpan &&
          (rightMinimum < 0 || height[index] < height[rightMinimum]))
        rightMinimum = index;
    }
    if (leftMinimum < 0 || rightMinimum < 0)
      return;

    const float leftProminence = height[peak] - height[leftMinimum];
    const float rightProminence = height[peak] - height[rightMinimum];
    const float prominence = std::min(leftProminence, rightProminence);
    if (prominence < config.centralBeadMinHeight)
      return;

    // Measure the width at half of the smaller two-sided prominence. This is
    // insensitive to a tilted floor line and rejects isolated spike samples.
    const float halfHeight = height[peak] - .5f * prominence;
    int halfBegin = peak, halfEnd = peak;
    while (halfBegin > 0 && curve[halfBegin - 1].valid &&
           height[halfBegin - 1] >= halfHeight)
      --halfBegin;
    while (halfEnd + 1 < static_cast<int>(curve.size()) &&
           curve[halfEnd + 1].valid && height[halfEnd + 1] >= halfHeight)
      ++halfEnd;
    int supportedLobes = 0;
    bool insideLobe = false;
    for (int index : usable) {
      const bool high = height[index] >= halfHeight;
      if (high && !insideLobe)
        ++supportedLobes;
      insideLobe = high;
    }
    const float beadWidth =
        curve[halfEnd].point.x() - curve[halfBegin].point.x();
    if (supportedLobes != 1 ||
        beadWidth < config.centralBeadMinWidth ||
        beadWidth > .55f * floorWidth)
      return;

    auto continuousValleyMedian = [&](int minimum, int begin, int end,
                                      float localProminence) {
      const float tolerance = std::max(
          .05f, std::min(.25f * localProminence,
                         .5f * config.profileWallTolerance));
      int lowBegin = minimum, lowEnd = minimum;
      while (lowBegin > begin && curve[lowBegin - 1].valid &&
             height[lowBegin - 1] <= height[minimum] + tolerance)
        --lowBegin;
      while (lowEnd < end && curve[lowEnd + 1].valid &&
             height[lowEnd + 1] <= height[minimum] + tolerance)
        ++lowEnd;
      const float totalArc = [&] {
        float value = 0.f;
        for (int i = lowBegin + 1; i <= lowEnd; ++i)
          value += (curve[i].point - curve[i - 1].point).norm();
        return value;
      }();
      int selected = lowBegin;
      float traveled = 0.f;
      for (int i = lowBegin + 1; i <= lowEnd; ++i) {
        traveled += (curve[i].point - curve[i - 1].point).norm();
        if (traveled >= .5f * totalArc) {
          selected = i;
          break;
        }
      }
      return selected;
    };
    const int leftValley = continuousValleyMedian(
        leftMinimum, usable.front(), halfBegin, leftProminence);
    const int rightValley = continuousValleyMedian(
        rightMinimum, halfEnd, usable.back(), rightProminence);
    if (rightValley <= leftValley ||
        curve[rightValley].point.x() - curve[leftValley].point.x() <
            std::max(1.f, config.centralBeadMinWidth))
      return;

    leftIntersection = curve[leftValley].point;
    rightIntersection = curve[rightValley].point;
    leftCurveTangent = curve[leftValley].tangent;
    rightCurveTangent = curve[rightValley].tangent;
    bottomShape = ProfileBottomShape::CentralBead;
  }

  void scoreGeometry() {
    residualSum = leftSurface.rms + broadBottom.rms + rightSurface.rms;
    const float angleQuality =
        clamp01(std::min(leftAngle, rightAngle) / (0.5f * kPi));
    const float residualQuality =
        std::exp(-residualSum / std::max(0.2f, 3.0f * config.surfaceRmsScale));
    const float bottomQuality =
        std::exp(-broadBottomSlope / std::max(0.1f, config.maxBottomSlope));
    const float boundaryQuality =
        std::exp(-(leftBoundaryError + rightBoundaryError) /
                 std::max(1.0f, 2.0f * config.surfaceIntersectionTolerance));
    feature = clamp01(0.40f * angleQuality + 0.35f * residualQuality +
                      0.15f * bottomQuality + 0.10f * boundaryQuality);
  }

  void recordSurfaceFits() {
    if (debug) {
      proposal.pair.debug = std::make_shared<DebugPair>();
      DebugPair &d = *proposal.pair.debug;
      auto line = [](const Line2D &l, int begin, int end) {
        DebugLine value;
        value.u = l.center.x();
        value.v = l.center.y();
        value.du = l.direction.x();
        value.dv = l.direction.y();
        value.rms = l.rms;
        value.begin = begin;
        value.end = end;
        return value;
      };
      d.lines = {
          {line(leftSurface, leftBoundary - guard - sideWindow,
                leftBoundary - guard),
           line(broadBottom, leftBoundary + guard, rightBoundary - guard),
           line(rightSurface, rightBoundary + guard,
                rightBoundary + guard + sideWindow)}};
      d.boundaries = {{leftBoundary, rightBoundary}};
      d.angles = {{leftAngle * 180 / kPi, rightAngle * 180 / kPi}};
      d.localAngles = {
          {leftLocalAngle * 180 / kPi, rightLocalAngle * 180 / kPi}};
      d.clearance = {{leftClearance, rightClearance}};
    }
  }

  void growRegionCandidates() {
    if (config.useRegionGrowingEdges) {
      // 底面从候选槽底中央的稳定小段向两侧生长；侧面则从原有大面
      // 拟合段向交界处生长。允许范围略越过原理论边界，使边界位置由
      // 局部法向和距离决定，而不是被枚举下标预先固定。
      const int groundBegin = leftBoundary + guard;
      const int groundEnd = rightBoundary - guard;
      const int seedWidth = std::max(
          3, std::min(scaledProfileCount(config, config.regionGrowingWindow,
                                         3) *
                          (config.useContourRegionGrowing ? 2 : 1),
                      std::max(3, groundEnd - groundBegin)));
      const int groundMiddle = (groundBegin + groundEnd) / 2;
      int groundSeedBegin = groundMiddle - seedWidth / 2;
      int groundSeedEnd = groundSeedBegin + seedWidth;
      if (groundSeedBegin < groundBegin) {
        groundSeedBegin = groundBegin;
        groundSeedEnd = groundSeedBegin + seedWidth;
      }
      if (groundSeedEnd > groundEnd) {
        groundSeedEnd = groundEnd;
        groundSeedBegin = groundSeedEnd - seedWidth;
      }

      const int leftSeedBegin = leftBoundary - guard - sideWindow;
      const int leftSeedEnd = leftBoundary - guard;
      const int rightSeedBegin = rightBoundary + guard;
      const int rightSeedEnd = rightBoundary + guard + sideWindow;
      const ProfileRegion groundRegion = growProfileRegion(
          profile, groundSeedBegin, groundSeedEnd,
          std::max(0, leftBoundary - sideWindow),
          std::min(pointCount, rightBoundary + sideWindow + 1), config,
          config.useContourRegionGrowing);
      const ProfileRegion leftRegion = growProfileRegion(
          profile, leftSeedBegin, leftSeedEnd, leftSeedBegin,
          std::min(pointCount, leftBoundary + guard + localWindow), config);
      const ProfileRegion rightRegion =
          growProfileRegion(profile, rightSeedBegin, rightSeedEnd,
                            std::max(0, rightBoundary - guard - localWindow),
                            rightSeedEnd, config);
      const RegionEdgePair leftEdges =
          matchRegionEdges(profile, leftRegion, groundRegion, true, config);
      const RegionEdgePair rightEdges =
          matchRegionEdges(profile, rightRegion, groundRegion, false, config);
      if (debug) {
        DebugPair &d = *proposal.pair.debug;
        d.regionAttempted = true;
        d.seeds = {{leftSeedBegin, leftSeedEnd, groundSeedBegin, groundSeedEnd,
                    rightSeedBegin, rightSeedEnd}};
        d.regions = {{leftRegion.begin, leftRegion.end, groundRegion.begin,
                      groundRegion.end, rightRegion.begin, rightRegion.end}};
        d.edges = {{leftEdges.sidePoint.x(), leftEdges.sidePoint.y(),
                    leftEdges.groundPoint.x(), leftEdges.groundPoint.y(),
                    rightEdges.sidePoint.x(), rightEdges.sidePoint.y(),
                    rightEdges.groundPoint.x(), rightEdges.groundPoint.y()}};
        d.edgeGaps = {{leftEdges.distance, rightEdges.distance}};
      }
      const bool allowIndependentEdges = config.useContourRegionGrowing;
      if ((leftEdges.valid && rightEdges.valid) ||
          (allowIndependentEdges && (leftEdges.valid || rightEdges.valid))) {
        const float grownLeftU =
            leftEdges.valid
                ? 0.5f * (leftEdges.sidePoint.x() + leftEdges.groundPoint.x())
                : leftIntersection.x();
        const float grownRightU =
            rightEdges.valid
                ? 0.5f * (rightEdges.sidePoint.x() + rightEdges.groundPoint.x())
                : rightIntersection.x();
        if (grownRightU - grownLeftU >= config.minSeamSeparation) {
          const float regionQuality =
              (groundRegion.quality + leftRegion.quality +
               rightRegion.quality) /
              3.0f;
          const float gapSum = (leftEdges.valid ? leftEdges.distance : 0.0f) +
                               (rightEdges.valid ? rightEdges.distance : 0.0f);
          const int edgeCount = static_cast<int>(leftEdges.valid) +
                                static_cast<int>(rightEdges.valid);
          const float gapQuality = std::exp(
              -gapSum / std::max(0.5f, static_cast<float>(edgeCount) *
                                           config.regionGrowingMaxEdgeGap));
          const float grownFeature = clamp01(
              0.70f * feature + 0.20f * regionQuality + 0.10f * gapQuality);
          // 配对路径仍联合约束左右间距，但可靠的一侧不应被另一侧的孔洞
          // 或生长失败牵连。只对失败的一侧回退到理论交点。
          proposal.pair.left =
              leftEdges.valid
                  ? candidateFromRegionEdges(profile, leftEdges,
                                             config.longitudinalAxis,
                                             grownFeature)
                  : candidateFromIntersection(profile, leftIntersection,
                                              config.longitudinalAxis, feature,
                                              profileGeometryScale(config));
          proposal.pair.right =
              rightEdges.valid
                  ? candidateFromRegionEdges(profile, rightEdges,
                                             config.longitudinalAxis,
                                             grownFeature)
                  : candidateFromIntersection(profile, rightIntersection,
                                              config.longitudinalAxis, feature,
                                              profileGeometryScale(config));
          regionCost = 0.15f * gapSum - 0.35f * (regionQuality + gapQuality) +
                       0.25f * static_cast<float>(2 - edgeCount);
          const int uAxis = transverseAxes(config.longitudinalAxis)[0];
          usedRegionEdges = component(proposal.pair.right.point, uAxis) -
                                component(proposal.pair.left.point, uAxis) >=
                            config.minSeamSeparation;
          usedLeftEdges = usedRegionEdges && leftEdges.valid;
          usedRightEdges = usedRegionEdges && rightEdges.valid;
          if (!usedRegionEdges)
            regionCost = 0.0f;
        }
      }
    }
    // 区域生长在遮挡、孔洞或极稀疏截面中无法形成可靠边缘时，保留
    // 原理论交点作为降级路径，避免整段轨迹因少数坏截面中断。
    if (!usedRegionEdges) {
      proposal.pair.left = candidateFromIntersection(
          profile, leftIntersection, config.longitudinalAxis, feature,
          profileGeometryScale(config));
      proposal.pair.right = candidateFromIntersection(
          profile, rightIntersection, config.longitudinalAxis, feature,
          profileGeometryScale(config));
      regionCost += config.useRegionGrowingEdges ? 0.5f : 0.0f;
    }
  }

  void refineContourCandidates() {
    if (config.useProfileCurvature || (config.useRegionGrowingEdges && config.useContourRegionGrowing)) {
      const auto leftCorner =
          findLocalContourCorner(profile, leftSurface, broadBottom,
                                 leftIntersection, leftBoundary, true, config, curvatures, turns, walls);
      const auto rightCorner = findLocalContourCorner(
          profile, rightSurface, broadBottom, rightIntersection, rightBoundary,
          false, config, curvatures, turns, walls);
      Candidate localLeft = proposal.pair.left;
      Candidate localRight = proposal.pair.right;
      auto useCorner = [&](const LocalContourCorner &corner,
                           Candidate *candidate) {
        if (corner.index < 0)
          return;
        candidate->point =
            profile[static_cast<std::size_t>(corner.index)].point;
        const auto axes = transverseAxes(config.longitudinalAxis);
        candidate->point.data[axes[0]] = corner.point.x();
        candidate->point.data[axes[1]] = corner.point.y();
        candidate->feature = clamp01(0.85f * feature + 0.15f * corner.quality);
        candidate->profileCurvature = corner.curvature;
        if (corner.wall.valid) {
          candidate->hasWall=true;
          candidate->wallNormal=corner.wall.normal;
          candidate->wallOffset=corner.wall.offset;
          candidate->wallRms=corner.wall.rms;
          candidate->trackedWall=config.useProfileCurvature;
          candidate->trackedNormal=corner.wall.normal;
          candidate->trackedOffset=corner.wall.offset;
        }
        candidate->fromRegionEdges = false;
        candidate->fromLocalCorner = true;
      };
      useCorner(leftCorner, &localLeft);
      useCorner(rightCorner, &localRight);
      const int uAxis = transverseAxes(config.longitudinalAxis)[0];
      if (component(localRight.point, uAxis) -
              component(localLeft.point, uAxis) >=
          config.minSeamSeparation) {
        proposal.pair.left = localLeft;
        proposal.pair.right = localRight;
        usedLocalCorner = leftCorner.index >= 0 || rightCorner.index >= 0;
        cornerIndices = {{leftCorner.index, rightCorner.index}};
        if (leftCorner.index >= 0)
          usedLeftEdges = false;
        if (rightCorner.index >= 0)
          usedRightEdges = false;
        if (debug) {
          proposal.pair.debug->localCornerIndices = {
              {leftCorner.index, rightCorner.index}};
          proposal.pair.debug->localCornerAngles = {
              {leftCorner.angle, rightCorner.angle}};
          proposal.pair.debug->profileCurvatures = {{leftCorner.curvature, rightCorner.curvature}};
          for (int side=0;side<2;++side) {
            const auto &wall=side ? rightCorner.wall : leftCorner.wall;
            auto &d=proposal.pair.debug->candidateWalls[side];
            d.u=wall.center.x(); d.v=wall.center.y();
            d.du=wall.direction.x(); d.dv=wall.direction.y(); d.rms=wall.rms;
          }
        }
      }
    }
  }

  void attachSurfaceConstraints() {
    // 将侧壁约束与候选一起传递到路径末端，不能在平滑后丢掉几何归属。
    auto attachWall = [&](Candidate &candidate, const Line2D &surface,
                          int corner, bool left) {
      if (config.useProfileCurvature && candidate.fromLocalCorner && candidate.hasWall) return;
      if (!surface.valid)
        return;
      Line2D wall = surface;
      const int w = scaledProfileCount(
          config, std::min(5, config.localAngleWindow), 3);
      if (corner >= w && corner + w < pointCount) {
        const auto local =
            left ? fitOrthogonalLine(profile, corner - w, corner, true)
                 : fitOrthogonalLine(profile, corner + 1, corner + w + 1, true);
        if (local.valid && local.rms <= config.regionGrowingDistance &&
            acuteAngle(local, surface) <= 12.f * kPi / 180.f)
          wall = local;
      }
      if (wall.rms > config.regionGrowingDistance)
        return;
      candidate.wallNormal = wall.normal;
      candidate.wallOffset = wall.offset;
      candidate.wallRms = wall.rms;
      candidate.hasWall = true;
    };
    attachWall(proposal.pair.left, leftSurface, cornerIndices[0], true);
    attachWall(proposal.pair.right, rightSurface, cornerIndices[1], false);
    proposal.pair.left.floorDirection = proposal.pair.right.floorDirection =
        broadBottom.direction;
  }

  void finishProposal() {
    proposal.centralBead = bottomShape == ProfileBottomShape::CentralBead;
    proposal.pair.centralBead = proposal.centralBead;
    // 左右独立残差和底面方向一致性负责候选排序；全局优化随后再根据
    // 整条轨迹的连续性选择，不会让一侧的低质量观测污染另一侧置信度。
    proposal.cost =
        1.2f * residualSum + 0.08f * (leftBoundaryError + rightBoundaryError) -
        0.8f * (leftAngle + rightAngle) -
        0.02f * std::log1p(leftSurface.support * rightSurface.support) +
                       regionCost;
    countGate("accepted");
    if (debug) {
      proposal.pair.debug->proposalCost = proposal.cost;
      proposal.pair.debug->profileShape =
          bottomShape == ProfileBottomShape::CentralBead
              ? "central_bead"
              : bottomShape == ProfileBottomShape::Wavy ? "wavy_floor"
                                                        : "flat_floor";
      proposal.pair.debug->source =
          bottomShape == ProfileBottomShape::CentralBead
              ? "central_bead_valleys"
              : config.profileLocalizationMethod ==
                  ProfileLocalizationMethod::CurveLineIntersection
              ? (bottomShape == ProfileBottomShape::Wavy
                           ? "curve_line_wavy_floor"
                           : "curve_line_flat_floor")
              : usedLocalCorner
              ? (config.useProfileCurvature ? "spline_curvature" : "local_contour")
              : (usedLeftEdges && usedRightEdges
                     ? "region_edges"
                     : (usedRegionEdges ? "mixed_edges" : "intersection"));
      proposal.pair.debug->regionUsed = {{usedLeftEdges, usedRightEdges}};
    }
  }

  const std::vector<ProfilePoint> &profile;
  const std::vector<ProfilePoint> &localizationProfile;
  const DetectorConfig &config;
  const int leftBoundary, rightBoundary;
  DebugSlice *debug;
  const std::vector<ProfileCurvatureSample> &curvatures;
  const std::vector<ProfileTurnSample> &turns;
  const std::array<ObservedWall,2> &walls;
  const int sideWindow, guard, localWindow, pointCount;
  Line2D leftSurface, broadBottom, rightSurface, leftBottom, rightBottom;
  Eigen::Vector2f leftCurveTangent = Eigen::Vector2f::Zero();
  Eigen::Vector2f rightCurveTangent = Eigen::Vector2f::Zero();
  bool leftLocalUsable = false, rightLocalUsable = false;
  float broadBottomSlope = 0, leftAngle = 0, rightAngle = 0;
  float leftLocalAngle = 0, rightLocalAngle = 0;
  Eigen::Vector2f leftIntersection, rightIntersection;
  ProfileBottomShape bottomShape = ProfileBottomShape::Flat;
  float leftBoundaryError = 0, rightBoundaryError = 0;
  float leftClearance = 0, rightClearance = 0, residualSum = 0, feature = 0;
  float regionCost = 0;
  bool usedRegionEdges = false, usedLeftEdges = false, usedRightEdges = false;
  bool usedLocalCorner = false;
  std::array<int, 2> cornerIndices{{-1, -1}};
  SurfaceProposal proposal;
};

std::vector<SurfacePair>
findSurfacePairs(const std::vector<ProfilePoint> &profile,
                 const std::vector<ProfilePoint> &localizationProfile,
                 const DetectorConfig &config,
                 const std::array<HistoryPrediction, 2> &predictions,
                 std::array<ObservedWall,2> &wallTrack,
                 bool *hasGeometry, DebugSlice *debug, bool fallbackPass = false) {
  const auto previousWallTrack = wallTrack;
  const bool stableWallTurn = config.profileLocalizationMethod ==
      ProfileLocalizationMethod::CurvatureCandidates &&
      config.profileTurnLocation == ProfileTurnLocation::WallTurnCentroid;
  const bool multiCandidates = config.useProfileCurvature &&
      config.profileLocalizationMethod == ProfileLocalizationMethod::CurvatureCandidates &&
      config.curvatureSearch.multiCandidates && !stableWallTurn;
  std::vector<SurfaceProposal> proposals;
  std::vector<Eigen::Vector2f> bottomProfile;
  const auto &referenceProfile = localizationProfile.empty() ? profile : localizationProfile;
  for (const auto &p : referenceProfile) bottomProfile.emplace_back(p.u,p.v);
  // A shallow feature on a tilted workpiece must be measured relative to
  // independent outer supports, not to absolute height at the scan edge.
  Line2D shoulderBaseline;
  if (referenceProfile.size() >= 9) {
    const int count = std::max(3, int(referenceProfile.size()) / 5);
    std::vector<ProfilePoint> shoulders(referenceProfile.begin(), referenceProfile.begin()+count);
    shoulders.insert(shoulders.end(), referenceProfile.end()-count, referenceProfile.end());
    shoulderBaseline = fitOrthogonalLine(shoulders, 0, int(shoulders.size()), true);
    if (!shoulderBaseline.valid || shoulderBaseline.rms > config.profileWallTolerance ||
        std::fabs(shoulderBaseline.direction.y()) /
            std::max(1.e-3f, std::fabs(shoulderBaseline.direction.x())) > config.maxBottomSlope) {
      shoulderBaseline.valid = false;
    } else {
      const float slope = shoulderBaseline.direction.y() / shoulderBaseline.direction.x();
      for (auto &p : bottomProfile) p.y() -= slope * p.x();
    }
  }
  float bottomMidpoint = profileBottomMidpoint(bottomProfile,config.interiorOpeningDirection,
      std::max(config.profileWallTolerance,.5f*profileGeometryScale(config)));
  // Cache once per frame, rather than fitting again for every boundary pair.
  std::vector<ProfileCurvatureSample> curvatures(profile.size());
  std::vector<ProfileTurnSample> turns;
  // Short walls can retain visible lower turns without passing main-wall gates.
  std::vector<std::vector<ProfileCurvatureSample>> recoverySplines;
  std::array<ObservedWall,2> walls;
  if (config.useProfileCurvature) {
    std::size_t begin = 0;
    while (begin<localizationProfile.size()) {
      std::size_t end=begin+1;
      while (end<localizationProfile.size() &&
             localizationProfile[end].u-localizationProfile[end-1].u<=maximumProfileGap(config)) ++end;
      std::vector<Eigen::Vector2f> points;
      for (auto i=begin;i<end;++i) points.emplace_back(localizationProfile[i].u,localizationProfile[i].v);
      const float span=std::min(config.profileCurvatureRadius,3.f*profileGeometryScale(config));
      const float step=std::max(.05f,profileGeometryScale(config)/8);
      const auto fitted=fitParametricProfileSpline(points,config.profileCurvatureRadius/4,step);
      observeWalls(fitted,config,walls,bottomMidpoint);
      std::vector<Eigen::Vector2f> curve;
      for (const auto &sample:fitted) {
        curve.push_back(sample.point);
        if (debug) {
          const auto axes=transverseAxes(config.longitudinalAxis);
          PointT p=profile.front().point;
          p.data[axes[0]]=sample.point.x(); p.data[axes[1]]=sample.point.y();
          debug->fittedProfile.push_back(p); debug->fittedCurvature.push_back(sample.curvature);
        }
      }
      // Tangents only validate which side is the wall. Peak positions and
      // magnitudes come exclusively from analytic B-spline derivatives.
      std::vector<bool> peaks(fitted.size(),false);
      for (std::size_t k=1;k+1<fitted.size();++k)
        peaks[k]=fitted[k].valid &&
            std::fabs(fitted[k].curvature)>=std::fabs(fitted[k-1].curvature) &&
            std::fabs(fitted[k].curvature)>=std::fabs(fitted[k+1].curvature);
      const bool regionMode = config.profileTurnLocation != ProfileTurnLocation::Peak;
      if (stableWallTurn) {
        std::fill(peaks.begin(),peaks.end(),false);
        for (int side=0;side<2;++side) if (walls[side].line.valid) {
          const auto &wall=walls[side];
          const int index=wallSupportedTurnCentroid(fitted,wall.line.center,
              wall.line.normal,wall.line.offset,config.profileWallTolerance,
              std::min(config.profileCurvatureRadius,4*profileGeometryScale(config)),
              config.interiorOpeningDirection,bottomMidpoint,side);
          if (index>=0) peaks[index]=true;
        }
      } else if (regionMode) {
        std::fill(peaks.begin(), peaks.end(), false);
        const auto regions = locateProfileTurnRegions(fitted,
            config.profileTurnLocation == ProfileTurnLocation::TurnMedian,
            config.profileTurnRegionFraction, config.profileTurnRegionMaxLength);
        for (auto region : regions) {
          if (config.profileTurnLocation == ProfileTurnLocation::FloorSideBoundary) {
            const int floorBoundary = profileTurnFloorBoundary(
                fitted, region, config.interiorOpeningDirection);
            const int side = fitted[floorBoundary].point.x() < bottomMidpoint ? 0 : 1;
            const auto &wall = walls[side];
            region.representative = floorBoundary;
            if (wall.line.valid) {
              const int stepToWall = floorBoundary == region.begin ? 1 : -1;
              const int beyondRegion = floorBoundary == region.begin
                                           ? region.end + 1
                                           : region.begin - 1;
              // Start at the floor-side endpoint and move only as far toward
              // the curvature peak as needed to retain measured wall support.
              // Raising profileWallTolerance therefore permits a deeper point
              // without changing the B-spline itself.
              for (int index = floorBoundary; index != beyondRegion;
                   index += stepToWall) {
                const auto &candidate = fitted[index].point;
                const float wallDistance = std::fabs(
                    wall.line.normal.dot(candidate) - wall.line.offset);
                const float footDistance = (candidate - wall.foot).norm();
                if (wallDistance <= config.profileWallTolerance &&
                    footDistance <= std::min(config.profileCurvatureRadius,
                                             5 * profileGeometryScale(config))) {
                  region.representative = index;
                  break;
                }
              }
            }
          }
          if (config.profileTurnLocation != ProfileTurnLocation::FloorSideBoundary) {
            const int side = fitted[region.representative].point.x() < bottomMidpoint ? 0 : 1;
            const auto &wall = walls[side];
            if (wall.line.valid) {
              const int supported = supportedProfileTurnRepresentative(fitted, region,
                  wall.line.normal, wall.line.offset, wall.foot,
                  config.profileWallTolerance,
                  std::min(config.profileCurvatureRadius, 5*profileGeometryScale(config)));
              if (supported >= 0) {
                if (debug && supported != region.representative)
                  ++debug->gates["turn_representative_inside_wall_support"];
                region.representative = supported;
              }
            }
          }
          peaks[region.representative] = true;
          if (debug) {
            const auto axes = transverseAxes(config.longitudinalAxis);
            auto point = [&](int index) {
              PointT p = profile.front().point;
              p.data[axes[0]] = fitted[index].point.x(); p.data[axes[1]] = fitted[index].point.y();
              return p;
            };
            DebugTurnRegion d;
            d.begin=point(region.begin); d.end=point(region.end);
            d.peak=point(region.peak); d.representative=point(region.representative);
            d.arcLength=region.arcLength; d.turningAngle=region.turningAngle;
            debug->turnRegions.push_back(d);
          }
        }
      }
      std::vector<unsigned> kinds(fitted.size(), 0);
      for (std::size_t k=0;k<peaks.size();++k) if (peaks[k]) kinds[k]=1;
      if (multiCandidates) {
        const auto extra = supplementalProfileCandidates(fitted,
            config.curvatureSearch.windowLength,
            config.curvatureSearch.valleyProminence, config.interiorOpeningDirection);
        for (std::size_t k=0;k<peaks.size();++k) {
          kinds[k] |= extra[k]; peaks[k] = kinds[k] != 0;
        }
      }
      auto measured=measureProfileTurns(curve,span,step,
          config.useWallCandidatePruning ? &peaks : nullptr);
      std::vector<bool> used(fitted.size(),false);
      for (auto &turn:measured) {
        const int k=turn.index;
        if (k<1 || k+1>=int(fitted.size()) || used[k] || !fitted[k].valid) continue;
        const float magnitude=std::fabs(fitted[k].curvature);
        if ((regionMode || multiCandidates ? !peaks[k] :
             (magnitude<std::fabs(fitted[k-1].curvature) || magnitude<std::fabs(fitted[k+1].curvature))) ||
            fitted[k].curvature*turn.angle<=0) continue;
        used[k]=true;
        turn.candidateKinds=kinds[k];
        turn.point=fitted[k].point;
        turn.curvature=fitted[k].curvature;
        float supportDistance=std::numeric_limits<float>::infinity();
        for (std::size_t j=1;j<points.size();++j) {
          const Eigen::Vector2f d=points[j]-points[j-1];
          const float t=std::max(0.f,std::min(1.f,d.dot(turn.point-points[j-1])/std::max(1.e-9f,d.squaredNorm())));
          supportDistance=std::min(supportDistance,(turn.point-points[j-1]-t*d).norm());
        }
        if (supportDistance>config.profileWallTolerance) continue;
        auto found=std::lower_bound(profile.begin(),profile.end(),turn.point.x(),
            [](const ProfilePoint &p,float u){return p.u<u;});
        int index=std::min(int(profile.size())-1,int(found-profile.begin()));
        if (index>0 && std::fabs(profile[index-1].u-turn.point.x())<std::fabs(profile[index].u-turn.point.x())) --index;
        turn.index=index;
        turns.push_back(turn);
      }
      recoverySplines.push_back(fitted);
      begin=end;
    }
  }
  // Track the observed wall itself, before candidate scoring or history.
  const float station=coordinate(profile.front().point,config.longitudinalAxis);
  for (int side=0;side<2 && config.useProfileCurvature;++side) {
    auto &wall=walls[side]; const auto &previous=wallTrack[side];
    if (!wall.line.valid) continue;
    const float gap=std::fabs(station-previous.station);
    if (previous.line.valid && gap<config.profileCurvatureRadius &&
        (wall.foot-previous.foot).norm()>std::max(3*profileGeometryScale(config),2*gap)) {
      wall.line.valid=false;
      if (debug) ++debug->gates["wall_identity_jump"];
    } else { wall.station=station; wallTrack[side]=wall; }
  }
  // Refine the dividing line from independently observed wall feet only when
  // they surround the measured bottom. Never derive it from the candidate
  // pair being tested (which would make a side check tautological).
  if (walls[0].line.valid && walls[1].line.valid &&
      profilePairStraddlesBottom(walls[0].foot.x(),walls[1].foot.x(),bottomMidpoint))
    bottomMidpoint=.5f*(walls[0].foot.x()+walls[1].foot.x());
  const int transverseAxis=transverseAxes(config.longitudinalAxis)[0];
  auto acceptBottomSides = [&](SurfacePair &pair) {
    if (!profilePairStraddlesBottom(component(pair.left.point,transverseAxis),
                                   component(pair.right.point,transverseAxis),bottomMidpoint)) {
      if (debug) ++debug->gates["candidate_crosses_bottom_midpoint"];
      return false;
    }
    if (pair.debug) {
      pair.debug->hasBottomMidpoint=true;
      pair.debug->bottomMidpointU=bottomMidpoint;
    }
    pair.bottomMidpointU=bottomMidpoint;
    return true;
  };
  const int sideWindow =
      scaledProfileCount(config, config.fitHalfWindow, 5);
  const int guard =
      scaledProfileCount(config, config.fitGuardPoints, 1, 1);
  const int pointCount = static_cast<int>(profile.size());
  const bool directCurveLine =
      config.profileLocalizationMethod ==
      ProfileLocalizationMethod::CurveLineIntersection;
  const int minimumBoundaryGap =
      directCurveLine && config.adaptiveSeamCount
          ? std::max(3, int(std::floor(
                            1.5f * profileGeometryScale(config) /
                            std::max(.05f, config.profileBinSize))))
          : std::max(8, static_cast<int>(0.70f * config.minSeamSeparation /
                                         std::max(0.05f, config.profileBinSize)));
  const bool enoughForSurfaceModel =
      pointCount >= 2 * (sideWindow + guard) + minimumBoundaryGap;
  // Coarse profiles already reduce raw samples to 1.2 mm bins; fine spline
  // peaks retain 0.3 mm input resolution. Only fit broad surface hypotheses
  // whose boundary can reach a supported peak at the independently found wall.
  // These are necessary conditions of findLocalContourCorner, using the SAME
  // index radius, not an extra geometric threshold or an approximate crop.
  const bool prune=config.useProfileCurvature && config.useWallCandidatePruning;
  std::array<std::vector<bool>,2> eligible{{
      std::vector<bool>(pointCount,!prune && !directCurveLine),
      std::vector<bool>(pointCount,!prune && !directCurveLine)}};
  if (debug && config.useProfileCurvature && turns.empty())
    ++debug->gates["curvature_turns_missing"];
  if (prune || directCurveLine) {
    const int radius=contourBoundaryRadius(config);
    for (int side=0;side<2;++side) {
      if (directCurveLine && walls[side].line.valid) {
        auto nearest = std::lower_bound(
            profile.begin(), profile.end(), walls[side].foot.x(),
            [](const ProfilePoint &point, float u) { return point.u < u; });
        int index = std::min(pointCount - 1, int(nearest - profile.begin()));
        if (index > 0 && std::fabs(profile[index-1].u-walls[side].foot.x()) <
                             std::fabs(profile[index].u-walls[side].foot.x()))
          --index;
        const int directRadius = std::max(
            1, int(std::ceil(.5f * profileGeometryScale(config) /
                             std::max(.05f, config.profileBinSize))));
        for (int i=std::max(0,index-directRadius);
             i<=std::min(pointCount-1,index+directRadius);++i)
          eligible[side][i]=true;
        continue;
      }
      const std::string prefix=side ? "right_" : "left_";
      std::size_t supportedTurns=0;
      if (debug && !walls[side].line.valid)
        ++debug->gates[prefix+"observed_wall_missing"];
      if (debug && walls[side].line.valid &&
          walls[side].line.rms>config.profileWallTolerance)
        ++debug->gates[prefix+"observed_wall_rms_too_large"];
      for (const auto &turn:turns) {
        if (debug && walls[side].line.valid) {
          if ((turn.point-walls[side].foot).norm()>
              std::min(config.profileCurvatureRadius,
                       5*profileGeometryScale(config)))
            ++debug->gates[prefix+"turn_far_from_wall_foot"];
          else if (std::fabs(walls[side].line.normal.dot(turn.point)-
                             walls[side].line.offset)>
                   config.profileWallTolerance)
            ++debug->gates[prefix+"turn_outside_wall_tolerance"];
          else if (std::fabs(turn.angle)<
                   config.minLocalSurfaceAngleDegrees*kPi/180.f)
            ++debug->gates[prefix+"turn_angle_too_small"];
        }
        if (!nearObservedWall(turn,walls[side],config) ||
            std::fabs(turn.angle)<config.minLocalSurfaceAngleDegrees*kPi/180.f ||
            walls[side].line.rms>config.profileWallTolerance) continue;
        ++supportedTurns;
        for (int i=std::max(0,turn.index-radius);i<=std::min(pointCount-1,turn.index+radius);++i)
          eligible[side][i]=true;
      }
      if (debug && config.useProfileCurvature && supportedTurns==0)
        ++debug->gates[prefix+"no_supported_wall_turn"];
    }
  }
  if (!enoughForSurfaceModel && debug)
    ++debug->gates["profile_insufficient"];
  if (enoughForSurfaceModel) {
    for (int leftBoundary = sideWindow + guard;
         leftBoundary + minimumBoundaryGap + sideWindow + guard < pointCount;
         ++leftBoundary) {
      if (!eligible[0][leftBoundary]) {
        if (debug) debug->gates["wall_roi_pruned"]+=
            pointCount-guard-sideWindow-leftBoundary-minimumBoundaryGap+1;
        continue;
      }
      for (int rightBoundary = leftBoundary + minimumBoundaryGap;
           rightBoundary + guard + sideWindow <= pointCount; ++rightBoundary) {
        if (!eligible[1][rightBoundary]) {
          if (debug) ++debug->gates["wall_roi_pruned"];
          continue;
        }
        if (debug)
          ++debug->gates["enumerated"];
        SurfaceProposal proposal;
        if (SurfacePairBuilder(profile, localizationProfile, config,
                               leftBoundary, rightBoundary,
                               debug, curvatures, turns, walls)
                .build(proposal) && acceptBottomSides(proposal.pair)) {
          // Near the bead threshold, preserve the independently fitted outer
          // intersections as an alternative. A single shape-label change must
          // not force both tracks onto interior valleys. Both locations share
          // the same measured surface evidence: a repeated shape bonus would
          // accumulate over frames and eventually pay for an artificial jump.
          // Keep only the proposal-order tie break; path costs use geometry
          // and continuity without a recurring preference for the inner branch.
          if (directCurveLine && proposal.centralBead) {
            SurfaceProposal outer;
            if (SurfacePairBuilder(profile,localizationProfile,config,
                    leftBoundary,rightBoundary,debug,curvatures,turns,walls)
                    .build(outer,true) && acceptBottomSides(outer.pair)) {
              outer.cost+=.06f;
              if (outer.pair.debug) {
                outer.pair.debug->source="curve_line_outer_alternative";
                outer.pair.debug->profileShape="ambiguous_bead";
                outer.pair.debug->proposalCost=outer.cost;
              }
              proposals.push_back(std::move(outer));
            }
          }
          proposals.push_back(std::move(proposal));
        }
      }
    }
  }
  // Rough floors and narrow grooves may not support a three-plane model.
  // Independently established walls can still validate two measured peaks.
  const bool hasBeadProposal = std::any_of(proposals.begin(), proposals.end(),
      [](const SurfaceProposal &p) { return p.centralBead; });
  const std::string supportedShape = debug && !proposals.empty() && proposals.front().pair.debug
      ? proposals.front().pair.debug->profileShape : "unknown";
  const bool broadSurfaceSupported = !proposals.empty();
  const std::size_t broadProposalCount = proposals.size();
  if ((proposals.empty() || (multiCandidates && !hasBeadProposal)) && config.useProfileCurvature && !directCurveLine) {
    const auto axes=transverseAxes(config.longitudinalAxis);
    for (const auto &left:turns) {
      if (!nearObservedWall(left,walls[0],config)) continue;
      const int opening=config.interiorOpeningDirection ? config.interiorOpeningDirection : (left.angle>0 ? 1 : -1);
      if (left.angle*opening<config.minLocalSurfaceAngleDegrees*kPi/180.f ||
          walls[0].line.direction.y()*opening>-.35f ||
          walls[0].line.rms>config.profileWallTolerance) continue;
      for (const auto &right:turns) {
        if (!nearObservedWall(right,walls[1],config)) continue;
        const float width=right.point.x()-left.point.x();
        if (width<config.minSeamSeparation ||
            right.angle*opening<config.minLocalSurfaceAngleDegrees*kPi/180.f ||
            walls[1].line.direction.y()*opening<.35f ||
            walls[1].line.rms>config.profileWallTolerance) continue;
        SurfaceProposal proposal;
        if (debug) proposal.pair.debug=std::make_shared<DebugPair>();
        bool supported=true;
        for (int side=0;side<2;++side) {
          const auto &turn=side ? right : left;
          const auto &direction=walls[side].line.direction;
          const auto &center=walls[side].line.center;
          const Eigen::Vector2f normal(-direction.y(),direction.x());
          if (std::fabs(normal.dot(turn.point-center))>config.profileWallTolerance) { supported=false; break; }
          Candidate &candidate=side ? proposal.pair.right : proposal.pair.left;
          candidate.point=profile[turn.index].point;
          candidate.point.data[axes[0]]=turn.point.x(); candidate.point.data[axes[1]]=turn.point.y();
          candidate.hasWall=true; candidate.wallNormal=normal; candidate.wallOffset=normal.dot(center);
          candidate.wallRms=walls[side].line.rms;
          candidate.trackedWall=true; candidate.trackedNormal=walls[side].line.normal;
          candidate.trackedOffset=walls[side].line.offset;
          candidate.fromLocalCorner=true; candidate.profileCurvature=turn.curvature;
          candidate.feature=clamp01(.5f+.5f*std::fabs(turn.angle)/(0.5f*kPi));
          candidate.floorDirection=(right.point-left.point).normalized();
          if (debug) {
            auto &d=*proposal.pair.debug;
            d.profileShape=supportedShape;
            d.source=!broadSurfaceSupported && width<config.doubleSeamWidth ? "narrow_spline_curvature" : "spline_curvature";
            d.profileCurvatures[side]=turn.curvature;
            d.candidateKinds[side]=turn.candidateKinds;
            auto &line=d.candidateWalls[side]; line.u=center.x(); line.v=center.y();
            line.du=direction.x(); line.dv=direction.y(); line.rms=candidate.wallRms;
          }
        }
        if (supported && broadSurfaceSupported) {
          // Independent wall support defines the same search domain whether
          // a coarse three-surface hypothesis succeeds or fails this frame.
          // Its proximity is only a soft preference, never a moving crop.
          float distance=std::numeric_limits<float>::infinity();
          for (std::size_t i=0;i<broadProposalCount;++i) {
            const float d=std::max(
                transverseDistance(proposal.pair.left.point,proposals[i].pair.left.point,config.longitudinalAxis),
                transverseDistance(proposal.pair.right.point,proposals[i].pair.right.point,config.longitudinalAxis));
            distance=std::min(distance,d);
          }
          proposal.cost=-proposal.pair.left.feature-proposal.pair.right.feature+
              .05f*std::min(2.f,distance/config.curvatureSearch.windowLength);
        } else if (supported) {
          proposal.cost=-proposal.pair.left.feature-proposal.pair.right.feature;
        }
        if (supported && acceptBottomSides(proposal.pair)) {
          if (proposal.pair.debug) proposal.pair.debug->proposalCost=proposal.cost;
          proposals.push_back(std::move(proposal));
        }
      }
    }
  }
  // At the beginning of a run the seam can first appear as a shallow notch,
  // a low central bead, or a broad ripple.  There is then a real measured
  // feature but not yet enough vertical extent to establish two stable walls.
  // Recover only such central, baseline-supported features after the normal
  // wall-constrained search failed.  These deliberately receive a lower
  // feature score, so the global path uses them only when they form a smooth
  // continuation into the independently observed wall segment.
  if (proposals.empty() && config.useProfileCurvature && !directCurveLine &&
      localizationProfile.size() >= 9) {
    std::vector<ProfileTurnSample> recoveryTurns;
    for (const auto &fitted:recoverySplines) {
      std::vector<Eigen::Vector2f> curve;
      for (const auto &p:fitted) curve.push_back(p.point);
      std::vector<bool> mask(fitted.size(),false), used(fitted.size(),false);
      for (const auto &region:locateProfileTurnRegions(fitted,true,
              config.profileTurnRegionFraction,config.profileTurnRegionMaxLength))
        mask[region.representative]=true;
      for (auto turn:measureProfileTurns(curve,
              std::min(config.profileCurvatureRadius,3.f*profileGeometryScale(config)),
              std::max(.05f,profileGeometryScale(config)/8),&mask)) {
        const int k=turn.index;
        if (k<1 || k+1>=int(fitted.size()) || used[k] || !fitted[k].valid ||
            fitted[k].curvature*turn.angle<=0) continue;
        used[k]=true;
        turn.point=fitted[k].point; turn.curvature=fitted[k].curvature;
        float distance=std::numeric_limits<float>::infinity();
        for (std::size_t j=1;j<localizationProfile.size();++j) {
          const Eigen::Vector2f a(localizationProfile[j-1].u,localizationProfile[j-1].v);
          const Eigen::Vector2f b(localizationProfile[j].u,localizationProfile[j].v);
          if ((b-a).norm()>maximumProfileGap(config)) continue;
          const Eigen::Vector2f d=b-a;
          const float t=std::max(0.f,std::min(1.f,d.dot(turn.point-a)/std::max(1.e-9f,d.squaredNorm())));
          distance=std::min(distance,(turn.point-a-t*d).norm());
        }
        if (distance<=config.profileWallTolerance) recoveryTurns.push_back(turn);
      }
    }
    std::sort(recoveryTurns.begin(),recoveryTurns.end(),
        [](const ProfileTurnSample &a,const ProfileTurnSample &b){return a.point.x()<b.point.x();});
    std::vector<Line2D> baselines{fitRobustBottomLine(
        localizationProfile, 0, int(localizationProfile.size()), config)};
    // A shallow trough can be explained by a diagonal through one shoulder
    // and half of the floor. If that creates only same-side features, retry
    // against the two independently measured outer support strips.
    if (shoulderBaseline.valid)
      baselines.push_back(shoulderBaseline);
    for (const auto &baseline:baselines) {
    if (!proposals.empty()) break;
    if (baseline.valid &&
        std::fabs(baseline.direction.y()) /
                std::max(1.e-3f, std::fabs(baseline.direction.x())) <=
            config.maxBottomSlope) {
      const float scale = profileGeometryScale(config);
      const float residualThreshold = std::max(.12f, .15f * scale);
      const float minimumAmplitude =
          std::max(.25f, .3f * config.centralBeadMinHeight);
      const float minimumWidth = std::max(.5f, .5f * scale);
      const float maximumWidth = std::max(
          2.f * config.profileCurvatureRadius,
          1.5f * config.doubleSeamWidth);
      const float profileCenter =
          .5f * (localizationProfile.front().u +
                 localizationProfile.back().u);
      const float centerTolerance = std::min(
          config.profileCurvatureRadius,
          .25f * (localizationProfile.back().u -
                  localizationProfile.front().u));
      std::vector<float> residual(localizationProfile.size(), 0.f);
      for (std::size_t i = 0; i < localizationProfile.size(); ++i)
        residual[i] = baseline.normal.dot(Eigen::Vector2f(
                          localizationProfile[i].u,
                          localizationProfile[i].v)) -
                      baseline.offset;

      const std::size_t margin = std::max<std::size_t>(
          2, localizationProfile.size() / 10);
      for (std::size_t seed = margin;
           seed + margin < localizationProfile.size();) {
        if (std::fabs(residual[seed]) < residualThreshold) {
          ++seed;
          continue;
        }
        const float sign = residual[seed] >= 0.f ? 1.f : -1.f;
        std::size_t begin = seed, end = seed;
        while (begin > margin &&
               sign * residual[begin - 1] >= residualThreshold)
          --begin;
        while (end + margin + 1 < localizationProfile.size() &&
               sign * residual[end + 1] >= residualThreshold)
          ++end;
        seed = end + 1;
        // A low groove viewed with the wrong opening must not be relabelled
        // as a raised bead by the emergency no-wall recovery branch.
        if (fallbackPass && config.interiorOpeningDirection &&
            sign * config.interiorOpeningDirection > 0) continue;

        std::size_t peak = begin;
        for (std::size_t i = begin + 1; i <= end; ++i)
          if (sign * residual[i] > sign * residual[peak]) peak = i;
        const float amplitude = sign * residual[peak];
        const std::size_t outerBegin = begin > margin ? begin - 1 : begin;
        const std::size_t outerEnd =
            end + margin + 1 < localizationProfile.size() ? end + 1 : end;
        const float leftU = localizationProfile[outerBegin].u;
        const float rightU = localizationProfile[outerEnd].u;
        const float width = rightU - leftU;
        const float featureCenter = .5f * (leftU + rightU);
        if (amplitude < minimumAmplitude || width < minimumWidth ||
            std::fabs(featureCenter - profileCenter) > centerTolerance)
          continue;
        // A one-sided residual is not a complete groove. A turn's search
        // radius must not let it borrow a point just beyond the dividing line.
        if (!profilePairStraddlesBottom(leftU,rightU,bottomMidpoint)) {
          if (debug) ++debug->gates["candidate_crosses_bottom_midpoint"];
          continue;
        }

        // Interval ends bound the feature; they are not observed lower turns.
        Eigen::Vector2f leftPoint=Eigen::Vector2f::Zero();
        Eigen::Vector2f rightPoint=Eigen::Vector2f::Zero();
        std::array<float, 2> selectedCurvature{{0.f, 0.f}};
        const int opening = config.interiorOpeningDirection
                                ? config.interiorOpeningDirection
                                : 1;
        const float minimumTurn = .25f *
            config.minLocalSurfaceAngleDegrees * kPi / 180.f;
        const bool groove=sign*opening<0;
        auto featureTurns=recoveryTurns;
        if (groove) {
          featureTurns.clear();
          for (const auto &fitted:recoverySplines) {
            std::vector<double> arc(fitted.size(),0);
            std::vector<Eigen::Vector2f> curve;
            for (std::size_t i=0;i<fitted.size();++i) {
              curve.push_back(fitted[i].point);
              if (i) arc[i]=arc[i-1]+(fitted[i].point-fitted[i-1].point).norm();
            }
            std::vector<bool> mask(fitted.size(),false),used(fitted.size(),false);
            for (int side=0;side<2;++side) {
              int anchor=-1;
              auto eligible=[&](int i) {
                const auto &p=fitted[i];
                return p.valid && p.point.x()>=leftU-scale && p.point.x()<=rightU+scale &&
                    (side ? p.point.x()>bottomMidpoint : p.point.x()<bottomMidpoint) &&
                    opening*(baseline.normal.dot(p.point)-baseline.offset)<=-.55f*amplitude;
              };
              // First entry into the measured lower band locates each short
              // wall independently of which curvature peak is strongest.
              for (int i=1;i+1<int(fitted.size());++i) if (eligible(i)) {
                anchor=i;if (!side) break;
              }
              if (anchor<0) continue;
              const float radius=std::min(config.profileCurvatureRadius,3*scale);
              double mass=0,moment=0;
              std::vector<int> supported;
              for (int i=1;i+1<int(fitted.size());++i) {
                const auto &p=fitted[i];
                const float distance=(p.point-fitted[anchor].point).norm();
                if (!eligible(i) || p.curvature*opening<=0 || distance>=radius) continue;
                const double taper=1-distance*distance/(radius*radius);
                const double weight=p.curvature*opening*.5*(arc[i+1]-arc[i-1])*taper*taper;
                mass+=weight;moment+=weight*arc[i];supported.push_back(i);
              }
              if (mass<minimumTurn) continue;
              int best=-1;double error=std::numeric_limits<double>::infinity();
              for (int i:supported) if (std::fabs(arc[i]-moment/mass)<error) {
                error=std::fabs(arc[i]-moment/mass);best=i;
              }
              if (best>=0) mask[best]=true;
            }
            for (auto turn:measureProfileTurns(curve,std::min(config.profileCurvatureRadius,3*scale),
                                                std::max(.05f,scale/8),&mask)) {
              const int k=turn.index;
              if (k<1 || k+1>=int(fitted.size()) || used[k] || !mask[k]) continue;
              used[k]=true;turn.point=fitted[k].point;turn.curvature=fitted[k].curvature;
              float support=std::numeric_limits<float>::infinity();
              for (std::size_t j=1;j<localizationProfile.size();++j) {
                const Eigen::Vector2f a(localizationProfile[j-1].u,localizationProfile[j-1].v);
                const Eigen::Vector2f b(localizationProfile[j].u,localizationProfile[j].v),d=b-a;
                if (d.norm()>maximumProfileGap(config)) continue;
                const float t=std::max(0.f,std::min(1.f,d.dot(turn.point-a)/std::max(1.e-9f,d.squaredNorm())));
                support=std::min(support,(turn.point-a-t*d).norm());
              }
              if (support<=config.profileWallTolerance) featureTurns.push_back(turn);
            }
          }
          std::sort(featureTurns.begin(),featureTurns.end(),
              [](const ProfileTurnSample &a,const ProfileTurnSample &b){return a.point.x()<b.point.x();});
        }
        auto supportedFoot=[&](const ProfileTurnSample &turn,int side) {
          const float h=opening*(baseline.normal.dot(turn.point)-baseline.offset);
          const float before=opening*baseline.normal.dot(turn.before);
          const float after=opening*baseline.normal.dot(turn.after);
          const float outer=side ? after : before;
          const float inner=side ? before : after;
          if (turn.curvature*opening<=0 || turn.angle*opening<minimumTurn) return false;
          if (groove) {
            // Relative to the local plate, the point must be in the lower
            // depression, with a steeper outside and a flatter inside floor.
            return h<=-.55f*amplitude &&
                (side ? outer>std::sin(minimumTurn) : outer<-std::sin(minimumTurn)) &&
                std::fabs(inner)<std::fabs(outer);
          }
          // A low bead has outward floor support and inward rising walls.
          return std::fabs(h)<=.4f*amplitude &&
              (side ? inner<-std::sin(minimumTurn) : inner>std::sin(minimumTurn)) &&
              std::fabs(outer)<std::fabs(inner);
        };
        float bestTurnScore = -std::numeric_limits<float>::infinity();
        for (std::size_t leftIndex = 0; leftIndex < featureTurns.size(); ++leftIndex) {
          const auto &leftTurn = featureTurns[leftIndex];
          if (!supportedFoot(leftTurn,0) ||
              leftTurn.point.x() < leftU - scale ||
              leftTurn.point.x() > rightU + scale)
            continue;
          for (std::size_t rightIndex = leftIndex + 1;
               rightIndex < featureTurns.size(); ++rightIndex) {
            const auto &rightTurn = featureTurns[rightIndex];
            const float turnWidth = rightTurn.point.x() - leftTurn.point.x();
            if (!supportedFoot(rightTurn,1) ||
                turnWidth < std::max(minimumWidth,config.minSeamSeparation) || turnWidth > maximumWidth ||
                rightTurn.point.x() > rightU + scale)
              continue;
            // The turns must surround the measured feature's midpoint and
            // stay on their own side of the independently determined floor.
            if (!profilePairStraddlesBottom(leftTurn.point.x(),rightTurn.point.x(),featureCenter) ||
                !profilePairStraddlesBottom(leftTurn.point.x(),rightTurn.point.x(),bottomMidpoint))
              continue;
            bool connected=true;
            for (std::size_t j=1;j<localizationProfile.size();++j)
              if (localizationProfile[j].u>leftTurn.point.x() &&
                  localizationProfile[j-1].u<rightTurn.point.x() &&
                  localizationProfile[j].u-localizationProfile[j-1].u>maximumProfileGap(config)) {
                connected=false;break;
              }
            if (!connected) {
              if (debug) ++debug->gates["emerging_floor_support_gap"];
              continue;
            }
            const float centerError = std::fabs(
                .5f * (leftTurn.point.x() + rightTurn.point.x()) -
                localizationProfile[peak].u);
            const float score = std::fabs(leftTurn.angle) +
                                std::fabs(rightTurn.angle) -
                                .05f * centerError;
            if (score <= bestTurnScore)
              continue;
            bestTurnScore = score;
            leftPoint = leftTurn.point;
            rightPoint = rightTurn.point;
            selectedCurvature =
                {{leftTurn.curvature, rightTurn.curvature}};
          }
        }

        if (!std::isfinite(bestTurnScore)) {
          if (debug) ++debug->gates["emerging_lower_turns_missing"];
          continue;
        }

        const auto axes = transverseAxes(config.longitudinalAxis);
        auto candidateAt = [&](const Eigen::Vector2f &fine,
                               float curvature) {
          auto nearest = std::lower_bound(
              profile.begin(), profile.end(), fine.x(),
              [](const ProfilePoint &point, float u) { return point.u < u; });
          int index = std::min(pointCount - 1,
                               int(nearest - profile.begin()));
          if (index > 0 &&
              std::fabs(profile[index - 1].u - fine.x()) <
                  std::fabs(profile[index].u - fine.x()))
            --index;
          Candidate candidate;
          candidate.point = profile[index].point;
          candidate.point.data[axes[0]] = fine.x();
          candidate.point.data[axes[1]] = fine.y();
          candidate.fromLocalCorner = true;
          candidate.profileCurvature = curvature;
          candidate.feature = clamp01(
              .12f + .32f * amplitude / std::max(scale, amplitude));
          return candidate;
        };

        SurfaceProposal proposal;
        proposal.pair.left = candidateAt(leftPoint, selectedCurvature[0]);
        proposal.pair.right = candidateAt(rightPoint, selectedCurvature[1]);
        proposal.pair.left.emergingProfile = true;
        proposal.pair.right.emergingProfile = true;
        const Eigen::Vector2f floorDirection = rightPoint - leftPoint;
        if (floorDirection.norm() > 1.e-4f)
          proposal.pair.left.floorDirection =
              proposal.pair.right.floorDirection = floorDirection.normalized();
        proposal.cost = -proposal.pair.left.feature -
                        proposal.pair.right.feature + .25f;
        if (debug) {
          proposal.pair.debug = std::make_shared<DebugPair>();
          proposal.pair.debug->source = "emerging_profile";
          proposal.pair.debug->profileCurvatures = selectedCurvature;
          proposal.pair.debug->profileShape =
              groove
                  ? "shallow_groove"
                  : "central_bead";
          ++debug->gates["emerging_profile_fallback"];
          ++debug->gates["emerging_lower_turns_accepted"];
        }
        if (acceptBottomSides(proposal.pair)) proposals.push_back(std::move(proposal));
      }
    }
    }
  }
  if (config.useProfileCurvature && !directCurveLine && !proposals.empty()) {
    float strongestLeft=0, strongestRight=0, bestPair=0;
    for (const auto &p:proposals) {
      strongestLeft=std::max(strongestLeft,std::fabs(p.pair.left.profileCurvature));
      strongestRight=std::max(strongestRight,std::fabs(p.pair.right.profileCurvature));
    }
    auto peakSupport=[&](const SurfaceProposal &p) {
      const float left=std::fabs(p.pair.left.profileCurvature)/std::max(1.e-6f,strongestLeft);
      const float right=std::fabs(p.pair.right.profileCurvature)/std::max(1.e-6f,strongestRight);
      return multiCandidates ? left*right : std::min(left,right);
    };
    for (const auto &p:proposals) bestPair=std::max(bestPair,peakSupport(p));
    // Boundary hypotheses must compete against the same frame's strongest
    // supported turns. Otherwise each window can label a tiny wall ripple
    // its own "maximum". Keep comparable peaks for cross-frame selection.
    // Joint evidence removes combinations dominated by stronger turns on
    // either side, while retaining similarly supported alternatives for DP.
    const float minimumSupport = .8f;
    proposals.erase(std::remove_if(proposals.begin(),proposals.end(),
        [&](const SurfaceProposal &p){return peakSupport(p)<minimumSupport*bestPair;}),proposals.end());
  }
  std::sort(proposals.begin(), proposals.end(),
            [](const SurfaceProposal &lhs, const SurfaceProposal &rhs) {
              return lhs.cost < rhs.cost;
            });

  std::vector<SurfacePair> result;
  const int maximumPairs = multiCandidates ? config.curvatureSearch.maximumPairs :
      std::max(2, config.maxCandidatesPerSlice / 2);
  const float candidateSpacing = multiCandidates ? config.curvatureSearch.candidateSpacing : config.minCandidateSpacing;
  for (const SurfaceProposal &proposal : proposals) {
    bool duplicate = false;
    for (const SurfacePair &existing : result) {
      if (directCurveLine && proposal.pair.centralBead != existing.centralBead)
        continue;
      if (transverseDistance(proposal.pair.left.point, existing.left.point,
                             config.longitudinalAxis) <
              candidateSpacing &&
          transverseDistance(proposal.pair.right.point, existing.right.point,
                             config.longitudinalAxis) <
              candidateSpacing) {
        duplicate = true;
        break;
      }
    }
    if (!duplicate) {
      result.push_back(proposal.pair);
    }
    if (static_cast<int>(result.size()) >= maximumPairs) {
      break;
    }
  }

  // A failed curve/line fit must not discard a measured short wall. The
  // fallback owns its settings; it never reads the standalone curvature bank.
  if (result.empty() && directCurveLine && config.curveLineFallback.enabled && config.useProfileCurvature) {
    DetectorConfig fallback = config;
    fallback.profileLocalizationMethod = ProfileLocalizationMethod::CurvatureCandidates;
    fallback.profileCurvatureRadius = config.curveLineFallback.splineRadius;
    fallback.profileWallTolerance = config.curveLineFallback.wallTolerance;
    fallback.profileTurnLocation = ProfileTurnLocation::RegionMedian;
    fallback.profileTurnRegionMaxLength = config.curveLineFallback.regionLength;
    fallback.profileTurnRegionFraction = config.curveLineFallback.regionFraction;
    fallback.curvatureSearch = config.curveLineFallback.search;
    wallTrack = previousWallTrack;
    if (debug) {
      ++debug->gates["curve_line_fallback_attempted"];
      debug->fittedProfile.clear(); debug->fittedCurvature.clear(); debug->turnRegions.clear();
    }
    result = findSurfacePairs(profile, localizationProfile, fallback, predictions,
                              wallTrack, hasGeometry, debug, true);
    for (auto &pair : result) {
      pair.left.curvatureFallback = pair.right.curvatureFallback = true;
      if (pair.debug) pair.debug->source = "curve_line_fallback/" + pair.debug->source;
    }
    if (debug && !result.empty()) ++debug->gates["curve_line_fallback_accepted"];
  }
  if (hasGeometry) {
    *hasGeometry = !result.empty();
  }

  // 焊到高层后，三面交角可能已经被焊道完全填圆。此时加入一对历史投影
  // 作为低特征候选，后续全局优化会在“当前三面模型”和“历史趋势”之间选择。
  if (!config.useProfileCurvature && predictions[0].valid && predictions[1].valid) {
    SurfacePair historyPair;
    historyPair.left = nearestHistoryCandidate(profile, predictions[0], config);
    historyPair.right =
        nearestHistoryCandidate(profile, predictions[1], config);
    const int sideAxis = transverseAxes(config.longitudinalAxis)[0];
    if (component(historyPair.right.point, sideAxis) -
            component(historyPair.left.point, sideAxis) >=
        config.minSeamSeparation) {
      bool duplicate = false;
      for (const SurfacePair &existing : result) {
        if (transverseDistance(historyPair.left.point, existing.left.point,
                               config.longitudinalAxis) <
                config.minCandidateSpacing &&
            transverseDistance(historyPair.right.point, existing.right.point,
                               config.longitudinalAxis) <
                config.minCandidateSpacing) {
          duplicate = true;
          break;
        }
      }
      if (!duplicate && acceptBottomSides(historyPair)) {
        if (debug) {
          historyPair.debug = std::make_shared<DebugPair>();
          historyPair.debug->source = "history";
          historyPair.debug->hasBottomMidpoint = true;
          historyPair.debug->bottomMidpointU = bottomMidpoint;
        }
        result.push_back(historyPair);
      }
    }
  }
  return result;
}

bool validateDetectionInput(const CloudT::ConstPtr &input,
                            const DetectorConfig &config, std::string &error) {
  if (!std::isfinite(config.minimumFeatureDepth) || config.minimumFeatureDepth < 0) {
    error = "Minimum feature depth must be finite and nonnegative.";
    return false;
  }
  auto validSearch = [](const CurvatureSearchConfig &s) {
    return std::isfinite(s.windowLength) && s.windowLength > 0 &&
           std::isfinite(s.valleyProminence) && s.valleyProminence > 0 &&
           std::isfinite(s.candidateSpacing) && s.candidateSpacing > 0 &&
           s.maximumPairs >= 2 && s.maximumPairs <= 32 &&
           std::isfinite(s.continuityWeight) && s.continuityWeight >= 0;
  };
  const bool curvature = config.profileLocalizationMethod == ProfileLocalizationMethod::CurvatureCandidates;
  const auto &fallback = config.curveLineFallback;
  if ((curvature && !validSearch(config.curvatureSearch)) ||
      (!curvature && fallback.enabled &&
       (!validSearch(fallback.search) || !std::isfinite(fallback.splineRadius) || fallback.splineRadius <= 0 ||
        !std::isfinite(fallback.wallTolerance) || fallback.wallTolerance <= 0 ||
        !std::isfinite(fallback.regionLength) || fallback.regionLength <= 0 ||
        !std::isfinite(fallback.regionFraction) || fallback.regionFraction <= 0 || fallback.regionFraction > 1))) {
    error = "Invalid active algorithm candidate search or fallback parameters.";
    return false;
  }
  if (config.profileLocalizationMethod !=
          ProfileLocalizationMethod::CurveLineIntersection &&
      config.profileLocalizationMethod !=
          ProfileLocalizationMethod::CurvatureCandidates) {
    error = "Invalid profile localization method.";
    return false;
  }
  if ((config.profileTurnLocation != ProfileTurnLocation::Peak &&
       config.profileTurnLocation != ProfileTurnLocation::RegionMedian &&
       config.profileTurnLocation != ProfileTurnLocation::TurnMedian &&
       config.profileTurnLocation != ProfileTurnLocation::FloorSideBoundary &&
       config.profileTurnLocation != ProfileTurnLocation::WallTurnCentroid) ||
      !std::isfinite(config.profileTurnRegionFraction) || config.profileTurnRegionFraction <= 0 ||
      config.profileTurnRegionFraction > 1 || !std::isfinite(config.profileTurnRegionMaxLength) ||
      config.profileTurnRegionMaxLength <= 0) {
    error = "Invalid profile turn location mode or region limits.";
    return false;
  }
  if (config.profileMergeFrameCount < 1 || config.profileMergeFrameCount > 101 ||
      config.profileMergeFrameCount % 2 == 0) {
    error = "Profile merge frame count must be an odd integer from 1 to 101.";
    return false;
  }
  if (!std::isfinite(config.profileCurvatureRadius) || config.profileCurvatureRadius <= 0) {
    error = "Profile curvature radius must be finite and positive.";
    return false;
  }
  if (!std::isfinite(config.centralBeadMinHeight) ||
      config.centralBeadMinHeight <= 0 ||
      !std::isfinite(config.centralBeadMinWidth) ||
      config.centralBeadMinWidth <= 0) {
    error = "Central bead height and width must be finite and positive.";
    return false;
  }
  if (!std::isfinite(config.profileWallTolerance) || config.profileWallTolerance<=0) {
    error = "Profile wall tolerance must be finite and positive.";
    return false;
  }
  if (!std::isfinite(config.doubleSeamWidth) || config.doubleSeamWidth<=0) {
    error="Double-seam width threshold must be finite and positive.";
    return false;
  }
  if(config.uniformWidthMetric!=UniformWidthMetric::BottomArc && config.uniformWidthMetric!=UniformWidthMetric::Chord) {
    error="Uniform width metric must be bottom arc or chord.";
    return false;
  }
  if(!std::isfinite(config.uniformWidthFilterLength) || config.uniformWidthFilterLength<=0) {
    error="Uniform-width filter length must be finite and positive.";return false;
  }
  if (!std::isfinite(config.endpointTrimWindow) || config.endpointTrimWindow<=0 ||
      !std::isfinite(config.endpointTrimDeviation) || config.endpointTrimDeviation<=0 ||
      !std::isfinite(config.endpointTrimAngleDegrees) || config.endpointTrimAngleDegrees<=0 || config.endpointTrimAngleDegrees>=90) {
    error="Endpoint trim window/deviation must be positive and angle must be between 0 and 90 degrees.";
    return false;
  }
  if (!std::isfinite(config.finalSmoothingLength) || config.finalSmoothingLength <= 0 ||
      !std::isfinite(config.finalSmoothingRobustScale) || config.finalSmoothingRobustScale <= 0) {
    error = "Final smoothing length and robust scale must be finite and positive.";
    return false;
  }
  if (config.interiorOpeningDirection < -1 || config.interiorOpeningDirection > 1 ||
      !std::isfinite(config.maxCandidateCorrection) || config.maxCandidateCorrection < 0) {
    error = "Opening direction or candidate correction limit is invalid.";
    return false;
  }
  // S00 输入与参数检查。失败通过 error 返回，不产生虚构轨迹。
  if (!input || input->size() < 30) {
    error = "Point cloud has fewer than 30 valid points.";
    return false;
  }
  if (config.sliceStep <= 0.0f || config.sliceThickness <= 0.0f ||
      config.profileBinSize <= 0.0f ||
      (stableProfileSampling(config) &&
       (!std::isfinite(config.profileReferenceScale) ||
        config.profileReferenceScale <= 0.0f))) {
    error = "Slice and profile dimensions must be positive.";
    return false;
  }
  if (!std::isfinite(config.floorPatchSize) || config.floorPatchSize <= 0 ||
      !std::isfinite(config.floorPatchGuard) || config.floorPatchGuard < 0 ||
      !std::isfinite(config.floorPatchMaxRms) || config.floorPatchMaxRms <= 0 ||
      !std::isfinite(config.floorPatchMaxCorrection) ||
      config.floorPatchMaxCorrection < 0) {
    error = "Local floor patch parameters are invalid.";
    return false;
  }
  if (config.useRegionGrowingEdges &&
      (config.regionGrowingWindow < 3 ||
       config.regionGrowingSmoothnessDegrees <= 0.0f ||
       config.regionGrowingMaxNormalDriftDegrees <= 0.0f ||
       config.regionGrowingDistance <= 0.0f ||
       config.regionGrowingMaxEdgeGap <= 0.0f)) {
    error = "Region-growing parameters must be positive.";
    return false;
  }
  return true;
}

struct StationRange {
  float minimum;
  float maximum;
};

// A lookup index, not a downsampler: all local predicates still see every
// original point in its original order. Include a floating-point guard band.
class LongitudinalIndex {
public:
  LongitudinalIndex(const CloudT::ConstPtr &cloud, Axis axis) : cloud_(cloud), axis_(axis) {
    order_.resize(cloud_->size());
    std::iota(order_.begin(), order_.end(), 0);
    auto less = [&](std::size_t a, std::size_t b) { return coordinate((*cloud_)[a], axis_) < coordinate((*cloud_)[b], axis_); };
    if (!std::is_sorted(order_.begin(), order_.end(), less))
      std::stable_sort(order_.begin(), order_.end(), less);
  }
  CloudT::ConstPtr around(float station, float radius) const {
    radius += 8 * std::numeric_limits<float>::epsilon() * (std::abs(station) + radius + 1);
    auto first = std::lower_bound(order_.begin(), order_.end(), station - radius,
        [&](std::size_t i, float value) { return coordinate((*cloud_)[i], axis_) < value; });
    auto last = std::upper_bound(first, order_.end(), station + radius,
        [&](float value, std::size_t i) { return value < coordinate((*cloud_)[i], axis_); });
    std::vector<std::size_t> selected(first, last);
    if (!std::is_sorted(selected.begin(), selected.end())) std::sort(selected.begin(), selected.end());
    CloudT::Ptr local(new CloudT);
    local->reserve(selected.size());
    for (auto i : selected) local->push_back((*cloud_)[i]);
    return local;
  }
private:
  CloudT::ConstPtr cloud_;
  Axis axis_;
  std::vector<std::size_t> order_;
};

StationRange findStationRange(const CloudT &cloud, Axis axis) {
  float minimum = std::numeric_limits<float>::max();
  float maximum = std::numeric_limits<float>::lowest();
  for (const PointT &point : cloud) {
    const float station = coordinate(point, axis);
    minimum = std::min(minimum, station);
    maximum = std::max(maximum, station);
  }
  return {minimum, maximum};
}

float findProfileGridOrigin(const CloudT &cloud, Axis longitudinal) {
  const int uAxis = transverseAxes(longitudinal)[0];
  float minimum = std::numeric_limits<float>::max();
  for (const PointT &point : cloud)
    if (pcl::isFinite(point))
      minimum = std::min(minimum, component(point, uAxis));
  return std::isfinite(minimum) ? minimum : 0.0f;
}

struct FrameProjection {
  int frameIndex = -1;
  float station = 0.0f;
  float fallbackLongitudinal = 0.0f;
  Eigen::Vector3f center = Eigen::Vector3f::Zero();
  Eigen::Vector3f normal = Eigen::Vector3f::UnitX();
};

PointT restoreFramePoint(PointT point, const FrameProjection &plane,
                         Axis axis) {
  const int a = static_cast<int>(axis);
  float rhs = plane.normal.dot(plane.center);
  for (int component = 0; component < 3; ++component)
    if (component != a)
      rhs -= plane.normal[component] * point.data[component];
  point.data[a] = std::fabs(plane.normal[a]) > 1.0e-4f
                      ? rhs / plane.normal[a]
                      : plane.fallbackLongitudinal;
  return point;
}

CloudT::Ptr preprocessOneFrame(const CloudT::Ptr &frame,
                               const DetectorConfig &config) {
#ifdef WELD_NO_PCL_FILTERS
  (void)config;
  return frame;
#else
  CloudT::Ptr denoised(new CloudT);
  if (config.outlierMeanK >= 3 &&
      frame->size() > static_cast<std::size_t>(config.outlierMeanK)) {
    pcl::StatisticalOutlierRemoval<PointT> filter;
    filter.setInputCloud(frame);
    filter.setMeanK(config.outlierMeanK);
    filter.setStddevMulThresh(config.outlierStddev);
    filter.filter(*denoised);
  } else {
    *denoised = *frame;
  }
  if (config.voxelSize <= 0.0f)
    return denoised;
  CloudT::Ptr downsampled(new CloudT);
  pcl::VoxelGrid<PointT> voxel;
  voxel.setInputCloud(denoised);
  voxel.setLeafSize(config.voxelSize, config.voxelSize, config.voxelSize);
  voxel.filter(*downsampled);
  return downsampled;
#endif
}

bool prepareFrameSections(const CloudT::ConstPtr &input,
                          const std::vector<ScanFrameRange> &frames,
                          const DetectorConfig &config,
                          CloudT::Ptr &flattened,
                          CloudT::Ptr &measured,
                          std::vector<ScanFrameRange> &processedFrames,
                          std::vector<FrameProjection> &projections,
                          float &frameStep, std::string &error) {
  if (!input || frames.empty()) {
    error = "No original frames were supplied.";
    return false;
  }
  std::size_t expected = 0;
  const Axis axis = config.longitudinalAxis;
  int previousId = -1;
  std::vector<float> rawStations;
  rawStations.reserve(frames.size());
  for (const auto &frame : frames) {
    if (frame.begin != expected || frame.count == 0 ||
        frame.count > input->size() - expected ||
        frame.frameIndex <= previousId) {
      error = "Original frame ranges are incomplete or unordered.";
      return false;
    }
    std::vector<float> values;
    values.reserve(frame.count);
    for (std::size_t i = frame.begin; i < frame.begin + frame.count; ++i) {
      if (!pcl::isFinite((*input)[i])) {
        error = "An original frame contains a non-finite point.";
        return false;
      }
      values.push_back(coordinate((*input)[i], axis));
    }
    rawStations.push_back(median(std::move(values)));
    expected += frame.count;
    previousId = frame.frameIndex;
  }
  if (expected != input->size()) {
    error = "Original frames do not cover the complete point cloud.";
    return false;
  }
  std::vector<float> pitches;
  for (std::size_t i = 1; i < frames.size(); ++i) {
    const int idGap = frames[i].frameIndex - frames[i - 1].frameIndex;
    const float pitch =
        std::fabs(rawStations[i] - rawStations[i - 1]) / float(idGap);
    if (std::isfinite(pitch) && pitch > 1.0e-4f)
      pitches.push_back(pitch);
  }
  frameStep = pitches.empty() ? 1.0f : median(std::move(pitches));
  flattened.reset(new CloudT);
  flattened->reserve(input->size());
  measured.reset(new CloudT);
  measured->reserve(input->size());
  processedFrames.reserve(frames.size());
  projections.reserve(frames.size());
  const int a = static_cast<int>(axis);
  const int firstId = frames.front().frameIndex;
  for (std::size_t f = 0; f < frames.size(); ++f) {
    const auto &range = frames[f];
    FrameProjection projection;
    projection.frameIndex = range.frameIndex;
    projection.station = rawStations.front() +
                         frameStep * float(range.frameIndex - firstId);
    projection.fallbackLongitudinal = rawStations[f];
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    for (std::size_t i = range.begin; i < range.begin + range.count; ++i)
      center += toVector((*input)[i]).cast<double>();
    center /= double(range.count);
    projection.center = center.cast<float>();
    Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
    CloudT::Ptr rawFrame(new CloudT);
    rawFrame->reserve(range.count);
    for (std::size_t i = range.begin; i < range.begin + range.count; ++i) {
      const auto delta =
          (toVector((*input)[i]).cast<double>() - center).eval();
      covariance += delta * delta.transpose();
      rawFrame->push_back((*input)[i]);
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(
        covariance / double(range.count));
    if (solver.info() == Eigen::Success)
      projection.normal = solver.eigenvectors().col(0).cast<float>();
    CloudT::Ptr cleaned = preprocessOneFrame(rawFrame, config);
    if (!cleaned || cleaned->empty()) {
      error = "A frame became empty during preprocessing.";
      return false;
    }
    ScanFrameRange processed;
    processed.begin = flattened->size();
    processed.count = cleaned->size();
    processed.frameIndex = range.frameIndex;
    for (PointT point : *cleaned) {
      measured->push_back(point);
      point.data[a] = projection.station;
      flattened->push_back(point);
    }
    processedFrames.push_back(processed);
    projections.push_back(projection);
  }
  flattened->width = static_cast<std::uint32_t>(flattened->size());
  flattened->height = 1;
  flattened->is_dense = true;
  return true;
}

void restoreOriginalFrameGeometry(SeamResult &result,
                                  const std::vector<int> &frameIds,
                                  const std::vector<FrameProjection> &planes,
                                  Axis axis) {
  result.leftFrameIndices.clear();
  result.rightFrameIndices.clear();
  auto restore = [&](CloudT &cloud, std::vector<int> &ids) {
    const std::size_t count = std::min(cloud.size(), frameIds.size());
    ids.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
      const int frameId = frameIds[i];
      ids.push_back(frameId);
      const auto found = std::find_if(
          planes.begin(), planes.end(),
          [frameId](const FrameProjection &p) { return p.frameIndex == frameId; });
      if (found == planes.end())
        continue;
      cloud[i] = restoreFramePoint(cloud[i], *found, axis);
    }
  };
  restore(*result.leftSeam, result.leftFrameIndices);
  restore(*result.rightSeam, result.rightFrameIndices);
}

void scoreSurfacePairs(std::vector<SurfacePair> &pairs,
                       const std::array<HistoryPrediction, 2> &predictions,
                       const DetectorConfig &config) {
  for (SurfacePair &pair : pairs) {
    Candidate *sideCandidates[2] = {&pair.left, &pair.right};
    for (int side = 0; side < 2; ++side) {
      Candidate &candidate = *sideCandidates[side];
      // 连续性使用真实横截面坐标，不减去每片的中位中心；扫描覆盖
      // 范围变化时，中位中心本身会跳动并制造假的连续轨迹。
      candidate.localOffset = toVector(candidate.point);
      candidate.localOffset[static_cast<int>(config.longitudinalAxis)] = 0.0f;
      const HistoryPrediction &prediction =
          predictions[static_cast<std::size_t>(side)];
      candidate.historyAvailable[static_cast<std::size_t>(side)] =
          prediction.valid;
      if (prediction.valid) {
        candidate.historyDistance[static_cast<std::size_t>(side)] =
            transverseDistance(candidate.point, prediction.point,
                               config.longitudinalAxis);
      }

      // 三面交线是主约束；当交角变弱时，再提高对应一侧历史的权重。
      const float adaptiveHistoryWeight =
          config.historyWeight * prediction.stability *
          (1.0f + config.weakFeatureHistoryBoost * (1.0f - candidate.feature));
      const float historyCost =
          prediction.valid
              ? adaptiveHistoryWeight *
                    std::pow(candidate.historyDistance[static_cast<std::size_t>(
                                 side)] /
                                 std::max(1.0f, config.historySearchRadius),
                             2.0f)
              : 0.0f;
      candidate.nodeCost[static_cast<std::size_t>(side)] =
          -config.featureWeight * candidate.feature + historyCost;
    }
  }
}

std::string describeNoCandidateReason(const DebugSlice &slice,
                                      const DetectorConfig &config) {
  const auto gate = [&](const char *name) {
    const auto found = slice.gates.find(name);
    return found == slice.gates.end() ? std::size_t(0) : found->second;
  };
  const auto number = [](float value) {
    std::ostringstream out;
    out.setf(std::ios::fixed);
    out.precision(2);
    out << value;
    return out.str();
  };

  if (gate("profile_insufficient"))
    return "分箱轮廓点数不足，无法同时拟合左右侧壁和槽底。";

  if (gate("emerging_floor_support_gap"))
    return "浅槽两侧之间存在超限的实测轮廓缺口，未输出跨缺口的左右壁脚点对。";

  if (gate("emerging_lower_turns_missing"))
    return "检测到浅槽或低凸起，但未找到同时满足实测支撑、转折方向和槽底位置的左右壁脚；不使用槽口上沿替代。";

  if (config.useProfileCurvature && gate("curvature_turns_missing")) {
    return "B 样条未形成满足条件的连续高曲率区间：相对曲率阈值 " +
           number(config.profileTurnRegionFraction) + "，区间最大长度 " +
           number(config.profileTurnRegionMaxLength) + "。";
  }

  const bool leftMissing = gate("left_observed_wall_missing") != 0;
  const bool rightMissing = gate("right_observed_wall_missing") != 0;
  if (leftMissing || rightMissing) {
    const float minimumHeight = 2 * profileGeometryScale(config);
    const std::string side = leftMissing && rightMissing
                                 ? "左右主侧壁"
                                 : (leftMissing ? "左主侧壁" : "右主侧壁");
    return "未识别到" + side +
           "；稳定斜壁的切线竖向分量需至少为 0.45，且连续高度当前至少为 " +
           number(minimumHeight) + "。";
  }

  const auto sideReason = [&](const char *prefix,
                              const char *sideName) -> std::string {
    const std::string key(prefix);
    if (!gate((key + "no_supported_wall_turn").c_str()))
      return {};
    std::vector<std::string> details;
    const auto far = gate((key + "turn_far_from_wall_foot").c_str());
    const auto outside = gate((key + "turn_outside_wall_tolerance").c_str());
    const auto angle = gate((key + "turn_angle_too_small").c_str());
    if (far)
      details.push_back(std::to_string(far) + " 个转折距侧壁脚点超过 " +
                        number(std::min(config.profileCurvatureRadius,
                                        5 * profileGeometryScale(config))));
    if (outside)
      details.push_back(std::to_string(outside) + " 个转折到侧壁的垂直距离超过 " +
                        number(config.profileWallTolerance));
    if (angle)
      details.push_back(std::to_string(angle) + " 个转折角小于 " +
                        number(config.minLocalSurfaceAngleDegrees) + "°");
    if (gate((key + "observed_wall_rms_too_large").c_str()))
      details.push_back("侧壁拟合 RMS 超过 " + number(config.profileWallTolerance));
    std::ostringstream out;
    out << sideName << "未找到贴合主侧壁的有效转折";
    if (!details.empty()) {
      out << "（";
      for (std::size_t i = 0; i < details.size(); ++i) {
        if (i) out << "；";
        out << details[i];
      }
      out << "）";
    }
    return out.str();
  };

  const std::string leftReason = sideReason("left_", "左侧");
  const std::string rightReason = sideReason("right_", "右侧");
  if (!leftReason.empty() || !rightReason.empty()) {
    std::string reason = leftReason;
    if (!reason.empty() && !rightReason.empty()) reason += "；";
    reason += rightReason;
    if (leftReason.empty() != rightReason.empty())
      reason += "；另一侧虽有候选，但左右候选必须成对";
    return reason + "。";
  }

  struct GateDescription {
    const char *gate;
    const char *description;
  };
  static const GateDescription descriptions[] = {
      {"candidate_crosses_bottom_midpoint", "左右候选未分别位于底部中点两侧"},
      {"invalid_fit", "侧壁或槽底直线拟合失败"},
      {"slope_or_broad_angle", "槽底斜率或侧壁与槽底夹角不满足阈值"},
      {"invalid_local_fit", "交界附近的局部直线拟合失败"},
      {"local_angle", "局部转折角小于阈值"},
      {"turn_direction", "左右转折方向不能构成槽形"},
      {"wrong_opening_direction", "槽口朝向与设置不一致"},
      {"parallel_lines", "侧壁与槽底接近平行，无法求稳定交点"},
      {"intersection_or_separation", "交点偏离边界或左右间距不足"},
      {"clearance", "侧壁与槽底的有效高度不足"},
      {"curve_line_no_intersection", "斜边曲线与槽底直线没有有效交点"},
      {"curve_line_geometry", "曲线交点的方向或间距不满足约束"},
      {"no_spline_curvature", "交界未匹配到 B 样条曲率候选"},
  };
  std::vector<std::pair<std::size_t, const char *>> rejected;
  for (const auto &item : descriptions) {
    const auto count = gate(item.gate);
    if (count) rejected.emplace_back(count, item.description);
  }
  std::sort(rejected.begin(), rejected.end(),
            [](const auto &a, const auto &b) { return a.first > b.first; });
  if (!rejected.empty()) {
    std::ostringstream out;
    out << "候选组合已生成，但被几何条件过滤：";
    const std::size_t count = std::min<std::size_t>(3, rejected.size());
    for (std::size_t i = 0; i < count; ++i) {
      if (i) out << "；";
      out << rejected[i].second << "（" << rejected[i].first << " 次）";
    }
    out << "。";
    return out.str();
  }
  return "候选未通过当前截面的侧壁贴合、几何或左右配对约束。";
}

// Plate-normal relief is invariant to scan tilt and groove opening. Outer
// support strips establish the plate when available; otherwise use the whole
// section's dominant direction. Median filtering and trimmed extremes keep a
// single spike from making an almost planar section appear deep.
float measuredFeatureDepth(const std::vector<ProfilePoint> &profile,
                           const DetectorConfig &config) {
  if (profile.size() < 9) return -1.f;
  const int count = std::max(3, int(profile.size()) / 5);
  std::vector<ProfilePoint> shoulders(profile.begin(), profile.begin()+count);
  shoulders.insert(shoulders.end(), profile.end()-count, profile.end());
  auto reference = fitOrthogonalLine(shoulders, 0, int(shoulders.size()), true);
  if (!reference.valid || reference.rms > config.profileWallTolerance)
    reference = fitOrthogonalLine(profile, 0, int(profile.size()), true);
  if (!reference.valid) return -1.f;
  std::vector<float> residual, filtered;
  for (const auto &p : profile)
    residual.push_back(reference.normal.dot(Eigen::Vector2f(p.u,p.v))-reference.offset);
  for (int i=0; i<int(residual.size()); ++i) {
    std::vector<float> window(residual.begin()+std::max(0,i-2),
                              residual.begin()+std::min(int(residual.size()),i+3));
    std::sort(window.begin(),window.end());
    filtered.push_back(window[window.size()/2]);
  }
  std::sort(filtered.begin(),filtered.end());
  const std::size_t trim = filtered.size()/50;
  return filtered[filtered.size()-1-trim]-filtered[trim];
}

void appendSlice(std::vector<Slice> &slices, const std::vector<PointT> &points,
                 float station, int frameIndex, const DetectorConfig &config,
                 const SeamHistory *leftHistory,
                 const SeamHistory *rightHistory, DebugTrace *debug,
                 std::array<ObservedWall,2> &wallTrack,
                 const std::vector<PointT> *unmergedPoints = nullptr) {
  const int sideAxis = transverseAxes(config.longitudinalAxis)[0];
  const int requiredProfilePoints =
      2 * (scaledProfileCount(config, config.fitHalfWindow, 5) +
           scaledProfileCount(config, config.fitGuardPoints, 1, 1)) +
      8;
  DebugSlice *debugSlice = nullptr;
  if (debug) {
    debug->slices.emplace_back();
    debugSlice = &debug->slices.back();
    debugSlice->station = station;
    debugSlice->frameIndex = frameIndex;
    debugSlice->raw = points;
  }
  if (points.size() < static_cast<std::size_t>(requiredProfilePoints)) {
    if (debugSlice) {
      debugSlice->noCandidateReason =
          "原始点数不足：当前 " + std::to_string(points.size()) +
          " 点，至少需要 " + std::to_string(requiredProfilePoints) + " 点。";
    }
    return;
  }

  const std::vector<ProfilePoint> profile =
      makeProfile(points, config.longitudinalAxis, config, station,
                  config.longitudinalDetrendMinSpan);
  std::vector<ProfilePoint> localizationProfile;
  if (config.useProfileCurvature) {
    auto fineConfig=config;
    fineConfig.profileBinSize=std::min(config.profileBinSize,.25f*profileGeometryScale(config));
    localizationProfile=makeProfile(points,config.longitudinalAxis,fineConfig,station,
                                     config.longitudinalDetrendMinSpan);
    std::vector<ProfilePoint> cleaned;
    for (std::size_t b=0;b<localizationProfile.size();) {
      std::size_t e=b+1;
      while (e<localizationProfile.size() && localizationProfile[e].u-localizationProfile[e-1].u<=maximumProfileGap(config)) ++e;
      std::vector<Eigen::Vector2f> contour;
      for (auto j=b;j<e;++j) contour.emplace_back(localizationProfile[j].u,localizationProfile[j].v);
      const auto removed=wallProfileSpikeMask(contour,4*profileGeometryScale(config),
                                              2*profileGeometryScale(config),config.profileWallTolerance);
      for (auto j=b;j<e;++j) {
        if (!removed[j-b]) cleaned.push_back(localizationProfile[j]);
        else if (debugSlice) ++debugSlice->gates["profile_spikes_removed"];
      }
      b=e;
    }
    localizationProfile=std::move(cleaned);
    if (debugSlice) for (const auto &p:localizationProfile) debugSlice->localizationProfile.push_back(p.point);
  }
  publishProfileStep(debug, profile, station, config);
  const std::array<HistoryPrediction, 2> predictions{
      {leftHistory ? leftHistory->predict(station, config.longitudinalAxis,
                                          config.maxHistoryExtrapolation,
                                          sideAxis, +1)
                   : HistoryPrediction(),
       rightHistory ? rightHistory->predict(station, config.longitudinalAxis,
                                            config.maxHistoryExtrapolation,
                                            sideAxis, -1)
                    : HistoryPrediction()}};
  if (debugSlice) {
    debugSlice->predictions = predictions;
    for (const auto &p : profile)
      debugSlice->profile.push_back(
          {p.u,
           p.v,
           p.support,
           {{p.fitMean.x(), p.fitMean.y()}},
           {{p.fitScatter(0, 0), p.fitScatter(0, 1), p.fitScatter(1, 1)}}});
    debugSlice->status = "no_candidates";
  }
  bool hasGeometry = false;
  if (config.minimumFeatureDepth > 0) {
    auto depthConfig = config;
    depthConfig.profileBinSize = std::min(config.profileBinSize, .25f*profileGeometryScale(config));
    const auto depthProfile = unmergedPoints
        ? makeProfile(*unmergedPoints, config.longitudinalAxis, depthConfig, station,
                      config.longitudinalDetrendMinSpan)
        : (localizationProfile.empty() ? profile : localizationProfile);
    const float depth = measuredFeatureDepth(depthProfile, config);
    if (debugSlice) debugSlice->featureDepth = depth;
    if (depth >= 0 && depth < config.minimumFeatureDepth) {
      Slice boundary;
      boundary.station = station;
      boundary.frameIndex = frameIndex;
      boundary.featureDepthBoundary = true;
      if (debugSlice) {
        boundary.debugIndex = debug->slices.size()-1;
        debugSlice->status = "feature_depth_cutoff";
        ++debugSlice->gates["feature_depth_cutoff"];
        debugSlice->noCandidateReason = "有效槽深 " + std::to_string(depth) +
            " 小于截断阈值 " + std::to_string(config.minimumFeatureDepth) +
            "：底部接近板面，停止本截面候选计算，禁止轨迹和插补跨越。";
      }
      slices.push_back(std::move(boundary));
      wallTrack = {};
      return;
    }
  }
  std::vector<SurfacePair> pairs =
      findSurfacePairs(profile, localizationProfile, config, predictions, wallTrack, &hasGeometry, debugSlice);
  publishPairStep(debug, pairs, station);
  if (pairs.empty()) {
    if (debugSlice)
      debugSlice->noCandidateReason =
          describeNoCandidateReason(*debugSlice, config);
    return;
  }

  scoreSurfacePairs(pairs, predictions, config);
  Slice slice;
  if (config.uniformWidthEnabled) {
    slice.widthSection.station=station;
    slice.widthSection.bottomCenterU=pairs.front().bottomMidpointU;
    const auto &source=localizationProfile.empty()?profile:localizationProfile;
    // Retain only the connected measured contour containing the bottom;
    // neither missing samples nor a shallow-feature boundary are extrapolated.
    for(std::size_t begin=0;begin<source.size();) {
      auto end=begin+1;
      while(end<source.size() && source[end].u-source[end-1].u<=maximumProfileGap(config))++end;
      if(source[begin].u<slice.widthSection.bottomCenterU && source[end-1].u>slice.widthSection.bottomCenterU) {
        std::vector<Eigen::Vector2f> points;
        for(auto i=begin;i<end;++i)points.emplace_back(source[i].u,source[i].v);
        const auto curve=fitParametricProfileSpline(points,config.profileCurvatureRadius/4,
                                                     std::max(.05f,profileGeometryScale(config)/8));
        const auto axes=transverseAxes(config.longitudinalAxis);
        for(const auto &sample:curve) {
          PointT p=source[begin].point;
          p.data[int(config.longitudinalAxis)]=station;
          p.data[axes[0]]=sample.point.x();p.data[axes[1]]=sample.point.y();
          slice.widthSection.contour.push_back(p);
        }
        break;
      }
      begin=end;
    }
  }
  if (debugSlice) {
    slice.debugIndex = debug->slices.size() - 1;
    debugSlice->hasGeometry = hasGeometry;
    debugSlice->status = "outside_main_segment";
    for (const auto &pair : pairs) {
      DebugPair d = pair.debug ? *pair.debug : DebugPair();
      d.points = {{pair.left.point, pair.right.point}};
      d.feature = {{pair.left.feature, pair.right.feature}};
      d.nodeCost = {{pair.left.nodeCost[0], pair.right.nodeCost[1]}};
      d.historyDistance =
          {{pair.left.historyDistance[0], pair.right.historyDistance[1]}};
      debugSlice->pairs.push_back(std::move(d));
      debugSlice->sourceCandidates[0].push_back(pair.left.point);
      debugSlice->sourceCandidates[1].push_back(pair.right.point);
    }
  }
  slice.station = station;
  slice.frameIndex = frameIndex;
  slice.hasGeometry = hasGeometry;
  slice.pairs = std::move(pairs);
  slices.push_back(std::move(slice));
}

// Pause only after a complete batch has been analysed. This never merges
// frames for fitting and never changes candidate/path computation.
void publishFrameBatch(DebugTrace *debug, std::size_t &begin, bool last,
                       const DetectorConfig &config, std::size_t base,
                       std::size_t total) {
  if (!stepDebugEnabled(debug) || config.debugFrameBatchSize <= 0)
    return;
  const auto end = debug->slices.size();
  if (begin >= end || (!last && end - begin < std::size_t(config.debugFrameBatchSize)))
    return;
  CloudT::Ptr raw(new CloudT), position(new CloudT), profile(new CloudT), left(new CloudT), right(new CloudT), spline(new CloudT);
  CloudT::Ptr overviewLeft(new CloudT), overviewRight(new CloudT);
  const auto axes = transverseAxes(config.longitudinalAxis);
  std::size_t withoutCandidates = 0;
  std::map<std::string, std::size_t> noCandidateReasons;
  for (auto i = begin; i < end; ++i) {
    const auto &slice = debug->slices[i];
    for (const auto &p : slice.raw) raw->push_back(p);
    const auto &source = slice.sourceRaw.empty() ? slice.raw : slice.sourceRaw;
    for (const auto &p : source) position->push_back(p);
    for (const auto &p : slice.fittedProfile) spline->push_back(p);
    for (const auto &p : slice.profile) {
      Eigen::Vector3f xyz = Eigen::Vector3f::Zero();
      xyz[int(config.longitudinalAxis)] = slice.station;
      xyz[axes[0]] = p.u; xyz[axes[1]] = p.v;
      profile->push_back(toPoint(xyz));
    }
    withoutCandidates += slice.pairs.empty();
    if (slice.pairs.empty() && !slice.noCandidateReason.empty())
      ++noCandidateReasons[slice.noCandidateReason];
    for (const auto &pair : slice.pairs) {
      left->push_back(pair.points[0]); right->push_back(pair.points[1]);
    }
  }
  for (std::size_t i = base; i < end; ++i) {
    const auto &slice = debug->slices[i];
    for (const auto &candidate : slice.sourceCandidates[0])
      overviewLeft->push_back(candidate);
    for (const auto &candidate : slice.sourceCandidates[1])
      overviewRight->push_back(candidate);
  }
  const bool framed = debug->slices[begin].frameIndex >= 0;
  const auto first = framed ? debug->slices[begin].frameIndex : int(begin + 1);
  const auto final = framed ? debug->slices[end - 1].frameIndex : int(end);
  const auto firstOrdinal = begin - base + 1;
  const auto lastOrdinal = end - base;
  const double progress = total ? 100.0 * double(lastOrdinal) / double(total) : 100.0;
  const float firstStation = debug->slices[begin].station;
  const float lastStation = debug->slices[end - 1].station;
  std::array<std::size_t, 3> shapeCounts{{0, 0, 0}};
  for (auto i = begin; i < end; ++i)
    for (const auto &pair : debug->slices[i].pairs) {
      if (pair.profileShape == "flat_floor") ++shapeCounts[0];
      else if (pair.profileShape == "wavy_floor") ++shapeCounts[1];
      else if (pair.profileShape == "central_bead") ++shapeCounts[2];
    }
  std::ostringstream title;
  title << "F 按" << (framed ? "帧" : "截面") << "单步 · 第 "
        << firstOrdinal << "–" << lastOrdinal << " / " << total;
  if (framed)
    title << " · 帧号 " << first << "–" << final;
  std::ostringstream message;
  message.setf(std::ios::fixed);
  message.precision(2);
  message << "当前已执行到整个点云的第 " << lastOrdinal << "/" << total
          << (framed ? " 帧" : " 个截面") << "（" << progress
          << "%），纵向位置 " << firstStation;
  if (firstStation != lastStation)
    message << "–" << lastStation;
  message << "。本组 " << (end - begin) << " 个截面，无候选 "
          << withoutCandidates << " 个。候选尚未经过全局路径优化。";
  if (!noCandidateReasons.empty()) {
    if (end - begin == 1) {
      message << " 【无候选原因】" << noCandidateReasons.begin()->first;
    } else {
      std::vector<std::pair<std::string, std::size_t>> reasons(
          noCandidateReasons.begin(), noCandidateReasons.end());
      std::sort(reasons.begin(), reasons.end(),
                [](const auto &a, const auto &b) { return a.second > b.second; });
      message << " 【无候选原因汇总】";
      const std::size_t count = std::min<std::size_t>(3, reasons.size());
      for (std::size_t i = 0; i < count; ++i) {
        if (i) message << "；";
        message << reasons[i].first << "（" << reasons[i].second << " 帧）";
      }
    }
  }
  if (shapeCounts[0] + shapeCounts[1] + shapeCounts[2] > 0)
    message << " 轮廓候选判定：平直/圆滑底 " << shapeCounts[0]
            << "，波动底 " << shapeCounts[1] << "，中央凸起 "
            << shapeCounts[2] << "。";
  if (framed)
    message << " 显示内部展开坐标；最终结果恢复真实扫描平面。";
  CloudT::ConstPtr overview;
  if (!debug->clouds.empty())
    overview = debug->clouds.front().cloud;
  publishDebugStep(debug, title.str(),
      {{"主体截面", raw}, {"分箱轮廓", profile}, {"左候选", left}, {"右候选", right},
       {"B 样条拟合", spline}},
      message.str(),
      {{"first_frame", double(first)}, {"last_frame", double(final)},
       {"first_ordinal", double(firstOrdinal)}, {"last_ordinal", double(lastOrdinal)},
       {"total_slices", double(total)}, {"progress_percent", progress},
       {"frame_based", framed ? 1.0 : 0.0},
       {"first_station", double(firstStation)}, {"last_station", double(lastStation)},
       {"batch_count", double(end - begin)}, {"without_candidates", double(withoutCandidates)}},
      true, overview, position, {{overviewLeft, overviewRight}});
  begin = end;
}

std::vector<Slice> collectSlices(const CloudT &cloud, const StationRange &range,
                                 const DetectorConfig &config,
                                 const SeamHistory *leftHistory,
                                 const SeamHistory *rightHistory,
                                 DebugTrace *debug) {
  // S04 从 min(A)+厚度/2 开始逐站位收集 abs(A-station)<=厚度/2 的点。
  // 这是一组几何截面，不是传感器原始帧；步长小于厚度时截面允许重叠。
  std::vector<Slice> slices;
  const float minimum = range.minimum, maximum = range.maximum;
  std::array<ObservedWall,2> wallTrack;
  const float halfThickness = 0.5f * config.sliceThickness;
  std::size_t total = 0;
  for (float station = minimum + halfThickness;
       station <= maximum - halfThickness; station += config.sliceStep)
    ++total;
  const std::size_t debugBase = debug ? debug->slices.size() : 0;
  const std::size_t startOffset = total
      ? std::min<std::size_t>(std::size_t(std::max(1, config.debugStartSlice) - 1), total - 1)
      : 0;
  std::size_t batchBegin = debugBase + startOffset;
  for (float station = minimum + halfThickness;
       station <= maximum - halfThickness; station += config.sliceStep) {
    std::vector<PointT> points;
    for (const PointT &point : cloud) {
      if (std::fabs(coordinate(point, config.longitudinalAxis) - station) <=
          halfThickness) {
        points.push_back(point);
      }
    }
    if (stepDebugEnabled(debug))
      publishDebugStep(debug, "S04 纵向切片 · " + std::to_string(station),
                       {{"截面原始点", pointsCloud(points)}});
    appendSlice(slices, points, station, -1, config, leftHistory,
                rightHistory, debug, wallTrack);
    publishFrameBatch(debug, batchBegin, false, config, debugBase, total);
  }
  publishFrameBatch(debug, batchBegin, true, config, debugBase, total);
  return slices;
}

std::vector<Slice>
collectFrameSlices(const CloudT &cloud,
                   const std::vector<ScanFrameRange> &frames,
                   const CloudT &measured,
                   const std::vector<FrameProjection> &planes,
                   const DetectorConfig &config,
                   const SeamHistory *leftHistory,
                   const SeamHistory *rightHistory, DebugTrace *debug) {
  std::vector<Slice> slices;
  std::array<ObservedWall,2> wallTrack;
  const std::size_t total = frames.size();
  const std::size_t debugBase = debug ? debug->slices.size() : 0;
  std::size_t startOffset = total
      ? std::min<std::size_t>(std::size_t(std::max(1, config.debugStartSlice) - 1), total - 1)
      : 0;
  if (total && config.debugSelectedFrameIndex >= 0) {
    const auto selected = std::lower_bound(
        frames.begin(), frames.end(), config.debugSelectedFrameIndex,
        [](const ScanFrameRange &frame, int id) { return frame.frameIndex < id; });
    startOffset = selected == frames.end() ? total - 1
                                           : std::size_t(selected - frames.begin());
  }
  std::size_t batchBegin = debugBase + startOffset;
  const int halfWindow = config.mergeProfileFrames ? config.profileMergeFrameCount / 2 : 0;
  const int axis = static_cast<int>(config.longitudinalAxis);
  for (std::size_t f = 0; f < frames.size(); ++f) {
    const auto &frame = frames[f];
    const auto &plane = planes[f];
    const float station = plane.station;
    std::size_t begin = f, end = f;
    // Do not borrow observations across missing/cropped frames.
    while (begin > 0 && f - begin < std::size_t(halfWindow) &&
           frames[begin].frameIndex - frames[begin - 1].frameIndex == 1)
      --begin;
    while (end + 1 < frames.size() && end - f < std::size_t(halfWindow) &&
           frames[end + 1].frameIndex - frames[end].frameIndex == 1)
      ++end;
    std::vector<PointT> points;
    if (halfWindow == 0) {
      // Preserve the single-frame path exactly when merging is disabled.
      points.assign(cloud.begin() + frame.begin, cloud.begin() + frame.begin + frame.count);
    } else {
      points.reserve(frames[end].begin + frames[end].count - frames[begin].begin);
      for (std::size_t j = frames[begin].begin; j < frames[end].begin + frames[end].count; ++j) {
        const Eigen::Vector3f p = toVector(measured[j]);
        const Eigen::Vector3f projected = p - plane.normal * plane.normal.dot(p - plane.center);
        PointT point(projected.x(), projected.y(), projected.z());
        // Store U/V of the real orthogonal projection in the existing
        // section coordinates; restore the target scan plane at output.
        point.data[axis] = station;
        points.push_back(point);
      }
    }
    if (stepDebugEnabled(debug))
      publishDebugStep(debug,
                       "S04 原始帧截面 · " + std::to_string(frame.frameIndex) +
                           " · 投影合并 " + std::to_string(end - begin + 1) + " 帧",
                       {{"截面拟合输入", pointsCloud(points)}});
    std::vector<PointT> unmerged;
    if (halfWindow > 0 && config.minimumFeatureDepth > 0)
      unmerged.assign(cloud.begin()+frame.begin, cloud.begin()+frame.begin+frame.count);
    appendSlice(slices, points, station, frame.frameIndex, config, leftHistory,
                rightHistory, debug, wallTrack, unmerged.empty() ? nullptr : &unmerged);
    if (debug) {
      auto &debugSlice = debug->slices.back();
      debugSlice.sourceRaw.assign(measured.begin() + frame.begin,
                                  measured.begin() + frame.begin + frame.count);
      for (auto &side : debugSlice.sourceCandidates)
        for (auto &candidate : side)
          candidate = restoreFramePoint(candidate, plane,
                                        config.longitudinalAxis);
      auto &gates = debugSlice.gates;
      gates["profile_merged_frames"] = int(end - begin + 1);
      gates["profile_first_frame"] = frames[begin].frameIndex;
      gates["profile_last_frame"] = frames[end].frameIndex;
    }
    publishFrameBatch(debug, batchBegin, false, config, debugBase, total);
  }
  publishFrameBatch(debug, batchBegin, true, config, debugBase, total);
  return slices;
}

// Dense scans can miss more than three frames over a very short physical
// interval. Reconnect only independently observed, compatible wall pairs.
bool supportedPairGap(const Candidate &al,const Candidate &ar,
                      const Candidate &bl,const Candidate &br,float gap,
                      const DetectorConfig &config) {
  const float scale=profileGeometryScale(config);
  if (!config.useProfileCurvature || gap<=0 || gap>2*scale) return false;
  const int u=transverseAxes(config.longitudinalAxis)[0];
  if (std::fabs((br.point.data[u]-bl.point.data[u])-(ar.point.data[u]-al.point.data[u]))>2*scale) return false;
  for (const auto &pair:{std::make_pair(&al,&bl),std::make_pair(&ar,&br)}) {
    if (!pair.first->trackedWall || !pair.second->trackedWall ||
        // Curve tangents at a wavy floor junction can rotate noticeably over
        // a handful of missing dense-scan frames even when the two measured
        // junctions and the groove width remain continuous. Keep a direction
        // check, but let the position/width gates carry the short-gap test.
        std::fabs(pair.first->trackedNormal.dot(pair.second->trackedNormal))<.75f ||
        transverseDistance(pair.first->point,pair.second->point,config.longitudinalAxis)>2*scale+2*gap) return false;
  }
  return true;
}

std::size_t retainMainSegment(std::vector<Slice> &slices,
                              const DetectorConfig &config,
                              bool hasBothHistories, DebugTrace *debug) {
  // 三面模型偶尔会在光滑上表面产生一两个孤立候选。焊缝应沿纵向形成
  // 连续主体，因此按“当前几何锚点”分段，只保留锚点最多的连续段。历史可
  // 在末端提供少量只参与优化的前视点，但最终输出仍严格截在几何锚点内。
  // S08 根据有效几何锚点选取主体连续段，避免把首尾圆滑上表面误当成焊缝。
  // 输出范围由可靠结构决定；本独立版无历史输入，不使用历史前视。
  std::size_t outputSliceCount = 0;
  if (!slices.empty()) {
    const float maximumSegmentGap =
        static_cast<float>(std::max(1, config.maxInterpolatedGapStations + 1)) *
        config.sliceStep;
    std::vector<std::size_t> anchors;
    std::vector<std::size_t> depthBoundaries(slices.size()+1, 0);
    for (std::size_t i = 0; i < slices.size(); ++i) {
      depthBoundaries[i+1] = depthBoundaries[i] + (slices[i].featureDepthBoundary ? 1 : 0);
      if (slices[i].hasGeometry) {
        anchors.push_back(i);
      }
    }

    std::size_t bestBegin = 0;
    std::size_t bestEnd = slices.size();
    if (!anchors.empty()) {
      std::size_t bestAnchorBegin = 0;
      std::size_t bestAnchorEnd = 1;
      std::size_t anchorBegin = 0;
      for (std::size_t i = 1; i <= anchors.size(); ++i) {
        const bool crossesDepthBoundary = i < anchors.size() &&
            depthBoundaries[anchors[i]] > depthBoundaries[anchors[i-1]];
        bool segmentEnds =
            i == anchors.size() || crossesDepthBoundary ||
            slices[anchors[i]].station - slices[anchors[i - 1]].station >
                maximumSegmentGap + 0.25f * config.sliceStep;
        if (segmentEnds && !crossesDepthBoundary && i<anchors.size()) {
          const auto &a=slices[anchors[i-1]],&b=slices[anchors[i]];
          for (const auto &pa:a.pairs) for (const auto &pb:b.pairs)
            if (supportedPairGap(pa.left,pa.right,pb.left,pb.right,b.station-a.station,config)) segmentEnds=false;
        }
        if (!segmentEnds) {
          continue;
        }
        if (i - anchorBegin > bestAnchorEnd - bestAnchorBegin) {
          bestAnchorBegin = anchorBegin;
          bestAnchorEnd = i;
        }
        anchorBegin = i;
      }
      bestBegin = anchors[bestAnchorBegin];
      bestEnd = anchors[bestAnchorEnd - 1] + 1;
    }
    std::size_t retainedEnd = bestEnd;
    if (!anchors.empty() && hasBothHistories) {
      retainedEnd =
          std::min(slices.size(),
                   bestEnd + static_cast<std::size_t>(
                                 std::max(0, config.historyLookaheadStations)));
      for (std::size_t i=bestEnd; i<retainedEnd; ++i)
        if (slices[i].featureDepthBoundary) { retainedEnd=i; break; }
      outputSliceCount = bestEnd - bestBegin;
    }
    if (bestBegin != 0 || retainedEnd != slices.size()) {
      std::vector<Slice> continuous;
      continuous.reserve(retainedEnd - bestBegin);
      for (std::size_t i = bestBegin; i < retainedEnd; ++i) {
        continuous.push_back(std::move(slices[i]));
      }
      slices = std::move(continuous);
    }
    // With no geometric anchors an entirely planar input only contains
    // boundary markers; they must never reach path optimization as nodes.
    if (anchors.empty() && depthBoundaries.back() > 0) slices.clear();
  }

  if (debug) {
    for (std::size_t i = 0; i < slices.size(); ++i) {
      auto &d = debug->slices[slices[i].debugIndex];
      d.retained = true;
      d.lookahead = outputSliceCount > 0 && i >= outputSliceCount;
      d.status = d.lookahead ? "lookahead_only" : "retained";
    }
  }
  publishRetainedSlices(debug, slices);
  return outputSliceCount;
}

float nodeCost(const SurfacePair &pair) {
  return pair.left.nodeCost[0] + pair.right.nodeCost[1];
}

bool isEmergingProfilePair(const SurfacePair &pair) {
  return pair.left.emergingProfile && pair.right.emergingProfile;
}

float pairAllowedJump(const SurfacePair &previousPair,
                      const SurfacePair &currentPair, float stationGap,
                      const DetectorConfig &config) {
  float allowed = config.maxTransitionDistance *
                  std::max(1.0f, stationGap / config.sliceStep);
  const bool previousEmerging = isEmergingProfilePair(previousPair);
  const bool currentEmerging = isEmergingProfilePair(currentPair);
  if (previousEmerging || currentEmerging) {
    const float scale = profileGeometryScale(config);
    const float emergingLimit =
        (previousEmerging && currentEmerging)
            ? std::max(2.f * scale, 5.f * stationGap)
            : std::max(3.f * scale, 6.f * stationGap);
    allowed = std::min(allowed, emergingLimit);
  }
  return allowed;
}

float pairTransitionCost(const SurfacePair &previousPair,
                         const SurfacePair &currentPair, float stationGap,
                         const DetectorConfig &config) {
  const int sideAxis = transverseAxes(config.longitudinalAxis)[0];
  const float allowedJump = pairAllowedJump(
      previousPair, currentPair, stationGap, config);
  const float leftDistance = transverseOffsetDistance(
      currentPair.left, previousPair.left, config.longitudinalAxis);
  const float rightDistance = transverseOffsetDistance(
      currentPair.right, previousPair.right, config.longitudinalAxis);
  const float leftNormalized = leftDistance / std::max(1.0f, allowedJump);
  const float rightNormalized = rightDistance / std::max(1.0f, allowedJump);
  float transition =
      config.smoothnessWeight *
      (leftNormalized * leftNormalized + rightNormalized * rightNormalized);
  if (leftDistance > allowedJump) {
    transition += 20.0f * (leftDistance - allowedJump) / allowedJump;
  }
  if (rightDistance > allowedJump) {
    transition += 20.0f * (rightDistance - allowedJump) / allowedJump;
  }
  const float currentWidth = component(currentPair.right.point, sideAxis) -
                             component(currentPair.left.point, sideAxis);
  const float previousWidth = component(previousPair.right.point, sideAxis) -
                              component(previousPair.left.point, sideAxis);
  const float widthChange =
      (currentWidth - previousWidth) / std::max(1.0f, allowedJump);
  return transition +
         0.75f * config.smoothnessWeight * widthChange * widthChange;
}

float tangentPredictionCost(const SurfacePair &beforePrevious,
                            const SurfacePair &previous,
                            const SurfacePair &current, float previousGap,
                            float currentGap, const DetectorConfig &config) {
  const float ratio = currentGap / std::max(1.0e-4f, previousGap);
  auto predictionError = [&](const Candidate &first, const Candidate &second,
                             const Candidate &third) {
    Eigen::Vector3f predicted =
        second.localOffset + (second.localOffset - first.localOffset) * ratio;
    Eigen::Vector3f difference = third.localOffset - predicted;
    difference[static_cast<int>(config.longitudinalAxis)] = 0.0f;
    return difference.norm();
  };
  const float leftError =
      predictionError(beforePrevious.left, previous.left, current.left);
  const float rightError =
      predictionError(beforePrevious.right, previous.right, current.right);
  const float allowedJump = std::min(
      pairAllowedJump(previous, current, currentGap, config),
      pairAllowedJump(beforePrevious, previous, previousGap, config));
  // 预测误差尺度小于单步最大位移：斜率可以很大，但斜率不应突然变化。
  const float curvatureScale = std::max(1.0f, 0.35f * allowedJump);
  const float leftNormalized = leftError / curvatureScale;
  const float rightNormalized = rightError / curvatureScale;
  float cost = config.curvatureWeight * (leftNormalized * leftNormalized +
                                         rightNormalized * rightNormalized);
  const bool curvature = config.profileLocalizationMethod == ProfileLocalizationMethod::CurvatureCandidates;
  const auto &search = curvature ? config.curvatureSearch : config.curveLineFallback.search;
  if (config.useProfileCurvature && search.multiCandidates &&
      (curvature || current.left.curvatureFallback || previous.left.curvatureFallback ||
       beforePrevious.left.curvatureFallback)) {
    // Compare geometrically supported alternatives at the localization scale,
    // rather than the much larger maximum permitted inter-frame jump.
    const float scale = std::max(.25f, curvature ? config.profileWallTolerance :
                                               config.curveLineFallback.wallTolerance);
    auto huber = [](float x) { return x <= 1.f ? x*x : 2.f*x-1.f; };
    cost += search.continuityWeight * (huber(leftError/scale) + huber(rightError/scale));
  }
  if (leftError > allowedJump) {
    cost += 30.0f * (leftError - allowedJump) / std::max(1.0f, allowedJump);
  }
  if (rightError > allowedJump) {
    cost += 30.0f * (rightError - allowedJump) / std::max(1.0f, allowedJump);
  }
  return cost;
}

bool optimizePairedPaths(const std::vector<Slice> &slices,
                         const DetectorConfig &config, DebugTrace *debug,
                         std::vector<Candidate> &leftPath,
                         std::vector<Candidate> &rightPath,
                         std::string &error) {
  // 状态保存连续两个截面的候选编号。这样不仅能限制相邻距离，还能用
  // 前两点的切线预测当前点；连续陡坡可以通过，突然换到另一块面会产生
  // 很大的预测残差。左右交点仍必须来自同一套截面结构候选。
  // S09 左右成对的二阶动态规划。
  // 状态=(前一截面候选编号, 当前截面候选编号)，从更早一片转移。
  // 总代价=特征节点代价+左右位移/宽度变化代价+前两点切线预测误差。
  // stateParents 记录最优前驱；回溯 selected 后得到 leftPath/rightPath。
  // 成对优化防止左右焊缝来自互不一致的截面结构。
  std::vector<std::vector<int>> stateParents(slices.size());
  const std::size_t firstCount = slices[0].pairs.size();
  const std::size_t secondCount = slices[1].pairs.size();
  std::vector<float> stateCosts(firstCount * secondCount,
                                std::numeric_limits<float>::infinity());
  stateParents[1].assign(firstCount * secondCount, -1);
  const float firstGap = slices[1].station - slices[0].station;
  for (std::size_t first = 0; first < firstCount; ++first) {
    for (std::size_t second = 0; second < secondCount; ++second) {
      stateCosts[first * secondCount + second] =
          nodeCost(slices[0].pairs[first]) +
          pairTransitionCost(slices[0].pairs[first], slices[1].pairs[second],
                             firstGap, config) +
          nodeCost(slices[1].pairs[second]);
    }
  }

  for (std::size_t i = 2; i < slices.size(); ++i) {
    const std::size_t beforeCount = slices[i - 2].pairs.size();
    const std::size_t previousCount = slices[i - 1].pairs.size();
    const std::size_t currentCount = slices[i].pairs.size();
    std::vector<float> nextCosts(previousCount * currentCount,
                                 std::numeric_limits<float>::infinity());
    stateParents[i].assign(previousCount * currentCount, -1);
    const float previousGap = slices[i - 1].station - slices[i - 2].station;
    const float currentGap = slices[i].station - slices[i - 1].station;
    for (std::size_t previous = 0; previous < previousCount; ++previous) {
      for (std::size_t current = 0; current < currentCount; ++current) {
        const float base =
            pairTransitionCost(slices[i - 1].pairs[previous],
                               slices[i].pairs[current], currentGap, config) +
            nodeCost(slices[i].pairs[current]);
        const std::size_t state = previous * currentCount + current;
        for (std::size_t before = 0; before < beforeCount; ++before) {
          const float preceding = stateCosts[before * previousCount + previous];
          if (!std::isfinite(preceding)) {
            continue;
          }
          const float total =
              preceding + base +
              tangentPredictionCost(
                  slices[i - 2].pairs[before], slices[i - 1].pairs[previous],
                  slices[i].pairs[current], previousGap, currentGap, config);
          if (total < nextCosts[state]) {
            nextCosts[state] = total;
            stateParents[i][state] = static_cast<int>(before);
          }
        }
      }
    }
    stateCosts = std::move(nextCosts);
  }

  const auto last = std::min_element(stateCosts.begin(), stateCosts.end());
  if (last == stateCosts.end() || !std::isfinite(*last)) {
    error = "No continuous left/right seam pair could be formed.";
    return false;
  }

  leftPath.resize(slices.size());
  if (debug)
    debug->optimalCost = *last;
  rightPath.resize(slices.size());
  std::vector<int> selected(slices.size(), -1);
  const std::size_t lastCurrentCount = slices.back().pairs.size();
  const std::size_t lastState =
      static_cast<std::size_t>(std::distance(stateCosts.begin(), last));
  selected[slices.size() - 2] = static_cast<int>(lastState / lastCurrentCount);
  selected[slices.size() - 1] = static_cast<int>(lastState % lastCurrentCount);
  for (std::size_t i = slices.size() - 1; i >= 2; --i) {
    const std::size_t currentCount = slices[i].pairs.size();
    const std::size_t state =
        static_cast<std::size_t>(selected[i - 1]) * currentCount +
        static_cast<std::size_t>(selected[i]);
    selected[i - 2] = stateParents[i][state];
  }
  for (std::size_t i = 0; i < slices.size(); ++i) {
    if (selected[i] < 0 ||
        selected[i] >= static_cast<int>(slices[i].pairs.size())) {
      error = "Internal paired-path backtracking failed.";
      return false;
    }
    const SurfacePair &pair =
        slices[i].pairs[static_cast<std::size_t>(selected[i])];
    if (debug) {
      auto &d = debug->slices[slices[i].debugIndex];
      d.selected = selected[i];
      d.selectedNodeCost = nodeCost(pair);
      if (i > 0)
        d.transitionCost = pairTransitionCost(
            slices[i - 1].pairs[selected[i - 1]], pair,
            slices[i].station - slices[i - 1].station, config);
      if (i > 1)
        d.tangentCost = tangentPredictionCost(
            slices[i - 2].pairs[selected[i - 2]],
            slices[i - 1].pairs[selected[i - 1]], pair,
            slices[i - 1].station - slices[i - 2].station,
            slices[i].station - slices[i - 1].station, config);
    }
    leftPath[i] = pair.left;
    rightPath[i] = pair.right;
    leftPath[i].point.data[static_cast<int>(config.longitudinalAxis)] =
        slices[i].station;
    rightPath[i].point.data[static_cast<int>(config.longitudinalAxis)] =
        slices[i].station;
  }
  return true;
}

float predictTransverseCoordinate(const std::vector<Candidate> &path,
                                  float station, const DetectorConfig &config) {
  const int sideAxis = transverseAxes(config.longitudinalAxis)[0];
  const int axis = static_cast<int>(config.longitudinalAxis);
  const auto upper =
      std::lower_bound(path.begin(), path.end(), station,
                       [&](const Candidate &candidate, float value) {
                         return candidate.point.data[axis] < value;
                       });
  if (upper == path.begin()) {
    return component(path.front().point, sideAxis);
  }
  if (upper == path.end()) {
    return component(path.back().point, sideAxis);
  }
  const Candidate &after = *upper;
  const Candidate &before = *(upper - 1);
  const float beforeStation = before.point.data[axis];
  const float afterStation = after.point.data[axis];
  const float ratio = (station - beforeStation) /
                      std::max(1.0e-4f, afterStation - beforeStation);
  return component(before.point, sideAxis) +
         ratio * (component(after.point, sideAxis) -
                  component(before.point, sideAxis));
}

Candidate continuityCandidate(const std::vector<ProfilePoint> &profile,
                              float station, float predictedTransverse,
                              const DetectorConfig &config) {
  const int sideAxis = transverseAxes(config.longitudinalAxis)[0];
  Candidate candidate;
  const auto nearest = std::min_element(
      profile.begin(), profile.end(),
      [predictedTransverse](const ProfilePoint &lhs, const ProfilePoint &rhs) {
        return std::fabs(lhs.u - predictedTransverse) <
               std::fabs(rhs.u - predictedTransverse);
      });
  candidate.point = nearest->point;
  candidate.point.data[static_cast<int>(config.longitudinalAxis)] = station;
  // 横向保持连续模型的位置，高度取本次点云中该横向附近的实测值。
  candidate.point.data[sideAxis] = predictedTransverse;
  candidate.feature = 0.08f;
  return candidate;
}

void interpolateInternalGaps(const CloudT &cloud, const StationRange &range,
                             const std::vector<Slice> &slices,
                             const DetectorConfig &config,
                             std::vector<Candidate> &leftPath,
                             std::vector<Candidate> &rightPath,
                             DebugTrace *debug,
                             const std::vector<ScanFrameRange> *frames = nullptr,
                             std::vector<int> *outputFrameIds = nullptr) {
  // S10 只在前后可靠锚点之间填补少量内部缺口。
  // U 由前后路径插值，V 取本次截面附近实测高度；首尾和长缺口不补。
  // 默认最多补 3 个站位，并继续保持左右最小间距。
  const float minimum = range.minimum, maximum = range.maximum;
  const float halfThickness = .5f * config.sliceThickness;
  auto gapSupported=[&](std::size_t anchor) {
    if (anchor==0 || anchor>=slices.size()) return false;
    const float gap=slices[anchor].station-slices[anchor-1].station;
    return gap<=float(config.maxInterpolatedGapStations+1)*config.sliceStep ||
        supportedPairGap(leftPath[anchor-1],rightPath[anchor-1],leftPath[anchor],rightPath[anchor],gap,config);
  };
  auto supportedMeasurements=[&](Candidate &left,Candidate &right,float station,std::size_t anchor) {
    const float first=slices[anchor-1].station,last=slices[anchor].station;
    if (last-first<=float(config.maxInterpolatedGapStations+1)*config.sliceStep) return true;
    const float t=(station-first)/(last-first);
    const auto axes=transverseAxes(config.longitudinalAxis);
    for (int side=0;side<2;++side) {
      auto &c=side ? right : left;
      const auto &path=side ? rightPath : leftPath;
      const auto &a=path[anchor-1],&b=path[anchor];
      const float sign=a.trackedNormal.dot(b.trackedNormal)<0 ? -1.f : 1.f;
      const Eigen::Vector2f normal=(1-t)*a.trackedNormal+t*sign*b.trackedNormal;
      if (normal.norm()<.9f) return false;
      const float offset=((1-t)*a.trackedOffset+t*sign*b.trackedOffset)/normal.norm();
      const Eigen::Vector2f n=normal.normalized(),p(c.point.data[axes[0]],c.point.data[axes[1]]);
      // The new point's height comes from this frame, not an extrapolated tail.
      if (std::fabs(n.dot(p)-offset)>config.profileWallTolerance) return false;
      c.trackedWall=true;c.trackedNormal=n;c.trackedOffset=offset;
    }
    return true;
  };
  if (frames) {
    std::vector<Candidate> completeLeft, completeRight;
    std::size_t anchor = 0;
    if (outputFrameIds)
      outputFrameIds->clear();
    for (const auto &frame : *frames) {
      const float station =
          coordinate(cloud[frame.begin], config.longitudinalAxis);
      if (anchor < slices.size() &&
          std::fabs(slices[anchor].station - station) <=
              .25f * config.sliceStep) {
        completeLeft.push_back(leftPath[anchor]);
        completeRight.push_back(rightPath[anchor]);
        if (outputFrameIds)
          outputFrameIds->push_back(frame.frameIndex);
        ++anchor;
        continue;
      }
      if (!gapSupported(anchor))
        continue;
      std::vector<PointT> points(cloud.begin() + frame.begin,
                                 cloud.begin() + frame.begin + frame.count);
      const auto profile =
          makeProfile(points, config.longitudinalAxis, config, station,
                      config.longitudinalDetrendMinSpan);
      if (profile.empty())
        continue;
      float predictedLeft =
          predictTransverseCoordinate(leftPath, station, config);
      float predictedRight =
          predictTransverseCoordinate(rightPath, station, config);
      if (predictedRight - predictedLeft < config.minSeamSeparation) {
        const float center = .5f * (predictedLeft + predictedRight);
        predictedLeft = center - .5f * config.minSeamSeparation;
        predictedRight = center + .5f * config.minSeamSeparation;
      }
      auto left=continuityCandidate(profile,station,predictedLeft,config);
      auto right=continuityCandidate(profile,station,predictedRight,config);
      if (!supportedMeasurements(left,right,station,anchor)) continue;
      completeLeft.push_back(left);completeRight.push_back(right);
      if (outputFrameIds)
        outputFrameIds->push_back(frame.frameIndex);
    }
    leftPath = std::move(completeLeft);
    rightPath = std::move(completeRight);
    if (debug) {
      for (const auto &c : leftPath)
        debug->interpolatedPath[0].push_back(c.point);
      for (const auto &c : rightPath)
        debug->interpolatedPath[1].push_back(c.point);
      debug->pathComplete = true;
      publishPathStep(debug, "S10 内部短缺口插补", leftPath, rightPath);
    }
    return;
  }
  // 只填补两个可靠锚点之间的少量漏检站位。首尾贴着弧形上表面的区段
  // 保持为空，不再为了视觉连续而生成实际上并不存在的焊缝。
  if (!leftPath.empty() && !rightPath.empty()) {
    std::vector<Candidate> completeLeft;
    std::vector<Candidate> completeRight;
    std::size_t anchor = 0;
    for (float station = minimum + halfThickness;
         station <= maximum - halfThickness + 0.25f * config.sliceStep;
         station += config.sliceStep) {
      if (anchor < slices.size() &&
          std::fabs(slices[anchor].station - station) <=
              0.25f * config.sliceStep) {
        completeLeft.push_back(leftPath[anchor]);
        completeRight.push_back(rightPath[anchor]);
        ++anchor;
        continue;
      }

      // 没有前后两个三面交线锚点时属于首尾弱特征区，直接不输出。内部
      // 缺口过长也不跨越，防止轨迹再次沿着弧形上表面延伸。
      if (!gapSupported(anchor)) {
        continue;
      }

      std::vector<PointT> stationPoints;
      for (const PointT &point : cloud) {
        if (std::fabs(coordinate(point, config.longitudinalAxis) - station) <=
            halfThickness) {
          stationPoints.push_back(point);
        }
      }
      const std::vector<ProfilePoint> profile = makeProfile(
          stationPoints, config.longitudinalAxis, config, station,
          config.longitudinalDetrendMinSpan);
      if (profile.empty()) {
        continue;
      }
      float predictedLeft =
          predictTransverseCoordinate(leftPath, station, config);
      float predictedRight =
          predictTransverseCoordinate(rightPath, station, config);
      if (predictedRight - predictedLeft < config.minSeamSeparation) {
        // 短缺口插值仍保持左右最小间距，避免两条焊缝在缺口中合并。
        const float center = 0.5f * (predictedLeft + predictedRight);
        predictedLeft = center - 0.5f * config.minSeamSeparation;
        predictedRight = center + 0.5f * config.minSeamSeparation;
      }
      auto left=continuityCandidate(profile,station,predictedLeft,config);
      auto right=continuityCandidate(profile,station,predictedRight,config);
      if (!supportedMeasurements(left,right,station,anchor)) continue;
      completeLeft.push_back(left);completeRight.push_back(right);
    }
    leftPath = std::move(completeLeft);
    rightPath = std::move(completeRight);
  }

  if (debug) {
    for (const auto &c : leftPath)
      debug->interpolatedPath[0].push_back(c.point);
    for (const auto &c : rightPath)
      debug->interpolatedPath[1].push_back(c.point);
    debug->pathComplete = true;
    publishPathStep(debug, "S10 内部短缺口插补", leftPath, rightPath);
  }
}

std::vector<Candidate> regularizePath(const std::vector<Candidate> &path,
                                      const DetectorConfig &config) {
  // S11 平滑定位：分别对 U/V 解加权最小二乘，A 保持不变。
  // 目标=sum(w_i*(q_i-p_i)^2)+lambda*sum(两段斜率不一致的残差^2)。
  // 低置信候选更服从邻点连续性；这是供贴面取邻域使用的前置平滑。
  std::vector<Candidate> regularizedPath = path;
  // 对每个横向坐标解一个置信度加权的平滑最小二乘问题：
  //   数据残差 + lambda * 相邻两段斜率差。
  // 二阶项不会把连续陡坡当成跳变，端点也不会像缩短窗口均值那样突然变平。
  if (regularizedPath.size() >= 3 && config.curveRegularizationWeight > 0.0f) {
    const int longitudinalAxis = static_cast<int>(config.longitudinalAxis);
    const int pointCount = static_cast<int>(regularizedPath.size());
    for (int axis = 0; axis < 3; ++axis) {
      if (axis == longitudinalAxis) {
        continue;
      }
      // The second-difference penalty is pentadiagonal, not a dense system.
      Eigen::SparseMatrix<float> system(pointCount, pointCount);
      std::vector<Eigen::Triplet<float>> entries;
      entries.reserve(10 * pointCount);
      Eigen::VectorXf rhs = Eigen::VectorXf::Zero(pointCount);
      for (int i = 0; i < pointCount; ++i) {
        const Candidate &candidate =
            regularizedPath[static_cast<std::size_t>(i)];
        const float dataWeight =
            std::max(0.05f, 4.0f * candidate.feature * candidate.feature);
        entries.emplace_back(i, i, dataWeight);
        rhs(i) += dataWeight * candidate.point.data[axis];
      }
      for (int i = 1; i + 1 < pointCount; ++i) {
        const float previousStation =
            regularizedPath[static_cast<std::size_t>(i - 1)]
                .point.data[longitudinalAxis];
        const float currentStation =
            regularizedPath[static_cast<std::size_t>(i)]
                .point.data[longitudinalAxis];
        const float nextStation =
            regularizedPath[static_cast<std::size_t>(i + 1)]
                .point.data[longitudinalAxis];
        const float previousGap =
            std::max(1.0e-4f, currentStation - previousStation);
        const float nextGap = std::max(1.0e-4f, nextStation - currentStation);
        const float totalGap = previousGap + nextGap;
        Eigen::Vector3f coefficients;
        coefficients << -nextGap / totalGap, 1.0f, -previousGap / totalGap;
        const int indices[3] = {i - 1, i, i + 1};
        for (int row = 0; row < 3; ++row) {
          for (int column = 0; column < 3; ++column) {
            entries.emplace_back(indices[row], indices[column],
                config.curveRegularizationWeight * coefficients[row] * coefficients[column]);
          }
        }
      }
      system.setFromTriplets(entries.begin(), entries.end());
      Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> solver;
      solver.compute(system);
      if (solver.info() != Eigen::Success) continue;
      const Eigen::VectorXf solution = solver.solve(rhs);
      if (!solution.allFinite()) {
        continue;
      }
      for (int i = 0; i < pointCount; ++i) {
        regularizedPath[static_cast<std::size_t>(i)].point.data[axis] =
            solution[i];
      }
    }
  }
  return regularizedPath;
}

struct WallEvidence {
  std::vector<PointT> patch, inliers;
  float beforeRms = -1;
  std::string method = "profile_fallback";
};

bool isSteepWallTrack(const std::vector<Candidate> &path) {
  // “贴垂直侧壁”按整段可靠侧壁的中值方向启用，避免单片角度在 45°
  // 附近摆动时反复切换约束。浅夹角焊道继续采用原有定位流程。
  std::vector<float> wallSteepness;
  for (const auto &c : path)
    if (c.hasWall)
      wallSteepness.push_back(std::abs(c.wallNormal.x()));
  return wallSteepness.size() >= 3 && median(wallSteepness) >= std::sqrt(.5f);
}

void fitLocalWallModel(Candidate &candidate,
                       const std::vector<Candidate> &regularizedPath,
                       std::size_t i, int side, const CloudT::ConstPtr &input,
                       const DetectorConfig &config, WallEvidence &evidence,
                       DebugTrace *debug) {
  // 用相邻站位的原始侧壁小片拟合三维局部面，避免逐片直线噪声导致贴壁后抖动。
  // 只收集候选外侧、靠近既有侧壁的点，不把槽内底面加入此拟合。
  if (candidate.hasWall) {
    const auto axes = transverseAxes(config.longitudinalAxis);
    const float profileScale = profileGeometryScale(config);
    const int a = int(config.longitudinalAxis);
    const auto anchor = toVector(candidate.point);
    Eigen::Vector3f trend = Eigen::Vector3f::Zero();
    if (regularizedPath.size() > 1) {
      const auto
          first = toVector(regularizedPath[i ? i - 1 : i].point),
          last = toVector(
              regularizedPath[std::min(i + 1, regularizedPath.size() - 1)]
                  .point);
      if (last[a] > first[a])
        trend = (last - first) / (last[a] - first[a]);
    }
    std::vector<Eigen::Vector3d> patch;
    for (const auto &sample : *input) {
      if (!pcl::isFinite(sample))
        continue;
      const Eigen::Vector3f delta = toVector(sample) - anchor;
      if (std::abs(delta[a]) >
          std::max(config.sliceThickness, 2 * config.sliceStep))
        continue;
      const auto d = (delta - trend * delta[a]).eval();
      const float du = d[axes[0]], dv = d[axes[1]];
      if ((side == 0 ? du >= -.5f * profileScale
                     : du <= .5f * profileScale) ||
          std::abs(du) > 6 * profileScale ||
          std::abs(dv) > 12 * profileScale)
        continue;
      const Eigen::Vector2f uv(sample.data[axes[0]] - trend[axes[0]] * delta[a],
                               sample.data[axes[1]] -
                                   trend[axes[1]] * delta[a]);
      if (std::abs(candidate.wallNormal.dot(uv) - candidate.wallOffset) >
          std::max(.5f, profileScale))
        continue;
      patch.push_back(toVector(sample).cast<double>());
    }
    if (debug)
      for (const auto &p : patch)
        evidence.patch.push_back(toPoint(p.cast<float>()));
    for (int pass = 0; pass < 2 && patch.size() >= 20; ++pass) {
      Eigen::Vector3d center = Eigen::Vector3d::Zero();
      for (const auto &p : patch)
        center += p;
      center /= double(patch.size());
      Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
      for (const auto &p : patch) {
        const auto d = (p - center).eval();
        cov += d * d.transpose();
      }
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(cov /
                                                         double(patch.size()));
      if (eig.info() != Eigen::Success || eig.eigenvalues()[1] < .01)
        break;
      const Eigen::Vector3d n = eig.eigenvectors().col(0);
      const float rms = float(std::sqrt(std::max(0., eig.eigenvalues()[0])));
      if (pass == 0) {
        std::vector<Eigen::Vector3d> kept;
        for (const auto &p : patch)
          if (std::abs(n.dot(p - center)) <= std::max(.15f, 2.5f * rms))
            kept.push_back(p);
        patch.swap(kept);
        continue;
      }
      Eigen::Vector2f normal{float(n[axes[0]]), float(n[axes[1]])};
      const float length = normal.norm();
      if (length < .1f || rms > std::min(.6f, config.regionGrowingDistance))
        break;
      normal /= length;
      if (std::abs(normal.dot(candidate.wallNormal)) <
          std::cos(12.f * kPi / 180.f))
        break;
      candidate.wallNormal = normal;
      candidate.wallOffset = float((n.dot(center) - n[a] * anchor[a]) / length);
      candidate.wallRms = rms;
      evidence.method = "local_plane";
      if (debug)
        for (const auto &p : patch)
          evidence.inliers.push_back(toPoint(p.cast<float>()));
    }
  }
}

void refitWallModel(Candidate &candidate,
                    const std::vector<Candidate> &regularizedPath,
                    const std::vector<Candidate> &wallReferences, std::size_t i,
                    int side, const CloudT::ConstPtr &input,
                    const DetectorConfig &config, WallEvidence &evidence,
                    DebugTrace *debug) {
  // S13 高残差侧壁重拟合，解决底面/外侧平板混入造成的侧壁偏斜。
  // 只处理 RMS>min(0.25,regionGrowingDistance) 的点；可靠段保持原模型。
  // 缩小纵向窗口、放宽旧面距离门槛，以固定种子采样 240 个三点平面。
  // 用可靠邻站位方向约束选平面，再对距离内点做两次 PCA 精修。
  // 至少 20 个内点且面内支撑不退化才采用；否则仍使用 S12 的模型。
  // 断点观察：reference、patch、best、rms、wallEvidence[i].method。
  const float reliableRms = std::min(.25f, config.regionGrowingDistance);
  if (!candidate.hasWall || candidate.wallRms <= reliableRms)
    return;
  const auto axes = transverseAxes(config.longitudinalAxis);
  const float profileScale = profileGeometryScale(config);
  const int a = int(config.longitudinalAxis);
  const Eigen::Vector3f anchor = toVector(candidate.point);
  std::vector<float> angles;
  for (std::size_t j = i > 12 ? i - 12 : 0;
       j < std::min(i + 13, wallReferences.size()); ++j) {
    const auto &ref = wallReferences[j];
    if (!ref.hasWall || ref.wallRms > reliableRms ||
        std::abs(coordinate(ref.point, config.longitudinalAxis) - anchor[a]) >
            12 * config.sliceStep)
      continue;
    Eigen::Vector2f n = ref.wallNormal;
    if (n.x() < 0)
      n = -n;
    angles.push_back(std::atan2(n.y(), n.x()));
  }
  Eigen::Vector2f reference = candidate.wallNormal;
  const float maxAngle = angles.size() >= 3 ? 12.f : 35.f;
  if (angles.size() >= 3) {
    const float angle = median(angles);
    reference = {std::cos(angle), std::sin(angle)};
  }
  const float orientationLimit = std::cos(maxAngle * kPi / 180.f);
  Eigen::Vector3f trend = Eigen::Vector3f::Zero();
  if (regularizedPath.size() > 1) {
    const Eigen::Vector3f delta =
        toVector(regularizedPath[std::min(i + 1, regularizedPath.size() - 1)]
                     .point) -
        toVector(regularizedPath[i ? i - 1 : i].point);
    if (delta[a] > 0)
      trend = delta / delta[a];
  }
  // 缩小纵向窗口，放宽旧平面的距离门槛，让实测侧壁能纠正旧方向。
  std::vector<Eigen::Vector3d> patch;
  for (const auto &sample : *input) {
    if (!pcl::isFinite(sample))
      continue;
    const Eigen::Vector3f delta = toVector(sample) - anchor;
    if (std::abs(delta[a]) >
        std::max(config.sliceStep, .5f * config.sliceThickness))
      continue;
    const Eigen::Vector3f d = delta - trend * delta[a];
    if (d[axes[0]] * (side == 0 ? 1.f : -1.f) >= -.5f * profileScale ||
        std::abs(d[axes[0]]) > 6 * profileScale ||
        std::abs(d[axes[1]]) > 12 * profileScale)
      continue;
    const Eigen::Vector2f uv(anchor[axes[0]] + d[axes[0]],
                             anchor[axes[1]] + d[axes[1]]);
    if (std::abs(candidate.wallNormal.dot(uv) - candidate.wallOffset) >
        3 * profileScale)
      continue;
    patch.push_back(delta.cast<double>());
  }
  if (patch.size() < 20)
    return;
  const double tolerance = std::min(
      double(reliableRms), std::max(.08, double(profileScale) * .15));
  auto oriented = [&](const Eigen::Vector3d &n) {
    Eigen::Vector2f uv{float(n[axes[0]]), float(n[axes[1]])};
    if (uv.norm() < .1f)
      return false;
    uv.normalize();
    return std::abs(uv.x()) >= std::sqrt(.5f) &&
           std::abs(uv.dot(reference)) >= orientationLimit;
  };
  std::mt19937 random(47);
  std::uniform_int_distribution<std::size_t> pick(0, patch.size() - 1);
  std::vector<Eigen::Vector3d> best;
  for (int trial = 0; trial < 240; ++trial) {
    const Eigen::Vector3d p = patch[pick(random)], q = patch[pick(random)],
                          r = patch[pick(random)];
    Eigen::Vector3d n = (q - p).cross(r - p);
    if (n.norm() < 1e-6)
      continue;
    n.normalize();
    if (!oriented(n))
      continue;
    std::vector<Eigen::Vector3d> inliers;
    for (const auto &s : patch)
      if (std::abs(n.dot(s - p)) <= tolerance)
        inliers.push_back(s);
    if (inliers.size() > best.size())
      best.swap(inliers);
  }
  for (int pass = 0; pass < 2 && best.size() >= 20; ++pass) {
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    for (const auto &p : best)
      center += p;
    center /= double(best.size());
    Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
    for (const auto &p : best) {
      const Eigen::Vector3d delta = p - center;
      covariance += delta * delta.transpose();
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(covariance /
                                                       double(best.size()));
    if (eig.info() != Eigen::Success || eig.eigenvalues()[1] < .01)
      break;
    const Eigen::Vector3d normal = eig.eigenvectors().col(0);
    const float rms = float(std::sqrt(std::max(0., eig.eigenvalues()[0])));
    if (!oriented(normal) || rms > reliableRms)
      break;
    if (pass == 0) {
      best.clear();
      for (const auto &p : patch)
        if (std::abs(normal.dot(p - center)) <= tolerance)
          best.push_back(p);
      continue;
    }
    Eigen::Vector2f uv{float(normal[axes[0]]), float(normal[axes[1]])};
    const float length = uv.norm();
    uv /= length;
    candidate.wallNormal = uv;
    candidate.wallOffset =
        float(normal.dot(center) / length) +
        uv.dot(Eigen::Vector2f(anchor[axes[0]], anchor[axes[1]]));
    candidate.wallRms = rms;
    evidence.method = "robust_neighborhood";
    if (debug) {
      evidence.patch.clear();
      evidence.inliers.clear();
      for (const auto &p : patch)
        evidence.patch.push_back(
            toPoint((p + anchor.cast<double>()).cast<float>()));
      for (const auto &p : best)
        evidence.inliers.push_back(
            toPoint((p + anchor.cast<double>()).cast<float>()));
    }
  }
}

std::vector<Candidate>
fitWallModels(const std::vector<Candidate> &path,
              const std::vector<Candidate> &regularizedPath, int side,
              const CloudT::ConstPtr &input, const DetectorConfig &config,
              std::vector<WallEvidence> &wallEvidence, DebugTrace *debug,
              const std::function<CloudT::ConstPtr(std::size_t)> &neighborhood) {
  const bool steepWallTrack = isSteepWallTrack(path);
  auto wallModels = regularizedPath;
  wallEvidence.resize(regularizedPath.size());
  // S12 先完成全部局部侧壁拟合；S13 再用这份固定快照辅助高残差重拟合。
  for (std::size_t i = 0; i < wallModels.size(); ++i) {
    wallModels[i].hasWall = wallModels[i].hasWall && steepWallTrack;
    fitLocalWallModel(wallModels[i], regularizedPath, i, side, neighborhood(i), config,
                      wallEvidence[i], debug);
    wallEvidence[i].beforeRms = wallModels[i].wallRms;
    if (stepDebugEnabled(debug))
      publishDebugStep(debug,
                       sideStepTitle("S12 局部侧壁拟合", side,
                                     wallModels[i].point,
                                     config.longitudinalAxis),
                       {{"侧壁取点", pointsCloud(wallEvidence[i].patch)},
                        {"侧壁内点", pointsCloud(wallEvidence[i].inliers)}},
                       wallEvidence[i].method);
  }
  const auto wallReferences = wallModels;
  for (std::size_t i = 0; i < wallModels.size(); ++i) {
    refitWallModel(wallModels[i], regularizedPath, wallReferences, i, side,
                   neighborhood(i), config, wallEvidence[i], debug);
    if (stepDebugEnabled(debug))
      publishDebugStep(debug,
                       sideStepTitle("S13 高残差侧壁重拟合", side,
                                     wallModels[i].point,
                                     config.longitudinalAxis),
                       {{"侧壁取点", pointsCloud(wallEvidence[i].patch)},
                        {"侧壁内点", pointsCloud(wallEvidence[i].inliers)}},
                       wallEvidence[i].method);
  }
  return wallModels;
}

DebugWallPoint projectOntoWall(Candidate &candidate, int side,
                               const CloudT::ConstPtr &input,
                               const DetectorConfig &config,
                               WallEvidence &evidence, DebugTrace *debug) {
  const float halfThickness = .5f * config.sliceThickness;
  // S14 将正则后的候选投影到当前站位的侧壁线：p'=p-n*(n·p-offset)。
  // 位移、侧壁 RMS、到当前原始截面最近点距离都通过才采用。
  // beforeFloor 是贴壁结果，也是下一阶段失败时保留的位置。
  DebugWallPoint wall;
  wall.side = side;
  wall.before = candidate.point;
  wall.after = candidate.point;
  wall.wallFitMethod = evidence.method;
  wall.wallRmsBeforeRefit = evidence.beforeRms;
  if (debug) {
    wall.wallPatch = std::move(evidence.patch);
    wall.wallInliers = std::move(evidence.inliers);
  }
  wall.normal = candidate.wallNormal;
  wall.offset = candidate.wallOffset;
  wall.rms = candidate.wallRms;
  wall.status = candidate.hasWall ? "wall_fit_rejected" : "no_reliable_wall";
  if (debug)
    debug->beforeWallPath[side].push_back(candidate.point);
  if (candidate.hasWall &&
      candidate.wallRms <= std::min(.6f, config.regionGrowingDistance)) {
    const auto axes = transverseAxes(config.longitudinalAxis);
    Eigen::Vector2f p(candidate.point.data[axes[0]],
                      candidate.point.data[axes[1]]);
    const float signedDistance =
        candidate.wallNormal.dot(p) - candidate.wallOffset;
    wall.distanceBefore = wall.distanceAfter = std::abs(signedDistance);
    const Eigen::Vector2f projected = p - signedDistance * candidate.wallNormal;
    const float limit = std::min(config.surfaceIntersectionTolerance,
                                 3.f * profileGeometryScale(config));
    wall.status = "displacement_rejected";
    if (projected.allFinite() && std::abs(signedDistance) <= limit) {
      float nearest = std::numeric_limits<float>::infinity();
      for (const auto &sample : *input) {
        if (!pcl::isFinite(sample) ||
            std::abs(coordinate(sample, config.longitudinalAxis) -
                     coordinate(candidate.point, config.longitudinalAxis)) >
                halfThickness)
          continue;
        nearest = std::min(nearest, (Eigen::Vector2f(sample.data[axes[0]],
                                                     sample.data[axes[1]]) -
                                     projected)
                                        .norm());
      }
      wall.supportDistance = nearest;
      wall.status = "support_rejected";
      if (nearest <= std::max(.5f, profileGeometryScale(config))) {
        candidate.point.data[axes[0]] = projected.x();
        candidate.point.data[axes[1]] = projected.y();
        wall.after = candidate.point;
        wall.distanceAfter = std::abs(candidate.wallNormal.dot(projected) -
                                      candidate.wallOffset);
        wall.status = "projected";
      }
    }
  }
  if (stepDebugEnabled(debug))
    publishDebugStep(debug,
                     sideStepTitle("S14 贴壁投影", side, candidate.point,
                                   config.longitudinalAxis),
                     {{"侧壁取点", pointsCloud(wall.wallPatch)},
                      {"修正前", pointsCloud({wall.before})},
                      {"贴壁结果", pointsCloud({wall.after})}},
                     wall.status);
  return wall;
}

void intersectWallAndFloor(Candidate &candidate, DebugWallPoint &wall,
                           const CloudT::ConstPtr &input,
                           const DetectorConfig &config) {
  const auto axes = transverseAxes(config.longitudinalAxis);
  const int a = int(config.longitudinalAxis);
  const auto anchor = toVector(candidate.point);
  const float halfThickness = .5f * config.sliceThickness;
  const Eigen::Vector2f n(wall.floorNormal[axes[0]], wall.floorNormal[axes[1]]);
  const float determinant =
      candidate.wallNormal.x() * n.y() - candidate.wallNormal.y() * n.x();
  wall.floorStatus = "floor_wall_parallel";
  if (std::abs(determinant) < .2f * n.norm())
    return;
  // S16 固定 A，联立侧壁线与底面在此截面内的直线，解 2×2 方程。
  // wallNu*U+wallNv*V=wallOffset
  // floorNu*U+floorNv*V=floorOffset-floorNa*A
  // 行列式过小表示近乎平行，拒绝；位移过大或没有原始点支持也拒绝。
  // 此阶段成功点同时落在两个面上；独立的最终轨迹平滑可离开局部拟合面。
  const float rhs = wall.floorOffset - wall.floorNormal[a] * anchor[a];
  Eigen::Vector3f point = anchor;
  point[axes[0]] =
      (candidate.wallOffset * n.y() - candidate.wallNormal.y() * rhs) /
      determinant;
  point[axes[1]] =
      (candidate.wallNormal.x() * rhs - candidate.wallOffset * n.x()) /
      determinant;
  wall.floorShift = (point - anchor).norm();
  wall.floorStatus = "floor_shift_rejected";
  if (!point.allFinite() || wall.floorShift > config.floorPatchMaxCorrection)
    return;
  float nearest = std::numeric_limits<float>::infinity();
  for (const auto &sample : *input) {
    if (!pcl::isFinite(sample) ||
        std::abs(sample.data[a] - anchor[a]) > halfThickness)
      continue;
    const auto delta = (toVector(sample) - point).eval();
    nearest = std::min(nearest, std::hypot(delta[axes[0]], delta[axes[1]]));
  }
  wall.floorStatus = "floor_intersection_unsupported";
  if (nearest > std::max(.5f, profileGeometryScale(config)))
    return;
  candidate.point = toPoint(point);
  wall.after = candidate.point;
  wall.supportDistance = nearest;
  wall.distanceAfter = std::abs(candidate.wallNormal.dot(Eigen::Vector2f(
                                    point[axes[0]], point[axes[1]])) -
                                candidate.wallOffset);
  wall.floorResidual = std::abs(wall.floorNormal.dot(point) - wall.floorOffset);
  wall.floorStatus = "floor_wall_intersection";
  wall.status = "projected";
}

void constrainToLocalFloor(Candidate &candidate,
                           const std::vector<Candidate> &regularizedPath,
                           std::size_t i, int side,
                           const CloudT::ConstPtr &input,
                           const DetectorConfig &config, DebugWallPoint &wall,
                           DebugTrace *debug) {
  // S15 局部底面：以轨迹切线 along 和向槽内 across 建立正方形取点区域。
  // 边长 floorPatchSize，靠壁边向内退让 floorPatchGuard；仅采集原始 input。
  // 两次 PCA 剔除离群点，至少 12 个内点，第二特征值>=0.1，RMS 不超上限。
  // floorNormal·p=floorOffset 是局部拟合面，不强制底面平坦或全局水平。
  // 断点观察：patchCenter、patch.size()、floorRms、floorStatus。
  wall.beforeFloor = candidate.point;
  wall.floorStatus =
      config.useLocalFloorConstraint ? "wall_unreliable" : "disabled";
  if (config.useLocalFloorConstraint && candidate.hasWall &&
      candidate.wallRms <= std::min(.6f, config.regionGrowingDistance)) {
    const auto axes = transverseAxes(config.longitudinalAxis);
    const int a = int(config.longitudinalAxis);
    const auto anchor = toVector(candidate.point);
    const float size = config.floorPatchSize, half = .5f * size;
    Eigen::Vector3f across = Eigen::Vector3f::Zero();
    across[axes[0]] = candidate.floorDirection.x();
    across[axes[1]] = candidate.floorDirection.y();
    if (across[axes[0]] * (side == 0 ? 1.f : -1.f) < 0)
      across = -across;
    Eigen::Vector3f along = Eigen::Vector3f::Zero();
    along[a] = 1;
    if (regularizedPath.size() > 1) {
      const auto
          first = toVector(regularizedPath[i ? i - 1 : i].point),
          last = toVector(
              regularizedPath[std::min(i + 1, regularizedPath.size() - 1)]
                  .point);
      if ((last - first).norm() > 1e-4f)
        along = (last - first).normalized();
    }
    along -= across * along.dot(across);
    along.normalize();
    const Eigen::Vector3f initialNormal = along.cross(across).normalized();
    const Eigen::Vector3f patchCenter =
        anchor + across * (half + config.floorPatchGuard);
    std::vector<Eigen::Vector3d> patch;
    for (const auto &sample : *input) {
      if (!pcl::isFinite(sample))
        continue;
      const Eigen::Vector3f delta = toVector(sample) - patchCenter;
      if (std::abs(delta.dot(along)) > half ||
          std::abs(delta.dot(across)) > half ||
          std::abs(delta.dot(initialNormal)) > std::max(2.f, .5f * size))
        continue;
      patch.push_back(toVector(sample).cast<double>());
      if (debug)
        wall.floorPatch.push_back(sample);
    }
    if (stepDebugEnabled(debug))
      publishDebugStep(debug,
                       sideStepTitle("S15 正方形底面取点", side,
                                     candidate.point, config.longitudinalAxis),
                       {{"底面取点", pointsCloud(wall.floorPatch)},
                        {"当前焊缝点", pointsCloud({candidate.point})}});
    wall.floorStatus = "insufficient_floor_points";
    for (int pass = 0; pass < 2 && patch.size() >= 12; ++pass) {
      Eigen::Vector3d center = Eigen::Vector3d::Zero();
      for (const auto &p : patch)
        center += p;
      center /= double(patch.size());
      Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
      for (const auto &p : patch) {
        const auto d = (p - center).eval();
        cov += d * d.transpose();
      }
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(cov /
                                                         double(patch.size()));
      wall.floorStatus = "floor_degenerate";
      if (eig.info() != Eigen::Success || eig.eigenvalues()[1] < .1)
        break;
      const Eigen::Vector3d normal = eig.eigenvectors().col(0);
      const float rms = float(std::sqrt(std::max(0., eig.eigenvalues()[0])));
      if (pass == 0) {
        std::vector<Eigen::Vector3d> kept;
        for (const auto &p : patch)
          if (std::abs(normal.dot(p - center)) <= std::max(.15f, 2.5f * rms))
            kept.push_back(p);
        patch.swap(kept);
        continue;
      }
      wall.floorNormal = normal.cast<float>();
      wall.floorOffset = float(normal.dot(center));
      wall.floorRms = rms;
      wall.floorInliers = int(patch.size());
      for (const auto &q : std::array<Eigen::Vector2f, 4>{
               {{-half, -half}, {half, -half}, {half, half}, {-half, half}}}) {
        Eigen::Vector3f p = patchCenter + along * q.x() + across * q.y();
        p -= wall.floorNormal * (wall.floorNormal.dot(p) - wall.floorOffset);
        wall.floorCorners.push_back(toPoint(p));
      }
      if (stepDebugEnabled(debug))
        publishDebugStep(debug,
                         sideStepTitle("S15 局部底面拟合", side,
                                       candidate.point,
                                       config.longitudinalAxis),
                         {{"底面取点", pointsCloud(wall.floorPatch)},
                          {"拟合面四角", pointsCloud(wall.floorCorners)}});
      wall.floorStatus = "floor_fit_rejected";
      // 局部底面可以与跨槽参考面方向不同；由残差、夹角和实测支持判定。
      if (rms > config.floorPatchMaxRms)
        break;
      intersectWallAndFloor(candidate, wall, input, config);
    }
  }
}

void recordWallRefinement(const Candidate &candidate, int side,
                          const DetectorConfig &config,
                          const DebugWallPoint &wall, DebugTrace *debug) {
  if (debug) {
    DebugStep step;
    step.title =
        std::string("局部底面 ") + (side == 0 ? "L" : "R") + " · " +
        std::to_string(coordinate(candidate.point, config.longitudinalAxis));
    step.status = "completed";
    step.message = wall.floorStatus;
    step.metrics = {{"wall_rms_before_refit", wall.wallRmsBeforeRefit},
                    {"wall_rms", wall.rms},
                    {"wall_patch_points", double(wall.wallPatch.size())},
                    {"wall_inliers", double(wall.wallInliers.size())},
                    {"patch_size", config.floorPatchSize},
                    {"patch_points", double(wall.floorPatch.size())},
                    {"inliers", double(wall.floorInliers)},
                    {"floor_rms", wall.floorRms},
                    {"floor_shift", wall.floorShift},
                    {"floor_residual", wall.floorResidual},
                    {"wall_residual", wall.distanceAfter}};
    CloudT::Ptr patch(new CloudT), plane(new CloudT), before(new CloudT),
        after(new CloudT);
    for (const auto &p : wall.floorPatch)
      patch->push_back(p);
    if (wall.floorCorners.size() == 4)
      for (int row = 0; row <= 10; ++row)
        for (int col = 0; col <= 10; ++col) {
          const Eigen::Vector3f p = toVector(wall.floorCorners[0]) +
                                    (toVector(wall.floorCorners[1]) -
                                     toVector(wall.floorCorners[0])) *
                                        (row / 10.f) +
                                    (toVector(wall.floorCorners[3]) -
                                     toVector(wall.floorCorners[0])) *
                                        (col / 10.f);
          plane->push_back(toPoint(p));
        }
    before->push_back(wall.beforeFloor);
    after->push_back(wall.after);
    CloudT::Ptr wallPatch(new CloudT), wallInliers(new CloudT);
    for (const auto &p : wall.wallPatch)
      wallPatch->push_back(p);
    for (const auto &p : wall.wallInliers)
      wallInliers->push_back(p);
    step.layers = {{"侧壁取点", wallPatch},   {"侧壁内点", wallInliers},
                   {"正方形底面取点", patch}, {"拟合底面", plane},
                   {"修正前", before},        {"最终点", after}};
    if (stepDebugEnabled(debug))
      publishDebugStep(debug,
                       sideStepTitle("S16 底面/侧壁/站位求交", side,
                                     candidate.point, config.longitudinalAxis),
                       step.layers, wall.floorStatus);
    debug->steps.push_back(std::move(step));
    debug->wallPoints.push_back(wall);
  }
}

float candidateConfidence(const Candidate &sourceCandidate, int side,
                          const DetectorConfig &config,
                          std::size_t &historyDominated) {
  float pointScore = sourceCandidate.feature;
  if (sourceCandidate.historyAvailable[static_cast<std::size_t>(side)]) {
    const float historyAgreement = std::exp(
        -sourceCandidate.historyDistance[static_cast<std::size_t>(side)] /
        std::max(1.0f, config.historySearchRadius));
    pointScore = 0.65f * sourceCandidate.feature + 0.35f * historyAgreement;
    if (sourceCandidate.feature < 0.35f && historyAgreement > 0.55f) {
      ++historyDominated;
    }
  }
  return clamp01(pointScore);
}

void fillSideResult(const std::vector<Candidate> &path, int side,
                    const CloudT::ConstPtr &input, const DetectorConfig &config,
                    DebugTrace *debug, CloudT::Ptr seam,
                    std::vector<float> *pointConfidence,
                    std::vector<SeamPointQuality> *pointQuality,
                    float *confidence, float *historyUsageRatio) {
  const auto regularizedPath = regularizePath(path, config);
  const LongitudinalIndex inputIndex(input, config.longitudinalAxis);
  const float halfFloor = .5f * config.floorPatchSize;
  const float localRadius = std::max({config.sliceThickness, 2 * config.sliceStep,
                                     halfFloor + std::max(2.f, halfFloor)});
  auto neighborhood = [&](std::size_t i) {
    return inputIndex.around(coordinate(regularizedPath[i].point, config.longitudinalAxis), localRadius);
  };
  if (stepDebugEnabled(debug))
    publishDebugStep(debug,
                     std::string("S11 曲线正则 · ") + (side == 0 ? "L" : "R"),
                     {{"原候选", candidateCloud(path)},
                      {"平滑路径", candidateCloud(regularizedPath)}});
  seam->reserve(regularizedPath.size());
  pointConfidence->reserve(path.size());
  pointQuality->reserve(path.size());
  std::vector<WallEvidence> wallEvidence;
  const auto wallModels = config.useSurfaceRefinement
      ? fitWallModels(path, regularizedPath, side, input,
                      config, wallEvidence, debug, neighborhood)
      : regularizedPath;
  std::size_t historyDominated = 0;
  for (std::size_t i = 0; i < wallModels.size(); ++i) {
    Candidate candidate = wallModels[i];
    DebugWallPoint wall;
    wall.before = wall.after = wall.beforeFloor = candidate.point;
    wall.rms = -1.0f;
    wall.status = "disabled";
    if (config.useSurfaceRefinement) {
      const auto local = neighborhood(i);
      wall = projectOntoWall(candidate, side, local, config, wallEvidence[i], debug);
      constrainToLocalFloor(candidate, regularizedPath, i, side, local, config,
                            wall, debug);
      recordWallRefinement(candidate, side, config, wall, debug);
    } else if (debug) {
      debug->beforeWallPath[side].push_back(candidate.point);
    }
    seam->push_back(candidate.point);
    const float sourceConfidence =
        candidateConfidence(path[i], side, config, historyDominated);
    pointConfidence->push_back(sourceConfidence);
    SeamPointQuality quality;
    if (path[i].curvatureFallback) quality.flags |= QualityCurvatureFallback;
    quality.sourceConfidence = sourceConfidence;
    quality.wallRms = wall.rms;
    quality.floorRms = wall.floorRms;
    quality.supportDistance = wall.supportDistance;
    quality.geometryCorrection =
        (toVector(wall.after) - toVector(wall.before)).norm();
    quality.wallAccepted = wall.status == "projected";
    quality.wallNormal = wall.normal;
    quality.wallOffset = wall.offset;
    quality.trackedWall=path[i].trackedWall;
    quality.trackedNormal=path[i].trackedNormal;
    quality.trackedOffset=path[i].trackedOffset;
    quality.floorAccepted = wall.floorStatus == "floor_wall_intersection";
    quality.floorNormal = wall.floorNormal;
    quality.floorOffset = wall.floorOffset;
    if (path[i].feature == 0.08f)
      quality.flags |= QualityInterpolated;
    if (!quality.wallAccepted)
      quality.flags |= QualityWallUnavailable;
    if (config.useSurfaceRefinement && config.useLocalFloorConstraint && !quality.floorAccepted)
      quality.flags |= QualityFloorUnavailable;
    pointQuality->push_back(quality);
  }
  seam->width = static_cast<std::uint32_t>(seam->size());
  seam->height = 1;
  seam->is_dense = true;
  *confidence =
      std::accumulate(pointConfidence->begin(), pointConfidence->end(), 0.0f) /
      static_cast<float>(pointConfidence->size());
  *historyUsageRatio = static_cast<float>(historyDominated) /
                       static_cast<float>(regularizedPath.size());
}

// Short internal omissions inherit the wall between two measured anchors,
// never extrapolate it beyond either end or across an unsupported long run.
void bridgeTrackedWalls(std::vector<Candidate> &path,const DetectorConfig &config) {
  for (std::size_t begin=0;begin<path.size();) {
    if (path[begin].trackedWall) { ++begin; continue; }
    std::size_t end=begin;
    while (end<path.size() && !path[end].trackedWall) ++end;
    if (begin>0 && end<path.size()) {
      const auto &a=path[begin-1]; const auto &b=path[end];
      const float sa=coordinate(a.point,config.longitudinalAxis), sb=coordinate(b.point,config.longitudinalAxis);
      if (sb>sa && (sb-sa<=(config.maxInterpolatedGapStations+1)*config.sliceStep ||
          (config.useProfileCurvature && sb-sa<=2*profileGeometryScale(config) &&
           std::fabs(a.trackedNormal.dot(b.trackedNormal))>=.94f))) {
        const float sign=a.trackedNormal.dot(b.trackedNormal)<0 ? -1.f : 1.f;
        for (auto j=begin;j<end;++j) {
          const float t=(coordinate(path[j].point,config.longitudinalAxis)-sa)/(sb-sa);
          const Eigen::Vector2f normal=(1-t)*a.trackedNormal+t*sign*b.trackedNormal;
          if (normal.norm()<.9f) continue;
          path[j].trackedWall=true; path[j].trackedNormal=normal.normalized();
          path[j].trackedOffset=((1-t)*a.trackedOffset+t*sign*b.trackedOffset)/normal.norm();
        }
      }
    }
    begin=end;
  }
}

void summarizeDetection(SeamResult &result,
                        const std::vector<Candidate> &leftPath,
                        const std::vector<Candidate> &rightPath) {
  result.confidence = 0.5f * (result.leftConfidence + result.rightConfidence);

  std::size_t regionEdgeStations = 0;
  std::size_t localContourStations = 0;
  for (std::size_t i = 0; i < leftPath.size() && i < rightPath.size(); ++i) {
    if (leftPath[i].fromRegionEdges && rightPath[i].fromRegionEdges) {
      ++regionEdgeStations;
    }
    if (leftPath[i].fromLocalCorner || rightPath[i].fromLocalCorner)
      ++localContourStations;
  }
  std::ostringstream message;
  message << "Detected left=" << result.leftSeam->size()
          << " and right=" << result.rightSeam->size()
          << " seam stations; confidence=" << result.confidence
           << "; region-edge stations=" << regionEdgeStations
           << "; local-contour stations=" << localContourStations;
  if (result.leftSmoothing.applied || result.rightSmoothing.applied) {
    message << "; robust-rts mean correction L/R="
            << result.leftSmoothing.meanCorrection << "/"
            << result.rightSmoothing.meanCorrection
            << "; innovation outliers L/R="
            << result.leftSmoothing.innovationOutliers << "/"
            << result.rightSmoothing.innovationOutliers;
  }
  result.summary = message.str();
}

} // namespace

CloudT::Ptr SeamDetector::preprocess(const CloudT::ConstPtr &input,
                                     const DetectorConfig &config,
                                     DebugTrace *debug) {
  // S01 去除 NaN/Inf。输出 finite
  // 不改变有效点坐标；先检查这里是否丢掉大量输入。
  CloudT::Ptr finite(new CloudT);
  finite->reserve(input->size());
  for (const PointT &point : *input) {
    if (pcl::isFinite(point)) {
      finite->push_back(point);
    }
  }
  finite->width = static_cast<std::uint32_t>(finite->size());
  finite->height = 1;
  finite->is_dense = true;
  if (debug) {
    debug->clouds.push_back({"01_finite", finite});
    publishDebugStep(debug, "S01 去除无效点", {{"处理后点云", finite}});
  }
  if (finite->empty()) {
    return finite;
  }

#ifdef WELD_NO_PCL_FILTERS
  // 仅供没有 ABI 兼容 PCL 库时运行合成单元测试；正式 CMake 构建不会定义它。
  (void)config;
  if (debug) {
    debug->clouds.push_back({"02_sor_bypassed_test_build", finite});
    debug->clouds.push_back({"03_voxel_bypassed_test_build", finite});
  }
  return finite;
#else
  // S02 SOR：每点与近邻的平均距离超过统计门槛时视为离群。
  // 默认 K=12、标准差倍数=2；有限点不足时原样通过。
  CloudT::Ptr denoised(new CloudT);
  if (config.outlierMeanK >= 3 &&
      finite->size() > static_cast<std::size_t>(config.outlierMeanK)) {
    pcl::StatisticalOutlierRemoval<PointT> filter;
    filter.setInputCloud(finite);
    filter.setMeanK(config.outlierMeanK);
    filter.setStddevMulThresh(config.outlierStddev);
    filter.filter(*denoised);
  } else {
    *denoised = *finite;
  }

  if (debug) {
    debug->clouds.push_back({"02_sor", denoised});
    publishDebugStep(debug, config.extractMainBody ? "S02 主体保护（跳过二次离群滤波）" : "S02 离群滤波",
                     {{"处理后点云", denoised}});
  }
  if (config.voxelSize <= 0.0f) {
    if (debug) {
      debug->clouds.push_back({"03_voxel_disabled", denoised});
      publishDebugStep(debug, "S03 体素降采样（已关闭）",
                       {{"保留点云", denoised}});
    }
    return denoised;
  }
  // S03 体素降采样：降低截面候选搜索成本；voxelSize<=0 时跳过。
  // 注意：后面的 S12/S13/S15/S16 实测支持检查仍使用完整 input。
  CloudT::Ptr downsampled(new CloudT);
  pcl::VoxelGrid<PointT> voxel;
  voxel.setInputCloud(denoised);
  voxel.setLeafSize(config.voxelSize, config.voxelSize, config.voxelSize);
  voxel.filter(*downsampled);
  if (debug) {
    debug->clouds.push_back({"03_voxel", downsampled});
    publishDebugStep(debug, "S03 体素降采样", {{"处理后点云", downsampled}});
  }
  return downsampled;
#endif
}

SeamResult SeamDetector::detect(const CloudT::ConstPtr &input,
                                const DetectorConfig &config,
                                const SeamHistory *leftHistory,
                                const SeamHistory *rightHistory,
                                DebugTrace *debug, const std::vector<ScanFrameRange> *frames,
                                bool generateSurfaceModels) const {
  SeamResult result;
  auto stageStarted = std::chrono::steady_clock::now();
  auto finishStage = [&](const std::string &name) {
    const auto now = std::chrono::steady_clock::now();
    result.stageMilliseconds[name] = std::chrono::duration<double, std::milli>(now - stageStarted).count();
    stageStarted = now;
  };
  const bool frameBased = frames && !frames->empty();
  result.frameBased = frameBased;
  result.inputFrameCount = frameBased ? frames->size() : 0;
  DetectorConfig workingConfig = config;
  if (config.uniformWidthEnabled) workingConfig.adaptiveSeamCount=false;
  if (workingConfig.adaptiveSeamCount) workingConfig.minSeamSeparation=std::min(config.minSeamSeparation,.5f);
  if (debug) {
    auto onStep = std::move(debug->onStep);
    *debug = DebugTrace();
    debug->onStep = std::move(onStep);
    debug->config = config;
    debug->algorithm = frameBased ? "scheme1_original_frames" : "reference";
    debug->leftHistoryPasses = leftHistory ? leftHistory->size() : 0;
    debug->rightHistoryPasses = rightHistory ? rightHistory->size() : 0;
    debug->clouds.push_back({"00_input", input});
    publishDebugStep(debug, "S00 输入点云", {{"原始点云", input}});
  }

  // S00-S03 校验与预处理。
  if (!validateDetectionInput(input, config, result.summary))
    return result;
  if (config.debugFrameBatchSize < 0) {
    result.summary = "Debug frame batch size must be nonnegative.";
    return result;
  }
  if (config.debugStartSlice < 1) {
    result.summary = "Debug start frame/section position must be positive.";
    return result;
  }
  // Validate original ranges BEFORE masking; otherwise a malformed frame map
  // could silently assign retained points to the wrong scan plane.
  if (frameBased) {
    std::size_t expected = 0;
    int previous = -1;
    for (const auto &frame : *frames) {
      if (frame.begin != expected || frame.count == 0 ||
          expected > input->size() || frame.count > input->size() - expected ||
          frame.frameIndex <= previous) {
        result.summary = "Original frame ranges are incomplete or unordered.";
        return result;
      }
      expected += frame.count;
      previous = frame.frameIndex;
    }
    if (expected != input->size()) {
      result.summary = "Original frames do not cover the complete point cloud.";
      return result;
    }
    for (const auto &point : *input)
      if (!pcl::isFinite(point)) {
        result.summary = "An original frame contains a non-finite point.";
        return result;
      }
  }
  finishStage("validation");
  const auto sampled=sampleMeasuredInput(input,config.inputVoxelSize,frames);
  if(!sampled.error.empty()) {result.summary=sampled.error;return result;}
  const auto &detectionInput=sampled.cloud;
  const auto *detectionFrames=frameBased ? &sampled.frames : nullptr;
  CloudT::ConstPtr sourceInput = detectionInput;
  result.inputSamplingPoints=detectionInput->size();
  result.inputSamplingRemovedPoints=input->size()-detectionInput->size();
  if(config.inputVoxelSize>0) {
    // Avoid a second, centroid-generating voxel pass on the sampled cloud.
    workingConfig.voxelSize=0;
    publishDebugStep(debug,"S00 入口实测点降采样",{{"采样后点云",detectionInput}},
        "逐帧保留实测点，主体、候选和表面提取共用此点云。",
        {{"input_points",double(input->size())},{"sampled_points",double(detectionInput->size())}});
  }
  finishStage("input_sampling");
  std::vector<ScanFrameRange> bodyFrames;
  if (config.extractMainBody) {
    const auto body = extractHullSurfaces(detectionInput, config);
    if (!body.ok()) {
      result.summary = "Main body extraction failed: " + body.error;
      publishDebugStep(debug, "H 主体提取失败", {{"输入点云", detectionInput}}, result.summary);
      return result;
    }
    result.mainBodyCloud = body.cloud;
    result.filteredCloud = body.cloud;
    result.removedBodyPoints = detectionInput->size() - body.cloud->size();
    sourceInput = body.cloud;
    if (frameBased) {
      std::size_t retained = 0;
      for (const auto &frame : *detectionFrames) {
        const auto begin = retained;
        while (retained < body.sourceIndices.size() &&
               body.sourceIndices[retained] < frame.begin + frame.count)
          ++retained;
        if (retained > begin)
          bodyFrames.push_back({begin, retained - begin, frame.frameIndex});
      }
    }
    if (stepDebugEnabled(debug)) {
      CloudT::Ptr rejected(new CloudT);
      std::size_t retained = 0;
      for (std::size_t i = 0; i < detectionInput->size(); ++i) {
        if (retained < body.sourceIndices.size() && body.sourceIndices[retained] == i)
          ++retained;
        else if (pcl::isFinite((*detectionInput)[i])) rejected->push_back((*detectionInput)[i]);
      }
      publishDebugStep(debug, "H 第一步 · 船体主要轮廓",
          {{"剔除点", rejected}, {"船体主体原始点", body.cloud}},
          "最大连通主体 + 有界曲面保护；保留 " + std::to_string(body.cloud->size()) +
          " 点，剔除 " + std::to_string(result.removedBodyPoints) +
          " 点。未改动原始坐标；底面和交界处仍需对照复核。",
          {{"retained_points", double(body.cloud->size())},
           {"removed_points", double(result.removedBodyPoints)},
           {"protected_uncertain_points", double(body.protectedUncertainIndices.size())}});
    }
    workingConfig.voxelSize = 0;
    workingConfig.outlierMeanK = 0;
  }
  CloudT::ConstPtr workingInput = sourceInput;
  finishStage("main_body");
  std::vector<FrameProjection> frameProjections;
  std::vector<ScanFrameRange> processedFrames;
  CloudT::Ptr measuredFrames;
  if (frameBased) {
    CloudT::Ptr flattened;
    float frameStep = 1.0f;
    if (!prepareFrameSections(sourceInput, config.extractMainBody ? bodyFrames : *detectionFrames,
                              workingConfig, flattened, measuredFrames,
                              processedFrames, frameProjections, frameStep,
                              result.summary))
      return result;
    // 帧号只替代截面生成。内部用原始帧间距估计连续性尺度，并把每帧
    // 暂时展平为一个剖面；候选、配对、路径和局部面修正仍走同一套流程。
    workingConfig.sliceStep = frameStep;
    workingConfig.sliceThickness = std::max(frameStep, 7.0f * frameStep);
    workingInput = flattened;
    result.filteredCloud = flattened;
    if (debug) {
      debug->clouds.push_back({"01_original_frames", detectionInput});
      publishDebugStep(debug, "S01 保留原始帧点",
                       {{"原始帧点云", detectionInput}},
                       "帧号定义截面，预处理不会混合不同帧。");
      publishDebugStep(debug, "S02 截面来源：原始帧",
                       {{"逐帧预处理点云", flattened}});
      publishDebugStep(debug, config.extractMainBody ? "S03 主体保护（保留逐帧原始点）" : "S03 逐帧体素降采样",
                       {{"处理后点云", flattened}});
    }
  } else {
    result.filteredCloud = preprocess(sourceInput, workingConfig, debug);
    workingInput = result.filteredCloud;
  }
  if (!result.filteredCloud || result.filteredCloud->size() < 30) {
    result.summary = "Too few points remain after preprocessing.";
    return result;
  }
  if (stableProfileSampling(workingConfig))
    workingConfig.profileGridOrigin =
        findProfileGridOrigin(*workingInput, workingConfig.longitudinalAxis);

  // S04-S08 生成截面候选并保留主体连续段。
  finishStage("prepare_sections");
  const auto range =
      findStationRange(*workingInput, workingConfig.longitudinalAxis);
  auto slices = frameBased
                    ? collectFrameSlices(*workingInput, processedFrames, *measuredFrames, frameProjections,
                                         workingConfig, leftHistory,
                                         rightHistory, debug)
                    : collectSlices(*workingInput, range, workingConfig,
                                    leftHistory, rightHistory, debug);
  const auto outputSliceCount =
      retainMainSegment(slices, workingConfig, leftHistory && rightHistory,
                        debug);
  finishStage("candidates");
  if (slices.size() < 3) {
    result.summary = "Fewer than three usable cross-sections were found.";
    return result;
  }

  // S09 成对优化左右轨迹，随后移除只参与优化的历史前视站位。
  std::vector<Candidate> leftPath, rightPath;
  if (!optimizePairedPaths(slices, workingConfig, debug, leftPath, rightPath,
                           result.summary))
    return result;
  // 前视站位已经完成了对末端路径的约束，现在立即裁掉，保证保存和显示中
  // 不会出现任何超出当前可靠几何范围的历史外推点。
  if (outputSliceCount > 0 && outputSliceCount < slices.size()) {
    slices.resize(outputSliceCount);
    leftPath.resize(outputSliceCount);
    rightPath.resize(outputSliceCount);
  }

  if (debug) {
    for (const auto &c : leftPath)
      debug->selectedPath[0].push_back(c.point);
    for (const auto &c : rightPath)
      debug->selectedPath[1].push_back(c.point);
  }

  publishPathStep(debug, "S09 成对路径优化", leftPath, rightPath);

  // S10 插补内部短缺口；S11-S16 平滑、贴壁及局部底面求交。
  std::vector<int> outputFrameIds;
  interpolateInternalGaps(*workingInput, range, slices, workingConfig, leftPath,
                          rightPath, debug,
                          frameBased ? &processedFrames : nullptr,
                          frameBased ? &outputFrameIds : nullptr);
  bridgeTrackedWalls(leftPath,workingConfig);
  bridgeTrackedWalls(rightPath,workingConfig);
  finishStage("path_and_interpolation");
  fillSideResult(leftPath, 0, workingInput, workingConfig, debug,
                  result.leftSeam,
                  &result.leftPointConfidence, &result.leftPointQuality,
                  &result.leftConfidence,
                  &result.leftHistoryUsageRatio);
  finishStage("left_refinement");
  fillSideResult(rightPath, 1, workingInput, workingConfig, debug,
                  result.rightSeam,
                  &result.rightPointConfidence, &result.rightPointQuality,
                  &result.rightConfidence,
                  &result.rightHistoryUsageRatio);
  finishStage("right_refinement");
  applySeamPostprocess(result, workingConfig);
  // A correct observed corner must not be pulled onto another feature by a
  // few wrong neighbouring candidates. Reject the entire excessive correction
  // instead of placing a shortened correction between two physical surfaces.
  std::array<std::vector<bool>, 2> correctionRejected{
      {std::vector<bool>(leftPath.size(), false), std::vector<bool>(rightPath.size(), false)}};
  auto restoreCandidate = [&](int side, std::size_t i) {
    auto &cloud = side ? result.rightSeam : result.leftSeam;
    auto &quality = side ? result.rightPointQuality : result.leftPointQuality;
    const auto &source = side ? rightPath : leftPath;
    (*cloud)[i] = source[i].point;
    quality[i].flags |= QualityCorrectionLimited | QualityWallUnavailable | QualityFloorUnavailable;
    quality[i].wallAccepted = quality[i].floorAccepted = false;
    quality[i].supportDistance = -1;
    quality[i].geometryCorrection = 0;
    if (!correctionRejected[side][i]) ++result.candidateCorrectionRejected;
    correctionRejected[side][i] = true;
  };
  for (int side = 0; side < 2; ++side) {
    auto &cloud = side ? result.rightSeam : result.leftSeam;
    const auto &source = side ? rightPath : leftPath;
    for (std::size_t i = 0; i < cloud->size(); ++i) {
      if ((toVector((*cloud)[i]) - toVector(source[i].point)).norm() <= config.maxCandidateCorrection)
        continue;
      restoreCandidate(side, i);
    }
  }
  const int widthAxis = transverseAxes(config.longitudinalAxis)[0];
  for (std::size_t i = 0; i < result.leftSeam->size(); ++i)
    if ((*result.rightSeam)[i].data[widthAxis] - (*result.leftSeam)[i].data[widthAxis] < workingConfig.minSeamSeparation) {
      restoreCandidate(0, i);
      restoreCandidate(1, i);
    }
  if (result.leftSmoothing.applied || result.rightSmoothing.applied) {
    result.leftConfidence =
        result.leftPointConfidence.empty()
            ? 0.0f
            : std::accumulate(result.leftPointConfidence.begin(),
                              result.leftPointConfidence.end(), 0.0f) /
                  static_cast<float>(result.leftPointConfidence.size());
    result.rightConfidence =
        result.rightPointConfidence.empty()
            ? 0.0f
            : std::accumulate(result.rightPointConfidence.begin(),
                              result.rightPointConfidence.end(), 0.0f) /
                  static_cast<float>(result.rightPointConfidence.size());
  }
  finishStage("candidate_guard");
  // Last transverse-coordinate operation: neither noisy local surfaces nor
  // per-point candidate rollback may undo this continuous final trajectory.
  applyFinalCurveSmoothing(result, workingConfig);
  if (config.finalSmoothingEnabled)
    publishDebugStep(debug, "S16b 最终轨迹连续平滑",
                     {{"左焊缝", result.leftSeam}, {"右焊缝", result.rightSeam}});
  finishStage("final_curve");
  std::vector<float> endpointStations;
  for (const auto &p:*result.leftSeam) endpointStations.push_back(coordinate(p,config.longitudinalAxis));
  const auto widthDecision=classifySeamWidth(result,workingConfig);
  if (frameBased) {
    restoreOriginalFrameGeometry(result, outputFrameIds, frameProjections,
                                 config.longitudinalAxis);
    result.filteredCloud.reset(new CloudT(*sourceInput));
  }
  // Emerging candidates remain eligible for output, but small consecutive
  // steps cannot prove that their trajectory follows the correct feature.
  // Apply the same sustained angle/deviation checks as for established walls;
  // the trimmer already preserves short tracks, gentle bends and local bumps.
  trimSeamEndpoints(result,workingConfig,&endpointStations,&widthDecision);
  const auto keptBegin=result.endpointTrimmedStart;
  const auto keptEnd=leftPath.size()-result.endpointTrimmedEnd;
  if (keptBegin || result.endpointTrimmedEnd) {
    leftPath=std::vector<Candidate>(leftPath.begin()+keptBegin,leftPath.begin()+keptEnd);
    rightPath=std::vector<Candidate>(rightPath.begin()+keptBegin,rightPath.begin()+keptEnd);
    publishDebugStep(debug,"S16c 首尾异常截断",{{"左焊缝",result.leftSeam},{"右焊缝",result.rightSeam}},
        "首端移除 "+std::to_string(keptBegin)+" 帧/站位，尾端移除 "+std::to_string(result.endpointTrimmedEnd)+
        (widthDecision.mode==1 ? "。依据单道中线截断，两侧范围同步。" : "。依据双道边界截断，两侧同步停止。"));
  }
  finishStage("endpoint_trim");
  // Preserve physical side/floor boundaries before output tracks are moved
  // by constant-width processing or replaced by a single center trajectory.
  std::vector<SurfaceSection> surfaceSections;
  for (std::size_t i=0; i<leftPath.size(); ++i) {
    SurfaceSection section;
    section.station=coordinate(leftPath[i].point,config.longitudinalAxis);
    section.left=(*result.leftSeam)[i];
    section.right=(*result.rightSeam)[i];
    if(frameBased) section.frameIndex=result.leftFrameIndices[i];
    auto normal=[](const Candidate &candidate, const SeamPointQuality &quality) {
      if(candidate.trackedWall) return candidate.trackedNormal;
      if(candidate.hasWall) return candidate.wallNormal;
      if(quality.wallAccepted) return quality.wallNormal;
      return Eigen::Vector2f(Eigen::Vector2f::Zero());
    };
    section.leftNormal=normal(leftPath[i],result.leftPointQuality[i]);
    section.rightNormal=normal(rightPath[i],result.rightPointQuality[i]);
    surfaceSections.push_back(section);
  }
  extractSurfaceClouds(sourceInput,surfaceSections,workingConfig,result,
      frameBased ? (config.extractMainBody ? &bodyFrames : detectionFrames) : nullptr,
      generateSurfaceModels);
  finishStage("surface_extraction");
  if(config.uniformWidthEnabled) {
    std::vector<UniformWidthSection> sections(leftPath.size());
    for(std::size_t i=0;i<leftPath.size();++i) {
      const float station=coordinate(leftPath[i].point,config.longitudinalAxis);
      auto slice=std::lower_bound(slices.begin(),slices.end(),station,
          [](const Slice &s,float x){return s.station<x;});
      if(slice==slices.end() || std::fabs(slice->station-station)>1.e-4f)continue;
      sections[i]=slice->widthSection;
      sections[i].sourceLeft=leftPath[i].point;sections[i].sourceRight=rightPath[i].point;
      if(frameBased) {
        auto plane=std::lower_bound(frameProjections.begin(),frameProjections.end(),slice->frameIndex,
            [](const FrameProjection &p,int id){return p.frameIndex<id;});
        if(plane==frameProjections.end() || plane->frameIndex!=slice->frameIndex) {sections[i].contour.clear();continue;}
        for(auto &p:sections[i].contour)p=restoreFramePoint(p,*plane,config.longitudinalAxis);
        sections[i].sourceLeft=restoreFramePoint(sections[i].sourceLeft,*plane,config.longitudinalAxis);
        sections[i].sourceRight=restoreFramePoint(sections[i].sourceRight,*plane,config.longitudinalAxis);
      }
    }
    applyUniformSeamWidth(result,workingConfig,sections);
    publishDebugStep(debug,"S16d 槽底统一宽度",{{"左焊缝",result.leftSeam},{"右焊缝",result.rightSeam}},
        "有效宽度 "+std::to_string(result.uniformWidthSamples)+" 个，两端各去掉 "+
        std::to_string(result.uniformWidthTrimEach)+" 个；目标宽度 "+std::to_string(result.uniformWidthTarget)+
        (config.uniformWidthMetric==UniformWidthMetric::BottomArc ? "（槽底弧长）" : "（两点直线距离）")+
        "。轮廓不足，未输出 "+std::to_string(result.uniformWidthOmitted)+" 帧/站位。");
    finishStage("uniform_width");
  }
  // S17 汇总左右结果。
  summarizeDetection(result, leftPath, rightPath);
  if (result.leftFinalCurve.applied || result.rightFinalCurve.applied) {
    result.summary += " Final curve smoothing: length=" + std::to_string(config.finalSmoothingLength) +
        ", mean correction L/R=" + std::to_string(result.leftFinalCurve.meanCorrection) + "/" +
        std::to_string(result.rightFinalCurve.meanCorrection) + ", max correction L/R=" +
        std::to_string(result.leftFinalCurve.maxCorrection) + "/" +
        std::to_string(result.rightFinalCurve.maxCorrection) + ".";
  }
  result.summary += config.useSurfaceRefinement
      ? (config.useLocalFloorConstraint ? " Surface refinement: wall + floor."
                                        : " Surface refinement: wall only.")
      : " Surface refinement: disabled (candidate correction guard still active).";
  if (config.extractMainBody)
    result.summary += " Main body: " + std::to_string(result.mainBodyCloud->size()) +
                      " original points retained, " + std::to_string(result.removedBodyPoints) + " removed.";
  if(config.inputVoxelSize>0)
    result.summary += " Input sampling: voxel="+std::to_string(config.inputVoxelSize)+
        ", retained="+std::to_string(result.inputSamplingPoints)+
        ", removed="+std::to_string(result.inputSamplingRemovedPoints)+" measured points.";
  result.summary += " Excessive candidate corrections rejected: " + std::to_string(result.candidateCorrectionRejected) + ".";
  if (frameBased)
    result.summary += " Original-frame input: " +
                      std::to_string(result.inputFrameCount) +
                      " frames; the common detection pipeline was used without "
                      "spatial slicing.";
  if (config.trimUnstableEndpoints) result.summary+=" Endpoint trim: head="+std::to_string(result.endpointTrimmedStart)+
      ", tail="+std::to_string(result.endpointTrimmedEnd)+" stations removed.";
  if(config.uniformWidthEnabled)
    result.summary+=" Uniform bottom width: target="+std::to_string(result.uniformWidthTarget)+
        (config.uniformWidthMetric==UniformWidthMetric::BottomArc ? " (arc)" : " (chord)")+
        ", samples="+std::to_string(result.uniformWidthSamples)+", trimmed_each="+std::to_string(result.uniformWidthTrimEach)+
        ", adjusted="+std::to_string(result.uniformWidthAdjusted)+", omitted="+std::to_string(result.uniformWidthOmitted)+".";
  if(result.uniformWidthFilterApplied)
    result.summary+=" Uniform-width final filter: length="+std::to_string(config.uniformWidthFilterLength)+
        ", mean correction="+std::to_string(result.uniformWidthFilterMeanCorrection)+
        ", max correction="+std::to_string(result.uniformWidthFilterMaxCorrection)+".";
  if(!config.uniformWidthEnabled)applySeamWidthRule(result,config,&widthDecision);
  result.summary+=" Surface points L/R/weld="+std::to_string(result.leftSurface->size())+"/"+
      std::to_string(result.rightSurface->size())+"/"+std::to_string(result.weldSurface->size())+".";
  result.summary+=" B-spline patches L/weld/R="+std::to_string(result.leftSurfaceModel.patches.size())+"/"+
      std::to_string(result.weldSurfaceModel.patches.size())+"/"+std::to_string(result.rightSurfaceModel.patches.size())+".";
  if(!generateSurfaceModels) result.summary+=" B-spline fitting disabled.";
  publishDebugStep(debug, "S17 最终焊缝",
                   {{"处理后点云", result.filteredCloud},
                    {"左焊缝", result.leftSeam},
                    {"右焊缝", result.rightSeam}, {"中间焊缝", result.centerSeam},
                    {"左侧面", result.leftSurface}, {"右侧面", result.rightSurface},
                    {"焊道面", result.weldSurface}},
                   result.summary);
  finishStage("restore_and_output");
  return result;
}

} // namespace weld
