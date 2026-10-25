#pragma once

#include "SeamTypes.h"
#include <array>
#include <functional>
#include <map>

namespace weld {

// 仅在 detect(..., &trace) 时收集中间结果；不改变候选与轨迹的计算。
// 二维数据均使用固定剖面轴 (u,v)，三维数据保持原始 XYZ。
struct DebugProfilePoint {
  float u, v, support;
  std::array<float, 2> mean;
  std::array<float, 3> scatter; // xx, xy, yy：箱内离差和（未除以点数）
};

struct DebugLine {
  float u = 0, v = 0, du = 0, dv = 0, rms = 0;
  int begin = 0, end = 0; // 拟合窗口，半开区间
};

struct DebugPair {
  float bottomMidpointU = 0;
  bool hasBottomMidpoint = false;
  std::array<unsigned, 2> candidateKinds{{0, 0}};
  std::array<PointT, 2> points;
  std::array<float, 2> feature{{0, 0}}, nodeCost{{0, 0}};
  std::array<float, 2> historyDistance{{-1, -1}};
  std::string source = "intersection";
  std::string profileShape = "unknown";
  float proposalCost = 0;
  std::array<DebugLine, 3> lines; // 左侧、跨槽底面、右侧
  std::array<int, 2> boundaries{{-1, -1}};
  std::array<int, 6> regions{{0, 0, 0, 0, 0, 0}}; // 左、底、右 begin/end
  std::array<int, 6> seeds{{0, 0, 0, 0, 0, 0}};
  std::array<float, 8> edges{{0, 0, 0, 0, 0, 0, 0, 0}};
  std::array<float, 2> edgeGaps{{-1, -1}}, angles{{0, 0}}, localAngles{{0, 0}};
  std::array<float, 2> clearance{{0, 0}};
  bool regionAttempted = false;
  std::array<bool, 2> regionUsed{{false, false}};
  std::array<int, 2> localCornerIndices{{-1, -1}};
  std::array<float, 2> localCornerAngles{{0, 0}};
  std::array<float, 2> profileCurvatures{{0, 0}};
  std::array<DebugLine, 2> candidateWalls; // Observed one-sided local wall tangents.
};

struct DebugTurnRegion {
  PointT begin, end, peak, representative;
  float arcLength = 0, turningAngle = 0;
};

struct DebugSlice {
  float featureDepth = -1.f; // Robust measured plate-normal relief; -1 unavailable.
  int frameIndex = -1;
  float station = 0;
  std::vector<PointT> raw;
  // Original measured geometry for locating a frame inside the complete 3-D
  // cloud. raw may be projected into the detector's unfolded coordinates.
  std::vector<PointT> sourceRaw;
  std::vector<DebugProfilePoint> profile;
  std::vector<PointT> localizationProfile, fittedProfile;
  std::vector<float> fittedCurvature;
  std::vector<DebugTurnRegion> turnRegions;
  std::array<HistoryPrediction, 2> predictions;
  std::vector<DebugPair> pairs;
  // Candidates restored to the original scan plane for overlaying on the
  // complete cloud. Spatial slices already use original coordinates.
  std::array<std::vector<PointT>, 2> sourceCandidates;
  // 每个枚举模型只计入最先失败的门槛；accepted 是通过全部门槛的数量。
  std::map<std::string, std::size_t> gates;
  std::string status = "raw_points_insufficient";
  // Human-readable explanation when this slice produced no left/right pair.
  // Stored with the slice so frame replay shows the same reason as live debug.
  std::string noCandidateReason;
  bool hasGeometry = false, retained = false, lookahead = false;
  int selected = -1;
  float selectedNodeCost = 0, transitionCost = 0, tangentCost = 0;
};

struct DebugCloudStage {
  std::string name;
  CloudT::ConstPtr cloud;
};

// 通用算法步骤，不依赖参考算法的截面、PCA 或路径规划实现。
struct DebugStep {
  bool hasFrameView = false;
  Eigen::Vector3f viewOrigin = Eigen::Vector3f::Zero();
  Eigen::Vector3f viewU = Eigen::Vector3f::UnitY(), viewV = Eigen::Vector3f::UnitZ();
  std::string title, status, message;
  double elapsedMs = 0;
  std::map<std::string, double> metrics;
  std::vector<DebugCloudStage> layers; // 发布时深拷贝，后续修改不会污染已记录步骤。
  // 帧/截面调试时，左侧总览显示完整点云并高亮当前切片；右侧仍显示
  // layers 中的截面轮廓与候选细节。
  CloudT::ConstPtr overviewCloud, positionCloud;
  std::array<CloudT::ConstPtr, 2> overviewCandidates;
};

struct DebugWallPoint {
  int side=0;
  PointT before, after;
  PointT beforeFloor;
  Eigen::Vector2f normal=Eigen::Vector2f::Zero();
  float offset=0, rms=0, distanceBefore=0, distanceAfter=0, supportDistance=-1;
  std::string status;
  std::vector<PointT> wallPatch, wallInliers;
  std::string wallFitMethod;
  float wallRmsBeforeRefit=-1;
  std::vector<PointT> floorPatch, floorCorners;
  Eigen::Vector3f floorNormal=Eigen::Vector3f::Zero();
  float floorOffset=0, floorRms=-1, floorResidual=-1, floorShift=0;
  int floorInliers=0;
  std::string floorStatus="not_attempted";
};

/** @brief 可选的详细调试输出；detect 开始时重置内容，但保留 onStep。
 * @details 二维点使用剖面 (u,v)，sourceRaw/sourceCandidates 为原始 XYZ。
 * 数据量随点云和候选数量增长；只需最终结果时传 nullptr。 */
struct DebugTrace {
  // 返回 false 可跳过本次检测剩余弹窗；算法仍继续执行。
  /** @brief 同步步骤回调，输入为当前步骤的只读引用。
   * @return false 停止后续步骤回调，但算法继续；引用仅在回调期间有效。
   * @note 若需跨线程交给 UI，调用者负责复制数据及线程切换。 */
  std::function<bool(const DebugStep &)> onStep;
  std::string algorithm = "reference";
  std::vector<DebugStep> steps;
  DetectorConfig config;
  std::size_t leftHistoryPasses = 0, rightHistoryPasses = 0;
  std::vector<DebugCloudStage> clouds;
  std::vector<DebugSlice> slices;
  std::array<std::vector<PointT>, 2> selectedPath, interpolatedPath;
  std::array<std::vector<PointT>, 2> beforeWallPath;
  std::vector<DebugWallPoint> wallPoints;
  float optimalCost = 0;
  bool pathComplete = false;
};

} // namespace weld
