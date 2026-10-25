#pragma once

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include "BSplineSurface.h"

#include <Eigen/Core>

#include <cstddef>
#include <cstdint>
#include <string>
#include <map>
#include <vector>

namespace weld {

using PointT = pcl::PointXYZ;
using CloudT = pcl::PointCloud<PointT>;

enum class Axis { X = 0, Y = 1, Z = 2 };

// LegacyRolling preserves the original data-anchored rolling bins and all
// thresholds coupled to profileBinSize. StableGrid uses one grid origin for
// the whole cloud and keeps the downstream physical scale independent from
// the sampling bin width.
enum class ProfileSamplingMode { LegacyRolling = 0, StableGrid = 1 };
enum class ProfileLocalizationMethod {
  CurveLineIntersection = 0,
  CurvatureCandidates = 1
};
enum class ProfileTurnLocation {
  Peak = 0,
  RegionMedian = 1,
  TurnMedian = 2,
  FloorSideBoundary = 3,
  WallTurnCentroid = 4
};

// Optional RTS stage before the candidate guard and final continuous fit.
// Final-curve smoothing has its own independent switch.
enum class SeamPostprocessMode { Legacy = 0, RobustRts = 1 };
enum class UniformWidthMetric { BottomArc = 0, Chord = 1 };

enum SeamQualityFlag : std::uint32_t {
  QualityNone = 0,
  QualityInterpolated = 1u << 0,
  QualityWallUnavailable = 1u << 1,
  QualityFloorUnavailable = 1u << 2,
  QualityInnovationOutlier = 1u << 3,
  QualityCorrectionLimited = 1u << 4,
  QualityGeometryLimited = 1u << 5,
  QualityPairConstraintLimited = 1u << 6,
  QualityFinalCurveSmoothed = 1u << 7,
  QualityCurvatureFallback = 1u << 8,
  QualityUniformWidth = 1u << 9
};

// 原始输入中的连续帧范围；保留真实帧号，不根据坐标或点距猜帧。
struct ScanFrameRange {
  /// @brief input 中起始下标与连续点数；范围 [begin, begin+count)。
  std::size_t begin = 0, count = 0;
  /// @brief 原始非负帧号，必须严格递增，无需从 0 开始。
  int frameIndex = 0;
};

// Stored separately for the standalone curvature algorithm and the curve/line
// fallback. These are values, never references into the other algorithm's settings.
struct CurvatureSearchConfig {
  /// @brief 叠加波谷/曲率变化率候选，默认 true；壁脚转角重心模式不使用此扩展。
  bool multiCandidates = true;
  /// @brief 候选搜索弧长窗口，长度单位；有效范围 [1e-6,1e9]。
  float windowLength = 2.0f;
  /// @brief 波谷两侧最小突出深度，长度单位；有效范围 [1e-6,1e9]。
  float valleyProminence = 0.08f;
  /// @brief 候选去重间距，长度单位；有效范围 [1e-6,1e9]。
  float candidateSpacing = 0.4f;
  /// @brief 每帧最多保留的左右候选对数；整数 [2,32]。
  int maximumPairs = 8;
  /// @brief 二阶路径连续性权重；有效范围 [0,1e9]。
  float continuityWeight = 2.0f;
};

struct CurveLineFallbackConfig {
  /// @brief 曲线求交整帧无候选对时，是否启用独立曲率回退。
  bool enabled = true;
  /// @brief 回退样条拟合半径，长度单位；[1e-6,1e9]。
  float splineRadius = 8.0f;
  /// @brief 回退候选到侧壁的最大距离，长度单位；[1e-6,1e9]。
  float wallTolerance = 0.6f;
  /// @brief 回退转折区域最大弧长，长度单位；[1e-6,1e9]。
  float regionLength = 4.0f;
  /// @brief 回退区域相对峰值比例；[1e-6,1]。
  float regionFraction = 0.6f;
  /// @brief 回退的独立候选搜索参数，不引用主曲率算法参数。
  CurvatureSearchConfig search;
};

/**
 * @brief 单次检测的完整参数快照。
 * @details 所有长度使用点云单位，默认按毫米设计；角度为度；权重无量纲。
 * 字段初值是曲线求交默认值；曲率默认值通过 AlgorithmParameters::get 读取。
 * JSON 配置中的字段名与此结构相同，枚举用整数；内部 profileGridOrigin 除外。
 */
struct DetectorConfig {
  /// @brief 连续曲率算法的独立搜索配置。
  CurvatureSearchConfig curvatureSearch;
  /// @brief 曲线求交算法的独立回退配置。
  CurveLineFallbackConfig curveLineFallback;
  // First isolate the measured hull body, using the reference project's
  // bounded surface recovery. No generated/smoothed points enter this mask.
  /** @brief 先提取工件主体。 @details 默认 true；有效范围 true/false。 */
  bool extractMainBody = true;
  /** @brief 主体提取体素边长（长度）。 @details 默认 1.0f；有效范围 0.000001..1000000000。 */
  float hullSurfaceVoxelSize = 1.0f;
  /** @brief 主体表面曲率阈值（无量纲）。 @details 默认 0.035f；有效范围 0..1。 */
  float hullSurfaceCurvature = 0.035f;
  /** @brief 主体表面距离阈值（长度）。 @details 默认 0.6f；有效范围 0.000001..1000000000。 */
  float hullSurfaceDistance = 0.6f;
  /** @brief 主体保留边界余量（长度）。 @details 默认 3.0f；有效范围 0..1000000000。 */
  float hullPreservationMargin = 3.0f;
  // 0: each algorithm stage; >0: pause after this many frames/sections.
  /** @brief 调试步长（帧）；0 按算法阶段回调。 @details 默认 1；有效范围 0..100000。 */
  int debugFrameBatchSize = 1;
  // 1-based frame/section ordinal at which frame stepping starts. Earlier
  // sections are still computed so path selection and the final result stay identical.
  /** @brief 调试起始帧/截面序号（从 1 开始）。 @details 默认 1；有效范围 1..100000。 */
  int debugStartSlice = 1;
  // GUI frame selection keeps the real source frame id so main-body cropping
  // cannot shift the requested slice to a different original frame.
  /** @brief 选定的原始帧号；-1 不指定。 @details 默认 -1；有效范围 -1..2147483647。 */
  int debugSelectedFrameIndex = -1;
  // +1: a trough opening toward +V; -1: the reverse; 0: legacy either sign.
  /** @brief 槽口朝向：-1 负 V、0 任意、1 正 V。 @details 默认 1；有效范围 -1..1。 */
  int interiorOpeningDirection = 1;
  // Do not let smoothing/local plane fits move a measured corner arbitrarily.
  /** @brief 候选允许的最大位移（长度）。 @details 默认 1.0f；有效范围 0..1000000000。 */
  float maxCandidateCorrection = 1.0f;
  // 开启后由 GUI/CLI 的 VTK 调试查看器逐步显示中间点云。
  /** @brief 启用调试步骤回调；库本身不创建窗口。 @details 默认 false；有效范围 true/false。 */
  bool debugFlag = false;
  // 方案1：在候选内侧取局部底面正方形，以底面、侧壁、当前站位共同定位。
  /** @brief 局部底面正方形边长（长度）。 @details 默认 6.0f；有效范围 0.000001..1000000000。 */
  float floorPatchSize = 6.0f;
  // 总开关：关闭时跳过 S12-S16 侧壁拟合、贴壁投影与底面求交。
  // 路径平滑、候选修正上限及左右间距保护仍然生效。
  /** @brief 启用侧壁/底面几何细化。 @details 默认 false；有效范围 true/false。 */
  bool useSurfaceRefinement = false;
  /** @brief 启用局部底面约束；受几何细化总开关控制。 @details 默认 true；有效范围 true/false。 */
  bool useLocalFloorConstraint = true;
  /** @brief 底面采样避让距离（长度）。 @details 默认 0.6f；有效范围 0..1000000000。 */
  float floorPatchGuard = 0.6f;
  /** @brief 底面拟合最大 RMS（长度）。 @details 默认 0.4f；有效范围 0.000001..1000000000。 */
  float floorPatchMaxRms = 0.4f;
  /** @brief 底面约束最大修正（长度）。 @details 默认 5.0f；有效范围 0..1000000000。 */
  float floorPatchMaxCorrection = 5.0f;
  /** @brief 纵向轴：0=X、1=Y、2=Z。 @details 默认 Axis::X；有效范围 0..2。 */
  Axis longitudinalAxis = Axis::X;

  // 所有长度参数都使用点云自身单位；当前样例数据按毫米设置默认值。
  /** @brief 预处理体素边长（长度）；0 禁用体素降采样。 @details 默认 0.8f；有效范围 0..1000000000。 */
  float voxelSize = 0.8f;
  /** @brief 统计离群滤波邻居数（点）；0 禁用。 @details 默认 12；有效范围 0..100000。 */
  int outlierMeanK = 12;
  /** @brief 离群标准差倍数。 @details 默认 2.0f；有效范围 0..1000000000。 */
  float outlierStddev = 2.0f;
  /** @brief 无原始帧时的空间切片步长（长度）。 @details 默认 2.0f；有效范围 0.000001..1000000000。 */
  float sliceStep = 2.0f;
  /** @brief 无原始帧时的切片厚度（长度）。 @details 默认 3.5f；有效范围 0.000001..1000000000。 */
  float sliceThickness = 3.5f;
  // Project a centered window of consecutive original frames onto the
  // current scan plane before profile binning and B-spline fitting.
  /** @brief 启用相邻原始帧投影合并。 @details 默认 true；有效范围 true/false。 */
  bool mergeProfileFrames = true;
  /** @brief 原始帧合并窗口（帧）；1~101 的奇数。 @details 默认 5；有效范围 1..101。 */
  int profileMergeFrameCount = 5;
  /** @brief 轮廓横向分箱尺寸（长度）。 @details 默认 1.2f；有效范围 0.000001..1000000000。 */
  float profileBinSize = 1.2f;
  /** @brief 分箱方式：0 滚动分箱、1 固定网格。 @details 默认 ProfileSamplingMode::LegacyRolling；有效范围 0..1。 */
  ProfileSamplingMode profileSamplingMode =
      ProfileSamplingMode::LegacyRolling;
  // In StableGrid mode this is the physical reference used by fitting
  // windows, gap checks and wall-support gates. Changing profileBinSize then
  // changes sampling resolution without silently changing those tolerances.
  /** @brief 固定网格模式的物理拟合尺度（长度）。 @details 默认 1.2f；有效范围 0.000001..1000000000。 */
  float profileReferenceScale = 1.2f;
  // Internal per-run anchor calculated from the complete preprocessed cloud.
  // It is intentionally not exposed as a user parameter.
  /** @brief 内部网格原点，由算法计算，可为负；不从配置文件读取，调用者无需设置。 */
  float profileGridOrigin = 0.0f;
  /** @brief 结构拟合半窗口（采样点）。 @details 默认 5；有效范围 1..100000。 */
  int fitHalfWindow = 5;
  // 拟合三块结构面时跳过交接点附近的焊道冠部，避免焊道余高拉偏板面。
  /** @brief 交界处跳过的采样点数。 @details 默认 1；有效范围 0..100000。 */
  int fitGuardPoints = 1;
  // 底面相对剖面横向的最大斜率，以及侧面与底面的最小夹角。
  // 二者共同排除底面纹理、扫描边缘等不符合“侧面-底面”结构的折角。
  /** @brief 底面允许的最大绝对斜率。 @details 默认 0.9f；有效范围 0..1000000000。 */
  float maxBottomSlope = 0.9f;
  /** @brief 侧壁与底面最小结构夹角（度）。 @details 默认 20.0f；有效范围 0..180。 */
  float minSurfaceAngleDegrees = 20.0f;
  // 大窗口确认三块大面的总体方向，小窗口再确认交界处确实发生了局部法向
  // 突变。后者用于排除把一整块平滑弧面误拆成三条直线的情况。
  /** @brief 局部表面最小夹角（度）。 @details 默认 8.0f；有效范围 0..180。 */
  float minLocalSurfaceAngleDegrees = 8.0f;
  /** @brief 局部法向估计窗口（采样点）。 @details 默认 3；有效范围 1..100000。 */
  int localAngleWindow = 3;
  // Fit each connected profile with one parametric cubic B-spline. The
  // bottom remains a coarse structure cue, not a required local flat patch.
  /** @brief 启用 B 样条轮廓拟合。 @details 默认 true；有效范围 true/false。 */
  bool useProfileCurvature = true;
  // Default path: fit each sloped wall as a parametric curve, fit the groove
  // floor as a line, and use their two intersections as the candidates.
  /** @brief 定位算法：0 曲线求交、1 连续曲率；须与所属参数组一致。 @details 默认 ProfileLocalizationMethod::CurveLineIntersection；有效范围 0..1。 */
  ProfileLocalizationMethod profileLocalizationMethod =
      ProfileLocalizationMethod::CurveLineIntersection;
  // A central raised bead is treated as a separate shape: the two continuous
  // valleys beside it become the candidates instead of the outer wall/floor
  // intersections. Units match the point cloud.
  /** @brief 中央凸起最小高度（长度）。 @details 默认 0.8f；有效范围 0.000001..1000000000。 */
  float centralBeadMinHeight = 0.8f;
  // Minimum measured relief normal to the plate reference; zero disables.
  // Independent per algorithm. A rejected section is a hard path boundary.
  /** @brief 有效槽深下限（长度）；0 禁用浅槽截断。 @details 默认 1.0f；有效范围 0..1000000000。 */
  float minimumFeatureDepth = 1.0f;
  /** @brief 中央凸起最小宽度（长度）。 @details 默认 2.0f；有效范围 0.000001..1000000000。 */
  float centralBeadMinWidth = 2.0f;
  // Skip broad surface hypotheses that cannot reach any supported wall peak.
  // Disable only for exhaustive-search equivalence/performance comparisons.
  /** @brief 使用侧壁支持提前剪枝候选。 @details 默认 true；有效范围 true/false。 */
  bool useWallCandidatePruning = true;
  /** @brief 轮廓样条拟合半径（长度）。 @details 默认 8.0f；有效范围 0.000001..1000000000。 */
  float profileCurvatureRadius = 8.0f;
  /** @brief 转折代表点：0 峰值、1 区域中值、2 转角中值、3 底面侧边、4 壁脚转角重心。 @details 默认 ProfileTurnLocation::RegionMedian；有效范围 0..4。 */
  ProfileTurnLocation profileTurnLocation = ProfileTurnLocation::RegionMedian;
  /** @brief 高曲率区域相对峰值比例。 @details 默认 0.6f；有效范围 0.000001..1。 */
  float profileTurnRegionFraction = 0.6f;
  // Maximum total arc length of one high-curvature region, in cloud units.
  /** @brief 单个转折区域最大弧长（长度）。 @details 默认 4.0f；有效范围 0.000001..1000000000。 */
  float profileTurnRegionMaxLength = 4.0f;
  // Hard perpendicular distance to the observed stable wall for a local turn.
  /** @brief 候选到独立侧壁的距离上限（长度）。 @details 默认 0.6f；有效范围 0.000001..1000000000。 */
  float profileWallTolerance = 0.6f;
  /** @brief 按整条轨迹槽宽投票选择单道或双道。 @details 默认 false；有效范围 true/false。 */
  bool adaptiveSeamCount = false;
  /** @brief 单/双道槽宽判定阈值（长度）。 @details 默认 13.0f；有效范围 0.000001..1000000000。 */
  float doubleSeamWidth = 13.0f;
  // 底面在每个交界处分别做多尺度局部拟合。窗口逐步扩大时，如果法向
  // 变化超过该角度，就停止扩大，避免把焊接变形或另一处折角混入切平面。
  /** @brief 局部拟合窗口增长时法向变化上限（度）。 @details 默认 6.0f；有效范围 0..180。 */
  float maxLocalFitNormalChangeDegrees = 6.0f;
  // 左右两个局部底面允许有缓慢翘曲，但不能指向完全不同的结构面。
  /** @brief 左右局部底面法向差上限（度）。 @details 默认 18.0f；有效范围 0..180。 */
  float maxLocalBottomAngleDifferenceDegrees = 18.0f;
  // 在每个二维剖面内从底面中央、左右侧面稳定段分别开始区域生长。
  // 相邻局部切线变化小且点到当前切线足够近时归入同一区域；这样平滑
  // 弯曲的底面不会被强制压成直线，同时会在侧面法向突变处停止。
  /** @brief 启用区域生长边界定位。 @details 默认 true；有效范围 true/false。 */
  bool useRegionGrowingEdges = true;
  // 底面生长跟踪分箱中值形成的轮廓切线。厚截面内的纵向高度变化会
  // 增大箱内竖向离差，不能把这种厚度当成底面在剖面内的法向变化。
  // 侧面和三面结构确认仍使用原始分布的加权 PCA。
  /** @brief 区域生长使用轮廓切线。 @details 默认 true；有效范围 true/false。 */
  bool useContourRegionGrowing = true;
  // 只有侧面与底面方向充分分离时才以局部轮廓折点校正位置；浅夹角
  // 焊道冠部可能产生多个焊趾，继续使用区域配对，避免追到冠部另一侧。
  /** @brief 轮廓折点校正的最小转角（度）。 @details 默认 50.0f；有效范围 0..180。 */
  float minContourCornerAngleDegrees = 50.0f;
  /** @brief 区域生长局部窗口（采样点）。 @details 默认 4；有效范围 3..100000。 */
  int regionGrowingWindow = 4;
  /** @brief 区域生长相邻切线变化上限（度）。 @details 默认 10.0f；有效范围 0.000001..180。 */
  float regionGrowingSmoothnessDegrees = 10.0f;
  /** @brief 区域生长相对种子法向漂移上限（度）。 @details 默认 35.0f；有效范围 0.000001..180。 */
  float regionGrowingMaxNormalDriftDegrees = 35.0f;
  // 基础点线距离。底面中值轮廓有噪声时按局部 RMS 自适应，最多放宽到
  // 两倍；步进角和相对种子的法向漂移仍作为独立限制。
  /** @brief 区域生长基础点线距离（长度）。 @details 默认 2.0f；有效范围 0.000001..1000000000。 */
  float regionGrowingDistance = 2.0f;
  // 底面边缘与对应侧面边缘之间允许的最大缺口。候选位置取两个边缘
  // 最近点的中点，不再主要依赖两条无限延长直线的理论交点。
  /** @brief 侧壁/底面边界最大间隙（长度）。 @details 默认 8.0f；有效范围 0.000001..1000000000。 */
  float regionGrowingMaxEdgeGap = 8.0f;
  // 拟合面理论交点与分段边界允许的最大距离。
  /** @brief 理论交点到分段边界最大距离（长度）。 @details 默认 10.0f；有效范围 0..1000000000。 */
  float surfaceIntersectionTolerance = 10.0f;
  // 侧面样本到当前截面局部上表面切线的最小正交距离。它衡量侧面是否
  // 真正从弧形上表面分离，不能用全局 Z 高度代替。
  /** @brief 侧面离上表面的最小距离（长度）。 @details 默认 2.0f；有效范围 0..1000000000。 */
  float minSurfaceClearance = 2.0f;
  // 三个面的正交拟合残差转成置信度时使用的毫米尺度。
  /** @brief 结构残差转置信度的尺度（长度）。 @details 默认 2.5f；有效范围 0.000001..1000000000。 */
  float surfaceRmsScale = 2.5f;
  /** @brief 每截面最大候选数量。 @details 默认 7；有效范围 1..100000。 */
  int maxCandidatesPerSlice = 7;
  /** @brief 同侧候选最小间隔（长度）。 @details 默认 4.0f；有效范围 0..1000000000。 */
  float minCandidateSpacing = 4.0f;
  // 两侧板面拟合线的理论交点只有落在此宽度内才可信，防止近似平行的
  // 拟合线在远处相交后把候选点拉离真实焊缝。
  /** @brief 可信理论交点的最大宽度（长度）。 @details 默认 18.0f；有效范围 0..1000000000。 */
  float maxBeadWidth = 18.0f;
  // 同一截面内左右焊缝在横向坐标上的最小间距。配对约束可防止把同一
  // 焊缝附近的两个焊趾误识别为左右两道焊缝。
  /** @brief 左右焊缝最小横向间距（长度）；曲率默认 6，求交默认 8。 @details 默认 8.0f；有效范围 0..1000000000。 */
  float minSeamSeparation = 8.0f;
  /** @brief 相邻站位路径转移距离上限（长度）。 @details 默认 12.0f；有效范围 0.000001..1000000000。 */
  float maxTransitionDistance = 12.0f;
  /** @brief 路径一阶平滑代价权重。 @details 默认 2.0f；有效范围 0..1000000000。 */
  float smoothnessWeight = 2.0f;
  // 相比只限制相邻点距离，切线预测残差可以允许连续陡坡，同时抑制候选
  // 在相邻截面间切换到另一块表面。
  /** @brief 路径切线预测残差权重。 @details 默认 4.0f；有效范围 0..1000000000。 */
  float curvatureWeight = 4.0f;
  /** @brief 几何特征代价权重。 @details 默认 4.0f；有效范围 0..1000000000。 */
  float featureWeight = 4.0f;

  // 每个横向分箱在厚截面内独立估计纵向高度斜率，再投影回当前站位。
  // 真实点云的扫描覆盖变化可能破坏该估计，因此默认关闭；该值大于 0 时
  // 表示参与斜率计算的最小纵向点间距。
  /** @brief 纵向去趋势的最小估计跨度（长度）；0 禁用。 @details 默认 0.0f；有效范围 0..1000000000。 */
  float longitudinalDetrendMinSpan = 0.0f;

  // 历史项是软约束。特征越弱，它在总分中的占比越高。
  /** @brief 历史轨迹代价权重。 @details 默认 2.5f；有效范围 0..1000000000。 */
  float historyWeight = 2.5f;
  /** @brief 弱特征时历史权重增益。 @details 默认 3.0f；有效范围 0..1000000000。 */
  float weakFeatureHistoryBoost = 3.0f;
  /** @brief 历史候选搜索半径（长度）。 @details 默认 18.0f；有效范围 0..1000000000。 */
  float historySearchRadius = 18.0f;
  /** @brief 历史层间外推位移上限（长度）。 @details 默认 8.0f；有效范围 0..1000000000。 */
  float maxHistoryExtrapolation = 8.0f;
  /** @brief 历史轨迹容量（条）；调用者创建 SeamHistory 时传入。 @details 默认 5；有效范围 1..100000。 */
  std::size_t maxHistoryPasses = 5;
  // 历史可在可靠几何段末端之后提供少量“只参与优化、不输出”的前视站位，
  // 用来稳定最后几个已被焊道填圆的交点。
  /** @brief 历史仅用于优化的前视站位数。 @details 默认 8；有效范围 0..100000。 */
  int historyLookaheadStations = 8;

  // 只插补少量内部漏检站位；首尾或长距离贴近上表面的区段保持不输出。
  /** @brief 允许插补的内部缺口站位数。 @details 默认 3；有效范围 0..100000。 */
  int maxInterpolatedGapStations = 3;
  // 输出曲线的二阶正则权重。它使用候选置信度作为数据权重，低置信度
  // 缺口主要服从前后切线趋势，不再使用会压平末端斜率的移动平均。
  /** @brief 轨迹二阶正则权重。 @details 默认 12.0f；有效范围 0..1000000000。 */
  float curveRegularizationWeight = 12.0f;

  // 最终轨迹在贴面/候选回退之后统一拟合，避免逐点回退重新制造尖刺。
  // 长度及异常残差尺度使用点云单位，不使用恢复后的非单调 X 作为帧序。
  /** @brief 启用最终连续轨迹平滑。 @details 默认 true；有效范围 true/false。 */
  bool finalSmoothingEnabled = true;
  /** @brief 启用槽底统一宽度双道输出。 @details 默认 false；有效范围 true/false。 */
  bool uniformWidthEnabled = false;
  /** @brief 统一宽度度量：0 槽底弧长、1 直线距离。 @details 默认 UniformWidthMetric::Chord；有效范围 0..1。 */
  UniformWidthMetric uniformWidthMetric = UniformWidthMetric::Chord;
  /** @brief 统一宽度后启用整体滤波。 @details 默认 true；有效范围 true/false。 */
  bool uniformWidthFilterEnabled = true;
  /** @brief 统一宽度滤波窗口（长度）。 @details 默认 12.0f；有效范围 0.000001..1000000000。 */
  float uniformWidthFilterLength = 12.0f;
  /** @brief 最终轨迹平滑窗口（长度）。 @details 默认 12.0f；有效范围 0.000001..1000000000。 */
  float finalSmoothingLength = 12.0f;
  /** @brief 最终平滑稳健残差尺度（长度）。 @details 默认 0.6f；有效范围 0.000001..1000000000。 */
  float finalSmoothingRobustScale = 0.6f;
  // Walk outward from the central reference span and stop at sustained bends/jumps.
  /** @brief 截去两端持续不稳定的轨迹。 @details 默认 true；有效范围 true/false。 */
  bool trimUnstableEndpoints = true;
  /** @brief 端点截断检查窗口（长度）。 @details 默认 12.0f；有效范围 0.000001..1000000000。 */
  float endpointTrimWindow = 12.0f;
  /** @brief 端点持续弯折阈值（度）。 @details 默认 20.0f；有效范围 0.000001..89.9999。 */
  float endpointTrimAngleDegrees = 20.0f;
  /** @brief 端点偏离阈值（长度）。 @details 默认 1.5f；有效范围 0.000001..1000000000。 */
  float endpointTrimDeviation = 1.5f;

  // Optional geometry-constrained Kalman/RTS pass, before the candidate guard
  // and the independently enabled final continuous fit.
  /** @brief 后处理：0 原有流程、1 稳健 Kalman/RTS。 @details 默认 SeamPostprocessMode::Legacy；有效范围 0..1。 */
  SeamPostprocessMode postprocessMode = SeamPostprocessMode::Legacy;
  // All values use the point-cloud coordinate unit.  The state is
  // [position, d(position)/d(station)] for each transverse coordinate.
  /** @brief RTS 过程噪声强度。 @details 默认 0.04f；有效范围 0..1000000000。 */
  float rtsProcessNoise = 0.04f;
  /** @brief RTS 观测标准差（长度）。 @details 默认 0.6f；有效范围 0.000001..1000000000。 */
  float rtsMeasurementStddev = 0.6f;
  /** @brief RTS 归一化创新平方门限。 @details 默认 9.21f；有效范围 0.000001..1000000000。 */
  float rtsInnovationGate = 9.21f;
  /** @brief RTS 最大修正（长度）。 @details 默认 2.0f；有效范围 0..1000000000。 */
  float rtsMaxCorrection = 2.0f;
  // A smoothed point is allowed to leave an accepted local wall/floor model
  // only by this residual.  Larger corrections are shortened, not forced.
  /** @brief RTS 允许偏离几何模型的距离（长度）。 @details 默认 0.6f；有效范围 0..1000000000。 */
  float rtsGeometryTolerance = 0.6f;
};

/** @brief 单点质量及几何约束。
 * @details sourceConfidence 为原始候选置信度（0..1）；RMS、supportDistance、
 * 各 correction 和 positionSigma 使用点云长度单位，RMS/supportDistance=-1 表示不可用。
 * innovationNis 为归一化创新平方；flags 是 SeamQualityFlag 位掩码。
 * tracked/wall 是剖面直线，floor 是检测工作坐标中的平面，方程均为 n.dot(p)=offset。
 * 相应 trackedWall/wallAccepted/floorAccepted 为 true 时法向与偏置才有效。
 */
struct SeamPointQuality {
  // Independent full-profile wall evidence, retained through all refinements.
  bool trackedWall = false;
  Eigen::Vector2f trackedNormal = Eigen::Vector2f::Zero();
  float trackedOffset = 0;
  float sourceConfidence = 0.0f;
  float wallRms = -1.0f;
  float floorRms = -1.0f;
  float supportDistance = -1.0f;
  float geometryCorrection = 0.0f;
  float rtsCorrection = 0.0f;
  float finalCurveCorrection = 0.0f;
  float innovationNis = 0.0f;
  float positionSigma = 0.0f;
  std::uint32_t flags = QualityNone;

  // Accepted local constraints are retained so the optional smoother can be
  // guarded by the actual geometry instead of smoothing XYZ freely.
  bool wallAccepted = false;
  Eigen::Vector2f wallNormal = Eigen::Vector2f::Zero();
  float wallOffset = 0.0f;
  bool floorAccepted = false;
  Eigen::Vector3f floorNormal = Eigen::Vector3f::Zero();
  float floorOffset = 0.0f;
};

/** @brief RTS 统计：applied 表示是否执行，各计数为受对应约束的点数，修正量为长度。 */
struct SeamSmoothingSummary {
  bool applied = false;
  std::size_t innovationOutliers = 0;
  std::size_t correctionLimited = 0;
  std::size_t geometryLimited = 0;
  std::size_t pairConstraintLimited = 0;
  float meanCorrection = 0.0f;
  float maxCorrection = 0.0f;
};

/** @brief 最终曲线拟合统计：异常点/侧壁约束点计数、平均/最大修正和最大侧壁残差（长度）。 */
struct FinalCurveSummary {
  bool applied = false;
  std::size_t outlierPoints = 0;
  float meanCorrection = 0.0f, maxCorrection = 0.0f;
  std::size_t wallConstrainedPoints = 0;
  float maximumWallResidual = 0.0f;
};

/** @brief 历史预测输出；valid=true 时 point 为有效 XYZ，stability 为 0..1 的历史一致性。 */
struct HistoryPrediction {
  bool valid = false;
  Eigen::Vector3f point = Eigen::Vector3f::Zero();
  // 0 表示历史轨迹相互矛盾，1 表示历次轨迹十分稳定。
  float stability = 0.0f;
};

/**
 * @brief 检测输出；轨迹为原始 XYZ 坐标，共享指针拥有输出点云。
 * @details left/right/center 的逐点置信度、质量、帧号按对应轨迹下标关联。
 * frameIndices 仅在 frameBased=true 时有原始帧含义。空轨迹表示未输出该侧。
 * centerSeam 非空表示单道输出；分段起点是轨迹下标，不能跨段连线。
 * 结构的普通复制共享点云；需要可独立修改的副本时由调用方深拷贝点云。
 */
struct SeamResult {
  /// @brief 统一宽度目标，长度单位；未启用时为 0。
  float uniformWidthTarget = 0;
  /// @brief 统一宽度有效样本数和每端去除的样本数。
  std::size_t uniformWidthSamples = 0, uniformWidthTrimEach = 0;
  /// @brief 完成统一宽度调整/因轮廓不足省略的站位数。
  std::size_t uniformWidthAdjusted = 0, uniformWidthOmitted = 0;
  /// @brief 是否执行了统一宽度后的整体滤波。
  bool uniformWidthFilterApplied = false;
  /// @brief 统一宽度滤波的平均/最大修正（长度）。
  float uniformWidthFilterMeanCorrection = 0, uniformWidthFilterMaxCorrection = 0;
  /// @brief 阶段名称到耗时的映射，单位毫秒；调试回调等待时间可能计入。
  std::map<std::string, double> stageMilliseconds;
  /// @brief 超过候选修正约束而拒绝的点数。
  std::size_t candidateCorrectionRejected = 0;
  /// @brief 首部/尾部截去的站位数。
  std::size_t endpointTrimmedStart = 0, endpointTrimmedEnd = 0;
  /// @brief 主体提取后的实测点云。
  CloudT::Ptr mainBodyCloud{new CloudT};
  /// @brief 左侧壁实测 XYZ 点云；原始坐标和输入顺序，不含拟合生成点。
  CloudT::Ptr leftSurface{new CloudT};
  /// @brief 右侧壁实测 XYZ 点云；与左侧壁、焊道面的输入点索引互斥。
  CloudT::Ptr rightSurface{new CloudT};
  /// @brief 左右侧壁拟合面之间的焊道面实测 XYZ 点云（不是中心轨迹）。
  CloudT::Ptr weldSurface{new CloudT};
  /// @brief 三次 B 样条曲面；独立曲面片对应连续支撑区段，缺口不相连。
  BSplineSurfaceModel leftSurfaceModel, weldSurfaceModel, rightSurfaceModel;
  /// @brief 主体提取移除的输入点数。
  std::size_t removedBodyPoints = 0;
  /// @brief true 表示按输入原始帧检测，false 表示空间切片。
  bool frameBased = false;
  /// @brief 输入原始帧总数；空间切片模式为 0。
  std::size_t inputFrameCount = 0;
  /// @brief 左/右轨迹点对应的原始帧号。
  std::vector<int> leftFrameIndices, rightFrameIndices;
  /// @brief 预处理后的点云。
  CloudT::Ptr filteredCloud{new CloudT};
  /// @brief 左侧焊缝 XYZ 点云；保持算法输出顺序。
  CloudT::Ptr leftSeam{new CloudT};
  /// @brief 右侧焊缝 XYZ 点云。
  CloudT::Ptr rightSeam{new CloudT};
  /// @brief 单道模式的中心焊缝 XYZ 点云。
  CloudT::Ptr centerSeam{new CloudT};
  /// @brief 中心轨迹对应的原始帧号。
  std::vector<int> centerFrameIndices;
  /// @brief 中心轨迹逐点置信度，0..1。
  std::vector<float> centerPointConfidence;
  /// @brief 中心轨迹逐点质量诊断。
  std::vector<SeamPointQuality> centerPointQuality;
  /// @brief 中心轨迹各独立连续段的起点下标。
  std::vector<std::size_t> centerSegmentStarts;
  /// @brief 检测得到的最小/最大槽宽（长度）。
  float minimumGrooveWidth = 0, maximumGrooveWidth = 0;
  /// @brief 左轨迹逐点置信度，0..1。
  std::vector<float> leftPointConfidence;
  /// @brief 右轨迹逐点置信度，0..1。
  std::vector<float> rightPointConfidence;
  /// @brief 左轨迹逐点质量诊断。
  std::vector<SeamPointQuality> leftPointQuality;
  /// @brief 右轨迹逐点质量诊断。
  std::vector<SeamPointQuality> rightPointQuality;
  // 独立轨迹段的起始点位置；绘图不得跨越段边界连线。
  std::vector<std::size_t> leftSegmentStarts, rightSegmentStarts; ///< 左/右轨迹独立段起点下标。
  float leftConfidence = 0.0f; ///< 左轨迹整体置信度，0..1。
  float rightConfidence = 0.0f; ///< 右轨迹整体置信度，0..1。
  float confidence = 0.0f; ///< 总体置信度，0..1。
  float leftHistoryUsageRatio = 0.0f; ///< 左轨迹历史使用比例，0..1。
  float rightHistoryUsageRatio = 0.0f; ///< 右轨迹历史使用比例，0..1。
  SeamSmoothingSummary leftSmoothing, rightSmoothing; ///< 左/右 RTS 统计。
  FinalCurveSummary leftFinalCurve, rightFinalCurve; ///< 左/右最终曲线平滑统计。
  std::string summary; ///< 成功概述或输入无效/检测失败原因。

  /** @return 帧模式至少 1 个有效站位、空间切片至少 3 个有效站位时为 true。
   * 单道计中心点数，双道计两侧较小点数；不代表所有输入帧均有输出。 */
  bool ok() const {
    const std::size_t minimum = frameBased ? 1 : 3;
    return (centerSeam ? centerSeam->size() : 0) +
           (leftSeam && rightSeam ? std::min(leftSeam->size(),rightSeam->size()) : 0) >= minimum;
  }
};

inline float coordinate(const PointT &point, Axis axis) {
  switch (axis) {
  case Axis::X:
    return point.x;
  case Axis::Y:
    return point.y;
  case Axis::Z:
    return point.z;
  }
  return point.x;
}

inline float coordinate(const Eigen::Vector3f &point, Axis axis) {
  return point[static_cast<int>(axis)];
}

inline Eigen::Vector3f toVector(const PointT &point) {
  return Eigen::Vector3f(point.x, point.y, point.z);
}

inline PointT toPoint(const Eigen::Vector3f &point) {
  return PointT(point.x(), point.y(), point.z());
}

} // namespace weld
