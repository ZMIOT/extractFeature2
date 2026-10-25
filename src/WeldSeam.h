#pragma once
/**
 * @file WeldSeam.h
 * @brief 点云输入，返回焊缝轨迹及左侧壁、右侧壁、焊道面点云。
 * @details 链接 WeldSeam::WeldSeam；C++14，PCL 1.11 common/filters 及 Eigen。
 * 长度与输入点云单位一致（默认参数按毫米），角度为度。
 * 不依赖 Qt、VTK 或 pcl_io，不打开窗口，不读取/保存点云文件。
 * Windows 调用方须使用与 SDK 匹配的 x64 MSVC/v142、Release /MD 及依赖版本。
 */
#include "AlgorithmParameters.h"
#include "ConfigIO.h"

namespace weld {
/** @brief 输出焊缝轨迹及三个面；双道轨迹为 left/right，单道轨迹为 center。
 * @details 坐标与输入 XYZ 使用相同坐标系和单位。无对应轨迹时为空点云。
 * 共享指针拥有输出数据；默认复制共享点云，需独立修改时由调用方深拷贝。
 * 轨迹点数量可能小于输入帧数，不应跨缺口强行连线。
 * 三个面均为输入实测点的子集，保持原始 XYZ、单位及输入顺序，各面不重复使用
 * 同一个输入点；不会把拟合点或体素中心作为面点输出。不保证覆盖全部输入。
 * 仅输出有效检测截面内有局部几何支撑的区域，缺失/不支持的面可能为空。
 * 左右以检测横向坐标增大方向定义（默认 Y 较小为左），与观察视角无关。
 * 三个新增 SurfaceModel 为分片三次 B 样条函数，须单独检查 valid()。
 * SDK 5 扩展了参数与检测统计布局；调用方必须与新版头文件、LIB、DLL 一起重新编译。
 */
struct SeamClouds {
  CloudT::Ptr left{new CloudT}; ///< 左焊缝 XYZ 点云。
  CloudT::Ptr right{new CloudT}; ///< 右焊缝 XYZ 点云。
  CloudT::Ptr center{new CloudT}; ///< 自适应单道模式的中间焊缝 XYZ 点云。
  CloudT::Ptr leftSurface{new CloudT}; ///< 左侧壁实测 XYZ 面点云，排除外部上沿。
  CloudT::Ptr rightSurface{new CloudT}; ///< 右侧壁实测 XYZ 面点云，排除外部上沿。
  CloudT::Ptr weldSurface{new CloudT}; ///< 左右侧壁拟合面之间的焊道面实测点云，不是中心轨迹。
  BSplineSurfaceModel leftSurfaceModel; ///< 左侧曲面函数；分片参数域均为 [0,1]^2。
  BSplineSurfaceModel weldSurfaceModel; ///< 中间焊道曲面函数。
  BSplineSurfaceModel rightSurfaceModel; ///< 右侧曲面函数。
};

/** @brief 返回连续曲率算法的内置默认参数，不读取文件。
 * @return 可直接修改并传给 extractSeams 的独立参数副本。 */
DetectorConfig defaultDetectorConfig();

/**
 * @brief 从输入点云提取焊缝轨迹和三个面的实测点云。
 * @param[in] input 只读 pcl::PointCloud<pcl::PointXYZ>；至少 30 点。
 * @param[out] output 成功时替换全部轨迹、面点云和曲面模型；失败时清空全部输出。
 * 面提取使用单道转换/等宽调整之前的几何边界；单道模式仍可输出三个面。
 * @param[in] config 算法参数；省略时使用连续曲率默认值。长度与输入点云单位一致。
 * @param[out] error 可空；成功时清空，失败时返回配置/输入/检测错误。
 * @param[in] frames 可空；有原始帧时传入完整连续点索引范围，帧号非负且递增。
 * 不传时按 config.longitudinalAxis 空间切片；不会根据点距猜测原始帧。
 * @return 成功提取有效轨迹为 true；输入/配置错误或未找到焊缝为 false。
 * true 不代表三个面都有足够的实测支撑，需要完整三面时检查各面是否为空。
 * 曲面模型至少需要连续四个有效截面，拟合失败不影响原有轨迹和实测面点云。
 * @note 不改变输入、不访问点云文件、不触发任何窗口；多次调用可独立并发，
 * 但输入不得同时被修改，output/error 须由各调用独占。内存不足异常仍可能抛出。
 */
bool extractSeams(const CloudT::ConstPtr &input, SeamClouds &output,
                  const DetectorConfig &config = defaultDetectorConfig(),
                  std::string *error = nullptr,
                  const std::vector<ScanFrameRange> *frames = nullptr);
} // namespace weld
