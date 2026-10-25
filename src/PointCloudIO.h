#pragma once

#include "SeamTypes.h"

#include <string>
#include <vector>

namespace weld {

class PointCloudIO {
public:
  /**
   * @brief 加载 PCD/PLY/TXT/XYZ/CSV 点云。
   * @param[in] path 文件路径；Windows 使用本地代码页，建议 ASCII 路径。
   * @param[out] error 可空；失败原因，成功时不保证清空。
   * @param[out] frames 可空；开始时清空，识别到帧字段时填充连续原始帧范围。
   * @return 点云共享指针；失败时为空点云，调用方应检查 empty()。
   * @details 有 x/y/z 表头时按列读取，无表头时使用前三列；保持有效点顺序。
   */
  // TXT/XYZ 有 x/y/z 表头时按列名读取（支持 index x y z）；
  // 无表头时仍取前三列 XYZ。保持有效点的原始行顺序，不排序、不重采样。
  static CloudT::Ptr load(const std::string &path, std::string *error = nullptr,
                          std::vector<ScanFrameRange> *frames = nullptr);

  /** @brief 保存按原始帧标识的轨迹 CSV。
   * @param[in] path 输出路径；父目录须存在，覆盖已有文件。
   * @param[in] result 含轨迹及对应 frameIndices 的帧模式结果。
   * @param[out] error 可空，失败原因。 @return 写入成功为 true。 */
  static bool saveFrameSeamCsv(const std::string &path,
                               const SeamResult &result,
                               std::string *error = nullptr);
  /** @brief 根据 frameBased 自动选择帧 CSV 或站位 CSV。
   * @param[in] path 输出路径；父目录须存在。
   * @param[in] result 检测结果，轨迹指针须非空。
   * @param[in] axis 空间切片模式的站位轴。
   * @param[out] error 可空，失败原因。 @return 写入成功为 true。 */
  static bool saveResultCsv(const std::string &path, const SeamResult &result,
                            Axis axis, std::string *error = nullptr);

  // 将单条焊缝作为独立点云保存。PCD 使用二进制格式，PLY/XYZ 也可用于
  // 与其它点云软件交换；左右焊缝由调用方分别调用本函数。
  /** @brief 保存单条轨迹或普通点云；PCD/PLY 为二进制，TXT/XYZ 为文本。
   * @param[in] path 输出路径，由扩展名选择格式；父目录须存在。
   * @param[in] cloud 非空 XYZ 点云。
   * @param[out] error 可空，失败原因。 @return 写入成功为 true。 */
  static bool savePointCloud(const std::string &path, const CloudT &cloud,
                             std::string *error = nullptr);

  /** @brief 导出单条轨迹的 station,x,y,z,confidence CSV。
   * @param[in] path 输出文件路径。 @param[in] seam XYZ 轨迹。
   * @param[in] confidence 逐点置信度，不足的条目写 0。
   * @param[in] axis 纵向站位轴。 @param[out] error 可空，失败原因。
   * @return 写入成功为 true。 */
  static bool saveSeamCsv(const std::string &path, const CloudT &seam,
                          const std::vector<float> &confidence, Axis axis,
                          std::string *error = nullptr);

  // 左右焊缝保存在同一文件中，side 列明确区分 L/R，历史导入时不会混线。
  /** @brief 保存带 L/R 标记的双轨迹 CSV，避免历史导入时混线。
   * @param[in] path 输出文件路径。 @param[in] leftSeam 左轨迹。
   * @param[in] leftConfidence 左轨迹逐点置信度，不足的条目写 0。
   * @param[in] rightSeam 右轨迹。 @param[in] rightConfidence 右轨迹置信度。
   * @param[in] axis 纵向轴。 @param[out] error 可空，失败原因。
   * @return 写入成功为 true。 */
  static bool saveDualSeamCsv(const std::string &path, const CloudT &leftSeam,
                              const std::vector<float> &leftConfidence,
                              const CloudT &rightSeam,
                              const std::vector<float> &rightConfidence,
                              Axis axis, std::string *error = nullptr);

  /** @brief 读取带 side 列的左右历史 CSV。
   * @param[in] path 输入路径。 @param[out] leftSeam 左轨迹指针地址，必须非空。
   * @param[out] rightSeam 右轨迹指针地址，必须非空。
   * @param[out] error 可空，失败原因。 @return 读取成功为 true。 */
  static bool loadDualSeamCsv(const std::string &path, CloudT::Ptr *leftSeam,
                              CloudT::Ptr *rightSeam,
                              std::string *error = nullptr);

  /** @brief 读取单条历史 CSV，仅使用 XYZ，兼容首列 station。
   * @param[in] path 输入文件路径。 @param[out] error 可空，失败原因。
   * @return 读取到的 XYZ 轨迹；失败时返回空点云。 */
  static CloudT::Ptr loadSeamCsv(const std::string &path,
                                 std::string *error = nullptr);

private:
  static CloudT::Ptr loadText(const std::string &path, std::string *error,
                              std::vector<ScanFrameRange> *frames);
};

} // namespace weld
