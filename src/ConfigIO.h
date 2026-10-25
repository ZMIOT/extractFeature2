#pragma once
#include "AlgorithmParameters.h"
#include <string>

namespace weld {
/** @brief 可持久化配置；两种算法拥有独立参数，默认使用连续曲率。 */
struct LibraryConfig {
  ProfileLocalizationMethod activeAlgorithm = ProfileLocalizationMethod::CurvatureCandidates; ///< 默认运行算法。
  AlgorithmParameters algorithms; ///< 独立算法参数快照，见 DetectorConfig 各字段。
};

/**
 * @brief 从本地 UTF-8 JSON 读取两套算法参数。
 * @param[in] path UTF-8 文件路径；相对路径以进程工作目录为基准。
 * @param[out] config 成功时替换为内置默认值叠加文件配置；失败时保持原值。
 * @param[out] error 可空；失败原因含配置键名，成功时清空。
 * @return 文件、格式、键名、枚举及数值范围全部有效时为 true。
 * @details schemaVersion=1。缺失参数使用对应算法默认值；未知键及非法值报错。
 * 不自动搜索文件、不改写文件；结果导出的 *_params.json 不属于此配置格式。
 */
bool loadLibraryConfig(const std::string &path, LibraryConfig &config,
                       std::string *error = nullptr);

/**
 * @brief 原子保存完整配置，包括两种算法的所有可配置字段。
 * @param[in] path UTF-8 输出路径，父目录必须存在；成功后替换原文件。
 * @param[in] config 待保存参数；内部 profileGridOrigin 不持久化。
 * @param[out] error 可空；非法配置或文件写入失败的原因，成功时清空。
 * @return 参数有效且文件提交成功时为 true。
 */
bool saveLibraryConfig(const std::string &path, const LibraryConfig &config,
                       std::string *error = nullptr);

/**
 * @brief 检查单次调用参数；建议直接调用 detect 前使用。
 * @param[in] config 参数快照。
 * @param[out] error 可空；首个无效字段及允许范围，成功时清空。
 * @return 所有可配置字段有效时为 true；不校验点云适用性。
 */
bool validateDetectorConfig(const DetectorConfig &config, std::string *error = nullptr);
} // namespace weld
