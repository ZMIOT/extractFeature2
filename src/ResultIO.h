#pragma once
#include "SeamTypes.h"
#include <QJsonObject>
#include <QString>
namespace weld {
/**
 * @brief 生成结果文件的参数元数据（不是 loadLibraryConfig 的配置格式）。
 * @param[in] config 本次检测实际使用的参数快照。
 * @return 含算法标识和运行参数的 JSON 对象；不执行文件操作。
 */
QJsonObject configJson(const DetectorConfig &config);
// Portable full-knot B-spline representation, no PCL/OpenNURBS-specific objects.
QJsonObject surfaceModelJson(const BSplineSurfaceModel &model);
// The selected PCD name supplies the export prefix: _left/_right.pcd
// (or _center.pcd in single mode), plus CSV and parameter JSON. No merged PCD.
// An explicit CSV path remains supported with the same per-track PCD files.
/**
 * @brief 导出独立焊缝 PCD、CSV 和运行参数 JSON，必要时创建父目录。
 * @param[in] outputPath PCD 前缀路径或 CSV 路径；无扩展名按 PCD 处理。
 * @param[in] result 成功的检测结果，须满足 ok()，轨迹指针须非空。
 * @param[in] config 产生 result 时的参数快照，不应传入之后修改的界面参数。
 * @param[out] error 可空；失败时填入原因，成功时不保证清空原字符串。
 * @return 全部文件保存成功为 true；失败时可能已经写出部分文件。
 * @details 生成 *_left.pcd / *_right.pcd 或 *_center.pcd、*_seams.csv、*_params.json。
 * 同时生成三个 *_surface.pcd，以及 *_left_surface_model.json、
 * *_weld_surface_model.json、*_right_surface_model.json。空模型也写入状态及空 patches，
 * 覆盖旧模型，避免残留上次拟合结果。
 * 同名前缀会覆盖旧结果，并删除该前缀下当前无对应轨迹的旧 PCD 文件。
 */
bool saveResults(const QString &outputPath, const SeamResult &result,
                 const DetectorConfig &config, QString *error);
} // namespace weld
