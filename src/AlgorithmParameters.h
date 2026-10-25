#pragma once
#include "SeamTypes.h"
#include <array>

namespace weld {
// A complete, independent configuration per algorithm, including preprocessing,
// path selection and postprocessing. DetectorConfig remains the per-run snapshot.
class AlgorithmParameters {
public:
  /** @brief 初始化两套独立的内置默认参数；曲率算法最小间距 6，代表点为壁脚转角重心。 */
  AlgorithmParameters() { reset(ProfileLocalizationMethod::CurveLineIntersection);
                          reset(ProfileLocalizationMethod::CurvatureCandidates); }
  /** @brief 读取参数。 @param[in] method 有效算法枚举 0 或 1。
   * @return 当前对象内的只读引用，有效期不超过当前对象。 */
  const DetectorConfig &get(ProfileLocalizationMethod method) const {
    return configs_[static_cast<std::size_t>(method)];
  }
  /** @brief 保存指定算法的独立参数副本。
   * @param[in] method 有效算法枚举 0 或 1。
   * @param[in] config 待保存快照，profileLocalizationMethod 自动调整为 method。 */
  void set(ProfileLocalizationMethod method, DetectorConfig config) {
    config.profileLocalizationMethod = method;
    configs_[static_cast<std::size_t>(method)] = std::move(config);
  }
  /** @brief 仅恢复指定算法默认值。 @param[in] method 有效算法枚举 0 或 1。 */
  void reset(ProfileLocalizationMethod method) {
    DetectorConfig config;
    // Narrow rounded grooves need their lower turns, not weaker outer ones.
    if (method == ProfileLocalizationMethod::CurvatureCandidates) {
      config.minSeamSeparation = 6.f;
      config.profileTurnLocation = ProfileTurnLocation::WallTurnCentroid;
    }
    set(method, config);
  }
private:
  std::array<DetectorConfig, 2> configs_;
};
} // namespace weld
