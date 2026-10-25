#pragma once

#include "SeamDebug.h"
#include "SeamHistory.h"
#include "SeamTypes.h"

namespace weld {

/** @brief 无持久状态的焊缝检测器；每次调用返回独立结果。 */
class SeamDetector {
public:
  /**
   * @brief 提取焊缝轨迹及左侧壁、右侧壁、焊道面实测点云。
   * @param[in] input 原始 XYZ 点云，只读；至少 30 点。原始帧模式要求全部坐标有限。
   * @param[in] config 本次参数快照，长度单位与点云一致，见 DetectorConfig。
   * @param[in] leftHistory 左侧历史轨迹；nullptr 不使用历史，不在调用内更新。
   * @param[in] rightHistory 右侧历史轨迹；nullptr 不使用历史。
   * @param[in,out] debug 可空；重置并填充过程信息，保留 onStep 回调；回调在调用线程同步执行。
   * @param[in] frames 可空/空数组表示空间切片；否则范围必须连续覆盖 input，count>0，
   * 帧号非负且严格递增。原始帧模式忽略 sliceStep/sliceThickness。
   * @return 包含 XYZ 轨迹、三个面实测点云、置信度、质量标志、分段边界及耗时的 SeamResult。
   * 先检查 ok()；失败或未检测到焊缝时读取 summary，不使用空轨迹。
   * @note 输出坐标恢复到原始点云坐标系；不得跨 segmentStarts 的边界连线。
   * 可并发使用，但每次调用的 debug 必须独立，输入/历史在调用期间不得被修改。
   * 分配失败、依赖库及用户回调抛出的异常由调用方处理。
   */
  SeamResult detect(const CloudT::ConstPtr &input, const DetectorConfig &config,
                    const SeamHistory *leftHistory = nullptr,
                    const SeamHistory *rightHistory = nullptr,
                    DebugTrace *debug = nullptr,
                    const std::vector<ScanFrameRange> *frames = nullptr) const;

private:
  static CloudT::Ptr preprocess(const CloudT::ConstPtr &input,
                                const DetectorConfig &config,
                                DebugTrace *debug);
};

} // namespace weld
