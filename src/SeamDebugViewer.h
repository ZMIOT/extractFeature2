#pragma once

#include "SeamDetector.h"

class QWidget;

namespace weld {
// GUI/CLI 共用入口。debugFlag=false 时直接执行算法；开启时须在 GUI 线程调用。
SeamResult
detectWithDebugViewer(const CloudT::ConstPtr &input,
                      const DetectorConfig &config, QWidget *parent = nullptr,
                      const std::vector<ScanFrameRange> &frames = {},
                      bool generateSurfaceModels = true);
} // namespace weld
