#include "SeamDebugViewer.h"

#include <QApplication>
#include <QCheckBox>
#include <QDialog>
#include <QHBoxLayout>
#include <QLabel>
#include <QProgressBar>
#include <QPushButton>
#include <QSignalBlocker>
#include <QSpinBox>
#include <QThread>
#include <QTimer>
#include <QVBoxLayout>
#include <QVTKOpenGLNativeWidget.h>
#include <pcl/common/point_tests.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderer.h>

namespace weld {
namespace {
std::array<int, 2> debugTransverseAxes(Axis longitudinal) {
  if (longitudinal == Axis::Y) return {{0, 2}};
  if (longitudinal == Axis::Z) return {{0, 1}};
  return {{1, 2}};
}

CloudT::Ptr debugCloud(const std::vector<PointT> &points) {
  CloudT::Ptr cloud(new CloudT);
  cloud->reserve(points.size());
  for (const auto &point : points)
    cloud->push_back(point);
  return cloud;
}

const char *profileShapeName(const DebugPair &pair) {
  if (pair.profileShape == "shallow_groove")
    return "浅槽（取实测下方壁脚，不使用槽口上沿）";
  if (pair.profileShape == "ambiguous_bead")
    return "凸起判定存在歧义（保留外侧交点供连续路径选择）";
  if (pair.profileShape == "central_bead" ||
      pair.source == "central_bead_valleys")
    return "中央凸起（取两侧低谷）";
  if (pair.profileShape == "wavy_floor" ||
      pair.source == "curve_line_wavy_floor")
    return "波动槽底（取外侧候选）";
  if (pair.profileShape == "flat_floor" ||
      pair.source == "curve_line_flat_floor")
    return "平直/圆滑槽底（取外侧候选）";
  return nullptr;
}

DebugStep historicalFrameStep(const DebugTrace &trace, int ordinal, int total) {
  DebugStep step;
  if (ordinal < 1 || std::size_t(ordinal) > trace.slices.size())
    return step;
  const auto &slice = trace.slices[std::size_t(ordinal - 1)];
  const bool framed = slice.frameIndex >= 0;
  auto raw = debugCloud(slice.raw);
  auto position = debugCloud(slice.sourceRaw.empty() ? slice.raw : slice.sourceRaw);
  CloudT::Ptr profile(new CloudT), left(new CloudT), right(new CloudT);
  const auto axes = debugTransverseAxes(trace.config.longitudinalAxis);
  for (const auto &point : slice.profile) {
    Eigen::Vector3f xyz = Eigen::Vector3f::Zero();
    xyz[int(trace.config.longitudinalAxis)] = slice.station;
    xyz[axes[0]] = point.u;
    xyz[axes[1]] = point.v;
    profile->push_back(toPoint(xyz));
  }
  for (const auto &pair : slice.pairs) {
    left->push_back(pair.points[0]);
    right->push_back(pair.points[1]);
  }
  step.title = std::string("F 回看") + (framed ? "帧" : "截面") +
               " · 第 " + std::to_string(ordinal) + "–" +
               std::to_string(ordinal) + " / " + std::to_string(total) +
               (framed ? " · 帧号 " + std::to_string(slice.frameIndex) : "");
  step.message = "重新显示已执行的切片；左侧黄色点为原始三维位置，右侧为算法内部截面。";
  if (slice.pairs.empty() && !slice.noCandidateReason.empty())
    step.message += " 【无候选原因】" + slice.noCandidateReason;
  if (slice.selected >= 0 &&
      slice.selected < static_cast<int>(slice.pairs.size()))
    if (const char *shape = profileShapeName(slice.pairs[slice.selected]))
      step.message += std::string(" 轮廓判定：") + shape + "。";
  step.metrics = {{"first_frame", double(framed ? slice.frameIndex : ordinal)},
                  {"last_frame", double(framed ? slice.frameIndex : ordinal)},
                  {"first_ordinal", double(ordinal)},
                  {"last_ordinal", double(ordinal)},
                  {"total_slices", double(total)},
                  {"progress_percent", total ? 100.0 * ordinal / total : 100.0},
                  {"frame_based", framed ? 1.0 : 0.0},
                  {"first_station", double(slice.station)},
                  {"last_station", double(slice.station)},
                  {"batch_count", 1.0},
                  {"without_candidates", slice.pairs.empty() ? 1.0 : 0.0}};
  step.layers = {{"主体截面", raw},
                 {"分箱轮廓", profile},
                 {"左候选", left},
                 {"右候选", right},
                 {"B 样条拟合", debugCloud(slice.fittedProfile)}};
  CloudT::Ptr overviewLeft(new CloudT), overviewRight(new CloudT);
  for (int i = 0; i < ordinal; ++i) {
    const auto &executed = trace.slices[std::size_t(i)];
    for (const auto &candidate : executed.sourceCandidates[0])
      overviewLeft->push_back(candidate);
    for (const auto &candidate : executed.sourceCandidates[1])
      overviewRight->push_back(candidate);
  }
  step.overviewCandidates = {{overviewLeft, overviewRight}};
  step.overviewCloud = trace.clouds.empty() ? CloudT::ConstPtr()
                                             : trace.clouds.front().cloud;
  step.positionCloud = position;
  step.hasFrameView = true;
  const int axis = int(trace.config.longitudinalAxis);
  step.viewU = Eigen::Vector3f::Zero();
  step.viewV = Eigen::Vector3f::Zero();
  step.viewU[axis == 0 ? 1 : 0] = 1;
  step.viewV[axis == 2 ? 1 : 2] = 1;
  if (!raw->empty()) {
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    for (const auto &point : *raw)
      center += toVector(point).cast<double>();
    step.viewOrigin = (center / double(raw->size())).cast<float>();
  }
  return step;
}

// 复用同一 VTK 窗口，每一步等待用户确认；不会累积数百个窗口。
class DebugViewer : public QDialog {
public:
  explicit DebugViewer(QWidget *parent) : QDialog(parent) {
    setObjectName("seamDebugViewer");
    resize(1400, 800);
    auto *layout = new QVBoxLayout(this);
    description_ = new QLabel(this);
    description_->setWordWrap(true);
    layout->addWidget(description_);
    progressLabel_ = new QLabel(this);
    progressLabel_->setObjectName("debugPositionLabel");
    progressLabel_->setStyleSheet("font-size:15px;font-weight:600;");
    layout->addWidget(progressLabel_);
    progress_ = new QProgressBar(this);
    progress_->setObjectName("debugPositionProgress");
    progress_->setRange(0, 1000);
    progress_->setTextVisible(false);
    layout->addWidget(progress_);
    auto *views = new QHBoxLayout;
    overviewPanel_ = new QWidget(this);
    auto *overviewLayout = new QVBoxLayout(overviewPanel_);
    overviewLayout->setContentsMargins(0, 0, 0, 0);
    auto *overviewTitle = new QLabel(
        QStringLiteral("整云位置（黄色为当前切片，绿色/红色为累计候选轨迹）"),
        overviewPanel_);
    overviewTitle->setStyleSheet("font-weight:600;");
    overviewLayout->addWidget(overviewTitle);
    overviewWidget_ = new QVTKOpenGLNativeWidget(overviewPanel_);
    overviewWidget_->setMinimumSize(400, 480);
    overviewLayout->addWidget(overviewWidget_, 1);
    auto *detailPanel = new QWidget(this);
    auto *detailLayout = new QVBoxLayout(detailPanel);
    detailLayout->setContentsMargins(0, 0, 0, 0);
    auto *detailTitle = new QLabel(QStringLiteral("当前切片轮廓与候选"), detailPanel);
    detailTitle->setStyleSheet("font-weight:600;");
    detailLayout->addWidget(detailTitle);
    detailWidget_ = new QVTKOpenGLNativeWidget(detailPanel);
    detailWidget_->setMinimumSize(700, 480);
    detailLayout->addWidget(detailWidget_, 1);
    views->addWidget(overviewPanel_, 2);
    views->addWidget(detailPanel, 3);
    layout->addLayout(views, 1);
    legend_ = new QLabel(this);
    legend_->setWordWrap(true);
    layout->addWidget(legend_);
    auto *buttons = new QHBoxLayout;
    keepCamera_ = new QCheckBox(QStringLiteral("保持当前视角"), this);
    buttons->addWidget(keepCamera_);
    jumpLabel_ = new QLabel(QStringLiteral("跳转到第"), this);
    frameSelector_ = new QSpinBox(this);
    frameSelector_->setObjectName("debugFrameSelector");
    frameSelector_->setRange(1, 1);
    jump_ = new QPushButton(QStringLiteral("跳转帧"), this);
    jump_->setObjectName("jumpDebugFrame");
    jump_->setEnabled(false);
    buttons->addWidget(jumpLabel_);
    buttons->addWidget(frameSelector_);
    buttons->addWidget(new QLabel(QStringLiteral("帧"), this));
    buttons->addWidget(jump_);
    buttons->addStretch();
    auto *skip = new QPushButton(QStringLiteral("跳过后续弹窗"), this);
    skip->setObjectName("skipDebugSteps");
    next_ = new QPushButton(QStringLiteral("下一步"), this);
    auto *next = next_;
    next->setObjectName("nextDebugStep");
    next->setDefault(true);
    buttons->addWidget(skip);
    buttons->addWidget(next);
    layout->addLayout(buttons);
    connect(skip, &QPushButton::clicked, this, &QDialog::reject);
    connect(next, &QPushButton::clicked, this, [this] {
      if (currentOrdinal_ > 0 && currentOrdinal_ < availableOrdinal_)
        requestedOrdinal_ = currentOrdinal_ + 1;
      accept();
    });
    connect(frameSelector_, QOverload<int>::of(&QSpinBox::valueChanged), this,
            [this](int value) { jump_->setEnabled(value != currentOrdinal_); });
    connect(jump_, &QPushButton::clicked, this, [this] {
      requestedOrdinal_ = frameSelector_->value();
      accept();
    });

    auto initialize = [](QVTKOpenGLNativeWidget *widget,
                         pcl::visualization::PCLVisualizer::Ptr &viewer,
                         const char *name) {
      auto renderer = vtkSmartPointer<vtkRenderer>::New();
      auto window = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
      window->AddRenderer(renderer);
      widget->SetRenderWindow(window);
      viewer.reset(new pcl::visualization::PCLVisualizer(renderer, window,
                                                          name, false));
      viewer->setupInteractor(widget->GetInteractor(), window);
      viewer->setBackgroundColor(.025, .055, .075);
    };
    initialize(overviewWidget_, overviewViewer_, "seam_debug_overview");
    initialize(detailWidget_, detailViewer_, "seam_debug_detail");
  }

  void setHistoryProvider(std::function<DebugStep(int, int)> provider) {
    historyProvider_ = std::move(provider);
  }

  bool handleStep(const DebugStep &step) {
    if (step.title.compare(0, 1, "F") == 0) {
      const auto last = step.metrics.find("last_ordinal");
      if (last != step.metrics.end())
        availableOrdinal_ = std::max(availableOrdinal_, int(last->second));
      const auto total = step.metrics.find("total_slices");
      if (total != step.metrics.end())
        totalOrdinal_ = int(total->second);
    }
    if (requestedOrdinal_ > 0 && step.title.compare(0, 1, "F") == 0) {
      const auto last = step.metrics.find("last_ordinal");
      if (last != step.metrics.end() && int(last->second) < requestedOrdinal_)
        return true;
      requestedOrdinal_ = 0;
    }
    bool keep = showStep(step);
    while (keep && requestedOrdinal_ > 0 &&
           requestedOrdinal_ <= availableOrdinal_ && historyProvider_) {
      const int requested = requestedOrdinal_;
      requestedOrdinal_ = 0;
      const auto historical = historyProvider_(requested, totalOrdinal_);
      if (historical.title.empty())
        break;
      keep = showStep(historical);
    }
    return keep;
  }

  bool showStep(const DebugStep &step) {
    next_->setText(step.title.compare(0, 1, "F") == 0
                       ? QStringLiteral("下一组帧 / 截面") : QStringLiteral("下一步"));
    setWindowTitle(QStringLiteral("点云调试 · %1 · %2")
                       .arg(++stepNumber_)
                       .arg(QString::fromStdString(step.title)));
    description_->setText(
        QString::fromStdString(step.title + "\n" + step.message));
    const auto totalIt = step.metrics.find("total_slices");
    const auto currentIt = step.metrics.find("last_ordinal");
    if (step.title.compare(0, 1, "F") == 0 && totalIt != step.metrics.end() &&
        currentIt != step.metrics.end()) {
      const double total = totalIt->second;
      const double current = currentIt->second;
      const double percent = total > 0 ? 100.0 * current / total : 100.0;
      const auto framedIt = step.metrics.find("frame_based");
      const bool framed = framedIt != step.metrics.end() && framedIt->second > 0.5;
      QString position = framed
          ? QStringLiteral("当前位置：第 %1 / %2 帧（%3%）")
                .arg(int(current)).arg(int(total)).arg(percent, 0, 'f', 1)
          : QStringLiteral("当前位置：第 %1 / %2 个截面（%3%）")
                .arg(int(current)).arg(int(total)).arg(percent, 0, 'f', 1);
      const auto frameIt = step.metrics.find("last_frame");
      if (frameIt != step.metrics.end() && framedIt != step.metrics.end() &&
          framed)
        position += QStringLiteral(" · 当前帧号 %1").arg(int(frameIt->second));
      const auto stationIt = step.metrics.find("last_station");
      if (stationIt != step.metrics.end())
        position += QStringLiteral(" · 纵向坐标 %1").arg(stationIt->second, 0, 'f', 2);
      progressLabel_->setText(position);
      progress_->setValue(int(std::min(1000.0, std::max(0.0, percent * 10.0))));
      progressLabel_->show();
      progress_->show();
      currentOrdinal_ = int(current);
      const QSignalBlocker blocker(frameSelector_);
      frameSelector_->setRange(1, int(total));
      frameSelector_->setValue(currentOrdinal_);
      frameSelector_->setSuffix(QStringLiteral(" / %1").arg(int(total)));
      jumpLabel_->show();
      frameSelector_->show();
      jump_->show();
      jump_->setEnabled(false);
    } else {
      progressLabel_->hide();
      progress_->hide();
      jumpLabel_->hide();
      frameSelector_->hide();
      jump_->hide();
    }
    detailViewer_->removeAllPointClouds();
    detailViewer_->removeAllShapes();
    overviewViewer_->removeAllPointClouds();
    overviewViewer_->removeAllShapes();
    overviewPanel_->setVisible(bool(step.overviewCloud));
    if (step.overviewCloud && !step.overviewCloud->empty()) {
      pcl::visualization::PointCloudColorHandlerCustom<PointT> wholeColor(
          step.overviewCloud, 65, 115, 140);
      overviewViewer_->addPointCloud<PointT>(step.overviewCloud, wholeColor, "whole_cloud");
      overviewViewer_->setPointCloudRenderingProperties(
          pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1.0, "whole_cloud");
      if (step.positionCloud && !step.positionCloud->empty()) {
        pcl::visualization::PointCloudColorHandlerCustom<PointT> sliceColor(
            step.positionCloud, 255, 220, 35);
        overviewViewer_->addPointCloud<PointT>(step.positionCloud, sliceColor, "current_slice");
        overviewViewer_->setPointCloudRenderingProperties(
            pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6.0, "current_slice");
      }
      const int candidateColors[2][3] = {{70, 255, 105}, {255, 95, 85}};
      for (int side = 0; side < 2; ++side) {
        const auto &candidates = step.overviewCandidates[side];
        if (!candidates || candidates->empty())
          continue;
        const std::string id = side ? "overview_right_candidates"
                                    : "overview_left_candidates";
        pcl::visualization::PointCloudColorHandlerCustom<PointT> color(
            candidates, candidateColors[side][0], candidateColors[side][1],
            candidateColors[side][2]);
        overviewViewer_->addPointCloud<PointT>(candidates, color, id);
        overviewViewer_->setPointCloudRenderingProperties(
            pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 8.0, id);
      }
      if (!overviewCameraReady_) {
        overviewViewer_->resetCamera();
        overviewCameraReady_ = true;
      }
    }
    QStringList legend;
    const int colors[][3] = {{80, 170, 215},  {255, 220, 50},  {80, 255, 110},
                             {255, 110, 100}, {205, 130, 255}, {250, 250, 250}};
    std::size_t index = 0;
    float frameRadius = 1.f;
    for (const auto &layer : step.layers) {
      CloudT::Ptr finite(new CloudT);
      if (layer.cloud) {
        finite->reserve(layer.cloud->size());
        for (const auto &p : *layer.cloud)
          if (pcl::isFinite(p)) {
            finite->push_back(p);
            if (step.hasFrameView)
              frameRadius =
                  std::max(frameRadius, (toVector(p) - step.viewOrigin).norm());
          }
      }
      const int *rgb = colors[index % 6];
      legend << QStringLiteral(
                    "<span style='color:rgb(%1,%2,%3)'>● %4 (%5 点)</span>")
                    .arg(rgb[0])
                    .arg(rgb[1])
                    .arg(rgb[2])
                    .arg(QString::fromStdString(layer.name).toHtmlEscaped())
                    .arg(finite->size());
      if (!finite->empty()) {
        const std::string id = "layer_" + std::to_string(index);
        pcl::visualization::PointCloudColorHandlerCustom<PointT> color(
            finite, rgb[0], rgb[1], rgb[2]);
        detailViewer_->addPointCloud<PointT>(finite, color, id);
        const bool candidateLayer = layer.name.find("候选") != std::string::npos ||
                                    step.title.compare(0, 4, "S06 ") == 0;
        const double pointSize = candidateLayer ? 15.0
            : (layer.name == "B 样条拟合" ? 2.0 : (index ? 5.0 : 2.0));
        detailViewer_->setPointCloudRenderingProperties(
            pcl::visualization::PCL_VISUALIZER_POINT_SIZE, pointSize, id);
      }
      ++index;
    }
    legend_->setText(legend.join(" &nbsp; "));
    if (!keepCamera_->isChecked() || stepNumber_ == 1) {
      detailViewer_->resetCamera();
      if (step.hasFrameView && step.viewU.cross(step.viewV).norm() > .1f) {
        const Eigen::Vector3f eye =
            step.viewOrigin +
            step.viewU.cross(step.viewV).normalized() * (3.f * frameRadius);
        detailViewer_->setCameraPosition(eye.x(), eye.y(), eye.z(),
                                   step.viewOrigin.x(), step.viewOrigin.y(),
                                   step.viewOrigin.z(), step.viewV.x(),
                                   step.viewV.y(), step.viewV.z());
      }
    }
    detailWidget_->update();
    overviewWidget_->update();
    // exec() re-shows the reused dialog. Render with its visible GL context;
    // QWidget::update alone can keep the previous frame's VTK framebuffer.
    QTimer::singleShot(0, this, [this] {
      detailWidget_->GetRenderWindow()->Render();
      detailWidget_->update();
      if (overviewPanel_->isVisible()) {
        overviewWidget_->GetRenderWindow()->Render();
        overviewWidget_->update();
      }
    });
    // 关闭窗口或 Esc 与“跳过后续弹窗”一致，仍完成后续算法。
    return exec() == QDialog::Accepted;
  }

private:
  QWidget *overviewPanel_ = nullptr;
  QVTKOpenGLNativeWidget *overviewWidget_ = nullptr, *detailWidget_ = nullptr;
  QLabel *description_ = nullptr, *progressLabel_ = nullptr, *legend_ = nullptr;
  QLabel *jumpLabel_ = nullptr;
  QProgressBar *progress_ = nullptr;
  QSpinBox *frameSelector_ = nullptr;
  QCheckBox *keepCamera_ = nullptr;
  QPushButton *next_ = nullptr, *jump_ = nullptr;
  pcl::visualization::PCLVisualizer::Ptr overviewViewer_, detailViewer_;
  std::function<DebugStep(int, int)> historyProvider_;
  bool overviewCameraReady_ = false;
  int stepNumber_ = 0, currentOrdinal_ = 0, requestedOrdinal_ = 0,
      availableOrdinal_ = 0, totalOrdinal_ = 0;
};
} // namespace

SeamResult detectWithDebugViewer(const CloudT::ConstPtr &input,
                                 const DetectorConfig &config, QWidget *parent,
                                 const std::vector<ScanFrameRange> &frames) {
  if (!config.debugFlag)
    return SeamDetector().detect(input, config, nullptr, nullptr, nullptr,
                                 &frames);
  if (!qobject_cast<QApplication *>(QCoreApplication::instance()) ||
      QThread::currentThread() != QCoreApplication::instance()->thread())
    throw std::runtime_error(
        "VTK step debugging requires the QApplication GUI thread.");
  DebugViewer viewer(parent);
  DebugTrace trace;
  viewer.setHistoryProvider([&trace](int ordinal, int total) {
    return historicalFrameStep(trace, ordinal, total);
  });
  trace.onStep = [&viewer](const DebugStep &step) {
    return viewer.handleStep(step);
  };
  const auto result = SeamDetector().detect(input, config, nullptr, nullptr, &trace,
                                            &frames);
  if (!result.ok() && trace.onStep) {
    DebugStep failed;
    failed.title = "本次提取未生成可靠焊缝";
    failed.message = result.summary;
    failed.layers = {{"主体点云", result.mainBodyCloud}};
    viewer.showStep(failed);
  }
  return result;
}
} // namespace weld
