#include "MainWindow.h"
#include "PointCloudIO.h"
#include "ResultIO.h"
#include "SeamDebugViewer.h"
#include "SeamDetector.h"
#include "ConfigIO.h"
#include <QCoreApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDateTime>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGroupBox>
#include <QGridLayout>
#include <QLabel>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QScrollArea>
#include <QSignalBlocker>
#include <QSpinBox>
#include <QStatusBar>
#include <QTimer>
#include <QVBoxLayout>
#include <QVTKOpenGLNativeWidget.h>
#include <QtConcurrent/QtConcurrentRun>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderer.h>

namespace weld {
MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
  setWindowTitle(QStringLiteral("焊缝提取 · 独立算法参数"));
  resize(1280, 860);
  buildUi();
  initializeViewer();
  resetDefaults();
  const QString configPath = QCoreApplication::applicationDirPath() + "/config/weld_seam.json";
  if (QFileInfo::exists(configPath)) {
    LibraryConfig saved;
    std::string error;
    if (loadLibraryConfig(configPath.toUtf8().constData(), saved, &error)) {
      algorithmParameters_ = saved.algorithms;
      activeAlgorithm_ = saved.activeAlgorithm;
      applyConfig(algorithmParameters_.get(activeAlgorithm_));
      appendLog(QStringLiteral("已加载本地配置：") + configPath);
    } else {
      appendLog(QStringLiteral("本地配置无效，使用内置默认值：") + QString::fromUtf8(error.c_str()));
    }
  }
  connect(&watcher_, &QFutureWatcher<SeamResult>::finished, this,
          [this] { finishDetection(watcher_.result()); });
}
void MainWindow::finishDetection(const SeamResult &result) {
  result_ = result;
  busy_ = false;
  parameters_->setEnabled(true);
  updateInputMode();
  defaults_->setEnabled(true);
  load_->setEnabled(true);
  run_->setEnabled(true);
  save_->setEnabled(result_.ok());
  resultLabel_->setText(result_.ok()
                            ? QStringLiteral("左 %1 点 · 右 %2 点\n置信度 %3%")
                                  .arg(result_.leftSeam->size())
                                  .arg(result_.rightSeam->size())
                                  .arg(result_.confidence * 100, 0, 'f', 1)
                            : QStringLiteral("本次未提取到完整焊缝"));
  if (result_.ok() && !result_.centerSeam->empty())
    resultLabel_->setText(QStringLiteral("中间 %1 点 · 左 %2 点 · 右 %3 点\n置信度 %4%")
        .arg(result_.centerSeam->size()).arg(result_.leftSeam->size()).arg(result_.rightSeam->size()).arg(result_.confidence*100,0,'f',1));
  if (result_.ok()) {
    const auto faces=QStringLiteral("左侧面 %1 点 · 右侧面 %2 点\n焊道面 %3 点")
        .arg(result_.leftSurface->size()).arg(result_.rightSurface->size())
        .arg(result_.weldSurface->size());
    resultLabel_->setText(resultLabel_->text()+"\n"+faces);
    appendLog(faces);
    const auto models=QStringLiteral("B 样条曲面片：左 %1 · 中 %2 · 右 %3")
        .arg(result_.leftSurfaceModel.patches.size()).arg(result_.weldSurfaceModel.patches.size())
        .arg(result_.rightSurfaceModel.patches.size());
    resultLabel_->setText(resultLabel_->text()+"\n"+models);appendLog(models);
    if(result_.leftSurface->empty() || result_.rightSurface->empty() || result_.weldSurface->empty())
      appendLog(QStringLiteral("部分面缺少有效实测支撑，未输出的面为空。"));
  }
  appendLog(QString::fromStdString(result_.summary));
  if(resultConfig_.uniformWidthEnabled) {
    const auto detail=QStringLiteral("统一宽度 %1（%2） · 有效 %3，两端各去掉 %4\n调整 %5 帧/站位，轮廓不足省略 %6")
        .arg(result_.uniformWidthTarget,0,'f',3)
        .arg(resultConfig_.uniformWidthMetric==UniformWidthMetric::Chord?QStringLiteral("直线距离"):QStringLiteral("槽底弧长"))
        .arg(result_.uniformWidthSamples).arg(result_.uniformWidthTrimEach)
        .arg(result_.uniformWidthAdjusted).arg(result_.uniformWidthOmitted);
    resultLabel_->setText(resultLabel_->text()+"\n"+detail);
    appendLog(detail);
    if(result_.uniformWidthFilterApplied)
      appendLog(QStringLiteral("统一宽度后滤波：长度 %1，平均移动 %2，最大移动 %3；左右直线宽度保持不变。")
          .arg(resultConfig_.uniformWidthFilterLength).arg(result_.uniformWidthFilterMeanCorrection,0,'f',3)
          .arg(result_.uniformWidthFilterMaxCorrection,0,'f',3));
  }
  for (const auto &stage : result_.stageMilliseconds)
    appendLog(QStringLiteral("耗时 %1：%2 ms").arg(QString::fromStdString(stage.first)).arg(stage.second, 0, 'f', 1));
  statusBar()->showMessage(
      result_.ok() ? QStringLiteral("提取完成：蓝色左侧面、橙色右侧面、绿色焊道面；可在“结果显示”中分别开关。")
                   : QStringLiteral("提取失败，请检查日志与参数。"));
  refreshView(true);
  emit extractionFinished(result_.ok());
}
void MainWindow::setDebugFlag(bool enabled) {
  if (!busy_)
    debugFlag_->setChecked(enabled);
}

MainWindow::~MainWindow() { watcher_.waitForFinished(); }
void MainWindow::buildUi() {
  auto *central = new QWidget(this);
  auto *layout = new QHBoxLayout(central);
  auto *panel = new QWidget(central);
  auto *column = new QVBoxLayout(panel);
  auto *title = new QLabel(QStringLiteral("焊缝提取 · 算法工作台"), panel);
  title->setStyleSheet("font-size:20px;font-weight:600;padding:6px 0;");
  column->addWidget(title);
  inputLabel_ =
      new QLabel(QStringLiteral("请加载 TXT、XYZ、PCD 或 PLY 点云"), panel);
  inputLabel_->setWordWrap(true);
  column->addWidget(inputLabel_);
  load_ = new QPushButton(QStringLiteral("加载点云…"), panel);
  column->addWidget(load_);
  connect(load_, &QPushButton::clicked, this, [this] {
    const auto path = QFileDialog::getOpenFileName(
        this, QStringLiteral("加载点云"),
        inputPath_.isEmpty() ? QDir::currentPath() + "/sample"
                             : QFileInfo(inputPath_).absolutePath(),
        QStringLiteral("点云 (*.txt *.xyz *.pcd *.ply);;所有文件 (*)"));
    if (!path.isEmpty())
      loadPath(path);
  });
  parameters_ = new QGroupBox(QStringLiteral("关键参数"), panel);
  auto *form = new QFormLayout(parameters_);
  auto spin = [this](double low, double high, double increment,
                     const QString &tip) {
    auto *s = new QDoubleSpinBox(parameters_);
    s->setRange(low, high);
    s->setDecimals(2);
    s->setSingleStep(increment);
    s->setToolTip(tip);
    connect(s, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
            [this](double) { invalidateResult(); });
    return s;
  };
  axis_ = new QComboBox(parameters_);
  axis_->addItems({"X", "Y", "Z"});
  connect(axis_, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
          [this](int) { invalidateResult(); });
  step_ = spin(
      .1, 20, .5,
      QStringLiteral("相邻输出截面的纵向间隔。越小采样越密，计算量越大。"));
  thickness_ = spin(
      .2, 30, .5,
      QStringLiteral("每个截面使用的纵向点云厚度。过厚会混合弯曲处的几何。"));
  floorSize_ = spin(
      1, 30, .5,
      QStringLiteral("候选内侧正方形底面取点边长。弯曲明显时可适当缩小。"));
  floorSize_->setObjectName("floorPatchSize");
  floorGuard_ =
      spin(0, 10, .1,
           QStringLiteral("底面取点区域离侧壁的内侧退让量，减少侧壁点混入。"));
  floorGuard_->setObjectName("floorPatchGuard");
  maxCorrection_ =
      spin(0, 20, .5,
           QStringLiteral("底面和侧壁联合交点相对贴壁结果的最大允许位移。"));
  maxCorrection_->setObjectName("floorPatchMaxCorrection");
  form->addRow(QStringLiteral("焊缝延伸轴"), axis_);
  opening_ = new QComboBox(parameters_);
  opening_->setObjectName("interiorOpeningDirection");
  opening_->addItems({QStringLiteral("朝 +V（下凹槽）"), QStringLiteral("朝 -V（反向槽）"), QStringLiteral("不限方向（兼容）")});
  opening_->setToolTip(QStringLiteral("剖面纵轴 V：延伸轴为 X/Y 时是 Z，为 Z 时是 Y。下凹槽只接受底部两侧交界，避免选到反向结构。"));
  connect(opening_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this](int) { invalidateResult(); });
  form->addRow(QStringLiteral("槽口方向"), opening_);
  form->addRow(QStringLiteral("截面步长"), step_);
  form->addRow(QStringLiteral("截面厚度"), thickness_);
  mergeProfileFrames_ = new QCheckBox(QStringLiteral("多帧投影合并拟合"), parameters_);
  mergeProfileFrames_->setObjectName("mergeProfileFrames");
  mergeProfileFrames_->setToolTip(QStringLiteral("把前后相邻原始帧正交投影到当前扫描平面，再分箱并拟合 B 样条。关闭后逐帧独立拟合。仅用于有帧号的输入。"));
  profileMergeFrameCount_ = new QSpinBox(parameters_);
  profileMergeFrameCount_->setObjectName("profileMergeFrameCount");
  profileMergeFrameCount_->setRange(1, 101);
  profileMergeFrameCount_->setSingleStep(2);
  profileMergeFrameCount_->setToolTip(QStringLiteral("输入奇数。5 帧表示当前帧加前后各 2 帧；首尾不足时缩小窗口，缺帧处截断窗口。与调试每步帧数无关。"));
  connect(mergeProfileFrames_, &QCheckBox::toggled, this, [this](bool) {
    updateInputMode();
    invalidateResult();
  });
  connect(profileMergeFrameCount_, QOverload<int>::of(&QSpinBox::valueChanged), this,
          [this](int) { invalidateResult(); });
  form->addRow(mergeProfileFrames_);
  form->addRow(QStringLiteral("合并帧数（奇数）"), profileMergeFrameCount_);
  localizationMethod_ = new QComboBox(parameters_);
  localizationMethod_->setObjectName("profileLocalizationMethod");
  localizationMethod_->addItems({QStringLiteral("斜边曲线 × 底部直线"),
                                 QStringLiteral("连续曲率候选")});
  localizationMethod_->setToolTip(QStringLiteral("斜边曲线模式分别拟合左右斜边的参数化 B 样条，并与槽底鲁棒直线求交；曲率模式沿用连续转折区候选。"));
  profileCurvature_ = new QCheckBox(QStringLiteral("B 样条轮廓拟合"), parameters_);
  profileCurvature_->setObjectName("useProfileCurvature");
  profileCurvature_->setToolTip(QStringLiteral("斜边曲线模式用参数化三次 B 样条拟合左右斜边并与底部直线求交；曲率候选模式拟合完整轮廓并通过解析曲率选点。"));
  curvatureRadius_ = spin(.5, 50, .5, QStringLiteral("样条拟合尺度：节点弧长间距为该值的四分之一。越大越平滑，也可能合并邻近转折。"));
  curvatureRadius_->setObjectName("profileCurvatureRadius");
  centralBeadHeight_ = spin(.05, 20, .05, QStringLiteral("中央凸起相对鲁棒槽底线的最小突出高度。低于该值按波浪底或平底处理。单位与点云一致。"));
  centralBeadHeight_->setObjectName("centralBeadMinHeight");
  centralBeadWidth_ = spin(.1, 50, .1, QStringLiteral("中央凸起在半峰高处的最小横向宽度。窄于该值视为噪声或局部尖点。单位与点云一致。"));
  centralBeadWidth_->setObjectName("centralBeadMinWidth");
  turnLocation_ = new QComboBox(parameters_);
  turnLocation_->setObjectName("profileTurnLocation");
  turnLocation_->addItems({QStringLiteral("单点曲率峰值"), QStringLiteral("连续区间中值"), QStringLiteral("累计转角中点"), QStringLiteral("槽底侧边界"), QStringLiteral("壁脚转角重心（稳定）")});
  turnLocation_->setToolTip(QStringLiteral("壁脚转角重心：汇总同一侧壁支持范围内的曲率转角，在实测样条上选一个代表点，减轻邻近峰值强弱互换引起的跳点。其他选项保留原有峰值/区间定位。所有模式仍受贴壁与底部中点约束。"));
  turnRegionLength_ = spin(.1, 30, .2, QStringLiteral("连续转折区域的最大总弧长，单位与点云一致。峰顶两侧各最多取一半长度；到曲率低谷或变号处提前停止。过大会混合邻近转折。"));
  turnRegionLength_->setObjectName("profileTurnRegionMaxLength");
  turnRegionFraction_ = spin(.05, 1, .05, QStringLiteral("连续转折区间保留的最低相对曲率。0.6 表示从峰顶向两侧扩展到峰值的 60%；越低区间越宽。"));
  turnRegionFraction_->setObjectName("profileTurnRegionFraction");
  connect(turnLocation_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this](int index) {
    turnRegionLength_->setEnabled(profileCurvature_->isChecked() && index != 0 && index != 4);
    turnRegionFraction_->setEnabled(profileCurvature_->isChecked() && index != 0 && index != 4);
    updateAlgorithmControls();
    invalidateResult();
  });
  profileWallTolerance_ = spin(.05, 5, .05, QStringLiteral("先独立识别主侧壁，再限制候选和最终左右轨迹到侧壁的垂直距离。整条曲线拟合也保留这个限制；单道输出取两侧交界中点。"));
  profileWallTolerance_->setObjectName("profileWallTolerance");
  connect(profileCurvature_, &QCheckBox::toggled, this, [this](bool enabled) {
    curvatureRadius_->setEnabled(enabled);
    turnLocation_->setEnabled(enabled);
    turnRegionLength_->setEnabled(enabled && turnLocation_->currentIndex() != 0 && turnLocation_->currentIndex() != 4);
    turnRegionFraction_->setEnabled(enabled && turnLocation_->currentIndex() != 0 && turnLocation_->currentIndex() != 4);
    profileWallTolerance_->setEnabled(enabled);
    invalidateResult();
  });
  connect(localizationMethod_, QOverload<int>::of(&QComboBox::currentIndexChanged),
          this, &MainWindow::switchAlgorithm);
  adaptiveSeamCount_ = new QCheckBox(QStringLiteral("按整条槽宽多数判定道数"), parameters_);
  adaptiveSeamCount_->setObjectName("adaptiveSeamCount");
  adaptiveSeamCount_->setToolTip(QStringLiteral("关闭时固定输出左右两道；开启后按整条焊缝的槽宽多数票决定输出中间单道或左右双道。"));
  doubleSeamWidth_ = spin(.5, 100, .5, QStringLiteral("统计整条焊缝有效站位的左右壁底交界横向间距。超过半数达到阈值，整条输出双道；否则整条输出中间单道，票数相同取单道。单位与点云一致（毫米数据即 mm）。"));
  doubleSeamWidth_->setObjectName("doubleSeamWidth");
  connect(adaptiveSeamCount_, &QCheckBox::toggled, doubleSeamWidth_, &QWidget::setEnabled);
  connect(adaptiveSeamCount_, &QCheckBox::toggled, this, [this] { invalidateResult(); });
  surfaceEnabled_ = new QCheckBox(QStringLiteral("贴面修正（侧壁拟合与投影）"), parameters_);
  surfaceEnabled_->setObjectName("useSurfaceRefinement");
  surfaceEnabled_->setToolTip(QStringLiteral("关闭时跳过局部贴面细化和底面求交；仍保留截面的主侧壁识别、轨迹平滑中的贴壁约束、候选修正上限及左右间距保护。"));
  connect(surfaceEnabled_, &QCheckBox::toggled, this, [this](bool) {
    updateRefinementControls();
    invalidateResult();
  });
  floorEnabled_ = new QCheckBox(QStringLiteral("底面与侧壁联合定位"), parameters_);
  floorEnabled_->setObjectName("useLocalFloorConstraint");
  floorEnabled_->setToolTip(QStringLiteral("贴面修正开启时，可单独关闭底面求交，仅保留侧壁拟合与投影。"));
  connect(floorEnabled_, &QCheckBox::toggled, this, [this](bool) {
    updateRefinementControls();
    invalidateResult();
  });
  finalSmoothing_ = new QCheckBox(QStringLiteral("最终轨迹平滑（抑制帧间跳动）"), parameters_);
  finalSmoothing_->setObjectName("finalSmoothingEnabled");
  finalSmoothing_->setToolTip(QStringLiteral("在贴面和候选回退之后联合拟合左右连续轨迹，抑制前后帧跳动，同时满足主侧壁距离与左右间距限制。关闭可对比未平滑结果。"));
  finalLength_ = spin(.1, 100, 1, QStringLiteral("最终平滑的物理长度。越大越平滑，但会抑制更短的真实弯曲；单位与点云一致。"));
  finalLength_->setObjectName("finalSmoothingLength");
  connect(finalSmoothing_, &QCheckBox::toggled, this, [this](bool enabled) {
    finalLength_->setEnabled(enabled);
    invalidateResult();
  });
  form->addRow(finalSmoothing_);
  form->addRow(QStringLiteral("轨迹平滑长度"), finalLength_);
  uniformWidth_ = new QCheckBox(QStringLiteral("沿槽底统一宽度（两端各去掉 10% 后取均值）"), parameters_);
  uniformWidth_->setObjectName("uniformWidthEnabled");
  uniformWidth_->setToolTip(QStringLiteral("用当前算法初选的有效左右点统计宽度。最终截断后，以每帧测得的槽底中心为基准沿轮廓对称定位，输出左右两道；轮廓不足的帧不外推。此开关和宽度定义由每种算法独立保存。"));
  uniformWidthMetric_ = new QComboBox(parameters_);
  uniformWidthMetric_->setObjectName("uniformWidthMetric");
  uniformWidthMetric_->addItems({QStringLiteral("槽底曲线长度（左右各一半）"),QStringLiteral("左右两点直线距离")});
  uniformWidthFilter_=new QCheckBox(QStringLiteral("统一宽度后整体滤波（保持直线宽度）"),parameters_);
  uniformWidthFilter_->setObjectName("uniformWidthFilterEnabled");
  uniformWidthFilter_->setToolTip(QStringLiteral("沿槽底定位后，联合平滑左右轨迹，保持目标直线宽度与扫描平面，不跨缺口。允许少量偏离原槽底噪声轮廓。仅直线距离模式生效。参数由当前算法独立保存。"));
  uniformWidthFilterLength_=spin(.1,100,1,QStringLiteral("统一宽度后的独立滤波长度，单位与点云一致。越大越平滑，局部细节保留越少。"));
  uniformWidthFilterLength_->setObjectName("uniformWidthFilterLength");
  auto updateWidthFilter=[this] {
    const bool active=uniformWidth_->isChecked() && uniformWidthMetric_->currentIndex()==1;
    uniformWidthFilter_->setEnabled(active);
    uniformWidthFilterLength_->setEnabled(active && uniformWidthFilter_->isChecked());
    invalidateResult();
  };
  connect(uniformWidth_,&QCheckBox::toggled,this,updateWidthFilter);
  connect(uniformWidthMetric_,QOverload<int>::of(&QComboBox::currentIndexChanged),this,updateWidthFilter);
  connect(uniformWidthFilter_,&QCheckBox::toggled,this,updateWidthFilter);
  connect(uniformWidth_, &QCheckBox::toggled, this, [this](bool enabled) {
    uniformWidthMetric_->setEnabled(enabled);
    adaptiveSeamCount_->setEnabled(!enabled);
    doubleSeamWidth_->setEnabled(!enabled && adaptiveSeamCount_->isChecked());
    invalidateResult();
  });
  connect(uniformWidthMetric_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this](int){invalidateResult();});
  form->addRow(uniformWidth_);
  form->addRow(QStringLiteral("统一宽度定义"),uniformWidthMetric_);
  form->addRow(uniformWidthFilter_);
  form->addRow(QStringLiteral("统一宽度后滤波长度"),uniformWidthFilterLength_);
  trimEndpoints_ = new QCheckBox(QStringLiteral("首尾大幅变动处截断"),parameters_);
  trimEndpoints_->setObjectName("trimUnstableEndpoints");
  trimEndpoints_->setToolTip(QStringLiteral("从稳定段向首尾检查方向变化和偏离量；连续异常达到确认长度时，去掉外侧尾段。两侧同步截断，保留中间真实弯曲。"));
  auto *advancedButton =
      new QPushButton(QStringLiteral("展开高级参数"), parameters_);
  advancedButton->setObjectName("advancedParametersButton");
  advancedButton->setCheckable(true);
  form->addRow(advancedButton);
  auto *advanced = new QWidget(parameters_);
  advanced->setObjectName("advancedParameters");
  auto *extra = new QFormLayout(advanced);
  extra->setContentsMargins(0, 0, 0, 0);
  form->insertRow(0, QStringLiteral("当前算法"), localizationMethod_);
  algorithmDescription_ = new QLabel(parameters_);
  algorithmDescription_->setWordWrap(true);
  algorithmDescription_->setStyleSheet("color:#16718a;padding:4px 0;");
  form->insertRow(1, algorithmDescription_);
  algorithmGroup_ = new QGroupBox(parameters_);
  algorithmGroup_->setObjectName("algorithmParameters");
  auto *algorithmForm = new QFormLayout(algorithmGroup_);
  algorithmForm->setRowWrapPolicy(QFormLayout::WrapLongRows);
  minimumFeatureDepth_ = spin(0, 50, .1, QStringLiteral("截面相对板面参考线的最小有效起伏，单位与点云一致。低于此值停止候选计算并断开轨迹，禁止回退和插补跨越。多帧合并时单独检查当前原始帧。0 表示关闭；此值由当前算法独立保存。"));
  minimumFeatureDepth_->setObjectName("minimumFeatureDepth");
  algorithmForm->addRow(QStringLiteral("最小有效槽深（截断）"), minimumFeatureDepth_);
  algorithmForm->addRow(QStringLiteral("B 样条拟合尺度"), curvatureRadius_);
  algorithmForm->addRow(QStringLiteral("中央凸起最小高度"), centralBeadHeight_);
  algorithmForm->addRow(QStringLiteral("中央凸起最小宽度"), centralBeadWidth_);
  algorithmForm->addRow(QStringLiteral("贴壁距离上限"), profileWallTolerance_);
  curvatureGroup_ = new QGroupBox(QStringLiteral("曲率算法 · 候选与连续性"), algorithmGroup_);
  curvatureGroup_->setObjectName("curvatureParameters");
  auto *curvatureForm = new QFormLayout(curvatureGroup_);
  curvatureForm->setRowWrapPolicy(QFormLayout::WrapLongRows);
  curvatureForm->addRow(profileCurvature_);
  curvatureForm->addRow(QStringLiteral("转折代表点"), turnLocation_);
  curvatureForm->addRow(QStringLiteral("转折区间长度"), turnRegionLength_);
  curvatureForm->addRow(QStringLiteral("区间曲率比例"), turnRegionFraction_);
  multiCandidates_ = new QCheckBox(QStringLiteral("加入波谷 / 曲率变化候选"), curvatureGroup_);
  multiCandidates_->setObjectName("multiCandidates");
  multiCandidates_->setToolTip(QStringLiteral("通过侧壁、方向和实测支持检查后保留多个候选，由跨帧路径选择稳定的组合。关闭可对比原曲率候选。"));
  candidateWindow_ = spin(.1,20,.2,QStringLiteral("波谷与曲率变化搜索的总弧长窗口。"));
  valleyProminence_ = spin(.001,5,.01,QStringLiteral("波谷相对窗口两侧的最小高度差。"));
  valleyProminence_->setDecimals(3);
  candidateSpacing_ = spin(.05,5,.05,QStringLiteral("候选配对的去重距离；较小值保留更多邻近定位方案。"));
  continuityWeight_ = spin(0,50,.5,QStringLiteral("切线预测残差的惩罚权重；增大可减少候选在相邻帧之间跳动。"));
  candidatePairs_ = new QSpinBox(curvatureGroup_); candidatePairs_->setRange(2,32);
  curvatureForm->addRow(multiCandidates_);
  curvatureForm->addRow(QStringLiteral("候选搜索窗口"), candidateWindow_);
  curvatureForm->addRow(QStringLiteral("波谷最小深度"), valleyProminence_);
  curvatureForm->addRow(QStringLiteral("候选去重距离"), candidateSpacing_);
  curvatureForm->addRow(QStringLiteral("每帧最多候选对"), candidatePairs_);
  curvatureForm->addRow(QStringLiteral("连续性权重"), continuityWeight_);
  algorithmForm->addRow(curvatureGroup_);
  fallbackGroup_ = new QGroupBox(QStringLiteral("曲线算法 · 独立曲率回退参数"), algorithmGroup_);
  fallbackGroup_->setObjectName("curveLineFallbackParameters");
  auto *fallbackForm = new QFormLayout(fallbackGroup_);
  fallbackEnabled_ = new QCheckBox(QStringLiteral("求交失败时使用曲率回退"), fallbackGroup_);
  fallbackEnabled_->setObjectName("curveLineFallbackEnabled");
  fallbackEnabled_->setToolTip(QStringLiteral("只在本帧无有效曲线 / 底线交点时启用。以下参数属于曲线算法，不读取独立曲率算法的设置。"));
  fallbackForm->addRow(fallbackEnabled_);
  auto *fallbackDetailsButton = new QPushButton(QStringLiteral("展开回退参数"), fallbackGroup_);
  fallbackDetailsButton->setCheckable(true);
  fallbackForm->addRow(fallbackDetailsButton);
  auto *fallbackDetails = new QWidget(fallbackGroup_);
  fallbackForm->addRow(fallbackDetails);
  fallbackForm = new QFormLayout(fallbackDetails);
  fallbackForm->setRowWrapPolicy(QFormLayout::WrapLongRows);
  fallbackForm->setContentsMargins(0,0,0,0);
  fallbackDetails->hide();
  connect(fallbackDetailsButton,&QPushButton::toggled,this,[fallbackDetails,fallbackDetailsButton](bool show) {
    fallbackDetails->setVisible(show);
    fallbackDetailsButton->setText(show ? QStringLiteral("收起回退参数") : QStringLiteral("展开回退参数"));
  });
  fallbackRadius_=spin(.5,50,.5,QStringLiteral("回退分支自己的 B 样条尺度。"));
  fallbackRadius_->setObjectName("fallbackSplineRadius");
  fallbackTolerance_=spin(.05,5,.05,QStringLiteral("回退候选的贴壁距离限制。"));
  fallbackRegionLength_=spin(.1,30,.2,QStringLiteral("回退转折区域的最大弧长。"));
  fallbackRegionFraction_=spin(.05,1,.05,QStringLiteral("回退转折区的最低相对曲率。"));
  fallbackWindow_=spin(.1,20,.2,QStringLiteral("回退波谷与曲率变化搜索的总弧长窗口。"));
  fallbackProminence_=spin(.001,5,.01,QStringLiteral("回退波谷相对窗口两侧的最小高度差。"));
  fallbackProminence_->setDecimals(3);
  fallbackSpacing_=spin(.05,5,.05,QStringLiteral("回退候选配对的去重距离。"));
  fallbackContinuity_=spin(0,50,.5,QStringLiteral("包含回退候选的路径连续性权重。"));
  fallbackPairs_=new QSpinBox(fallbackGroup_); fallbackPairs_->setRange(2,32);
  fallbackForm->addRow(QStringLiteral("回退拟合尺度"),fallbackRadius_);
  fallbackForm->addRow(QStringLiteral("回退贴壁距离"),fallbackTolerance_);
  fallbackForm->addRow(QStringLiteral("回退转折区长度"),fallbackRegionLength_);
  fallbackForm->addRow(QStringLiteral("回退曲率比例"),fallbackRegionFraction_);
  fallbackForm->addRow(QStringLiteral("回退搜索窗口"),fallbackWindow_);
  fallbackForm->addRow(QStringLiteral("回退波谷深度"),fallbackProminence_);
  fallbackForm->addRow(QStringLiteral("回退去重距离"),fallbackSpacing_);
  fallbackForm->addRow(QStringLiteral("回退最多候选对"),fallbackPairs_);
  fallbackForm->addRow(QStringLiteral("回退连续性权重"),fallbackContinuity_);
  algorithmForm->addRow(fallbackGroup_);
  form->insertRow(2,algorithmGroup_);
  connect(multiCandidates_, &QCheckBox::toggled, this, [this] { updateAlgorithmControls(); invalidateResult(); });
  connect(profileCurvature_, &QCheckBox::toggled, this, [this] { updateAlgorithmControls(); });
  connect(fallbackEnabled_, &QCheckBox::toggled, this, [this] { updateAlgorithmControls(); invalidateResult(); });
  for (auto *count : {candidatePairs_,fallbackPairs_})
    connect(count,QOverload<int>::of(&QSpinBox::valueChanged),this,[this] { invalidateResult(); });
  extra->addRow(adaptiveSeamCount_);
  extra->addRow(QStringLiteral("双道宽度阈值"), doubleSeamWidth_);
  extra->addRow(trimEndpoints_);
  endpointAngle_=spin(1,89,1,QStringLiteral("首尾相邻两段轨迹的方向变化阈值（度），越小越早停止。"));
  endpointWindow_=spin(1,100,1,QStringLiteral("首尾检测的物理窗口长度，单位与点云一致。从中间参考段向两端查找持续异常，短轨迹不强行裁剪。"));
  endpointDeviation_=spin(.1,20,.1,QStringLiteral("允许首尾轨迹相对稳定趋势的变动范围；超出后还需满足转角和持续长度条件才会截断。单位与点云一致。"));
  endpointAngle_->setObjectName("endpointTrimAngleDegrees");
  endpointWindow_->setObjectName("endpointTrimWindow");
  endpointDeviation_->setObjectName("endpointTrimDeviation");
  extra->addRow(QStringLiteral("首尾转角阈值"),endpointAngle_);
  extra->addRow(QStringLiteral("首尾检测长度"),endpointWindow_);
  extra->addRow(QStringLiteral("首尾变动范围"),endpointDeviation_);
  connect(trimEndpoints_,&QCheckBox::toggled,this,[this](bool enabled) {
    endpointAngle_->setEnabled(enabled); endpointWindow_->setEnabled(enabled); endpointDeviation_->setEnabled(enabled);
    invalidateResult();
  });
  voxel_ = spin(0, 10, .1, QStringLiteral("预处理体素大小，0 关闭降采样。"));
  bin_ = spin(.1, 10, .1, QStringLiteral("横截面轮廓分箱宽度。"));
  profileMode_ = new QComboBox(advanced);
  profileMode_->addItem(QStringLiteral("原始滚动分箱（兼容）"));
  profileMode_->addItem(QStringLiteral("固定网格与物理尺度（推荐试用）"));
  profileMode_->setToolTip(QStringLiteral(
      "固定网格模式让整云共用分箱原点，并避免分箱宽度同时改变拟合和几何阈值。"));
  profileReferenceScale_ = spin(
      .1, 10, .1,
      QStringLiteral("稳定模式使用的几何参考尺度；通常保持 1.2，不随分箱宽度调整。"));
  connect(profileMode_, QOverload<int>::of(&QComboBox::currentIndexChanged),
          this, [this](int index) {
            profileReferenceScale_->setEnabled(index == 1);
            invalidateResult();
          });
  floorRms_ = spin(
      .02, 5, .05,
      QStringLiteral("局部底面拟合的 RMS 上限。放大阈值会接受更差的拟合。"));
  floorRms_->setObjectName("floorPatchMaxRms");
  separation_ =
      spin(1, 500, 1, QStringLiteral("左右焊缝在横截面中的最小分离距离。"));
  extra->addRow(surfaceEnabled_);
  extra->addRow(floorEnabled_);
  extra->addRow(QStringLiteral("底面取点边长"), floorSize_);
  extra->addRow(QStringLiteral("底面内侧退让"), floorGuard_);
  extra->addRow(QStringLiteral("最大联合修正"), maxCorrection_);
  extra->addRow(QStringLiteral("底面 RMS 上限"), floorRms_);
  extra->addRow(QStringLiteral("体素尺寸"), voxel_);
  extra->addRow(QStringLiteral("剖面分箱"), bin_);
  extra->addRow(QStringLiteral("分箱模式"), profileMode_);
  extra->addRow(QStringLiteral("几何参考尺度"), profileReferenceScale_);
  extra->addRow(QStringLiteral("左右最小间距"), separation_);
  candidateCorrection_ = spin(0, 20, .1, QStringLiteral("最终轨迹平滑之前，相对本帧候选的修正上限。最终平滑使用前后帧趋势，不再逐点回退。"));
  extra->addRow(QStringLiteral("候选修正上限"), candidateCorrection_);
  hullVoxel_ = spin(.1, 10, .1, QStringLiteral("仅用于分析主体，不改变输出原始点坐标。"));
  hullDistance_ = spin(.01, 10, .1, QStringLiteral("参考曲面拟合容差。"));
  hullMargin_ = spin(0, 20, .1, QStringLiteral("底面与交界原始点的清理保护带宽度。"));
  extra->addRow(QStringLiteral("主体分析网格"), hullVoxel_);
  extra->addRow(QStringLiteral("主体曲面容差"), hullDistance_);
  extra->addRow(QStringLiteral("交界保护余量"), hullMargin_);
  debugFlag_ =
      new QCheckBox(QStringLiteral("逐步查看点云（VTK）"), parameters_);
  debugFlag_->setObjectName("debugFlag");
  debugFlag_->setToolTip(
      QStringLiteral("对应 debugFlag。每一步弹出 VTK "
                     "查看器，点击下一步继续；可跳过后续弹窗。"));
  form->addRow(debugFlag_);
  debugMode_ = new QComboBox(parameters_);
  debugMode_->setObjectName("debugMode");
  debugMode_->addItems({QStringLiteral("按帧 / 截面单步"), QStringLiteral("按算法步骤单步")});
  debugBatch_ = new QSpinBox(parameters_);
  debugBatch_->setObjectName("debugFrameBatchSize");
  debugBatch_->setRange(1, 1000);
  debugBatch_->setToolTip(QStringLiteral("每次处理多少帧后暂停；不合并帧参与拟合。无帧号时按空间截面计数。"));
  debugStart_ = new QSpinBox(parameters_);
  debugStart_->setObjectName("debugStartSlice");
  debugStart_->setRange(1, 1);
  debugPositionLabel_ = new QLabel(parameters_);
  debugPositionLabel_->setObjectName("debugSlicePosition");
  debugPositionLabel_->setWordWrap(true);
  debugPositionLabel_->setStyleSheet("color:#16718a;font-weight:600;");
  form->addRow(QStringLiteral("调试方式"), debugMode_);
  form->addRow(QStringLiteral("每步帧数"), debugBatch_);
  form->addRow(QStringLiteral("选择切片 / 帧"), debugStart_);
  form->addRow(QStringLiteral("点云中的位置"), debugPositionLabel_);
  auto updateDebugControls = [this] {
    debugMode_->setEnabled(debugFlag_->isChecked());
    debugBatch_->setEnabled(debugFlag_->isChecked() && debugMode_->currentIndex() == 0);
    debugStart_->setEnabled(debugFlag_->isChecked() && debugMode_->currentIndex() == 0);
  };
  connect(debugFlag_, &QCheckBox::toggled, this, updateDebugControls);
  connect(debugMode_, QOverload<int>::of(&QComboBox::currentIndexChanged), this, updateDebugControls);
  connect(debugStart_, QOverload<int>::of(&QSpinBox::valueChanged), this,
          [this](int) { updateDebugSliceRange(); refreshView(false); });
  connect(debugFlag_, &QCheckBox::toggled, this,
          [this](bool) { refreshView(false); });
  connect(axis_, QOverload<int>::of(&QComboBox::currentIndexChanged), this,
          [this](int) { updateDebugSliceRange(); });
  connect(step_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
          [this](double) { updateDebugSliceRange(); });
  connect(thickness_, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this,
          [this](double) { updateDebugSliceRange(); });
  updateDebugControls();
  advanced->hide();
  form->addRow(advanced);
  connect(advancedButton, &QPushButton::toggled, this,
          [advanced, advancedButton](bool show) {
            advanced->setVisible(show);
            advancedButton->setText(show ? QStringLiteral("收起高级参数")
                                         : QStringLiteral("展开高级参数"));
          });
  auto *unit =
      new QLabel(QStringLiteral("所有长度均使用点云自身单位。"), parameters_);
  unit->setStyleSheet("color:#687787;");
  form->addRow(unit);
  column->addWidget(parameters_);
  defaults_ = new QPushButton(QStringLiteral("仅恢复当前算法默认参数"), panel);
  defaults_->setObjectName("resetAlgorithmDefaults");
  column->addWidget(defaults_);
  connect(defaults_, &QPushButton::clicked, this, &MainWindow::resetDefaults);
  run_ = new QPushButton(QStringLiteral("提取焊缝"), panel);
  run_->setMinimumHeight(40);
  run_->setStyleSheet("background:#146b86;color:white;font-weight:600;");
  run_->setEnabled(false);
  column->addWidget(run_);
  connect(run_, &QPushButton::clicked, this, &MainWindow::runDetection);
  save_ = new QPushButton(QStringLiteral("保存结果…"), panel);
  save_->setEnabled(false);
  column->addWidget(save_);
  connect(save_, &QPushButton::clicked, this, &MainWindow::saveResult);
  auto *displayGroup = new QGroupBox(QStringLiteral("结果显示"), panel);
  auto *view = new QGridLayout(displayGroup);
  auto layer=[&](const QString &title,const char *name,const char *color,bool checked,int row,int col) {
    auto *box=new QCheckBox(title,displayGroup);
    box->setObjectName(name);
    box->setChecked(checked);
    if(color)box->setStyleSheet(QString("color:%1;font-weight:600;").arg(color));
    view->addWidget(box,row,col);
    connect(box,&QCheckBox::toggled,this,[this]{refreshView();});
    return box;
  };
  showLeftSurface_=layer(QStringLiteral("左侧面 · 蓝"),"showLeftSurface","#1677b8",true,0,0);
  showRightSurface_=layer(QStringLiteral("右侧面 · 橙"),"showRightSurface","#ad5b13",true,0,1);
  showWeldSurface_=layer(QStringLiteral("焊道面 · 绿"),"showWeldSurface","#208244",true,1,0);
  showTracks_=layer(QStringLiteral("焊缝轨迹"),"showSeamTracks",nullptr,true,1,1);
  showCloud_=layer(QStringLiteral("叠加原始点云"),"showInputCloud",nullptr,false,2,0);
  showBody_=layer(QStringLiteral("叠加主体点云"),"showMainBody",nullptr,false,2,1);
  auto *reset = new QPushButton(QStringLiteral("重置视角"), displayGroup);
  view->addWidget(reset,3,0,1,2);
  connect(reset, &QPushButton::clicked, this, [this] { refreshView(true); });
  resultLabel_ = new QLabel(QStringLiteral("尚未提取"), panel);
  resultLabel_->setWordWrap(true);
  column->addWidget(resultLabel_);
  log_ = new QPlainTextEdit(panel);
  log_->setReadOnly(true);
  log_->setMaximumBlockCount(300);
  log_->setMinimumHeight(100);
  log_->setMaximumHeight(150);
  column->addWidget(log_);
  column->addStretch();
  auto *scroll = new QScrollArea(central);
  scroll->setWidgetResizable(true);
  scroll->setWidget(panel);
  // Keep extraction, reset and save reachable while either algorithm's
  // parameter panel is expanded and scrolled.
  auto *controls = new QWidget(central);
  controls->setMinimumWidth(430);
  controls->setMaximumWidth(460);
  auto *controlsLayout = new QVBoxLayout(controls);
  controlsLayout->setContentsMargins(0,0,0,0);
  controlsLayout->addWidget(scroll,1);
  controlsLayout->addWidget(displayGroup);
  for (QWidget *action : {static_cast<QWidget *>(defaults_),static_cast<QWidget *>(run_),
                          static_cast<QWidget *>(save_),static_cast<QWidget *>(resultLabel_)})
    column->removeWidget(action);
  controlsLayout->addWidget(defaults_);
  auto *actions = new QHBoxLayout;
  actions->addWidget(run_,1);
  actions->addWidget(save_,1);
  save_->setMinimumHeight(40);
  controlsLayout->addLayout(actions);
  controlsLayout->addWidget(resultLabel_);
  layout->addWidget(controls);
  vtkWidget_ = new QVTKOpenGLNativeWidget(central);
  vtkWidget_->setMinimumSize(600, 500);
  layout->addWidget(vtkWidget_, 1);
  setCentralWidget(central);
}
void MainWindow::initializeViewer() {
  auto renderer = vtkSmartPointer<vtkRenderer>::New();
  auto window = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
  window->AddRenderer(renderer);
  vtkWidget_->SetRenderWindow(window);
  viewer_.reset(new pcl::visualization::PCLVisualizer(renderer, window,
                                                      "scheme1", false));
  viewer_->setupInteractor(vtkWidget_->GetInteractor(),
                           vtkWidget_->GetRenderWindow());
  viewer_->setBackgroundColor(.025, .055, .075);
}
void MainWindow::updateInputMode() {
  const bool framed = !frames_.empty();
  step_->setEnabled(!framed);
  thickness_->setEnabled(!framed);
  mergeProfileFrames_->setEnabled(framed);
  profileMergeFrameCount_->setEnabled(framed && mergeProfileFrames_->isChecked());
  floorEnabled_->setText(QStringLiteral("底面与侧壁联合定位"));
  parameters_->setTitle(framed ? QStringLiteral("原始帧参数（不使用空间切片）")
                               : QStringLiteral("关键参数"));
  updateDebugSliceRange();
}
void MainWindow::updateDebugSliceRange() {
  if (!debugStart_)
    return;
  int total = 1;
  float minimum = 0.0f, maximum = 0.0f;
  bool hasRange = false;
  if (!frames_.empty()) {
    total = std::max(1, int(frames_.size()));
  } else if (cloud_ && !cloud_->empty()) {
    minimum = std::numeric_limits<float>::max();
    maximum = -std::numeric_limits<float>::max();
    const int axis = axis_->currentIndex();
    for (const auto &point : *cloud_) {
      const float value = point.data[axis];
      if (!std::isfinite(value))
        continue;
      minimum = std::min(minimum, value);
      maximum = std::max(maximum, value);
      hasRange = true;
    }
    if (hasRange) {
      const double first = minimum + 0.5 * thickness_->value();
      const double last = maximum - 0.5 * thickness_->value();
      total = last < first ? 1
                           : std::max(1, int(std::floor((last - first) /
                                                       step_->value() + 1e-6)) + 1);
    }
  }
  const int selected = std::min(debugStart_->value(), total);
  const QSignalBlocker blocker(debugStart_);
  debugStart_->setRange(1, total);
  debugStart_->setValue(selected);
  debugStart_->setSuffix(QStringLiteral(" / %1").arg(total));
  if (!frames_.empty()) {
    const int frameId = frames_[std::size_t(selected - 1)].frameIndex;
    const auto &frame = frames_[std::size_t(selected - 1)];
    double station = 0.0;
    std::size_t finite = 0;
    const int axis = axis_->currentIndex();
    for (std::size_t i = frame.begin; i < frame.begin + frame.count && i < cloud_->size(); ++i) {
      const float value = (*cloud_)[i].data[axis];
      if (std::isfinite(value)) { station += value; ++finite; }
    }
    if (finite) station /= double(finite);
    const double percent = total > 1 ? 100.0 * (selected - 1) / (total - 1) : 0.0;
    debugStart_->setToolTip(
        QStringLiteral("选择整个点云的第 %1/%2 帧（原始帧号 %3）。主视图会高亮该帧，执行时首先显示该帧。")
            .arg(selected).arg(total).arg(frameId));
    debugPositionLabel_->setText(
        QStringLiteral("整云 %1% · 原始帧号 %2 · %3≈%4")
            .arg(percent, 0, 'f', 1).arg(frameId)
            .arg(axis == 0 ? "X" : axis == 1 ? "Y" : "Z")
            .arg(station, 0, 'f', 2));
  } else if (hasRange) {
    const double station = minimum + 0.5 * thickness_->value() +
                           (selected - 1) * step_->value();
    const double percent = total > 1 ? 100.0 * (selected - 1) / (total - 1) : 0.0;
    debugStart_->setToolTip(
        QStringLiteral("选择整个点云的第 %1/%2 个空间截面，纵向坐标约为 %3。主视图会高亮该截面。")
            .arg(selected).arg(total).arg(station, 0, 'f', 2));
    debugPositionLabel_->setText(
        QStringLiteral("整云 %1% · %2≈%3")
            .arg(percent, 0, 'f', 1)
            .arg(axis_->currentIndex() == 0 ? "X" : axis_->currentIndex() == 1 ? "Y" : "Z")
            .arg(station, 0, 'f', 2));
  } else {
    debugStart_->setToolTip(QStringLiteral("加载点云后可选择并高亮一个帧或空间截面。"));
    debugPositionLabel_->setText(QStringLiteral("加载点云后显示切片位置"));
  }
}
void MainWindow::updateRefinementControls() {
  const bool enabled = surfaceEnabled_->isChecked();
  floorEnabled_->setEnabled(enabled);
  const bool floor = enabled && floorEnabled_->isChecked();
  floorSize_->setEnabled(floor);
  floorGuard_->setEnabled(floor);
  maxCorrection_->setEnabled(floor);
  floorRms_->setEnabled(floor);
}
void MainWindow::resetDefaults() {
  algorithmParameters_.reset(activeAlgorithm_);
  applyConfig(algorithmParameters_.get(activeAlgorithm_));
}
void MainWindow::switchAlgorithm(int index) {
  if (settingDefaults_ || busy_ || index < 0 || index > 1) return;
  algorithmParameters_.set(activeAlgorithm_, currentConfig());
  activeAlgorithm_ = static_cast<ProfileLocalizationMethod>(index);
  applyConfig(algorithmParameters_.get(activeAlgorithm_));
}
void MainWindow::updateAlgorithmControls() {
  if (!algorithmGroup_ || !fallbackPairs_) return;
  const bool curvature = activeAlgorithm_ == ProfileLocalizationMethod::CurvatureCandidates;
  algorithmGroup_->setTitle(curvature ? QStringLiteral("连续曲率 · 独立参数") : QStringLiteral("斜边曲线 × 底部直线 · 独立参数"));
  algorithmDescription_->setText(curvature
      ? QStringLiteral("默认以壁脚转角重心定位，减少相邻曲率峰之间跳点；可切回多候选定位。下方参数仅属于本算法。")
      : QStringLiteral("优先曲线与底线求交；无有效交点时回退。下方所有参数仅属于本算法，回退参数独立设置。"));
  curvatureGroup_->setVisible(curvature);
  fallbackGroup_->setVisible(!curvature);
  const bool enabled = curvature && profileCurvature_->isChecked();
  profileCurvature_->setEnabled(curvature);
  curvatureRadius_->setEnabled(!curvature || enabled);
  profileWallTolerance_->setEnabled(!curvature || enabled);
  turnLocation_->setEnabled(enabled);
  const bool centroid=turnLocation_->currentIndex()==int(ProfileTurnLocation::WallTurnCentroid);
  turnRegionLength_->setEnabled(enabled && turnLocation_->currentIndex()!=0 && !centroid);
  turnRegionFraction_->setEnabled(enabled && turnLocation_->currentIndex()!=0 && !centroid);
  multiCandidates_->setEnabled(enabled && !centroid);
  for (auto *s : {candidateWindow_,valleyProminence_,candidateSpacing_})
    s->setEnabled(enabled && !centroid && multiCandidates_->isChecked());
  continuityWeight_->setEnabled(enabled && multiCandidates_->isChecked());
  candidatePairs_->setEnabled(enabled && !centroid && multiCandidates_->isChecked());
  for (auto *s : {fallbackRadius_,fallbackTolerance_,fallbackRegionLength_,fallbackRegionFraction_,
                  fallbackWindow_,fallbackProminence_,fallbackSpacing_,fallbackContinuity_})
    s->setEnabled(fallbackEnabled_->isChecked());
  fallbackPairs_->setEnabled(fallbackEnabled_->isChecked());
}
void MainWindow::applyConfig(const DetectorConfig &c) {
  settingDefaults_ = true;
  axis_->setCurrentIndex(int(c.longitudinalAxis));
  step_->setValue(c.sliceStep);
  thickness_->setValue(c.sliceThickness);
  mergeProfileFrames_->setChecked(c.mergeProfileFrames);
  profileMergeFrameCount_->setValue(c.profileMergeFrameCount);
  floorSize_->setValue(c.floorPatchSize);
  floorGuard_->setValue(c.floorPatchGuard);
  maxCorrection_->setValue(c.floorPatchMaxCorrection);
  voxel_->setValue(c.voxelSize);
  bin_->setValue(c.profileBinSize);
  profileMode_->setCurrentIndex(
      c.profileSamplingMode == ProfileSamplingMode::StableGrid ? 1 : 0);
  profileReferenceScale_->setValue(c.profileReferenceScale);
  profileReferenceScale_->setEnabled(profileMode_->currentIndex() == 1);
  voxel_->setEnabled(false);
  floorRms_->setValue(c.floorPatchMaxRms);
  separation_->setValue(c.minSeamSeparation);
  opening_->setCurrentIndex(c.interiorOpeningDirection == 1 ? 0 : c.interiorOpeningDirection == -1 ? 1 : 2);
  candidateCorrection_->setValue(c.maxCandidateCorrection);
  floorEnabled_->setChecked(c.useLocalFloorConstraint);
  surfaceEnabled_->setChecked(c.useSurfaceRefinement);
  trimEndpoints_->setChecked(c.trimUnstableEndpoints);
  endpointAngle_->setValue(c.endpointTrimAngleDegrees); endpointWindow_->setValue(c.endpointTrimWindow);
  endpointDeviation_->setValue(c.endpointTrimDeviation);
  endpointAngle_->setEnabled(c.trimUnstableEndpoints); endpointWindow_->setEnabled(c.trimUnstableEndpoints); endpointDeviation_->setEnabled(c.trimUnstableEndpoints);
  finalSmoothing_->setChecked(c.finalSmoothingEnabled);
  uniformWidth_->setChecked(c.uniformWidthEnabled);
  uniformWidthMetric_->setCurrentIndex(int(c.uniformWidthMetric));
  uniformWidthMetric_->setEnabled(c.uniformWidthEnabled);
  uniformWidthFilter_->setChecked(c.uniformWidthFilterEnabled);
  uniformWidthFilterLength_->setValue(c.uniformWidthFilterLength);
  uniformWidthFilter_->setEnabled(c.uniformWidthEnabled && c.uniformWidthMetric==UniformWidthMetric::Chord);
  uniformWidthFilterLength_->setEnabled(uniformWidthFilter_->isEnabled() && c.uniformWidthFilterEnabled);
  adaptiveSeamCount_->setEnabled(!c.uniformWidthEnabled);
  doubleSeamWidth_->setEnabled(!c.uniformWidthEnabled && c.adaptiveSeamCount);
  finalLength_->setValue(c.finalSmoothingLength);
  finalLength_->setEnabled(c.finalSmoothingEnabled);
  adaptiveSeamCount_->setChecked(c.adaptiveSeamCount);
  doubleSeamWidth_->setValue(c.doubleSeamWidth);
  doubleSeamWidth_->setEnabled(c.adaptiveSeamCount && !c.uniformWidthEnabled);
  profileCurvature_->setChecked(c.useProfileCurvature);
  localizationMethod_->setCurrentIndex(int(c.profileLocalizationMethod));
  const bool curvatureCandidates =
      c.profileLocalizationMethod == ProfileLocalizationMethod::CurvatureCandidates;
  profileCurvature_->setEnabled(curvatureCandidates);
  curvatureRadius_->setValue(c.profileCurvatureRadius);
  centralBeadHeight_->setValue(c.centralBeadMinHeight);
  minimumFeatureDepth_->setValue(c.minimumFeatureDepth);
  centralBeadWidth_->setValue(c.centralBeadMinWidth);
  centralBeadHeight_->setEnabled(true);
  centralBeadWidth_->setEnabled(true);
  turnLocation_->setCurrentIndex(int(c.profileTurnLocation));
  turnLocation_->setEnabled(c.useProfileCurvature && curvatureCandidates);
  turnRegionLength_->setValue(c.profileTurnRegionMaxLength);
  turnRegionLength_->setEnabled(c.useProfileCurvature && curvatureCandidates && c.profileTurnLocation != ProfileTurnLocation::Peak && c.profileTurnLocation != ProfileTurnLocation::WallTurnCentroid);
  turnRegionFraction_->setValue(c.profileTurnRegionFraction);
  turnRegionFraction_->setEnabled(c.useProfileCurvature && curvatureCandidates && c.profileTurnLocation != ProfileTurnLocation::Peak && c.profileTurnLocation != ProfileTurnLocation::WallTurnCentroid);
  curvatureRadius_->setEnabled(c.useProfileCurvature);
  profileWallTolerance_->setValue(c.profileWallTolerance);
  profileWallTolerance_->setEnabled(c.useProfileCurvature);
  multiCandidates_->setChecked(c.curvatureSearch.multiCandidates);
  candidateWindow_->setValue(c.curvatureSearch.windowLength);
  valleyProminence_->setValue(c.curvatureSearch.valleyProminence);
  candidateSpacing_->setValue(c.curvatureSearch.candidateSpacing);
  candidatePairs_->setValue(c.curvatureSearch.maximumPairs);
  continuityWeight_->setValue(c.curvatureSearch.continuityWeight);
  fallbackEnabled_->setChecked(c.curveLineFallback.enabled);
  fallbackRadius_->setValue(c.curveLineFallback.splineRadius);
  fallbackTolerance_->setValue(c.curveLineFallback.wallTolerance);
  fallbackRegionLength_->setValue(c.curveLineFallback.regionLength);
  fallbackRegionFraction_->setValue(c.curveLineFallback.regionFraction);
  fallbackWindow_->setValue(c.curveLineFallback.search.windowLength);
  fallbackProminence_->setValue(c.curveLineFallback.search.valleyProminence);
  fallbackSpacing_->setValue(c.curveLineFallback.search.candidateSpacing);
  fallbackPairs_->setValue(c.curveLineFallback.search.maximumPairs);
  fallbackContinuity_->setValue(c.curveLineFallback.search.continuityWeight);
  updateAlgorithmControls();
  updateRefinementControls();
  debugFlag_->setChecked(c.debugFlag);
  hullVoxel_->setValue(c.hullSurfaceVoxelSize);
  hullDistance_->setValue(c.hullSurfaceDistance);
  hullMargin_->setValue(c.hullPreservationMargin);
  debugMode_->setCurrentIndex(c.debugFrameBatchSize == 0 ? 1 : 0);
  debugBatch_->setValue(std::max(1,c.debugFrameBatchSize));
  debugStart_->setValue(c.debugStartSlice);
  updateInputMode();
  settingDefaults_ = false;
  invalidateResult();
}
DetectorConfig MainWindow::currentConfig() const {
  DetectorConfig c = algorithmParameters_.get(activeAlgorithm_);
  c.longitudinalAxis = Axis(axis_->currentIndex());
  c.sliceStep = float(step_->value());
  c.sliceThickness = float(thickness_->value());
  c.mergeProfileFrames = mergeProfileFrames_->isChecked();
  c.profileMergeFrameCount = profileMergeFrameCount_->value();
  c.floorPatchSize = float(floorSize_->value());
  c.floorPatchGuard = float(floorGuard_->value());
  c.floorPatchMaxCorrection = float(maxCorrection_->value());
  c.voxelSize = float(voxel_->value());
  c.profileBinSize = float(bin_->value());
  c.profileSamplingMode = profileMode_->currentIndex() == 1
                              ? ProfileSamplingMode::StableGrid
                              : ProfileSamplingMode::LegacyRolling;
  c.profileReferenceScale = float(profileReferenceScale_->value());
  c.floorPatchMaxRms = float(floorRms_->value());
  c.minSeamSeparation = float(separation_->value());
  c.interiorOpeningDirection = opening_->currentIndex() == 0 ? 1 : opening_->currentIndex() == 1 ? -1 : 0;
  c.maxCandidateCorrection = float(candidateCorrection_->value());
  c.useLocalFloorConstraint = floorEnabled_->isChecked();
  c.useSurfaceRefinement = surfaceEnabled_->isChecked();
  c.trimUnstableEndpoints=trimEndpoints_->isChecked();
  c.endpointTrimAngleDegrees=float(endpointAngle_->value()); c.endpointTrimWindow=float(endpointWindow_->value());
  c.endpointTrimDeviation=float(endpointDeviation_->value());
  c.finalSmoothingEnabled = finalSmoothing_->isChecked();
  c.uniformWidthEnabled=uniformWidth_->isChecked();
  c.uniformWidthMetric=static_cast<UniformWidthMetric>(uniformWidthMetric_->currentIndex());
  c.uniformWidthFilterEnabled=uniformWidthFilter_->isChecked();
  c.uniformWidthFilterLength=float(uniformWidthFilterLength_->value());
  c.finalSmoothingLength = float(finalLength_->value());
  c.adaptiveSeamCount = adaptiveSeamCount_->isChecked();
  c.doubleSeamWidth = float(doubleSeamWidth_->value());
  c.profileLocalizationMethod = activeAlgorithm_;
  c.useProfileCurvature = c.profileLocalizationMethod == ProfileLocalizationMethod::CurveLineIntersection ||
                          profileCurvature_->isChecked();
  c.profileCurvatureRadius = float(curvatureRadius_->value());
  c.centralBeadMinHeight = float(centralBeadHeight_->value());
  c.minimumFeatureDepth = float(minimumFeatureDepth_->value());
  c.centralBeadMinWidth = float(centralBeadWidth_->value());
  c.profileTurnLocation = static_cast<ProfileTurnLocation>(turnLocation_->currentIndex());
  c.profileTurnRegionMaxLength = float(turnRegionLength_->value());
  c.profileTurnRegionFraction = float(turnRegionFraction_->value());
  c.profileWallTolerance = float(profileWallTolerance_->value());
  c.curvatureSearch.multiCandidates = multiCandidates_->isChecked();
  c.curvatureSearch.windowLength = float(candidateWindow_->value());
  c.curvatureSearch.valleyProminence = float(valleyProminence_->value());
  c.curvatureSearch.candidateSpacing = float(candidateSpacing_->value());
  c.curvatureSearch.maximumPairs = candidatePairs_->value();
  c.curvatureSearch.continuityWeight = float(continuityWeight_->value());
  c.curveLineFallback.enabled = fallbackEnabled_->isChecked();
  c.curveLineFallback.splineRadius = float(fallbackRadius_->value());
  c.curveLineFallback.wallTolerance = float(fallbackTolerance_->value());
  c.curveLineFallback.regionLength = float(fallbackRegionLength_->value());
  c.curveLineFallback.regionFraction = float(fallbackRegionFraction_->value());
  c.curveLineFallback.search.windowLength = float(fallbackWindow_->value());
  c.curveLineFallback.search.valleyProminence = float(fallbackProminence_->value());
  c.curveLineFallback.search.candidateSpacing = float(fallbackSpacing_->value());
  c.curveLineFallback.search.maximumPairs = fallbackPairs_->value();
  c.curveLineFallback.search.continuityWeight = float(fallbackContinuity_->value());
  c.debugFlag = debugFlag_->isChecked();
  // Preserve extractMainBody from the selected local algorithm configuration.
  c.hullSurfaceVoxelSize = float(hullVoxel_->value());
  c.hullSurfaceDistance = float(hullDistance_->value());
  c.hullPreservationMargin = float(hullMargin_->value());
  c.debugFrameBatchSize = debugMode_->currentIndex() == 0 ? debugBatch_->value() : 0;
  c.debugStartSlice = debugStart_->value();
  if (!frames_.empty()) {
    const std::size_t selected = std::min<std::size_t>(
        std::size_t(std::max(1, debugStart_->value()) - 1), frames_.size() - 1);
    c.debugSelectedFrameIndex = frames_[selected].frameIndex;
  }
  return c;
}
void MainWindow::invalidateResult() {
  if (settingDefaults_ || busy_)
    return;
  result_ = SeamResult();
  if (save_)
    save_->setEnabled(false);
  if (resultLabel_)
    resultLabel_->setText(cloud_->empty()
                              ? QStringLiteral("尚未提取")
                              : QStringLiteral("参数已更新，请重新提取。"));
  refreshView();
}
bool MainWindow::loadPath(const QString &path) {
  if (busy_)
    return false;
  std::string error;
  std::vector<ScanFrameRange> frames;
  auto input =
      PointCloudIO::load(path.toLocal8Bit().constData(), &error, &frames);
  if (!input || input->empty()) {
    appendLog(QString::fromLocal8Bit(error.c_str()));
    return false;
  }
  cloud_ = input;
  frames_ = std::move(frames);
  updateInputMode();
  inputPath_ = QFileInfo(path).absoluteFilePath();
  inputLabel_->setText(QStringLiteral("%1\n%2 个点")
                           .arg(QFileInfo(path).fileName())
                           .arg(cloud_->size()));
  if (!frames_.empty())
    inputLabel_->setText(
        inputLabel_->text() +
        QStringLiteral("\n%1 个原始轮廓帧（不重新切片）").arg(frames_.size()));
  inputLabel_->setToolTip(inputPath_);
  invalidateResult();
  run_->setEnabled(true);
  appendLog(
      QStringLiteral("加载 %1，点数 %2").arg(inputPath_).arg(cloud_->size()));
  refreshView(true);
  return true;
}
void MainWindow::runDetection() {
  if (busy_ || cloud_->empty())
    return;
  resultConfig_ = currentConfig();
  result_ = SeamResult();
  busy_ = true;
  parameters_->setEnabled(false);
  defaults_->setEnabled(false);
  load_->setEnabled(false);
  run_->setEnabled(false);
  save_->setEnabled(false);
  resultLabel_->setText(QStringLiteral("正在提取…"));
  if (frames_.empty())
    appendLog(QStringLiteral("开始方案一：步长 %1，厚度 %2，底面边长 %3")
                  .arg(resultConfig_.sliceStep)
                  .arg(resultConfig_.sliceThickness)
                  .arg(resultConfig_.floorPatchSize));
  // 工作线程只接触输入快照和算法；所有 VTK/Qt 控件更新都在 GUI 线程完成。
  const auto input = cloud_;
  const auto config = resultConfig_;
  const auto frames = frames_;
  if (!frames.empty())
    appendLog(QStringLiteral(
                  "按 %1 个原始帧提取：轮廓投影合并 %2 帧，保留当前帧号与扫描平面。")
                  .arg(frames.size())
                  .arg(config.mergeProfileFrames ? config.profileMergeFrameCount : 1));
  if (config.debugFlag) {
    // VTK 和模态步骤窗口只在主线程创建；调试模式逐步等待用户继续。
    QTimer::singleShot(0, this, [this, input, config, frames] {
      SeamResult result;
      try {
        result = detectWithDebugViewer(input, config, this, frames);
      } catch (const std::exception &e) {
        result.summary = e.what();
      } catch (...) {
        result.summary = "Unexpected algorithm failure";
      }
      finishDetection(result);
    });
    return;
  }
  watcher_.setFuture(QtConcurrent::run([input, config, frames] {
    try {
      return SeamDetector().detect(input, config, nullptr, nullptr, nullptr,
                                   &frames);
    } catch (const std::exception &e) {
      SeamResult r;
      r.summary = e.what();
      return r;
    } catch (...) {
      SeamResult r;
      r.summary = "Unexpected algorithm failure";
      return r;
    }
  }));
}
void MainWindow::saveResult() {
  if (!result_.ok())
    return;
  QString selectedFilter;
  auto path = QFileDialog::getSaveFileName(
      this, QStringLiteral("保存轨迹、三面点云及 B 样条曲面模型"),
      QFileInfo(inputPath_).absolutePath() + "/" +
          QFileInfo(inputPath_).completeBaseName() + "_seams.pcd",
      QStringLiteral("PCD 点云（轨迹和面分别保存） (*.pcd);;CSV 表格 (*.csv)"), &selectedFilter);
  if (path.isEmpty())
    return;
  if (QFileInfo(path).suffix().isEmpty())
    path += selectedFilter.startsWith("CSV") ? ".csv" : ".pcd";
  QString error;
  if (!saveResultToPath(path, &error)) {
    QMessageBox::critical(this, QStringLiteral("保存失败"), error);
    return;
  }
}
bool MainWindow::saveResultToPath(const QString &path, QString *error) {
  if (busy_ || !result_.ok()) {
    if (error) *error=QStringLiteral("没有可保存的完整提取结果。");
    return false;
  }
  if (!saveResults(path, result_, resultConfig_, error)) return false;
  appendLog(QStringLiteral("已保存到：%1；%2；三个面为 *_left_surface.pcd、*_right_surface.pcd、*_weld_surface.pcd，曲面函数为 *_left_surface_model.json、*_weld_surface_model.json、*_right_surface_model.json，另附 CSV 和参数 JSON。")
                .arg(QFileInfo(path).absolutePath())
                .arg(result_.centerSeam->empty()
                         ? QStringLiteral("左焊缝为 *_left.pcd，右焊缝为 *_right.pcd")
                         : QStringLiteral("单道焊缝为 *_center.pcd")));
  return true;
}
void MainWindow::refreshView(bool reset) {
  if (!viewer_)
    return;
  viewer_->removeAllPointClouds();
  viewer_->removeAllShapes();
  if (cloud_ && !cloud_->empty() &&
      (showCloud_->isChecked() || !result_.ok())) {
    pcl::visualization::PointCloudColorHandlerCustom<PointT> color(cloud_, 160,160,160);
    viewer_->addPointCloud<PointT>(cloud_, color, "input");
    viewer_->setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "input");
    viewer_->setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_OPACITY, result_.ok()?.2:1., "input");
  }
  if (debugFlag_->isChecked() && cloud_ && !cloud_->empty() && debugStart_) {
    CloudT::Ptr selected(new CloudT);
    if (!frames_.empty()) {
      const std::size_t index = std::min<std::size_t>(
          std::size_t(std::max(1, debugStart_->value()) - 1), frames_.size() - 1);
      const auto &frame = frames_[index];
      selected->reserve(frame.count);
      for (std::size_t i = frame.begin; i < frame.begin + frame.count && i < cloud_->size(); ++i)
        selected->push_back((*cloud_)[i]);
    } else {
      float minimum = std::numeric_limits<float>::max();
      const int axis = axis_->currentIndex();
      for (const auto &point : *cloud_)
        if (std::isfinite(point.data[axis])) minimum = std::min(minimum, point.data[axis]);
      const float station = minimum + float(0.5 * thickness_->value() +
                            (debugStart_->value() - 1) * step_->value());
      const float half = float(0.5 * thickness_->value());
      for (const auto &point : *cloud_)
        if (std::isfinite(point.data[axis]) && std::fabs(point.data[axis] - station) <= half)
          selected->push_back(point);
    }
    if (!selected->empty()) {
      pcl::visualization::PointCloudColorHandlerCustom<PointT> selectedColor(selected, 255, 220, 35);
      viewer_->addPointCloud<PointT>(selected, selectedColor, "selected_slice_preview");
      viewer_->setPointCloudRenderingProperties(
          pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 6, "selected_slice_preview");
    }
  }
  if (result_.ok()) {
    auto surface=[&](const CloudT::Ptr &cloud,QCheckBox *toggle,const char *id,int r,int g,int b) {
      if(!toggle->isChecked() || !cloud || cloud->empty())return;
      pcl::visualization::PointCloudColorHandlerCustom<PointT> color(cloud,r,g,b);
      viewer_->addPointCloud<PointT>(cloud,color,id);
      viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE,2,id);
    };
    surface(result_.leftSurface,showLeftSurface_,"left_surface",60,155,245);
    surface(result_.rightSurface,showRightSurface_,"right_surface",245,145,45);
    surface(result_.weldSurface,showWeldSurface_,"weld_surface",55,195,105);
  }
  if (result_.ok() && showTracks_->isChecked())
    for (int side = 0; side < 3; ++side) {
      auto seam = side==2 ? result_.centerSeam : side ? result_.rightSeam : result_.leftSeam;
      if (seam->empty()) continue;
      const auto &starts =
          side==2 ? result_.centerSegmentStarts : side ? result_.rightSegmentStarts : result_.leftSegmentStarts;
      const std::string id = side==2 ? "center" : side ? "right" : "left";
      pcl::visualization::PointCloudColorHandlerCustom<PointT> color(
          seam, side==2 ? 80 : 255, side==1 ? 90 : 235, side==2 ? 255 : side==1 ? 235 : 30);
      viewer_->addPointCloud<PointT>(seam, color, id);
      viewer_->setPointCloudRenderingProperties(
          pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, id);
      for (std::size_t i = 1; i < seam->size(); ++i) {
        if (std::find(starts.begin(), starts.end(), i) != starts.end())
          continue;
        viewer_->addLine(seam->points[i - 1], seam->points[i], side==2 ? .31 : 1.,
                         side==1 ? .35 : .92, side==2 ? 1. : side==1 ? .92 : .12,
                         id + std::to_string(i));
      }
    }
  if (showBody_->isChecked() && result_.mainBodyCloud && !result_.mainBodyCloud->empty()) {
    pcl::visualization::PointCloudColorHandlerCustom<PointT> bodyColor(result_.mainBodyCloud, 185,185,185);
    viewer_->addPointCloud<PointT>(result_.mainBodyCloud, bodyColor, "main_body");
    viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, "main_body");
    viewer_->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_OPACITY, .2, "main_body");
  }
  if (reset)
    viewer_->resetCamera();
  vtkWidget_->update();
}
void MainWindow::appendLog(const QString &message) {
  log_->appendPlainText(QDateTime::currentDateTime().toString("HH:mm:ss ") +
                        message);
}
} // namespace weld
