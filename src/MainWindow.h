#pragma once
#include "SeamTypes.h"
#include "AlgorithmParameters.h"
#include <QElapsedTimer>
#include <QFutureWatcher>
#include <QMainWindow>
#include <pcl/visualization/pcl_visualizer.h>
class QVTKOpenGLNativeWidget;
class QComboBox;
class QVBoxLayout;
class QDoubleSpinBox;
class QSpinBox;
class QCheckBox;
class QPushButton;
class QLabel;
class QPlainTextEdit;
class QGroupBox;
namespace weld {
class MainWindow : public QMainWindow {
  Q_OBJECT
public:
  explicit MainWindow(QWidget *parent = nullptr);
  ~MainWindow() override;
  bool loadPath(const QString &path);
  void runDetection();
  bool saveResultToPath(const QString &path, QString *error = nullptr);
  void setDebugFlag(bool enabled);
signals:
  void extractionFinished(bool ok);

private:
  struct DetectionRun {
    SeamResult result;
    double algorithmMilliseconds = 0;
  };
  void finishDetection(const DetectionRun &run);
  void updateInputMode();
  void updateDebugSliceRange();
  void updateRefinementControls();
  void buildUi();
  void initializeViewer();
  void refreshView(bool reset = false);
  void invalidateResult();
  void resetDefaults();
  void applyConfig(const DetectorConfig &config);
  void switchAlgorithm(int index);
  void updateAlgorithmControls();
  void saveResult();
  void appendLog(const QString &message);
  DetectorConfig currentConfig() const;
  QVTKOpenGLNativeWidget *vtkWidget_ = nullptr;
  QGroupBox *parameters_ = nullptr;
  QComboBox *axis_ = nullptr, *profileMode_ = nullptr;
  QDoubleSpinBox *step_ = nullptr, *thickness_ = nullptr, *floorSize_ = nullptr,
                 *floorGuard_ = nullptr, *maxCorrection_ = nullptr;
  QDoubleSpinBox *inputVoxel_ = nullptr, *voxel_ = nullptr, *bin_ = nullptr,
                 *profileReferenceScale_ = nullptr, *floorRms_ = nullptr,
                 *separation_ = nullptr;
  QCheckBox *surfaceEnabled_ = nullptr, *floorEnabled_ = nullptr,
            *showCloud_ = nullptr;
  QCheckBox *showLeftSurface_ = nullptr, *showRightSurface_ = nullptr,
            *showWeldSurface_ = nullptr, *showTracks_ = nullptr,
            *showBody_ = nullptr;
  QCheckBox *trimEndpoints_ = nullptr;
  QDoubleSpinBox *minimumFeatureDepth_ = nullptr;
  QDoubleSpinBox *endpointAngle_ = nullptr, *endpointWindow_ = nullptr, *endpointDeviation_ = nullptr;
  QCheckBox *finalSmoothing_ = nullptr;
  QCheckBox *uniformWidth_ = nullptr;
  QComboBox *uniformWidthMetric_ = nullptr;
  QCheckBox *uniformWidthFilter_ = nullptr;
  QDoubleSpinBox *uniformWidthFilterLength_ = nullptr;
  QDoubleSpinBox *finalLength_ = nullptr;
  QCheckBox *adaptiveSeamCount_ = nullptr;
  QDoubleSpinBox *doubleSeamWidth_ = nullptr;
  QCheckBox *profileCurvature_ = nullptr;
  QComboBox *localizationMethod_ = nullptr;
  QGroupBox *algorithmGroup_ = nullptr, *curvatureGroup_ = nullptr, *fallbackGroup_ = nullptr;
  QLabel *algorithmDescription_ = nullptr;
  QCheckBox *multiCandidates_ = nullptr, *fallbackEnabled_ = nullptr;
  QDoubleSpinBox *candidateWindow_ = nullptr, *valleyProminence_ = nullptr,
                 *candidateSpacing_ = nullptr, *continuityWeight_ = nullptr;
  QSpinBox *candidatePairs_ = nullptr;
  QDoubleSpinBox *fallbackRadius_ = nullptr, *fallbackTolerance_ = nullptr,
                 *fallbackRegionLength_ = nullptr, *fallbackRegionFraction_ = nullptr,
                 *fallbackWindow_ = nullptr, *fallbackProminence_ = nullptr,
                 *fallbackSpacing_ = nullptr, *fallbackContinuity_ = nullptr;
  QSpinBox *fallbackPairs_ = nullptr;
  QCheckBox *mergeProfileFrames_ = nullptr;
  QSpinBox *profileMergeFrameCount_ = nullptr;
  QDoubleSpinBox *curvatureRadius_ = nullptr;
  QDoubleSpinBox *centralBeadHeight_ = nullptr, *centralBeadWidth_ = nullptr;
  QComboBox *turnLocation_ = nullptr;
  QDoubleSpinBox *turnRegionLength_ = nullptr, *turnRegionFraction_ = nullptr;
  QDoubleSpinBox *profileWallTolerance_ = nullptr;
  QCheckBox *debugFlag_ = nullptr;
  QComboBox *debugMode_ = nullptr;
  QComboBox *opening_ = nullptr;
  QDoubleSpinBox *candidateCorrection_ = nullptr;
  QSpinBox *debugBatch_ = nullptr, *debugStart_ = nullptr;
  QDoubleSpinBox *hullVoxel_ = nullptr, *hullDistance_ = nullptr,
                 *hullMargin_ = nullptr;
  QPushButton *load_ = nullptr, *run_ = nullptr, *save_ = nullptr,
              *defaults_ = nullptr;
  QLabel *inputLabel_ = nullptr, *resultLabel_ = nullptr,
         *debugPositionLabel_ = nullptr;
  QPlainTextEdit *log_ = nullptr;
  pcl::visualization::PCLVisualizer::Ptr viewer_;
  std::vector<ScanFrameRange> frames_;
  CloudT::Ptr cloud_{new CloudT};
  SeamResult result_;
  DetectorConfig resultConfig_;
  AlgorithmParameters algorithmParameters_;
  ProfileLocalizationMethod activeAlgorithm_ = ProfileLocalizationMethod::CurvatureCandidates;
  QString inputPath_;
  QFutureWatcher<DetectionRun> watcher_;
  QElapsedTimer runTimer_;
  bool busy_ = false, settingDefaults_ = false;
};
} // namespace weld
