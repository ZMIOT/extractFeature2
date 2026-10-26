#include "MainWindow.h"
#include <QDialog>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QPushButton>
#include <QPlainTextEdit>
#include <QSpinBox>
#include <QScrollArea>
#include <QScrollBar>
#include <QApplication>
#include <QSurfaceFormat>
#include <QTimer>
#include <QVTKOpenGLNativeWidget.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkActorCollection.h>
#include <cstdio>
#include <cstring>

int main(int argc, char **argv) {
  bool smoke = false;
  for (int i=1;i<argc;++i) if(std::strcmp(argv[i],"--smoke-test")==0) smoke=true;
  auto stage=[smoke](const char *message){if(smoke){std::fprintf(stderr,"%s\n",message);std::fflush(stderr);}};
  stage("Before QApplication");
  QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
  QApplication app(argc, argv);
  stage("Before MainWindow");
  weld::MainWindow window;
  if (smoke) {
    auto *advanced = window.findChild<QWidget *>("advancedParameters");
    auto *surface = window.findChild<QCheckBox *>("useSurfaceRefinement");
    auto *floor = window.findChild<QCheckBox *>("useLocalFloorConstraint");
    auto *floorSize = window.findChild<QDoubleSpinBox *>("floorPatchSize");
    auto *floorGuard = window.findChild<QDoubleSpinBox *>("floorPatchGuard");
    auto *maxCorrection = window.findChild<QDoubleSpinBox *>("floorPatchMaxCorrection");
    auto *floorRms = window.findChild<QDoubleSpinBox *>("floorPatchMaxRms");
    auto *turn = window.findChild<QComboBox *>("profileTurnLocation");
    auto *localization = window.findChild<QComboBox *>("profileLocalizationMethod");
    auto *spline = window.findChild<QCheckBox *>("useProfileCurvature");
    auto *trim = window.findChild<QCheckBox *>("trimUnstableEndpoints");
    auto *radius = window.findChild<QDoubleSpinBox *>("profileCurvatureRadius");
    auto *beadHeight = window.findChild<QDoubleSpinBox *>("centralBeadMinHeight");
    auto *beadWidth = window.findChild<QDoubleSpinBox *>("centralBeadMinWidth");
    auto *regionLength = window.findChild<QDoubleSpinBox *>("profileTurnRegionMaxLength");
    auto *regionFraction = window.findChild<QDoubleSpinBox *>("profileTurnRegionFraction");
    auto *wallTolerance = window.findChild<QDoubleSpinBox *>("profileWallTolerance");
    auto *widthRule = window.findChild<QCheckBox *>("adaptiveSeamCount");
    auto *surfaceModels = window.findChild<QCheckBox *>("generateSurfaceModels");
    auto *doubleWidth = window.findChild<QDoubleSpinBox *>("doubleSeamWidth");
    auto *debugStart = window.findChild<QSpinBox *>("debugStartSlice");
    auto *debugPosition = window.findChild<QLabel *>("debugSlicePosition");
    const bool defaultsOk = advanced && surface && floor && floorSize && floorGuard &&
        maxCorrection && floorRms && turn && localization && spline && trim &&
        radius && beadHeight && beadWidth && regionLength && regionFraction && wallTolerance && widthRule && doubleWidth &&
        debugStart && debugPosition && surfaceModels && surfaceModels->isChecked() &&
        !surface->isChecked() && !widthRule->isChecked() &&
        turn->count() == 5 && turn->currentIndex() == 4 &&
        localization->count() == 2 && localization->currentIndex() == 1 &&
        std::fabs(beadHeight->value() - .8) < 1.e-6 &&
        std::fabs(beadWidth->value() - 2.) < 1.e-6 &&
        beadHeight->isEnabled() && beadWidth->isEnabled() &&
        debugStart->value() == 1 &&
        advanced->isAncestorOf(surface) && advanced->isAncestorOf(floor) &&
        advanced->isAncestorOf(floorSize) && advanced->isAncestorOf(floorGuard) &&
        advanced->isAncestorOf(maxCorrection) && advanced->isAncestorOf(floorRms) &&
        !advanced->isAncestorOf(localization) && !advanced->isAncestorOf(radius) &&
        advanced->isAncestorOf(widthRule) && advanced->isAncestorOf(doubleWidth) &&
        advanced->isAncestorOf(trim) &&
        !window.findChild<QCheckBox *>("extractMainBody");
    if (!defaultsOk) {
      std::fprintf(stderr, "UI defaults or advanced-parameter grouping are invalid.\n");
      return 10;
    }
    localization->setCurrentIndex(1);
    const bool curvatureControlsOk = beadHeight->isEnabled() &&
        beadWidth->isEnabled() && spline->isEnabled() && turn->isEnabled() &&
        turn->currentIndex()==4 && !regionLength->isEnabled() && !regionFraction->isEnabled();
    localization->setCurrentIndex(0);
    if (!curvatureControlsOk) {
      std::fprintf(stderr, "Continuous-curvature controls are invalid.\n");
      return 11;
    }
    auto *length = window.findChild<QDoubleSpinBox *>("finalSmoothingLength");
    auto *fallbackRadius = window.findChild<QDoubleSpinBox *>("fallbackSplineRadius");
    auto *reset = window.findChild<QPushButton *>("resetAlgorithmDefaults");
    auto *merge = window.findChild<QSpinBox *>("profileMergeFrameCount");
    auto *depth = window.findChild<QDoubleSpinBox *>("minimumFeatureDepth");
    auto *uniform = window.findChild<QCheckBox *>("uniformWidthEnabled");
    auto *metric = window.findChild<QComboBox *>("uniformWidthMetric");
    auto *widthFilter=window.findChild<QCheckBox *>("uniformWidthFilterEnabled");
    auto *widthFilterLength=window.findChild<QDoubleSpinBox *>("uniformWidthFilterLength");
    if(!widthFilter || !widthFilterLength || !widthFilter->isChecked() || widthFilterLength->value()!=12)return 27;
    widthFilter->setChecked(false);widthFilterLength->setValue(21);
    if(!uniform || !metric || uniform->isChecked() || metric->currentIndex()!=1)return 23;
    uniform->setChecked(true);metric->setCurrentIndex(0);
    if (!depth || depth->value()!=1) return 17;
    depth->setValue(2);
    if (!length || !fallbackRadius || !reset || !merge) return 12;
    radius->setValue(6.5); length->setValue(17); fallbackRadius->setValue(4.5); merge->setValue(7);
    localization->setCurrentIndex(1);
    if (radius->value()!=8 || length->value()!=12 || merge->value()!=5) return 13;
    if (depth->value()!=1) return 18;
    if(uniform->isChecked() || metric->currentIndex()!=1)return 24;
    if(!widthFilter->isChecked() || widthFilterLength->value()!=12)return 28;
    depth->setValue(3);
    radius->setValue(11); length->setValue(23); merge->setValue(9);
    spline->setChecked(false);
    localization->setCurrentIndex(0);
    if (radius->value()!=6.5 || length->value()!=17 || fallbackRadius->value()!=4.5 || merge->value()!=7) return 14;
    if (depth->value()!=2) return 19;
    if(!uniform->isChecked() || metric->currentIndex()!=0)return 25;
    if(widthFilter->isChecked() || widthFilterLength->value()!=21)return 29;
    localization->setCurrentIndex(1);
    if (radius->value()!=11 || length->value()!=23 || merge->value()!=9 || spline->isChecked()) return 15;
    if (depth->value()!=3) return 20;
    reset->click();
    if (depth->value()!=1) return 21;
    localization->setCurrentIndex(0);
    if (radius->value()!=6.5 || length->value()!=17 || fallbackRadius->value()!=4.5 || merge->value()!=7) return 16;
    if (depth->value()!=2) return 22;
    reset->click();
    std::fprintf(stderr,"Independent algorithm parameters and current-only reset passed.\n");
    if(uniform->isChecked() || metric->currentIndex()!=1)return 26;
    if(!widthFilter->isChecked() || widthFilterLength->value()!=12)return 30;
    localization->setCurrentIndex(1);
  }
  stage("Before show");
  window.show();
  stage("Window shown");
  const auto args = app.arguments();
  window.setDebugFlag(args.contains("--debug-flag"));
  if(smoke)for(const auto &arg:args){std::fprintf(stderr,"arg: %s\n",arg.toLocal8Bit().constData());std::fflush(stderr);}
  auto option = [&](const QString &name) {
    const int i = args.indexOf(name);
    return i >= 0 && i + 1 < args.size() ? args[i + 1] : QString();
  };
  // 打包自检入口：真实加载与提取，保存窗口截图后退出，不增加产品界面按钮。
  if (args.contains("--smoke-test")) {
    stage("Smoke branch entered");
    const auto input = option("--input"), shot = option("--screenshot"), output = option("--output");
    if (option("--localization") == "curvature")
      window.findChild<QComboBox *>("profileLocalizationMethod")->setCurrentIndex(1);
    else if (option("--localization") == "curve-line")
      window.findChild<QComboBox *>("profileLocalizationMethod")->setCurrentIndex(0);
    if(!option("--uniform-width").isEmpty()) {
      window.findChild<QCheckBox *>("uniformWidthEnabled")->setChecked(true);
      window.findChild<QComboBox *>("uniformWidthMetric")->setCurrentIndex(option("--uniform-width")=="arc"?0:1);
    }
    const auto debugStart = option("--debug-start");
    const auto debugJump = option("--debug-jump");
    const auto debugJumpBack = option("--debug-jump-back");
    // Exercise the real modal VTK loop in unattended integration checks.
    if (args.contains("--smoke-auto-step") && args.contains("--debug-flag")) {
      if (auto *batch = window.findChild<QSpinBox *>("debugFrameBatchSize"))
        batch->setValue(option("--debug-frames").isEmpty() ? 1 : option("--debug-frames").toInt());
      auto *advance = new QTimer(&app);
      QObject::connect(advance, &QTimer::timeout, &app,
                       [shot, debugJump, debugJumpBack, jumped = false,
                        jumpedBack = false]() mutable {
        auto *dialog = qobject_cast<QDialog *>(QApplication::activeModalWidget());
        if (!dialog || dialog->objectName() != "seamDebugViewer") return;
        if (!shot.isEmpty()) {
          const auto suffix = dialog->windowTitle().contains(QStringLiteral("H 第一步"))
                                  ? ".body.png"
                                  : dialog->windowTitle().contains(QStringLiteral("F 按"))
                                        ? ".frames.png" : ".step.png";
          dialog->grab().save(shot + suffix);
        }
        std::fprintf(stderr, "Debug pause: %s\n", dialog->windowTitle().toUtf8().constData());
        if (!jumped && !debugJump.isEmpty() &&
            dialog->windowTitle().contains(QStringLiteral("F 按"))) {
          if (auto *selector = dialog->findChild<QSpinBox *>("debugFrameSelector"))
            selector->setValue(debugJump.toInt());
          if (auto *jump = dialog->findChild<QPushButton *>("jumpDebugFrame")) {
            jumped = true;
            jump->click();
            return;
          }
        }
        if (jumped && !jumpedBack && !debugJumpBack.isEmpty() &&
            dialog->windowTitle().contains(QStringLiteral("F 按"))) {
          if (auto *selector = dialog->findChild<QSpinBox *>("debugFrameSelector"))
            selector->setValue(debugJumpBack.toInt());
          if (auto *jump = dialog->findChild<QPushButton *>("jumpDebugFrame")) {
            jumpedBack = true;
            jump->click();
            return;
          }
        }
        if (auto *next = dialog->findChild<QPushButton *>("nextDebugStep")) next->click();
      });
      advance->start(200);
    }
    QTimer::singleShot(180000, &app, [&app] { app.exit(9); });
    QObject::connect(&window, &weld::MainWindow::extractionFinished, &app,
                     [&window, &app, shot, output](bool ok) {
                       QTimer::singleShot(800, &app, [&window, &app, shot, output, ok] {
                         // Verify real renderer actors disappear/reappear when
                         // each independent face toggle changes after detection.
                         bool layersOk=ok;
                         auto *widget=window.findChild<QVTKOpenGLNativeWidget *>();
                         auto *renderer=widget?widget->GetRenderWindow()->GetRenderers()->GetFirstRenderer():nullptr;
                         if(!renderer)layersOk=false;
                         if(layersOk)for(const auto *name:{"showLeftSurface","showRightSurface","showWeldSurface"}) {
                           auto *box=window.findChild<QCheckBox *>(name);
                           if(!box || !box->isChecked()){layersOk=false;break;}
                           const auto before=renderer->GetActors()->GetNumberOfItems();
                           box->setChecked(false);
                           const auto hidden=renderer->GetActors()->GetNumberOfItems();
                           box->setChecked(true);
                           const auto restored=renderer->GetActors()->GetNumberOfItems();
                           if(hidden!=before-1 || restored!=before)layersOk=false;
                           std::fprintf(stderr,"Surface toggle %s: actors %d -> %d -> %d\n",name,before,hidden,restored);
                         }
                         for (auto *scroll : window.findChildren<QScrollArea *>()) {
                           scroll->verticalScrollBar()->setValue(0);
                           scroll->horizontalScrollBar()->setValue(0);
                         }
                         const bool saved =
                             shot.isEmpty() || window.grab().save(shot);
                         QString error;
                         const bool exported=output.isEmpty() || window.saveResultToPath(output,&error);
                         for(auto *log:window.findChildren<QPlainTextEdit *>())
                           std::fprintf(stderr,"%s\n",log->toPlainText().toUtf8().constData());
                         if(!exported)std::fprintf(stderr,"Export failed: %s\n",error.toUtf8().constData());
                         app.exit(ok && layersOk && saved && exported ? 0 : 8);
                       });
                     });
    const bool generateSurfaceModels = !args.contains("--no-surface-models");
    // Apply output options once the GUI event loop can service redraws.
    QTimer::singleShot(500, &app, [&window, &app, input, debugStart, generateSurfaceModels] {
      window.findChild<QCheckBox *>("generateSurfaceModels")->setChecked(generateSurfaceModels);
      std::fprintf(stderr,"Loading input\n");std::fflush(stderr);
      if (!window.loadPath(input)) {
        app.exit(7);
        return;
      }
      if (!debugStart.isEmpty())
        if (auto *start = window.findChild<QSpinBox *>("debugStartSlice"))
          start->setValue(debugStart.toInt());
      window.runDetection();
      std::fprintf(stderr,"Detection dispatched\n");std::fflush(stderr);
    });
  }
  return app.exec();
}
