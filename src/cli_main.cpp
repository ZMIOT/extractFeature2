#include "PointCloudIO.h"
#include "ResultIO.h"
#include "SeamDebugViewer.h"
#include "SeamDetector.h"
#include "AlgorithmParameters.h"
#include "ConfigIO.h"
#include <QFileInfo>
#include <QApplication>
#include <QCoreApplication>
#include <QElapsedTimer>
#include <QFile>
#include <QMap>
#include <QSurfaceFormat>
#include <QTextStream>
#include <QVTKOpenGLNativeWidget.h>
#include <cmath>
#include <cstring>
#include <exception>
#include <memory>

int main(int argc, char **argv) {
  bool debugFlag = false;
  for (int i = 1; i < argc; ++i)
    if (std::strcmp(argv[i], "--debug-flag") == 0)
      debugFlag = true;
  std::unique_ptr<QCoreApplication> app;
  if (debugFlag) {
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
    auto *gui = new QApplication(argc, argv);
    gui->setQuitOnLastWindowClosed(false);
    app.reset(gui);
  } else {
    app.reset(new QCoreApplication(argc, argv));
  }
  QTextStream out(stdout), err(stderr);
  const auto args = app->arguments();
  if (args.contains("--help") || args.size() == 1) {
    out << "Scheme 1 weld seam extraction\nUsage: weld_seam_cli --input "
           "cloud.txt --output results/seams.pcd\n"
        << "  --axis x|y|z  --slice-step 2  --slice-thickness 3.5\n"
        << "  --config path.json  Local SDK config (default: bin/config/weld_seam.json)\n"
        << "  Explicit CLI options override the selected algorithm's file settings.\n"
        << "  --profile-merge-frames 5  --no-profile-frame-merge (original frames only)\n"
        << "  --floor-patch-size 6  --floor-patch-guard 0.6  "
           "--floor-max-correction 5\n"
        << "  --voxel-size 0.8  --profile-bin 1.2  --floor-patch-rms 0.4\n"
        << "  --profile-mode legacy-rolling|stable-grid  "
           "--profile-reference-scale 1.2\n"
        << "  --min-separation N (default: curve-line 8, curvature 6)  --no-floor-constraint\n"
        << "  --surface-refinement  Enable wall fitting/projection and floor intersection\n"
        << "  --no-final-smoothing  Disable final continuous trajectory fit\n"
        << "  --uniform-width arc|chord  Trim each width tail by 10%, then place a uniform pair along the bottom\n"
        << "  --no-uniform-width-filter  --uniform-width-filter-length 12\n"
        << "  --no-endpoint-trim  --endpoint-window 12  --endpoint-angle 20  --endpoint-deviation 1.5\n"
        << "  --final-smoothing-length 12  --final-robust-scale 0.6\n"
        << "  --profile-curvature-radius 8  --no-profile-curvature\n"
        << "  --localization curve-line|curvature (default curvature)\n"
        << "  Curvature only: --no-multi-candidates --curvature-window 2\n"
        << "    --curvature-prominence 0.08 --curvature-spacing 0.4\n"
        << "    --curvature-continuity 2\n"
        << "  Curve-line only: --no-curve-line-fallback --fallback-radius 8\n"
        << "    --fallback-wall-tolerance 0.6 --fallback-region-length 4\n"
        << "    --fallback-region-fraction 0.6 --fallback-window 2\n"
        << "    --fallback-prominence 0.08 --fallback-spacing 0.4 --fallback-continuity 2\n"
        << "  --central-bead-height 0.8 --central-bead-width 2\n"
        << "  --minimum-feature-depth 1  Stop at shallow sections; 0 disables\n"
        << "  --profile-turn-location peak|region-median|turn-median|floor-side-boundary|wall-turn-centroid\n"
        << "    Curvature default: wall-turn-centroid; curve-line default: region-median.\n"
        << "  --profile-turn-region-fraction 0.6 --profile-turn-region-length 4\n"
        << "  --exhaustive-candidates  Disable wall-region pruning for comparison\n"
        << "  --profile-wall-tolerance 0.6  Maximum local candidate distance to wall\n"
        << "  --double-seam-width 13  Whole-seam majority vote (ties: single)\n"
        << "  --width-rule  Enable whole-seam width vote (default: always two tracks)\n"
        << "  --opening positive-v|negative-v|either --candidate-max-correction 1\n"
        << "  --postprocess legacy|robust-rts  --rts-process-noise 0.04\n"
        << "  --rts-measurement-stddev 0.6  --rts-gate 9.21\n"
        << "  --rts-max-correction 2  --rts-geometry-tolerance 0.6\n"
        << "  --debug-flag  Show each processing step in a VTK window\n"
        << "  --debug-frames N  Pause every N frames/sections (default 1)\n"
        << "  --debug-start N   Start frame/section popups at 1-based position N\n"
        << "  --debug-stages    Pause at every algorithm stage\n"
        << "  --hull-voxel 1 --hull-distance 0.6 --hull-margin 3\n"
        << "Frame-indexed TXT/CSV: each original frame is one profile "
           "(auto-detected).\n"
        << "  Frame and spatial-slice inputs share the same detection "
           "parameters.\n"
        << "  --slice-step and --slice-thickness are ignored for original "
           "frames.\n"
        << "Outputs: use .pcd for an export prefix or .csv for the table path.\n"
        << "  Separate _left.pcd / _right.pcd (or _center.pcd), CSV and parameter JSON; no merged PCD.\n";
    return 0;
  }
  const int methodArgument = args.indexOf("--localization");
  weld::LibraryConfig libraryConfig;
  const int configArgument = args.indexOf("--config");
  const QString configPath = configArgument >= 0 ? args.value(configArgument + 1)
      : QCoreApplication::applicationDirPath() + "/config/weld_seam.json";
  if (configArgument >= 0 && (configPath.isEmpty() || configPath.startsWith("--"))) {
    err << "--config requires a file path.\n"; return 2;
  }
  if (configArgument >= 0 || QFileInfo::exists(configPath)) {
    std::string detail;
    if (!weld::loadLibraryConfig(configPath.toUtf8().constData(), libraryConfig, &detail)) {
      err << "Invalid configuration: " << QString::fromUtf8(detail.c_str()) << '\n'; return 2;
    }
  }
  const auto initialMethod = methodArgument >= 0 && args.value(methodArgument+1) == "curve-line"
      ? weld::ProfileLocalizationMethod::CurveLineIntersection
      : methodArgument >= 0 ? weld::ProfileLocalizationMethod::CurvatureCandidates : libraryConfig.activeAlgorithm;
  weld::DetectorConfig config = libraryConfig.algorithms.get(initialMethod);
  // The CLI creates a GUI application only when --debug-flag was explicit.
  config.debugFlag = debugFlag;
  QString input, output;
  const QMap<QString, float *> numeric{
      {"--curvature-window", &config.curvatureSearch.windowLength},
      {"--curvature-prominence", &config.curvatureSearch.valleyProminence},
      {"--curvature-spacing", &config.curvatureSearch.candidateSpacing},
      {"--curvature-continuity", &config.curvatureSearch.continuityWeight},
      {"--fallback-radius", &config.curveLineFallback.splineRadius},
      {"--fallback-wall-tolerance", &config.curveLineFallback.wallTolerance},
      {"--fallback-region-length", &config.curveLineFallback.regionLength},
      {"--fallback-region-fraction", &config.curveLineFallback.regionFraction},
      {"--fallback-window", &config.curveLineFallback.search.windowLength},
      {"--fallback-prominence", &config.curveLineFallback.search.valleyProminence},
      {"--fallback-spacing", &config.curveLineFallback.search.candidateSpacing},
      {"--fallback-continuity", &config.curveLineFallback.search.continuityWeight},
      {"--endpoint-window", &config.endpointTrimWindow},
      {"--endpoint-angle", &config.endpointTrimAngleDegrees},
      {"--endpoint-deviation", &config.endpointTrimDeviation},
      {"--final-smoothing-length", &config.finalSmoothingLength},
      {"--uniform-width-filter-length", &config.uniformWidthFilterLength},
      {"--profile-curvature-radius", &config.profileCurvatureRadius},
      {"--central-bead-height", &config.centralBeadMinHeight},
      {"--minimum-feature-depth", &config.minimumFeatureDepth},
      {"--central-bead-width", &config.centralBeadMinWidth},
      {"--profile-turn-region-fraction", &config.profileTurnRegionFraction},
      {"--profile-turn-region-length", &config.profileTurnRegionMaxLength},
      {"--profile-wall-tolerance", &config.profileWallTolerance},
      {"--double-seam-width", &config.doubleSeamWidth},
      {"--final-robust-scale", &config.finalSmoothingRobustScale},
      {"--hull-voxel", &config.hullSurfaceVoxelSize},
      {"--hull-distance", &config.hullSurfaceDistance},
      {"--hull-margin", &config.hullPreservationMargin},
      {"--candidate-max-correction", &config.maxCandidateCorrection},
      {"--slice-step", &config.sliceStep},
      {"--slice-thickness", &config.sliceThickness},
      {"--floor-patch-size", &config.floorPatchSize},
      {"--floor-patch-guard", &config.floorPatchGuard},
      {"--floor-max-correction", &config.floorPatchMaxCorrection},
      {"--voxel-size", &config.voxelSize},
      {"--profile-bin", &config.profileBinSize},
      {"--profile-reference-scale", &config.profileReferenceScale},
      {"--floor-patch-rms", &config.floorPatchMaxRms},
      {"--min-separation", &config.minSeamSeparation},
      {"--rts-process-noise", &config.rtsProcessNoise},
      {"--rts-measurement-stddev", &config.rtsMeasurementStddev},
      {"--rts-gate", &config.rtsInnovationGate},
      {"--rts-max-correction", &config.rtsMaxCorrection},
      {"--rts-geometry-tolerance", &config.rtsGeometryTolerance}};
  for (int i = 1; i < args.size(); ++i) {
    const QString key = args[i];
    if (key == "--config") { ++i; continue; }
    if (key == "--no-multi-candidates") { config.curvatureSearch.multiCandidates=false; continue; }
    if (key == "--no-curve-line-fallback") { config.curveLineFallback.enabled=false; continue; }
    if (key == "--localization") {
      const QString value = i+1 < args.size() ? args[++i] : QString();
      if (value == "curve-line")
        config.profileLocalizationMethod = weld::ProfileLocalizationMethod::CurveLineIntersection;
      else if (value == "curvature")
        config.profileLocalizationMethod = weld::ProfileLocalizationMethod::CurvatureCandidates;
      else { err << "--localization must be curve-line or curvature.\n"; return 2; }
      continue;
    }
    if (key == "--profile-turn-location") {
      const QString value = i+1 < args.size() ? args[++i] : QString();
      if (value == "peak") config.profileTurnLocation = weld::ProfileTurnLocation::Peak;
      else if (value == "region-median") config.profileTurnLocation = weld::ProfileTurnLocation::RegionMedian;
      else if (value == "turn-median") config.profileTurnLocation = weld::ProfileTurnLocation::TurnMedian;
      else if (value == "floor-side-boundary") config.profileTurnLocation = weld::ProfileTurnLocation::FloorSideBoundary;
      else if (value == "wall-turn-centroid") config.profileTurnLocation = weld::ProfileTurnLocation::WallTurnCentroid;
      else { err << "Invalid --profile-turn-location. See --help.\n"; return 2; }
      continue;
    }
    if (key == "--no-profile-frame-merge") { config.mergeProfileFrames = false; continue; }
    if (key == "--profile-merge-frames") {
      bool ok = false;
      const int count = i + 1 < args.size() ? args[++i].toInt(&ok) : 0;
      if (!ok || count < 1 || count > 101 || count % 2 == 0) {
        err << "--profile-merge-frames requires an odd integer from 1 to 101.\n";
        return 2;
      }
      config.profileMergeFrameCount = count;
      continue;
    }
    if (key == "--exhaustive-candidates") { config.useWallCandidatePruning=false; continue; }
    if (key == "--no-endpoint-trim") { config.trimUnstableEndpoints=false; continue; }
    if (key == "--no-uniform-width-filter") {config.uniformWidthFilterEnabled=false;continue;}
    if (key == "--uniform-width") {
      const QString value=i+1<args.size()?args[++i].toLower():QString();
      if(value!="arc" && value!="chord") {err<<"--uniform-width must be arc or chord.\n";return 2;}
      config.uniformWidthEnabled=true;
      config.uniformWidthMetric=value=="arc"?weld::UniformWidthMetric::BottomArc:weld::UniformWidthMetric::Chord;
      continue;
    }
    if (key == "--debug-flag")
      continue;
    if (key == "--opening") {
      const QString value = i + 1 < args.size() ? args[++i].toLower() : QString();
      if (value == "positive-v") config.interiorOpeningDirection = 1;
      else if (value == "negative-v") config.interiorOpeningDirection = -1;
      else if (value == "either") config.interiorOpeningDirection = 0;
      else { err << "--opening must be positive-v, negative-v or either.\n"; return 2; }
      continue;
    }
    if (key == "--no-profile-curvature") {
      config.useProfileCurvature = false;
      continue;
    }
    if (key == "--no-width-rule") { config.adaptiveSeamCount=false; continue; }
    if (key == "--width-rule") { config.adaptiveSeamCount=true; continue; }
    if (key == "--debug-stages") {
      config.debugFrameBatchSize = 0;
      continue;
    }
    if (key == "--debug-frames") {
      bool ok = false;
      const int count = i + 1 < args.size() ? args[++i].toInt(&ok) : 0;
      if (!ok || count < 1 || count > 1000) {
        err << "--debug-frames requires an integer from 1 to 1000.\n";
        return 2;
      }
      config.debugFrameBatchSize = count;
      continue;
    }
    if (key == "--debug-start") {
      bool ok = false;
      const int start = i + 1 < args.size() ? args[++i].toInt(&ok) : 0;
      if (!ok || start < 1) {
        err << "--debug-start requires a positive 1-based frame/section position.\n";
        return 2;
      }
      config.debugStartSlice = start;
      continue;
    }
    if (key == "--no-final-smoothing") {
      config.finalSmoothingEnabled = false;
      continue;
    }
    if (key == "--no-surface-refinement") {
      config.useSurfaceRefinement = false;
      continue;
    }
    if (key == "--surface-refinement") {
      config.useSurfaceRefinement = true;
      continue;
    }
    if (key == "--no-floor-constraint") {
      config.useLocalFloorConstraint = false;
      continue;
    }
    if (key == "--postprocess") {
      if (i + 1 >= args.size() || args[i + 1].startsWith("--")) {
        err << "Missing value for --postprocess\n";
        return 2;
      }
      const QString value = args[++i].toLower();
      if (value == "legacy")
        config.postprocessMode = weld::SeamPostprocessMode::Legacy;
      else if (value == "robust-rts")
        config.postprocessMode = weld::SeamPostprocessMode::RobustRts;
      else {
        err << "Postprocess must be legacy or robust-rts.\n";
        return 2;
      }
      continue;
    }
    if (key == "--profile-mode") {
      if (i + 1 >= args.size() || args[i + 1].startsWith("--")) {
        err << "Missing value for --profile-mode\n";
        return 2;
      }
      const QString value = args[++i].toLower();
      if (value == "legacy-rolling")
        config.profileSamplingMode =
            weld::ProfileSamplingMode::LegacyRolling;
      else if (value == "stable-grid")
        config.profileSamplingMode = weld::ProfileSamplingMode::StableGrid;
      else {
        err << "Profile mode must be legacy-rolling or stable-grid.\n";
        return 2;
      }
      continue;
    }
    if (key != "--input" && key != "--output" && key != "--axis" &&
        !numeric.contains(key)) {
      err << "Unknown option: " << key << '\n';
      return 2;
    }
    if (i + 1 >= args.size() || args[i + 1].startsWith("--")) {
      err << "Missing value for " << key << '\n';
      return 2;
    }
    const QString value = args[++i];
    if (key == "--input")
      input = value;
    else if (key == "--output")
      output = value;
    else if (key == "--axis") {
      const int axis = QStringList{"x", "y", "z"}.indexOf(value.toLower());
      if (axis < 0) {
        err << "Axis must be x, y or z.\n";
        return 2;
      }
      config.longitudinalAxis = weld::Axis(axis);
    } else {
      bool ok = false;
      const float number = value.toFloat(&ok);
      if (!ok || !std::isfinite(number)) {
        err << "Invalid number for " << key << '\n';
        return 2;
      }
      *numeric[key] = number;
    }
  }
  for (const auto &arg : args) {
    const bool curveLine = config.profileLocalizationMethod == weld::ProfileLocalizationMethod::CurveLineIntersection;
    if ((curveLine && (arg.startsWith("--curvature-") || arg == "--no-multi-candidates")) ||
        (!curveLine && (arg.startsWith("--fallback-") || arg == "--no-curve-line-fallback"))) {
      err << "Option " << arg << " belongs to the other localization algorithm.\n";
      return 2;
    }
  }
  if (config.profileLocalizationMethod ==
      weld::ProfileLocalizationMethod::CurveLineIntersection)
    config.useProfileCurvature = true;
  if (input.isEmpty() || output.isEmpty()) {
    err << "--input and --output are required.\n";
    return 2;
  }
  std::string configError;
  if (!weld::validateDetectorConfig(config, &configError)) {
    err << QString::fromUtf8(configError.c_str()) << '\n'; return 2;
  }
  try {
    std::string error;
    std::vector<weld::ScanFrameRange> frames;
    const auto cloud = weld::PointCloudIO::load(input.toLocal8Bit().constData(),
                                                &error, &frames);
    if (!cloud || cloud->empty()) {
      err << QString::fromLocal8Bit(error.c_str()) << '\n';
      return 3;
    }
    if (!frames.empty())
      out << "Original-frame mode: " << frames.size()
          << " frames. Profile projection window: "
          << (config.mergeProfileFrames ? config.profileMergeFrameCount : 1)
          << " frames; one output cross-section per original frame.\n";
    QElapsedTimer timer;
    timer.start();
    // 与界面共用入口；仅 debugFlag 开启时同步展示 VTK 步骤窗口。
    const auto result =
        weld::detectWithDebugViewer(cloud, config, nullptr, frames);
    out << QString::fromStdString(result.summary)
        << "\nElapsed ms: " << timer.elapsed() << '\n';
    for (const auto &stage : result.stageMilliseconds)
      out << "  " << QString::fromStdString(stage.first) << ": " << stage.second << " ms\n";
    if (!result.ok())
      return 5;
    QString saveError;
    if (!weld::saveResults(output, result, config, &saveError)) {
      err << saveError << '\n';
      return 6;
    }
    out << "Saved tracks/surfaces PCD, three B-spline surface model JSON files, CSV and parameters: " << output << '\n';
    return 0;
  } catch (const std::exception &e) {
    err << e.what() << '\n';
    return 7;
  }
}
