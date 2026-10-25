#include "ResultIO.h"
#include "PointCloudIO.h"
#include <QDir>
#include <QFileInfo>
#include <QFile>
#include <QJsonDocument>
#include <QJsonArray>
#include <QSaveFile>

namespace weld {
QJsonObject surfaceModelJson(const BSplineSurfaceModel &model) {
  QJsonObject json;
  json["schemaVersion"]=1;
  json["type"]="non_rational_tensor_product_bspline";
  json["valid"]=model.valid();
  json["status"]=QString::fromStdString(model.message);
  json["rejectedPatchCount"]=double(model.rejectedPatchCount);
  json["coordinateSystem"]="input_xyz";
  json["units"]="same_as_input";
  json["controlPointOrder"]="u_major_index_i_times_countV_plus_j";
  json["knotConvention"]="full_clamped_with_endpoint_multiplicity_degree_plus_one";
  json["parameterization"]="u_scan_centroid_chord_length_v_section_contour_chord_length";
  json["errorMetric"]="euclidean_residual_at_assigned_parameters_not_closest_point_distance";
  QJsonArray patches;
  auto numbers=[](const std::vector<double> &values){QJsonArray a;for(double x:values)a.append(x);return a;};
  for(const auto &p:model.patches) {
    QJsonObject patch;
    patch["degreeU"]=p.degreeU;patch["degreeV"]=p.degreeV;
    patch["countU"]=p.countU;patch["countV"]=p.countV;
    patch["knotsU"]=numbers(p.knotsU);patch["knotsV"]=numbers(p.knotsV);
    patch["domainU"]=QJsonArray{0.,1.};patch["domainV"]=QJsonArray{0.,1.};
    QJsonArray points;for(const auto &cp:p.controlPoints)points.append(QJsonArray{cp[0],cp[1],cp[2]});
    patch["controlPoints"]=points;
    patch["firstFrame"]=p.firstFrame;patch["lastFrame"]=p.lastFrame;
    patch["firstStation"]=p.firstStation;patch["lastStation"]=p.lastStation;
    patch["sampleCount"]=double(p.sampleCount);
    patch["rmsError"]=p.rmsError;patch["maxError"]=p.maxError;
    patches.append(patch);
  }
  json["patches"]=patches;return json;
}
QJsonObject configJson(const DetectorConfig &c) {
  QJsonObject j;
  j["parameterSchemaVersion"] = 2;
  j["parameterScope"] = c.profileLocalizationMethod == ProfileLocalizationMethod::CurveLineIntersection
      ? "curve_line" : "curvature";
  auto searchJson = [](const CurvatureSearchConfig &s) {
    QJsonObject p;
    p["multiCandidates"] = s.multiCandidates;
    p["windowLength"] = s.windowLength;
    p["valleyProminence"] = s.valleyProminence;
    p["candidateSpacing"] = s.candidateSpacing;
    p["maximumPairs"] = s.maximumPairs;
    p["continuityWeight"] = s.continuityWeight;
    return p;
  };
  if (c.profileLocalizationMethod == ProfileLocalizationMethod::CurveLineIntersection) {
    QJsonObject fallback = searchJson(c.curveLineFallback.search);
    fallback["enabled"] = c.curveLineFallback.enabled;
    fallback["splineRadius"] = c.curveLineFallback.splineRadius;
    fallback["wallTolerance"] = c.curveLineFallback.wallTolerance;
    fallback["regionLength"] = c.curveLineFallback.regionLength;
    fallback["regionFraction"] = c.curveLineFallback.regionFraction;
    j["curveLineFallback"] = fallback;
  } else {
    j["curvatureSearch"] = searchJson(c.curvatureSearch);
  }
  j["algorithm"] = "scheme1";
  j["longitudinalAxis"] = int(c.longitudinalAxis);
  j["historyEnabled"] =
      false; // 独立版入口不导入历史；算法库保留历史接口供二次开发。
#define PARAM(name) j[#name] = c.name
  PARAM(voxelSize);
  PARAM(hullSurfaceVoxelSize);
  PARAM(hullSurfaceCurvature);
  PARAM(hullSurfaceDistance);
  PARAM(hullPreservationMargin);
  PARAM(debugFrameBatchSize);
  PARAM(debugStartSlice);
  PARAM(debugSelectedFrameIndex);
  PARAM(mergeProfileFrames);
  PARAM(profileMergeFrameCount);
  PARAM(interiorOpeningDirection);
  PARAM(maxCandidateCorrection);
  PARAM(sliceStep);
  PARAM(sliceThickness);
  PARAM(profileBinSize);
  j["profileSamplingMode"] =
      c.profileSamplingMode == ProfileSamplingMode::StableGrid
          ? "stable_grid"
          : "legacy_rolling";
  PARAM(profileReferenceScale);
  PARAM(useProfileCurvature);
  j["profileLocalizationMethod"] =
      c.profileLocalizationMethod == ProfileLocalizationMethod::CurveLineIntersection
          ? "curve_line_intersection"
          : "curvature_candidates";
  PARAM(useWallCandidatePruning);
  PARAM(profileCurvatureRadius);
  PARAM(centralBeadMinHeight);
  PARAM(centralBeadMinWidth);
  j["profileTurnLocation"] =
      c.profileTurnLocation == ProfileTurnLocation::Peak ? "peak" :
      c.profileTurnLocation == ProfileTurnLocation::RegionMedian ? "region_median" :
      c.profileTurnLocation == ProfileTurnLocation::TurnMedian ? "turn_median" :
      c.profileTurnLocation == ProfileTurnLocation::WallTurnCentroid ? "wall_turn_centroid" :
      "floor_side_boundary";
  PARAM(profileTurnRegionFraction);
  PARAM(profileTurnRegionMaxLength);
  PARAM(profileWallTolerance);
  PARAM(adaptiveSeamCount);
  PARAM(doubleSeamWidth);
  j["seamCountPolicy"] = c.uniformWidthEnabled ? "uniform_width_double" :
      c.adaptiveSeamCount ? "whole_seam_majority" : "always_double";
  j["seamCountTieBreak"] = "single";
  j["grooveWidthDefinition"]="transverse_distance_between_fitted_wall_bottom_junctions";
  PARAM(useSurfaceRefinement);
  PARAM(useLocalFloorConstraint);
  PARAM(floorPatchSize);
  PARAM(floorPatchGuard);
  PARAM(floorPatchMaxRms);
  PARAM(floorPatchMaxCorrection);
  PARAM(minSeamSeparation);
  PARAM(maxTransitionDistance);
  PARAM(curvatureWeight);
  PARAM(curveRegularizationWeight);
  PARAM(trimUnstableEndpoints);
  PARAM(minimumFeatureDepth);
  PARAM(endpointTrimAngleDegrees);
  PARAM(endpointTrimWindow);
  PARAM(endpointTrimDeviation);
  PARAM(finalSmoothingEnabled);
  PARAM(uniformWidthEnabled);
  PARAM(uniformWidthFilterEnabled);
  PARAM(uniformWidthFilterLength);
  j["uniformWidthMetric"] = c.uniformWidthMetric==UniformWidthMetric::BottomArc ? "bottom_arc" : "chord";
  j["uniformWidthTrimFractionEachTail"] = 0.1;
  PARAM(finalSmoothingLength);
  PARAM(finalSmoothingRobustScale);
  j["postprocessMode"] = c.postprocessMode == SeamPostprocessMode::RobustRts
                             ? "robust_rts"
                             : "legacy";
  PARAM(rtsProcessNoise);
  PARAM(rtsMeasurementStddev);
  PARAM(rtsInnovationGate);
  PARAM(rtsMaxCorrection);
  PARAM(rtsGeometryTolerance);
  PARAM(outlierMeanK);
  PARAM(outlierStddev);
  PARAM(fitHalfWindow);
  PARAM(regionGrowingDistance);
  PARAM(regionGrowingSmoothnessDegrees);
#undef PARAM
  // Preserve the remaining SDK parameters as well as the historical metadata
  // spelling of enum values above, so exports describe the actual run.
#define CONFIG_BOOL(name) if (!j.contains(#name)) j[#name] = c.name;
#define CONFIG_NUMBER(name, low, high) if (!j.contains(#name)) j[#name] = static_cast<double>(c.name);
#include "ConfigFields.inc"
#undef CONFIG_BOOL
#undef CONFIG_NUMBER
  return j;
}
bool saveResults(const QString &outputPath, const SeamResult &r,
                 const DetectorConfig &c, QString *error) {
  if (outputPath.trimmed().isEmpty()) {
    if (error) *error = QStringLiteral("请选择保存路径。");
    return false;
  }
  if (!r.ok()) {
    if (error)
      *error = QStringLiteral("没有可保存的焊缝轨迹。");
    return false;
  }
  const QString path = QFileInfo(outputPath).suffix().isEmpty() ? outputPath+".pcd" : outputPath;
  const QFileInfo info(path);
  const bool pcd = info.suffix().compare("pcd", Qt::CaseInsensitive)==0;
  if (!pcd && info.suffix().compare("csv", Qt::CaseInsensitive)!=0) {
    if (error) *error = QStringLiteral("请选择 PCD (*.pcd) 或 CSV (*.csv) 格式。");
    return false;
  }
  if (!QDir().mkpath(info.absolutePath())) {
    if (error)
      *error = QStringLiteral("无法创建输出目录。");
    return false;
  }
  QString name = info.completeBaseName();
  if (pcd && name.endsWith("_left_surface"))
    name.chop(13);
  else if (pcd && name.endsWith("_right_surface"))
    name.chop(14);
  else if (pcd && name.endsWith("_weld_surface"))
    name.chop(13);
  else if (name.endsWith("_seams"))
    name.chop(6);
  else if (name.endsWith("_seam"))
    name.chop(5);
  else if (pcd && name.endsWith("_left"))
    name.chop(5);
  else if (pcd && name.endsWith("_right"))
    name.chop(6);
  else if (pcd && name.endsWith("_center"))
    name.chop(7);
  const QString base = info.absolutePath() + "/" + name;
  std::string detail;
  const QString csvPath = pcd ? base+"_seams.csv" : path;
  const bool savedCsv = PointCloudIO::saveResultCsv(csvPath.toLocal8Bit().constData(),r,c.longitudinalAxis,&detail);
  auto saveSide=[&](const QString &suffix,const CloudT &cloud) {
    const QString path=base+suffix+".pcd";
    if (cloud.empty()) {
      if (QFile::exists(path) && !QFile::remove(path)) { detail="Cannot remove stale empty-mode output: "+path.toStdString(); return false; }
      return true;
    }
    return PointCloudIO::savePointCloud(path.toLocal8Bit().constData(),cloud,&detail);
  };
  if (!savedCsv || !saveSide("_left",*r.leftSeam) || !saveSide("_right",*r.rightSeam) ||
      !saveSide("_center",*r.centerSeam) ||
      !saveSide("_left_surface",*r.leftSurface) ||
      !saveSide("_right_surface",*r.rightSurface) ||
      !saveSide("_weld_surface",*r.weldSurface)) {
    if (error)
      *error = QString::fromLocal8Bit(detail.c_str());
    return false;
  }
  QJsonObject metadata = configJson(c);
  const std::pair<QString,const BSplineSurfaceModel *> models[]={
    {"left",&r.leftSurfaceModel},{"weld",&r.weldSurfaceModel},{"right",&r.rightSurfaceModel}};
  QJsonObject modelFiles;
  for(const auto &entry:models) {
    const QString modelPath=base+"_"+entry.first+"_surface_model.json";
    QSaveFile modelFile(modelPath);
    const auto data=QJsonDocument(surfaceModelJson(*entry.second)).toJson();
    if(!modelFile.open(QIODevice::WriteOnly) || modelFile.write(data)!=data.size() || !modelFile.commit()) {
      if(error)*error=modelPath+": "+modelFile.errorString();return false;
    }
    modelFiles[entry.first]=QFileInfo(modelPath).fileName();
  }
  metadata["surfaceModelFiles"]=modelFiles;
  if (r.frameBased) {
    metadata["algorithm"] = "scheme1_original_frames";
    metadata["inputFrameCount"] = static_cast<int>(r.inputFrameCount);
    metadata["spatialSlicing"] = false;
    metadata.remove("sliceStep");
    metadata.remove("sliceThickness");
  }
  metadata["leftPoints"] = int(r.leftSeam->size());
  metadata["leftSurfacePoints"] = int(r.leftSurface->size());
  metadata["rightSurfacePoints"] = int(r.rightSurface->size());
  metadata["weldSurfacePoints"] = int(r.weldSurface->size());
  metadata["mainBodyPoints"] = int(r.mainBodyCloud->size());
  metadata["removedBodyPoints"] = int(r.removedBodyPoints);
  metadata["endpointTrimmedStart"] = int(r.endpointTrimmedStart);
  metadata["endpointTrimmedEnd"] = int(r.endpointTrimmedEnd);
  metadata["candidateCorrectionRejected"] = int(r.candidateCorrectionRejected);
  metadata["rightPoints"] = int(r.rightSeam->size());
  auto fallbackCount = [](const std::vector<SeamPointQuality> &quality) {
    int count=0;
    for (const auto &q:quality) if (q.flags & QualityCurvatureFallback) ++count;
    return count;
  };
  metadata["leftFallbackPoints"] = fallbackCount(r.leftPointQuality);
  metadata["rightFallbackPoints"] = fallbackCount(r.rightPointQuality);
  metadata["centerFallbackPoints"] = fallbackCount(r.centerPointQuality);
  metadata["centerPoints"] = int(r.centerSeam->size());
  metadata["seamMode"] = r.centerSeam->empty() ? "double" : "single";
  metadata["minimumGrooveWidth"] = r.minimumGrooveWidth;
  metadata["maximumGrooveWidth"] = r.maximumGrooveWidth;
  metadata["confidence"] = r.confidence;
  metadata["uniformWidthTarget"] = r.uniformWidthTarget;
  metadata["uniformWidthSamples"] = int(r.uniformWidthSamples);
  metadata["uniformWidthTrimEach"] = int(r.uniformWidthTrimEach);
  metadata["uniformWidthAdjusted"] = int(r.uniformWidthAdjusted);
  metadata["uniformWidthOmitted"] = int(r.uniformWidthOmitted);
  metadata["uniformWidthFilterApplied"] = r.uniformWidthFilterApplied;
  metadata["uniformWidthFilterMeanCorrection"] = r.uniformWidthFilterMeanCorrection;
  metadata["uniformWidthFilterMaxCorrection"] = r.uniformWidthFilterMaxCorrection;
  metadata["summary"] = QString::fromStdString(r.summary);
  auto finalCurve = [](const FinalCurveSummary &s) {
    QJsonObject j;
    j["applied"] = s.applied;
    j["outlierPoints"] = int(s.outlierPoints);
    j["meanCorrection"] = s.meanCorrection;
    j["maxCorrection"] = s.maxCorrection;
    j["wallConstrainedPoints"] = int(s.wallConstrainedPoints);
    j["maximumWallResidual"] = s.maximumWallResidual;
    return j;
  };
  metadata["leftFinalCurve"] = finalCurve(r.leftFinalCurve);
  metadata["rightFinalCurve"] = finalCurve(r.rightFinalCurve);
  QJsonObject timings;
  for (const auto &stage : r.stageMilliseconds)
    timings[QString::fromStdString(stage.first)] = stage.second;
  metadata["stageMilliseconds"] = timings;
  QSaveFile file(base + "_params.json");
  const auto bytes = QJsonDocument(metadata).toJson();
  if (!file.open(QIODevice::WriteOnly) || file.write(bytes) != bytes.size() ||
      !file.commit()) {
    if (error)
      *error = file.errorString();
    return false;
  }
  return true;
}
} // namespace weld
