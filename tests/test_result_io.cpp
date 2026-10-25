#include "ResultIO.h"
#include "PointCloudIO.h"
#include <QCoreApplication>
#include <QFile>
#include <QTemporaryDir>
#include <QJsonDocument>
#include <QJsonArray>
#include <iostream>
#include <stdexcept>
#include <cmath>

namespace {
void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}
weld::CloudT::Ptr readPcd(const QString &path) {
  QFile file(path);
  require(file.open(QIODevice::ReadOnly), "Requested PCD path was not written");
  require(file.read(512).contains("DATA binary\n"), "PCD path does not contain a binary PCD");
  std::string error;
  auto cloud=weld::PointCloudIO::load(path.toLocal8Bit().constData(), &error);
  require(cloud && !cloud->empty(), "Saved PCD cannot be read by PCL");
  return cloud;
}
}

int main(int argc, char **argv) {
  QCoreApplication app(argc, argv);
  try {
    QTemporaryDir dir;
    require(dir.isValid(), "Cannot create export test directory");
    weld::SeamResult result;
    weld::DetectorConfig config;
    config.uniformWidthEnabled=true;
    result.uniformWidthTarget=10.05f;result.uniformWidthSamples=30;
    result.uniformWidthTrimEach=3;result.uniformWidthAdjusted=3;
    result.frameBased=true;
    for (int i=0;i<3;++i) {
      result.leftSeam->push_back({float(i),-4.25f,2.5f});
      result.rightSeam->push_back({float(i),5.75f,3.5f});
      result.leftSurface->push_back({float(i),-8.f,8.f});
      result.rightSurface->push_back({float(i),8.f,8.f});
      result.weldSurface->push_back({float(i),0.f,1.f});
      result.leftFrameIndices.push_back(10+i);
      result.rightFrameIndices.push_back(10+i);
      result.leftPointConfidence.push_back(.9f);
      result.rightPointConfidence.push_back(.9f);
    }
    QString error;
    weld::BSplineSurfacePatch patch;
    patch.countU=patch.countV=4;
    patch.knotsU=patch.knotsV={0,0,0,0,1,1,1,1};
    for(int i=0;i<4;++i)for(int j=0;j<4;++j)patch.controlPoints.push_back({double(i),double(j),.1*i*j});
    result.leftSurfaceModel.patches.push_back(patch);
    result.weldSurfaceModel=result.rightSurfaceModel=result.leftSurfaceModel;
    const auto path=dir.path()+"/test_seams.pcd";
    require(weld::saveResults(path,result,config,&error), "Double-track PCD save failed");
    require(!QFile::exists(path), "Export unexpectedly merged the two seam tracks");
    QFile metadata(dir.path()+"/test_params.json");
    require(metadata.open(QIODevice::ReadOnly),"Uniform-width metadata was not saved");
    const auto json=QJsonDocument::fromJson(metadata.readAll()).object();
    metadata.close();
    require(std::abs(json["uniformWidthTarget"].toDouble()-10.05)<1.e-4 &&
            json["uniformWidthSamples"].toInt()==30 && json["uniformWidthTrimEach"].toInt()==3,
            "Uniform-width statistic metadata changed");
    const auto left=readPcd(dir.path()+"/test_left.pcd");
    const auto right=readPcd(dir.path()+"/test_right.pcd");
    require(left->size()==3 && right->size()==3, "Per-track PCD point counts differ");
    const QString faceNames[]={"left_surface","right_surface","weld_surface"};
    const weld::CloudT::Ptr faces[]={result.leftSurface,result.rightSurface,result.weldSurface};
    for(int face=0;face<3;++face) {
      const auto saved=readPcd(dir.path()+"/test_"+faceNames[face]+".pcd");
      require(saved->size()==faces[face]->size(),"Surface PCD point count changed");
      for(std::size_t i=0;i<saved->size();++i)
        require((weld::toVector((*saved)[i])-weld::toVector((*faces[face])[i])).norm()==0,
                "Surface PCD coordinates changed");
    }
    require(json["leftSurfacePoints"].toInt()==3 && json["rightSurfacePoints"].toInt()==3 &&
            json["weldSurfacePoints"].toInt()==3,"Surface metadata counts differ");
    for(const auto &face:faceNames) {
      QFile modelFile(dir.path()+"/test_"+face+"_model.json");
      require(modelFile.open(QIODevice::ReadOnly),"B-spline model file missing");
      const auto modelJson=QJsonDocument::fromJson(modelFile.readAll()).object();
      require(modelJson["valid"].toBool(),"Saved model lost valid status");
      const auto saved=modelJson["patches"].toArray()[0].toObject();
      weld::BSplineSurfacePatch loaded;
      loaded.countU=saved["countU"].toInt();loaded.countV=saved["countV"].toInt();
      loaded.degreeU=saved["degreeU"].toInt();loaded.degreeV=saved["degreeV"].toInt();
      for(const auto &k:saved["knotsU"].toArray())loaded.knotsU.push_back(k.toDouble());
      for(const auto &k:saved["knotsV"].toArray())loaded.knotsV.push_back(k.toDouble());
      for(const auto &cp:saved["controlPoints"].toArray()) {
        const auto xyz=cp.toArray();loaded.controlPoints.push_back({xyz[0].toDouble(),xyz[1].toDouble(),xyz[2].toDouble()});
      }
      Eigen::Vector3d original,restored;
      for(double u:{0.,.13,.5,1.})for(double v:{0.,.37,.9,1.}) {
        require(patch.evaluate(u,v,original) && loaded.evaluate(u,v,restored),"Saved model cannot be evaluated");
        require((original-restored).norm()<1e-12,"JSON round-trip changed the surface function");
      }
    }
    for (int i=0;i<3;++i) {
      require((weld::toVector((*left)[i])-weld::toVector((*result.leftSeam)[i])).norm()==0,
              "Left PCD coordinates changed during export");
      require((weld::toVector((*right)[i])-weld::toVector((*result.rightSeam)[i])).norm()==0,
              "Right PCD coordinates changed during export");
    }
    require(readPcd(dir.path()+"/test_left.pcd")->size()==3 &&
            readPcd(dir.path()+"/test_right.pcd")->size()==3, "Per-track PCD files missing");
    require(QFile::exists(dir.path()+"/test_seams.csv") && QFile::exists(dir.path()+"/test_params.json"),
            "Frame and parameter companion files missing");
    require(weld::saveResults(dir.path()+"/legacy_seams.csv",result,config,&error),
            "Legacy CSV export failed");
    require(readPcd(dir.path()+"/legacy_left.pcd")->size()==3, "Legacy export lost its PCD");
    require(weld::saveResults(dir.path()+"/no_suffix",result,config,&error), "Default extension failed");
    require(readPcd(dir.path()+"/no_suffix_left.pcd")->size()==3 &&
            readPcd(dir.path()+"/no_suffix_right.pcd")->size()==3,
            "Extensionless export did not save separate PCD tracks");
    require(weld::saveResults(dir.path()+"/test_left.pcd",result,config,&error),
            "Selecting an existing side file failed");
    require(!QFile::exists(dir.path()+"/test_left_left.pcd"), "Repeated save duplicated side suffix");
    for(const auto &name:faceNames) {
      require(weld::saveResults(dir.path()+"/test_"+name+".pcd",result,config,&error),
              "Selecting an existing surface file failed");
      require(!QFile::exists(dir.path()+"/test_"+name+"_left.pcd"),"Repeated save duplicated surface suffix");
    }

    *result.centerSeam=*result.leftSeam;
    result.centerFrameIndices=result.leftFrameIndices;
    result.centerPointConfidence=result.leftPointConfidence;
    result.leftSeam->clear(); result.rightSeam->clear();
    result.leftSurface->clear(); result.rightSurface->clear(); result.weldSurface->clear();
    result.leftSurfaceModel={};result.rightSurfaceModel={};result.weldSurfaceModel={};
    result.leftFrameIndices.clear(); result.rightFrameIndices.clear();
    require(weld::saveResults(path,result,config,&error), "Single-track PCD save failed");
    require(!QFile::exists(path) && readPcd(dir.path()+"/test_center.pcd")->size()==3,
            "Single-track export did not contain exactly the center seam");
    require(!QFile::exists(dir.path()+"/test_left.pcd") && !QFile::exists(dir.path()+"/test_right.pcd"),
            "Single-track overwrite retained stale double-track files");
    for(const auto &name:faceNames)
      require(!QFile::exists(dir.path()+"/test_"+name+".pcd"),"Empty surface retained a stale PCD");
    for(const auto &name:faceNames) {
      QFile f(dir.path()+"/test_"+name+"_model.json");require(f.open(QIODevice::ReadOnly),"Empty model status missing");
      const auto j=QJsonDocument::fromJson(f.readAll()).object();
      require(!j["valid"].toBool() && j["patches"].toArray().isEmpty(),"Overwrite retained stale surface model");
    }
    require(!weld::saveResults(dir.path()+"/bad.xyz",result,config,&error), "Unsupported extension accepted");
    std::cout << "PCD double/single track round-trip and CSV compatibility passed.\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
