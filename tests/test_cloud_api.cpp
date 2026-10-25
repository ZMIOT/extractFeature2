#include "WeldSeam.h"
#include "SeamDetector.h"
#include "PointCloudIO.h"
#include <iostream>
#include <cmath>

#define CHECK(x) do { if (!(x)) { std::cerr << "Failed line " << __LINE__ << ": " << #x << '\n'; return 1; } } while (0)
bool same(const weld::CloudT &a, const weld::CloudT &b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i=0; i<a.size(); ++i)
    if (a[i].x != b[i].x || a[i].y != b[i].y || a[i].z != b[i].z) return false;
  return true;
}
int main() {
  std::string error;
  std::vector<weld::ScanFrameRange> frames;
  const auto input = weld::PointCloudIO::load(
      std::string(SEAM_TEST_DATA_DIR) + "/curvature_2_frames_2040_2110.txt", &error, &frames);
  CHECK(input && !input->empty() && frames.size() == 71);
  const weld::CloudT unchanged(*input);
  auto config = weld::defaultDetectorConfig();
  CHECK(config.profileLocalizationMethod == weld::ProfileLocalizationMethod::CurvatureCandidates);
  config.extractMainBody = false; config.voxelSize = 0; config.outlierMeanK = 0;
  config.trimUnstableEndpoints = false;
  for (bool single : {false, true}) {
    config.adaptiveSeamCount = single;
    config.doubleSeamWidth = single ? 1000.f : 13.f;
    const auto expected = weld::SeamDetector().detect(input, config, nullptr, nullptr, nullptr, &frames);
    CHECK(expected.ok());
    weld::SeamClouds output;
    error = "stale error";
    CHECK(weld::extractSeams(input, output, config, &error, &frames));
    CHECK(error.empty());
    CHECK(same(*output.left, *expected.leftSeam));
    CHECK(same(*output.right, *expected.rightSeam));
    CHECK(same(*output.center, *expected.centerSeam));
    CHECK(same(*output.leftSurface, *expected.leftSurface));
    CHECK(same(*output.rightSurface, *expected.rightSurface));
    CHECK(same(*output.weldSurface, *expected.weldSurface));
    CHECK(!output.leftSurface->empty() && !output.rightSurface->empty() && !output.weldSurface->empty());
    for(const auto *model:{&output.leftSurfaceModel,&output.weldSurfaceModel,&output.rightSurfaceModel})
      std::cout<<"model patches="<<model->patches.size()<<" rejected="<<model->rejectedPatchCount<<" status="<<model->message<<'\n';
    CHECK(output.leftSurfaceModel.valid() && output.weldSurfaceModel.valid() && output.rightSurfaceModel.valid());
    CHECK(output.leftSurfaceModel.patches.size()==expected.leftSurfaceModel.patches.size());
    Eigen::Vector3d fitted,reference;
    CHECK(output.weldSurfaceModel.evaluate(0,.5,.5,fitted));
    CHECK(expected.weldSurfaceModel.evaluate(0,.5,.5,reference));
    CHECK((fitted-reference).norm()<1e-10);
    CHECK(single ? !output.center->empty() : !output.left->empty() && !output.right->empty());
    CHECK(same(*input, unchanged));
    auto invalid = config; invalid.sliceStep = 0;
    CHECK(!weld::extractSeams(input, output, invalid, &error));
    CHECK(!error.empty() && output.left->empty() && output.right->empty() && output.center->empty());
    CHECK(output.leftSurface->empty() && output.rightSurface->empty() && output.weldSurface->empty());
    CHECK(!output.leftSurfaceModel.valid() && !output.rightSurfaceModel.valid() && !output.weldSurfaceModel.valid());
    auto badFrames = frames; ++badFrames[0].begin;
    CHECK(!weld::extractSeams(input, output, config, &error, &badFrames));
    CHECK(!error.empty());
  }
  weld::SeamClouds output;
  CHECK(!weld::extractSeams(weld::CloudT::ConstPtr(), output, config, &error));
  CHECK(!error.empty());
  // No frame metadata: the public API uses spatial sections of an ordinary cloud.
  weld::CloudT::Ptr spatial(new weld::CloudT);
  for (int x=0; x<80; ++x) for (int j=-250; j<=250; ++j) {
    const float y = .1f*j;
    spatial->push_back(weld::PointT(.8f*x, y, std::min(15.f, std::max(0.f, 2*(std::fabs(y)-8.f)))));
  }
  config = weld::AlgorithmParameters().get(weld::ProfileLocalizationMethod::CurveLineIntersection);
  config.extractMainBody = false; config.voxelSize = 0; config.outlierMeanK = 0;
  config.trimUnstableEndpoints = false;
  CHECK(weld::extractSeams(spatial, output, config, &error));
  CHECK(output.left->size() >= 3 && output.right->size() >= 3);
  CHECK(!output.leftSurface->empty() && !output.rightSurface->empty() && !output.weldSurface->empty());
  // Track-width conversion must not change physical surface segmentation.
  const auto physicalLeft=output.leftSurface, physicalRight=output.rightSurface, physicalWeld=output.weldSurface;
  config.uniformWidthEnabled=true;
  CHECK(weld::extractSeams(spatial, output, config, &error));
  CHECK(same(*physicalLeft,*output.leftSurface) && same(*physicalRight,*output.rightSurface) &&
        same(*physicalWeld,*output.weldSurface));
  // Output may own the next input. New surfaces must not clear or mutate it.
  config=weld::defaultDetectorConfig();
  output.weldSurface=spatial;
  const auto inputSize=spatial->size();
  CHECK(weld::extractSeams(output.weldSurface,output,config,&error));
  CHECK(spatial->size()==inputSize);
  CHECK(!output.leftSurface->empty() && !output.rightSurface->empty() && !output.weldSurface->empty());
  std::cout << "Cloud-only API: tracks, three surfaces, framed/unframed, single/dual, uniform width, invalid input and input ownership passed.\n";
  return 0;
}
