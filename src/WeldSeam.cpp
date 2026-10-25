#include "WeldSeam.h"
#include "SeamDetector.h"
#include <exception>
#include <new>

namespace weld {
DetectorConfig defaultDetectorConfig() {
  return AlgorithmParameters().get(ProfileLocalizationMethod::CurvatureCandidates);
}
bool extractSeams(const CloudT::ConstPtr &input, SeamClouds &output,
                  const DetectorConfig &config, std::string *error,
                  const std::vector<ScanFrameRange> *frames) {
  // Do not mutate pointed-to clouds: an output from an earlier call may also
  // be the input to this call. Read it completely before replacing output.
  std::string detail;
  if (validateDetectorConfig(config, &detail)) {
    try {
      const auto result = SeamDetector().detect(input, config, nullptr, nullptr, nullptr, frames);
      if (result.ok()) {
        output.left = result.leftSeam;
        output.right = result.rightSeam;
        output.center = result.centerSeam;
        output.leftSurface = result.leftSurface;
        output.rightSurface = result.rightSurface;
        output.weldSurface = result.weldSurface;
        output.leftSurfaceModel = result.leftSurfaceModel;
        output.weldSurfaceModel = result.weldSurfaceModel;
        output.rightSurfaceModel = result.rightSurfaceModel;
        if (error) error->clear();
        return true;
      }
      detail = result.summary;
    } catch (const std::bad_alloc &) { throw; }
      catch (const std::exception &e) { detail = e.what(); }
  }
  output = SeamClouds{};
  if (error) *error = detail;
  return false;
}
} // namespace weld
