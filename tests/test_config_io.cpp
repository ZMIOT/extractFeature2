#include "ConfigIO.h"
#include <fstream>
#include <iterator>
#include <cstdio>
#include <chrono>
#ifdef _WIN32
#include <windows.h>
#endif
#include <iostream>
#include <limits>

#define CHECK(x) do { if (!(x)) { std::cerr << "Failed line " << __LINE__ << ": " << #x << '\n'; return 1; } } while (0)
int main() {
  const auto prefix = std::string("config_test_") + std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  const auto path = prefix + ".json";
  weld::LibraryConfig original, loaded;
  auto a = original.algorithms.get(weld::ProfileLocalizationMethod::CurvatureCandidates);
  a.minSeamSeparation = 7.25f;
  a.curvatureSearch.maximumPairs = 13;
  a.curveLineFallback.search.continuityWeight = 4.5f;
  a.historyWeight = 8.5f;
  original.algorithms.set(weld::ProfileLocalizationMethod::CurvatureCandidates, a);
  std::string error;
  CHECK(weld::saveLibraryConfig(path, original, &error));
  CHECK(weld::loadLibraryConfig(path, loaded, &error));
  auto b = loaded.algorithms.get(loaded.activeAlgorithm);
  CHECK(b.minSeamSeparation == 7.25f && b.curvatureSearch.maximumPairs == 13);
  CHECK(b.curveLineFallback.search.continuityWeight == 4.5f && b.historyWeight == 8.5f);
  CHECK(loaded.algorithms.get(weld::ProfileLocalizationMethod::CurveLineIntersection).minSeamSeparation == 8);
  const auto second = prefix + "_second.json";
  CHECK(weld::saveLibraryConfig(second, loaded, &error));
  auto read = [](const std::string &name) {
    std::ifstream file(name, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>());
  };
  CHECK(read(path) == read(second));
  auto write = [&](const std::string &text) {
    std::ofstream file(path, std::ios::binary | std::ios::trunc);
    file << text;
    file.close();
    return bool(file);
  };
  for (const auto &bad : {
       "{", "[]", "{\"schemaVersion\":2}", "{\"schemaVersion\":1,\"unknown\":1}",
       "{\"schemaVersion\":1,\"algorithms\":{\"curvature\":{\"sliceStep\":0}}}",
       "{\"schemaVersion\":1,\"algorithms\":{\"curvature\":{\"profileMergeFrameCount\":4}}}",
       "{\"schemaVersion\":1,\"algorithms\":{\"curvature\":{\"outlierMeanK\":2.5}}}",
       "{\"schemaVersion\":1,\"algorithms\":{\"curvature\":{\"debugFlag\":1}}}",
       "{\"schemaVersion\":1,\"algorithms\":{\"curvature\":{\"profileLocalizationMethod\":0}}}",
       "{\"schemaVersion\":1,\"algorithms\":{\"curvature\":{\"curvatureSearch\":{\"typo\":1}}}}"}) {
    CHECK(write(bad)); CHECK(!weld::loadLibraryConfig(path, loaded, &error));
    CHECK(!error.empty());
    CHECK(loaded.algorithms.get(loaded.activeAlgorithm).minSeamSeparation == 7.25f);
  }
  CHECK(write("{\"schemaVersion\":1}"));
  CHECK(weld::loadLibraryConfig(path, loaded, &error)); CHECK(error.empty());
  CHECK(loaded.algorithms.get(loaded.activeAlgorithm).minSeamSeparation == 6);
  // The shipped file must remain equivalent to the library's built-in defaults.
  weld::LibraryConfig shipped;
  CHECK(weld::loadLibraryConfig(SEAM_DEFAULT_CONFIG, shipped, &error));
  CHECK(weld::saveLibraryConfig(path, shipped, &error));
  CHECK(weld::saveLibraryConfig(second, loaded, &error));
  CHECK(read(path) == read(second));
  // UTF-8 paths work without Qt; save validation must not overwrite a valid file.
  const auto before = read(path);
  auto invalid = loaded;
  auto invalidParams = invalid.algorithms.get(invalid.activeAlgorithm);
  invalidParams.sliceStep = -1;
  invalid.algorithms.set(invalid.activeAlgorithm, invalidParams);
  CHECK(!weld::saveLibraryConfig(path, invalid, &error));
  CHECK(read(path) == before);
#ifdef _WIN32
  const auto unicodePath = prefix + u8"_参数.json";
  CHECK(weld::saveLibraryConfig(unicodePath, loaded, &error));
  CHECK(weld::loadLibraryConfig(unicodePath, shipped, &error));
  std::wstring wide(prefix.begin(), prefix.end());
  wide += L"_参数.json";
  CHECK(DeleteFileW(wide.c_str()));
#endif
  a.sliceStep = std::numeric_limits<float>::quiet_NaN();
  CHECK(!weld::validateDetectorConfig(a, &error));
  CHECK(!weld::loadLibraryConfig(path + ".missing", loaded, &error));
  CHECK(!weld::saveLibraryConfig(path + "/impossible", loaded, &error));
  std::remove(path.c_str());
  std::remove(second.c_str());
  return 0;
}
