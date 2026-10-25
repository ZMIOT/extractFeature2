#include "ConfigIO.h"
#include <nlohmann/json.hpp>
#include <cmath>
#include <type_traits>
#include <fstream>
#include <sstream>
#include <set>
#include <atomic>
#include <cstdio>
#include <stdexcept>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace weld {
namespace {
using Json = nlohmann::json;
// One codec validates both file input and in-memory values before saving.
struct Codec {
  Json object;
  bool reading;
  std::string &error;
  std::string prefix;
  std::set<std::string> known;

  void fail(const std::string &key, const std::string &reason) {
    if (error.empty()) error = prefix + key + ": " + reason;
  }
  void boolean(const char *key, bool &value) {
    known.insert(key);
    if (reading && object.count(key)) {
      if (!object[key].is_boolean()) { fail(key, "expected boolean"); return; }
      value = object[key].get<bool>();
    } else if (!reading) object[key] = value;
  }
  template<class T> void number(const char *key, T &value, double low, double high) {
    known.insert(key);
    double n = static_cast<double>(value);
    if (reading && object.count(key)) {
      if (!object[key].is_number()) { fail(key, "expected number"); return; }
      n = object[key].get<double>();
    }
    if (!std::isfinite(n) || n < low || n > high ||
        ((std::is_integral<T>::value || std::is_enum<T>::value) && n != std::floor(n))) {
      fail(key, "expected " + std::string(std::is_integral<T>::value || std::is_enum<T>::value ? "integer" : "number")
          + " in [" + std::to_string(low) + ", " + std::to_string(high) + "]");
      return;
    }
    value = static_cast<T>(n);
    if (!reading) object[key] = n;
  }
  template<class T, class F> void nested(const char *key, T &value, F visit) {
    known.insert(key);
    if (reading && object.count(key) && !object[key].is_object()) {
      fail(key, "expected object"); return;
    }
    Codec child{object.count(key) ? object[key] : Json::object(), reading, error, prefix + key + ".", {}};
    visit(child, value);
    child.finish();
    if (!reading) object[key] = child.object;
  }
  void finish() {
    if (reading) for (auto it = object.begin(); it != object.end(); ++it)
      if (!known.count(it.key())) fail(it.key(), "unknown configuration key");
  }
};
void visitSearch(Codec &v, CurvatureSearchConfig &s) {
  v.boolean("multiCandidates", s.multiCandidates);
  v.number("windowLength", s.windowLength, 1e-6, 1e9);
  v.number("valleyProminence", s.valleyProminence, 1e-6, 1e9);
  v.number("candidateSpacing", s.candidateSpacing, 1e-6, 1e9);
  v.number("maximumPairs", s.maximumPairs, 2, 32);
  v.number("continuityWeight", s.continuityWeight, 0, 1e9);
}
void visitFallback(Codec &v, CurveLineFallbackConfig &c) {
  v.boolean("enabled", c.enabled);
  v.number("splineRadius", c.splineRadius, 1e-6, 1e9);
  v.number("wallTolerance", c.wallTolerance, 1e-6, 1e9);
  v.number("regionLength", c.regionLength, 1e-6, 1e9);
  v.number("regionFraction", c.regionFraction, 1e-6, 1);
  v.nested("search", c.search, visitSearch);
}
void visitDetector(Codec &v, DetectorConfig &c) {
#define CONFIG_BOOL(name) v.boolean(#name, c.name);
#define CONFIG_NUMBER(name, low, high) v.number(#name, c.name, low, high);
#include "ConfigFields.inc"
#undef CONFIG_BOOL
#undef CONFIG_NUMBER
  v.nested("curvatureSearch", c.curvatureSearch, visitSearch);
  v.nested("curveLineFallback", c.curveLineFallback, visitFallback);
  if (c.profileMergeFrameCount % 2 == 0)
    v.fail("profileMergeFrameCount", "expected odd integer from 1 to 101");
}
Json transfer(Json object, bool reading, LibraryConfig &c, std::string &error) {
  Codec root{object, reading, error, "", {"schemaVersion", "activeAlgorithm", "algorithms"}};
  int version = 1;
  if (reading && !object.count("schemaVersion")) root.fail("schemaVersion", "required");
  root.number("schemaVersion", version, 1, 1);
  std::string method = c.activeAlgorithm == ProfileLocalizationMethod::CurveLineIntersection ? "curve_line" : "curvature";
  if (reading && object.count("activeAlgorithm")) {
    if (!object["activeAlgorithm"].is_string()) root.fail("activeAlgorithm", "expected string");
    else method = object["activeAlgorithm"].get<std::string>();
  }
  if (method != "curve_line" && method != "curvature") root.fail("activeAlgorithm", "expected curve_line or curvature");
  if (c.activeAlgorithm != ProfileLocalizationMethod::CurveLineIntersection &&
      c.activeAlgorithm != ProfileLocalizationMethod::CurvatureCandidates)
    root.fail("activeAlgorithm", "invalid enum");
  c.activeAlgorithm = method == "curve_line" ? ProfileLocalizationMethod::CurveLineIntersection : ProfileLocalizationMethod::CurvatureCandidates;
  if (reading && object.count("algorithms") && !object["algorithms"].is_object()) root.fail("algorithms", "expected object");
  Codec algorithms{object.count("algorithms") && object["algorithms"].is_object() ? object["algorithms"] : Json::object(),
                   reading, error, "algorithms.", {}};
  for (int i = 0; i < 2; ++i) {
    const auto selected = static_cast<ProfileLocalizationMethod>(i);
    auto snapshot = c.algorithms.get(selected);
    const char *key = i == 0 ? "curve_line" : "curvature";
    algorithms.nested(key, snapshot, visitDetector);
    if (snapshot.profileLocalizationMethod != selected) algorithms.fail(key, "profileLocalizationMethod conflicts with algorithm group");
    c.algorithms.set(selected, snapshot);
  }
  algorithms.finish();
  root.finish();
  if (!reading) {
    root.object["activeAlgorithm"] = method;
    root.object["algorithms"] = algorithms.object;
  }
  return root.object;
}
} // namespace

bool validateDetectorConfig(const DetectorConfig &config, std::string *error) {
  std::string detail;
  auto copy = config;
  Codec codec{Json::object(), false, detail, "", {}};
  visitDetector(codec, copy);
  if (error) *error = detail;
  return detail.empty();
}
namespace {
#ifdef _WIN32
std::wstring widePath(const std::string &path) {
  if (path.find('\0') != std::string::npos) throw std::runtime_error("Path contains a null character");
  const int count = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path.c_str(), -1, nullptr, 0);
  if (!count) throw std::runtime_error("Path is not valid UTF-8");
  std::wstring result(static_cast<std::size_t>(count), L'\0');
  MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path.c_str(), -1, &result[0], count);
  result.pop_back();
  return result;
}
#endif
void atomicWrite(const std::string &path, const std::string &bytes) {
#ifdef _WIN32
  const auto destination = widePath(path);
  static std::atomic<unsigned long> sequence{0};
  std::wstring temporary;
  HANDLE file = INVALID_HANDLE_VALUE;
  for (int attempt = 0; attempt < 100 && file == INVALID_HANDLE_VALUE; ++attempt) {
    temporary = destination + L".tmp." + std::to_wstring(GetCurrentProcessId()) + L"." + std::to_wstring(sequence++);
    file = CreateFileW(temporary.c_str(), GENERIC_WRITE, 0, nullptr, CREATE_NEW, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (file == INVALID_HANDLE_VALUE && GetLastError() != ERROR_FILE_EXISTS)
      throw std::runtime_error("Cannot create configuration temporary file: " + path);
  }
  if (file == INVALID_HANDLE_VALUE) throw std::runtime_error("Cannot reserve configuration temporary file");
  DWORD written = 0;
  const bool saved = WriteFile(file, bytes.data(), static_cast<DWORD>(bytes.size()), &written, nullptr)
      && written == bytes.size() && FlushFileBuffers(file);
  CloseHandle(file);
  if (!saved || !MoveFileExW(temporary.c_str(), destination.c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH)) {
    DeleteFileW(temporary.c_str());
    throw std::runtime_error("Cannot write or replace configuration: " + path);
  }
#else
  std::string name = path + ".tmp.XXXXXX";
  std::vector<char> temporary(name.begin(), name.end());
  temporary.push_back('\0');
  const int descriptor = mkstemp(temporary.data());
  if (descriptor < 0) throw std::runtime_error("Cannot create configuration temporary file: " + path);
  FILE *file = fdopen(descriptor, "wb");
  if (!file) { close(descriptor); std::remove(temporary.data()); throw std::runtime_error("Cannot open temporary file"); }
  bool saved = std::fwrite(bytes.data(), 1, bytes.size(), file) == bytes.size();
  saved = std::fflush(file) == 0 && saved;
  saved = fsync(descriptor) == 0 && saved;
  saved = std::fclose(file) == 0 && saved;
  if (!saved || std::rename(temporary.data(), path.c_str()) != 0) {
    std::remove(temporary.data());
    throw std::runtime_error("Cannot write or replace configuration: " + path);
  }
#endif
}
} // namespace

bool loadLibraryConfig(const std::string &path, LibraryConfig &config, std::string *error) {
  std::string detail;
  try {
#ifdef _WIN32
    std::ifstream file(widePath(path).c_str(), std::ios::binary);
#else
    if (path.find('\0') != std::string::npos) throw std::runtime_error("Path contains a null character");
    std::ifstream file(path, std::ios::binary);
#endif
    if (!file) throw std::runtime_error("Cannot open configuration: " + path);
    const auto object = Json::parse(file);
    if (!object.is_object()) throw std::runtime_error("Configuration root must be an object");
    LibraryConfig next;
    transfer(object, true, next, detail);
    if (detail.empty()) config = next;
  } catch (const std::exception &e) { detail = e.what(); }
  if (error) *error = detail;
  return detail.empty();
}
bool saveLibraryConfig(const std::string &path, const LibraryConfig &config, std::string *error) {
  auto copy = config;
  std::string detail;
  try {
    const auto object = transfer(Json::object(), false, copy, detail);
    if (detail.empty()) atomicWrite(path, object.dump(4) + "\n");
  } catch (const std::exception &e) { detail = e.what(); }
  if (error) *error = detail;
  return detail.empty();
}
} // namespace weld
