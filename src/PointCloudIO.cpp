#include "PointCloudIO.h"

#include <pcl/common/point_tests.h>
#ifndef WELD_SEAM_TEXT_IO_ONLY
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#endif

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

namespace weld {
namespace {

std::string lowerExtension(const std::string &path) {
  const std::size_t dot = path.find_last_of('.');
  if (dot == std::string::npos) {
    return std::string();
  }
  std::string result = path.substr(dot);
  std::transform(
      result.begin(), result.end(), result.begin(),
      [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
  return result;
}

void finalizeCloud(CloudT &cloud) {
  cloud.width = static_cast<std::uint32_t>(cloud.size());
  cloud.height = 1;
  cloud.is_dense = false;
}

} // namespace

CloudT::Ptr PointCloudIO::load(const std::string &path, std::string *error,
                               std::vector<ScanFrameRange> *frames) {
  if (frames)
    frames->clear();
  const std::string extension = lowerExtension(path);
  CloudT::Ptr cloud(new CloudT);

  if (extension == ".txt" || extension == ".xyz" || extension == ".csv") {
    return loadText(path, error, frames);
  }

  int status = -1;
#ifndef WELD_SEAM_TEXT_IO_ONLY
  if (extension == ".pcd") {
    status = pcl::io::loadPCDFile<PointT>(path, *cloud);
  } else if (extension == ".ply") {
    status = pcl::io::loadPLYFile<PointT>(path, *cloud);
  } else
#endif
  {
    if (error) {
      *error = "Unsupported point-cloud format: " + extension;
    }
    return CloudT::Ptr(new CloudT);
  }

  if (status < 0 || cloud->empty()) {
    if (error) {
      *error = "Unable to load point cloud: " + path;
    }
    cloud->clear();
  }
  return cloud;
}

bool PointCloudIO::savePointCloud(const std::string &path, const CloudT &cloud,
                                  std::string *error) {
  if (cloud.empty()) {
    if (error) {
      *error = "Cannot save an empty point cloud: " + path;
    }
    return false;
  }

  const std::string extension = lowerExtension(path);
  int status = -1;
#ifndef WELD_SEAM_TEXT_IO_ONLY
  if (extension == ".pcd") {
    status = pcl::io::savePCDFileBinary(path, cloud);
  } else if (extension == ".ply") {
    status = pcl::io::savePLYFileBinary(path, cloud);
  } else
#endif
  if (extension == ".xyz" || extension == ".txt") {
    std::ofstream output(path.c_str());
    if (output) {
      output << std::fixed << std::setprecision(6);
      for (const PointT &point : cloud) {
        output << point.x << ' ' << point.y << ' ' << point.z << '\n';
      }
      status = output.good() ? 0 : -1;
    }
  } else {
    if (error) {
      *error = "Unsupported output point-cloud format: " + extension;
    }
    return false;
  }

  if (status < 0) {
    if (error) {
      *error = "Unable to save point cloud: " + path;
    }
    return false;
  }
  return true;
}

CloudT::Ptr PointCloudIO::loadText(const std::string &path, std::string *error,
                                   std::vector<ScanFrameRange> *frames) {
  CloudT::Ptr cloud(new CloudT);
  std::ifstream input(path);
  std::vector<ScanFrameRange> ranges;
  if (frames)
    frames->clear();
  if (!input) {
    if (error)
      *error = "Unable to open point cloud: " + path;
    return cloud;
  }
  int xyz[3] = {0, 1, 2}, frameColumn = -1;
  bool resolved = false;
  std::size_t lineNumber = 0;
  auto fail = [&](const std::string &reason) {
    cloud->clear();
    if (error)
      *error = "Line " + std::to_string(lineNumber) + ": " + reason;
    return cloud;
  };
  std::string line;
  while (std::getline(input, line)) {
    ++lineNumber;
    if (line.compare(0, 3, "\xEF\xBB\xBF") == 0)
      line.erase(0, 3);
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream parser(line);
    std::vector<std::string> fields;
    std::string token;
    while (parser >> token)
      fields.push_back(token);
    if (fields.empty())
      continue;
    if (!resolved) {
      auto columns = fields;
      for (auto &column : columns)
        std::transform(column.begin(), column.end(), column.begin(),
                       [](unsigned char ch) {
                         return static_cast<char>(std::tolower(ch));
                       });
      const auto x = std::find(columns.begin(), columns.end(), "x"),
                 y = std::find(columns.begin(), columns.end(), "y"),
                 z = std::find(columns.begin(), columns.end(), "z");
      auto frame = std::find(columns.begin(), columns.end(), "frameindex");
      if (frame == columns.end())
        frame = std::find(columns.begin(), columns.end(), "frame_index");
      if (frame == columns.end())
        frame = std::find(columns.begin(), columns.end(), "target");
      if (frame != columns.end() &&
          (x == columns.end() || y == columns.end() || z == columns.end()))
        return fail("A frame header must include x, y and z.");
      if (x != columns.end() && y != columns.end() && z != columns.end()) {
        xyz[0] = int(x - columns.begin());
        xyz[1] = int(y - columns.begin());
        xyz[2] = int(z - columns.begin());
        if (frame != columns.end())
          frameColumn = int(frame - columns.begin());
        resolved = true;
        continue;
      }
    }
    PointT point;
    int frameId = 0;
    bool valid = fields.size() >
                 static_cast<std::size_t>(std::max(
                     frameColumn, std::max(xyz[0], std::max(xyz[1], xyz[2]))));
    try {
      if (valid)
        for (int axis = 0; axis < 3; ++axis) {
          std::size_t used = 0;
          const double value = std::stod(fields[xyz[axis]], &used);
          if (used != fields[xyz[axis]].size() || !std::isfinite(value))
            valid = false;
          point.data[axis] = static_cast<float>(value);
        }
      if (valid && frameColumn >= 0) {
        std::size_t used = 0;
        const long long value = std::stoll(fields[frameColumn], &used);
        if (used != fields[frameColumn].size() || value < 0 ||
            value > std::numeric_limits<int>::max())
          valid = false;
        else
          frameId = static_cast<int>(value);
      }
    } catch (const std::exception &) {
      valid = false;
    }
    if (!valid || !pcl::isFinite(point)) {
      if (frameColumn >= 0)
        return fail("Invalid XYZ or frameIndex; framed input was not partially "
                    "loaded.");
      continue;
    }
    resolved = true;
    if (frameColumn >= 0) {
      if (!ranges.empty() && frameId < ranges.back().frameIndex)
        return fail("frameIndex must be nondecreasing; interleaved frames are "
                    "not supported.");
      if (ranges.empty() || frameId != ranges.back().frameIndex)
        ranges.push_back({cloud->size(), 0, frameId});
      ++ranges.back().count;
    }
    cloud->push_back(point);
  }
  finalizeCloud(*cloud);
  if (frames)
    *frames = std::move(ranges);
  if (cloud->empty() && error)
    *error = "No valid XYZ points were found in " + path;
  return cloud;
}

bool PointCloudIO::saveFrameSeamCsv(const std::string &path,
                                    const SeamResult &result,
                                    std::string *error) {
  if (result.leftFrameIndices.size() != result.leftSeam->size() ||
      result.rightFrameIndices.size() != result.rightSeam->size() ||
      result.centerFrameIndices.size() != result.centerSeam->size()) {
    if (error)
      *error = "Frame indices do not match seam points.";
    return false;
  }
  std::ofstream output(path);
  if (!output) {
    if (error)
      *error = "Unable to write " + path;
    return false;
  }
  output << "side,frameIndex,x,y,z,confidence,segment\n"
         << std::fixed << std::setprecision(6);
  for (int side = 0; side < 3; ++side) {
    const auto &points = side==2 ? *result.centerSeam : side ? *result.rightSeam : *result.leftSeam;
    const auto &ids = side==2 ? result.centerFrameIndices : side ? result.rightFrameIndices : result.leftFrameIndices;
    const auto &scores =
        side==2 ? result.centerPointConfidence : side ? result.rightPointConfidence : result.leftPointConfidence;
    const auto &segments =
        side==2 ? result.centerSegmentStarts : side ? result.rightSegmentStarts : result.leftSegmentStarts;
    std::size_t segment = 0;
    for (std::size_t i = 0; i < points.size(); ++i) {
      if (segment + 1 < segments.size() && i == segments[segment + 1])
        ++segment;
      const auto &p = points[i];
      output << (side==2 ? 'C' : side ? 'R' : 'L') << ',' << ids[i] << ',' << p.x << ',' << p.y
             << ',' << p.z << ',' << scores[i] << ',' << segment << '\n';
    }
  }
  if (!output) {
    if (error)
      *error = "Writing frame CSV failed: " + path;
    return false;
  }
  return true;
}

bool PointCloudIO::saveResultCsv(const std::string &path, const SeamResult &result,
                                Axis axis, std::string *error) {
  if (result.frameBased) return saveFrameSeamCsv(path,result,error);
  std::ofstream output(path);
  if (!output) { if (error) *error="Unable to write "+path; return false; }
  output << "side,station,x,y,z,confidence,segment\n" << std::fixed << std::setprecision(6);
  for (int side=0;side<3;++side) {
    const auto &points=side==2 ? *result.centerSeam : side ? *result.rightSeam : *result.leftSeam;
    const auto &scores=side==2 ? result.centerPointConfidence : side ? result.rightPointConfidence : result.leftPointConfidence;
    const auto &starts=side==2 ? result.centerSegmentStarts : side ? result.rightSegmentStarts : result.leftSegmentStarts;
    std::size_t segment=0;
    for (std::size_t i=0;i<points.size();++i) {
      if (segment+1<starts.size() && i==starts[segment+1]) ++segment;
      const auto &p=points[i];
      output << (side==2 ? 'C' : side ? 'R' : 'L') << ',' << coordinate(p,axis) << ','
             << p.x << ',' << p.y << ',' << p.z << ',' << (i<scores.size()?scores[i]:0) << ',' << segment << '\n';
    }
  }
  if (!output) { if (error) *error="Writing seam CSV failed: "+path; return false; }
  return true;
}

bool PointCloudIO::saveSeamCsv(const std::string &path, const CloudT &seam,
                               const std::vector<float> &confidence, Axis axis,
                               std::string *error) {
  std::ofstream output(path.c_str());
  if (!output) {
    if (error) {
      *error = "Unable to create output file: " + path;
    }
    return false;
  }

  output << "station,x,y,z,confidence\n";
  output << std::fixed << std::setprecision(6);
  for (std::size_t i = 0; i < seam.size(); ++i) {
    const float score = i < confidence.size() ? confidence[i] : 0.0f;
    output << coordinate(seam[i], axis) << ',' << seam[i].x << ',' << seam[i].y
           << ',' << seam[i].z << ',' << score << '\n';
  }
  return true;
}

bool PointCloudIO::saveDualSeamCsv(const std::string &path,
                                   const CloudT &leftSeam,
                                   const std::vector<float> &leftConfidence,
                                   const CloudT &rightSeam,
                                   const std::vector<float> &rightConfidence,
                                   Axis axis, std::string *error) {
  std::ofstream output(path.c_str());
  if (!output) {
    if (error) {
      *error = "Unable to create output file: " + path;
    }
    return false;
  }

  output << "side,station,x,y,z,confidence\n";
  output << std::fixed << std::setprecision(6);
  const auto writeSide = [&](char side, const CloudT &seam,
                             const std::vector<float> &confidence) {
    for (std::size_t i = 0; i < seam.size(); ++i) {
      const float score = i < confidence.size() ? confidence[i] : 0.0f;
      output << side << ',' << coordinate(seam[i], axis) << ',' << seam[i].x
             << ',' << seam[i].y << ',' << seam[i].z << ',' << score << '\n';
    }
  };
  writeSide('L', leftSeam, leftConfidence);
  writeSide('R', rightSeam, rightConfidence);
  return true;
}

bool PointCloudIO::loadDualSeamCsv(const std::string &path,
                                   CloudT::Ptr *leftSeam,
                                   CloudT::Ptr *rightSeam, std::string *error) {
  if (!leftSeam || !rightSeam) {
    if (error) {
      *error = "Dual seam output pointers are null.";
    }
    return false;
  }
  std::ifstream input(path.c_str());
  if (!input) {
    if (error) {
      *error = "Unable to open dual seam file: " + path;
    }
    return false;
  }

  CloudT::Ptr left(new CloudT);
  CloudT::Ptr right(new CloudT);
  std::string line;
  while (std::getline(input, line)) {
    const std::size_t comma = line.find(',');
    if (comma == std::string::npos) {
      continue;
    }
    std::string side = line.substr(0, comma);
    std::transform(
        side.begin(), side.end(), side.begin(),
        [](unsigned char ch) { return static_cast<char>(std::toupper(ch)); });
    if (side != "L" && side != "LEFT" && side != "R" && side != "RIGHT") {
      continue;
    }
    std::string valuesText = line.substr(comma + 1);
    std::replace(valuesText.begin(), valuesText.end(), ',', ' ');
    std::istringstream parser(valuesText);
    double station = 0.0;
    PointT point;
    if (!(parser >> station >> point.x >> point.y >> point.z) ||
        !pcl::isFinite(point)) {
      continue;
    }
    (void)station;
    if (side == "L" || side == "LEFT") {
      left->push_back(point);
    } else {
      right->push_back(point);
    }
  }
  finalizeCloud(*left);
  finalizeCloud(*right);
  if (left->size() < 2 || right->size() < 2) {
    if (error) {
      *error = "The file does not contain both L and R seam tracks: " + path;
    }
    return false;
  }
  *leftSeam = left;
  *rightSeam = right;
  return true;
}

CloudT::Ptr PointCloudIO::loadSeamCsv(const std::string &path,
                                      std::string *error) {
  // 历史文件也只依赖 XYZ。saveSeamCsv 的首列是 station，所以这里兼容
  // “station,x,y,z,...”和普通“x,y,z”两种格式。
  std::ifstream input(path.c_str());
  CloudT::Ptr cloud(new CloudT);
  if (!input) {
    if (error) {
      *error = "Unable to open seam file: " + path;
    }
    return cloud;
  }

  std::string line;
  bool generatedFormat = false;
  bool formatKnown = false;
  while (std::getline(input, line)) {
    if (!formatKnown) {
      std::string lower = line;
      std::transform(
          lower.begin(), lower.end(), lower.begin(),
          [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
      if (lower.find("station") != std::string::npos) {
        generatedFormat = true;
        formatKnown = true;
        continue;
      }
    }
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream parser(line);
    std::vector<double> values;
    double value = 0.0;
    while (parser >> value) {
      values.push_back(value);
    }
    if (values.empty()) {
      // 表头无法解析成数字，会自然走到这里。
      continue;
    }

    PointT point;
    if (generatedFormat && values.size() >= 4) {
      point.x = static_cast<float>(values[1]);
      point.y = static_cast<float>(values[2]);
      point.z = static_cast<float>(values[3]);
    } else if (values.size() >= 3) {
      point.x = static_cast<float>(values[0]);
      point.y = static_cast<float>(values[1]);
      point.z = static_cast<float>(values[2]);
    } else {
      continue;
    }
    formatKnown = true;
    if (pcl::isFinite(point)) {
      cloud->push_back(point);
    }
  }
  finalizeCloud(*cloud);
  if (cloud->empty() && error) {
    *error = "No valid XYZ seam points were found in: " + path;
  }
  return cloud;
}

} // namespace weld
