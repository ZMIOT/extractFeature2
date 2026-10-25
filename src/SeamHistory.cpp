#include "SeamHistory.h"

#include <algorithm>
#include <cmath>

namespace weld {

SeamHistory::SeamHistory(std::size_t capacity)
    : capacity_(std::max<std::size_t>(1, capacity)) {}

void SeamHistory::setCapacity(std::size_t capacity) {
    capacity_ = std::max<std::size_t>(1, capacity);
    while (passes_.size() > capacity_) {
        passes_.pop_front();
    }
}

void SeamHistory::add(const CloudT::ConstPtr& seam, Axis axis) {
    if (!seam || seam->size() < 2) {
        return;
    }
    CloudT::Ptr copy(new CloudT(*seam));
    std::sort(copy->begin(), copy->end(),
              [axis](const PointT& lhs, const PointT& rhs) {
                  return coordinate(lhs, axis) < coordinate(rhs, axis);
              });
    copy->width = static_cast<std::uint32_t>(copy->size());
    copy->height = 1;
    passes_.push_back(copy);
    while (passes_.size() > capacity_) {
        passes_.pop_front();
    }
}

void SeamHistory::clear() {
    passes_.clear();
}

bool SeamHistory::interpolate(const CloudT& seam,
                              float station,
                              Axis axis,
                              Eigen::Vector3f* point) {
    if (seam.size() < 2 || !point) {
        return false;
    }
    const float begin = coordinate(seam.front(), axis);
    const float end = coordinate(seam.back(), axis);
    if (station < begin || station > end) {
        return false;
    }

    const auto upper = std::lower_bound(
        seam.begin(), seam.end(), station,
        [axis](const PointT& sample, float value) {
            return coordinate(sample, axis) < value;
        });
    if (upper == seam.begin()) {
        *point = toVector(*upper);
        return true;
    }
    if (upper == seam.end()) {
        *point = toVector(seam.back());
        return true;
    }

    const PointT& right = *upper;
    const PointT& left = *(upper - 1);
    const float leftStation = coordinate(left, axis);
    const float span = coordinate(right, axis) - leftStation;
    const float t = std::fabs(span) > 1.0e-6f ? (station - leftStation) / span : 0.0f;
    *point = (1.0f - t) * toVector(left) + t * toVector(right);
    return true;
}

HistoryPrediction SeamHistory::predict(float station,
                                       Axis axis,
                                       float maxExtrapolation,
                                       int inwardAxis,
                                       int inwardDirection) const {
    HistoryPrediction result;
    std::vector<Eigen::Vector3f> samples;
    for (const CloudT::Ptr& pass : passes_) {
        Eigen::Vector3f point;
        if (pass && interpolate(*pass, station, axis, &point)) {
            samples.push_back(point);
        }
    }
    if (samples.empty()) {
        return result;
    }

    // 只有一道历史时直接沿用；有多道时，用最近两道的位移作有限外推。
    // 这对应多层焊逐层上移的物理过程，同时限幅可防止一次误检不断放大。
    Eigen::Vector3f predicted = samples.back();
    if (samples.size() >= 2) {
        Eigen::Vector3f delta = samples.back() - samples[samples.size() - 2];
        delta[static_cast<int>(axis)] = 0.0f;
        // 船体焊接变形只允许两侧向中间收缩。左侧传 +1，阻止 Y 负向外扩；
        // 右侧传 -1，阻止 Y 正向外扩。高度方向仍可随焊层正常变化。
        if (inwardAxis >= 0 && inwardAxis < 3) {
            if ((inwardDirection > 0 && delta[inwardAxis] < 0.0f) ||
                (inwardDirection < 0 && delta[inwardAxis] > 0.0f)) {
                delta[inwardAxis] = 0.0f;
            }
        }
        const float norm = delta.norm();
        if (norm > maxExtrapolation && norm > 1.0e-6f) {
            delta *= maxExtrapolation / norm;
        }
        predicted += delta;
    }
    predicted[static_cast<int>(axis)] = station;

    float spread = 0.0f;
    for (const Eigen::Vector3f& sample : samples) {
        Eigen::Vector3f difference = sample - samples.back();
        difference[static_cast<int>(axis)] = 0.0f;
        spread += difference.norm();
    }
    spread /= static_cast<float>(samples.size());

    result.valid = true;
    result.point = predicted;
    result.stability = 1.0f / (1.0f + spread / std::max(1.0f, maxExtrapolation));
    return result;
}

} // namespace weld
