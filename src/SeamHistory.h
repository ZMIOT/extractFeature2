#pragma once

#include "SeamTypes.h"

#include <deque>
#include <vector>

namespace weld {

class SeamHistory {
public:
    /** @brief 创建历史容器。 @param[in] capacity 最多保留轨迹条数，0 按 1 处理。 */
    explicit SeamHistory(std::size_t capacity = 5);

    /** @brief 修改容量并淘汰最早轨迹。 @param[in] capacity 最小按 1 处理。 */
    void setCapacity(std::size_t capacity);
    /** @brief 深拷贝轨迹并按纵向坐标排序；空指针或不足 2 点忽略。
     * @param[in] seam 有限坐标的历史 XYZ 轨迹。 @param[in] axis 纵向坐标轴。 */
    void add(const CloudT::ConstPtr& seam, Axis axis);
    /** @brief 清空历史轨迹，保持容量。 */
    void clear();

    /** @return 是否没有已保存的轨迹。 */
    bool empty() const { return passes_.empty(); }
    /** @return 当前历史轨迹条数。 */
    std::size_t size() const { return passes_.size(); }
    /** @return 从旧到新的轨迹容器引用；调用者不得修改其中点云。 */
    const std::deque<CloudT::Ptr>& passes() const { return passes_; }

    /** @brief 插值并基于最近两层历史估计下一层位置。
     * @param[in] station 纵向站位坐标，与点云同单位。
     * @param[in] axis 与 add 时相同的纵向轴。
     * @param[in] maxExtrapolation 非负的层间外推位移上限（长度）。
     * @param[in] inwardAxis 限制外扩的轴索引 0..2；-1 不限制。
     * @param[in] inwardDirection 1 仅正向、-1 仅负向、0 不限制。
     * @return valid 为 false 表示该站位没有历史支持，否则含预测点及稳定度。
     */
    HistoryPrediction predict(float station,
                              Axis axis,
                              float maxExtrapolation,
                              int inwardAxis = -1,
                              int inwardDirection = 0) const;

private:
    static bool interpolate(const CloudT& seam,
                            float station,
                            Axis axis,
                            Eigen::Vector3f* point);

    std::size_t capacity_;
    std::deque<CloudT::Ptr> passes_;
};

} // namespace weld
