#pragma once

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdlib>
#include <exception>
#include <functional>
#include <mutex>
#include <system_error>
#include <thread>
#include <vector>

namespace weld {
namespace detail {

// Per detection-call budget; configure the environment before starting calls.
// No public configuration/ABI change and no process-wide Eigen setting change.
inline unsigned candidateWorkerCount() {
  const unsigned available = std::max(1u, std::thread::hardware_concurrency());
  unsigned requested = 4;
  if (const char *value = std::getenv("WELD_SEAM_CANDIDATE_THREADS")) {
    char *end = nullptr;
    const long parsed = std::strtol(value, &end, 10);
    if (end != value && *end == '\0' && parsed >= 1 && parsed <= 8)
      requested = static_cast<unsigned>(parsed);
  }
  return std::min(requested, available);
}

// One caller plus a fixed set of reusable workers. Owned by one section sweep,
// never shared across concurrent detect() calls. run() is synchronous: every
// worker finishes before captured frame data is destroyed or an error escapes.
class CandidateExecutor {
public:
  explicit CandidateExecutor(unsigned workers) : requested_(std::max(1u, workers)) {}
  CandidateExecutor(const CandidateExecutor &) = delete;
  CandidateExecutor &operator=(const CandidateExecutor &) = delete;
  ~CandidateExecutor() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stopping_ = true;
    }
    ready_.notify_all();
    for (auto &worker : workers_) worker.join();
  }

  bool parallelFor(std::size_t count) const { return requested_ > 1 && count >= 16; }

  void run(std::size_t count, const std::function<void(std::size_t)> &work) {
    if (!parallelFor(count)) {
      for (std::size_t i = 0; i < count; ++i) work(i);
      return;
    }
    if (!started_) {
      workers_.reserve(requested_ - 1);
      started_ = true;
      for (unsigned i = 1; i < requested_; ++i) {
        try { workers_.emplace_back([this] { workerLoop(); }); }
        catch (const std::system_error &) { break; } // Use available workers.
      }
    }
    {
      std::lock_guard<std::mutex> lock(mutex_);
      task_ = work;
      count_ = count;
      next_.store(0, std::memory_order_relaxed);
      failed_.store(false, std::memory_order_relaxed);
      failure_ = nullptr;
      remaining_ = workers_.size();
      ++generation_;
    }
    ready_.notify_all();
    executeChunks(); // The calling thread participates too.
    std::unique_lock<std::mutex> lock(mutex_);
    done_.wait(lock, [this] { return remaining_ == 0; });
    task_ = {};
    const auto failure = failure_;
    lock.unlock();
    if (failure) std::rethrow_exception(failure);
  }

private:
  void executeChunks() noexcept {
    try {
      constexpr std::size_t grain = 4;
      while (!failed_.load(std::memory_order_relaxed)) {
        const auto begin = next_.fetch_add(grain, std::memory_order_relaxed);
        if (begin >= count_) break;
        for (auto i = begin; i < std::min(count_, begin + grain); ++i)
          task_(i);
      }
    } catch (...) {
      std::lock_guard<std::mutex> lock(mutex_);
      if (!failure_) failure_ = std::current_exception();
      failed_.store(true, std::memory_order_relaxed);
    }
  }

  void workerLoop() {
    std::size_t seen = 0;
    std::unique_lock<std::mutex> lock(mutex_);
    for (;;) {
      ready_.wait(lock, [&] { return stopping_ || generation_ != seen; });
      if (stopping_) return;
      seen = generation_;
      lock.unlock();
      executeChunks();
      lock.lock();
      if (--remaining_ == 0) done_.notify_one();
    }
  }

  const unsigned requested_;
  bool started_ = false, stopping_ = false;
  std::vector<std::thread> workers_;
  std::mutex mutex_;
  std::condition_variable ready_, done_;
  std::size_t generation_ = 0, remaining_ = 0, count_ = 0;
  std::atomic<std::size_t> next_{0};
  std::atomic<bool> failed_{false};
  std::function<void(std::size_t)> task_;
  std::exception_ptr failure_;
};
}
}
