#pragma once

#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/base/error.h"

#include <condition_variable>
#include <functional>
#include <memory>
#include <mutex>
#include <span>
#include <thread>
#include <vector>

namespace sirius::render {

// Independent camera probes share a bounded set of trace workers. Each worker
// exclusively owns its callback index (and thus its tracer); callers receive
// results in request order, regardless of execution order. Cancellation belongs
// to the trace callback, including work that was queued when cancellation began.
// As with other owned executors, callers must finish Execute before destruction.
class RayWorkQueue {
  public:
    using TraceFunction = std::function<backend::TraceResult(std::size_t, const core::CameraRay&)>;
    static constexpr std::size_t kMaximumBatch = 64;

    RayWorkQueue(std::size_t workers, TraceFunction trace);
    ~RayWorkQueue();
    RayWorkQueue(const RayWorkQueue&) = delete;
    RayWorkQueue& operator=(const RayWorkQueue&) = delete;

    [[nodiscard]] base::Expected<std::vector<backend::TraceResult>> Execute(
        std::span<const core::CameraRay> rays);

  private:
    struct Batch;
    struct Job {
        std::shared_ptr<Batch> batch;
        std::size_t index = 0;
    };
    void Run(std::size_t worker);
    void Stop();

    TraceFunction trace_;
    std::mutex mutex_;
    std::condition_variable available_, space_;
    std::vector<Job> queue_;
    std::size_t head_ = 0, size_ = 0;
    bool stopping_ = false;
    std::vector<std::thread> workers_;
};

}  // namespace sirius::render
