#include "sirius/render/session/ray_work_queue.h"

#include <exception>
#include <stdexcept>
#include <utility>

namespace sirius::render {

struct RayWorkQueue::Batch {
    explicit Batch(std::span<const core::CameraRay> input)
        : rays(input.begin(), input.end()), results(input.size()) {}
    std::vector<core::CameraRay> rays;
    std::vector<backend::TraceResult> results;
    std::mutex mutex;
    std::condition_variable completed;
    std::size_t remaining = 0;
    bool failed = false;
};

RayWorkQueue::RayWorkQueue(std::size_t workers, TraceFunction trace) : trace_(std::move(trace)) {
    if (workers == 0 || workers > 1024 || !trace_)
        throw std::invalid_argument("ray queue requires 1..1024 workers and a trace callback");
    queue_.resize(2 * workers);
    workers_.reserve(workers);
    try {
        for (std::size_t i = 0; i < workers; ++i) workers_.emplace_back([this, i] { Run(i); });
    } catch (...) {
        Stop();
        throw;
    }
}

RayWorkQueue::~RayWorkQueue() { Stop(); }

void RayWorkQueue::Stop() {
    {
        std::lock_guard lock(mutex_);
        stopping_ = true;
    }
    available_.notify_all();
    space_.notify_all();
    for (auto& worker : workers_)
        if (worker.joinable()) worker.join();
}

base::Expected<std::vector<backend::TraceResult>> RayWorkQueue::Execute(
    std::span<const core::CameraRay> rays) {
    if (rays.size() > kMaximumBatch)
        return base::Fail(base::ErrorDomain::kConfiguration, "trace probe batch",
                          "batch exceeds the bounded ray capacity");
    if (rays.empty()) return std::vector<backend::TraceResult>{};
    try {
        const auto batch = std::make_shared<Batch>(rays);
        batch->remaining = rays.size();
        for (std::size_t i = 0; i < rays.size(); ++i) {
            std::unique_lock lock(mutex_);
            space_.wait(lock, [&] { return stopping_ || size_ < queue_.size(); });
            if (stopping_)
                return base::Fail(base::ErrorDomain::kInternal, "trace probe batch",
                                  "ray queue is stopping");
            queue_[(head_ + size_) % queue_.size()] = {batch, i};
            ++size_;
            lock.unlock();
            available_.notify_one();
        }
        std::unique_lock lock(batch->mutex);
        batch->completed.wait(lock, [&] { return batch->remaining == 0; });
        if (batch->failed)
            return base::Fail(base::ErrorDomain::kInternal, "trace probe batch",
                              "a ray worker threw while tracing");
        return std::move(batch->results);
    } catch (const std::exception& error) {
        return base::Fail(base::ErrorDomain::kInternal, "trace probe batch", error.what());
    }
}

void RayWorkQueue::Run(std::size_t worker) {
    for (;;) {
        Job job;
        {
            std::unique_lock lock(mutex_);
            available_.wait(lock, [&] { return stopping_ || size_ > 0; });
            if (size_ == 0) return;
            job = std::move(queue_[head_]);
            head_ = (head_ + 1) % queue_.size();
            --size_;
        }
        space_.notify_one();
        bool failed = false;
        try {
            job.batch->results[job.index] = trace_(worker, job.batch->rays[job.index]);
        } catch (...) {
            failed = true;
        }
        {
            std::lock_guard lock(job.batch->mutex);
            job.batch->failed = job.batch->failed || failed;
            --job.batch->remaining;
        }
        job.batch->completed.notify_one();
    }
}

}  // namespace sirius::render
