#pragma once

#include "sirius/backend/retained_integrator.h"
#include "sirius/backend/trace_step_executor.h"

#include <condition_variable>
#include <deque>
#include <functional>
#include <map>
#include <mutex>
#include <thread>

namespace sirius::backend {

// Coalesces synchronous trace workers into bounded device batches. A worker's
// retained phase survives accepted intervals and a rollback after event checks.
// Public binary64 trace fields are views, not the next device phase record.
class RetainedTraceExecutor final : public TraceStepExecutor {
  public:
    explicit RetainedTraceExecutor(RetainedCompute& compute,
                                   std::function<bool()> should_cancel = {},
                                   double maximum_submission_ms = 0);
    ~RetainedTraceExecutor() override;
    void BeginTrace() override;
    void EndTrace() override;
    void RejectLastInterval() override;
    std::optional<core::CameraLaunch> Launch(core::IMetric& metric, double absolute_spin,
                                             const core::CameraRay& camera) override;
    bool Step(core::Lightray& ray, core::IMetric& metric, const core::IntegratorConfig& config,
              core::Rk45CoupledState& coupled, core::Rk45CoupledComparison& comparison) override;
    [[nodiscard]] std::optional<base::Error> Error() const;
    struct Stats {
        std::uint64_t interval_batches = 0;
        std::uint64_t camera_batches = 0;
        std::uint64_t batch_subdivisions = 0;
        std::uint64_t safety_fallbacks = 0;
        std::size_t maximum_batch_rows = 0;
        std::uint64_t reused_phases = 0;
        std::uint64_t accepted_intervals = 0;
        std::uint64_t rejected_intervals = 0;
    };
    [[nodiscard]] Stats Statistics() const;

  private:
    struct Snapshot {
        std::array<double, 40> physical{};
        std::array<double, 4> metric{};
        double chart = 0;
        float affine = 0;
        bool operator==(const Snapshot&) const = default;
    };
    struct Continuation {
        Snapshot before, after;
        RetainedEndpointOutput start, finish;
    };
    struct Request {
        bool camera = false;
        RetainedRayCameraInput camera_input;
        base::Expected<RetainedCameraOutput> camera_result = RetainedCameraOutput{};
        RetainedInitializeInput initialization;
        RetainedIntervalInput interval;
        base::Expected<RetainedIntervalOutput> result = RetainedIntervalOutput{};
        bool completed = false;
    };
    void Run();
    void Execute(std::span<Request*> requests);
    RetainedCompute& compute_;
    double maximum_submission_ms_;
    std::function<bool()> should_cancel_;
    mutable std::mutex mutex_;
    std::condition_variable available_, completed_;
    std::deque<Request*> requests_;
    std::map<std::thread::id, Continuation> continuations_;
    std::optional<base::Error> error_;
    Stats stats_;
    bool stopping_ = false;
    std::thread dispatcher_;
};

}  // namespace sirius::backend
