#pragma once

#include "sirius/backend/retained_integrator.h"
#include "sirius/backend/trace_step_executor.h"

#include <condition_variable>
#include <deque>
#include <functional>
#include <map>
#include <mutex>
#include <thread>
#include <vector>

namespace sirius::backend {

// Coalesces synchronous trace workers into bounded device batches. A worker's
// retained phase survives accepted intervals and a rollback after event checks.
// Public binary64 trace fields are views, not the next device phase record.
class RetainedTraceExecutor final : public TraceStepExecutor {
  public:
    explicit RetainedTraceExecutor(RetainedCompute& compute,
                                   std::function<bool()> should_cancel = {},
                                   double maximum_submission_ms = 0,
                                   std::size_t maximum_batch_rows = 0);
    ~RetainedTraceExecutor() override;
    void BeginTrace() override;
    void EndTrace() override;
    void RejectLastInterval() override;
    std::optional<core::CameraLaunch> Launch(core::IMetric& metric, double absolute_spin,
                                             const core::CameraRay& camera) override;
    bool Step(core::Lightray& ray, core::IMetric& metric, const core::IntegratorConfig& config,
              core::Rk45CoupledState& coupled, core::Rk45CoupledComparison& comparison) override;
    // A submission safety error includes its observed one-row stage and host
    // submit/wait duration; it does not attribute the driver's compilation work.
    // Sticky errors complete queued calls without further device work.
    [[nodiscard]] std::optional<base::Error> Error() const;
    struct Stats {
        std::uint64_t interval_batches = 0;
        std::uint64_t camera_batches = 0;
        std::uint64_t batch_subdivisions = 0;
        std::uint64_t safety_fallbacks = 0;
        // A private singleton attempt repeated with serialized projections
        // before publication; its discarded stage work remains charged.
        std::uint64_t paired_projection_retries = 0;
        std::size_t maximum_batch_rows = 0;
        std::uint64_t reused_phases = 0;
        std::uint64_t initialized_phases = 0;
        std::uint64_t accepted_intervals = 0;
        std::uint64_t rejected_intervals = 0;
        // Completed gathered batches; mixed requests and a private retry count
        // once. Histogram includes rejected interval rows.
        std::uint64_t batches = 0;
        std::uint64_t full_batches = 0;
        std::uint64_t interval_rows = 0;
        std::uint64_t camera_rows = 0;
        std::vector<std::uint64_t> batch_row_counts;
        // Predicate wait wall time includes lock reacquisition. Untimed idle
        // waits for the first request are outside this coalescing observation.
        std::uint64_t coalescing_timeouts = 0;
        std::uint64_t coalescing_underfilled = 0;
        std::uint64_t coalescing_stopped = 0;
        std::uint64_t coalescing_traces_ready = 0;
        double coalescing_wait_ms = 0;
        double maximum_coalescing_wait_ms = 0;
        double execute_ms = 0;  // Serialized, including packing and device calls.
        // Summed worker intervals can overlap each other and the dispatcher.
        // They are not an exclusive component of render wall time.
        std::uint64_t acceleration_calls = 0;
        double acceleration_ms = 0;
    };
    [[nodiscard]] Stats Statistics() const;

  private:
    friend struct RetainedTraceExecutorTestPeer;
    struct Snapshot {
        std::array<double, 40> physical{};
        std::array<double, 4> metric{};
        double chart = 0;
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
        // Known completed attempts remain charged after retry, error or cancellation.
        std::uint32_t completed_stages = 0;
        bool completed = false;
        bool cancelled = false;
        bool registered = false;
    };
    void Run();
    void Execute(std::span<Request*> requests, std::size_t projection_row_budget);
    RetainedCompute& compute_;
    // A lower logical limit reserves both endpoint orders of each queued ray.
    // Zero in the constructor selects the complete storage capacity.
    std::size_t maximum_batch_rows_;
    double maximum_submission_ms_;
    std::function<bool()> should_cancel_;
    mutable std::mutex mutex_;
    std::condition_variable available_, completed_;
    std::deque<Request*> requests_;
    std::map<std::thread::id, Continuation> continuations_;
    // Synchronous Launch/Step permits at most one queued request per thread.
    // Nesting depth preserves distinct membership across nested trace scopes.
    std::map<std::thread::id, std::size_t> active_traces_;
    std::size_t queued_registered_ = 0;
    std::optional<base::Error> error_;
    Stats stats_;
    bool stopping_ = false;
    std::thread dispatcher_;
};

}  // namespace sirius::backend
