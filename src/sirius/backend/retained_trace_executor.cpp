#include "sirius/backend/retained_trace_executor.h"

#include "sirius/core/metrics/outgoing_kerr_schild.h"
#include "sirius/core/twofold.h"

#include <algorithm>
#include <chrono>
#include <cmath>

namespace sirius::backend {
namespace {
double Rounded(const RetainedValue& value) {
    return (core::Twofold(value.high) + core::Twofold(value.low) + core::Twofold(value.tail))
        .Rounded();
}

void Physical(const RetainedEndpointOutput& input, core::Lightray& ray,
              core::GeodesicVariations& variations) {
    for (std::size_t i = 0; i < 4; ++i) {
        ray.position(static_cast<int>(i)) = Rounded(input.physical[i]);
        ray.velocity(static_cast<int>(i)) = Rounded(input.physical[4 + i]);
        for (std::size_t column = 0; column < 4; ++column) {
            variations[column].displacement(static_cast<int>(i)) =
                Rounded(input.physical[8 + 8 * column + i]);
            variations[column].derivative(static_cast<int>(i)) =
                Rounded(input.physical[12 + 8 * column + i]);
        }
    }
}

core::CoupledSegmentIncrement Increment(const std::array<RetainedValue, 20>& input) {
    core::CoupledSegmentIncrement result;
    for (std::size_t i = 0; i < 4; ++i) {
        result.position(static_cast<int>(i)) = Rounded(input[i]);
        for (std::size_t c = 0; c < 4; ++c)
            result.displacement[c](static_cast<int>(i)) = Rounded(input[4 + 4 * c + i]);
    }
    return result;
}
}  // namespace

RetainedTraceExecutor::RetainedTraceExecutor(RetainedCompute& compute,
                                             std::function<bool()> should_cancel,
                                             double maximum_submission_ms)
    : compute_(compute),
      maximum_submission_ms_(maximum_submission_ms),
      should_cancel_(std::move(should_cancel)),
      dispatcher_(&RetainedTraceExecutor::Run, this) {}

RetainedTraceExecutor::~RetainedTraceExecutor() {
    {
        std::lock_guard lock(mutex_);
        stopping_ = true;
    }
    available_.notify_one();
    dispatcher_.join();
}

std::optional<base::Error> RetainedTraceExecutor::Error() const {
    std::lock_guard lock(mutex_);
    return error_;
}

void RetainedTraceExecutor::BeginTrace() {
    std::lock_guard lock(mutex_);
    continuations_.erase(std::this_thread::get_id());
}

void RetainedTraceExecutor::EndTrace() { BeginTrace(); }

void RetainedTraceExecutor::RejectLastInterval() {
    std::lock_guard lock(mutex_);
    const auto found = continuations_.find(std::this_thread::get_id());
    if (found != continuations_.end()) found->second.finish = {};
}

RetainedTraceExecutor::Stats RetainedTraceExecutor::Statistics() const {
    std::lock_guard lock(mutex_);
    return stats_;
}

void RetainedTraceExecutor::Run() {
    std::size_t batch_limit = compute_.Capacity();
    std::size_t safety_cap = batch_limit;
    for (;;) {
        std::vector<Request*> batch;
        {
            std::unique_lock lock(mutex_);
            available_.wait(lock, [&] { return stopping_ || !requests_.empty(); });
            if (stopping_ && requests_.empty()) return;
            // One bounded coalescing window; no worker count or full batch is
            // required for progress when rays finish at different times.
            available_.wait_for(lock, std::chrono::milliseconds(1),
                                [&] { return stopping_ || requests_.size() >= batch_limit; });
            while (!requests_.empty() && batch.size() < batch_limit) {
                batch.push_back(requests_.front());
                requests_.pop_front();
            }
        }
        Execute(batch);
        const double peak = compute_.TakeSubmissionPeakMs();
        {
            std::lock_guard lock(mutex_);
            if (std::any_of(batch.begin(), batch.end(),
                            [](const Request* r) { return !r->camera; }))
                ++stats_.interval_batches;
            if (std::any_of(batch.begin(), batch.end(), [](const Request* r) { return r->camera; }))
                ++stats_.camera_batches;
            const bool safety = maximum_submission_ms_ > 0 && peak > maximum_submission_ms_;
            if (safety && batch.size() == 1 && !error_)
                error_.emplace(base::ErrorDomain::kDevice, "dispatch retained renderer",
                               "single-ray submission exceeded the safety duration");
            if ((safety || peak == 0) && batch_limit > 1) {
                safety_cap = std::min(safety_cap,
                                      peak == 0 ? 1 : std::max<std::size_t>(1, batch.size() / 2));
                batch_limit = std::min(batch_limit, safety_cap);
                ++stats_.safety_fallbacks;
            }
            if (compute_.DispatchTargetMs() > 0 && peak > compute_.DispatchTargetMs() &&
                batch_limit > 1) {
                batch_limit = std::max<std::size_t>(
                    1, std::min(batch_limit - 1,
                                static_cast<std::size_t>(batch.size() * .8 *
                                                         compute_.DispatchTargetMs() / peak)));
                ++stats_.batch_subdivisions;
            } else if (compute_.DispatchTargetMs() > 0 && peak > 0 &&
                       peak < .5 * compute_.DispatchTargetMs() && batch.size() == batch_limit) {
                // Recover throughput after a temporary soft-target overshoot,
                // while preserving any irreversible safety ceiling.
                batch_limit = std::min(safety_cap, batch_limit * 2);
            }
            stats_.maximum_batch_rows = std::max(stats_.maximum_batch_rows, batch.size());
            for (auto* request : batch) {
                if (!request->camera && request->result) {
                    if (request->result->admissible)
                        ++stats_.accepted_intervals;
                    else
                        ++stats_.rejected_intervals;
                }
                if (request->camera && !request->camera_result && !error_)
                    error_ = request->camera_result.error();
                if (!request->result && !error_) error_ = request->result.error();
                if (error_) {
                    request->camera_result = std::unexpected(*error_);
                    request->result = std::unexpected(*error_);
                }
                request->completed = true;
            }
        }
        completed_.notify_all();
    }
}

void RetainedTraceExecutor::Execute(std::span<Request*> requests) {
    if (std::any_of(requests.begin(), requests.end(),
                    [](const Request* request) { return request->camera; })) {
        std::vector<RetainedRayCameraInput> packets;
        std::vector<Request*> cameras, steps;
        for (auto* request : requests) {
            if (request->camera) {
                cameras.push_back(request);
                packets.push_back(request->camera_input);
            } else
                steps.push_back(request);
        }
        const auto outputs = compute_.RayCamera(packets);
        for (std::size_t row = 0; row < cameras.size(); ++row)
            cameras[row]->camera_result =
                outputs ? base::Expected<RetainedCameraOutput>((*outputs)[row])
                        : std::unexpected(outputs.error());
        if (!steps.empty()) Execute(steps);
        return;
    }
    const auto fail = [&](const base::Error& error) {
        for (auto* request : requests) request->result = std::unexpected(error);
    };
    bool initialize = false;
    std::vector<RetainedInitializeInput> packets(requests.size());
    for (std::size_t row = 0; row < requests.size(); ++row)
        if (!requests[row]->interval.start.valid) {
            initialize = true;
            packets[row] = requests[row]->initialization;
        }
    if (initialize) {
        const auto phase = compute_.Initialize(packets);
        if (!phase) {
            fail(phase.error());
            return;
        }
        std::vector<RetainedEndpointInput> endpoints(requests.size());
        for (std::size_t row = 0; row < requests.size(); ++row)
            if ((*phase)[row].valid) {
                std::copy(requests[row]->interval.metric.begin(),
                          requests[row]->interval.metric.end(), endpoints[row].values.begin());
                std::copy((*phase)[row].phase.begin(), (*phase)[row].phase.end(),
                          endpoints[row].values.begin() + 4);
                endpoints[row].values[44] =
                    RetainedValue::FromDouble(requests[row]->interval.chart);
            }
        const auto projected = compute_.Endpoint(endpoints);
        if (!projected) {
            fail(projected.error());
            return;
        }
        for (std::size_t row = 0; row < requests.size(); ++row)
            if (!requests[row]->interval.start.valid)
                requests[row]->interval.start = (*projected)[row];
    }
    std::vector<RetainedIntervalInput> intervals;
    for (const auto* request : requests) intervals.push_back(request->interval);
    const auto results = AttemptRetainedIntervals(compute_, intervals);
    if (!results) {
        fail(results.error());
        return;
    }
    for (std::size_t row = 0; row < requests.size(); ++row) requests[row]->result = (*results)[row];
}

std::optional<core::CameraLaunch> RetainedTraceExecutor::Launch(core::IMetric& metric,
                                                                double absolute_spin,
                                                                const core::CameraRay& camera) {
    const auto* family = dynamic_cast<core::KerrSchildFamily*>(&metric);
    if (!family || !core::IsRepresentedCameraRay(camera) || !camera.active ||
        !std::isfinite(absolute_spin))
        return std::nullopt;
    if (should_cancel_ && should_cancel_()) return std::nullopt;
    const auto p = family->GetParams();
    const core::coordinates::Vec4Bl origin(camera.origin(0), camera.origin(1), camera.origin(2),
                                           camera.origin(3));
    const auto cart = core::coordinates::BlToKerrSchildCart(origin, absolute_spin);
    const double theta = camera.origin(2), phi = std::atan2(cart.y, cart.x);
    const double st = std::sin(theta), ct = std::cos(theta), sp = std::sin(phi), cp = std::cos(phi);
    std::array<double, 45> values{};
    const std::array<double, 25> base{p.M,
                                      p.a,
                                      p.Q,
                                      p.Lambda,
                                      cart.t,
                                      cart.x,
                                      cart.y,
                                      cart.z,
                                      st * cp,
                                      st * sp,
                                      ct,
                                      -ct * cp,
                                      -ct * sp,
                                      st,
                                      -sp,
                                      cp,
                                      0,
                                      -camera.beta_forward,
                                      -camera.beta_up,
                                      camera.beta_right,
                                      camera.aperture_right,
                                      -camera.aperture_up,
                                      camera.direction(1),
                                      camera.direction(2),
                                      camera.direction(3)};
    std::copy(base.begin(), base.end(), values.begin());
    if (camera.phase_space)
        for (std::size_t column = 0; column < 4; ++column) {
            for (std::size_t axis = 0; axis < 3; ++axis)
                values[25 + 4 * axis + column] = camera.phase_space->direction[axis][column];
            values[37 + column] = camera.phase_space->pupil_right[column];
            values[41 + column] = -camera.phase_space->pupil_up[column];
        }
    Request request;
    request.camera = true;
    for (std::size_t i = 0; i < values.size(); ++i)
        request.camera_input.values[i] = RetainedValue::FromDouble(values[i]);
    {
        std::unique_lock lock(mutex_);
        if (error_ || stopping_) return std::nullopt;
        requests_.push_back(&request);
        available_.notify_one();
        completed_.wait(lock, [&] { return request.completed; });
    }
    if (!request.camera_result || !request.camera_result->valid) return std::nullopt;
    const auto& output = request.camera_result->values;
    core::CameraLaunch launch;
    for (int mu = 0; mu < 4; ++mu) {
        launch.position(mu) = Rounded(output[mu]);
        launch.tangent(mu) = Rounded(output[4 + mu]);
        launch.observer.time(mu) = Rounded(output[8 + mu]);
        for (int axis = 0; axis < 3; ++axis)
            launch.observer.spatial[axis](mu) = Rounded(output[12 + 4 * axis + mu]);
        for (int column = 0; column < 4; ++column) {
            launch.variations[column].displacement(mu) = Rounded(output[24 + 4 * column + mu]);
            launch.variations[column].derivative(mu) = Rounded(output[56 + 4 * column + mu]);
        }
    }
    return launch;
}

bool RetainedTraceExecutor::Step(core::Lightray& ray, core::IMetric& metric,
                                 const core::IntegratorConfig& config,
                                 core::Rk45CoupledState& coupled,
                                 core::Rk45CoupledComparison& comparison) {
    if (should_cancel_ && should_cancel_()) {
        ray.terminated = 3;
        coupled.failure = core::CoupledStepFailure::InvalidState;
        return false;
    }
    auto* outgoing = dynamic_cast<core::OutgoingKerrSchild*>(&metric);
    auto* family = outgoing ? &outgoing->Source() : dynamic_cast<core::KerrSchildFamily*>(&metric);
    if (!family) {
        ray.terminated = 3;
        coupled.failure = core::CoupledStepFailure::InvalidState;
        return false;
    }
    const auto parameters = family->GetParams();
    Snapshot snapshot;
    snapshot.metric = {parameters.M, parameters.a, parameters.Q, parameters.Lambda};
    snapshot.chart = outgoing ? -1 : 1;
    snapshot.affine = ray.proper_time;
    const auto snapshot_physical = [&](Snapshot& target) {
        for (std::size_t i = 0; i < 4; ++i) {
            target.physical[i] = ray.position(static_cast<int>(i));
            target.physical[4 + i] = ray.velocity(static_cast<int>(i));
            for (std::size_t c = 0; c < 4; ++c) {
                target.physical[8 + 8 * c + i] =
                    coupled.variations[c].displacement(static_cast<int>(i));
                target.physical[12 + 8 * c + i] =
                    coupled.variations[c].derivative(static_cast<int>(i));
            }
        }
    };
    snapshot_physical(snapshot);
    Request request;
    for (std::size_t i = 0; i < 4; ++i)
        request.interval.metric[i] = request.initialization.values[i] =
            RetainedValue::FromDouble(snapshot.metric[i]);
    for (std::size_t i = 0; i < 40; ++i)
        request.initialization.values[4 + i] = RetainedValue::FromDouble(snapshot.physical[i]);
    request.initialization.values[44] = RetainedValue::FromDouble(snapshot.chart);
    request.interval.chart = snapshot.chart;
    request.interval.interval = ray.step_size;
    request.interval.control = {config, coupled.length_scale, coupled.frequency_scale,
                                coupled.tolerance, coupled.column_scale};
    const auto thread = std::this_thread::get_id();
    {
        std::unique_lock lock(mutex_);
        if (error_ || stopping_) {
            ray.terminated = 3;
            coupled.failure = core::CoupledStepFailure::InvalidState;
            return false;
        }
        const auto cached = continuations_.find(thread);
        if (cached != continuations_.end()) {
            if (snapshot == cached->second.after && cached->second.finish.valid)
                request.interval.start = cached->second.finish;
            else if (snapshot == cached->second.before)
                request.interval.start = cached->second.start;
            if (request.interval.start.valid) ++stats_.reused_phases;
        }
        requests_.push_back(&request);
        available_.notify_one();
        completed_.wait(lock, [&] { return request.completed; });
    }
    if (!request.result) {
        ray.terminated = 3;
        coupled.failure = core::CoupledStepFailure::InvalidState;
        return false;
    }
    const auto& result = *request.result;
    coupled.central_stages += result.attempted_stages;
    coupled.variation_stages += result.attempted_stages;
    coupled.failure = result.failure;
    Continuation continuation;
    continuation.before = snapshot;
    continuation.start = request.interval.start;
    if (!result.admissible) {
        if (ray.step_size <= config.min_step)
            ray.terminated = 5;
        else
            ray.step_size = std::max(config.min_step, ray.step_size * .5f);
        std::lock_guard lock(mutex_);
        continuations_[thread] = continuation;
        return false;
    }
    const float interval = ray.step_size;
    comparison = {};
    comparison.lower_order = comparison.midpoint = comparison.refined_endpoint = ray;
    Physical(result.lower, comparison.lower_order, comparison.lower_variations);
    Physical(result.midpoint, comparison.midpoint, comparison.midpoint_variations);
    Physical(result.refined, comparison.refined_endpoint, comparison.refined_variations);
    comparison.lower_order.proper_time += interval;
    comparison.midpoint.proper_time += interval * .5f;
    comparison.refined_endpoint.proper_time += interval;
    comparison.full_increment = Increment(result.full_increment);
    comparison.lower_increment = Increment(result.lower_increment);
    comparison.midpoint_increment = Increment(result.midpoint_increment);
    comparison.refined_increment = Increment(result.refined_increment);
    comparison.error_ratio = result.error_ratio;
    Physical(result.full, ray, coupled.variations);
    ray.acceleration = core::Geodesic::CalculateAcceleration(ray.velocity, ray.position, &metric);
    ray.proper_time += interval;
    ray.coordinate_time += static_cast<float>(interval * std::abs(ray.velocity(0)));
    ray.step_size = core::Geodesic::ComputeOptimalStep(
        interval, static_cast<float>(result.error_ratio), 1, config);
    continuation.after = snapshot;
    continuation.after.affine = ray.proper_time;
    snapshot_physical(continuation.after);
    continuation.finish = result.full;
    {
        std::lock_guard lock(mutex_);
        continuations_[thread] = continuation;
    }
    return true;
}

}  // namespace sirius::backend
