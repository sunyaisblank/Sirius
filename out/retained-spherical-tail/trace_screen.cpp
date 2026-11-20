#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera_sampling.h"
#include "sirius/core/spherical_continuation/coupled_admission.h"
#include "sirius/render/trace_domain.h"

#include "exterior_turn.h"

#ifdef SIRIUS_HAS_RETAINED_COMPUTE
#include "sirius/backend/device.h"
#include "sirius/backend/retained_trace_executor.h"
#endif

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <optional>
#include <span>
#include <string_view>

#ifdef SIRIUS_HAS_RETAINED_COMPUTE
namespace sirius::backend {
struct RetainedTraceExecutorTestPeer {
    static std::optional<std::array<RetainedValue, 40>> CommittedState(
        const RetainedTraceExecutor& executor, const core::Lightray& ray,
        const core::GeodesicVariations& variations, const core::KerrSchildParams& params) {
        const auto thread = std::this_thread::get_id();
        RetainedTraceExecutor::Snapshot observed;
        observed.chart = -1;
        observed.metric = {params.M, params.a, params.Q, params.Lambda};
        for (unsigned i = 0; i < 4; ++i) {
            observed.physical[i] = ray.position(i);
            observed.physical[4 + i] = ray.velocity(i);
            for (unsigned column = 0; column < variations.size(); ++column) {
                observed.physical[8 + 8 * column + i] = variations[column].displacement(i);
                observed.physical[12 + 8 * column + i] = variations[column].derivative(i);
            }
        }
        std::lock_guard lock(executor.mutex_);
        const auto state = executor.continuations_.find(thread);
        if (state == executor.continuations_.end() || !executor.active_traces_.contains(thread) ||
            !state->second.finish.valid || !state->second.dopri || state->second.after != observed)
            return {};
        return state->second.finish.physical;
    }
};
}  // namespace sirius::backend
#endif

namespace {
using namespace sirius::core;
using namespace sirius::core::relativity;
using sirius::backend::GeodesicTracer;
using sirius::backend::TracerConfig;
using sirius::backend::TraceResult;
using sirius::backend::TraceStepExecutor;
using sirius::render::BuildTraceDomainParameters;
using sirius::render::kRenderTraceMaximumAttempts;

void Check(bool passed, const char* name) {
    if (passed) return;
    std::fprintf(stderr, "Failed: %s\n", name);
    std::abort();
}

struct Proposal {
    Vec4 position;
    Vec4 tangent;
    SourceSkyAngularMap finite_sky;
    SourceSkyAngularMap infinity_sky;
    double full_affine;
    double normalized_bound;
    int prefix_attempts;
    unsigned panels;
    bool through_turn;
};

std::optional<Proposal> ProjectToPublicChart(const OutgoingKerrSchild& metric,
                                             const CameraRay& camera,
                                             const SphericalCoupledTail& tail, double affine_prefix,
                                             int attempts, double bound) {
    const auto incoming = metric.ToIngoing(tail.position);
    const auto angular = CameraAngularVariations(tail.variations, camera.phase_space);
    if (!incoming || !angular) return {};
    const auto tangent = incoming->Apply(tail.tangent);
    std::array<Vec4, 2> positions, derivatives;
    for (unsigned i = 0; i < positions.size(); ++i) {
        const auto event =
            VarySphericalEndpoint(tail.position, tail.tangent, (*angular)[i].displacement);
        if (!event) return {};
        positions[i] = incoming->Apply(event->displacement);
        derivatives[i] = incoming->Apply((*angular)[i].derivative);
    }
    Metric4d values, inverse;
    Tensor<Dual<double>, 4, 4, 4> gradient;
    metric.Source().Evaluate(incoming->position, values, gradient);
    if (!metric.Source().InverseMetric(incoming->position, inverse)) return {};
    const auto sky =
        MeasureSourceSkyAngularMap(values, inverse, gradient, tangent, positions, derivatives, 1.0);
    if (!sky) return {};
    const auto connection = TensorOps::Christoffel(values, gradient);
    for (unsigned i = 0; i < derivatives.size(); ++i)
        derivatives[i] += NegConnectionContraction(connection, tangent, positions[i]);
    const auto params = metric.Source().GetParams();
    const auto infinity =
        TraceKerrInfinity(params.M, params.a, values, gradient, incoming->position, tangent,
                          positions, derivatives, 1.0);
    if (!infinity) return {};
    return Proposal{incoming->position,
                    tangent,
                    *sky,
                    infinity->map,
                    affine_prefix + tail.affine_length,
                    bound,
                    attempts,
                    tail.panels,
                    tail.periapsis.has_value()};
}

// Run the ordinary coordinator to its usual terminal event. The next Step
// supplies the previous interval's committed endpoint and four columns. A
// rejected interval cannot seed a proposal and no proposal changes the trace.
class TraceScreen final : public TraceStepExecutor {
  public:
    explicit TraceScreen(const TracerConfig& config, TraceStepExecutor* executor = nullptr)
        : config_(config), executor_(executor) {}

    void BeginTrace() override {
        if (executor_) executor_->BeginTrace();
        attempts_ = 0;
        next_radius_ = source_support_;
        affine_prefix_ = 0;
        sampled_first_ = false;
        exterior_turn.reset();
        pending_step_.reset();
        pending_device_curve_ = false;
        proposal.reset();
    }

    void EndTrace() override {
        if (executor_) executor_->EndTrace();
    }

    std::optional<CameraLaunch> Launch(IMetric& metric, double spin,
                                       const CameraRay& camera) override {
        camera_ = camera;
        return executor_ ? executor_->Launch(metric, spin, camera)
                         : LaunchCameraRay(metric, spin, camera);
    }

    bool Step(Lightray& ray, IMetric& metric, const IntegratorConfig& config,
              Rk45CoupledState& coupled, Rk45CoupledComparison& comparison) override {
        if (pending_step_ && cancel_after_checks == 0) {
            affine_prefix_ += *pending_step_;
            if (pending_device_curve_) ++committed_device_intervals;
            Check(!executor_ || pending_device_curve_,
                  "a device prefix uses retained dense curves");
            if (!sampled_first_) {
                sampled_first_ = true;
                const auto* chart = dynamic_cast<const OutgoingKerrSchild*>(&metric);
                if (chart) {
                    const auto p = chart->Source().GetParams();
                    if (p.a == 0 && p.Q == 0 && p.Lambda == 0)
                        exterior_turn = BracketVacuumTurn(ray.position, ray.velocity, p.M,
                                                          SphericalKerrSchildChart::Outgoing,
                                                          source_support_, config_.escape_radius);
                }
                if (exterior_turn)
                    TryProposal(ray, metric, config, coupled, SphericalTailPath::ExteriorTurn);
            }
            if (!proposal) TryProposal(ray, metric, config, coupled);
        }
        ++attempts_;
        const double attempted_step = ray.step_size;
        const bool success =
            executor_ ? executor_->Step(ray, metric, config, coupled, comparison)
                      : Geodesic::IntegrateStepRk45(ray, &metric, config, &coupled, &comparison);
        pending_step_ = success ? std::optional(attempted_step) : std::nullopt;
        pending_device_curve_ =
            executor_ && success &&
            std::all_of(comparison.dopri_positions.begin(), comparison.dopri_positions.end(),
                        [](const auto& curve) { return curve.has_value(); });
        return success;
    }

    std::optional<CoupledSegmentSample> Sample(CoupledTrial trial, double fraction,
                                               const Vec4* normal) override {
        return executor_ ? executor_->Sample(trial, fraction, normal) : std::nullopt;
    }

    void RejectLastInterval() override {
        pending_step_.reset();
        pending_device_curve_ = false;
        if (executor_) executor_->RejectLastInterval();
    }

    bool SupportsSphericalEscape() const override {
        return executor_ && executor_->SupportsSphericalEscape();
    }

    std::optional<sirius::backend::SphericalEscapeTrial> TrySphericalEscape(
        IMetric& metric, const Lightray& ray, const IntegratorConfig& step,
        const Rk45CoupledState& coupled, const sirius::backend::SphericalEscapeRequest& request,
        const std::function<bool()>& should_cancel) override {
        const std::function<bool()> cancel = [&] {
            if (cancel_after_checks > 0 && ++tail_checks >= cancel_after_checks) cancelled = true;
            return cancelled || (should_cancel && should_cancel());
        };
        return executor_->TrySphericalEscape(metric, ray, step, coupled, request, cancel);
    }

    std::optional<Proposal> proposal;
    std::optional<SphericalTailBound> exterior_turn;
    unsigned trials = 0;
    unsigned budget_refusals = 0;
    unsigned accuracy_refusals = 0;
    unsigned other_refusals = 0;
    unsigned committed_device_intervals = 0;
    unsigned expansion_seeded_trials = 0;
    unsigned cancel_after_checks = 0;
    unsigned tail_checks = 0;
    bool cancelled = false;
    double outward_ms = 0;
    double turn_ms = 0;

  private:
    std::expected<SphericalCoupledTail, FiniteTailFailure> Solve(
        const Lightray& ray, OutgoingKerrSchild& metric, const Rk45CoupledState& coupled,
        const SphericalTailControl& tail_control, SphericalTailPath path) {
#ifdef SIRIUS_HAS_RETAINED_COMPUTE
        if (executor_) {
            using namespace finite_tail_detail;
            using namespace sirius::backend;
            const auto* retained = dynamic_cast<const RetainedTraceExecutor*>(executor_);
            Check(retained != nullptr, "retained test executor");
            const auto params = metric.Source().GetParams();
            auto altered = ray;
            altered.position(0) = std::nextafter(ray.position(0), kInfinity);
            Check(!RetainedTraceExecutorTestPeer::CommittedState(*retained, altered,
                                                                 coupled.variations, params),
                  "different physical view cannot claim retained state");
            const auto packet = RetainedTraceExecutorTestPeer::CommittedState(
                *retained, ray, coupled.variations, params);
            if (!packet) return std::unexpected(FiniteTailFailure::Input);
            const auto expansion = [&](unsigned i) {
                const auto& value = (*packet)[i];
                // An accepted expansion defines the initial point of the next
                // local problem. Bound the full three-word sum using outward
                // binary64 interval operations, as opposed to a rounded view.
                return value.IsRepresented()
                           ? Range(value.high) + Range(value.low) + Range(value.tail)
                           : InvalidRange();
            };
            CoupledLaunch launch;
            for (unsigned i = 0; i < 4; ++i) {
                launch.position[i] = expansion(i);
                launch.tangent[i] = expansion(4 + i);
                for (unsigned column = 0; column < coupled.variations.size(); ++column) {
                    launch.displacement[column][i] = expansion(8 + 8 * column + i);
                    launch.derivative[column][i] = expansion(12 + 8 * column + i);
                }
            }
            ++expansion_seeded_trials;
            return ContinueCoupled(metric, launch, source_support_, config_.escape_radius,
                                   tail_control, path);
        }
#endif
        return PrepareSphericalCoupledTail(metric, ray.position, ray.velocity, coupled.variations,
                                           source_support_, config_.escape_radius, tail_control,
                                           path);
    }

    void TryProposal(const Lightray& ray, IMetric& metric, const IntegratorConfig& config,
                     const Rk45CoupledState& coupled,
                     SphericalTailPath path = SphericalTailPath::Outward) {
        auto* outgoing = dynamic_cast<OutgoingKerrSchild*>(&metric);
        if (!outgoing) return;
        const double r = std::hypot(ray.position(1), ray.position(2), ray.position(3));
        const double dot = ray.position(1) * ray.velocity(1) + ray.position(2) * ray.velocity(2) +
                           ray.position(3) * ray.velocity(3);
        const double speed = std::hypot(ray.velocity(1), ray.velocity(2), ray.velocity(3));
        if (path == SphericalTailPath::Outward) {
            if (!(r > next_radius_ && dot > .25 * r * speed)) return;
            next_radius_ = r * 1.25;
        }
        ++trials;
        SphericalTailControl tail_control;
        tail_control.maximum_panels = std::min(
            tail_control.maximum_panels, static_cast<unsigned>(config_.max_steps - attempts_));
        const auto start = std::chrono::steady_clock::now();
        const auto tail = Solve(ray, *outgoing, coupled, tail_control, path);
        const double elapsed =
            std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start)
                .count();
        (path == SphericalTailPath::Outward ? outward_ms : turn_ms) += elapsed;
        if (!tail) {
            if (tail.error() == FiniteTailFailure::EventAccuracy)
                ++accuracy_refusals;
            else
                ++other_refusals;
            return;
        }
        const auto bound = SphericalCoupledStateErrorBound(*tail, config, coupled);
        if (!(bound <= 1.0)) {
            ++budget_refusals;
            return;
        }
        proposal =
            ProjectToPublicChart(*outgoing, camera_, *tail, affine_prefix_, attempts_, bound);
        if (!proposal) ++other_refusals;
    }

    const TracerConfig& config_;
    TraceStepExecutor* executor_;
    const double source_support_ = std::nextafter(static_cast<float>(config_.disk_outer),
                                                  std::numeric_limits<float>::infinity());
    CameraRay camera_;
    int attempts_ = 0;
    bool sampled_first_ = false;
    double next_radius_ = 0;
    double affine_prefix_ = 0;
    std::optional<double> pending_step_;
    bool pending_device_curve_ = false;
};

double ScaledError(double a, double b) {
    return std::abs(a - b) / std::max({1.0, std::abs(a), std::abs(b)});
}

double SkyError(const SourceSkyAngularMap& a, const SourceSkyAngularMap& b) {
    const auto basis_a = MakeCelestialTangentBasis(a.direction);
    const auto basis_b = MakeCelestialTangentBasis(b.direction);
    Check(basis_a.has_value() && basis_b.has_value(), "represented celestial bases");
    double error = 0;
    for (unsigned i = 0; i < a.direction.size(); ++i) {
        error = std::max(error, ScaledError(a.direction[i], b.direction[i]));
        for (unsigned j = 0; j < a.jacobian[0].size(); ++j) {
            // Symmetric sky directions can straddle the least-aligned-axis
            // tie. Compare the physical derivatives in Cartesian components.
            const double da =
                basis_a->first[i] * a.jacobian[0][j] + basis_a->second[i] * a.jacobian[1][j];
            const double db =
                basis_b->first[i] * b.jacobian[0][j] + basis_b->second[i] * b.jacobian[1][j];
            error = std::max(error, ScaledError(da, db));
        }
    }
    return error;
}

void CheckVacuumTurn() {
    Vec4 x, k;
    x(1) = 25;
    k(0) = -std::sqrt(2.0);
    k(1) = -1;
    k(2) = 1;
    constexpr auto chart = SphericalKerrSchildChart::Ingoing;
    const auto turn = BracketVacuumTurn(x, k, 0, chart, 10, 100);
    const long double expected = -25.L / k(0);
    Check(turn && turn->lower <= expected && expected <= turn->upper, "flat closest point");
    Check(!BracketVacuumTurn(x, k, 0, chart, 20, 100), "inward ray enters source support");
    k(1) = 1;
    Check(!BracketVacuumTurn(x, k, 0, chart, 10, 100), "outward ray has no future turn");
}

void CheckBudgetEnclosesOriginalComparison() {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1));
    Vec4 x, k;
    x(1) = 25.;
    // The ingoing-chart outward radial null branch has k^t=-k^r.
    k(0) = -1.;
    k(1) = 1.;
    GeodesicVariations variations;
    variations[0].derivative(3) = .01;
    const auto result = PrepareSphericalCoupledTail(metric, x, k, variations, 20, 200);
    Check(result.has_value(), "central-scale fixture");
    IntegratorConfig step;
    Rk45CoupledState coupled;
    coupled.length_scale = coupled.frequency_scale = 1.0;
    coupled.tolerance = 1.0e-4 / (4.0 * kRenderTraceMaximumAttempts);
    coupled.column_scale = {.008, .008, 1.0, 1.0};
    const auto conservative = SphericalCoupledStateErrorBound(*result, step, coupled);
    Check(conservative > 0.0 && conservative <= 1.0, "certified original scale bound");
    Lightray center;
    center.position = result->position;
    center.velocity = result->tangent;
    for (unsigned corner = 0; corner < 256; ++corner) {
        const auto choose = [&](SphericalTailBound bound, unsigned i) {
            return corner & (1U << (i % 8)) ? bound.upper : bound.lower;
        };
        Lightray alternate = center;
        GeodesicVariations alternate_variations;
        for (unsigned i = 0; i < 4; ++i) {
            alternate.position(i) = choose(result->event_bounds.position[i], i);
            alternate.velocity(i) = choose(result->event_bounds.tangent[i], i + 4);
            for (unsigned j = 0; j < 4; ++j) {
                const auto& variation = result->variation_bounds[j];
                alternate_variations[j].displacement(i) = choose(variation.displacement[i], i + j);
                alternate_variations[j].derivative(i) = choose(variation.derivative[i], i + j + 4);
            }
        }
        Check(Geodesic::CoupledStateError(center, result->variations, alternate,
                                          alternate_variations, step, coupled) <= conservative,
              "live comparison dominated on a represented corner");
    }
    auto malformed = *result;
    malformed.event_bounds.position[0] = {1., -1.};
    Check(std::isinf(SphericalCoupledStateErrorBound(malformed, step, coupled)),
          "invalid central enclosure fails closed");
    malformed = *result;
    malformed.variation_bounds[0].displacement[0] = {1., -1.};
    Check(std::isinf(SphericalCoupledStateErrorBound(malformed, step, coupled)),
          "invalid variation enclosure fails closed");
    coupled.variations[0].displacement(0) = std::numeric_limits<double>::quiet_NaN();
    Check(std::isinf(SphericalCoupledStateErrorBound(*result, step, coupled)),
          "invalid launch variation fails closed");
}

struct FilmPixel {
    int x, y;
};

CameraConfig ScreenCameraConfig() {
    CameraConfig config;
    config.r = 50;
    config.width = config.height = 128;
    config.fov = 60;
    return config;
}

TracerConfig ScreenTracerConfig(const CameraConfig& camera_config) {
    const auto domain = BuildTraceDomainParameters(
        {MetricId::Kerr, 1., 0., camera_config.r, kDefaultMorrisThorneThroatRadius,
         kDefaultAlcubierreBubbleRadius, kDefaultAlcubierreBubbleSigma});
    TracerConfig config;
    config.escape_radius = domain.escape_radius;
    config.max_steps = kRenderTraceMaximumAttempts;
    config.integrator.initial_step = domain.cpu_initial_step;
    config.integrator.min_step = domain.cpu_min_step;
    config.integrator.max_step = domain.max_step;
    config.integrator.abs_tolerance = config.integrator.rel_tolerance = 5.0e-6f;
    return config;
}

void ScreenOriginalGeometry(std::span<const FilmPixel> pixels,
                            TraceStepExecutor* executor = nullptr, bool print_each = false) {
    PinholeCamera camera(ScreenCameraConfig());
    const auto config = ScreenTracerConfig(camera.GetConfig());
    KerrSchildFamily metric(KerrSchildParams::Kerr(1., 0.));
    TraceScreen screen(config, executor);
    GeodesicTracer tracer(&metric, config);
    tracer.SetStepExecutor(&screen);
    CameraSample sample;
    Check(ForEachCameraSample(1, [&](CameraSample next) { sample = next; }) == 1,
          "exact camera sample count");
    unsigned escapes = 0, disk = 0, horizon = 0, proposals = 0, panels = 0;
    unsigned steps = 0, proposed_steps = 0, exterior_turns = 0, turn_proposals = 0, turn_panels = 0;
    double event_error = 0, map_error = 0, infinity_error = 0, worst_bound = 0, widest_turn = 0;
    for (const auto [x, y] : pixels) {
        const auto camera_ray = camera.ProjectFilmForObserver(
            x + sample.image_u, y + sample.image_v, sample.pupil_u, sample.pupil_v);
        Check(camera_ray.has_value(), "camera projection");
        const auto traced = tracer.Trace(camera_ray->ray);
        if (print_each) {
            std::printf(
                "pixel=%d,%d outcome=%d attempts=%d failure=%d termination=%d proposal=%d "
                "turn=%d\n",
                x, y, static_cast<int>(traced.outcome), traced.steps_taken,
                static_cast<int>(traced.coupled_failure), traced.integrator_termination,
                static_cast<bool>(screen.proposal),
                screen.proposal && screen.proposal->through_turn);
            std::fflush(stdout);
        }
        Check(!traced.cancelled && !traced.numerical_failure, "ordinary trace admission");
        steps += traced.steps_taken;
        if (traced.outcome != TraceResult::Outcome::Escaped) {
            Check(!screen.proposal && !screen.exterior_turn,
                  "no continuation bypasses an opaque event");
            Check(traced.outcome == TraceResult::Outcome::DiskHit ||
                      traced.outcome == TraceResult::Outcome::Horizon,
                  "opaque terminal event");
            if (traced.outcome == TraceResult::Outcome::DiskHit)
                ++disk;
            else
                ++horizon;
            proposed_steps += traced.steps_taken;
            continue;
        }
        ++escapes;
        if (screen.exterior_turn) {
            const auto periapsis = *screen.exterior_turn;
            Check(
                periapsis.lower > config.disk_outer && traced.min_radius + 1e-5 >= periapsis.lower,
                "bounded exterior periapsis agrees with full trace clearance");
            ++exterior_turns;
            widest_turn = std::max(widest_turn, periapsis.upper - periapsis.lower);
        }
        if (!screen.proposal) {
            proposed_steps += traced.steps_taken;
            continue;
        }
        const auto& tail = *screen.proposal;
        ++proposals;
        proposed_steps += tail.prefix_attempts;
        panels += tail.panels;
        if (tail.through_turn) {
            ++turn_proposals;
            turn_panels += tail.panels;
        }
        worst_bound = std::max(worst_bound, tail.normalized_bound);
        Check(traced.final_tangent && traced.beam.finite_source_map &&
                  traced.beam.infinity_source_map,
              "ordinary finite-boundary publication");
        event_error = std::max(event_error, ScaledError(tail.full_affine, traced.affine_length));
        for (unsigned i = 0; i < 4; ++i) {
            event_error =
                std::max({event_error, ScaledError(tail.position(i), traced.final_position(i)),
                          ScaledError(tail.tangent(i), (*traced.final_tangent)(i))});
        }
        map_error = std::max(map_error, SkyError(tail.finite_sky, *traced.beam.finite_source_map));
        infinity_error = std::max(
            infinity_error, SkyError(tail.infinity_sky, traced.beam.infinity_source_map->map));
    }
    std::printf(
        "rays=%zu escape=%u disk=%u capture=%u proposals=%u trials=%u budget_refusals=%u "
        "accuracy_refusals=%u other_refusals=%u\n",
        pixels.size(), escapes, disk, horizon, proposals, screen.trials, screen.budget_refusals,
        screen.accuracy_refusals, screen.other_refusals);
    std::printf(
        "baseline_steps=%u projected_retained_attempts=%u quadrature_panels=%u "
        "max_coupled_bound=%.9g max_scaled_event_error=%.9g max_scaled_finite_sky_error=%.9g "
        "max_scaled_infinity_sky_error=%.9g committed_device_intervals=%u "
        "expansion_seeded_trials=%u\n",
        steps, proposed_steps, panels, worst_bound, event_error, map_error, infinity_error,
        screen.committed_device_intervals, screen.expansion_seeded_trials);
    std::printf(
        "certified_exterior_turns=%u turn_proposals=%u widest_turn_bracket=%.9g turn_panels=%u "
        "turn_ms=%.6f outward_ms=%.6f\n",
        exterior_turns, turn_proposals, widest_turn, turn_panels, screen.turn_ms,
        screen.outward_ms);
    std::fflush(stdout);
    Check(proposals > 0 && disk > 0 && horizon > 0, "joint escape, disk, and capture coverage");
    Check(exterior_turns > 0 && exterior_turns < escapes, "mixed inward vacuum turn eligibility");
    Check(turn_proposals > 0, "at least one full turn meets the live budget");
    Check(!executor || (screen.committed_device_intervals > 0 &&
                        screen.expansion_seeded_trials == screen.trials),
          "device-owned inputs seed every Vulkan proposal");
    Check(event_error <= 1e-4 && map_error <= 1e-4 && infinity_error <= 1e-4,
          "diagnostic agreement with the complete trace");
}

void ScreenCpu(bool print_each) {
    constexpr std::array coordinates{4, 25, 49, 64, 79, 103, 123};
    std::array<FilmPixel, coordinates.size() * coordinates.size()> pixels;
    auto next = pixels.begin();
    for (const int x : coordinates)
        for (const int y : coordinates) *next++ = {x, y};
    ScreenOriginalGeometry(pixels, nullptr, print_each);
}

#ifdef SIRIUS_HAS_RETAINED_COMPUTE
template <typename T>
T Require(sirius::base::Expected<T> result) {
    if (result) return std::move(*result);
    std::fprintf(stderr, "%s\n", result.error().Description().c_str());
    std::exit(EXIT_FAILURE);
}

constexpr std::array kVulkanPixels{FilmPixel{4, 4}, FilmPixel{25, 49}, FilmPixel{64, 64},
                                   FilmPixel{49, 64}};

void CheckPublishedVulkan(sirius::backend::RetainedCompute& compute) {
    using namespace sirius::backend;
    RetainedTraceExecutor ordinary(compute, {}, 1000.0, 4);
    RetainedTraceExecutor candidate(compute, {}, 1000.0, 4,
                                    RetainedTraceExecutor::SphericalEscape::Certified);
    const PinholeCamera camera(ScreenCameraConfig());
    auto config = ScreenTracerConfig(camera.GetConfig());
    config.enable_ray_bundles = true;
    config.bundle_point_source = true;
    KerrSchildFamily metric(KerrSchildParams::Kerr(1., 0.));
    GeodesicTracer tracer(&metric, config);
    CameraSample sample;
    Check(ForEachCameraSample(1, [&](CameraSample next) { sample = next; }) == 1,
          "exact camera sample count");
    const auto project = [&](FilmPixel pixel) {
        const auto camera_ray = camera.ProjectFilmForObserver(
            pixel.x + sample.image_u, pixel.y + sample.image_v, sample.pupil_u, sample.pupil_v);
        Check(camera_ray.has_value(), "publication camera projection");
        return camera_ray->ray;
    };
    double event_error = 0, sky_error = 0, beam_error = 0;
    unsigned ordinary_attempts = 0, attempts = 0, panels = 0;
    for (const auto pixel : kVulkanPixels) {
        const auto input = project(pixel);
        tracer.SetStepExecutor(&ordinary);
        const auto baseline = tracer.Trace(input);
        tracer.SetStepExecutor(&candidate);
        const auto trial = tracer.Trace(input);
        Check(!baseline.numerical_failure && !trial.numerical_failure && !trial.cancelled &&
                  baseline.final_tangent && trial.final_tangent &&
                  baseline.outcome == trial.outcome &&
                  baseline.terminal_chart == trial.terminal_chart &&
                  baseline.beam.valid == trial.beam.valid,
              "same published terminal contracts");
        const bool escaped = trial.outcome == TraceResult::Outcome::Escaped;
        Check((trial.spherical_terminal_panels > 0) == escaped &&
                  trial.spherical_terminal_turn == (pixel.x == 4),
              "expected live terminal route");
        Check(trial.num_disk_crossings == baseline.num_disk_crossings &&
                  trial.disk_radius == baseline.disk_radius &&
                  trial.redshift == baseline.redshift &&
                  trial.optical_depth == baseline.optical_depth &&
                  trial.volumetric_emission == baseline.volumetric_emission,
              "same accumulated transfer");
        Check(escaped ? trial.steps_taken < baseline.steps_taken
                      : trial.steps_taken == baseline.steps_taken,
              "actual omitted retained attempts");
        ordinary_attempts += baseline.steps_taken;
        attempts += trial.steps_taken;
        panels += trial.spherical_terminal_panels;
        event_error =
            std::max(event_error, ScaledError(trial.affine_length, baseline.affine_length));
        for (int i = 0; i < 4; ++i) {
            event_error = std::max(
                {event_error, ScaledError(trial.final_position(i), baseline.final_position(i)),
                 ScaledError((*trial.final_tangent)(i), (*baseline.final_tangent)(i))});
        }
        if (escaped) {
            Check(trial.beam.finite_source_map && baseline.beam.finite_source_map &&
                      trial.beam.infinity_source_map && baseline.beam.infinity_source_map &&
                      !trial.beam.infinity_source_failure,
                  "both published source maps");
            sky_error =
                std::max({sky_error,
                          SkyError(*trial.beam.finite_source_map, *baseline.beam.finite_source_map),
                          SkyError(trial.beam.infinity_source_map->map,
                                   baseline.beam.infinity_source_map->map)});
            beam_error =
                std::max({beam_error, ScaledError(trial.beam.semi_major, baseline.beam.semi_major),
                          ScaledError(trial.beam.semi_minor, baseline.beam.semi_minor),
                          ScaledError(trial.beam.area_ratio, baseline.beam.area_ratio)});
        }
        std::printf(
            "pixel=%d,%d outcome=%d ordinary_attempts=%d certified_attempts=%d panels=%u turn=%d\n",
            pixel.x, pixel.y, static_cast<int>(trial.outcome), baseline.steps_taken,
            trial.steps_taken, trial.spherical_terminal_panels, trial.spherical_terminal_turn);
        std::fflush(stdout);
    }
    Check(event_error <= 1e-7 && sky_error <= 1e-7 && beam_error <= 1e-6,
          "published terminal numerical comparison");

    TraceScreen interrupted(config, &candidate);
    interrupted.cancel_after_checks = 8;
    tracer.SetStepExecutor(&interrupted);
    tracer.SetCancellationCallback([&] { return interrupted.cancelled; });
    const auto input = project(kVulkanPixels.front());
    const auto cancelled = tracer.Trace(input);
    Check(cancelled.cancelled && interrupted.tail_checks == 8 && !cancelled.numerical_failure &&
              !cancelled.final_tangent && !cancelled.beam.valid &&
              !cancelled.beam.finite_source_map && !cancelled.beam.infinity_source_map &&
              cancelled.spherical_terminal_panels == 0 && cancelled.affine_length == 0 &&
              cancelled.num_disk_crossings == 0,
          "cancelled continuation publishes no physical state");
    tracer.SetStepExecutor(&candidate);
    tracer.SetCancellationCallback({});
    const auto retry = tracer.Trace(input);
    Check(retry.outcome == TraceResult::Outcome::Escaped && !retry.cancelled &&
              !retry.numerical_failure && retry.spherical_terminal_turn &&
              retry.beam.finite_source_map,
          "retained executor resumes after cancellation");
    Check(!ordinary.Error() && !candidate.Error(), "no sticky device errors");
    std::printf(
        "published_ordinary_attempts=%u published_certified_attempts=%u panels=%u "
        "event_error=%.9g sky_error=%.9g beam_error=%.9g cancelled_panels_check=%u recovered=%d\n",
        ordinary_attempts, attempts, panels, event_error, sky_error, beam_error,
        interrupted.tail_checks, retry.spherical_terminal_turn);
}

void ScreenVulkan(bool publish) {
    using namespace sirius::backend;
    const auto inventory = Require(EnumerateVulkanDevices());
    const auto index = Require(ResolveVulkanDeviceIndex(inventory));
    const auto device = Require(CreateVulkanDevice(index));
    const auto& info = device->Info();
    std::printf("device=%s driver=%s fp32_portable=%d\n", info.name.c_str(),
                info.driver_info.c_str(), RetainedUsesPortableArithmetic(info));
    std::fflush(stdout);
    auto compute = Require(RetainedCompute::Create(*device, 8));
    RetainedCompute::PreparationStats preparation;
    const auto ready = compute->PrepareSoftwareRendererStages(preparation);
    if (!ready) {
        std::fprintf(stderr, "Vulkan preparation failed after %.3f ms: %s\n", preparation.wall_ms,
                     ready.error().Description().c_str());
        std::exit(EXIT_FAILURE);
    }
    std::printf("prepared_ms=%.3f allocated_bytes=%llu\n", preparation.wall_ms,
                static_cast<unsigned long long>(device->BufferAllocationBytes()));
    std::fflush(stdout);
    if (publish) {
        CheckPublishedVulkan(*compute);
        return;
    }
    RetainedTraceExecutor executor(*compute, {}, 1000.0, 4);
    ScreenOriginalGeometry(kVulkanPixels, &executor, true);
    const auto stats = executor.Statistics();
    std::printf("live_initialized=%llu live_reused=%llu accepted=%llu tracer_rollbacks=%llu\n",
                static_cast<unsigned long long>(stats.initialized_phases),
                static_cast<unsigned long long>(stats.reused_phases),
                static_cast<unsigned long long>(stats.accepted_intervals),
                static_cast<unsigned long long>(stats.tracer_rollbacks));
    Check(!executor.Error() && stats.accepted_intervals > 0 && stats.reused_phases > 0,
          "actual retained execution and reuse");
}
#endif
}  // namespace

int main(int argc, char** argv) {
    CheckBudgetEnclosesOriginalComparison();
    CheckVacuumTurn();
    const std::string_view mode = argc > 1 ? argv[1] : "--cpu";
#ifdef SIRIUS_HAS_RETAINED_COMPUTE
    if (argc == 2 && (mode == "--vulkan" || mode == "--vulkan-publish")) {
        ScreenVulkan(mode == "--vulkan-publish");
        return 0;
    }
#endif
    Check(argc <= 2 && (mode == "--cpu" || mode == "--cpu-details"), "trace screen mode");
    ScreenCpu(mode == "--cpu-details");
}
