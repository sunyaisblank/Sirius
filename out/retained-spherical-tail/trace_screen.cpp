#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera_sampling.h"
#include "sirius/render/trace_domain.h"

#include "coupled_admission.h"
#include "exterior_turn.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <optional>

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
class CpuTraceScreen final : public TraceStepExecutor {
  public:
    explicit CpuTraceScreen(const TracerConfig& config) : config_(config) {}

    void BeginTrace() override {
        attempts_ = 0;
        next_radius_ = source_support_;
        affine_prefix_ = 0;
        sampled_first_ = false;
        exterior_turn.reset();
        pending_step_.reset();
        proposal.reset();
    }

    std::optional<CameraLaunch> Launch(IMetric& metric, double spin,
                                       const CameraRay& camera) override {
        camera_ = camera;
        return LaunchCameraRay(metric, spin, camera);
    }

    bool Step(Lightray& ray, IMetric& metric, const IntegratorConfig& config,
              Rk45CoupledState& coupled, Rk45CoupledComparison& comparison) override {
        if (pending_step_) {
            affine_prefix_ += *pending_step_;
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
            Geodesic::IntegrateStepRk45(ray, &metric, config, &coupled, &comparison);
        pending_step_ = success ? std::optional(attempted_step) : std::nullopt;
        return success;
    }

    void RejectLastInterval() override { pending_step_.reset(); }

    std::optional<Proposal> proposal;
    std::optional<SphericalTailBound> exterior_turn;
    unsigned trials = 0;
    unsigned budget_refusals = 0;
    unsigned accuracy_refusals = 0;
    unsigned other_refusals = 0;
    double outward_ms = 0;
    double turn_ms = 0;

  private:
    void TryProposal(const Lightray& ray, IMetric& metric, const IntegratorConfig& config,
                     const Rk45CoupledState& coupled,
                     SphericalTailPath path = SphericalTailPath::Outward) {
        const auto* outgoing = dynamic_cast<const OutgoingKerrSchild*>(&metric);
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
        const auto tail =
            PrepareSphericalCoupledTail(metric, ray.position, ray.velocity, coupled.variations,
                                        source_support_, config_.escape_radius, tail_control, path);
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
    const double source_support_ = std::nextafter(static_cast<float>(config_.disk_outer),
                                                  std::numeric_limits<float>::infinity());
    CameraRay camera_;
    int attempts_ = 0;
    bool sampled_first_ = false;
    double next_radius_ = 0;
    double affine_prefix_ = 0;
    std::optional<double> pending_step_;
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

void ScreenOriginalGeometry() {
    CameraConfig camera_config;
    camera_config.r = 50;
    camera_config.width = camera_config.height = 128;
    camera_config.fov = 60;
    PinholeCamera camera(camera_config);
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
    KerrSchildFamily metric(KerrSchildParams::Kerr(1., 0.));
    CpuTraceScreen screen(config);
    GeodesicTracer tracer(&metric, config);
    tracer.SetStepExecutor(&screen);
    CameraSample sample;
    Check(ForEachCameraSample(1, [&](CameraSample next) { sample = next; }) == 1,
          "exact camera sample count");
    constexpr std::array coordinates{4, 25, 49, 64, 79, 103, 123};
    unsigned escapes = 0, disk = 0, horizon = 0, proposals = 0, panels = 0;
    unsigned steps = 0, proposed_steps = 0, exterior_turns = 0, turn_proposals = 0, turn_panels = 0;
    double event_error = 0, map_error = 0, infinity_error = 0, worst_bound = 0, widest_turn = 0;
    for (const int x : coordinates)
        for (const int y : coordinates) {
            const auto camera_ray = camera.ProjectFilmForObserver(
                x + sample.image_u, y + sample.image_v, sample.pupil_u, sample.pupil_v);
            Check(camera_ray.has_value(), "camera projection");
            const auto traced = tracer.Trace(camera_ray->ray);
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
                Check(periapsis.lower > config.disk_outer &&
                          traced.min_radius + 1e-5 >= periapsis.lower,
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
            event_error =
                std::max(event_error, ScaledError(tail.full_affine, traced.affine_length));
            for (unsigned i = 0; i < 4; ++i) {
                event_error =
                    std::max({event_error, ScaledError(tail.position(i), traced.final_position(i)),
                              ScaledError(tail.tangent(i), (*traced.final_tangent)(i))});
            }
            map_error =
                std::max(map_error, SkyError(tail.finite_sky, *traced.beam.finite_source_map));
            infinity_error = std::max(
                infinity_error, SkyError(tail.infinity_sky, traced.beam.infinity_source_map->map));
        }
    std::printf(
        "rays=%zu escape=%u disk=%u capture=%u proposals=%u trials=%u budget_refusals=%u "
        "accuracy_refusals=%u other_refusals=%u\n",
        coordinates.size() * coordinates.size(), escapes, disk, horizon, proposals, screen.trials,
        screen.budget_refusals, screen.accuracy_refusals, screen.other_refusals);
    std::printf(
        "baseline_steps=%u projected_retained_attempts=%u quadrature_panels=%u "
        "max_coupled_bound=%.9g max_scaled_event_error=%.9g max_scaled_finite_sky_error=%.9g "
        "max_scaled_infinity_sky_error=%.9g\n",
        steps, proposed_steps, panels, worst_bound, event_error, map_error, infinity_error);
    std::printf(
        "certified_exterior_turns=%u turn_proposals=%u widest_turn_bracket=%.9g turn_panels=%u "
        "turn_ms=%.6f outward_ms=%.6f\n",
        exterior_turns, turn_proposals, widest_turn, turn_panels, screen.turn_ms,
        screen.outward_ms);
    std::fflush(stdout);
    Check(proposals > 0 && disk > 0 && horizon > 0, "joint escape, disk, and capture coverage");
    Check(exterior_turns > 0 && exterior_turns < escapes, "mixed inward vacuum turn eligibility");
    Check(turn_proposals > 0, "at least one full turn meets the live budget");
    Check(event_error <= 1e-4 && map_error <= 1e-4 && infinity_error <= 1e-4,
          "diagnostic agreement with the complete trace");
}
}  // namespace

int main() {
    CheckBudgetEnclosesOriginalComparison();
    CheckVacuumTurn();
    ScreenOriginalGeometry();
}
