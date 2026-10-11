// CPU accepted-segment event witnesses. These are backend computations only:
// no render session, image dispatch, or presentation path is exercised.

#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera.h"
#include "sirius/core/metrics/cpu_metric_factory.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/metrics/morris_thorne_family.h"
#include "sirius/core/metrics/registry.h"
#include "sirius/render/trace_domain.h"

#include <gtest/gtest.h>

#include "support/ellis_geodesic_reference.h"
#include "support/kottler_geodesic_reference.h"

#include <cmath>
#include <iomanip>
#include <numbers>
#include <sstream>
#include <vector>

namespace {

using sirius::backend::GeodesicTracer;
using sirius::backend::TracerConfig;
using sirius::backend::TraceResult;
using sirius::core::CameraConfig;
using sirius::core::CameraRay;
using sirius::core::KerrSchildFamily;
using sirius::core::KerrSchildParams;
using sirius::core::MetricId;
using sirius::core::PinholeCamera;

constexpr double kLambda = 1.0e-2;

struct BoundaryTrace {
    TraceResult result;
    float boundary_radius;
};

BoundaryTrace TraceDeSitterBoundary(bool enable_bundle, float maximum_step) {
    KerrSchildFamily metric(KerrSchildParams::DeSitter(kLambda));
    const double exact_boundary =
        sirius::core::MetricCosmologicalHorizonRadius(MetricId::DeSitter, 0.0, kLambda).value();
    float stored_boundary = static_cast<float>(exact_boundary);
    if (static_cast<double>(stored_boundary) > exact_boundary) {
        stored_boundary = std::nextafter(stored_boundary, 0.0f);
    }

    TracerConfig config;
    config.escape_radius = stored_boundary;
    config.finite_causal_boundary = true;
    config.max_steps = 4000;
    config.enable_disk = false;
    config.enable_ray_bundles = enable_bundle;
    config.bundle_point_source = enable_bundle;
    config.bundle_angular_size = 1.0e-3f;
    config.integrator.initial_step = maximum_step;
    config.integrator.max_step = maximum_step;
    config.integrator.min_step = 1.0e-5f;
    config.integrator.abs_tolerance = 1.0e-7f;
    config.integrator.rel_tolerance = 1.0e-7f;

    CameraConfig camera_config;
    camera_config.r = 5.0;
    camera_config.theta = std::numbers::pi / 2.0;
    camera_config.phi = 0.0;
    camera_config.fov = 60.0f;
    camera_config.width = 7;
    camera_config.height = 7;
    PinholeCamera camera(camera_config);
    const CameraRay ray = camera.GenerateRay(5, 3, 0.5f, 0.5f);

    GeodesicTracer tracer(&metric, config);
    return {tracer.Trace(ray), stored_boundary};
}

BoundaryTrace TraceFlatBoundaryBundle() {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    constexpr float kBoundaryRadius = 12.0f;

    TracerConfig config;
    config.escape_radius = kBoundaryRadius;
    config.finite_causal_boundary = true;
    config.max_steps = 20;
    config.enable_disk = false;
    config.enable_ray_bundles = true;
    config.bundle_point_source = true;
    config.bundle_angular_size = 1.0e-3f;
    config.integrator.initial_step = 7.0f;
    config.integrator.max_step = 7.0f;
    config.integrator.min_step = 1.0e-5f;
    config.integrator.abs_tolerance = 1.0e-7f;
    config.integrator.rel_tolerance = 1.0e-7f;

    CameraConfig camera_config;
    camera_config.r = 5.0;
    camera_config.theta = std::numbers::pi / 2.0;
    camera_config.phi = 0.0;
    camera_config.fov = 60.0f;
    camera_config.width = 3;
    camera_config.height = 3;
    PinholeCamera camera(camera_config);

    GeodesicTracer tracer(&metric, config);
    return {tracer.Trace(camera.GenerateRay(1, 1, 0.5f, 0.5f)), kBoundaryRadius};
}

TraceResult TraceCentralEllisRay(sirius::core::WormholeTopology topology, double b0,
                                 double observer_radius, double escape_radius, float maximum_step) {
    sirius::core::MorrisThorneCartesian metric(sirius::core::MorrisThorneParams::Ellis(b0));

    TracerConfig config;
    config.escape_radius = static_cast<float>(escape_radius);
    config.max_steps = 20000;
    config.wormhole_topology = topology;
    config.enable_disk = false;
    config.integrator.initial_step = maximum_step;
    config.integrator.max_step = maximum_step;
    config.integrator.min_step = static_cast<float>(1.0e-6 * b0);
    config.integrator.abs_tolerance = 1.0e-7f;
    config.integrator.rel_tolerance = 1.0e-7f;

    CameraConfig camera_config;
    camera_config.r = observer_radius;
    camera_config.theta = std::numbers::pi / 2.0;
    camera_config.phi = 0.0;
    camera_config.fov = 60.0f;
    camera_config.width = 3;
    camera_config.height = 3;
    PinholeCamera camera(camera_config);

    GeodesicTracer tracer(&metric, config);
    return tracer.Trace(camera.GenerateRay(1, 1, 0.5f, 0.5f));
}

double TerminalRadius(const TraceResult& trace) {
    return std::sqrt(trace.final_position(1) * trace.final_position(1) +
                     trace.final_position(2) * trace.final_position(2) +
                     trace.final_position(3) * trace.final_position(3));
}

TracerConfig OrdinaryBoundaryConfig(float maximum_step = 7.0f) {
    TracerConfig config;
    config.escape_radius = 12.0f;
    config.enable_disk = false;
    config.max_steps = 2000;
    config.integrator.initial_step = maximum_step;
    config.integrator.max_step = maximum_step;
    config.integrator.min_step = 1.0e-6f;
    config.integrator.abs_tolerance = 1.0e-7f;
    config.integrator.rel_tolerance = 1.0e-7f;
    return config;
}

TraceResult TraceRadialBoundaryRay(sirius::core::IMetric& metric, const TracerConfig& config,
                                   double launch_radius = 5.0, double radial_direction = -1.0) {
    CameraRay ray;
    ray.origin(1) = launch_radius;
    ray.origin(2) = std::numbers::pi / 2.0;
    ray.direction(1) = radial_direction;
    GeodesicTracer tracer(&metric, config);
    return tracer.Trace(ray);
}

sirius::core::MetricConstructionParameters RepresentedParametersFor(MetricId id) {
    sirius::core::MetricConstructionParameters parameters;
    switch (id) {
        case MetricId::Minkowski:
        case MetricId::DeSitter:
        case MetricId::MorrisThorne:
        case MetricId::Alcubierre:
            parameters.mass = 0.0;
            break;
        case MetricId::Kerr:
            parameters.mass = 2.0;
            parameters.dimensionless_spin = 0.5;
            break;
        case MetricId::ReissnerNordstrom:
            parameters.mass = 2.0;
            parameters.dimensionless_charge = 0.3;
            break;
        case MetricId::KerrNewman:
            parameters.mass = 2.0;
            parameters.dimensionless_spin = 0.3;
            parameters.dimensionless_charge = 0.3;
            break;
        case MetricId::Schwarzschild:
        case MetricId::SchwarzschildDeSitter:
            parameters.mass = 2.0;
            break;
    }
    if (id == MetricId::DeSitter || id == MetricId::SchwarzschildDeSitter) {
        parameters.cosmological_constant = 1.0e-3;
    }
    return parameters;
}

void ExpectConstructedParameters(MetricId id,
                                 const sirius::core::MetricConstructionParameters& expected,
                                 sirius::core::IMetric* metric) {
    switch (id) {
        case MetricId::Minkowski:
        case MetricId::Schwarzschild:
        case MetricId::Kerr:
        case MetricId::ReissnerNordstrom:
        case MetricId::KerrNewman:
        case MetricId::DeSitter:
        case MetricId::SchwarzschildDeSitter: {
            const auto* family = dynamic_cast<const KerrSchildFamily*>(metric);
            ASSERT_NE(family, nullptr);
            const KerrSchildParams actual = family->GetParams();
            EXPECT_DOUBLE_EQ(actual.M, expected.mass);
            EXPECT_DOUBLE_EQ(actual.a, expected.dimensionless_spin * expected.mass);
            EXPECT_DOUBLE_EQ(actual.Q, expected.dimensionless_charge * expected.mass);
            EXPECT_DOUBLE_EQ(actual.Lambda, expected.cosmological_constant);
            break;
        }
        case MetricId::MorrisThorne: {
            const auto* family = dynamic_cast<const sirius::core::MorrisThorneCartesian*>(metric);
            ASSERT_NE(family, nullptr);
            EXPECT_DOUBLE_EQ(family->SphericalFamily().GetParams().b0, expected.throat_radius);
            break;
        }
        case MetricId::Alcubierre: {
            const auto* family = dynamic_cast<const sirius::core::WarpDriveFamily*>(metric);
            ASSERT_NE(family, nullptr);
            const sirius::core::WarpDriveParams actual = family->GetParams();
            EXPECT_DOUBLE_EQ(actual.vs, expected.warp_velocity);
            EXPECT_DOUBLE_EQ(actual.R, expected.bubble_radius);
            EXPECT_DOUBLE_EQ(actual.sigma, expected.bubble_sigma);
            break;
        }
    }
}

// Controlled executor samples exercise coordinator ownership, not geodesic
// accuracy. Position locators and physical tangents are deliberately distinct.
struct QuarticCoordinatorExecutor final : sirius::backend::TraceStepExecutor {
    struct VolumeSample {
        double fraction;
        sirius::core::Vec4 position;
    };
    bool fail_volume = false;
    int rejected = 0, ended = 0;
    std::vector<VolumeSample> volume_samples;
    std::array<sirius::core::DopriPositionSegment, sirius::core::kCoupledTrialCount> curves;
    sirius::core::Lightray origin{};
    sirius::core::GeodesicVariations initial;
    sirius::core::IMetric* metric = nullptr;
    double interval = 0.0;

    void BeginTrace() override { volume_samples.clear(); }
    void EndTrace() override { ++ended; }
    void RejectLastInterval() override { ++rejected; }
    std::optional<sirius::core::CameraLaunch> Launch(sirius::core::IMetric& authority, double spin,
                                                     const CameraRay& camera) override {
        return sirius::core::LaunchCameraRay(authority, spin, camera);
    }
    sirius::core::CoupledSegmentSample At(sirius::core::CoupledTrial trial, double fraction,
                                          const sirius::core::Vec4* normal = nullptr) {
        const auto index = static_cast<std::size_t>(trial);
        const auto locator = curves[index].Sample(fraction);
        const double global_fraction =
            index < 2 ? fraction : (fraction + (index == 3 ? 1.0 : 0.0)) * 0.5;
        sirius::core::CoupledSegmentSample sample;
        sample.ray = origin;
        sample.ray.position = locator.position;
        // Unit outward radial affine velocity, with its past null time
        // component derived from the current metric at this supplied event.
        sirius::core::Metric4d values;
        sirius::core::Tensor<sirius::core::Dual<double>, 4, 4, 4> derivatives;
        metric->Evaluate(sample.ray.position, values, derivatives);
        sample.ray.velocity = {};
        sample.ray.velocity(1) = 1.0;
        const double g00 = values(0, 0).real, g01 = values(0, 1).real, g11 = values(1, 1).real;
        sample.ray.velocity(0) = (-g01 + std::sqrt(g01 * g01 - g00 * g11)) / g00;
        sample.polynomial_tangent = locator.tangent;
        sample.variations = initial;
        for (auto& column : sample.variations) {
            column.displacement += column.derivative * (interval * global_fraction);
            if (normal) {
                double numerator = 0.0, denominator = 0.0;
                for (int axis = 0; axis < 4; ++axis) {
                    numerator += (*normal)(axis)*column.displacement(axis);
                    denominator += (*normal)(axis)*sample.ray.velocity(axis);
                }
                column.displacement -= sample.ray.velocity * (numerator / denominator);
            }
        }
        return sample;
    }
    bool Step(sirius::core::Lightray& ray, sirius::core::IMetric& authority,
              const sirius::core::IntegratorConfig&, sirius::core::Rk45CoupledState& coupled,
              sirius::core::Rk45CoupledComparison& comparison) override {
        metric = &authority;
        origin = ray;
        initial = coupled.variations;
        interval = ray.step_size;
        // x(s)=x0+40s(1-s): both full endpoints lie inside R=12,
        // while the interior crosses that boundary and the later disk annulus.
        curves[0] = {};
        curves[0].origin = origin.position;
        curves[0].increment(0) = -interval;
        curves[0].a(1) = 40.0;
        curves[0].interval = interval;
        curves[1] = curves[0];
        curves[2] = {};
        curves[2].origin = origin.position;
        curves[2].increment(0) = -interval * 0.5;
        curves[2].increment(1) = curves[2].a(1) = 10.0;
        curves[2].interval = interval * 0.5;
        curves[3] = curves[2];
        curves[3].origin = curves[0].Sample(0.5).position;
        curves[3].increment(1) = -10.0;
        for (std::size_t trial = 0; trial < curves.size(); ++trial)
            comparison.dopri_positions[trial] = curves[trial];
        const auto full = At(sirius::core::CoupledTrial::Full, 1.0);
        const auto half = At(sirius::core::CoupledTrial::FirstHalf, 1.0);
        ray = comparison.lower_order = comparison.refined_endpoint = full.ray;
        coupled.variations = comparison.lower_variations = comparison.refined_variations =
            full.variations;
        comparison.midpoint = half.ray;
        comparison.midpoint_variations = half.variations;
        return true;
    }
    std::optional<sirius::core::CoupledSegmentSample> Sample(
        sirius::core::CoupledTrial trial, double fraction,
        const sirius::core::Vec4* normal = nullptr) override {
        const auto sample = At(trial, fraction, normal);
        if (trial == sirius::core::CoupledTrial::Full && !normal && fraction > 0.0 &&
            fraction < 1.0) {
            volume_samples.push_back({fraction, sample.ray.position});
            if (fail_volume && volume_samples.size() == 3) return std::nullopt;
        }
        return sample;
    }
};

TEST(CpuTraceBoundary, RetainedQuarticEventsKeepOriginalFractionsAndExcludeLaterDisk) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    auto config = OrdinaryBoundaryConfig(10.0f);
    config.finite_causal_boundary = true;
    config.enable_disk = true;
    config.disk_inner = 13.0;
    config.disk_outer = 20.0;
    config.max_steps = 1;
    CameraRay camera;
    camera.origin(1) = 5.0;
    camera.origin(2) = std::numbers::pi / 2.0;
    camera.direction(1) = 1.0;
    QuarticCoordinatorExecutor executor;
    GeodesicTracer tracer(&metric, config);
    tracer.SetStepExecutor(&executor);
    const auto result = tracer.Trace(camera);
    ASSERT_FALSE(result.numerical_failure);
    ASSERT_EQ(result.outcome, TraceResult::Outcome::Escaped);
    EXPECT_NEAR(TerminalRadius(result), 12.0, 2.0e-10);
    EXPECT_NEAR(result.affine_length, 5.0 * (1.0 - std::sqrt(0.3)), 2.0e-10);
    EXPECT_EQ(result.num_disk_crossings, 0);
    ASSERT_TRUE(result.final_tangent);
    EXPECT_NEAR((*result.final_tangent)(1), 1.0, 1.0e-12);
    EXPECT_GT(executor.curves[0].Sample(result.affine_length / 10.0).tangent(1), 2.0);
    EXPECT_EQ(executor.rejected, 0);
}

TEST(CpuTraceBoundary, RetainedQuarticVolumeFailurePublishesNothingAndRecovers) {
    // A charged metric avoids an infinity-handoff requirement in this finite
    // coordinator witness. These supplied samples do not qualify charged flow.
    KerrSchildFamily metric(KerrSchildParams::ReissnerNordstrom(1.0, 0.5));
    auto config = OrdinaryBoundaryConfig(10.0f);
    config.finite_causal_boundary = true;
    config.enable_disk = config.enable_volumetric = true;
    config.volumetric_samples = 4;
    config.volumetric_tau_midplane = 0.01f;
    config.max_steps = 1;
    CameraRay camera;
    camera.origin(1) = 8.0;
    camera.origin(2) = std::numbers::pi / 2.0;
    camera.direction(1) = 1.0;
    QuarticCoordinatorExecutor executor;
    GeodesicTracer tracer(&metric, config);
    tracer.SetStepExecutor(&executor);
    executor.fail_volume = true;
    const auto failed = tracer.Trace(camera);
    EXPECT_TRUE(failed.numerical_failure);
    ASSERT_EQ(executor.volume_samples.size(), 3);
    EXPECT_FALSE(failed.volumetric_hit);
    EXPECT_FLOAT_EQ(failed.optical_depth, 0.0f);
    EXPECT_FLOAT_EQ(failed.volumetric_affine_length, 0.0f);
    for (float channel : failed.volumetric_emission) EXPECT_FLOAT_EQ(channel, 0.0f);
    EXPECT_DOUBLE_EQ(failed.affine_length, 0.0);
    EXPECT_GE(executor.rejected, 1);

    executor.fail_volume = false;
    const auto recovered = tracer.Trace(camera);
    ASSERT_FALSE(recovered.numerical_failure);
    ASSERT_EQ(recovered.outcome, TraceResult::Outcome::Escaped);
    ASSERT_EQ(executor.volume_samples.size(), 4);
    ASSERT_TRUE(recovered.volumetric_hit);
    EXPECT_GT(recovered.optical_depth, 0.0f);
    EXPECT_GT(recovered.volumetric_affine_length, 0.0f);
    const double limit = 0.5 * (1.0 - std::sqrt(0.6));
    for (std::size_t index = 0; index < executor.volume_samples.size(); ++index) {
        const double fraction = limit * (static_cast<double>(index) + 0.5) / 4.0;
        EXPECT_NEAR(executor.volume_samples[index].fraction, fraction, 2.0e-12);
        EXPECT_NEAR(executor.volume_samples[index].position(1),
                    8.0 + 40.0 * fraction * (1.0 - fraction), 2.0e-10);
        EXPECT_LT(executor.volume_samples[index].position(1), 12.0);
    }
    EXPECT_NEAR(recovered.affine_length, 10.0 * limit, 2.0e-10);
    EXPECT_EQ(executor.ended, 2);
}

TEST(CpuTraceBoundary, EveryAdvertisedCpuMetricConstructsAndTracesOneRay) {
    auto absent_parameter = RepresentedParametersFor(MetricId::Minkowski);
    absent_parameter.dimensionless_spin = 0.5;
    EXPECT_EQ(sirius::core::CreateCpuMetric(MetricId::Minkowski, absent_parameter), nullptr);

    auto two_sheet_topology = RepresentedParametersFor(MetricId::MorrisThorne);
    two_sheet_topology.wormhole_topology = sirius::core::WormholeTopology::TwoSheet;
    EXPECT_NE(sirius::core::CreateCpuMetric(MetricId::MorrisThorne, two_sheet_topology), nullptr);

    std::size_t advertised_cpu_metrics = 0;
    for (const auto& info : sirius::core::MetricRegistry()) {
        if (!info.cpu_supported) continue;
        ++advertised_cpu_metrics;
        SCOPED_TRACE(info.canonical_name);

        const auto parameters = RepresentedParametersFor(info.id);
        auto metric = sirius::core::CreateCpuMetric(info.id, parameters);
        ASSERT_NE(metric, nullptr);
        EXPECT_EQ(sirius::core::ParseMetricName(metric->GetName()), info.id);
        ExpectConstructedParameters(info.id, parameters, metric.get());

        const double scale =
            sirius::core::MetricSceneLengthScale(info.id, parameters.mass, parameters.throat_radius,
                                                 parameters.bubble_radius, parameters.bubble_sigma);
        CameraConfig camera_config;
        camera_config.r = 10.0 * scale;
        camera_config.theta = std::numbers::pi / 2.0;
        camera_config.phi = 0.0;
        camera_config.fov = 60.0f;
        camera_config.width = 3;
        camera_config.height = 3;
        PinholeCamera camera(camera_config);

        TracerConfig config;
        config.enable_disk = false;
        config.escape_radius = static_cast<float>(100.0 * scale);
        config.max_steps = 5000;
        config.integrator.initial_step = 0.1f;
        config.integrator.max_step = 0.25f;
        if (const auto horizon = sirius::core::MetricCosmologicalHorizonRadius(
                info.id, parameters.mass, parameters.cosmological_constant);
            horizon.has_value()) {
            const double interior_boundary =
                sirius::core::kMaxCosmologicalObserverFraction * *horizon;
            config.escape_radius = static_cast<float>(interior_boundary);
            if (static_cast<double>(config.escape_radius) > interior_boundary) {
                config.escape_radius = std::nextafter(config.escape_radius, 0.0f);
            }
            config.finite_causal_boundary = true;
        }

        GeodesicTracer tracer(metric.get(), config);
        CameraRay ray = camera.GenerateRay(1, 1, 0.5f, 0.5f);
        ray.direction(1) = 1.0;
        ray.direction(2) = 0.0;
        ray.direction(3) = 0.0;
        const TraceResult result = tracer.Trace(ray);
        EXPECT_FALSE(result.numerical_failure);
        EXPECT_NE(result.outcome, TraceResult::Outcome::MaxSteps);
        EXPECT_GT(result.steps_taken, 0);
        EXPECT_GT(result.affine_length, 0.0f);
        for (int component = 0; component < 4; ++component) {
            EXPECT_TRUE(std::isfinite(result.final_position(component)));
        }
    }
    EXPECT_EQ(advertised_cpu_metrics, 9u);
}

TEST(CpuTraceBoundary, FinitePupilOffsetMovesTheLiveCpuLaunchEvent) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());

    TracerConfig config;
    config.escape_radius = 12.0f;
    config.finite_causal_boundary = true;
    config.max_steps = 20;
    config.enable_disk = false;
    config.integrator.initial_step = 7.0f;
    config.integrator.max_step = 7.0f;
    config.integrator.min_step = 1.0e-5f;
    config.integrator.abs_tolerance = 1.0e-7f;
    config.integrator.rel_tolerance = 1.0e-7f;

    CameraConfig camera_config;
    camera_config.r = 5.0;
    camera_config.theta = std::numbers::pi / 2.0;
    camera_config.phi = 0.0;
    camera_config.width = 3;
    camera_config.height = 3;
    PinholeCamera camera(camera_config);

    CameraRay central = camera.GenerateRay(1, 1, 0.5f, 0.5f);
    central.direction(1) = 1.0;
    central.direction(2) = 0.0;
    central.direction(3) = 0.0;
    CameraRay pupil = central;
    pupil.aperture_right = 0.25;

    GeodesicTracer tracer(&metric, config);
    const TraceResult central_trace = tracer.Trace(central);
    const TraceResult pupil_trace = tracer.Trace(pupil);

    ASSERT_EQ(central_trace.outcome, TraceResult::Outcome::Escaped)
        << "steps=" << central_trace.steps_taken << " numerical=" << central_trace.numerical_failure
        << " coupled=" << static_cast<int>(central_trace.coupled_failure)
        << " x=" << central_trace.final_position(1);
    ASSERT_EQ(pupil_trace.outcome, TraceResult::Outcome::Escaped)
        << "steps=" << pupil_trace.steps_taken << " numerical=" << pupil_trace.numerical_failure
        << " coupled=" << static_cast<int>(pupil_trace.coupled_failure)
        << " x=" << pupil_trace.final_position(1);
    // The central interval ends exactly on the sphere. That endpoint owns the
    // escape; advancing another interval would start outside the locator's domain.
    EXPECT_EQ(central_trace.steps_taken, 1);
    EXPECT_FALSE(central_trace.numerical_failure);
    EXPECT_FALSE(pupil_trace.numerical_failure);
    EXPECT_NEAR(central_trace.final_position(2), 0.0, 2.0e-5);
    EXPECT_NEAR(pupil_trace.final_position(2), pupil.aperture_right, 2.0e-5);
    EXPECT_NEAR(pupil_trace.final_position(2) - central_trace.final_position(2),
                pupil.aperture_right, 2.0e-5)
        << "the live tracer ignored the finite-pupil launch-event displacement";
}

TEST(CpuTraceBoundary, CancellationDiscardsPrivateRayDataAndAllowsTracerReuse) {
    struct CancellingExecutor final : sirius::backend::TraceStepExecutor {
        bool cancel_after_step = true, cancelled = false;
        int accepted = 0, rejected = 0, ended = 0;
        std::optional<sirius::core::CameraLaunch> Launch(sirius::core::IMetric& metric, double spin,
                                                         const CameraRay& camera) override {
            return sirius::core::LaunchCameraRay(metric, spin, camera);
        }
        bool Step(sirius::core::Lightray& ray, sirius::core::IMetric& metric,
                  const sirius::core::IntegratorConfig& config,
                  sirius::core::Rk45CoupledState& coupled,
                  sirius::core::Rk45CoupledComparison& comparison) override {
            const bool success = sirius::core::Geodesic::IntegrateStepRk45(ray, &metric, config,
                                                                           &coupled, &comparison);
            if (success) {
                ++accepted;
                if (cancel_after_step) cancelled = true;
            }
            return success;
        }
        void RejectLastInterval() override { ++rejected; }
        void EndTrace() override { ++ended; }
    };
    TracerConfig config;
    config.escape_radius = 7;
    config.finite_causal_boundary = true;
    config.max_steps = 20;
    config.enable_disk = false;
    config.integrator.initial_step = config.integrator.max_step = 1;
    CameraConfig camera_config;
    camera_config.r = 5;
    camera_config.theta = std::numbers::pi / 2;
    camera_config.width = camera_config.height = 3;
    PinholeCamera camera(camera_config);
    auto ray = camera.GenerateRay(1, 1, .5f, .5f);
    ray.direction(1) = 1;
    ray.direction(2) = ray.direction(3) = 0;
    for (const auto parameters :
         {KerrSchildParams::Minkowski(), KerrSchildParams::Schwarzschild(1)}) {
        SCOPED_TRACE(parameters.M);
        KerrSchildFamily metric(parameters);
        CancellingExecutor executor;
        GeodesicTracer tracer(&metric, config);
        tracer.SetStepExecutor(&executor);
        tracer.SetCancellationCallback([&] { return executor.cancelled; });
        const auto stopped = tracer.Trace(ray);
        EXPECT_TRUE(stopped.cancelled);
        EXPECT_EQ(executor.accepted, 1);
        EXPECT_GE(executor.rejected, 1);
        EXPECT_EQ(executor.ended, 1);
        EXPECT_FALSE(stopped.numerical_failure);
        EXPECT_FALSE(stopped.final_tangent);
        EXPECT_FALSE(stopped.beam.valid);
        EXPECT_EQ(stopped.num_disk_crossings, 0);

        executor.cancel_after_step = executor.cancelled = false;
        const auto complete = tracer.Trace(ray);
        EXPECT_FALSE(complete.cancelled);
        EXPECT_FALSE(complete.numerical_failure);
        EXPECT_EQ(complete.outcome, TraceResult::Outcome::Escaped);
        EXPECT_EQ(executor.ended, 2);
    }
}

TEST(CpuTraceBoundary, TruncatedPageThorneLiveProfileUsesDeclaredZeroTorqueEdge) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1.0));

    TracerConfig config;
    config.escape_radius = 100.0f;
    config.horizon_factor = 1.05f;
    config.max_steps = 5000;
    config.enable_disk = true;
    config.disk_inner = 6.0f;
    config.disk_outer = 20.0f;
    config.disk_temperature_inner = 10000.0f;
    config.integrator.abs_tolerance = 1.0e-6f;
    config.integrator.rel_tolerance = 1.0e-6f;
    GeodesicTracer tracer(&metric, config);

    CameraConfig camera_config;
    camera_config.r = 50.0;
    camera_config.theta = std::numbers::pi / 2.0;
    camera_config.phi = 0.0;
    camera_config.fov = 60.0f;
    camera_config.width = 64;
    camera_config.height = 64;
    PinholeCamera camera(camera_config);

    // Populate the live tracer's cached Page-Thorne authority before changing
    // the declared edge. SetConfig must not retain the former ISCO profile.
    (void)tracer.Trace(camera.GenerateRay(32, 32, 0.5f, 0.5f));
    config.disk_inner = 10.0f;
    tracer.SetConfig(config);

    sirius::core::AccretionDiskD::Config oracle_config;
    oracle_config.M = 1.0;
    oracle_config.a_star = 0.0;
    oracle_config.r_inner = config.disk_inner;
    oracle_config.r_outer = 100.0;
    sirius::core::AccretionDiskD oracle(oracle_config);
    const double reference = oracle.Temperature(
        sirius::core::constants::disk::kTemperatureReferenceRadiusRatio * oracle.InnerRadius());
    ASSERT_GT(reference, 0.0);

    int compared = 0;
    for (int y = 8; y < 56 && compared < 12; y += 2) {
        for (int x = 8; x < 56 && compared < 12; x += 2) {
            const TraceResult result = tracer.Trace(camera.GenerateRay(x, y, 0.5f, 0.5f));
            for (int crossing_index = 0; crossing_index < result.num_disk_crossings;
                 ++crossing_index) {
                const auto& crossing = result.disk_crossings[crossing_index];
                if (!crossing.valid) continue;
                const double expected =
                    config.disk_temperature_inner * oracle.Temperature(crossing.r) / reference;
                EXPECT_NEAR(crossing.temperature, expected, 2.0e-3)
                    << "live crossing at r=" << crossing.r
                    << " did not consume the declared truncated Page-Thorne edge";
                ++compared;
            }
        }
    }
    EXPECT_GE(compared, 6) << "insufficient live disk crossings for the truncated Page-Thorne gate";
}

TEST(CpuTraceBoundary, CentralEventIsInvariantUnderBundleFeatureToggle) {
    const BoundaryTrace without_bundle = TraceDeSitterBoundary(false, 4.0f);
    const BoundaryTrace with_bundle = TraceDeSitterBoundary(true, 4.0f);

    ASSERT_EQ(without_bundle.result.outcome, TraceResult::Outcome::Escaped);
    ASSERT_EQ(with_bundle.result.outcome, TraceResult::Outcome::Escaped);
    ASSERT_FALSE(without_bundle.result.numerical_failure);
    ASSERT_FALSE(with_bundle.result.numerical_failure);
    EXPECT_NEAR(TerminalRadius(without_bundle.result), without_bundle.boundary_radius, 2.0e-5);
    EXPECT_NEAR(TerminalRadius(with_bundle.result), with_bundle.boundary_radius, 2.0e-5);

    for (int component = 0; component < 4; ++component) {
        EXPECT_NEAR(without_bundle.result.final_position(component),
                    with_bundle.result.final_position(component), 1.0e-10);
        EXPECT_NEAR(without_bundle.result.final_direction(component),
                    with_bundle.result.final_direction(component), 1.0e-10);
    }
}

TEST(CpuTraceBoundary, JacobiBundleTerminatesAtTheSameCausalEvent) {
    const BoundaryTrace trace = TraceFlatBoundaryBundle();

    ASSERT_EQ(trace.result.outcome, TraceResult::Outcome::Escaped);
    ASSERT_FALSE(trace.result.numerical_failure);
    ASSERT_TRUE(trace.result.beam.valid);
    EXPECT_NEAR(TerminalRadius(trace.result), trace.boundary_radius, 2.0e-6);

    // In Minkowski space the point-source Jacobi map is xi=lambda*xi_dot.
    // The centre ray travels 5 units to the origin and 12 to the far boundary,
    // so a 1e-3 angular seed has a 17e-3 transverse semi-axis there. Advancing
    // through the deliberately seven-unit overshoot would instead report 21e-3.
    constexpr double kExpectedSemiAxis = 17.0e-3;
    constexpr double kExpectedFootprint = kExpectedSemiAxis / 12.0;
    EXPECT_NEAR(trace.result.beam.semi_major, kExpectedSemiAxis, 2.0e-7);
    EXPECT_NEAR(trace.result.beam.semi_minor, kExpectedSemiAxis, 2.0e-7);
    EXPECT_NEAR(trace.result.beam.footprint_major, kExpectedFootprint, 2.0e-8);
    EXPECT_NEAR(trace.result.beam.footprint_minor, kExpectedFootprint, 2.0e-8);
}

TEST(CpuTraceBoundary, OrdinaryEscapeClipsCentralAndJacobiAffineIntervals) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    // Both a multi-step overshoot and a single segment through the origin
    // must select the first outward R=12 event, at lambda=5+12.
    for (float maximum_step : {7.0f, 17.0f, 30.0f}) {
        for (bool bundle : {false, true}) {
            SCOPED_TRACE(maximum_step);
            SCOPED_TRACE(bundle);
            auto config = OrdinaryBoundaryConfig(maximum_step);
            config.enable_ray_bundles = bundle;
            config.bundle_point_source = bundle;
            config.bundle_angular_size = bundle ? 1.0e-3f : config.bundle_angular_size;
            const auto result = TraceRadialBoundaryRay(metric, config);
            ASSERT_EQ(result.outcome, TraceResult::Outcome::Escaped);
            ASSERT_FALSE(result.numerical_failure);
            EXPECT_NEAR(TerminalRadius(result), 12.0, 2.0e-9);
            EXPECT_NEAR(result.affine_length, 17.0, 2.0e-9);
            EXPECT_NEAR(result.final_position(0), -17.0, 2.0e-9);
            EXPECT_NEAR(result.final_position(1), -12.0, 2.0e-9);
            EXPECT_NEAR(result.final_direction(1), -1.0, 2.0e-9);
            if (bundle) {
                ASSERT_TRUE(result.beam.valid);
                EXPECT_NEAR(result.beam.semi_major, 17.0e-3, 2.0e-7);
                EXPECT_NEAR(result.beam.semi_minor, 17.0e-3, 2.0e-7);
            }
        }
    }
}

TEST(CpuTraceBoundary, OrdinaryEscapeBoundaryLaunchUsesCrossingDirection) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    const auto config = OrdinaryBoundaryConfig();
    for (double direction : {-1.0, 1.0}) {
        const auto result = TraceRadialBoundaryRay(metric, config, 12.0, direction);
        ASSERT_EQ(result.outcome, TraceResult::Outcome::Escaped);
        ASSERT_FALSE(result.numerical_failure);
        EXPECT_NEAR(TerminalRadius(result), 12.0, 2.0e-9);
        EXPECT_NEAR(result.affine_length, direction > 0.0 ? 0.0 : 24.0, 2.0e-9);
        EXPECT_NEAR(result.final_position(1), direction * 12.0, 2.0e-9);
    }
}

TEST(CpuTraceBoundary, OrdinaryEscapeDeclinesExteriorOutwardLaunchWithoutSourceEvent) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    const auto config = OrdinaryBoundaryConfig();
    const auto outward = TraceRadialBoundaryRay(metric, config, 15.0, 1.0);
    EXPECT_EQ(outward.outcome, TraceResult::Outcome::MaxSteps);
    EXPECT_TRUE(outward.numerical_failure);
    EXPECT_DOUBLE_EQ(outward.affine_length, 0.0);
    EXPECT_NEAR(TerminalRadius(outward), 15.0, 2.0e-9);
    EXPECT_DOUBLE_EQ(outward.final_position(0), 0.0);

    // Exterior origin alone is not a refusal: an inward ray can traverse
    // the sphere and reach its future outward boundary at lambda=15+12.
    const auto inward = TraceRadialBoundaryRay(metric, config, 15.0, -1.0);
    ASSERT_EQ(inward.outcome, TraceResult::Outcome::Escaped);
    EXPECT_FALSE(inward.numerical_failure);
    EXPECT_NEAR(inward.affine_length, 27.0, 2.0e-9);
    EXPECT_NEAR(TerminalRadius(inward), 12.0, 2.0e-9);
}

TEST(CpuTraceBoundary, OrdinaryEscapeExcludesDiskAndVolumeBeyondBoundary) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1.0));
    for (bool volume : {true, false}) {
        auto config = OrdinaryBoundaryConfig();
        config.enable_disk = true;
        config.enable_volumetric = volume;
        config.enable_polarisation = !volume;
        config.disk_inner = 12.00001;
        config.disk_outer = 20.0;
        const auto result = TraceRadialBoundaryRay(metric, config, 5.0, 1.0);
        ASSERT_EQ(result.outcome, TraceResult::Outcome::Escaped);
        ASSERT_FALSE(result.numerical_failure);
        EXPECT_NEAR(TerminalRadius(result), 12.0, 2.0e-8);
        EXPECT_EQ(result.num_disk_crossings, 0);
        EXPECT_FALSE(result.volumetric_hit);
        EXPECT_FLOAT_EQ(result.optical_depth, 0.0f);
        EXPECT_FLOAT_EQ(result.volumetric_affine_length, 0.0f);
        for (float channel : result.volumetric_emission) EXPECT_FLOAT_EQ(channel, 0.0f);

        // The same consumer is active when its source is inside the boundary.
        config.disk_inner = 8.0;
        const auto inside = TraceRadialBoundaryRay(metric, config, 5.0, 1.0);
        ASSERT_FALSE(inside.numerical_failure);
        if (volume) {
            ASSERT_EQ(inside.outcome, TraceResult::Outcome::Escaped);
            EXPECT_TRUE(inside.volumetric_hit);
            EXPECT_GT(inside.optical_depth, 0.0f);
            EXPECT_GT(inside.volumetric_affine_length, 0.0f);
            EXPECT_NEAR(TerminalRadius(inside), 12.0, 2.0e-8);
        } else {
            ASSERT_EQ(inside.outcome, TraceResult::Outcome::DiskHit);
            ASSERT_EQ(inside.num_disk_crossings, 1);
            EXPECT_NEAR(TerminalRadius(inside), 8.0, 2.0e-6);
            EXPECT_TRUE(inside.disk_crossings[0].polarisation_valid);
            config.escape_radius = 40.0f;
            const auto farther = TraceRadialBoundaryRay(metric, config, 5.0, 1.0);
            ASSERT_EQ(farther.outcome, TraceResult::Outcome::DiskHit);
            EXPECT_NEAR(inside.affine_length, farther.affine_length, 2.0e-8);
            EXPECT_NEAR(inside.disk_crossings[0].polarisation_evpa,
                        farther.disk_crossings[0].polarisation_evpa, 2.0e-7);
        }
    }
}

TEST(CpuTraceBoundary, OneSheetEllisNamesTheRegularThroatBoundary) {
    constexpr double kB0 = 1.0;
    const TraceResult trace = TraceCentralEllisRay(sirius::core::WormholeTopology::OneSheetCapture,
                                                   kB0, 10.0, 20.0, 0.25f);

    ASSERT_EQ(trace.outcome, TraceResult::Outcome::Throat);
    ASSERT_FALSE(trace.numerical_failure);
    EXPECT_EQ(trace.asymptotic_sheet, TraceResult::AsymptoticSheet::Observer);
    EXPECT_NEAR(TerminalRadius(trace), 0.5 * kB0, 2.0e-5);
}

TEST(CpuTraceBoundary, TwoSheetEllisCrossesThroatAndReachesInversionMatchedInfinity) {
    for (const double b0 : {0.25, 1.0, 10.0}) {
        const double observer_radius = 10.0 * b0;
        const double escape_radius = 20.0 * b0;
        const double opposite_radius = b0 * b0 / (4.0 * escape_radius);
        const TraceResult trace =
            TraceCentralEllisRay(sirius::core::WormholeTopology::TwoSheet, b0, observer_radius,
                                 escape_radius, static_cast<float>(0.25 * b0));

        SCOPED_TRACE(b0);
        ASSERT_EQ(trace.outcome, TraceResult::Outcome::Escaped);
        ASSERT_FALSE(trace.numerical_failure);
        EXPECT_EQ(trace.asymptotic_sheet, TraceResult::AsymptoticSheet::Opposite);
        EXPECT_LT(trace.min_radius, 0.5 * b0);
        EXPECT_NEAR(TerminalRadius(trace), opposite_radius, 5.0e-5 * b0);
        EXPECT_GT(trace.final_direction(1), 0.999);
        EXPECT_NEAR(trace.final_direction(2), 0.0, 2.0e-5);
        EXPECT_NEAR(trace.final_direction(3), 0.0, 2.0e-5);

        const double launch_l =
            sirius::core::EllisProperRadialDistanceFromIsotropic(b0, observer_radius);
        const double terminal_l =
            sirius::core::EllisProperRadialDistanceFromIsotropic(b0, opposite_radius);
        EXPECT_NEAR(trace.affine_length, launch_l - terminal_l, 3.0e-3 * b0);
    }
}

// Keep reference errors in long double; EXPECT_NEAR takes binary64 arguments.
TEST(CpuTraceBoundary, EllisReferenceRetainsExactRadialLengthAndIndependentRefinement) {
    namespace reference = sirius::test::ellis_reference;
    for (const long double terminal_rho : {0.5L, 1.0L / 48}) {
        const auto radial = reference::Orbit(8, terminal_rho, 0, terminal_rho < 0.5L, 8192);
        const long double exact = 8 - 1.0L / 32 - terminal_rho + 1 / (4 * terminal_rho);
        EXPECT_LE(std::abs(radial.affine - exact), 1.0e-8L);
        EXPECT_EQ(radial.angle, 0);
        EXPECT_LE(std::abs(radial.sky[1] - (terminal_rho < 0.5L ? 1 : -1)), 1.0e-14L);
    }
    EXPECT_THROW(reference::Orbit(8, 12, 1, false, 8192), std::invalid_argument);
}

TEST(CpuTraceBoundary, NonradialEllisFiniteEventsMatchIndependentQuadrature) {
    namespace reference = sirius::test::ellis_reference;
    using sirius::core::WormholeTopology;
    constexpr double kEnvelope = 1.0e-4;
    for (const double b0 : {0.1, 1.0, 1000.0}) {
        sirius::core::MetricConstructionParameters parameters;
        parameters.mass = 0;
        parameters.throat_radius = b0;
        auto metric = sirius::core::CreateCpuMetric(MetricId::MorrisThorne, parameters);
        ASSERT_NE(metric, nullptr);
        for (const double impact : {0.5, -0.95, 1.05, -2.0}) {
            for (const auto topology :
                 {WormholeTopology::OneSheetCapture, WormholeTopology::TwoSheet}) {
                const bool reflection = std::abs(impact) > 1;
                const bool opposite = !reflection && topology == WormholeTopology::TwoSheet;
                const float stored_boundary = static_cast<float>(12 * b0);
                const long double boundary = static_cast<long double>(stored_boundary) / b0;
                const long double terminal_rho =
                    reflection ? boundary : (opposite ? 1 / (4 * boundary) : 0.5L);
                const auto coarse = reference::Orbit(8, terminal_rho, impact, opposite, 4096);
                const auto expected = reference::Orbit(8, terminal_rho, impact, opposite, 8192);
                ASSERT_LT(std::abs(coarse.angle - expected.angle), 1.0e-6L);
                ASSERT_LT(std::abs(coarse.affine - expected.affine), 1.0e-6L);
                for (int axis = 0; axis < 4; ++axis) {
                    ASSERT_LT(std::abs(coarse.position[axis] - expected.position[axis]), 1.0e-6L);
                    ASSERT_LT(std::abs(coarse.tangent[axis] - expected.tangent[axis]), 1.0e-6L);
                    ASSERT_LT(std::abs(coarse.sky[axis] - expected.sky[axis]), 1.0e-6L);
                }
                for (const double maximum_step : {0.25, 0.125, 0.0625}) {
                    SCOPED_TRACE(::testing::Message()
                                 << "b0=" << b0 << " J/b0=" << impact << " opposite=" << opposite
                                 << " topology=" << static_cast<int>(topology)
                                 << " max_step/b0=" << maximum_step);
                    TracerConfig config;
                    config.escape_radius = stored_boundary;
                    config.horizon_factor = 1;
                    config.wormhole_topology = topology;
                    config.enable_disk = false;
                    config.max_steps = 20000;
                    config.integrator.initial_step = static_cast<float>(maximum_step * b0);
                    config.integrator.max_step = config.integrator.initial_step;
                    config.integrator.min_step = static_cast<float>(1.0e-6 * b0);
                    config.integrator.abs_tolerance = 1.0e-7f;
                    config.integrator.rel_tolerance = 1.0e-7f;
                    CameraRay ray;
                    ray.origin(1) = 8 * b0;
                    ray.origin(2) = std::numbers::pi / 2;
                    ray.direction(3) = impact / (8 + 1.0 / 32);
                    ray.direction(1) = -std::sqrt(1 - ray.direction(3) * ray.direction(3));
                    GeodesicTracer tracer(metric.get(), config);
                    const auto actual = tracer.Trace(ray);
                    ASSERT_EQ(actual.outcome, reflection || opposite
                                                  ? TraceResult::Outcome::Escaped
                                                  : TraceResult::Outcome::Throat);
                    ASSERT_FALSE(actual.numerical_failure);
                    ASSERT_FALSE(actual.cancelled);
                    ASSERT_EQ(actual.integrator_termination, 0);
                    ASSERT_EQ(actual.terminal_chart, TraceResult::TerminalChart::MetricNative);
                    ASSERT_EQ(actual.asymptotic_sheet,
                              opposite ? TraceResult::AsymptoticSheet::Opposite
                                       : TraceResult::AsymptoticSheet::Observer);
                    ASSERT_TRUE(actual.final_tangent.has_value());
                    EXPECT_LE(std::abs(actual.affine_length / b0 - expected.affine), kEnvelope);
                    for (int axis = 0; axis < 4; ++axis) {
                        EXPECT_LE(
                            std::abs(actual.final_position(axis) / b0 - expected.position[axis]),
                            kEnvelope);
                        EXPECT_LE(std::abs((*actual.final_tangent)(axis)-expected.tangent[axis]),
                                  kEnvelope);
                        if (reflection || opposite) {
                            EXPECT_LE(std::abs(actual.final_direction(axis) - expected.sky[axis]),
                                      kEnvelope);
                        }
                    }
                }
            }
        }
    }
}

TEST(CpuTraceBoundary, KottlerReferenceRetainsExactDeSitterAndRadialIdentities) {
    namespace reference = sirius::test::kottler_reference;
    // Exact de Sitter primitives: A^2=C^2+Lambda*L^2/3,
    // affine primitive sqrt(A^2*r^2-L^2)/A^2,
    // angle primitive sign(L)*acos(|L|/(A*r)). Neither uses quadrature.
    for (const auto input : {std::array<long double, 3>{5, -0.6L, 0.8L},
                             std::array<long double, 3>{15, 0.6L, -0.8L}}) {
        const long double r0 = input[0], nr = input[1], nphi = input[2];
        const long double boundary = 17;
        const long double h = 0.01L * r0 * r0 / 3;
        const long double c = (1 + h * nr) / std::sqrt(1 + h), l = r0 * nphi;
        const long double a2 = c * c + 0.01L * l * l / 3;
        const auto length = [=](long double r) { return std::sqrt(a2 * r * r - l * l) / a2; };
        const auto angle = [=](long double r) {
            return std::copysign(std::acos(std::abs(l) / (std::sqrt(a2) * r)), l);
        };
        const bool inward = h + nr < 0;
        const long double affine = length(boundary) + (inward ? length(r0) : -length(r0));
        const long double phi = angle(boundary) + (inward ? angle(r0) : -angle(r0));
        const auto actual = reference::Orbit(0, 0.01L, r0, boundary, nr, nphi, 8192);
        ASSERT_EQ(actual.fate, reference::Fate::Escape);
        EXPECT_LE(std::abs(actual.affine - affine), 1.0e-8L);
        EXPECT_LE(std::abs(actual.angle - phi), 1.0e-10L);
        EXPECT_LE(std::abs(actual.position[1] - boundary * std::cos(phi)), 1.0e-8L);
        EXPECT_LE(std::abs(actual.position[2] - boundary * std::sin(phi)), 1.0e-8L);
    }
    const auto roots = reference::Roots(1, 0.001L);
    EXPECT_LE(std::abs(reference::Lapse(roots.capture, 1, 0.001L)), 1.0e-14L);
    EXPECT_LE(std::abs(reference::Lapse(roots.cosmological, 1, 0.001L)), 1.0e-14L);
    // Radial outgoing-chart capture: r'=-C and t_out'=-C, including the
    // declared initial time gauge. Check its quadrature independently from a
    // derivative-free exact cubic partial-fraction integral of H/f.
    constexpr long double radius = 8, lambda = 0.001L;
    const long double h = 2 / radius + lambda * radius * radius / 3;
    const long double c = (1 - h) / std::sqrt(1 + h);
    const long double anchor = (roots.capture + roots.cosmological) / 2;
    long double primitive = -(radius - anchor);
    for (const auto root :
         {roots.capture, roots.cosmological, -roots.capture - roots.cosmological}) {
        const long double derivative = 2 / (root * root) - 2 * lambda * root / 3;
        primitive += std::log(std::abs((radius - root) / (anchor - root))) / derivative;
    }
    const auto radial = reference::Orbit(1, lambda, radius, 50, -1, 0, 8192);
    ASSERT_EQ(radial.fate, reference::Fate::Capture);
    EXPECT_LE(std::abs(radial.affine - (radius - roots.capture) / c), 1.0e-9L);
    EXPECT_LE(std::abs(radial.position[0] - (-2 * primitive - (radius - roots.capture))), 1.0e-8L);
    EXPECT_LE(std::abs(radial.tangent[0] + c), 1.0e-14L);
    EXPECT_LE(std::abs(radial.tangent[1] + c), 1.0e-14L);
    EXPECT_EQ(radial.angular_momentum, 0);
}

TEST(CpuTraceBoundary, NonradialKottlerFiniteEventsMatchIndependentQuadrature) {
    namespace reference = sirius::test::kottler_reference;
    struct Case {
        double mass, lambda, radius, radial, azimuth;
        bool capture;
    };
    const std::array<Case, 6> cases{{
        {0, 0.01, 5, -0.6, 0.8, false},
        {0, 0.01, 15, 0.6, -0.8, false},
        {1, 0.001, 8, -std::sqrt(0.96), 0.2, true},
        {1, 0.01, 8, -0.8, -0.6, false},
        {0.1, 0.1, 0.8, -0.8, 0.6, false},
        {100, 1.0e-7, 800, -std::sqrt(0.96), -0.2, true},
    }};
    constexpr double kEnvelope = 1.0e-4;
    long double maximum_reference_gap = 0, maximum_position_error = 0;
    long double maximum_tangent_error = 0, maximum_affine_error = 0, maximum_sky_error = 0;
    for (const auto& input : cases) {
        const auto id = input.mass == 0 ? MetricId::DeSitter : MetricId::SchwarzschildDeSitter;
        sirius::core::MetricConstructionParameters parameters;
        parameters.mass = input.mass;
        parameters.cosmological_constant = input.lambda;
        auto metric = sirius::core::CreateCpuMetric(id, parameters);
        ASSERT_NE(metric, nullptr);
        const auto domain = sirius::render::BuildTraceDomainParameters(
            {id, input.mass, input.lambda, input.radius, parameters.throat_radius,
             parameters.bubble_radius, parameters.bubble_sigma});
        ASSERT_TRUE(domain.finite_causal_boundary);
        const double scale = input.mass == 0 ? 1 : input.mass;
        CameraConfig camera_config;
        camera_config.r = input.radius;
        camera_config.theta = std::numbers::pi / 2;
        camera_config.phi = 0;
        camera_config.yaw = static_cast<float>(std::atan2(-input.azimuth, -input.radial));
        camera_config.width = camera_config.height = 3;
        PinholeCamera camera(camera_config);
        // One centre sample per family: six camera seeds, each traced at three
        // refinements below. This is not a nine-pixel image acceptance check.
        const auto ray = camera.GenerateRay(1, 1, 0.5f, 0.5f);
        ASSERT_EQ(ray.direction(2), 0);
        // The local ray slots are (radial,theta,azimuth), not Cartesian xyz.
        // theta=pi/2, phi=0 embeds the orbit in XY and +azimuth in +y.
        // Binary64 cos(pi/2) introduces only <1e-14 scaled plane error.
        const auto orbit = [&](unsigned panels) {
            return reference::Orbit(input.mass / scale, input.lambda * scale * scale,
                                    input.radius / scale, domain.escape_radius / scale,
                                    ray.direction(1), ray.direction(3), panels);
        };
        const auto coarse = orbit(4096), expected = orbit(8192);
        // Check both references before std::max: a NaN coarse field can
        // otherwise disappear while reducing the refinement uncertainty.
        for (const auto* state : {&coarse, &expected}) {
            ASSERT_TRUE(std::isfinite(state->affine));
            ASSERT_TRUE(std::isfinite(state->angle));
            for (int axis = 0; axis < 4; ++axis) {
                ASSERT_TRUE(std::isfinite(state->position[axis]));
                ASSERT_TRUE(std::isfinite(state->tangent[axis]));
                ASSERT_TRUE(std::isfinite(state->sky[axis]));
            }
        }
        ASSERT_EQ(expected.fate,
                  input.capture ? reference::Fate::Capture : reference::Fate::Escape);
        ASSERT_EQ(std::signbit(expected.angular_momentum), std::signbit(input.azimuth));
        ASSERT_EQ(expected.turning_radius > 0, !input.capture && input.radial < 0);
        long double reference_gap = std::max(std::abs(coarse.affine - expected.affine),
                                             std::abs(coarse.angle - expected.angle));
        for (int axis = 0; axis < 4; ++axis) {
            reference_gap =
                std::max(reference_gap, std::abs(coarse.position[axis] - expected.position[axis]));
            reference_gap =
                std::max(reference_gap, std::abs(coarse.tangent[axis] - expected.tangent[axis]));
            reference_gap =
                std::max(reference_gap, std::abs(coarse.sky[axis] - expected.sky[axis]));
        }
        ASSERT_LT(reference_gap, 1.0e-6L);
        maximum_reference_gap = std::max(maximum_reference_gap, reference_gap);
        for (int refinement = 0; refinement < 3; ++refinement) {
            SCOPED_TRACE(::testing::Message()
                         << "M=" << input.mass << " Lambda=" << input.lambda
                         << " r=" << input.radius << " L=" << expected.angular_momentum
                         << " capture=" << input.capture << " refinement=" << refinement);
            TracerConfig config;
            config.escape_radius = domain.escape_radius;
            config.finite_causal_boundary = domain.finite_causal_boundary;
            config.horizon_factor = sirius::render::kRenderTraceCaptureFactor;
            config.max_steps = sirius::render::kRenderTraceMaximumAttempts;
            config.enable_disk = false;
            config.integrator.initial_step = domain.cpu_initial_step;
            config.integrator.min_step = domain.cpu_min_step;
            config.integrator.max_step = std::ldexp(domain.max_step, -refinement);
            config.integrator.abs_tolerance = std::ldexp(5.0e-6f, -refinement);
            config.integrator.rel_tolerance = config.integrator.abs_tolerance;
            GeodesicTracer tracer(metric.get(), config);
            const auto actual = tracer.Trace(ray);
            ASSERT_EQ(actual.outcome, input.capture ? TraceResult::Outcome::Horizon
                                                    : TraceResult::Outcome::Escaped);
            ASSERT_FALSE(actual.numerical_failure);
            ASSERT_FALSE(actual.cancelled);
            ASSERT_EQ(actual.integrator_termination, 0);
            ASSERT_EQ(actual.terminal_chart, input.capture
                                                 ? TraceResult::TerminalChart::OutgoingKerrSchild
                                                 : TraceResult::TerminalChart::MetricNative);
            ASSERT_EQ(actual.asymptotic_sheet, TraceResult::AsymptoticSheet::Observer);
            ASSERT_TRUE(actual.final_tangent.has_value());
            const long double affine_error =
                std::abs(actual.affine_length / scale - expected.affine);
            EXPECT_LT(affine_error, kEnvelope);
            maximum_affine_error = std::max(maximum_affine_error, affine_error);
            for (int axis = 0; axis < 4; ++axis) {
                const long double position_error =
                    std::abs(actual.final_position(axis) / scale - expected.position[axis]);
                const long double tangent_error =
                    std::abs((*actual.final_tangent)(axis)-expected.tangent[axis]);
                EXPECT_LT(position_error, kEnvelope);
                EXPECT_LT(tangent_error, kEnvelope);
                maximum_position_error = std::max(maximum_position_error, position_error);
                maximum_tangent_error = std::max(maximum_tangent_error, tangent_error);
                if (!input.capture) {
                    const long double sky_error =
                        std::abs(actual.final_direction(axis) - expected.sky[axis]);
                    EXPECT_LT(sky_error, kEnvelope);
                    maximum_sky_error = std::max(maximum_sky_error, sky_error);
                }
            }
            if (!input.capture) {
                EXPECT_LE(TerminalRadius(actual),
                          static_cast<double>(domain.escape_radius) + 1.0e-7 * scale);
            }
        }
    }
    const auto record = [](const char* name, long double value) {
        std::ostringstream text;
        text << std::setprecision(17) << value;
        ::testing::Test::RecordProperty(name, text.str());
    };
    RecordProperty("physical_case_count", 6);
    RecordProperty("actual_trace_count", 18);
    RecordProperty("reference_coarse_panels", 4096);
    RecordProperty("reference_fine_panels", 8192);
    record("reference_max_normalized_refinement_gap", maximum_reference_gap);
    record("max_normalized_position_error", maximum_position_error);
    record("max_tangent_error", maximum_tangent_error);
    record("max_normalized_affine_error", maximum_affine_error);
    record("max_finite_sky_error", maximum_sky_error);
    RecordProperty("reference_scope",
                   "finite_CPU_XY_spherical_positive_Lambda_no_beam_or_transfer");
}

// Preserve the exact analytic metric while exercising the existing native
// IMetric trace route as an independent chart reference for outward rays.
class NativeKottlerReference final : public sirius::core::IMetric {
  public:
    explicit NativeKottlerReference(KerrSchildFamily& source) : source_(source) {}
    void Evaluate(const sirius::core::Vec4& position, sirius::core::Metric4d& metric,
                  sirius::core::Tensor<sirius::core::Dual<double>, 4, 4, 4>& derivative) override {
        source_.Evaluate(position, metric, derivative);
    }
    bool InverseMetric(const sirius::core::Vec4& position,
                       sirius::core::Metric4d& inverse) const override {
        return source_.InverseMetric(position, inverse);
    }
    bool IsValidEvent(const sirius::core::Vec4& position) const override {
        return source_.IsValidEvent(position);
    }
    bool InsideCaptureSurface(const sirius::core::Vec4& position, double margin) const override {
        return source_.InsideCaptureSurface(position, margin);
    }
    const sirius::core::Config& GetParameters() const override { return source_.GetParameters(); }
    void SetParameter(const std::string& key, double value) override {
        source_.SetParameter(key, value);
    }
    const char* GetName() const override { return source_.GetName(); }

  private:
    KerrSchildFamily& source_;
};

TracerConfig KottlerBoundaryConfig(const KerrSchildFamily& metric, float maximum_step) {
    TracerConfig config;
    const double horizon = metric.CosmologicalHorizonRadius();
    config.escape_radius = static_cast<float>(horizon);
    if (config.escape_radius > horizon) {
        config.escape_radius = std::nextafter(config.escape_radius, 0.0f);
    }
    config.finite_causal_boundary = true;
    config.horizon_factor = 1.0f;
    config.enable_disk = false;
    config.max_steps = 20000;
    config.integrator.initial_step = 0.1f;
    config.integrator.max_step = maximum_step;
    config.integrator.min_step = 1.0e-5f;
    config.integrator.abs_tolerance = 5.0e-6f;
    config.integrator.rel_tolerance = 5.0e-6f;
    return config;
}

TEST(CpuTraceBoundary, KottlerBothPastHorizonsRetainRadialAffineOracleAndTracerReuse) {
    KerrSchildFamily metric({1.0, 0.0, 0.0, 0.001});
    const auto config = KottlerBoundaryConfig(metric, 2.0f);
    GeodesicTracer tracer(&metric, config);
    for (double radius : {8.0, 50.0, 8.0}) {
        for (double direction : {1.0, -1.0}) {
            SCOPED_TRACE(radius);
            SCOPED_TRACE(direction);
            CameraRay ray;
            ray.origin(1) = radius;
            ray.origin(2) = std::numbers::pi / 2.0;
            ray.direction(1) = direction;
            const auto result = tracer.Trace(ray);
            ASSERT_FALSE(result.numerical_failure);
            ASSERT_EQ(result.outcome, direction > 0.0 ? TraceResult::Outcome::Escaped
                                                      : TraceResult::Outcome::Horizon);
            const double boundary =
                direction > 0.0 ? config.escape_radius : metric.OuterHorizonRadius();
            // Radial Kottler null rays have constant dr/dlambda. The public
            // ingoing Eulerian observer fixes its value at launch to
            // (H+n)/sqrt(1+H), where H=2M/r+Lambda*r^2/3 and n=+/-1.
            const double h = 2.0 / radius + 0.001 * radius * radius / 3.0;
            const double radial_rate = (h + direction) / std::sqrt(1.0 + h);
            EXPECT_NEAR(result.affine_length, (boundary - radius) / radial_rate, 2.0e-3);
            EXPECT_NEAR(TerminalRadius(result), boundary, 2.0e-7);
            EXPECT_EQ(result.integrator_termination, 0);
            for (int component = 0; component < 4; ++component) {
                EXPECT_TRUE(std::isfinite(result.final_position(component)));
            }
            if (direction > 0.0) {
                EXPECT_NEAR(result.final_direction(1), 1.0, 2.0e-7);
            }
        }
    }
}

TEST(CpuTraceBoundary, KottlerOverlapPreservesCoupledSourceMapAgainstNativeChart) {
    KerrSchildFamily metric({1.0, 0.0, 0.0, 0.001});
    NativeKottlerReference native(metric);
    CameraRay ray;
    ray.origin(1) = 8.0;
    ray.origin(2) = std::numbers::pi / 2.0;
    ray.direction(1) = 0.8;
    ray.direction(2) = 0.6;
    for (float maximum_step : {1.0f, 0.5f}) {
        SCOPED_TRACE(maximum_step);
        auto config = KottlerBoundaryConfig(metric, maximum_step);
        config.enable_ray_bundles = true;
        config.bundle_point_source = true;
        config.bundle_angular_size = 1.0e-3f;
        GeodesicTracer switched(&metric, config);
        GeodesicTracer reference(&native, config);
        const auto actual = switched.Trace(ray);
        const auto expected = reference.Trace(ray);
        ASSERT_FALSE(actual.numerical_failure);
        ASSERT_FALSE(expected.numerical_failure);
        ASSERT_EQ(actual.outcome, TraceResult::Outcome::Escaped);
        ASSERT_EQ(expected.outcome, TraceResult::Outcome::Escaped);
        ASSERT_TRUE(actual.beam.finite_source_map);
        ASSERT_TRUE(expected.beam.finite_source_map);
        EXPECT_NEAR(actual.affine_length, expected.affine_length, 2.0e-4);
        for (int component = 0; component < 4; ++component) {
            EXPECT_NEAR(actual.final_position(component), expected.final_position(component),
                        2.0e-4);
            EXPECT_NEAR(actual.final_direction(component), expected.final_direction(component),
                        2.0e-5);
        }
        EXPECT_NEAR(actual.beam.semi_major, expected.beam.semi_major, 2.0e-5);
        EXPECT_NEAR(actual.beam.semi_minor, expected.beam.semi_minor, 2.0e-5);
        for (int row = 0; row < 2; ++row) {
            for (int column = 0; column < 2; ++column) {
                EXPECT_NEAR(actual.beam.finite_source_map->jacobian[row][column],
                            expected.beam.finite_source_map->jacobian[row][column], 5.0e-4);
            }
        }
    }
}

}  // namespace
