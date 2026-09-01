// CPU accepted-segment event witnesses. These are backend computations only:
// no render session, image dispatch, or presentation path is exercised.

#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera.h"
#include "sirius/core/metrics/cpu_metric_factory.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/metrics/morris_thorne_family.h"
#include "sirius/core/metrics/registry.h"

#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

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

}  // namespace
