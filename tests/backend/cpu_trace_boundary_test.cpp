// CPU accepted-segment event witnesses. These are backend computations only:
// no render session, image dispatch, or presentation path is exercised.

#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera.h"
#include "sirius/core/metrics/kerr_schild_family.h"
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

double TerminalRadius(const TraceResult& trace) {
    return std::sqrt(trace.final_position(1) * trace.final_position(1) +
                     trace.final_position(2) * trace.final_position(2) +
                     trace.final_position(3) * trace.final_position(3));
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

}  // namespace
