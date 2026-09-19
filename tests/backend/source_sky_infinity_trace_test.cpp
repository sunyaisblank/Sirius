// Actual CPU trace handoffs to radial infinity. Analytic flat Lorentz geometry
// and independently retraced, fixed-pupil neighbours are the angular oracles.
#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/metrics/kerr_schild_family.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <iostream>
#include <numbers>
#include <optional>

namespace {
using sirius::backend::GeodesicTracer;
using sirius::backend::TracerConfig;
using sirius::backend::TraceResult;
using sirius::core::CameraRay;
using sirius::core::KerrSchildFamily;
using sirius::core::KerrSchildParams;
using Vector = std::array<double, 3>;
using Matrix = std::array<std::array<double, 2>, 2>;

double Dot(const Vector& a, const Vector& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vector Unit(Vector value) {
    const double length = std::hypot(value[0], value[1], value[2]);
    for (double& component : value) component /= length;
    return value;
}
Vector Cross(const Vector& a, const Vector& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
Vector Direction(const sirius::core::Vec4& value) { return Unit({value(1), value(2), value(3)}); }
std::array<Vector, 2> TangentBasis(const Vector& direction) {
    // Reproduce the declared coordinate convention, independently of the
    // production observer/frame and angular-map implementations.
    std::size_t least = 0;
    for (std::size_t axis = 1; axis < 3; ++axis)
        if (std::abs(direction[axis]) < std::abs(direction[least])) least = axis;
    Vector first{};
    for (std::size_t axis = 0; axis < 3; ++axis)
        first[axis] = (axis == least ? 1.0 : 0.0) - direction[least] * direction[axis];
    first = Unit(first);
    return {first, Cross(direction, first)};
}
double Difference(const Matrix& a, const Matrix& b) {
    double maximum = 0.0;
    for (std::size_t row = 0; row < 2; ++row)
        for (std::size_t column = 0; column < 2; ++column)
            maximum = std::max(maximum, std::abs(a[row][column] - b[row][column]));
    return maximum;
}

CameraRay PupilRay(double radius, double theta, double phi, Vector direction) {
    CameraRay ray;
    ray.origin(1) = radius;
    ray.origin(2) = theta;
    ray.origin(3) = phi;
    direction = Unit(direction);
    for (int axis = 0; axis < 3; ++axis) ray.direction(axis + 1) = direction[axis];
    ray.aperture_up = 0.021;
    ray.aperture_right = -0.017;
    return ray;
}
TracerConfig MapConfig(float radius, float maximum_step, float tolerance) {
    TracerConfig config;
    config.escape_radius = radius;
    config.max_steps = 20000;
    config.enable_disk = false;
    config.enable_ray_bundles = true;
    config.bundle_point_source = true;
    config.bundle_angular_size = 1.0e-3f;
    config.integrator.initial_step = maximum_step;
    config.integrator.max_step = maximum_step;
    config.integrator.min_step = 1.0e-6f;
    config.integrator.abs_tolerance = tolerance;
    config.integrator.rel_tolerance = tolerance;
    return config;
}
TraceResult Trace(KerrSchildFamily& metric, const CameraRay& ray, const TracerConfig& config) {
    GeodesicTracer tracer(&metric, config);
    return tracer.Trace(ray);
}

void ExpectSameInfinity(const TraceResult& actual, const TraceResult& reference) {
    ASSERT_FALSE(actual.numerical_failure);
    ASSERT_FALSE(reference.numerical_failure);
    ASSERT_EQ(actual.outcome, TraceResult::Outcome::Escaped);
    ASSERT_EQ(reference.outcome, TraceResult::Outcome::Escaped);
    ASSERT_TRUE(actual.beam.infinity_source_map);
    ASSERT_TRUE(reference.beam.infinity_source_map);
    const auto& a = *actual.beam.infinity_source_map;
    const auto& b = *reference.beam.infinity_source_map;
    EXPECT_LE(Difference(a.map.jacobian, b.map.jacobian), 2.0e-6);
    EXPECT_NEAR(a.frequency, b.frequency, 2.0e-8);
    for (std::size_t axis = 0; axis < 3; ++axis)
        EXPECT_NEAR(a.map.direction[axis], b.map.direction[axis], 2.0e-7);
    for (std::size_t column = 0; column < 2; ++column)
        EXPECT_NEAR(a.frequency_derivative[column], b.frequency_derivative[column], 2.0e-6);
    EXPECT_LE(a.maximum_local_error_ratio, 1.0);
}

struct FlatMap {
    Vector direction;
    Matrix jacobian;
    double frequency;
    std::array<double, 2> frequency_derivative{};
};
FlatMap AnalyticFlatMap(const CameraRay& ray) {
    const Vector n = Direction(ray.direction);
    const auto launch_basis = TangentBasis(n);
    const Vector beta{-ray.beta_forward, -ray.beta_up, ray.beta_right};
    const double beta_squared = Dot(beta, beta);
    const double gamma = 1.0 / std::sqrt(1.0 - beta_squared);
    const double boost_coefficient = beta_squared == 0.0 ? 0.0 : (gamma - 1.0) / beta_squared;
    const double frequency = gamma * (1.0 - Dot(beta, n));
    Vector local_direction{};
    for (std::size_t axis = 0; axis < 3; ++axis)
        local_direction[axis] =
            (n[axis] + (boost_coefficient * Dot(beta, n) - gamma) * beta[axis]) / frequency;

    const double theta = ray.origin(2), phi = ray.origin(3);
    const std::array<Vector, 3> triad{
        Vector{std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)},
        Vector{std::cos(theta) * std::cos(phi), std::cos(theta) * std::sin(phi), -std::sin(theta)},
        Vector{-std::sin(phi), std::cos(phi), 0.0}};
    const auto world = [&](const Vector& local) {
        Vector result{};
        for (std::size_t axis = 0; axis < 3; ++axis)
            for (std::size_t component = 0; component < 3; ++component)
                result[component] += triad[axis][component] * local[axis];
        return result;
    };
    FlatMap result{Unit(world(local_direction)), {}, frequency};
    const auto source_basis = TangentBasis(result.direction);
    for (std::size_t column = 0; column < 2; ++column) {
        const auto& e = launch_basis[column];
        const double frequency_derivative = -gamma * Dot(beta, e);
        result.frequency_derivative[column] = frequency_derivative;
        Vector derivative{};
        for (std::size_t axis = 0; axis < 3; ++axis)
            derivative[axis] = (e[axis] + boost_coefficient * Dot(beta, e) * beta[axis] -
                                local_direction[axis] * frequency_derivative) /
                               frequency;
        derivative = world(derivative);
        for (std::size_t row = 0; row < 2; ++row)
            result.jacobian[row][column] = Dot(source_basis[row], derivative);
    }
    return result;
}

std::optional<Matrix> RetracedJacobian(KerrSchildFamily& metric, const CameraRay& ray,
                                       TracerConfig config, const Vector& central_source,
                                       double angular_step) {
    const Vector direction = Direction(ray.direction);
    const auto launch_basis = TangentBasis(direction);
    const auto source_basis = TangentBasis(central_source);
    Matrix result{};
    for (std::size_t column = 0; column < 2; ++column) {
        std::array<Vector, 2> endpoints;
        for (int side = 0; side < 2; ++side) {
            // The full camera state and actual pupil are fixed; perturb only
            // its rest-frame launch angle along a great circle.
            CameraRay neighbour = ray;
            for (std::size_t axis = 0; axis < 3; ++axis)
                neighbour.direction(static_cast<int>(axis) + 1) =
                    std::cos(angular_step) * direction[axis] +
                    (side == 0 ? -1.0 : 1.0) * std::sin(angular_step) * launch_basis[column][axis];
            const auto trace = Trace(metric, neighbour, config);
            if (trace.numerical_failure || trace.outcome != TraceResult::Outcome::Escaped)
                return std::nullopt;
            if (!trace.beam.infinity_source_map) return std::nullopt;
            endpoints[side] = trace.beam.infinity_source_map->map.direction;
        }
        Vector derivative{};
        for (std::size_t axis = 0; axis < 3; ++axis)
            derivative[axis] = (endpoints[1][axis] - endpoints[0][axis]) / (2.0 * angular_step);
        for (std::size_t row = 0; row < 2; ++row)
            result[row][column] = Dot(source_basis[row], derivative);
    }
    return result;
}

TEST(SourceSkyInfinityTrace, FlatPupilMapAndFrequencyAreHandoffIndependent) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    for (bool moving : {false, true}) {
        auto ray = PupilRay(7.0, 1.3, -0.41, {0.6, 0.3, -0.74});
        if (moving) {
            ray.beta_forward = 0.2;
            ray.beta_up = -0.1;
            ray.beta_right = 0.05;
        }
        const auto expected = AnalyticFlatMap(ray);
        for (float handoff : {20.0f, 60.0f, 200.0f}) {
            SCOPED_TRACE(moving);
            SCOPED_TRACE(handoff);
            const auto result = Trace(metric, ray, MapConfig(handoff, 3.0f, 1.0e-9f));
            ASSERT_FALSE(result.numerical_failure);
            ASSERT_EQ(result.outcome, TraceResult::Outcome::Escaped);
            ASSERT_TRUE(result.beam.infinity_source_map);
            const auto& infinity = *result.beam.infinity_source_map;
            EXPECT_LE(Difference(infinity.map.jacobian, expected.jacobian), 2.0e-8);
            EXPECT_NEAR(infinity.frequency, expected.frequency, 2.0e-9);
            for (std::size_t column = 0; column < 2; ++column) {
                EXPECT_NEAR(infinity.frequency_derivative[column],
                            expected.frequency_derivative[column], 2.0e-8);
            }
            EXPECT_NEAR(infinity.map.determinant, 1.0 / (expected.frequency * expected.frequency),
                        2.0e-8);
            for (std::size_t axis = 0; axis < 3; ++axis) {
                EXPECT_NEAR(infinity.map.direction[axis], expected.direction[axis], 2.0e-9);
            }
            EXPECT_LE(infinity.maximum_local_error_ratio, 1.0);
        }
    }
}

TEST(SourceSkyInfinityTrace, KerrMapConvergesAcrossHandoffRadiiAndFixedPupilNeighbours) {
    for (double spin : {-0.7, 0.7}) {
        KerrSchildFamily metric(KerrSchildParams::Kerr(1.0, spin));
        auto ray = PupilRay(8.0, 1.1, 0.31, {0.5, 0.65, 0.57});
        ray.beta_forward = 0.08;
        ray.beta_up = 0.025;
        ray.beta_right = -0.015;
        std::optional<sirius::core::relativity::KerrInfinityResult> first;
        for (float handoff : {20.0f, 40.0f, 80.0f}) {
            SCOPED_TRACE(spin);
            SCOPED_TRACE(handoff);
            const auto config = MapConfig(handoff, 0.125f, 1.0e-10f);
            const auto result = Trace(metric, ray, config);
            ASSERT_FALSE(result.numerical_failure);
            ASSERT_EQ(result.outcome, TraceResult::Outcome::Escaped);
            ASSERT_TRUE(result.beam.infinity_source_map);
            const auto& infinity = *result.beam.infinity_source_map;
            const auto wide = RetracedJacobian(metric, ray, config, infinity.map.direction, 2.0e-4);
            const auto narrow =
                RetracedJacobian(metric, ray, config, infinity.map.direction, 1.0e-4);
            ASSERT_TRUE(wide);
            ASSERT_TRUE(narrow);
            EXPECT_LE(Difference(infinity.map.jacobian, *wide), 2.0e-6);
            EXPECT_LE(Difference(infinity.map.jacobian, *narrow), 2.0e-6);
            EXPECT_LE(Difference(*wide, *narrow), 2.0e-7);
            if (first) {
                EXPECT_LE(Difference(infinity.map.jacobian, first->map.jacobian), 2.0e-6);
                EXPECT_NEAR(infinity.frequency, first->frequency, 2.0e-8);
                for (std::size_t axis = 0; axis < 3; ++axis) {
                    EXPECT_NEAR(infinity.map.direction[axis], first->map.direction[axis], 2.0e-7);
                }
            } else {
                first = infinity;
            }
        }
    }
}

TEST(SourceSkyInfinityTrace, PointSourceHandoffPreservesTheFlatObserverMap) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    auto ray = PupilRay(7.0, 1.3, -0.41, {0.6, 0.3, -0.74});
    ray.beta_forward = 0.2;
    ray.beta_up = -0.1;
    ray.beta_right = 0.05;
    const auto expected = AnalyticFlatMap(ray);
    const auto config = MapConfig(200.0f, 3.0f, 1.0e-9f);
    GeodesicTracer tracer(&metric, config);
    const auto result = tracer.TracePointSource(ray);
    ASSERT_FALSE(result.numerical_failure);
    ASSERT_EQ(result.outcome, TraceResult::Outcome::Escaped);
    ASSERT_TRUE(result.beam.infinity_source_map);
    const auto& infinity = *result.beam.infinity_source_map;
    EXPECT_LE(Difference(infinity.map.jacobian, expected.jacobian), 2.0e-8);
    EXPECT_NEAR(infinity.frequency, expected.frequency, 2.0e-9);
    for (std::size_t axis = 0; axis < 3; ++axis)
        EXPECT_NEAR(infinity.map.direction[axis], expected.direction[axis], 2.0e-9);
    for (std::size_t column = 0; column < 2; ++column)
        EXPECT_NEAR(infinity.frequency_derivative[column], expected.frequency_derivative[column],
                    2.0e-8);
    EXPECT_LT(result.steps_taken, 10);
    EXPECT_LT(
        std::hypot(result.final_position(1), result.final_position(2), result.final_position(3)),
        20.0);
}

TEST(SourceSkyInfinityTrace, PointSourceHandoffAvoidsVacuumTravelForMovingPhysicalPupils) {
    sirius::core::CameraConfig camera_config;
    camera_config.r = 50;
    camera_config.theta = 60 * std::numbers::pi / 180;
    camera_config.width = 4;
    camera_config.height = 2;
    camera_config.fov = 2;
    camera_config.beta_x = .1;
    camera_config.beta_y = .8;
    camera_config.focus_distance = 50;
    sirius::core::ThinLensCamera camera(camera_config);
    int ordinary_attempts = 0, probe_attempts = 0;
    double direction_gap = 0, jacobian_gap = 0;
    for (double spin : {-0.7, 0.7}) {
        KerrSchildFamily metric(KerrSchildParams::Kerr(1.0, spin));
        auto config = MapConfig(200.0f, 2.0f, 5.0e-6f);
        config.max_steps = 20000;
        config.integrator.initial_step = .1f;
        GeodesicTracer tracer(&metric, config);
        for (const auto offset : {std::array{0.0, 0.0}, std::array{-1.0, .5}}) {
            SCOPED_TRACE(spin);
            SCOPED_TRACE(offset[0]);
            const auto projection =
                camera.ProjectFilmOffsetForObserver(1.5, .5, offset[0], offset[1], .2, 1.0 / 7.0);
            ASSERT_TRUE(projection);
            ASSERT_TRUE(projection->ray.phase_space);
            const auto reference = tracer.Trace(projection->ray);
            const auto result = tracer.TracePointSource(projection->ray);
            ExpectSameInfinity(result, reference);
            ordinary_attempts += reference.steps_taken;
            probe_attempts += result.steps_taken;
            if (result.beam.infinity_source_map && reference.beam.infinity_source_map) {
                const auto& actual_map = result.beam.infinity_source_map->map;
                const auto& reference_map = reference.beam.infinity_source_map->map;
                for (int axis = 0; axis < 3; ++axis)
                    direction_gap = std::max(
                        direction_gap,
                        std::abs(actual_map.direction[axis] - reference_map.direction[axis]));
                jacobian_gap =
                    std::max(jacobian_gap, Difference(actual_map.jacobian, reference_map.jacobian));
            }
            EXPECT_LT(result.steps_taken * 2, reference.steps_taken);
            EXPECT_LT(result.central_stages * 2, reference.central_stages);
            EXPECT_LT(std::hypot(result.final_position(1), result.final_position(2),
                                 result.final_position(3)),
                      100.0);
            EXPECT_EQ(result.optical_depth, reference.optical_depth);
            auto without_beam = config;
            without_beam.enable_ray_bundles = false;
            without_beam.bundle_point_source = false;
            without_beam.bundle_angular_size = TracerConfig{}.bundle_angular_size;
            GeodesicTracer unbundled(&metric, without_beam);
            const auto plain = unbundled.TracePointSource(projection->ray);
            ExpectSameInfinity(plain, result);
            EXPECT_EQ(plain.steps_taken, result.steps_taken);
            EXPECT_EQ(plain.central_stages, result.central_stages);
        }
    }
    RecordProperty("ordinary_inner_attempts", ordinary_attempts);
    RecordProperty("point_source_inner_attempts", probe_attempts);
    RecordProperty("maximum_direction_difference", std::format("{:.17g}", direction_gap));
    RecordProperty("maximum_jacobian_difference", std::format("{:.17g}", jacobian_gap));
}

TEST(SourceSkyInfinityTrace, PointSourceHandoffRetainsOpaqueAndVolumetricSources) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1.0));
    for (bool volume : {false, true}) {
        SCOPED_TRACE(volume);
        const auto ray = PupilRay(8.0, std::numbers::pi / 2 - .05, 0.0, {.8, .6, 0});
        auto config = MapConfig(60.0f, .25f, 1.0e-9f);
        config.enable_disk = true;
        config.enable_volumetric = volume;
        if (volume) config.volumetric_tau_midplane = .01f;
        GeodesicTracer tracer(&metric, config);
        const auto reference = tracer.Trace(ray);
        const auto result = tracer.TracePointSource(ray);
        ASSERT_FALSE(reference.numerical_failure);
        ASSERT_FALSE(result.numerical_failure);
        if (volume) {
            ExpectSameInfinity(result, reference);
            EXPECT_TRUE(result.volumetric_hit);
            EXPECT_GT(result.optical_depth, 0);
            EXPECT_EQ(result.optical_depth, reference.optical_depth);
            for (int channel = 0; channel < 3; ++channel)
                EXPECT_EQ(result.volumetric_emission[channel],
                          reference.volumetric_emission[channel]);
            EXPECT_LT(result.steps_taken, reference.steps_taken);
            EXPECT_GT(std::hypot(result.final_position(1), result.final_position(2),
                                 result.final_position(3)),
                      config.disk_outer);
        } else {
            EXPECT_EQ(result.outcome, TraceResult::Outcome::DiskHit);
            EXPECT_EQ(result.outcome, reference.outcome);
            EXPECT_EQ(result.steps_taken, reference.steps_taken);
            EXPECT_EQ(result.disk_radius, reference.disk_radius);
            EXPECT_FALSE(result.beam.infinity_source_map);
        }
    }
}

TEST(SourceSkyInfinityTrace, PointSourceWorkExhaustionAndCancellationHaveNoSourceMap) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1.0));
    auto config = MapConfig(60.0f, .125f, 1.0e-9f);
    config.max_steps = 1;
    const auto ray = PupilRay(8.0, 1.1, 0, {-1, 0, 0});
    GeodesicTracer tracer(&metric, config);
    const auto exhausted = tracer.TracePointSource(ray);
    EXPECT_TRUE(exhausted.numerical_failure);
    EXPECT_EQ(exhausted.outcome, TraceResult::Outcome::MaxSteps);
    EXPECT_EQ(exhausted.coupled_failure, sirius::core::CoupledStepFailure::WorkLimit);
    EXPECT_FALSE(exhausted.beam.infinity_source_map);
    tracer.SetCancellationCallback([] { return true; });
    const auto cancelled = tracer.TracePointSource(ray);
    EXPECT_TRUE(cancelled.cancelled);
    EXPECT_EQ(cancelled.steps_taken, 0);
    EXPECT_FALSE(cancelled.beam.infinity_source_map);
}

TEST(SourceSkyInfinityTrace, NonVacuumAndCapturedRaysDoNotClaimVacuumInfinity) {
    const auto outward = PupilRay(8.0, 1.1, 0.31, {0.8, 0.3, 0.1});
    for (const auto parameters :
         {KerrSchildParams::ReissnerNordstrom(1.0, 0.2), KerrSchildParams::DeSitter(0.001)}) {
        KerrSchildFamily metric(parameters);
        auto config = MapConfig(20.0f, 0.125f, 1.0e-10f);
        config.finite_causal_boundary = parameters.Lambda > 0.0;
        const auto result = Trace(metric, outward, config);
        ASSERT_FALSE(result.numerical_failure);
        ASSERT_EQ(result.outcome, TraceResult::Outcome::Escaped);
        ASSERT_TRUE(result.beam.finite_source_map);
        EXPECT_FALSE(result.beam.infinity_source_map);
        EXPECT_FALSE(result.beam.infinity_source_failure);
        GeodesicTracer tracer(&metric, config);
        const auto probe = tracer.TracePointSource(outward);
        ASSERT_FALSE(probe.numerical_failure);
        EXPECT_EQ(probe.outcome, result.outcome);
        EXPECT_EQ(probe.steps_taken, result.steps_taken);
        EXPECT_FALSE(probe.beam.infinity_source_map);
    }
    KerrSchildFamily schwarzschild(KerrSchildParams::Schwarzschild(1.0));
    CameraRay captured;
    captured.origin(1) = 8.0;
    captured.origin(2) = 1.1;
    captured.direction(1) = -1.0;
    const auto result = Trace(schwarzschild, captured, MapConfig(20.0f, 0.125f, 1.0e-10f));
    ASSERT_FALSE(result.numerical_failure);
    ASSERT_EQ(result.outcome, TraceResult::Outcome::Horizon);
    EXPECT_FALSE(result.beam.infinity_source_map);
    EXPECT_FALSE(result.beam.infinity_source_failure);
    GeodesicTracer tracer(&schwarzschild, MapConfig(20.0f, 0.125f, 1.0e-10f));
    const auto probe = tracer.TracePointSource(captured);
    ASSERT_FALSE(probe.numerical_failure);
    EXPECT_EQ(probe.outcome, TraceResult::Outcome::Horizon);
    EXPECT_EQ(probe.steps_taken, result.steps_taken);
    EXPECT_FALSE(probe.beam.infinity_source_map);
}

}  // namespace
