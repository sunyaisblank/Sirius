// Actual CPU trace handoffs to radial infinity. Analytic flat Lorentz geometry
// and independently retraced, fixed-pupil neighbours are the angular oracles.
#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/metrics/kerr_schild_family.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
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
}

}  // namespace
