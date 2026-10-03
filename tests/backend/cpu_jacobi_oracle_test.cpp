// Production CPU Jacobi witnesses against independent Schwarzschild geometry.
// These are backend calculations only: no render session, image, dispatch, or
// presentation path is created.

#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/oracle/kerr_boyer_lindquist.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>

namespace {

using sirius::backend::GeodesicTracer;
using sirius::backend::TracerConfig;
using sirius::backend::TraceResult;
using sirius::core::CameraRay;
using sirius::core::KerrSchildFamily;
using sirius::core::KerrSchildParams;
using sirius::core::Vec4;
using sirius::oracle::KerrMetricD;
using sirius::oracle::Vec4d;

Vec4 ToLiveVector(const Vec4d& vector_bl, const Vec4d& event_bl, double mass) {
    const double sin_theta = std::sin(event_bl.theta);
    const double cos_theta = std::cos(event_bl.theta);
    const double sin_phi = std::sin(event_bl.phi);
    const double cos_phi = std::cos(event_bl.phi);
    Vec4 result;
    result(0) = vector_bl.t + (2.0 * mass / (event_bl.r - 2.0 * mass)) * vector_bl.r;
    result(1) = sin_theta * cos_phi * vector_bl.r +
                event_bl.r * cos_theta * cos_phi * vector_bl.theta -
                event_bl.r * sin_theta * sin_phi * vector_bl.phi;
    result(2) = sin_theta * sin_phi * vector_bl.r +
                event_bl.r * cos_theta * sin_phi * vector_bl.theta +
                event_bl.r * sin_theta * cos_phi * vector_bl.phi;
    result(3) = cos_theta * vector_bl.r - event_bl.r * sin_theta * vector_bl.theta;
    return result;
}

Vec4d ToOracleVector(const Vec4& vector_cart, const Vec4d& event_bl, double mass) {
    const double sin_theta = std::sin(event_bl.theta);
    const double cos_theta = std::cos(event_bl.theta);
    const double sin_phi = std::sin(event_bl.phi);
    const double cos_phi = std::cos(event_bl.phi);
    const double radial = sin_theta * cos_phi * vector_cart(1) +
                          sin_theta * sin_phi * vector_cart(2) + cos_theta * vector_cart(3);
    const double polar = (cos_theta * cos_phi * vector_cart(1) +
                          cos_theta * sin_phi * vector_cart(2) - sin_theta * vector_cart(3)) /
                         event_bl.r;
    const double azimuthal =
        (-sin_phi * vector_cart(1) + cos_phi * vector_cart(2)) / (event_bl.r * sin_theta);
    return Vec4d(vector_cart(0) - (2.0 * mass / (event_bl.r - 2.0 * mass)) * radial, radial, polar,
                 azimuthal);
}

Vec4d OracleTidalAcceleration(KerrMetricD& metric, const Vec4d& event, const Vec4d& tangent,
                              const Vec4d& deviation) {
    double riemann[4][4][4][4];
    metric.Riemann(event, riemann);
    Vec4d acceleration;
    for (int mu = 0; mu < 4; ++mu) {
        double contraction = 0.0;
        for (int nu = 0; nu < 4; ++nu)
            for (int rho = 0; rho < 4; ++rho)
                for (int sigma = 0; sigma < 4; ++sigma)
                    contraction +=
                        riemann[mu][nu][rho][sigma] * tangent[nu] * deviation[rho] * tangent[sigma];
        acceleration[mu] = -contraction;
    }
    return acceleration;
}

Vec4 LiveEvent(double radius, double theta, double phi, double spin) {
    Vec4 position;
    position(1) = (radius * std::cos(phi) - spin * std::sin(phi)) * std::sin(theta);
    position(2) = (radius * std::sin(phi) + spin * std::cos(phi)) * std::sin(theta);
    position(3) = radius * std::cos(theta);
    return position;
}

TEST(CpuJacobiOracle, TidalContractionMatchesAnalyticSchwarzschildAtMatchedEvents) {
    constexpr double mass = 1.0;
    KerrSchildFamily live_metric(KerrSchildParams::Schwarzschild(mass));
    GeodesicTracer tracer(&live_metric, TracerConfig{});
    KerrMetricD oracle_metric(mass, 0.0);

    for (const double radius : {4.0, 6.0, 10.0, 20.0}) {
        constexpr double theta = 1.1;
        const Vec4d event(0.0, radius, theta, 0.37);
        const double lapse = 1.0 - 2.0 * mass / radius;
        constexpr double energy = 1.0;
        constexpr double polar_momentum = 1.3;
        constexpr double angular_momentum = 2.0;
        const double sin_theta = std::sin(theta);
        const double angular_norm =
            (polar_momentum * polar_momentum +
             angular_momentum * angular_momentum / (sin_theta * sin_theta)) /
            (radius * radius);
        const double radial = std::sqrt(energy * energy - lapse * angular_norm);
        const Vec4d tangent(energy / lapse, radial, polar_momentum / (radius * radius),
                            angular_momentum / (radius * radius * sin_theta * sin_theta));
        const Vec4d deviation(0.2, 0.15, 0.07 / radius, -0.11 / (radius * sin_theta));

        // Events are points, not vectors: use the exact Schwarzschild spatial
        // chart map and exploit stationarity for the arbitrary time origin.
        const Vec4 position = LiveEvent(event.r, event.theta, event.phi, 0.0);

        const Vec4 live_acceleration = tracer.TidalAcceleration(
            position, ToLiveVector(tangent, event, mass), ToLiveVector(deviation, event, mass));
        const Vec4d actual = ToOracleVector(live_acceleration, event, mass);
        const Vec4d expected = OracleTidalAcceleration(oracle_metric, event, tangent, deviation);

        double scale = 0.0;
        double maximum_error = 0.0;
        for (int component = 0; component < 4; ++component) {
            scale = std::max(scale, std::abs(expected[component]));
            maximum_error =
                std::max(maximum_error, std::abs(actual[component] - expected[component]));
        }
        EXPECT_LT(maximum_error / std::max(scale, 1.0e-30), 5.0e-7)
            << "matched-event tidal contraction at r=" << radius;
    }
}

TEST(CpuJacobiOracle, CurvatureScalarMatchesAnalyticKerrOffEquator) {
    constexpr double mass = 1.0;
    for (const double spin : {0.0, 0.9}) {
        KerrSchildFamily live_metric(KerrSchildParams::Kerr(mass, spin));
        GeodesicTracer tracer(&live_metric, TracerConfig{});
        KerrMetricD oracle_metric(mass, spin);
        for (const double radius : {3.0, 4.0, 8.0, 20.0}) {
            for (const double theta : {0.6, 1.1, std::numbers::pi / 2.0}) {
                const double actual =
                    tracer.KretschmannScalar(LiveEvent(radius, theta, 0.37, spin));
                const double expected = oracle_metric.Kretschmann(Vec4d(0.0, radius, theta, 0.37));
                const double relative_error =
                    std::abs(actual - expected) / std::max(std::abs(expected), 1.0e-30);
                EXPECT_LT(relative_error, 1.0e-6) << "Kerr curvature scalar at a=" << spin
                                                  << " r=" << radius << " theta=" << theta;
            }
        }
    }
}

TEST(CpuJacobiOracle, RadialPointSourceCongruenceMatchesClosedForm) {
    constexpr double mass = 1.0;
    constexpr float boundary_radius = 12.0f;
    constexpr float angular_derivative = 1.0e-3f;
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(mass));

    TracerConfig config;
    config.escape_radius = boundary_radius;
    config.finite_causal_boundary = true;
    config.max_steps = 1000;
    config.enable_disk = false;
    config.enable_ray_bundles = true;
    config.bundle_point_source = true;
    config.bundle_angular_size = angular_derivative;
    config.integrator.initial_step = 2.0e-2f;
    config.integrator.max_step = 2.0e-2f;
    config.integrator.min_step = 1.0e-6f;
    config.integrator.abs_tolerance = 1.0e-8f;
    config.integrator.rel_tolerance = 1.0e-8f;

    CameraRay ray;
    ray.origin(1) = 10.0;
    ray.origin(2) = std::numbers::pi / 2.0;
    ray.direction(1) = 1.0;  // Past-directed ray leaves the hole radially.

    GeodesicTracer tracer(&metric, config);
    const TraceResult result = tracer.Trace(ray);

    ASSERT_EQ(result.outcome, TraceResult::Outcome::Escaped);
    ASSERT_FALSE(result.numerical_failure);
    ASSERT_TRUE(result.beam.valid);
    ASSERT_GT(result.affine_length, 0.0f);

    // Morales-Ruiz and Raposo (2023), eqs. 19-24: for point-source screen data
    // xi(0)=0 with unit physical Dxi/dlambda, both radial Schwarzschild screen
    // axes equal lambda exactly. The live seed scales that derivative by eps.
    const double expected_axis =
        static_cast<double>(angular_derivative) * static_cast<double>(result.affine_length);
    EXPECT_NEAR(result.beam.semi_major, expected_axis, 2.0e-6 * expected_axis);
    EXPECT_NEAR(result.beam.semi_minor, expected_axis, 2.0e-6 * expected_axis);
    EXPECT_NEAR(result.beam.transverse_area, expected_axis * expected_axis,
                4.0e-6 * expected_axis * expected_axis);
}

TEST(CpuJacobiOracle, FlatScreenEllipseRetainsAnisotropyAndScale) {
    KerrSchildFamily metric(KerrSchildParams::Minkowski());
    TracerConfig config;
    config.escape_radius = 12.0f;
    config.enable_disk = false;
    config.enable_ray_bundles = true;
    config.integrator.initial_step = config.integrator.max_step = 2.0f;
    config.max_steps = 100;
    GeodesicTracer tracer(&metric, config);
    // A flat, outward central ray with a separately specified constant pupil
    // displacement family. Its transverse map is R(angle) diag(major,minor);
    // flat parallel transport preserves that map exactly at the source event.
    // No production geometry calculation supplies the expected singular axes.
    for (const auto shape : {std::array{1.0, 1.0e-9, 0.0}, std::array{3.0, 1.0, 0.4},
                             std::array{1.0, 0.0, 0.0}, std::array{0.0, 0.0, 0.0}}) {
        for (const double scale : {1.0e-10, 1.0, 1.0e10}) {
            SCOPED_TRACE(shape[0]);
            SCOPED_TRACE(shape[1]);
            SCOPED_TRACE(scale);
            CameraRay ray;
            constexpr double launch_radius = 5.0;
            ray.origin(1) = launch_radius;
            ray.origin(2) = std::numbers::pi / 2;
            ray.direction(1) = 1.0;
            ray.phase_space.emplace();
            auto& family = *ray.phase_space;
            family.direction[1][0] = family.direction[2][1] = 1;
            family.film_to_angle[0][0] = family.film_to_angle[1][1] = 1;
            const double c = std::cos(shape[2]), s = std::sin(shape[2]);
            family.pupil_right[2] = launch_radius * scale * shape[0] * c;
            family.pupil_up[2] = launch_radius * scale * shape[0] * s;
            family.pupil_right[3] = -launch_radius * scale * shape[1] * s;
            family.pupil_up[3] = launch_radius * scale * shape[1] * c;
            const auto result = tracer.Trace(ray);
            ASSERT_FALSE(result.numerical_failure);
            ASSERT_EQ(result.outcome, TraceResult::Outcome::Escaped);
            ASSERT_TRUE(result.beam.valid);
            const double seed = config.bundle_angular_size;
            const double expected_major = seed * scale * shape[0];
            const double expected_minor = seed * scale * shape[1];
            const double expected_area = expected_major * expected_minor;
            EXPECT_NEAR(result.beam.semi_major, expected_major, 5.0e-6 * expected_major);
            EXPECT_NEAR(result.beam.semi_minor, expected_minor, 5.0e-6 * expected_minor);
            EXPECT_NEAR(result.beam.transverse_area, expected_area, 5.0e-6 * expected_area);
            EXPECT_NEAR(static_cast<double>(result.beam.semi_major) * result.beam.semi_minor,
                        expected_area, 1.0e-5 * expected_area);
            if (shape[0] > shape[1]) {
                EXPECT_NEAR(std::cos(2 * result.beam.orientation), std::cos(2 * shape[2]), 1.0e-6);
                EXPECT_NEAR(std::sin(2 * result.beam.orientation), std::sin(2 * shape[2]), 1.0e-6);
            }
        }
    }

    // Dyadic seed/radius and a flat launch whose polar corrections round below
    // binary64 resolution ensure the pupil map reaches extraction unchanged.
    // The camera domain excludes the exact pole. In this regular chart the
    // screen is (+y,-x); (3,1)^T(1,3) has det=0, nextafter(d) has det=3*(d-3).
    config.bundle_angular_size = 1.f / 64;
    GeodesicTracer rank_tracer(&metric, config);
    for (const double d : {3.0, std::nextafter(3.0, 4.0)}) {
        SCOPED_TRACE(d);
        CameraRay ray;
        ray.origin(1) = 4;
        ray.origin(2) = 1e-100;
        ray.direction(1) = 1;
        ray.phase_space.emplace();
        auto& family = *ray.phase_space;
        family.direction[1][0] = family.direction[2][1] = 1;
        family.film_to_angle[0][0] = family.film_to_angle[1][1] = 1;
        constexpr double inverse_publication_scale = 256;
        family.pupil_right[2] = 3 * inverse_publication_scale;
        family.pupil_up[2] = inverse_publication_scale;
        family.pupil_right[3] = 9 * inverse_publication_scale;
        family.pupil_up[3] = d * inverse_publication_scale;
        const auto result = rank_tracer.Trace(ray);
        ASSERT_FALSE(result.numerical_failure);
        ASSERT_EQ(result.outcome, TraceResult::Outcome::Escaped);
        ASSERT_TRUE(result.beam.valid);
        const double expected_area = 3 * (d - 3);
        const double expected_minor = expected_area / 10;
        EXPECT_NEAR(result.beam.semi_major, 10, 5e-6 * 10);
        EXPECT_NEAR(result.beam.semi_minor, expected_minor, 5e-6 * expected_minor);
        EXPECT_NEAR(result.beam.transverse_area, expected_area, 5e-6 * expected_area);
        EXPECT_NEAR(result.beam.orientation, std::atan2(1.0, 3.0), 1e-6);
    }
}

}  // namespace
