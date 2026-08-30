// Production CPU Jacobi witnesses against independent Schwarzschild geometry.
// These are backend calculations only: no render session, image, dispatch, or
// presentation path is created.

#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/oracle/kerr_boyer_lindquist.h"

#include <gtest/gtest.h>

#include <algorithm>
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

}  // namespace
