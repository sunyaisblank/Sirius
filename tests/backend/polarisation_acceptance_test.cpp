// The actual tracer transports screens on accepted segments. The independent
// algebraic reference reconstructs their disk angle from conserved
// Walker-Penrose constants, using the test-only Boyer-Lindquist formulation.
// References: Walker & Penrose (1970); James et al. (2015), arXiv:1502.03808.
#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera_launch.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/oracle/polarisation_transport.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <complex>
#include <format>
#include <limits>
#include <numbers>

namespace {
using namespace sirius::core;
using namespace sirius::backend;
using Bl = sirius::oracle::Vec4d;

struct BlEvent {
    Bl position;
    Bl tangent;
};

// Differentiate the defining oblate-coordinate quartic and the azimuthal
// twist explicitly. No production coordinate adapter supplies the reference.
BlEvent ToBl(const Vec4& x, const Vec4& v, double a) {
    const double X = x(1), Y = x(2), Z = x(3);
    const double s = X * X + Y * Y + Z * Z - a * a;
    const double u = (s + std::sqrt(s * s + 4 * a * a * Z * Z)) / 2;
    const double r = std::sqrt(u);
    const double theta = std::acos(Z / r);
    const double ds = 2 * (X * v(1) + Y * v(2) + Z * v(3));
    const double dr = (u * ds + 2 * a * a * Z * v(3)) / ((2 * u - s) * 2 * r);
    const double dtheta = (std::cos(theta) * dr - v(3)) / (r * std::sin(theta));
    const double dphi = (X * v(2) - Y * v(1)) / (X * X + Y * Y) + a * dr / (r * r + a * a) -
                        a * dr / (r * r - 2 * r + a * a);
    return {{0, r, theta, std::atan2(Y, X) - std::atan2(a, r)},
            {v(0) - 2 * r * dr / (r * r - 2 * r + a * a), dr, dtheta, dphi}};
}

struct DiskReference {
    double angle, degree, intensity, redshift;
};

DiskReference DiskAngle(const CameraLaunch& launch, const CameraRay& camera,
                        const TraceResult& result, double a) {
    const std::array<double, 3> rest{camera.direction(1), camera.direction(2), camera.direction(3)};
    const auto screens = relativity::ObserverScreenBasis(launch.observer, rest);
    const auto initial = ToBl(launch.position, launch.tangent, a);
    std::array<std::complex<double>, 2> constants;
    for (int column = 0; column < 2; ++column) {
        const auto f = ToBl(launch.position, (*screens)[column], a).tangent;
        constants[column] = sirius::oracle::WalkerPenroseConstant(
            sirius::oracle::PolarisedStateD(initial.position, initial.tangent, f), a);
    }
    const auto terminal = ToBl(result.final_position, *result.final_tangent, a);
    sirius::oracle::KerrMetricD metric(1, a);
    double g[4][4], inv[4][4];
    metric.Evaluate(terminal.position, g, inv);
    const auto dot = [&](const Bl& lhs, const Bl& rhs) {
        return sirius::oracle::InnerProductD(g, lhs, rhs);
    };
    // Circular emitter, directly normalized with the independent BL metric.
    const double omega = 1 / (std::pow(terminal.position.r, 1.5) + a);
    Bl emitter(1, 0, 0, omega);
    emitter = emitter * (1 / std::sqrt(-dot(emitter, emitter)));
    const double frequency = dot(terminal.tangent, emitter);
    Bl n = terminal.tangent * (-1 / frequency) - emitter;
    Bl normal(0, 0, -1 / std::sqrt(g[2][2]), 0);
    const double mu = std::abs(dot(n, normal));
    Bl meridian = normal - n * dot(normal, n);
    meridian = meridian * (1 / std::sqrt(dot(meridian, meridian)));
    const auto endpoint = sirius::oracle::WalkerPenroseConstant(
        sirius::oracle::PolarisedStateD(terminal.position, terminal.tangent, meridian), a);
    const double det =
        constants[0].real() * constants[1].imag() - constants[1].real() * constants[0].imag();
    const double first =
        (endpoint.real() * constants[1].imag() - constants[1].real() * endpoint.imag()) / det;
    const double second =
        (constants[0].real() * endpoint.imag() - endpoint.real() * constants[0].imag()) / det;
    return {std::remainder(std::atan2(second, first) + std::numbers::pi / 2, std::numbers::pi),
            .1171 * (1 - mu) / (1 + 3.582 * mu), (1 + 2.06 * mu) / (1 + (2.0 / 3) * 2.06),
            1 / frequency};
}

TEST(PolarisationAcceptance, NormalEmissionHasFiniteUnpolarisedStokes) {
    KerrSchildFamily metric(KerrSchildParams::Schwarzschild(1));
    TracerConfig config;
    config.enable_disk = true;
    config.enable_polarisation = true;
    config.disk_inner = 6;
    config.disk_outer = 20;
    for (double tilt : {0., 5e-4, 1e-3}) {
        CameraRay ray;
        ray.origin(1) = 10;
        ray.origin(2) = std::numbers::pi / 2;
        ray.direction(1) = tilt;
        ray.direction(2) = std::sqrt(1 - tilt * tilt);
        ray.beta_forward = -.2;
        // Match the represented circular emitter, including its fp32 omega.
        ray.beta_right =
            static_cast<double>(static_cast<float>(std::pow(10., -1.5))) * 10 * std::sqrt(1.2);
        GeodesicTracer tracer(&metric, config);
        const auto result = tracer.Trace(ray);
        ASSERT_EQ(result.outcome, TraceResult::Outcome::DiskHit);
        ASSERT_FALSE(result.numerical_failure);
        ASSERT_EQ(result.num_disk_crossings, 1);
        const auto& crossing = result.disk_crossings[0];
        ASSERT_TRUE(crossing.polarisation_valid);
        EXPECT_TRUE(std::isfinite(crossing.polarisation_evpa));
        // The atmosphere law gives p=O(tilt^2); EVPA is unobservable at p=0.
        EXPECT_LE(crossing.polarisation_degree, tilt * tilt * .1);
        if (tilt == 0) {
            EXPECT_FLOAT_EQ(crossing.polarisation_degree, 0);
            EXPECT_FLOAT_EQ(crossing.polarisation_evpa, 0);
        }
    }
}

TEST(PolarisationAcceptance, ActualDiskStokesMatchWalkerPenroseReconstruction) {
    // Existing launch tests separately judge the measured input frame. This
    // comparison covers real accepted-segment transport and disk extraction;
    // independent full-path trajectory tests judge the central event itself.
    for (double a : {0., .7, .998}) {
        KerrSchildFamily metric(KerrSchildParams::Kerr(1, a));
        CameraRay ray;
        ray.origin(1) = 50;
        ray.origin(2) = std::numbers::pi / 3;
        ray.direction(1) = -std::sqrt(1 - .2 * .2 - .08 * .08);
        ray.direction(2) = .2;
        ray.direction(3) = .08;
        const auto launch = LaunchCameraRay(metric, a, ray);
        ASSERT_TRUE(launch);
        double previous_error = 1;
        double first_error = 0;
        for (int refinement = 0; refinement < 3; ++refinement) {
            SCOPED_TRACE(std::format("a={} refinement={}", a, refinement));
            TracerConfig config;
            config.enable_disk = true;
            config.enable_polarisation = true;
            config.disk_inner = 6;
            config.disk_outer = 20;
            config.integrator.initial_step = .02f;
            config.integrator.max_step = .4f / static_cast<float>(1 << refinement);
            config.integrator.abs_tolerance = 1e-7f / static_cast<float>(1 << (3 * refinement));
            config.integrator.rel_tolerance = config.integrator.abs_tolerance;
            GeodesicTracer tracer(&metric, config);
            const auto result = tracer.Trace(ray);
            ASSERT_EQ(result.outcome, TraceResult::Outcome::DiskHit);
            ASSERT_FALSE(result.numerical_failure);
            ASSERT_TRUE(result.final_tangent);
            const auto reference = DiskAngle(*launch, ray, result, a);
            const auto& crossing = result.disk_crossings[0];
            ASSERT_TRUE(crossing.polarisation_valid);
            const double error = std::abs(
                std::remainder(crossing.polarisation_evpa - reference.angle, std::numbers::pi));
            // Same 1e-4 live transport budget as the existing WP gate. The
            // thermal angle/degree public fields have binary32 precision.
            EXPECT_LT(error, 1e-4);
            if (refinement == 0) first_error = error;
            EXPECT_LE(error, previous_error + 2e-7);
            EXPECT_NEAR(crossing.polarisation_degree, reference.degree, 2e-7);
            EXPECT_NEAR(crossing.polarisation_intensity_scale, reference.intensity, 2e-6);
            EXPECT_NEAR(crossing.redshift, reference.redshift, 2e-6);
            for (bool sine : {false, true}) {
                const double observed = crossing.polarisation_degree *
                                        (sine ? std::sin(2 * crossing.polarisation_evpa)
                                              : std::cos(2 * crossing.polarisation_evpa));
                const double expected = reference.degree * (sine ? std::sin(2 * reference.angle)
                                                                 : std::cos(2 * reference.angle));
                EXPECT_NEAR(observed, expected, 2e-5);
            }
            const auto prefix = std::format("spin_{:.3f}_r{}", a, refinement);
            RecordProperty(prefix + "_evpa_error_rad", std::format("{:.17g}", error));
            RecordProperty(prefix + "_disk_radius", std::format("{:.17g}", crossing.r));
            previous_error = error;
        }
        // An error already below the public binary32 angle's rounding scale
        // cannot exhibit a reliable convergence slope. Above that scale,
        // refinement must reduce the independently measured error.
        const double angle_rounding = std::numeric_limits<float>::epsilon() * std::numbers::pi / 2;
        EXPECT_LE(previous_error, std::max(first_error / 2, angle_rounding));
        RecordProperty(std::format("spin_{:.3f}_refinement_behavior", a),
                       first_error <= angle_rounding ? "binary32_angle_plateau" : "error_reduced");
    }
}
}  // namespace
