// TSVAL001A.cpp - DNGR Validation and Analytic Tests
// Component ID: TSVAL001A
// Purpose: Validate against analytic solutions from DNGR paper
//
// REFERENCE: James et al. (2015) "Gravitational Lensing by Spinning Black Holes"
// arXiv: 1502.03808

#include "sirius/core/disk/novikov_thorne_disk.h"
#include "sirius/oracle/beam_integrator.h"
#include "sirius/oracle/kerr_boyer_lindquist.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <numbers>

using namespace sirius::core;
using namespace sirius::oracle;

namespace {

double IndependentSpecificEnergy(double radius, double spin) {
    const double sqrt_radius = std::sqrt(radius);
    const double root = std::sqrt(radius * radius - 3.0 * radius + 2.0 * spin * sqrt_radius);
    return (radius * radius - 2.0 * radius + spin * sqrt_radius) / (radius * root);
}

double IndependentSpecificAngularMomentum(double radius, double spin) {
    const double sqrt_radius = std::sqrt(radius);
    const double root = std::sqrt(radius * radius - 3.0 * radius + 2.0 * spin * sqrt_radius);
    return (radius * radius - 2.0 * spin * sqrt_radius + spin * spin) / (sqrt_radius * root);
}

double IndependentAngularVelocity(double radius, double spin) {
    return 1.0 / (std::pow(radius, 1.5) + spin);
}

double IndependentAngularMomentumDerivative(double radius, double spin) {
    const double step = 1.0e-4 * radius;
    return (-IndependentSpecificAngularMomentum(radius + 2.0 * step, spin) +
            8.0 * IndependentSpecificAngularMomentum(radius + step, spin) -
            8.0 * IndependentSpecificAngularMomentum(radius - step, spin) +
            IndependentSpecificAngularMomentum(radius - 2.0 * step, spin)) /
           (12.0 * step);
}

double IndependentAngularVelocityDerivative(double radius, double spin) {
    const double step = 1.0e-4 * radius;
    return (-IndependentAngularVelocity(radius + 2.0 * step, spin) +
            8.0 * IndependentAngularVelocity(radius + step, spin) -
            8.0 * IndependentAngularVelocity(radius - step, spin) +
            IndependentAngularVelocity(radius - 2.0 * step, spin)) /
           (12.0 * step);
}

// Page & Thorne (1974), Eq. 15n. This uses a composite midpoint rule and
// finite-difference derivatives, independent of Core's 16-point Gauss-Legendre
// rule and analytic derivatives. The returned shape is Q(r)/r^3; its physical
// prefactor cancels when temperature profiles are normalised.
double IndependentPageThorneFluxShape(double radius, double spin) {
    const double isco = AccretionDiskD::ComputeIsco(spin);
    if (radius <= isco) return 0.0;

    constexpr int kPanels = 32768;
    const double width = (radius - isco) / kPanels;
    double integral = 0.0;
    for (int panel = 0; panel < kPanels; ++panel) {
        const double sample_radius = isco + (static_cast<double>(panel) + 0.5) * width;
        const double energy = IndependentSpecificEnergy(sample_radius, spin);
        const double angular_momentum = IndependentSpecificAngularMomentum(sample_radius, spin);
        const double angular_velocity = IndependentAngularVelocity(sample_radius, spin);
        integral += (energy - angular_velocity * angular_momentum) *
                    IndependentAngularMomentumDerivative(sample_radius, spin);
    }
    integral *= width;

    const double energy = IndependentSpecificEnergy(radius, spin);
    const double angular_momentum = IndependentSpecificAngularMomentum(radius, spin);
    const double angular_velocity = IndependentAngularVelocity(radius, spin);
    const double invariant = energy - angular_velocity * angular_momentum;
    const double correction =
        -IndependentAngularVelocityDerivative(radius, spin) * integral / (invariant * invariant);
    return correction / (radius * radius * radius);
}

double NormalisedOracleTemperature(double radius, double reference_radius, double spin) {
    return std::pow(IndependentPageThorneFluxShape(radius, spin) /
                        IndependentPageThorneFluxShape(reference_radius, spin),
                    0.25);
}

void ExpectCircularNullOrbit(KerrMetricD& metric, double radius, double impact_parameter) {
    const Vec4d position(0.0, radius, std::numbers::pi / 2.0, 0.0);
    const Vec4d momentum(-1.0, 0.0, 0.0, impact_parameter);
    EXPECT_NEAR(metric.Hamiltonian(position, momentum), 0.0, 2.0e-13)
        << "null condition at r/M=" << radius;
    EXPECT_NEAR(metric.dHdq(position, momentum).r, 0.0, 2.0e-13)
        << "radial canonical acceleration at r/M=" << radius;
}

//==============================================================================
// Photon Sphere Radius Tests (DNGR Eq. 5)
//==============================================================================

TEST(AnalyticValidationTest, PhotonSphereSchwarzschildExact) {
    KerrMetricD metric(1.0, 0.0);
    ASSERT_NO_FATAL_FAILURE(ExpectCircularNullOrbit(metric, 3.0, 3.0 * std::sqrt(3.0)));
}

TEST(AnalyticValidationTest, PhotonSphereKerrPrograde) {
    constexpr double spin = 0.9;
    const double radius = 2.0 * (1.0 + std::cos((2.0 / 3.0) * std::acos(-spin)));
    const double impact_parameter =
        (radius * radius * (radius - 3.0) + spin * spin * (radius + 1.0)) / (spin * (1.0 - radius));
    EXPECT_NEAR(radius, 1.5578546274233829, 1.0e-14);
    KerrMetricD metric(1.0, spin);
    ExpectCircularNullOrbit(metric, radius, impact_parameter);
}

TEST(AnalyticValidationTest, PhotonSphereKerrRetrograde) {
    constexpr double spin = 0.9;
    const double radius = 2.0 * (1.0 + std::cos((2.0 / 3.0) * std::acos(spin)));
    const double impact_parameter =
        (radius * radius * (radius - 3.0) + spin * spin * (radius + 1.0)) / (spin * (1.0 - radius));
    EXPECT_NEAR(radius, 3.910267939103037, 1.0e-14);
    KerrMetricD metric(1.0, spin);
    ExpectCircularNullOrbit(metric, radius, impact_parameter);
}

//==============================================================================
// ISCO Location Tests (verified in TSDK001A, repeated for completeness)
//==============================================================================

TEST(AnalyticValidationTest, ISCOSchwarzschildExact) {
    // For Schwarzschild: r_ISCO = 6M exactly
    AccretionDiskD::Config config;
    config.M = 1.0;
    config.a_star = 0.0;
    AccretionDiskD disk(config);

    EXPECT_NEAR(disk.IscoRadius(), 6.0, 1.0e-14) << "Schwarzschild ISCO should be at r = 6M";
}

TEST(AnalyticValidationTest, ISCOKerrPrograde) {
    // For a/M = 0.9: r_ISCO ≈ 2.32M (prograde)
    AccretionDiskD::Config config;
    config.M = 1.0;
    config.a_star = 0.9;
    AccretionDiskD disk(config);

    double r_isco = disk.IscoRadius();

    EXPECT_NEAR(r_isco, 2.320883041761887, 1.0e-14);
}

TEST(AnalyticValidationTest, ISCONearExtremal) {
    // For a/M → 1: r_ISCO → M (prograde)
    AccretionDiskD::Config config;
    config.M = 1.0;
    config.a_star = 0.998;
    AccretionDiskD disk(config);

    double r_isco = disk.IscoRadius();

    EXPECT_NEAR(r_isco, 1.2369706551751847, 1.0e-14);
}

TEST(AnalyticValidationTest, PageThorneFluxMatchesIndependentQuadrature) {
    for (const double spin : {-0.7, 0.0, 0.9}) {
        AccretionDiskD::Config config;
        config.a_star = spin;
        AccretionDiskD disk(config);
        const double reference_radius = 1.5 * disk.IscoRadius();
        const double reference_temperature = disk.Temperature(reference_radius);
        ASSERT_GT(reference_temperature, 0.0);
        for (const double radius_scale : {1.05, 1.2, 2.0, 4.0}) {
            const double radius = radius_scale * disk.IscoRadius();
            const double production = disk.Temperature(radius) / reference_temperature;
            const double oracle = NormalisedOracleTemperature(radius, reference_radius, spin);
            EXPECT_NEAR(production / oracle, 1.0, 2.0e-5) << "spin=" << spin << ", r/M=" << radius;
        }
    }
}

TEST(AnalyticValidationTest, PageThorneTemperatureHasZeroTorqueInnerEdge) {
    AccretionDiskD disk;
    const double isco = disk.IscoRadius();
    EXPECT_EQ(disk.Temperature(isco), 0.0);
    EXPECT_EQ(IndependentPageThorneFluxShape(isco, 0.0), 0.0);
    EXPECT_GT(disk.Temperature(1.05 * isco), 0.0);
    EXPECT_GT(disk.Temperature(1.5 * isco), disk.Temperature(1.05 * isco));
}

//==============================================================================
// DNGR Configuration Tests (a/M = 0.999, r_cam = 6.03M)
//==============================================================================

TEST(DNGRParityTest, ExtremalKerrConfiguration) {
    // DNGR paper uses a/M = 0.999
    double a = 0.999;
    double M = 1.0;

    KerrMetricD metric(M, a);

    // Horizon should be very close to r = M
    // r_+ = M + √(M² - a²) = 1 + √(1 - 0.998001) = 1 + 0.0447 ≈ 1.045
    const double analytic_horizon = M + std::sqrt(M * M - a * a);
    EXPECT_NEAR(metric.HorizonRadius(), analytic_horizon, 1.0e-14);
}

TEST(DNGRParityTest, CameraDistance) {
    // DNGR uses camera at r = 6.03M
    double r_camera = 6.03;
    double M = 1.0;

    // Camera should be outside ISCO for a=0.999
    AccretionDiskD::Config config;
    config.M = M;
    config.a_star = 0.999;
    AccretionDiskD disk(config);

    double r_isco = disk.IscoRadius();

    EXPECT_GT(r_camera, r_isco) << "Camera at 6.03M should be outside ISCO for a=0.999";
}

//==============================================================================
// Numerical Stability Tests
//==============================================================================

TEST(NumericalStabilityTest, NoNaNInMetric) {
    KerrMetricD metric(1.0, 0.999);

    // Test at various locations
    std::vector<Vec4d> positions = {
        Vec4d(0, 100, std::numbers::pi / 2, 0),   // Far field
        Vec4d(0, 10, std::numbers::pi / 4, 1.0),  // Moderate distance
        Vec4d(0, 2.0, std::numbers::pi / 2, 0),   // Near horizon
        Vec4d(0, 6.03, 1.449, 0),                 // DNGR camera position
    };

    for (const auto& pos : positions) {
        double g[4][4], g_inv[4][4];
        metric.Evaluate(pos, g, g_inv);

        // Check for NaN
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                EXPECT_FALSE(std::isnan(g[i][j])) << "Metric should not contain NaN at r=" << pos.r;
                EXPECT_FALSE(std::isinf(g[i][j])) << "Metric should not contain Inf at r=" << pos.r;
            }
        }
    }
}

TEST(NumericalStabilityTest, DeterministicIntegration) {
    KerrMetricD metric(1.0, 0.9);
    BeamIntegratorD::Config config;
    BeamIntegratorD integrator(&metric, config);

    // Run integration twice with same initial conditions
    auto runIntegration = [&]() {
        BeamStateD beam;
        beam.x = Vec4d(0, 100, std::numbers::pi / 2, 0);
        beam.k.t = -1.0;
        beam.k.r = -0.1;
        beam.k.theta = 0;
        beam.k.phi = 0.05;
        beam.Initialise();

        for (int i = 0; i < 100 && !beam.terminated; ++i) {
            integrator.Step(beam, 0.1);
        }

        return beam.x.r;
    };

    double r1 = runIntegration();
    double r2 = runIntegration();

    EXPECT_DOUBLE_EQ(r1, r2) << "Integration should be deterministic";
}

}  // namespace
