// Mandatory Killing Vector Conservation Tests
// Purpose: MANDATORY tests for Killing vector conservation (moved from soft to mandatory)
// Ported from TSDG006A.cpp; assertions and tolerances unchanged.
//
// MATHEMATICAL BASIS:
// For stationary spacetimes (∂/∂t Killing vector):
//   E = -g_tμ k^μ = const (energy conservation)
//
// For axisymmetric spacetimes (∂/∂φ Killing vector):
//   L = g_φμ k^μ = const (angular momentum conservation)
//
// SPECIFICATION TARGETS (sirius/core/constants.h):
// - Energy conservation: |ΔE/E| < 10^-4
// - Angular momentum: |ΔL/L| < 10^-4
//
// LABEL: Mandatory (build gate)
// TESTS: Schwarzschild, Kerr, Reissner-Nordström

#include "sirius/core/geodesic_integrator.h"
#include "sirius/core/metrics/kerr_schild_family.h"  // Unified Kerr-Schild Family
#include "sirius/core/metrics/metric.h"
#include "sirius/core/tensor.h"

#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <vector>

namespace sirius::test {
using namespace sirius::core;

// =============================================================================
// MANDATORY Conservation Tolerances
// =============================================================================
// These tolerances define pass/fail criteria for build gating.
// Failures indicate fundamental integrator issues that must be fixed.

constexpr double MANDATORY_ENERGY_TOLERANCE = 1e-4;  // |ΔE/E| < 10^-4
constexpr double MANDATORY_L_TOLERANCE = 1e-4;       // |ΔL/L| < 10^-4
constexpr int INTEGRATION_STEPS = 500;               // Sufficient for meaningful test

// =============================================================================
// Test Fixture
// =============================================================================

class MandatoryKillingTests : public ::testing::Test {
  protected:
    IntegratorConfig config;

    void SetUp() override {
        config = Geodesic::GetDefaultConfig();
        config.abs_tolerance = 1e-8f;
        config.rel_tolerance = 1e-8f;
        config.min_step = 1e-8f;
        config.max_step = 0.05f;
        config.initial_step = 0.005f;
        config.use_rk45 = true;
    }

    // Compute Killing energy: E = -g_tμ k^μ
    double computeKillingEnergy(const Vec4& k, const Metric4d& g) {
        double E = 0.0;
        for (int mu = 0; mu < 4; ++mu) {
            E -= g(0, mu).real * k(mu);  // ξ^t = (1,0,0,0)
        }
        return E;
    }

    // Compute Killing angular momentum: L = ξ_φ^μ k_μ
    // In Cartesian: ξ^φ = x ∂_y - y ∂_x = (0, -y, x, 0)
    double computeKillingAngularMomentum(const Vec4& k, const Metric4d& g, const Vec4& pos) {
        double L = 0.0;
        Vec4 xi;
        xi(0) = 0.0;
        xi(1) = -pos(2);  // -y
        xi(2) = pos(1);   // x
        xi(3) = 0.0;

        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                L += g(mu, nu).real * k(mu) * xi(nu);
            }
        }
        return L;
    }

    // Initialize lightray from spherical to Cartesian
    Lightray createLightray(const Vec4& spherical_pos, double v_r, double v_theta, double v_phi,
                            IMetric* metric) {
        Lightray ray;

        double r = spherical_pos(1);
        double th = spherical_pos(2);
        double ph = spherical_pos(3);

        double sin_th = std::sin(th);
        double cos_th = std::cos(th);
        double sin_ph = std::sin(ph);
        double cos_ph = std::cos(ph);

        ray.position(0) = spherical_pos(0);
        ray.position(1) = r * sin_th * cos_ph;
        ray.position(2) = r * sin_th * sin_ph;
        ray.position(3) = r * cos_th;

        ray.proper_time = 0.0f;
        ray.coordinate_time = 0.0f;
        ray.terminated = 0;
        ray.step_size = config.initial_step;

        ray.velocity(1) =
            v_r * sin_th * cos_ph + r * v_theta * cos_th * cos_ph - r * v_phi * sin_th * sin_ph;
        ray.velocity(2) =
            v_r * sin_th * sin_ph + r * v_theta * cos_th * sin_ph + r * v_phi * sin_th * cos_ph;
        ray.velocity(3) = v_r * cos_th - r * v_theta * sin_th;

        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric->Evaluate(ray.position, g, dg);
        ray.velocity = TensorOps::NormalizeNull(ray.velocity, g);
        ray.acceleration = Geodesic::CalculateAcceleration(ray.velocity, ray.position, metric);

        return ray;
    }
};

// =============================================================================
// MANDATORY: Schwarzschild Energy Conservation
// =============================================================================

TEST_F(MandatoryKillingTests, SchwarzschildEnergyConservation) {
    // MANDATORY: Schwarzschild is stationary → E must be conserved
    sirius::core::KerrSchildFamily metric{sirius::core::KerrSchildParams::Schwarzschild(1.0)};

    Vec4 pos;
    pos(0) = 0.0;
    pos(1) = 20.0;  // Well outside horizon
    pos(2) = std::numbers::pi / 2.0;
    pos(3) = 0.0;

    Lightray ray = createLightray(pos, -0.5, 0.0, 0.5, &metric);

    Metric4d g_init;
    Tensor<Dual<double>, 4, 4, 4> dg_init;
    metric.Evaluate(ray.position, g_init, dg_init);
    double E_initial = computeKillingEnergy(ray.velocity, g_init);
    ASSERT_GT(std::abs(E_initial), 1e-10)
        << "deterministic witness produced degenerate initial energy";

    double max_drift = 0.0;
    int accepted_steps = 0;

    for (int step = 0; step < INTEGRATION_STEPS; ++step) {
        bool success = Geodesic::IntegrateStepRk45(ray, &metric, config);
        if (!success || ray.terminated) break;
        ++accepted_steps;

        double r = std::sqrt(ray.position(1) * ray.position(1) + ray.position(2) * ray.position(2) +
                             ray.position(3) * ray.position(3));
        if (r < 2.5 || r > 100.0) break;

        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.Evaluate(ray.position, g, dg);
        double E_current = computeKillingEnergy(ray.velocity, g);

        double drift = std::abs((E_current - E_initial) / E_initial);
        max_drift = std::max(max_drift, drift);
    }

    ASSERT_GT(accepted_steps, 0) << "conservation witness accepted no integration step";
    EXPECT_LT(max_drift, MANDATORY_ENERGY_TOLERANCE)
        << "MANDATORY: Schwarzschild energy drift " << max_drift << " exceeds tolerance "
        << MANDATORY_ENERGY_TOLERANCE;
}

// =============================================================================
// MANDATORY: Schwarzschild Angular Momentum Conservation
// =============================================================================

TEST_F(MandatoryKillingTests, SchwarzschildAngularMomentumConservation) {
    // MANDATORY: Schwarzschild is spherically symmetric → L must be conserved
    sirius::core::KerrSchildFamily metric{sirius::core::KerrSchildParams::Schwarzschild(1.0)};

    Vec4 pos;
    pos(0) = 0.0;
    pos(1) = 15.0;
    pos(2) = std::numbers::pi / 2.0;  // Equatorial
    pos(3) = 0.0;

    Lightray ray = createLightray(pos, -0.3, 0.0, 0.8, &metric);

    Metric4d g_init;
    Tensor<Dual<double>, 4, 4, 4> dg_init;
    metric.Evaluate(ray.position, g_init, dg_init);
    double L_initial = computeKillingAngularMomentum(ray.velocity, g_init, ray.position);

    ASSERT_GT(std::abs(L_initial), 1e-10)
        << "deterministic witness produced degenerate initial angular momentum";

    double max_drift = 0.0;
    int accepted_steps = 0;

    for (int step = 0; step < INTEGRATION_STEPS; ++step) {
        bool success = Geodesic::IntegrateStepRk45(ray, &metric, config);
        if (!success || ray.terminated) break;
        ++accepted_steps;

        double r = std::sqrt(ray.position(1) * ray.position(1) + ray.position(2) * ray.position(2) +
                             ray.position(3) * ray.position(3));
        if (r < 2.5 || r > 100.0) break;

        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.Evaluate(ray.position, g, dg);
        double L_current = computeKillingAngularMomentum(ray.velocity, g, ray.position);

        double drift = std::abs((L_current - L_initial) / L_initial);
        max_drift = std::max(max_drift, drift);
    }

    ASSERT_GT(accepted_steps, 0) << "conservation witness accepted no integration step";
    EXPECT_LT(max_drift, MANDATORY_L_TOLERANCE) << "MANDATORY: Schwarzschild L drift " << max_drift
                                                << " exceeds tolerance " << MANDATORY_L_TOLERANCE;
}

// =============================================================================
// MANDATORY: Kerr Energy Conservation
// =============================================================================

TEST_F(MandatoryKillingTests, KerrEnergyConservation) {
    // MANDATORY: Kerr is stationary → E must be conserved
    sirius::core::KerrSchildFamily metric{sirius::core::KerrSchildParams::Kerr(1.0, 0.9)};

    Vec4 pos;
    pos(0) = 0.0;
    pos(1) = 15.0;
    pos(2) = std::numbers::pi / 3.0;  // Off-equatorial
    pos(3) = 0.0;

    Lightray ray = createLightray(pos, -0.4, 0.1, 0.5, &metric);

    Metric4d g_init;
    Tensor<Dual<double>, 4, 4, 4> dg_init;
    metric.Evaluate(ray.position, g_init, dg_init);
    double E_initial = computeKillingEnergy(ray.velocity, g_init);
    ASSERT_GT(std::abs(E_initial), 1e-10)
        << "deterministic witness produced degenerate initial energy";

    double max_drift = 0.0;
    int accepted_steps = 0;

    for (int step = 0; step < INTEGRATION_STEPS; ++step) {
        bool success = Geodesic::IntegrateStepRk45(ray, &metric, config);
        if (!success || ray.terminated) break;
        ++accepted_steps;

        double r = std::sqrt(ray.position(1) * ray.position(1) + ray.position(2) * ray.position(2) +
                             ray.position(3) * ray.position(3));
        if (r < 2.0 || r > 100.0) break;

        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.Evaluate(ray.position, g, dg);
        double E_current = computeKillingEnergy(ray.velocity, g);

        double drift = std::abs((E_current - E_initial) / E_initial);
        max_drift = std::max(max_drift, drift);
    }

    ASSERT_GT(accepted_steps, 0) << "conservation witness accepted no integration step";
    EXPECT_LT(max_drift, MANDATORY_ENERGY_TOLERANCE)
        << "MANDATORY: Kerr energy drift " << max_drift << " exceeds tolerance "
        << MANDATORY_ENERGY_TOLERANCE;
}

// =============================================================================
// MANDATORY: Kerr Angular Momentum Conservation
// =============================================================================

TEST_F(MandatoryKillingTests, KerrAngularMomentumConservation) {
    // MANDATORY: Kerr is axisymmetric → L_z must be conserved
    sirius::core::KerrSchildFamily metric{sirius::core::KerrSchildParams::Kerr(1.0, 0.7)};

    Vec4 pos;
    pos(0) = 0.0;
    pos(1) = 12.0;
    pos(2) = std::numbers::pi / 2.0;
    pos(3) = 0.0;

    Lightray ray = createLightray(pos, -0.3, 0.0, 0.7, &metric);

    Metric4d g_init;
    Tensor<Dual<double>, 4, 4, 4> dg_init;
    metric.Evaluate(ray.position, g_init, dg_init);
    double L_initial = computeKillingAngularMomentum(ray.velocity, g_init, ray.position);

    ASSERT_GT(std::abs(L_initial), 1e-10)
        << "deterministic witness produced degenerate initial angular momentum";

    double max_drift = 0.0;
    int accepted_steps = 0;

    for (int step = 0; step < INTEGRATION_STEPS; ++step) {
        bool success = Geodesic::IntegrateStepRk45(ray, &metric, config);
        if (!success || ray.terminated) break;
        ++accepted_steps;

        double r = std::sqrt(ray.position(1) * ray.position(1) + ray.position(2) * ray.position(2) +
                             ray.position(3) * ray.position(3));
        if (r < 2.0 || r > 100.0) break;

        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.Evaluate(ray.position, g, dg);
        double L_current = computeKillingAngularMomentum(ray.velocity, g, ray.position);

        double drift = std::abs((L_current - L_initial) / L_initial);
        max_drift = std::max(max_drift, drift);
    }

    ASSERT_GT(accepted_steps, 0) << "conservation witness accepted no integration step";
    EXPECT_LT(max_drift, MANDATORY_L_TOLERANCE) << "MANDATORY: Kerr L drift " << max_drift
                                                << " exceeds tolerance " << MANDATORY_L_TOLERANCE;
}

// =============================================================================
// MANDATORY: Reissner-Nordström Energy Conservation
// =============================================================================

TEST_F(MandatoryKillingTests, ReissnerNordstromEnergyConservation) {
    // MANDATORY: RN is stationary → E must be conserved
    sirius::core::KerrSchildFamily metric{
        sirius::core::KerrSchildParams::ReissnerNordstrom(1.0, 0.5)};

    Vec4 pos;
    pos(0) = 0.0;
    pos(1) = 15.0;
    pos(2) = std::numbers::pi / 2.0;
    pos(3) = 0.0;

    Lightray ray = createLightray(pos, -0.5, 0.0, 0.4, &metric);

    Metric4d g_init;
    Tensor<Dual<double>, 4, 4, 4> dg_init;
    metric.Evaluate(ray.position, g_init, dg_init);
    double E_initial = computeKillingEnergy(ray.velocity, g_init);
    ASSERT_GT(std::abs(E_initial), 1e-10)
        << "deterministic witness produced degenerate initial energy";

    double max_drift = 0.0;
    int accepted_steps = 0;

    for (int step = 0; step < INTEGRATION_STEPS; ++step) {
        bool success = Geodesic::IntegrateStepRk45(ray, &metric, config);
        if (!success || ray.terminated) break;
        ++accepted_steps;

        double r = std::sqrt(ray.position(1) * ray.position(1) + ray.position(2) * ray.position(2) +
                             ray.position(3) * ray.position(3));
        if (r < 2.0 || r > 100.0) break;

        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.Evaluate(ray.position, g, dg);
        double E_current = computeKillingEnergy(ray.velocity, g);

        double drift = std::abs((E_current - E_initial) / E_initial);
        max_drift = std::max(max_drift, drift);
    }

    ASSERT_GT(accepted_steps, 0) << "conservation witness accepted no integration step";
    EXPECT_LT(max_drift, MANDATORY_ENERGY_TOLERANCE)
        << "MANDATORY: RN energy drift " << max_drift << " exceeds tolerance "
        << MANDATORY_ENERGY_TOLERANCE;
}

}  // namespace sirius::test
