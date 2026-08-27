// Symplectic integrator oracle: canonical two-form, Hamiltonian, energy, and
// angular-momentum conservation; horizon/escape termination; Yoshida order
// comparison. The fixed-step structural gate disables null projection because
// projection is intentionally outside the symplectic map.

#include "sirius/oracle/symplectic_integrator.h"

#include "sirius/oracle/adaptive_hamiltonian_integrator.h"
#include "sirius/oracle/kerr_boyer_lindquist.h"
#include "sirius/oracle/metric_interface.h"
#include "sirius/oracle/transport_types.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>

using namespace sirius::oracle;

namespace {

// Fixture for symplectic integrator tests
class SymplecticIntegratorTest : public ::testing::Test {
  protected:
    void SetUp() override {
        schwarzschild = std::make_unique<KerrMetricD>(1.0, 0.0);
        kerr = std::make_unique<KerrMetricD>(1.0, 0.5);

        SymplecticIntegratorD::Config config;
        config.initial_step_size = 0.05;
        config.tolerance = 1e-10;

        integrator_sch = std::make_unique<SymplecticIntegratorD>(schwarzschild.get(), config);
        integrator_kerr = std::make_unique<SymplecticIntegratorD>(kerr.get(), config);
    }

    // Create a null ray at given position - use NEGATIVE p_r for ingoing
    GeodesicStateD createIngoingNullRay(const IMetricD* metric, double r, double theta,
                                        double angular_momentum = 0.1) {
        Vec4d x(0, r, theta, 0);

        double g[4][4], g_inv[4][4];
        metric->Evaluate(x, g, g_inv);

        Vec4d p;
        p.t = -1.0;  // E = 1 (E = -p_t)
        p.phi = angular_momentum;
        p.theta = 0;

        // Solve for |p_r| from null condition: g^μν p_μ p_ν = 0
        double A = g_inv[1][1];
        double C =
            g_inv[0][0] * p.t * p.t + 2 * g_inv[0][3] * p.t * p.phi + g_inv[3][3] * p.phi * p.phi;

        if (A > 0 && C < 0) {
            // Use negative p_r for ingoing (dr/dλ = g^rr * p_r < 0)
            p.r = -std::sqrt(-C / A);
        } else {
            p.r = 0;
        }

        GeodesicStateD state(x, p);
        state.ComputeConservedQuantities(0);  // a=0 for Schwarzschild
        return state;
    }

    std::unique_ptr<KerrMetricD> schwarzschild;
    std::unique_ptr<KerrMetricD> kerr;
    std::unique_ptr<SymplecticIntegratorD> integrator_sch;
    std::unique_ptr<SymplecticIntegratorD> integrator_kerr;
};

//==============================================================================
// Test: Hamiltonian Conservation Over 10^4 Steps
// Verifies: |H| < 10^-10 after 10000 steps
//==============================================================================

TEST_F(SymplecticIntegratorTest, HamiltonianConservation_10kSteps) {
    constexpr double tolerance = 1.0e-10;
    const int numSteps = 10000;
    const double step_size = 0.05;

    // Create a ray from r=100M with high angular momentum (scattered ray)
    // This will orbit/scatter rather than plunge into the horizon
    Vec4d x(0, 100.0, std::numbers::pi / 2, 0);

    double g[4][4], g_inv[4][4];
    schwarzschild->Evaluate(x, g, g_inv);

    Vec4d p;
    p.t = -1.0;    // E = 1
    p.phi = 10.0;  // High angular momentum to prevent plunge
    p.theta = 0;

    // Solve for p_r from null condition (ingoing)
    double A = g_inv[1][1];
    double C =
        g_inv[0][0] * p.t * p.t + 2 * g_inv[0][3] * p.t * p.phi + g_inv[3][3] * p.phi * p.phi;
    if (A > 0 && C < 0) {
        p.r = -std::sqrt(-C / A);  // Negative = ingoing
    }

    GeodesicStateD initial(x, p);

    auto stats = integrator_sch->TestConservation(initial, numSteps, step_size);

    EXPECT_LT(stats.max_h_error, tolerance)
        << "Hamiltonian error exceeded tolerance after " << stats.steps << " steps";
    EXPECT_GT(stats.steps, numSteps * 0.5)
        << "Integration terminated early after " << stats.steps << " steps";
}

//==============================================================================
// Test: Energy Conservation (E = -p_t)
// Verifies: |ΔE|/E < 10^-10 after integration
//==============================================================================

TEST_F(SymplecticIntegratorTest, EnergyConservation) {
    constexpr double tolerance = 1.0e-10;
    const int numSteps = 5000;
    const double step_size = 0.05;

    // b=L/E=6 exceeds the Schwarzschild capture threshold 3 sqrt(3), so this
    // is a long scattering ray rather than a vacuous horizon termination.
    GeodesicStateD initial =
        createIngoingNullRay(schwarzschild.get(), 15.0, std::numbers::pi / 3, 6.0);

    auto stats = integrator_sch->TestConservation(initial, numSteps, step_size);

    double deltaE = std::abs(stats.final_e - stats.initial_e);
    double relError = deltaE / std::abs(stats.initial_e);

    EXPECT_LT(relError, tolerance)
        << "Relative energy drift: " << relError << " (ΔE = " << deltaE << ")";
    EXPECT_GT(stats.steps, numSteps * 0.5)
        << "energy witness terminated before a meaningful trajectory";
}

//==============================================================================
// Test: Angular Momentum Conservation (Lz = p_φ)
// Verifies: |ΔLz|/Lz < 10^-10 after integration
//==============================================================================

TEST_F(SymplecticIntegratorTest, AngularMomentumConservation) {
    constexpr double tolerance = 1.0e-10;
    const int numSteps = 5000;
    const double step_size = 0.05;

    GeodesicStateD initial =
        createIngoingNullRay(schwarzschild.get(), 12.0, std::numbers::pi / 2, 6.0);

    auto stats = integrator_sch->TestConservation(initial, numSteps, step_size);

    double deltaLz = std::abs(stats.final_lz - stats.initial_lz);
    double relError = deltaLz / std::abs(stats.initial_lz);

    EXPECT_LT(relError, tolerance)
        << "Relative Lz drift: " << relError << " (ΔLz = " << deltaLz << ")";
    EXPECT_GT(stats.steps, numSteps * 0.5)
        << "angular-momentum witness terminated before a meaningful trajectory";
}

//==============================================================================
// Test: Symplectic Structure (Phase Space Volume)
// Verifies: Integration preserves symplectic 2-form (area in phase space)
//==============================================================================

TEST_F(SymplecticIntegratorTest, SymplecticStructurePreserved) {
    SymplecticIntegratorD::Config config;
    config.order = IntegratorOrder::kYoshida4;
    config.enforce_null_condition = false;
    SymplecticIntegratorD integrator(schwarzschild.get(), config);
    const HamiltonianStateD initial(createIngoingNullRay(schwarzschild.get(), 10.0, 1.1, 6.0));

    constexpr double step_size = 0.01;
    constexpr double perturbation = 1.0e-6;
    std::array<std::array<double, 8>, 8> jacobian{};
    for (int column = 0; column < 8; ++column) {
        HamiltonianStateD plus = initial;
        HamiltonianStateD minus = initial;
        if (column < 4) {
            plus.q[column] += perturbation;
            minus.q[column] -= perturbation;
        } else {
            plus.p[column - 4] += perturbation;
            minus.p[column - 4] -= perturbation;
        }
        const HamiltonianStateD plus_result = integrator.Step(plus, step_size).state;
        const HamiltonianStateD minus_result = integrator.Step(minus, step_size).state;
        for (int row = 0; row < 4; ++row) {
            jacobian[row][column] =
                (plus_result.q[row] - minus_result.q[row]) / (2.0 * perturbation);
            jacobian[row + 4][column] =
                (plus_result.p[row] - minus_result.p[row]) / (2.0 * perturbation);
        }
    }

    // A phase map Phi is symplectic exactly when D(Phi)^T J D(Phi) = J.
    double worst_residual = 0.0;
    for (int row = 0; row < 8; ++row) {
        for (int column = 0; column < 8; ++column) {
            double transported = 0.0;
            for (int axis = 0; axis < 4; ++axis) {
                transported += jacobian[axis][row] * jacobian[axis + 4][column] -
                               jacobian[axis + 4][row] * jacobian[axis][column];
            }
            double canonical = 0.0;
            if (column == row + 4) canonical = 1.0;
            if (row == column + 4) canonical = -1.0;
            worst_residual = std::max(worst_residual, std::abs(transported - canonical));
        }
    }
    EXPECT_LT(worst_residual, 5.0e-7)
        << "fixed-step phase map does not preserve the canonical two-form";
}

//==============================================================================
// Test: Horizon Termination
// Verifies: Integration correctly detects and terminates at horizon
//==============================================================================

TEST_F(SymplecticIntegratorTest, HorizonTermination) {
    // Start at r=3M with ingoing null ray (p_r < 0 for ingoing)
    GeodesicStateD initial = createIngoingNullRay(schwarzschild.get(), 3.0, std::numbers::pi / 2);

    auto result = integrator_sch->Integrate(initial, 100.0);

    // Should hit horizon or move inward
    EXPECT_TRUE(result.hit_horizon || result.final_state.x.r < 3.0)
        << "Should have hit horizon or moved inward, r=" << result.final_state.x.r;
}

//==============================================================================
// Test: Escape Detection
// Verifies: Integration correctly detects rays escaping to infinity
//==============================================================================

TEST_F(SymplecticIntegratorTest, EscapeDetection) {
    // Start at r=50M with OUTgoing null ray (positive p_r)
    Vec4d x(0, 50.0, std::numbers::pi / 2, 0);

    double g[4][4], g_inv[4][4];
    schwarzschild->Evaluate(x, g, g_inv);

    Vec4d p;
    p.t = -1.0;  // E = 1
    p.phi = 0;   // Radial ray
    p.theta = 0;

    // Solve for p_r, use POSITIVE for outgoing
    double A = g_inv[1][1];
    double C = g_inv[0][0] * p.t * p.t;
    if (A > 0 && C < 0) {
        p.r = std::sqrt(-C / A);  // Positive = outgoing
    }

    GeodesicStateD initial(x, p);

    auto result = integrator_sch->Integrate(initial, 10000.0, 100.0);

    EXPECT_TRUE(result.escaped || result.final_state.x.r > 50.0)
        << "Should have escaped or moved outward, r=" << result.final_state.x.r;
}

//==============================================================================
// Test: Kerr Metric Integration
// Verifies: Integrator works correctly with spinning black hole
//==============================================================================

TEST_F(SymplecticIntegratorTest, KerrMetricIntegration) {
    // Tolerance for Kerr metric with fully analytic dHdq
    // Measured: ~1.17e-3 over 2k steps
    const double tolerance = 2e-3;
    const int numSteps = 2000;
    const double step_size = 0.05;

    Vec4d x(0, 8.0, std::numbers::pi / 3, 0);

    double g[4][4], g_inv[4][4];
    kerr->Evaluate(x, g, g_inv);

    Vec4d p;
    p.t = -1.0;
    p.phi = 0.3;
    p.theta = 0;

    // Null condition
    double A = g_inv[1][1];
    double C =
        g_inv[0][0] * p.t * p.t + 2 * g_inv[0][3] * p.t * p.phi + g_inv[3][3] * p.phi * p.phi;
    if (A > 0 && C < 0) {
        p.r = std::sqrt(-C / A);
    }

    GeodesicStateD initial(x, p);

    auto stats = integrator_kerr->TestConservation(initial, numSteps, step_size);

    EXPECT_LT(stats.max_h_error, tolerance) << "Kerr Hamiltonian error: " << stats.max_h_error;
}

TEST(AnalyticValidationTest, NearExtremalKerrConservesEnergyAngularMomentumAndCarter) {
    KerrMetricD near_extremal(1.0, 0.998);
    AdaptiveHamiltonianIntegratorD::Config config;
    config.absolute_tolerance = 1.0e-13;
    config.relative_tolerance = 1.0e-13;
    AdaptiveHamiltonianIntegratorD integrator(&near_extremal, config);

    Vec4d position(0.0, 20.0, 1.1, 0.0);
    double metric[4][4], inverse[4][4];
    near_extremal.Evaluate(position, metric, inverse);
    Vec4d momentum;
    momentum.t = -1.0;
    momentum.theta = 0.4;
    momentum.phi = 3.0;
    const double non_radial = inverse[0][0] * momentum.t * momentum.t +
                              2.0 * inverse[0][3] * momentum.t * momentum.phi +
                              inverse[2][2] * momentum.theta * momentum.theta +
                              inverse[3][3] * momentum.phi * momentum.phi;
    ASSERT_LT(non_radial, 0.0);
    momentum.r = std::sqrt(-non_radial / inverse[1][1]);  // Outgoing complete ray.

    GeodesicStateD initial(position, momentum);
    initial.ComputeConservedQuantities(0.998);
    HamiltonianStateD state(initial);
    const double initial_energy = initial.E;
    const double initial_angular_momentum = initial.Lz;
    const double initial_carter = initial.Q;
    double max_energy_drift = 0.0;
    double max_angular_momentum_drift = 0.0;
    double max_carter_drift = 0.0;

    int accepted_steps = 0;
    constexpr double escape_radius = 200.0;
    for (int segment = 0; segment < 250 && state.q.r < escape_radius; ++segment) {
        const auto result = integrator.Integrate(state, 1.0);
        ASSERT_TRUE(result.complete);
        state = result.state;
        accepted_steps += result.accepted_steps;
        GeodesicStateD current = state.ToGeodesicState();
        current.ComputeConservedQuantities(0.998);
        max_energy_drift = std::max(max_energy_drift, std::abs(current.E - initial_energy));
        max_angular_momentum_drift =
            std::max(max_angular_momentum_drift, std::abs(current.Lz - initial_angular_momentum) /
                                                     std::abs(initial_angular_momentum));
        max_carter_drift =
            std::max(max_carter_drift, std::abs(current.Q - initial_carter) /
                                           std::max(std::abs(initial_carter), 1.0e-30));
    }

    std::cout << "near-extremal oracle drift E=" << max_energy_drift
              << " Lz=" << max_angular_momentum_drift << " Q=" << max_carter_drift
              << " final_r=" << state.q.r << " steps=" << accepted_steps << '\n';
    EXPECT_GE(state.q.r, escape_radius)
        << "the oracle witness did not span the complete ray to the escape surface";
    EXPECT_LT(max_energy_drift, 1.0e-10);
    EXPECT_LT(max_angular_momentum_drift, 1.0e-10);
    EXPECT_LT(max_carter_drift, 1.0e-10);
}

//==============================================================================
// Test: Single Step Accuracy
// Verifies: Individual Step preserves null condition
//==============================================================================

TEST_F(SymplecticIntegratorTest, SingleStepNullCondition) {
    // Tolerance for single Step with fully analytic dHdq
    // Measured: ~3.93e-5 per Step
    const double tolerance = 1e-4;

    GeodesicStateD initial = createIngoingNullRay(schwarzschild.get(), 10.0, std::numbers::pi / 2);
    HamiltonianStateD hs(initial);

    // Initial Hamiltonian should be ~0
    double H0 = schwarzschild->Hamiltonian(hs.q, hs.p);
    EXPECT_NEAR(H0, 0.0, 1e-10) << "Initial state not null";

    // Take one Step
    auto result = integrator_sch->Step(hs, 0.1);

    // Hamiltonian should still be ~0
    EXPECT_NEAR(result.hamiltonian_error, 0.0, tolerance) << "Single Step broke null condition";
}

//==============================================================================
// Test: Integrator Order Comparison
// Verifies: Higher-order compositions reduce phase-state truncation error
//==============================================================================

TEST_F(SymplecticIntegratorTest, IntegratorOrderComparison) {
    const HamiltonianStateD initial(createIngoingNullRay(kerr.get(), 8.0, 1.1, 2.5));
    constexpr double affine_span = 0.5;

    AdaptiveHamiltonianIntegratorD::Config reference_config;
    reference_config.absolute_tolerance = 1.0e-14;
    reference_config.relative_tolerance = 1.0e-14;
    reference_config.maximum_step = 1.0e-3;
    AdaptiveHamiltonianIntegratorD reference_integrator(kerr.get(), reference_config);
    const auto reference = reference_integrator.Integrate(initial, affine_span);
    ASSERT_TRUE(reference.complete);

    SymplecticIntegratorD::Config config4;
    config4.order = IntegratorOrder::kYoshida4;
    config4.enforce_null_condition = false;
    SymplecticIntegratorD integrator4(kerr.get(), config4);
    const HamiltonianStateD result4 = integrator4.Step(initial, affine_span).state;

    SymplecticIntegratorD::Config config6;
    config6.order = IntegratorOrder::kYoshida6;
    config6.enforce_null_condition = false;
    SymplecticIntegratorD integrator6(kerr.get(), config6);
    const HamiltonianStateD result6 = integrator6.Step(initial, affine_span).state;

    SymplecticIntegratorD::Config config8;
    config8.order = IntegratorOrder::kYoshida8;
    config8.enforce_null_condition = false;
    SymplecticIntegratorD integrator8(kerr.get(), config8);
    const HamiltonianStateD result8 = integrator8.Step(initial, affine_span).state;

    const auto phase_error = [&reference](const HamiltonianStateD& result) {
        double error = 0.0;
        for (int component = 0; component < 4; ++component) {
            error = std::max(error, std::abs(result.q[component] - reference.state.q[component]));
            error = std::max(error, std::abs(result.p[component] - reference.state.p[component]));
        }
        return error;
    };
    const double error4 = phase_error(result4);
    const double error6 = phase_error(result6);
    const double error8 = phase_error(result8);

    std::cout << "Yoshida phase-state error: order4=" << error4 << " order6=" << error6
              << " order8=" << error8 << '\n';
    EXPECT_LT(error6, error4);
    EXPECT_LT(error8, error6);
}

}  // namespace
