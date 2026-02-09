// TSPH010A.cpp - Numerical Geodesic Integrator Tests (PHGD002A)
// Tests: RK4 for numerical spacetimes, termination, coordinate handling.

#include <gtest/gtest.h>
#include <cmath>
#include "PHGD002A.h"
#include "PHMT000A.h"
#include "PHMT100A.h" // Unified Kerr-Schild family
#include "MTTN001A.h"

namespace sirius::test {
using namespace Sirius;

using namespace Sirius::Physics;

// =============================================================================
// Constants
// =============================================================================

constexpr double PI = 3.14159265358979323846;
constexpr double M = 1.0;  // Mass in geometric units

// =============================================================================
// Test Fixture
// =============================================================================

class NumericalGeodesicTests : public ::testing::Test {
protected:
    Sirius::KerrSchildFamily minkowski{Sirius::KerrSchildParams::Minkowski()};
    Sirius::KerrSchildFamily schwarzschild{Sirius::KerrSchildParams::Schwarzschild(M)};
    NumericalGeodesicIntegrator integrator;
    
    void SetUp() override {
        integrator.setStepSize(0.01);
        integrator.setMaxSteps(1000);
        integrator.setMaxDistance(100.0);
        integrator.setHorizonRadius(2.0 * M);
    }
    
    void TearDown() override {}
    
    // Create initial state at given position moving outward
    GeodesicState createOutwardState(double r, double theta = PI/2) {
        GeodesicState state;
        state.x(0) = 0.0;      // t
        state.x(1) = r;        // r
        state.x(2) = theta;    // theta
        state.x(3) = 0.0;      // phi
        
        state.u(0) = 1.0;      // dt/dλ
        state.u(1) = 0.5;      // dr/dλ (outward)
        state.u(2) = 0.0;      // dθ/dλ
        state.u(3) = 0.3;      // dφ/dλ
        
        return state;
    }
    
    // Create initial state moving inward
    GeodesicState createInwardState(double r, double theta = PI/2) {
        GeodesicState state;
        state.x(0) = 0.0;
        state.x(1) = r;
        state.x(2) = theta;
        state.x(3) = 0.0;
        
        state.u(0) = 1.0;
        state.u(1) = -0.5;   // dr/dλ (inward)
        state.u(2) = 0.0;
        state.u(3) = 0.3;
        
        return state;
    }
};

// =============================================================================
// Constructor and Parameter Tests
// =============================================================================

// Test: Default construction
TEST_F(NumericalGeodesicTests, DefaultConstruction) {
    NumericalGeodesicIntegrator integ;
    // Should not crash
    EXPECT_TRUE(true);
}

// Test: Parameter setters work
TEST_F(NumericalGeodesicTests, ParameterSetters) {
    NumericalGeodesicIntegrator integ;
    
    // These should not throw
    EXPECT_NO_THROW(integ.setMaxSteps(500));
    EXPECT_NO_THROW(integ.setStepSize(0.05));
    EXPECT_NO_THROW(integ.setTolerance(1e-8));
    EXPECT_NO_THROW(integ.setMaxDistance(200.0));
    EXPECT_NO_THROW(integ.setHorizonRadius(3.0));
    EXPECT_NO_THROW(integ.setMetric(&minkowski));
}

// Test: Integration without metric returns error
TEST_F(NumericalGeodesicTests, IntegrateWithoutMetric) {
    NumericalGeodesicIntegrator integ;
    // Don't set metric
    
    GeodesicState initial = createOutwardState(10.0);
    IntegrationResult result = integ.integrate(initial);
    
    EXPECT_TRUE(result.error) << "Should return error when no metric is set";
}

// =============================================================================
// Derivative Computation Tests
// =============================================================================

// Test: computeDerivative returns finite values
TEST_F(NumericalGeodesicTests, DerivativeFinite) {
    integrator.setMetric(&schwarzschild);
    
    GeodesicState state = createOutwardState(10.0);
    GeodesicState deriv = integrator.computeDerivative(state);
    
    for (int i = 0; i < 4; ++i) {
        EXPECT_FALSE(std::isnan(deriv.x(i))) << "deriv.x(" << i << ") is NaN";
        EXPECT_FALSE(std::isinf(deriv.x(i))) << "deriv.x(" << i << ") is Inf";
        EXPECT_FALSE(std::isnan(deriv.u(i))) << "deriv.u(" << i << ") is NaN";
        EXPECT_FALSE(std::isinf(deriv.u(i))) << "deriv.u(" << i << ") is Inf";
    }
}

// Test: dx/dλ = u (position derivative equals velocity)
TEST_F(NumericalGeodesicTests, DerivativePositionEqualsVelocity) {
    integrator.setMetric(&schwarzschild);
    
    GeodesicState state = createOutwardState(10.0);
    GeodesicState deriv = integrator.computeDerivative(state);
    
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(deriv.x(i), state.u(i), 1e-10)
            << "dx^" << i << "/dλ should equal u^" << i;
    }
}

// Test: Minkowski derivative has zero acceleration
TEST_F(NumericalGeodesicTests, MinkowskiZeroAcceleration) {
    integrator.setMetric(&minkowski);
    
    GeodesicState state = createOutwardState(10.0);
    GeodesicState deriv = integrator.computeDerivative(state);
    
    // In Minkowski space, Christoffel symbols are zero, so du/dλ = 0
    double accel_norm = 0.0;
    for (int i = 0; i < 4; ++i) {
        accel_norm += deriv.u(i) * deriv.u(i);
    }
    accel_norm = std::sqrt(accel_norm);
    
    EXPECT_LT(accel_norm, 1e-6) << "Minkowski should have near-zero acceleration";
}

// =============================================================================
// RK4 Step Tests
// =============================================================================

// Test: rk4Step returns valid state
TEST_F(NumericalGeodesicTests, RK4StepValid) {
    integrator.setMetric(&schwarzschild);
    
    GeodesicState state = createOutwardState(10.0);
    GeodesicState next = integrator.rk4Step(state, 0.01);
    
    for (int i = 0; i < 4; ++i) {
        EXPECT_FALSE(std::isnan(next.x(i))) << "next.x(" << i << ") is NaN";
        EXPECT_FALSE(std::isinf(next.x(i))) << "next.x(" << i << ") is Inf";
    }
}

// Test: rk4Step moves position
TEST_F(NumericalGeodesicTests, RK4StepMovesPosition) {
    integrator.setMetric(&schwarzschild);
    
    GeodesicState state = createOutwardState(10.0);
    GeodesicState next = integrator.rk4Step(state, 0.1);
    
    // Position should have changed
    double pos_change = 0.0;
    for (int i = 0; i < 4; ++i) {
        double diff = next.x(i) - state.x(i);
        pos_change += diff * diff;
    }
    pos_change = std::sqrt(pos_change);
    
    EXPECT_GT(pos_change, 0.0) << "Position should change after RK4 step";
}

// Test: Minkowski straight line
TEST_F(NumericalGeodesicTests, MinkowskiStraightLine) {
    integrator.setMetric(&minkowski);
    
    GeodesicState state = createOutwardState(10.0);
    double initial_r = state.x(1);
    double dr_dlambda = state.u(1);
    
    // Take 10 steps
    double h = 0.1;
    for (int i = 0; i < 10; ++i) {
        state = integrator.rk4Step(state, h);
    }
    
    // In flat space, r should increase linearly
    double expected_r = initial_r + dr_dlambda * (10 * h);
    EXPECT_NEAR(state.x(1), expected_r, 0.01)
        << "r should increase linearly in flat space";
}

// =============================================================================
// Integration Tests
// =============================================================================

// Test: Integration terminates properly (escape or horizon)
TEST_F(NumericalGeodesicTests, IntegrationTerminates) {
    integrator.setMetric(&schwarzschild);
    integrator.setMaxDistance(50.0);
    integrator.setMaxSteps(500);
    
    GeodesicState initial = createOutwardState(10.0);
    IntegrationResult result = integrator.integrate(initial);
    
    // Should terminate via escape, horizon, or max steps
    EXPECT_TRUE(result.escaped || result.hitHorizon || result.steps >= 500)
        << "Integration should terminate";
    EXPECT_FALSE(result.error) << "Should not have numerical error";
    EXPECT_GT(result.steps, 0);
    EXPECT_GT(result.affineDist, 0.0);
    
    // Log what happened for debugging
    std::cout << "\n[NumericalGeodesic] Integration result:\n";
    std::cout << "  Steps: " << result.steps << "\n";
    std::cout << "  Affine distance: " << result.affineDist << "\n";
    std::cout << "  Escaped: " << result.escaped << "\n";
    std::cout << "  Hit horizon: " << result.hitHorizon << "\n";
    std::cout << "  Final r: " << result.finalState.x(1) << "\n";
}

// Test: Integration of inward ray terminates
TEST_F(NumericalGeodesicTests, InwardRayTerminates) {
    integrator.setMetric(&schwarzschild);
    integrator.setHorizonRadius(2.0);
    integrator.setMaxSteps(500);
    
    GeodesicState initial = createInwardState(5.0);  // Start close, moving inward
    IntegrationResult result = integrator.integrate(initial);
    
    // Should terminate via horizon, escape (if deflected), or max steps
    EXPECT_TRUE(result.escaped || result.hitHorizon || result.steps >= 500)
        << "Should terminate eventually";
    EXPECT_GT(result.steps, 0);
    
    std::cout << "\n[NumericalGeodesic] Inward ray result:\n";
    std::cout << "  Steps: " << result.steps << "\n";
    std::cout << "  Escaped: " << result.escaped << "\n";
    std::cout << "  Hit horizon: " << result.hitHorizon << "\n";
    std::cout << "  Final r: " << result.finalState.x(1) << "\n";
}

// Test: Integration accumulates affine distance
TEST_F(NumericalGeodesicTests, IntegrationAffineDistance) {
    integrator.setMetric(&minkowski);
    integrator.setMaxDistance(50.0);
    
    GeodesicState initial = createOutwardState(10.0);
    IntegrationResult result = integrator.integrate(initial);
    
    EXPECT_GT(result.affineDist, 0.0) << "Affine distance should accumulate";
}

// Test: Step count varies with step size
TEST_F(NumericalGeodesicTests, StepCountVariesWithStepSize) {
    integrator.setMetric(&schwarzschild);
    integrator.setMaxDistance(30.0);
    
    GeodesicState initial = createOutwardState(10.0);
    
    // Small step size
    integrator.setStepSize(0.01);
    IntegrationResult result1 = integrator.integrate(initial);
    
    // Large step size
    integrator.setStepSize(0.1);
    IntegrationResult result2 = integrator.integrate(initial);
    
    // Smaller step size should require more steps
    EXPECT_GT(result1.steps, result2.steps)
        << "Smaller step size should require more steps";
}

// =============================================================================
// Termination Condition Tests
// =============================================================================

// Test: Max steps limit enforced
TEST_F(NumericalGeodesicTests, MaxStepsEnforced) {
    integrator.setMetric(&schwarzschild);
    integrator.setMaxSteps(10);
    integrator.setMaxDistance(1000.0);  // Very large distance
    
    // Create a ray that will take many steps to escape
    GeodesicState initial;
    initial.x(0) = 0; initial.x(1) = 5.0; initial.x(2) = PI/2; initial.x(3) = 0;
    initial.u(0) = 1; initial.u(1) = 0.01; initial.u(2) = 0; initial.u(3) = 0.5;
    
    IntegrationResult result = integrator.integrate(initial);
    
    EXPECT_LE(result.steps, 10) << "Should not exceed max steps";
}

// Test: Phi wrapping works
TEST_F(NumericalGeodesicTests, PhiWrapping) {
    integrator.setMetric(&schwarzschild);
    integrator.setMaxSteps(100);
    integrator.setMaxDistance(20.0);
    
    GeodesicState initial = createOutwardState(15.0);
    initial.u(3) = 2.0;  // Fast phi rotation
    
    IntegrationResult result = integrator.integrate(initial);
    
    // Final phi should be in [0, 2π)
    EXPECT_GE(result.finalState.x(3), 0.0) << "phi should be >= 0";
    EXPECT_LT(result.finalState.x(3), 2 * PI) << "phi should be < 2π";
}

// Test: Theta clamping works
TEST_F(NumericalGeodesicTests, ThetaClamping) {
    integrator.setMetric(&schwarzschild);
    integrator.setMaxSteps(100);
    
    // Start near pole with theta velocity
    GeodesicState initial;
    initial.x(0) = 0; initial.x(1) = 10.0; initial.x(2) = 0.1; initial.x(3) = 0;
    initial.u(0) = 1; initial.u(1) = 0.5; initial.u(2) = -0.5; initial.u(3) = 0;
    
    IntegrationResult result = integrator.integrate(initial);
    
    // Theta should be clamped
    EXPECT_GE(result.finalState.x(2), 0.01) << "theta should be clamped above 0.01";
    EXPECT_LE(result.finalState.x(2), PI - 0.01) << "theta should be clamped below π-0.01";
}

// =============================================================================
// Numerical Stability Tests
// =============================================================================

// Test: NaN detection works
TEST_F(NumericalGeodesicTests, NaNDetection) {
    integrator.setMetric(&schwarzschild);
    
    // Create invalid state with NaN
    GeodesicState invalid;
    invalid.x(0) = 0; invalid.x(1) = NAN; invalid.x(2) = PI/2; invalid.x(3) = 0;
    invalid.u(0) = 1; invalid.u(1) = 0; invalid.u(2) = 0; invalid.u(3) = 0;
    
    IntegrationResult result = integrator.integrate(invalid);
    
    EXPECT_TRUE(result.error) << "Should detect NaN and return error";
}

// Test: Multiple integrations are independent
TEST_F(NumericalGeodesicTests, IndependentIntegrations) {
    integrator.setMetric(&schwarzschild);
    integrator.setMaxDistance(30.0);
    
    GeodesicState initial = createOutwardState(10.0);
    
    IntegrationResult result1 = integrator.integrate(initial);
    IntegrationResult result2 = integrator.integrate(initial);
    
    // Same initial conditions should give same result
    EXPECT_EQ(result1.steps, result2.steps);
    EXPECT_DOUBLE_EQ(result1.affineDist, result2.affineDist);
    EXPECT_DOUBLE_EQ(result1.finalState.x(1), result2.finalState.x(1));
}

} // namespace sirius::test
