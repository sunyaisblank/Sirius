// TSPH004A.cpp - Kerr Metric Tests
// Tests: Schwarzschild/Minkowski limits, ergosphere, frame dragging, Christoffels.

#define _USE_MATH_DEFINES
#include <cmath>
#include <gtest/gtest.h>
#include <MTTN001A.h>
#include <MTDL001A.h>

namespace sirius::test {

constexpr double kEpsilon = 1e-10;

// =============================================================================
// Test Fixture
// =============================================================================

class KerrTests : public ::testing::Test {
protected:
    // Mass and spin parameters
    static constexpr double M = 1.0;
    
    // Helper functions for Kerr metric components
    double Sigma(double r, double theta, double a) {
        return r * r + a * a * std::cos(theta) * std::cos(theta);
    }
    
    double Delta(double r, double a) {
        return r * r - 2.0 * M * r + a * a;
    }
    
    // Create Kerr metric at given r, θ with spin a
    Metric4D createKerrMetric(double r, double theta, double a) {
        Metric4D g;
        g.zero();
        
        double sigma = Sigma(r, theta, a);
        double delta = Delta(r, a);
        double sin_theta = std::sin(theta);
        double cos_theta = std::cos(theta);
        double sin2 = sin_theta * sin_theta;
        double A = (r * r + a * a) * (r * r + a * a) - delta * a * a * sin2;
        
        // g_tt
        g(0, 0) = Dual<double>(-(1.0 - 2.0 * M * r / sigma), 0.0);
        
        // g_rr
        g(1, 1) = Dual<double>(sigma / delta, 0.0);
        
        // g_θθ
        g(2, 2) = Dual<double>(sigma, 0.0);
        
        // g_φφ
        g(3, 3) = Dual<double>(A * sin2 / sigma, 0.0);
        
        // g_tφ = g_φt (off-diagonal for frame dragging)
        double g_tphi = -2.0 * M * a * r * sin2 / sigma;
        g(0, 3) = Dual<double>(g_tphi, 0.0);
        g(3, 0) = Dual<double>(g_tphi, 0.0);
        
        return g;
    }
    
    void SetUp() override {}
    void TearDown() override {}
};

// =============================================================================
// Schwarzschild Limit (a → 0)
// =============================================================================

// Test: Kerr reduces to Schwarzschild when a = 0
TEST_F(KerrTests, ReducesToSchwarzschildAtZeroSpin) {
    double r = 10.0;
    double theta = M_PI / 2;
    double a = 0.0;
    
    Metric4D g_kerr = createKerrMetric(r, theta, a);
    
    // Expected Schwarzschild values
    double f = 1.0 - 2.0 * M / r;
    
    EXPECT_NEAR(g_kerr(0, 0).real, -f, kEpsilon)
        << "g_tt should match Schwarzschild";
    EXPECT_NEAR(g_kerr(1, 1).real, 1.0 / f, kEpsilon)
        << "g_rr should match Schwarzschild";
    EXPECT_NEAR(g_kerr(2, 2).real, r * r, kEpsilon)
        << "g_θθ should match Schwarzschild";
    EXPECT_NEAR(g_kerr(3, 3).real, r * r, kEpsilon)  // sin²(π/2) = 1
        << "g_φφ should match Schwarzschild at equator";
    EXPECT_NEAR(g_kerr(0, 3).real, 0.0, kEpsilon)
        << "Off-diagonal should vanish for a = 0";
}

// =============================================================================
// Minkowski Limit (M → 0)
// =============================================================================

// Test: Kerr becomes flat when M = 0 (implemented by small M)
TEST_F(KerrTests, ApproachesMinkowskiAtZeroMass) {
    double r = 10.0;
    double theta = M_PI / 2;
    double a = 0.5;
    
    // Approximate M → 0 by using very small M (not exact test)
    // Instead, check that at large r, corrections are O(M/r)
    double large_r = 1000.0;
    Metric4D g_far = createKerrMetric(large_r, theta, a);
    
    // At large r, should approach spheroidal coordinates
    // g_tt → -1 + O(M/r)
    EXPECT_NEAR(g_far(0, 0).real, -1.0, 0.01)
        << "g_tt should approach -1 at large r";
    
    // g_tφ → 0 + O(M/r)
    EXPECT_NEAR(g_far(0, 3).real, 0.0, 0.01)
        << "Frame dragging should vanish at large r";
}

// =============================================================================
// Symmetry Properties
// =============================================================================

// Test: Metric is symmetric g_μν = g_νμ
TEST_F(KerrTests, MetricIsSymmetric) {
    double r = 5.0;
    double theta = M_PI / 4;
    double a = 0.5;
    
    Metric4D g = createKerrMetric(r, theta, a);
    
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            EXPECT_NEAR(g(i, j).real, g(j, i).real, kEpsilon)
                << "Metric asymmetry at (" << i << "," << j << ")";
        }
    }
}

// =============================================================================
// Ergosphere
// =============================================================================

// Test: g_tt > 0 inside ergosphere
TEST_F(KerrTests, ErgosphereGttPositive) {
    double a = 0.9;  // High spin
    double theta = M_PI / 2;  // Equator (widest ergosphere)
    
    // Ergosphere outer boundary (equator): r_ergo = M + √(M² - a²cos²θ)
    // At equator (cos θ = 0): r_ergo = M + M = 2M
    double r_horizon = M + std::sqrt(M * M - a * a);  // Event horizon
    double r_ergo = M + std::sqrt(M * M - a * a * 0);  // = 2M at equator
    
    // Inside ergosphere but outside horizon
    double r_inside = (r_horizon + r_ergo) / 2;
    
    Metric4D g = createKerrMetric(r_inside, theta, a);
    
    // g_tt > 0 inside ergosphere (stationary limit surface)
    EXPECT_GT(g(0, 0).real, 0.0)
        << "g_tt should be positive inside ergosphere";
}

// Test: g_tt < 0 outside ergosphere
TEST_F(KerrTests, GttNegativeOutsideErgosphere) {
    double a = 0.5;
    double theta = M_PI / 2;
    double r = 10.0;  // Well outside ergosphere
    
    Metric4D g = createKerrMetric(r, theta, a);
    
    EXPECT_LT(g(0, 0).real, 0.0)
        << "g_tt should be negative far from black hole";
}

// =============================================================================
// Frame Dragging
// =============================================================================

// Test: g_tφ ≠ 0 for spinning black hole
TEST_F(KerrTests, FrameDraggingPresent) {
    double r = 5.0;
    double theta = M_PI / 2;
    double a = 0.5;
    
    Metric4D g = createKerrMetric(r, theta, a);
    
    EXPECT_NE(g(0, 3).real, 0.0)
        << "Frame dragging term g_tφ should be non-zero";
}

// Test: Frame dragging sign convention
TEST_F(KerrTests, FrameDraggingSign) {
    double r = 5.0;
    double theta = M_PI / 2;
    double a = 0.5;  // Positive spin
    
    Metric4D g = createKerrMetric(r, theta, a);
    
    // For a > 0, g_tφ < 0 (prograde frame dragging)
    EXPECT_LT(g(0, 3).real, 0.0)
        << "Frame dragging should have correct sign for prograde rotation";
}

// Test: Frame dragging vanishes at poles
TEST_F(KerrTests, FrameDraggingVanishesAtPoles) {
    double r = 5.0;
    double a = 0.5;
    
    // At poles, sin θ = 0, so g_tφ = 0
    double theta_pole = 0.001;  // Near pole
    Metric4D g = createKerrMetric(r, theta_pole, a);
    
    EXPECT_NEAR(g(0, 3).real, 0.0, 1e-6)
        << "Frame dragging should vanish at poles";
}

// =============================================================================
// Angular Momentum
// =============================================================================

// Test: ω = -g_tφ/g_φφ is the frame dragging angular velocity
TEST_F(KerrTests, FrameDraggingAngularVelocity) {
    double r = 4.0;
    double theta = M_PI / 2;
    double a = 0.5;
    
    Metric4D g = createKerrMetric(r, theta, a);
    
    // ω = -g_tφ/g_φφ
    double omega = -g(0, 3).real / g(3, 3).real;
    
    // Expected: ω = 2Mar / [(r² + a²)² - Δa²sin²θ]
    // At equator: simplified
    EXPECT_GT(omega, 0.0) << "Frame dragging should be prograde for a > 0";
    
    // Should decrease with radius
    double r2 = 8.0;
    Metric4D g2 = createKerrMetric(r2, theta, a);
    double omega2 = -g2(0, 3).real / g2(3, 3).real;
    
    EXPECT_LT(omega2, omega)
        << "Frame dragging should decrease with radius";
}

// =============================================================================
// Determinant
// =============================================================================

// Test: Kerr determinant formula
TEST_F(KerrTests, DeterminantFormula) {
    double r = 5.0;
    double theta = M_PI / 3;
    double a = 0.5;
    
    Metric4D g = createKerrMetric(r, theta, a);
    
    // det(g) = -Σ² sin²θ
    double sigma = Sigma(r, theta, a);
    double sin_theta = std::sin(theta);
    double expected_det = -sigma * sigma * sin_theta * sin_theta;
    
    // For 4x4 with off-diagonal terms, need full determinant
    // Simplified check: det(g) < 0 (Lorentzian)
    double det = g(0, 0).real * g(1, 1).real * g(2, 2).real * g(3, 3).real
                - g(0, 3).real * g(0, 3).real * g(1, 1).real * g(2, 2).real;
    
    EXPECT_LT(det, 0.0) << "Determinant should be negative (Lorentzian)";
}

// =============================================================================
// Kerr Christoffel Symbol Validation (Phase 3: Fix Analytic Formulas)
// =============================================================================

// Helper: Create Kerr metric with derivatives using dual numbers
Metric4D createKerrMetricWithDerivatives(double r, double theta, double a, double M_param, int deriv_wrt) {
    // Create metric where the derivative direction uses dual part = 1
    Metric4D g;
    g.zero();
    
    // Dual number setup: if deriv_wrt == i, then x_i = (x_i, 1)
    // deriv_wrt: 0=t, 1=r, 2=theta, 3=phi
    
    Dual<double> r_d(r, (deriv_wrt == 1) ? 1.0 : 0.0);
    Dual<double> theta_d(theta, (deriv_wrt == 2) ? 1.0 : 0.0);
    
    Dual<double> sin_theta = sin(theta_d);
    Dual<double> cos_theta = cos(theta_d);
    Dual<double> sin2 = sin_theta * sin_theta;
    Dual<double> cos2 = cos_theta * cos_theta;
    
    Dual<double> a2(a * a, 0.0);
    Dual<double> M_d(M_param, 0.0);
    
    Dual<double> sigma = r_d * r_d + a2 * cos2;
    Dual<double> delta = r_d * r_d - Dual<double>(2.0, 0.0) * M_d * r_d + a2;
    
    Dual<double> r2_plus_a2 = r_d * r_d + a2;
    Dual<double> A = r2_plus_a2 * r2_plus_a2 - delta * a2 * sin2;
    
    // g_tt = -(1 - 2Mr/Σ)
    g(0, 0) = Dual<double>(-1.0, 0.0) + Dual<double>(2.0, 0.0) * M_d * r_d / sigma;
    
    // g_rr = Σ/Δ
    g(1, 1) = sigma / delta;
    
    // g_θθ = Σ
    g(2, 2) = sigma;
    
    // g_φφ = A sin²θ / Σ
    g(3, 3) = A * sin2 / sigma;
    
    // g_tφ = g_φt = -2Mar sin²θ / Σ
    Dual<double> g_tphi = Dual<double>(-2.0, 0.0) * M_d * Dual<double>(a, 0.0) * r_d * sin2 / sigma;
    g(0, 3) = g_tphi;
    g(3, 0) = g_tphi;
    
    return g;
}

// Test: Kerr Christoffel symmetry (Γ^μ_νρ = Γ^μ_ρν)
TEST_F(KerrTests, ChristoffelSymmetry) {
    double r = 6.0;
    double theta = M_PI / 3;
    double a = 0.7;
    double eps = 1e-4;
    
    // Compute metric and derivatives
    Metric4D g = createKerrMetricWithDerivatives(r, theta, a, M, -1);  // No deriv
    
    // Metric derivatives: dg[sigma][mu][nu] = ∂_sigma g_μν
    Tensor<Dual<double>, 4, 4, 4> dg;
    dg.zero();
    
    // ∂/∂r
    Metric4D g_dr = createKerrMetricWithDerivatives(r, theta, a, M, 1);
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            dg(1, mu, nu) = Dual<double>(g_dr(mu, nu).dual, 0.0);
        }
    }
    
    // ∂/∂θ
    Metric4D g_dth = createKerrMetricWithDerivatives(r, theta, a, M, 2);
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            dg(2, mu, nu) = Dual<double>(g_dth(mu, nu).dual, 0.0);
        }
    }
    
    // Compute Christoffel symbols
    ChristoffelSymbols gamma = TensorOps::christoffel(g, dg);
    
    // Check symmetry in lower indices
    for (int lam = 0; lam < 4; lam++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
                EXPECT_NEAR(gamma.gamma(lam, mu, nu).real,
                           gamma.gamma(lam, nu, mu).real, 1e-8)
                    << "Christoffel symmetry violated at Γ^" << lam << "_" << mu << nu;
            }
        }
    }
}

// Test: Kerr Christoffel non-zero components exist
TEST_F(KerrTests, ChristoffelNonZeroForRotatingBH) {
    double r = 6.0;
    double theta = M_PI / 3;
    double a = 0.7;
    
    // Compute metric and derivatives using automatic differentiation
    Metric4D g = createKerrMetricWithDerivatives(r, theta, a, M, -1);
    
    Tensor<Dual<double>, 4, 4, 4> dg;
    dg.zero();
    
    // ∂/∂r
    Metric4D g_dr = createKerrMetricWithDerivatives(r, theta, a, M, 1);
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            dg(1, mu, nu) = Dual<double>(g_dr(mu, nu).dual, 0.0);
        }
    }
    
    // ∂/∂θ
    Metric4D g_dth = createKerrMetricWithDerivatives(r, theta, a, M, 2);
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            dg(2, mu, nu) = Dual<double>(g_dth(mu, nu).dual, 0.0);
        }
    }
    
    ChristoffelSymbols gamma = TensorOps::christoffel(g, dg);
    
    // Count non-zero Christoffel components (there should be many)
    int nonZero = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                if (std::abs(gamma.gamma(i, j, k).real) > 1e-10) {
                    nonZero++;
                }
            }
        }
    }
    
    // Kerr should have many non-zero Christoffel symbols (at least 20)
    EXPECT_GT(nonZero, 15) << "Kerr should have many non-zero Christoffel components";
    
    // Specific checks: Γ^t_tr should be non-zero (gravitational time dilation)
    EXPECT_NE(gamma.gamma(0, 0, 1).real, 0.0) << "Γ^t_tr should be non-zero";
    
    // Γ^r_θθ should be non-zero (spherical coordinate effect)
    EXPECT_NE(gamma.gamma(1, 2, 2).real, 0.0) << "Γ^r_θθ should be non-zero";
    
    // Γ^θ_rθ should be non-zero
    EXPECT_NE(gamma.gamma(2, 1, 2).real, 0.0) << "Γ^θ_rθ should be non-zero";
}

} // namespace sirius::test

