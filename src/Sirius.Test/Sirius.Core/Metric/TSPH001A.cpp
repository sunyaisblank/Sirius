// TSPH001A.cpp - Schwarzschild Metric Tests
// ds² = -(1-2M/r)dt² + dr²/(1-2M/r) + r²dΩ²
// Tests: signature, asymptotic flatness, horizon, proper time, redshift.

#define _USE_MATH_DEFINES
#include <cmath>
#include <gtest/gtest.h>
#include <MTTN001A.h>
#include <MTDL001A.h>

namespace sirius::test {
using namespace Sirius;

constexpr double kEpsilon = 1e-10;

// =============================================================================
// Test Fixture
// =============================================================================

class SchwarzschildTests : public ::testing::Test {
protected:
    // Mass parameter (in geometric units where G = c = 1)
    static constexpr double M = 1.0;
    
    // Schwarzschild radius
    static constexpr double rs = 2.0 * M;
    
    // Create Schwarzschild metric at given r, θ
    Metric4D createSchwarzschildMetric(double r, double theta = M_PI/2) {
        Metric4D g;
        g.zero();
        
        double f = 1.0 - rs / r;
        double sin_theta = std::sin(theta);
        
        g(0, 0) = Dual<double>(-f, 0.0);              // g_tt
        g(1, 1) = Dual<double>(1.0 / f, 0.0);         // g_rr
        g(2, 2) = Dual<double>(r * r, 0.0);           // g_θθ
        g(3, 3) = Dual<double>(r * r * sin_theta * sin_theta, 0.0);  // g_φφ
        
        return g;
    }
    
    void SetUp() override {}
    void TearDown() override {}
};

// =============================================================================
// Metric Structure Tests
// =============================================================================

// Test: Schwarzschild metric is diagonal
TEST_F(SchwarzschildTests, MetricIsDiagonal) {
    double r = 10.0;
    Metric4D g = createSchwarzschildMetric(r);
    
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i != j) {
                EXPECT_NEAR(g(i, j).real, 0.0, kEpsilon)
                    << "Off-diagonal non-zero at (" << i << "," << j << ")";
            }
        }
    }
}

// Test: Schwarzschild metric is symmetric
TEST_F(SchwarzschildTests, MetricIsSymmetric) {
    double r = 5.0;
    Metric4D g = createSchwarzschildMetric(r);
    
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            EXPECT_NEAR(g(i, j).real, g(j, i).real, kEpsilon)
                << "Asymmetry at (" << i << "," << j << ")";
        }
    }
}

// Test: Lorentzian signature (-,+,+,+)
TEST_F(SchwarzschildTests, MetricHasLorentzianSignature) {
    double r = 10.0;
    Metric4D g = createSchwarzschildMetric(r);
    
    EXPECT_LT(g(0, 0).real, 0.0) << "g_tt should be negative outside horizon";
    EXPECT_GT(g(1, 1).real, 0.0) << "g_rr should be positive outside horizon";
    EXPECT_GT(g(2, 2).real, 0.0) << "g_θθ should be positive";
    EXPECT_GT(g(3, 3).real, 0.0) << "g_φφ should be positive";
}

// =============================================================================
// Asymptotic Flatness
// =============================================================================

// Test: Metric approaches Minkowski as r → ∞
TEST_F(SchwarzschildTests, AsymptoticFlatness) {
    std::vector<double> radii = {100.0, 1000.0, 10000.0};
    
    for (double r : radii) {
        Metric4D g = createSchwarzschildMetric(r);
        
        double expected_gtt = -(1.0 - rs/r);
        double expected_grr = 1.0 / (1.0 - rs/r);
        
        EXPECT_NEAR(g(0, 0).real, expected_gtt, kEpsilon)
            << "g_tt incorrect at r=" << r;
        EXPECT_NEAR(g(1, 1).real, expected_grr, kEpsilon)
            << "g_rr incorrect at r=" << r;
        
        // As r → ∞, should approach Minkowski
        double deviation_tt = std::abs(g(0, 0).real - (-1.0));
        double deviation_rr = std::abs(g(1, 1).real - 1.0);
        
        EXPECT_LT(deviation_tt, rs/r + kEpsilon)
            << "g_tt deviates too much from -1 at large r=" << r;
        EXPECT_LT(deviation_rr, 2*rs/r + kEpsilon)  // First-order expansion
            << "g_rr deviates too much from 1 at large r=" << r;
    }
}

// =============================================================================
// Schwarzschild Radius (Event Horizon)
// =============================================================================

// Test: g_tt → 0 as r → r_s (coordinate singularity)
TEST_F(SchwarzschildTests, EventHorizonGtt) {
    // Approach from outside
    std::vector<double> radii = {2.1, 2.05, 2.01, 2.001};
    
    for (double r : radii) {
        Metric4D g = createSchwarzschildMetric(r);
        
        double expected_gtt = -(1.0 - rs/r);
        EXPECT_NEAR(g(0, 0).real, expected_gtt, kEpsilon);
        
        // g_tt should approach 0 as r → rs
        EXPECT_GT(g(0, 0).real, -1.0) << "g_tt should approach 0 as r → rs";
    }
}

// Test: g_rr → ∞ as r → r_s
TEST_F(SchwarzschildTests, EventHorizonGrr) {
    double r = 2.001;  // Just outside horizon
    Metric4D g = createSchwarzschildMetric(r);
    
    // g_rr should be large near horizon
    EXPECT_GT(g(1, 1).real, 1000.0)
        << "g_rr should blow up near horizon";
}

// =============================================================================
// Angular Components
// =============================================================================

// Test: g_θθ = r² (exact)
TEST_F(SchwarzschildTests, MetricThetaTheta) {
    std::vector<double> radii = {3.0, 5.0, 10.0, 100.0};
    
    for (double r : radii) {
        Metric4D g = createSchwarzschildMetric(r);
        
        EXPECT_NEAR(g(2, 2).real, r * r, kEpsilon)
            << "g_θθ should equal r² at r=" << r;
    }
}

// Test: g_φφ = r² sin²θ
TEST_F(SchwarzschildTests, MetricPhiPhi) {
    double r = 10.0;
    std::vector<double> thetas = {M_PI/6, M_PI/4, M_PI/3, M_PI/2};
    
    for (double theta : thetas) {
        Metric4D g = createSchwarzschildMetric(r, theta);
        double sin_theta = std::sin(theta);
        
        EXPECT_NEAR(g(3, 3).real, r * r * sin_theta * sin_theta, kEpsilon)
            << "g_φφ incorrect at θ=" << theta;
    }
}

// Test: g_φφ → 0 at poles (θ = 0, π)
TEST_F(SchwarzschildTests, MetricPhiPhiAtPoles) {
    double r = 10.0;
    std::vector<double> thetas = {0.001, M_PI - 0.001};
    
    for (double theta : thetas) {
        Metric4D g = createSchwarzschildMetric(r, theta);
        
        // g_φφ should be small near poles
        EXPECT_LT(g(3, 3).real, r * r * 0.001 * 0.001 + kEpsilon)
            << "g_φφ should vanish at poles";
    }
}

// =============================================================================
// Physical Invariants
// =============================================================================

// Test: Proper time for stationary observer
// dτ² = -g_tt dt² = (1 - 2M/r) dt²
TEST_F(SchwarzschildTests, ProperTimeStationaryObserver) {
    double r = 10.0;
    double dt = 1.0;  // Coordinate time interval
    
    Metric4D g = createSchwarzschildMetric(r);
    
    // Proper time: dτ² = -g_tt dt²
    double dtau_squared = -g(0, 0).real * dt * dt;
    double dtau = std::sqrt(dtau_squared);
    
    // Expected: dτ = √(1 - rs/r) dt
    double expected_dtau = std::sqrt(1.0 - rs/r) * dt;
    
    EXPECT_NEAR(dtau, expected_dtau, kEpsilon)
        << "Proper time gravitational redshift incorrect";
}

// Test: Gravitational redshift factor
// z = √(g_tt(∞)/g_tt(r)) - 1 = 1/√(1 - rs/r) - 1
TEST_F(SchwarzschildTests, GravitationalRedshift) {
    double r = 3.0;  // Just outside photon sphere
    
    double f = 1.0 - rs / r;
    double redshift = 1.0 / std::sqrt(f) - 1.0;
    
    // At r = 3M (rs = 2M), f = 1/3, so z = √3 - 1 ≈ 0.732
    double expected_redshift = std::sqrt(3.0) - 1.0;
    
    EXPECT_NEAR(redshift, expected_redshift, kEpsilon)
        << "Gravitational redshift incorrect at r = 3M";
}

// =============================================================================
// Determinant
// =============================================================================

// Test: det(g) = -r⁴ sin²θ
TEST_F(SchwarzschildTests, MetricDeterminant) {
    double r = 10.0;
    double theta = M_PI / 3;
    Metric4D g = createSchwarzschildMetric(r, theta);
    
    // For diagonal metric: det(g) = g_tt * g_rr * g_θθ * g_φφ
    double det = g(0, 0).real * g(1, 1).real * g(2, 2).real * g(3, 3).real;
    
    // Expected: -r⁴ sin²θ (the f factors cancel: -f * (1/f) = -1)
    double sin_theta = std::sin(theta);
    double expected_det = -r * r * r * r * sin_theta * sin_theta;
    
    EXPECT_NEAR(det, expected_det, kEpsilon)
        << "Metric determinant incorrect";
}

// Test: Determinant is negative (Lorentzian)
TEST_F(SchwarzschildTests, DeterminantIsNegative) {
    std::vector<double> radii = {3.0, 5.0, 10.0, 100.0};
    
    for (double r : radii) {
        Metric4D g = createSchwarzschildMetric(r, M_PI/2);
        double det = g(0, 0).real * g(1, 1).real * g(2, 2).real * g(3, 3).real;
        
        EXPECT_LT(det, 0.0) << "Determinant should be negative at r=" << r;
    }
}

} // namespace sirius::test
