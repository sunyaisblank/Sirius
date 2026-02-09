// TSPH007A.cpp - Reissner-Nordström Metric Tests
// f(r)=1-2M/r+Q²/r². Tests: horizons, Schwarzschild limit, signature, asymptotic.

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

class ReissnerNordstromTests : public ::testing::Test {
protected:
    // Mass parameter (in geometric units where G = c = 1)
    static constexpr double M = 1.0;
    
    // Schwarzschild radius (for comparison)
    static constexpr double rs = 2.0 * M;
    
    // Create Reissner-Nordström metric at given r, θ with charge Q
    Metric4D createRNMetric(double r, double Q, double theta = M_PI/2) {
        Metric4D g;
        g.zero();
        
        // Q is given as fraction of M for sub-extremal condition
        double Q_actual = Q * M;  // Actual charge
        
        // f(r) = 1 - 2M/r + Q²/r²
        double f = 1.0 - 2.0 * M / r + (Q_actual * Q_actual) / (r * r);
        double sin_theta = std::sin(theta);
        
        g(0, 0) = Dual<double>(-f, 0.0);              // g_tt
        g(1, 1) = Dual<double>(1.0 / f, 0.0);         // g_rr
        g(2, 2) = Dual<double>(r * r, 0.0);           // g_θθ
        g(3, 3) = Dual<double>(r * r * sin_theta * sin_theta, 0.0);  // g_φφ
        
        return g;
    }
    
    // Calculate outer horizon radius r₊ = M + √(M² - Q²)
    double outerHorizon(double Q) const {
        double Q_actual = Q * M;
        double discriminant = M * M - Q_actual * Q_actual;
        return M + std::sqrt(std::max(discriminant, 0.0));
    }
    
    // Calculate inner horizon radius r₋ = M - √(M² - Q²)
    double innerHorizon(double Q) const {
        double Q_actual = Q * M;
        double discriminant = M * M - Q_actual * Q_actual;
        return M - std::sqrt(std::max(discriminant, 0.0));
    }
    
    // Calculate f(r) for a given Q
    double f_of_r(double r, double Q) const {
        double Q_actual = Q * M;
        return 1.0 - 2.0 * M / r + (Q_actual * Q_actual) / (r * r);
    }
    
    void SetUp() override {}
    void TearDown() override {}
};

// =============================================================================
// Metric Structure Tests
// =============================================================================

// Test: Reissner-Nordström metric is diagonal (like Schwarzschild)
TEST_F(ReissnerNordstromTests, MetricIsDiagonal) {
    double r = 10.0;
    double Q = 0.5;  // Sub-extremal
    Metric4D g = createRNMetric(r, Q);
    
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i != j) {
                EXPECT_NEAR(g(i, j).real, 0.0, kEpsilon)
                    << "Off-diagonal non-zero at (" << i << "," << j << ")";
            }
        }
    }
}

// Test: Reissner-Nordström metric is symmetric
TEST_F(ReissnerNordstromTests, MetricIsSymmetric) {
    double r = 5.0;
    double Q = 0.8;
    Metric4D g = createRNMetric(r, Q);
    
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            EXPECT_NEAR(g(i, j).real, g(j, i).real, kEpsilon)
                << "Asymmetry at (" << i << "," << j << ")";
        }
    }
}

// Test: Lorentzian signature (-,+,+,+) outside horizon
TEST_F(ReissnerNordstromTests, MetricHasLorentzianSignature) {
    double r = 10.0;  // Well outside horizon
    double Q = 0.5;
    Metric4D g = createRNMetric(r, Q);
    
    EXPECT_LT(g(0, 0).real, 0.0) << "g_tt should be negative outside horizon";
    EXPECT_GT(g(1, 1).real, 0.0) << "g_rr should be positive outside horizon";
    EXPECT_GT(g(2, 2).real, 0.0) << "g_θθ should be positive";
    EXPECT_GT(g(3, 3).real, 0.0) << "g_φφ should be positive";
}

// Test: Signature preserved for various charge values
TEST_F(ReissnerNordstromTests, SignatureForVariousCharges) {
    std::vector<double> charges = {0.0, 0.3, 0.5, 0.7, 0.9, 0.99};
    double r = 10.0;
    
    for (double Q : charges) {
        Metric4D g = createRNMetric(r, Q);
        
        EXPECT_LT(g(0, 0).real, 0.0) << "g_tt negative failed at Q=" << Q;
        EXPECT_GT(g(1, 1).real, 0.0) << "g_rr positive failed at Q=" << Q;
        EXPECT_GT(g(2, 2).real, 0.0) << "g_θθ positive failed at Q=" << Q;
        EXPECT_GT(g(3, 3).real, 0.0) << "g_φφ positive failed at Q=" << Q;
    }
}

// =============================================================================
// Reduction to Schwarzschild (Q → 0)
// =============================================================================

// Test: When Q=0, RN reduces to Schwarzschild exactly
TEST_F(ReissnerNordstromTests, ReducesToSchwarzschildZeroCharge) {
    double r = 10.0;
    double Q = 0.0;  // Zero charge
    Metric4D g_rn = createRNMetric(r, Q);
    
    // Expected Schwarzschild values
    double f_schw = 1.0 - rs / r;
    
    EXPECT_NEAR(g_rn(0, 0).real, -f_schw, kEpsilon)
        << "g_tt should match Schwarzschild at Q=0";
    EXPECT_NEAR(g_rn(1, 1).real, 1.0 / f_schw, kEpsilon)
        << "g_rr should match Schwarzschild at Q=0";
    EXPECT_NEAR(g_rn(2, 2).real, r * r, kEpsilon)
        << "g_θθ should match Schwarzschild at Q=0";
    
    double sin_theta = std::sin(M_PI / 2);
    EXPECT_NEAR(g_rn(3, 3).real, r * r * sin_theta * sin_theta, kEpsilon)
        << "g_φφ should match Schwarzschild at Q=0";
}

// Test: At multiple radii, Q=0 equals Schwarzschild
TEST_F(ReissnerNordstromTests, SchwarzschildLimitMultipleRadii) {
    std::vector<double> radii = {3.0, 5.0, 10.0, 50.0, 100.0};
    double Q = 0.0;
    
    for (double r : radii) {
        Metric4D g = createRNMetric(r, Q);
        double f_schw = 1.0 - rs / r;
        
        EXPECT_NEAR(g(0, 0).real, -f_schw, kEpsilon) << "at r=" << r;
        EXPECT_NEAR(g(1, 1).real, 1.0 / f_schw, kEpsilon) << "at r=" << r;
    }
}

// =============================================================================
// Asymptotic Flatness
// =============================================================================

// Test: Metric approaches Minkowski as r → ∞
TEST_F(ReissnerNordstromTests, AsymptoticFlatness) {
    std::vector<double> radii = {100.0, 1000.0, 10000.0};
    double Q = 0.5;
    
    for (double r : radii) {
        Metric4D g = createRNMetric(r, Q);
        
        // As r → ∞: g_tt → -1, g_rr → 1
        double deviation_tt = std::abs(g(0, 0).real - (-1.0));
        double deviation_rr = std::abs(g(1, 1).real - 1.0);
        
        // Leading order: deviation ~ 2M/r
        double expected_deviation = 3.0 * rs / r;  // Conservative bound
        
        EXPECT_LT(deviation_tt, expected_deviation)
            << "g_tt deviates too much from -1 at large r=" << r;
        EXPECT_LT(deviation_rr, expected_deviation)
            << "g_rr deviates too much from 1 at large r=" << r;
    }
}

// Test: Charge effects diminish faster than mass (Q²/r² vs 2M/r)
TEST_F(ReissnerNordstromTests, ChargeEffectsFallOffFaster) {
    double r = 100.0;
    double Q = 0.99;  // Near-extremal charge
    
    // At large r, f(r) ≈ 1 - 2M/r + Q²/r²
    // The Q² term falls off as 1/r², faster than 2M/r term
    double mass_correction = 2.0 * M / r;
    double Q_actual = Q * M;
    double charge_correction = (Q_actual * Q_actual) / (r * r);
    
    EXPECT_GT(mass_correction, charge_correction * 10)
        << "Mass correction should dominate at large r";
}

// =============================================================================
// Horizon Structure
// =============================================================================

// Test: Outer horizon at r₊ = M + √(M² - Q²)
TEST_F(ReissnerNordstromTests, OuterHorizonFormula) {
    std::vector<double> charges = {0.0, 0.3, 0.5, 0.7, 0.9};
    
    for (double Q : charges) {
        double r_plus = outerHorizon(Q);
        double Q_actual = Q * M;
        double expected = M + std::sqrt(M * M - Q_actual * Q_actual);
        
        EXPECT_NEAR(r_plus, expected, kEpsilon) << "at Q=" << Q;
        
        // Verify f(r₊) = 0
        double f_at_horizon = f_of_r(r_plus, Q);
        EXPECT_NEAR(f_at_horizon, 0.0, 1e-9) 
            << "f(r₊) should equal 0 at Q=" << Q;
    }
}

// Test: Inner horizon at r₋ = M - √(M² - Q²)
TEST_F(ReissnerNordstromTests, InnerHorizonFormula) {
    std::vector<double> charges = {0.3, 0.5, 0.7, 0.9};  // Not 0 (gives r₋=0)
    
    for (double Q : charges) {
        double r_minus = innerHorizon(Q);
        double Q_actual = Q * M;
        double expected = M - std::sqrt(M * M - Q_actual * Q_actual);
        
        EXPECT_NEAR(r_minus, expected, kEpsilon) << "at Q=" << Q;
        
        // r₋ should be positive for Q > 0
        EXPECT_GT(r_minus, 0.0) << "Inner horizon should be positive at Q=" << Q;
    }
}

// Test: Q=0 gives r₊ = 2M (Schwarzschild), r₋ = 0
TEST_F(ReissnerNordstromTests, SchwarzschildHorizon) {
    double Q = 0.0;
    double r_plus = outerHorizon(Q);
    double r_minus = innerHorizon(Q);
    
    EXPECT_NEAR(r_plus, rs, kEpsilon) << "r₊ should equal rs for Q=0";
    EXPECT_NEAR(r_minus, 0.0, kEpsilon) << "r₋ should equal 0 for Q=0";
}

// Test: Near-extremal case Q ≈ M, horizons approach each other
TEST_F(ReissnerNordstromTests, NearExtremalHorizons) {
    double Q = 0.999;  // Very close to extremal
    double r_plus = outerHorizon(Q);
    double r_minus = innerHorizon(Q);
    
    // Both should approach M as Q → M
    double gap = r_plus - r_minus;
    EXPECT_LT(gap, 0.1) << "Horizons should be close at near-extremal charge";
    
    // Average should be approximately M
    double average = (r_plus + r_minus) / 2.0;
    EXPECT_NEAR(average, M, 0.01) << "Horizon average should approach M";
}

// Test: g_tt → 0 as r → r₊ (coordinate singularity)
TEST_F(ReissnerNordstromTests, EventHorizonGtt) {
    double Q = 0.5;
    double r_plus = outerHorizon(Q);
    
    // Approach from outside
    std::vector<double> factors = {1.2, 1.1, 1.05, 1.01, 1.001};
    
    for (double factor : factors) {
        double r = r_plus * factor;
        Metric4D g = createRNMetric(r, Q);
        
        // g_tt should be small and negative, approaching 0
        EXPECT_LT(g(0, 0).real, 0.0) << "g_tt should be negative at r=" << r;
        EXPECT_GT(g(0, 0).real, -1.0) << "g_tt should approach 0 as r → r₊";
    }
}

// Test: g_rr diverges as r → r₊
TEST_F(ReissnerNordstromTests, EventHorizonGrr) {
    double Q = 0.5;
    double r_plus = outerHorizon(Q);
    double r = r_plus * 1.001;  // Just outside horizon
    
    Metric4D g = createRNMetric(r, Q);
    
    // g_rr = 1/f(r) should be very large near horizon
    EXPECT_GT(g(1, 1).real, 100.0) << "g_rr should blow up near horizon";
}

// =============================================================================
// Angular Components
// =============================================================================

// Test: g_θθ = r² (independent of charge, same as Schwarzschild)
TEST_F(ReissnerNordstromTests, MetricThetaTheta) {
    std::vector<double> radii = {3.0, 5.0, 10.0, 100.0};
    std::vector<double> charges = {0.0, 0.3, 0.5, 0.9};
    
    for (double r : radii) {
        for (double Q : charges) {
            Metric4D g = createRNMetric(r, Q);
            
            EXPECT_NEAR(g(2, 2).real, r * r, kEpsilon)
                << "g_θθ should equal r² at r=" << r << ", Q=" << Q;
        }
    }
}

// Test: g_φφ = r² sin²θ (independent of charge)
TEST_F(ReissnerNordstromTests, MetricPhiPhi) {
    double r = 10.0;
    double Q = 0.7;
    std::vector<double> thetas = {M_PI/6, M_PI/4, M_PI/3, M_PI/2};
    
    for (double theta : thetas) {
        Metric4D g = createRNMetric(r, Q, theta);
        double sin_theta = std::sin(theta);
        
        EXPECT_NEAR(g(3, 3).real, r * r * sin_theta * sin_theta, kEpsilon)
            << "g_φφ incorrect at θ=" << theta;
    }
}

// =============================================================================
// Determinant
// =============================================================================

// Test: det(g) = -r⁴ sin²θ (same as Schwarzschild, charge cancels)
TEST_F(ReissnerNordstromTests, MetricDeterminant) {
    double r = 10.0;
    double theta = M_PI / 3;
    double Q = 0.6;
    Metric4D g = createRNMetric(r, Q, theta);
    
    // For diagonal metric: det(g) = g_tt * g_rr * g_θθ * g_φφ
    double det = g(0, 0).real * g(1, 1).real * g(2, 2).real * g(3, 3).real;
    
    // Expected: -r⁴ sin²θ (the f factors cancel: -f * (1/f) = -1)
    double sin_theta = std::sin(theta);
    double expected_det = -r * r * r * r * sin_theta * sin_theta;
    
    EXPECT_NEAR(det, expected_det, kEpsilon)
        << "Metric determinant incorrect";
}

// Test: Determinant is negative (Lorentzian) for all charges
TEST_F(ReissnerNordstromTests, DeterminantIsNegative) {
    std::vector<double> radii = {3.0, 5.0, 10.0, 100.0};
    std::vector<double> charges = {0.0, 0.3, 0.5, 0.9};
    
    for (double r : radii) {
        for (double Q : charges) {
            Metric4D g = createRNMetric(r, Q, M_PI/2);
            double det = g(0, 0).real * g(1, 1).real * g(2, 2).real * g(3, 3).real;
            
            EXPECT_LT(det, 0.0) 
                << "Determinant should be negative at r=" << r << ", Q=" << Q;
        }
    }
}

// =============================================================================
// Physical Invariants
// =============================================================================

// Test: Proper time for stationary observer
// dτ² = -g_tt dt² = f(r) dt²
TEST_F(ReissnerNordstromTests, ProperTimeStationaryObserver) {
    double r = 10.0;
    double Q = 0.5;
    double dt = 1.0;  // Coordinate time interval
    
    Metric4D g = createRNMetric(r, Q);
    
    // Proper time: dτ² = -g_tt dt²
    double dtau_squared = -g(0, 0).real * dt * dt;
    double dtau = std::sqrt(dtau_squared);
    
    // Expected: dτ = √f(r) dt
    double expected_dtau = std::sqrt(f_of_r(r, Q)) * dt;
    
    EXPECT_NEAR(dtau, expected_dtau, kEpsilon)
        << "Proper time gravitational redshift incorrect";
}

// Test: Charge increases redshift compared to Schwarzschild (for same r)
TEST_F(ReissnerNordstromTests, ChargeIncreasesRedshift) {
    double r = 5.0;  // Fixed radius
    
    // Compare f(r) for different charges
    double f_schw = f_of_r(r, 0.0);  // Schwarzschild
    double f_rn = f_of_r(r, 0.5);    // With charge
    
    // For RN, Q²/r² term is positive, so f(r) is larger than Schwarzschild
    // This means LESS redshift (observer is less deep in potential well)
    EXPECT_GT(f_rn, f_schw)
        << "Charge should increase f(r), reducing gravitational potential";
}

// Test: Effective potential structure (photon sphere differs from Schwarzschild)
TEST_F(ReissnerNordstromTests, PhotonSphereExists) {
    double Q = 0.5;
    
    // For RN, photon sphere is at r = 3M/2 * (1 + √(1 - 8Q²/9M²))
    // For Q=0, this gives r = 3M (Schwarzschild photon sphere)
    double Q_actual = Q * M;
    double discriminant = 1.0 - 8.0 * Q_actual * Q_actual / (9.0 * M * M);
    
    if (discriminant > 0) {
        double r_photon = 1.5 * M * (1.0 + std::sqrt(discriminant));
        
        // Verify it's outside outer horizon
        double r_plus = outerHorizon(Q);
        EXPECT_GT(r_photon, r_plus) 
            << "Photon sphere should be outside outer horizon";
    }
}

// =============================================================================
// Charge Comparison Tests
// =============================================================================

// Test: Increasing charge moves outer horizon inward
TEST_F(ReissnerNordstromTests, ChargeShrinksOuterHorizon) {
    std::vector<double> charges = {0.0, 0.3, 0.5, 0.7, 0.9};
    double prev_r_plus = std::numeric_limits<double>::max();
    
    for (double Q : charges) {
        double r_plus = outerHorizon(Q);
        
        if (Q > 0.0) {
            EXPECT_LT(r_plus, prev_r_plus)
                << "Outer horizon should shrink with increasing charge";
        }
        prev_r_plus = r_plus;
    }
}

// Test: Increasing charge moves inner horizon outward
TEST_F(ReissnerNordstromTests, ChargeExpandsInnerHorizon) {
    std::vector<double> charges = {0.1, 0.3, 0.5, 0.7, 0.9};
    double prev_r_minus = 0.0;
    
    for (double Q : charges) {
        double r_minus = innerHorizon(Q);
        
        EXPECT_GT(r_minus, prev_r_minus)
            << "Inner horizon should expand with increasing charge";
        prev_r_minus = r_minus;
    }
}

// =============================================================================
// f(r) Function Tests
// =============================================================================

// Test: f(r) formula correctness
TEST_F(ReissnerNordstromTests, FunctionFCorrectness) {
    std::vector<std::pair<double, double>> test_cases = {
        {10.0, 0.0},   // Schwarzschild limit
        {10.0, 0.5},   // Medium charge
        {5.0, 0.3},    // Different r
        {20.0, 0.9},   // Large r, high charge
    };
    
    for (auto [r, Q] : test_cases) {
        double f_computed = f_of_r(r, Q);
        double Q_actual = Q * M;
        double f_expected = 1.0 - 2.0 * M / r + (Q_actual * Q_actual) / (r * r);
        
        EXPECT_NEAR(f_computed, f_expected, kEpsilon)
            << "f(r) incorrect at r=" << r << ", Q=" << Q;
    }
}

// Test: f(r) is positive outside outer horizon
TEST_F(ReissnerNordstromTests, FPositiveOutsideHorizon) {
    std::vector<double> charges = {0.0, 0.3, 0.5, 0.7, 0.9};
    
    for (double Q : charges) {
        double r_plus = outerHorizon(Q);
        std::vector<double> distances = {1.01, 1.1, 1.5, 2.0, 5.0};
        
        for (double factor : distances) {
            double r = r_plus * factor;
            double f = f_of_r(r, Q);
            
            EXPECT_GT(f, 0.0) 
                << "f(r) should be positive at r=" << r << ", Q=" << Q;
        }
    }
}

} // namespace sirius::test
