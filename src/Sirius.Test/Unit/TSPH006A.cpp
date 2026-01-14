// TSPH006A.cpp - Ellis Drainhole (Wormhole) Metric Tests
// Tests: throat geometry, redshift function, asymptotic flatness.

#define _USE_MATH_DEFINES
#include <cmath>
#include <gtest/gtest.h>

constexpr double kEpsilon = 1e-6;
constexpr double M_PI_VAL = 3.14159265358979323846;

// =============================================================================
// Test Fixture
// =============================================================================

class EllisDrainholeTests : public ::testing::Test {
protected:
    // Throat factor
    double alpha(double m, double n) {
        if (n <= m) return 0.01;
        return std::sqrt(n*n - m*m);
    }
    
    // Pseudo-angle
    double pseudophi(double r, double m, double n) {
        double a = alpha(m, n);
        return (n / a) * (M_PI_VAL/2.0 - std::atan2(r - m, a));
    }
    
    // Redshift function
    double Fp(double r, double m, double n) {
        double psi = pseudophi(r, m, n);
        double exp_term = std::exp(-(2.0 * m / n) * psi);
        double Fp_sq = 1.0 - exp_term;
        if (Fp_sq < 1e-10) return -1e-5;
        return -std::sqrt(Fp_sq);
    }
    
    // Areal radius
    double Rp(double r, double m, double n) {
        double a = alpha(m, n);
        double fp = Fp(r, m, n);
        double numerator = (r - m) * (r - m) + a * a;
        double denominator = 1.0 - fp*fp;
        if (denominator < 1e-10) denominator = 1e-10;
        return std::sqrt(numerator / denominator);
    }
};

// =============================================================================
// Throat Geometry Tests
// =============================================================================

TEST_F(EllisDrainholeTests, ThroatFactorPositive) {
    // alpha = sqrt(n^2 - m^2) should be positive for n > m
    double m = 0.5, n = 1.0;
    double a = alpha(m, n);
    
    EXPECT_GT(a, 0.0) << "Throat factor should be positive";
    EXPECT_NEAR(a, std::sqrt(0.75), kEpsilon) << "alpha = sqrt(1 - 0.25) = sqrt(0.75)";
}

TEST_F(EllisDrainholeTests, ThroatMinimumRadius) {
    // At r = m, the areal radius should be at or near minimum
    double m = 0.5, n = 1.0;
    
    double Rp_at_throat = Rp(m, m, n);
    double Rp_far = Rp(m + 10.0, m, n);
    
    EXPECT_LT(Rp_at_throat, Rp_far) 
        << "Areal radius at throat should be smaller than far from throat";
}

TEST_F(EllisDrainholeTests, ArealRadiusSymmetry) {
    // Wormhole connects two asymptotic regions
    // Rp(r) should have minimum at r ≈ m
    double m = 0.5, n = 1.0;
    
    double Rp_neg = Rp(m - 5.0, m, n);  // "Other side" of wormhole
    double Rp_throat = Rp(m, m, n);
    double Rp_pos = Rp(m + 5.0, m, n);
    
    // Both far regions should have similar large radii
    EXPECT_GT(Rp_neg, Rp_throat);
    EXPECT_GT(Rp_pos, Rp_throat);
}

// =============================================================================
// Redshift Function Tests
// =============================================================================

TEST_F(EllisDrainholeTests, RedshiftFunctionNegative) {
    // Fp should be negative (by construction)
    double m = 0.5, n = 1.0;
    
    EXPECT_LT(Fp(0.0, m, n), 0.0);
    EXPECT_LT(Fp(m, m, n), 0.0);
    EXPECT_LT(Fp(10.0, m, n), 0.0);
}

TEST_F(EllisDrainholeTests, RedshiftFunctionBounded) {
    // |Fp| should be bounded by 1
    double m = 0.5, n = 1.0;
    
    for (double r = -5.0; r <= 15.0; r += 1.0) {
        double fp = Fp(r, m, n);
        EXPECT_LE(std::abs(fp), 1.0) 
            << "Redshift function should have magnitude <= 1 at r=" << r;
    }
}

// =============================================================================
// Asymptotic Flatness Tests
// =============================================================================

TEST_F(EllisDrainholeTests, AsymptoticFlatness) {
    // Far from wormhole, metric should approach flat spacetime
    double m = 0.5, n = 1.0;
    double r_far = 100.0;
    
    double fp = Fp(r_far, m, n);
    double g_tt = -(1.0 + fp*fp);
    
    // As r -> infinity, Fp -> 0, so g_tt -> -1 (Minkowski)
    EXPECT_NEAR(g_tt, -1.0, 0.1) 
        << "g_tt should approach -1 far from wormhole";
}

// =============================================================================
// Parameter Sensitivity Tests
// =============================================================================

TEST_F(EllisDrainholeTests, ThroatSizeIncreases) {
    // Increasing n (with fixed m) should increase throat size
    double m = 0.5;
    
    double a1 = alpha(m, 1.0);
    double a2 = alpha(m, 2.0);
    double a3 = alpha(m, 3.0);
    
    EXPECT_LT(a1, a2);
    EXPECT_LT(a2, a3);
}

TEST_F(EllisDrainholeTests, NMustBeGreaterThanM) {
    // For valid wormhole, n > m is required
    double m = 1.0, n = 0.9;  // Invalid: n < m
    
    double a = alpha(m, n);
    
    // With n < m, our implementation returns small positive value as fallback
    EXPECT_GT(a, 0.0) << "Implementation should handle n <= m gracefully";
}

// main() is provided by GTest::gtest_main
