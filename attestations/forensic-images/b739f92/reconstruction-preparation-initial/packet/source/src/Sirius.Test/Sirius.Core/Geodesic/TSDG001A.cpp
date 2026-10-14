// TSDG001A.cpp - Numerical Stability Diagnostics
// Tests edge cases: near-singularities, small/large values, precision limits.

#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <MTTN001A.h>
#include <MTDL001A.h>

namespace sirius::test {
using namespace Sirius;

constexpr double kEpsilon = 1e-10;

// =============================================================================
// Test Fixture
// =============================================================================

class NumericalStabilityTests : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

// =============================================================================
// Dual Number Edge Cases
// =============================================================================

// Test: Dual numbers handle zero correctly
TEST_F(NumericalStabilityTests, DualZeroHandling) {
    Dual<double> zero(0.0, 0.0);
    Dual<double> x(3.0, 2.0);
    
    Dual<double> sum = zero + x;
    Dual<double> prod = zero * x;
    
    EXPECT_NEAR(sum.real, 3.0, kEpsilon);
    EXPECT_NEAR(sum.dual, 2.0, kEpsilon);
    EXPECT_NEAR(prod.real, 0.0, kEpsilon);
    EXPECT_NEAR(prod.dual, 0.0, kEpsilon);
}

// Test: Dual division by small numbers
TEST_F(NumericalStabilityTests, DualDivisionSmallDenominator) {
    Dual<double> numerator(1.0, 1.0);
    Dual<double> small_denom(1e-6, 0.0);
    
    Dual<double> result = numerator / small_denom;
    
    // Should get large but finite result
    EXPECT_FALSE(std::isnan(result.real));
    EXPECT_FALSE(std::isnan(result.dual));
    EXPECT_FALSE(std::isinf(result.real));
    EXPECT_NEAR(result.real, 1e6, 1.0);  // Approximate
}

// Test: sqrt near zero (derivative blows up)
TEST_F(NumericalStabilityTests, SqrtNearZero) {
    double small = 1e-10;
    Dual<double> x(small, 1.0);
    
    Dual<double> result = sqrt(x);
    
    // sqrt(small) should be finite
    EXPECT_FALSE(std::isnan(result.real));
    EXPECT_NEAR(result.real, std::sqrt(small), kEpsilon);
    
    // Derivative 1/(2*sqrt(small)) will be large
    // Just check it's not NaN
    EXPECT_FALSE(std::isnan(result.dual));
}

// Test: sin/cos at large angles
TEST_F(NumericalStabilityTests, TrigLargeAngles) {
    double large_angle = 1e6;  // Many rotations
    Dual<double> x(large_angle, 1.0);
    
    Dual<double> s = sin(x);
    Dual<double> c = cos(x);
    
    // Results should be bounded [-1, 1]
    EXPECT_GE(s.real, -1.0 - kEpsilon);
    EXPECT_LE(s.real, 1.0 + kEpsilon);
    EXPECT_GE(c.real, -1.0 - kEpsilon);
    EXPECT_LE(c.real, 1.0 + kEpsilon);
    
    // Derivatives should also be bounded
    EXPECT_GE(s.dual, -1.0 - kEpsilon);  // cos
    EXPECT_LE(s.dual, 1.0 + kEpsilon);
    EXPECT_GE(c.dual, -1.0 - kEpsilon);  // -sin
    EXPECT_LE(c.dual, 1.0 + kEpsilon);
}

// =============================================================================
// Metric Singularities
// =============================================================================

// Test: Schwarzschild near horizon (r → 2M)
TEST_F(NumericalStabilityTests, SchwarzschildNearHorizon) {
    double M = 1.0;
    double rs = 2.0 * M;
    
    // Approach horizon from outside
    std::vector<double> radii = {2.1, 2.01, 2.001};
    
    for (double r : radii) {
        double f = 1.0 - rs / r;
        
        // g_tt should be small but not NaN
        double g_tt = -f;
        EXPECT_FALSE(std::isnan(g_tt)) << "g_tt NaN at r=" << r;
        EXPECT_GT(g_tt, -1.0) << "g_tt should approach 0 near horizon";
        
        // g_rr should be large but finite
        double g_rr = 1.0 / f;
        EXPECT_FALSE(std::isnan(g_rr)) << "g_rr NaN at r=" << r;
        EXPECT_FALSE(std::isinf(g_rr)) << "g_rr Inf at r=" << r;
    }
}

// Test: Polar axis (sin θ → 0)
TEST_F(NumericalStabilityTests, PolarAxisSingularity) {
    double r = 10.0;
    double M = 1.0;
    
    // Near pole
    std::vector<double> thetas = {0.01, 0.001, 0.0001};
    
    for (double theta : thetas) {
        double sin_theta = std::sin(theta);
        double g_phiphi = r * r * sin_theta * sin_theta;
        
        // Should be small but not negative
        EXPECT_GE(g_phiphi, 0.0) << "g_φφ should be non-negative";
        EXPECT_FALSE(std::isnan(g_phiphi)) << "g_φφ NaN near pole";
    }
}

// =============================================================================
// Vector Operations
// =============================================================================

// Test: Vector normalization with small magnitude
TEST_F(NumericalStabilityTests, VectorSmallMagnitude) {
    Vec4 v;
    v(0) = 1e-8;
    v(1) = 1e-8;
    v(2) = 1e-8;
    v(3) = 1e-8;
    
    double len = v.length();
    
    // Should be small but not zero or NaN
    EXPECT_GT(len, 0.0);
    EXPECT_FALSE(std::isnan(len));
    EXPECT_FALSE(std::isinf(len));
}

// Test: Vector operations with large components
TEST_F(NumericalStabilityTests, VectorLargeComponents) {
    Vec4 u, v;
    u(0) = 1e10; u(1) = 1e10; u(2) = 1e10; u(3) = 1e10;
    v(0) = 1e10; v(1) = 1e10; v(2) = 1e10; v(3) = 1e10;
    
    Vec4 w = u + v;
    
    // Should handle large numbers
    EXPECT_FALSE(std::isnan(w(0)));
    EXPECT_FALSE(std::isinf(w(0)));
    EXPECT_NEAR(w(0), 2e10, 1e5);
}

// =============================================================================
// Inner Product Stability
// =============================================================================

// Test: Inner product near null (almost cancellation)
TEST_F(NumericalStabilityTests, InnerProductNearNull) {
    Metric4D eta;
    eta.zero();
    eta(0, 0) = Dual<double>(-1.0, 0.0);
    eta(1, 1) = Dual<double>(1.0, 0.0);
    eta(2, 2) = Dual<double>(1.0, 0.0);
    eta(3, 3) = Dual<double>(1.0, 0.0);
    
    // Almost null vector
    Vec4 k;
    k(0) = 1.0;
    k(1) = 1.0 + 1e-15;  // Small deviation from null
    k(2) = 0.0;
    k(3) = 0.0;
    
    double norm_sq = TensorOps::innerProduct(k, k, eta);
    
    // Should be very small, not exactly zero due to perturbation
    EXPECT_NEAR(norm_sq, 0.0, 1e-10);
    EXPECT_FALSE(std::isnan(norm_sq));
}

// =============================================================================
// Christoffel Symbol Edge Cases
// =============================================================================

// Test: Christoffel with near-singular metric inverse
TEST_F(NumericalStabilityTests, ChristoffelNearSingular) {
    // Create metric approaching singularity
    Metric4D g;
    g.zero();
    double small = 1e-8;
    g(0, 0) = Dual<double>(-small, 0.0);  // Small but non-zero
    g(1, 1) = Dual<double>(1.0 / small, 0.0);  // Large
    g(2, 2) = Dual<double>(1.0, 0.0);
    g(3, 3) = Dual<double>(1.0, 0.0);
    
    Tensor<Dual<double>, 4, 4, 4> dg;
    dg.zero();
    dg(1, 0, 0) = Dual<double>(1e-9, 0.0);  // Small derivative
    
    ChristoffelSymbols gamma = TensorOps::christoffel(g, dg);
    
    // Results should be finite
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                EXPECT_FALSE(std::isnan(gamma.gamma(i, j, k).real))
                    << "Christoffel NaN at (" << i << "," << j << "," << k << ")";
            }
        }
    }
}

// =============================================================================
// Precision Tests
// =============================================================================

// Test: Machine epsilon handling
TEST_F(NumericalStabilityTests, MachineEpsilonAddition) {
    double x = 1.0;
    double eps = std::numeric_limits<double>::epsilon();
    
    // x + eps/2 should equal x (catastrophic cancellation)
    double y = x + eps / 2.0;
    
    // This tests understanding of floating point limits
    EXPECT_EQ(x, y) << "Addition below machine epsilon should not change value";
    
    // x + eps should differ from x
    double z = x + eps;
    EXPECT_NE(x, z) << "Addition at machine epsilon should change value";
}

// Test: Compensated summation concept (Kahan)
TEST_F(NumericalStabilityTests, AccumulationError) {
    // Accumulate many small numbers
    double naive_sum = 0.0;
    double small = 1e-15;
    int N = 1000000;
    
    for (int i = 0; i < N; ++i) {
        naive_sum += small;
    }
    
    double expected = small * N;  // 1e-9
    double relative_error = std::abs(naive_sum - expected) / expected;
    
    // Naive summation has some error, but should be reasonable
    EXPECT_LT(relative_error, 1e-5)
        << "Accumulation error too large: " << relative_error;
}

} // namespace sirius::test
