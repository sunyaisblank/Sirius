// Numerical Stability Diagnostics
// Tests edge cases: near-singularities, small/large values, precision limits.
// Ported from TSDG001A.cpp; assertions and tolerances unchanged.

#include "sirius/core/dual_number.h"
#include "sirius/core/tensor.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

namespace sirius::test {
using namespace sirius::core;

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
// Vector Operations
// =============================================================================

// Test: Vector normalization with small magnitude
TEST_F(NumericalStabilityTests, VectorSmallMagnitude) {
    Vec4 v;
    v(0) = 1e-8;
    v(1) = 1e-8;
    v(2) = 1e-8;
    v(3) = 1e-8;

    double len = v.Length();

    // Should be small but not zero or NaN
    EXPECT_GT(len, 0.0);
    EXPECT_FALSE(std::isnan(len));
    EXPECT_FALSE(std::isinf(len));
}

// Test: Vector operations with large components
TEST_F(NumericalStabilityTests, VectorLargeComponents) {
    Vec4 u, v;
    u(0) = 1e10;
    u(1) = 1e10;
    u(2) = 1e10;
    u(3) = 1e10;
    v(0) = 1e10;
    v(1) = 1e10;
    v(2) = 1e10;
    v(3) = 1e10;

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
    Metric4d eta;
    eta.Zero();
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

    double norm_sq = TensorOps::InnerProduct(k, k, eta);

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
    Metric4d g;
    g.Zero();
    double small = 1e-8;
    g(0, 0) = Dual<double>(-small, 0.0);       // Small but non-zero
    g(1, 1) = Dual<double>(1.0 / small, 0.0);  // Large
    g(2, 2) = Dual<double>(1.0, 0.0);
    g(3, 3) = Dual<double>(1.0, 0.0);

    Tensor<Dual<double>, 4, 4, 4> dg;
    dg.Zero();
    dg(1, 0, 0) = Dual<double>(1e-9, 0.0);  // Small derivative

    ChristoffelSymbols gamma = TensorOps::Christoffel(g, dg);

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

}  // namespace sirius::test
