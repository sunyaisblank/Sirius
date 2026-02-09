// =============================================================================
// TSMT001A.cpp - Dual Number Arithmetic Tests
// Component ID: TSMT001A (Test/Unit/DualNumber)
// =============================================================================
//
// PURPOSE:
// Validates dual number arithmetic for automatic differentiation.
// Tests ring axioms (nilpotency ε²=0), chain rule, elementary functions.
//
// MATHEMATICAL BASIS:
// Dual numbers: D = {a + bε : ε² = 0}
// Forward-mode autodiff: f(a + bε) = f(a) + b·f'(a)·ε
//
// LABEL: Mandatory;Correctness
// =============================================================================

#define _USE_MATH_DEFINES
#include <cmath>
#include <gtest/gtest.h>
#include <MTDL001A.h>
#include <PHCN001A.h>  // Centralized constants

namespace sirius::test {
using namespace Sirius;

// Numerical tolerance for dual number comparisons
// Dual arithmetic is exact to machine epsilon when inputs are exact
constexpr double kEpsilon = 1e-12;  // Reasonable for double precision operations

// =============================================================================
// Test Fixture
// =============================================================================

class DualNumberTests : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

// =============================================================================
// Ring Axioms
// =============================================================================

// Test: Additive identity (a + bε) + 0 = (a + bε)
TEST_F(DualNumberTests, AdditiveIdentity) {
    Dual<double> x(3.5, 2.0);
    Dual<double> zero(0.0, 0.0);
    Dual<double> result = x + zero;
    
    EXPECT_NEAR(result.real, 3.5, kEpsilon);
    EXPECT_NEAR(result.dual, 2.0, kEpsilon);
}

// Test: Multiplicative identity (a + bε) * 1 = (a + bε)
TEST_F(DualNumberTests, MultiplicativeIdentity) {
    Dual<double> x(3.5, 2.0);
    Dual<double> one(1.0, 0.0);
    Dual<double> result = x * one;
    
    EXPECT_NEAR(result.real, 3.5, kEpsilon);
    EXPECT_NEAR(result.dual, 2.0, kEpsilon);
}

// Test: Additive inverse (a + bε) + (-a - bε) = 0
TEST_F(DualNumberTests, AdditiveInverse) {
    Dual<double> x(3.5, 2.0);
    Dual<double> neg_x = -x;
    Dual<double> result = x + neg_x;
    
    EXPECT_NEAR(result.real, 0.0, kEpsilon);
    EXPECT_NEAR(result.dual, 0.0, kEpsilon);
}

// Test: Commutativity of addition
TEST_F(DualNumberTests, AdditionCommutativity) {
    Dual<double> x(3.5, 2.0);
    Dual<double> y(1.2, 4.5);
    
    Dual<double> xy = x + y;
    Dual<double> yx = y + x;
    
    EXPECT_NEAR(xy.real, yx.real, kEpsilon);
    EXPECT_NEAR(xy.dual, yx.dual, kEpsilon);
}

// Test: Commutativity of multiplication
TEST_F(DualNumberTests, MultiplicationCommutativity) {
    Dual<double> x(3.5, 2.0);
    Dual<double> y(1.2, 4.5);
    
    Dual<double> xy = x * y;
    Dual<double> yx = y * x;
    
    EXPECT_NEAR(xy.real, yx.real, kEpsilon);
    EXPECT_NEAR(xy.dual, yx.dual, kEpsilon);
}

// Test: Associativity of addition
TEST_F(DualNumberTests, AdditionAssociativity) {
    Dual<double> x(1.0, 2.0);
    Dual<double> y(3.0, 4.0);
    Dual<double> z(5.0, 6.0);
    
    Dual<double> left = (x + y) + z;
    Dual<double> right = x + (y + z);
    
    EXPECT_NEAR(left.real, right.real, kEpsilon);
    EXPECT_NEAR(left.dual, right.dual, kEpsilon);
}

// Test: Distributivity a*(b+c) = a*b + a*c
TEST_F(DualNumberTests, Distributivity) {
    Dual<double> a(2.0, 1.0);
    Dual<double> b(3.0, 2.0);
    Dual<double> c(4.0, 3.0);
    
    Dual<double> left = a * (b + c);
    Dual<double> right = a * b + a * c;
    
    EXPECT_NEAR(left.real, right.real, kEpsilon);
    EXPECT_NEAR(left.dual, right.dual, kEpsilon);
}

// =============================================================================
// Nilpotency Property: ε² = 0
// =============================================================================

// Test: (0 + ε)² = 0 (pure dual squared is zero)
TEST_F(DualNumberTests, Nilpotency) {
    Dual<double> eps(0.0, 1.0);  // Pure dual: 0 + 1·ε
    Dual<double> result = eps * eps;
    
    // ε * ε = (0 + 1ε) * (0 + 1ε) = 0 + (0*1 + 1*0)ε = 0
    EXPECT_NEAR(result.real, 0.0, kEpsilon);
    EXPECT_NEAR(result.dual, 0.0, kEpsilon);
}

// =============================================================================
// Product Rule: (a + bε)(c + dε) = ac + (ad + bc)ε
// =============================================================================

// Test: Explicit product formula
TEST_F(DualNumberTests, ProductFormula) {
    double a = 3.0, b = 2.0, c = 5.0, d = 4.0;
    Dual<double> x(a, b);
    Dual<double> y(c, d);
    Dual<double> result = x * y;
    
    // Expected: ac + (ad + bc)ε
    double expected_real = a * c;          // 15
    double expected_dual = a * d + b * c;  // 12 + 10 = 22
    
    EXPECT_NEAR(result.real, expected_real, kEpsilon);
    EXPECT_NEAR(result.dual, expected_dual, kEpsilon);
}

// =============================================================================
// Division: (a + bε)/(c + dε) = (a/c) + ((bc - ad)/c²)ε
// =============================================================================

// Test: Division formula
TEST_F(DualNumberTests, DivisionFormula) {
    double a = 6.0, b = 8.0, c = 2.0, d = 1.0;
    Dual<double> x(a, b);
    Dual<double> y(c, d);
    Dual<double> result = x / y;
    
    // Expected: (a/c) + ((bc - ad)/c²)ε
    double expected_real = a / c;                        // 3
    double expected_dual = (b * c - a * d) / (c * c);   // (16 - 6) / 4 = 2.5
    
    EXPECT_NEAR(result.real, expected_real, kEpsilon);
    EXPECT_NEAR(result.dual, expected_dual, kEpsilon);
}

// =============================================================================
// Derivative of sin: d/dx[sin(x)] = cos(x)
// =============================================================================

// Test: sin(a + ε) = sin(a) + cos(a)·ε
TEST_F(DualNumberTests, SinDerivative) {
    double a = 0.5;  // Test point
    Dual<double> x(a, 1.0);  // x with dx = 1
    Dual<double> result = sin(x);
    
    EXPECT_NEAR(result.real, std::sin(a), kEpsilon);
    EXPECT_NEAR(result.dual, std::cos(a), kEpsilon);  // d/dx sin(x) = cos(x)
}

// Test: sin at multiple angles
TEST_F(DualNumberTests, SinDerivativeMultipleAngles) {
    std::vector<double> angles = {0.0, M_PI/6, M_PI/4, M_PI/3, M_PI/2, M_PI};
    
    for (double theta : angles) {
        Dual<double> x(theta, 1.0);
        Dual<double> result = sin(x);
        
        EXPECT_NEAR(result.real, std::sin(theta), kEpsilon)
            << "sin(" << theta << ") failed";
        EXPECT_NEAR(result.dual, std::cos(theta), kEpsilon)
            << "d/dx sin(" << theta << ") failed";
    }
}

// =============================================================================
// Derivative of cos: d/dx[cos(x)] = -sin(x)
// =============================================================================

// Test: cos(a + ε) = cos(a) - sin(a)·ε
TEST_F(DualNumberTests, CosDerivative) {
    double a = 0.7;
    Dual<double> x(a, 1.0);
    Dual<double> result = cos(x);
    
    EXPECT_NEAR(result.real, std::cos(a), kEpsilon);
    EXPECT_NEAR(result.dual, -std::sin(a), kEpsilon);  // d/dx cos(x) = -sin(x)
}

// =============================================================================
// Derivative of sqrt: d/dx[√x] = 1/(2√x)
// =============================================================================

// Test: sqrt(a + ε) = √a + 1/(2√a)·ε
TEST_F(DualNumberTests, SqrtDerivative) {
    double a = 4.0;
    Dual<double> x(a, 1.0);
    Dual<double> result = sqrt(x);
    
    EXPECT_NEAR(result.real, std::sqrt(a), kEpsilon);  // √4 = 2
    EXPECT_NEAR(result.dual, 1.0 / (2.0 * std::sqrt(a)), kEpsilon);  // 1/(2*2) = 0.25
}

// Test: sqrt derivative at various points
TEST_F(DualNumberTests, SqrtDerivativeMultiplePoints) {
    std::vector<double> values = {0.25, 1.0, 2.0, 9.0, 100.0};
    
    for (double v : values) {
        Dual<double> x(v, 1.0);
        Dual<double> result = sqrt(x);
        
        EXPECT_NEAR(result.real, std::sqrt(v), kEpsilon);
        EXPECT_NEAR(result.dual, 0.5 / std::sqrt(v), kEpsilon)
            << "d/dx sqrt(" << v << ") failed";
    }
}

// =============================================================================
// Chain Rule: d/dx[f(g(x))] = f'(g(x)) · g'(x)
// =============================================================================

// Test: d/dx[sin(x²)] = cos(x²) · 2x
TEST_F(DualNumberTests, ChainRuleSinSquare) {
    double x_val = 1.5;
    Dual<double> x(x_val, 1.0);
    
    // Compute x²
    Dual<double> x_squared = x * x;  // x² with derivative 2x
    
    // Compute sin(x²)
    Dual<double> result = sin(x_squared);
    
    // Expected: sin(x²), derivative = cos(x²) * 2x
    double expected_real = std::sin(x_val * x_val);
    double expected_dual = std::cos(x_val * x_val) * 2.0 * x_val;
    
    EXPECT_NEAR(result.real, expected_real, kEpsilon);
    EXPECT_NEAR(result.dual, expected_dual, kEpsilon);
}

// Test: d/dx[√(sin(x))] = cos(x) / (2√sin(x))
TEST_F(DualNumberTests, ChainRuleSqrtSin) {
    double x_val = 1.0;  // sin(1) ≈ 0.841 > 0
    Dual<double> x(x_val, 1.0);
    
    // Compute sqrt(sin(x))
    Dual<double> sin_x = sin(x);
    Dual<double> result = sqrt(sin_x);
    
    // Expected derivative: cos(x) / (2*sqrt(sin(x)))
    double sin_val = std::sin(x_val);
    double cos_val = std::cos(x_val);
    double expected_real = std::sqrt(sin_val);
    double expected_dual = cos_val / (2.0 * std::sqrt(sin_val));
    
    EXPECT_NEAR(result.real, expected_real, kEpsilon);
    EXPECT_NEAR(result.dual, expected_dual, kEpsilon);
}

// =============================================================================
// Scalar Operations
// =============================================================================

// Test: Scalar multiplication
TEST_F(DualNumberTests, ScalarMultiplication) {
    Dual<double> x(3.0, 2.0);
    double scalar = 5.0;
    
    Dual<double> result1 = x * scalar;
    Dual<double> result2 = scalar * x;
    
    EXPECT_NEAR(result1.real, 15.0, kEpsilon);
    EXPECT_NEAR(result1.dual, 10.0, kEpsilon);
    EXPECT_NEAR(result2.real, 15.0, kEpsilon);
    EXPECT_NEAR(result2.dual, 10.0, kEpsilon);
}

// Test: Scalar division
TEST_F(DualNumberTests, ScalarDivision) {
    Dual<double> x(10.0, 6.0);
    double scalar = 2.0;
    
    Dual<double> result = x / scalar;
    
    EXPECT_NEAR(result.real, 5.0, kEpsilon);
    EXPECT_NEAR(result.dual, 3.0, kEpsilon);
}

} // namespace sirius::test
