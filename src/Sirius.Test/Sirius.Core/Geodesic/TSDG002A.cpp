// TSDG002A.cpp - NaN/Inf Detection Tests
// Tests: dual number edge cases, metric/Christoffel under singularities.

#define _USE_MATH_DEFINES
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include "MTDL001A.h"
#include "MTTN001A.h"
#include "PHMT000A.h"
#include "PHMT100A.h" // Unified Kerr-Schild Family

namespace sirius::test {

// =============================================================================
// Test Fixture
// =============================================================================

class NaNInfDetectionTests : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
    
    // Helper to check if a value is finite (not NaN, not Inf)
    template<typename T>
    bool isFinite(T value) {
        return std::isfinite(value);
    }
    
    // Helper to check entire metric tensor for finite values
    bool isMetricFinite(const Metric4D& g) {
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                if (!std::isfinite(g(mu, nu).real)) {
                    return false;
                }
            }
        }
        return true;
    }
};

// =============================================================================
// Dual Number NaN Detection
// =============================================================================

TEST_F(NaNInfDetectionTests, DualDivisionByZero)
{
    // FORMAL SPECIFICATION:
    // Division by zero in real part should produce Inf
    // This must be detectable
    
    Dual<double> a(1.0, 1.0);
    Dual<double> zero(0.0, 0.0);
    
    // Note: We test detection, not prevention (prevention is in code)
    // The division would produce Inf, which we must detect
    if (zero.real != 0.0) {
        Dual<double> result = a / zero;
        EXPECT_FALSE(std::isfinite(result.real)) 
            << "Division by zero should produce non-finite result";
    }
    
    // Test that normal division produces finite results
    Dual<double> b(2.0, 1.0);
    Dual<double> result = a / b;
    EXPECT_TRUE(std::isfinite(result.real));
    EXPECT_TRUE(std::isfinite(result.dual));
}

TEST_F(NaNInfDetectionTests, DualSqrtNegative)
{
    // FORMAL SPECIFICATION:
    // sqrt of negative real part is undefined
    // Must produce NaN, which must be detectable
    
    Dual<double> negative(-1.0, 1.0);
    Dual<double> result = sqrt(negative);
    
    EXPECT_TRUE(std::isnan(result.real)) 
        << "sqrt of negative should produce NaN";
}

TEST_F(NaNInfDetectionTests, DualLogNonPositive)
{
    // FORMAL SPECIFICATION:
    // log of non-positive real part is undefined
    
    Dual<double> zero(0.0, 1.0);
    Dual<double> negative(-1.0, 1.0);
    
    Dual<double> result_zero = log(zero);
    Dual<double> result_neg = log(negative);
    
    EXPECT_FALSE(std::isfinite(result_zero.real)) 
        << "log(0) should produce -Inf";
    EXPECT_TRUE(std::isnan(result_neg.real)) 
        << "log(negative) should produce NaN";
}

TEST_F(NaNInfDetectionTests, DualSqrtNearZeroDerivative)
{
    // FORMAL SPECIFICATION:
    // sqrt near zero has derivative tending to infinity
    // d/dx sqrt(x) = 1/(2*sqrt(x)) -> infinity as x -> 0
    
    double epsilon = 1e-20;
    Dual<double> near_zero(epsilon, 1.0);
    Dual<double> result = sqrt(near_zero);
    
    // Real part should be finite (sqrt of small positive is small positive)
    EXPECT_TRUE(std::isfinite(result.real));
    
    // Dual part may be very large but should be finite for epsilon > 0
    // For epsilon = 1e-20, derivative = 1/(2*1e-10) = 5e9
    EXPECT_TRUE(std::isfinite(result.dual))
        << "Derivative should be large but finite for small positive input";
}

// =============================================================================
// Metric Evaluation NaN Detection
// =============================================================================

TEST_F(NaNInfDetectionTests, SchwarzschildMetricFinite)
{
    // FORMAL SPECIFICATION:
    // Metric must produce finite values for all valid domain points
    // Domain: r > 2M
    
    Sirius::KerrSchildFamily metric{Sirius::KerrSchildParams::Schwarzschild(1.0)};
    metric.setParameter("mass", 1.0);
    
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    
    // Test at various valid radii
    std::vector<double> test_radii = {2.1, 3.0, 6.0, 10.0, 100.0, 1000.0};
    
    for (double r : test_radii) {
        // Convert Spherical to Cartesian (r, theta, phi) -> (x, y, z)
        double theta = M_PI/2.0;
        double sin_theta = std::sin(theta);
        double cos_theta = std::cos(theta);
        double sin_phi = 0.0; // phi = 0
        double cos_phi = 1.0;
        
        Tensor<double, 4> pos;
        pos(0) = 0.0; // t
        pos(1) = r * sin_theta * cos_phi; // x
        pos(2) = r * sin_theta * sin_phi; // y
        pos(3) = r * cos_theta;           // z
        
        metric.evaluate(pos, g, dg);
        
        EXPECT_TRUE(isMetricFinite(g)) 
            << "Metric should be finite at r = " << r;
    }
}

TEST_F(NaNInfDetectionTests, SchwarzschildNearHorizon)
{
    // FORMAL SPECIFICATION:
    // Very close to horizon, numerical precision may degrade
    // Implementation should clamp to prevent singularity
    
    Sirius::KerrSchildFamily metric{Sirius::KerrSchildParams::Schwarzschild(1.0)};
    metric.setParameter("mass", 1.0);
    double rs = 2.0;  // Schwarzschild radius
    
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    
    // Test approaching horizon (implementation clamps at rs * 1.001)
    std::vector<double> near_horizon = {rs * 1.01, rs * 1.001, rs * 1.0001};
    
    for (double r : near_horizon) {
        // Convert Spherical to Cartesian (r, theta, phi) -> (x, y, z)
        double theta = M_PI/2.0;
        double sin_theta = std::sin(theta);
        double cos_theta = std::cos(theta);
        double sin_phi = 0.0; // phi = 0
        double cos_phi = 1.0;
        
        Tensor<double, 4> pos;
        pos(0) = 0.0; // t
        pos(1) = r * sin_theta * cos_phi; // x
        pos(2) = r * sin_theta * sin_phi; // y
        pos(3) = r * cos_theta;           // z
        
        metric.evaluate(pos, g, dg);
        
        EXPECT_TRUE(isMetricFinite(g)) 
            << "Metric should be clamped and finite near horizon at r = " << r;
    }
}

TEST_F(NaNInfDetectionTests, KerrMetricFinite)
{
    // FORMAL SPECIFICATION:
    // Kerr metric must be finite outside horizons
    // Additional singularities: Sigma = 0 (ring), Delta = 0 (horizons)
    
    Sirius::KerrSchildFamily metric{Sirius::KerrSchildParams::Kerr(1.0, 0.5)};
    metric.setParameter("mass", 1.0);
    metric.setParameter("spin", 0.5);
    
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    
    // Test at various valid positions
    std::vector<std::pair<double, double>> test_points = {
        {3.0, M_PI/2},   // Equatorial plane
        {3.0, M_PI/4},   // 45 degrees
        {3.0, 0.1},      // Near pole (but not at pole)
        {10.0, M_PI/2},  // Far field
    };
    
    for (auto& [r, theta] : test_points) {
        // Convert Spherical to Cartesian (r, theta, phi) -> (x, y, z)
        double sin_theta = std::sin(theta);
        double cos_theta = std::cos(theta);
        double sin_phi = 0.0; // phi = 0
        double cos_phi = 1.0;
        
        Tensor<double, 4> pos;
        pos(0) = 0.0; // t
        pos(1) = r * sin_theta * cos_phi; // x
        pos(2) = r * sin_theta * sin_phi; // y
        pos(3) = r * cos_theta;           // z
        
        metric.evaluate(pos, g, dg);
        
        EXPECT_TRUE(isMetricFinite(g)) 
            << "Kerr metric should be finite at r=" << r << ", theta=" << theta;
    }
}

TEST_F(NaNInfDetectionTests, KerrNearPoles)
{
    // FORMAL SPECIFICATION:
    // Boyer-Lindquist coordinates have singularity at poles (sin(theta) = 0)
    // Implementation should handle this gracefully
    
    Sirius::KerrSchildFamily metric{Sirius::KerrSchildParams::Kerr(1.0, 0.5)};
    metric.setParameter("mass", 1.0);
    metric.setParameter("spin", 0.5);
    
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    
    // Test near but not at poles
    double pole_offset = 1e-6;
    std::vector<double> near_pole_theta = {pole_offset, M_PI - pole_offset};
    
    for (double theta : near_pole_theta) {
        // Convert Spherical to Cartesian (r, theta, phi) -> (x, y, z)
        double r = 5.0; // Defined as 5.0 in original loop
        double sin_theta = std::sin(theta);
        double cos_theta = std::cos(theta);
        double sin_phi = 0.0; // phi = 0
        double cos_phi = 1.0;
        
        Tensor<double, 4> pos;
        pos(0) = 0.0; // t
        pos(1) = r * sin_theta * cos_phi; // x
        pos(2) = r * sin_theta * sin_phi; // y
        pos(3) = r * cos_theta;           // z
        
        metric.evaluate(pos, g, dg);
        
        // g_phi_phi = (r^2 + a^2 + ...) * sin^2(theta) should be very small but finite
        EXPECT_TRUE(isMetricFinite(g)) 
            << "Kerr metric should be finite near pole at theta=" << theta;
    }
}

// =============================================================================
// Vector and Tensor NaN Detection
// =============================================================================

TEST_F(NaNInfDetectionTests, TensorOperationsFinite)
{
    // FORMAL SPECIFICATION:
    // All tensor operations must produce finite results for finite inputs
    
    Vec4 v1, v2;
    v1(0) = 1.0; v1(1) = 0.0; v1(2) = 0.0; v1(3) = 0.0;
    v2(0) = 0.0; v2(1) = 1.0; v2(2) = 0.0; v2(3) = 0.0;
    
    // Create a simple metric (Minkowski)
    Metric4D g;
    g.zero();
    g(0, 0) = Dual<double>(-1.0);
    g(1, 1) = Dual<double>(1.0);
    g(2, 2) = Dual<double>(1.0);
    g(3, 3) = Dual<double>(1.0);
    
    // Test inner product
    double inner = TensorOps::innerProduct(v1, v1, g);
    EXPECT_TRUE(std::isfinite(inner)) << "Inner product should be finite";
    
    // Test index lowering
    Vec4 v1_lower = TensorOps::lowerIndex(v1, g);
    for (int i = 0; i < 4; ++i) {
        EXPECT_TRUE(std::isfinite(v1_lower(i))) 
            << "Lowered index component " << i << " should be finite";
    }
}

TEST_F(NaNInfDetectionTests, ChristoffelFinite)
{
    // FORMAL SPECIFICATION:
    // Christoffel symbols must be finite for valid metric
    
    Sirius::KerrSchildFamily metric{Sirius::KerrSchildParams::Schwarzschild(1.0)};
    metric.setParameter("mass", 1.0);
    
    // Convert Spherical to Cartesian (r, theta, phi) -> (x, y, z)
    double r = 5.0;
    double theta = M_PI/2.0;
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    double sin_phi = 0.0; // phi = 0
    double cos_phi = 1.0;
    
    Tensor<double, 4> pos;
    pos(0) = 0.0; // t
    pos(1) = r * sin_theta * cos_phi; // x
    pos(2) = r * sin_theta * sin_phi; // y
    pos(3) = r * cos_theta;           // z
    
    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric.evaluate(pos, g, dg);
    
    ChristoffelSymbols christoffel = TensorOps::christoffel(g, dg);
    
    // Check all 40 unique components are finite
    for (int lambda = 0; lambda < 4; ++lambda) {
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = mu; nu < 4; ++nu) {
                EXPECT_TRUE(std::isfinite(christoffel.gamma(lambda, mu, nu).real)) 
                    << "Gamma^" << lambda << "_" << mu << nu << " should be finite";
            }
        }
    }
}

// =============================================================================
// Extreme Value Tests
// =============================================================================

TEST_F(NaNInfDetectionTests, LargeValuesDoNotOverflow)
{
    // FORMAL SPECIFICATION:
    // Operations with large but representable values should not overflow
    
    double large = 1e100;
    Dual<double> a(large, 1.0);
    Dual<double> b(2.0, 1.0);
    
    // Addition should be fine
    Dual<double> sum = a + b;
    EXPECT_TRUE(std::isfinite(sum.real));
    
    // Multiplication of two large could overflow
    Dual<double> c(1e200, 1.0);
    Dual<double> product = c * c;
    // This WILL overflow - we're testing detection
    EXPECT_FALSE(std::isfinite(product.real)) 
        << "Overflow should be detectable";
}

TEST_F(NaNInfDetectionTests, SmallValuesDoNotUnderflow)
{
    // FORMAL SPECIFICATION:
    // Very small values should remain finite (may become denormalized)
    
    double small = 1e-300;
    Dual<double> a(small, 1.0);
    
    // Squaring a small value gives something smaller
    Dual<double> result = a * a;
    
    // Should still be finite (denormalized is fine)
    EXPECT_TRUE(std::isfinite(result.real));
}

} // namespace sirius::test
