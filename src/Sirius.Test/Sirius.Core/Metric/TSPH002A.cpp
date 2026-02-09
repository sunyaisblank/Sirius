// TSPH002A.cpp - Christoffel Symbol Tests
// Tests: torsion-free symmetry, flat space zeros, spherical coordinate values.

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

class ChristoffelTests : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

// =============================================================================
// Minkowski Spacetime (All Christoffel = 0)
// =============================================================================

// Test: Flat space has zero Christoffel symbols
TEST_F(ChristoffelTests, FlatSpaceChristoffelAllZero) {
    // Create Minkowski metric
    Metric4D eta;
    eta.zero();
    eta(0, 0) = Dual<double>(-1.0, 0.0);
    eta(1, 1) = Dual<double>(1.0, 0.0);
    eta(2, 2) = Dual<double>(1.0, 0.0);
    eta(3, 3) = Dual<double>(1.0, 0.0);
    
    // Zero derivatives (constant metric)
    Tensor<Dual<double>, 4, 4, 4> dg;
    dg.zero();
    
    ChristoffelSymbols gamma = TensorOps::christoffel(eta, dg);
    
    // Count non-zero elements
    int nonZeroCount = 0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                if (std::abs(gamma.gamma(i, j, k).real) > kEpsilon) {
                    nonZeroCount++;
                }
            }
        }
    }
    
    EXPECT_EQ(nonZeroCount, 0) 
        << "Flat space should have all zero Christoffel symbols";
}

// =============================================================================
// Torsion-Free Property
// =============================================================================

// Test: Γᵘ_αβ = Γᵘ_βα (symmetry in lower indices)
TEST_F(ChristoffelTests, TorsionFreeSymmetry) {
    // Create non-trivial metric (spherical coordinates in flat space)
    // ds² = -dt² + dr² + r²dθ² + r²sin²θ dφ²
    double r = 5.0;
    double theta = M_PI / 4;
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    
    Metric4D g;
    g.zero();
    g(0, 0) = Dual<double>(-1.0, 0.0);
    g(1, 1) = Dual<double>(1.0, 0.0);
    g(2, 2) = Dual<double>(r * r, 0.0);
    g(3, 3) = Dual<double>(r * r * sin_theta * sin_theta, 0.0);
    
    // Create metric derivatives dg[σ,μ,ν] = ∂_σ g_μν
    Tensor<Dual<double>, 4, 4, 4> dg;
    dg.zero();
    // ∂_r g_θθ = 2r
    dg(1, 2, 2) = Dual<double>(2.0 * r, 0.0);
    // ∂_r g_φφ = 2r sin²θ
    dg(1, 3, 3) = Dual<double>(2.0 * r * sin_theta * sin_theta, 0.0);
    // ∂_θ g_φφ = 2r² sinθ cosθ
    dg(2, 3, 3) = Dual<double>(2.0 * r * r * sin_theta * cos_theta, 0.0);
    
    ChristoffelSymbols gamma = TensorOps::christoffel(g, dg);
    
    // Check symmetry for all combinations
    for (int lam = 0; lam < 4; ++lam) {
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                EXPECT_NEAR(gamma.gamma(lam, mu, nu).real,
                           gamma.gamma(lam, nu, mu).real, kEpsilon)
                    << "Torsion-free symmetry violated at Γ^" 
                    << lam << "_" << mu << nu;
            }
        }
    }
}

// =============================================================================
// Known Analytic Results for Spherical Coordinates
// =============================================================================

// Test: Γʳ_θθ = -r in flat spherical coordinates
TEST_F(ChristoffelTests, SphericalGammaRThetaTheta) {
    double r = 5.0;
    double theta = M_PI / 3;
    double sin_theta = std::sin(theta);
    
    Metric4D g;
    g.zero();
    g(0, 0) = Dual<double>(-1.0, 0.0);
    g(1, 1) = Dual<double>(1.0, 0.0);
    g(2, 2) = Dual<double>(r * r, 0.0);
    g(3, 3) = Dual<double>(r * r * sin_theta * sin_theta, 0.0);
    
    Tensor<Dual<double>, 4, 4, 4> dg;
    dg.zero();
    dg(1, 2, 2) = Dual<double>(2.0 * r, 0.0);
    dg(1, 3, 3) = Dual<double>(2.0 * r * sin_theta * sin_theta, 0.0);
    dg(2, 3, 3) = Dual<double>(2.0 * r * r * sin_theta * std::cos(theta), 0.0);
    
    ChristoffelSymbols gamma = TensorOps::christoffel(g, dg);
    
    // Γʳ_θθ = -r (r = index 1, θ = index 2)
    EXPECT_NEAR(gamma.gamma(1, 2, 2).real, -r, kEpsilon)
        << "Γʳ_θθ should equal -r";
}

// Test: Γᶿ_rθ = 1/r in flat spherical coordinates
TEST_F(ChristoffelTests, SphericalGammaThetaRTheta) {
    double r = 5.0;
    double theta = M_PI / 3;
    double sin_theta = std::sin(theta);
    
    Metric4D g;
    g.zero();
    g(0, 0) = Dual<double>(-1.0, 0.0);
    g(1, 1) = Dual<double>(1.0, 0.0);
    g(2, 2) = Dual<double>(r * r, 0.0);
    g(3, 3) = Dual<double>(r * r * sin_theta * sin_theta, 0.0);
    
    Tensor<Dual<double>, 4, 4, 4> dg;
    dg.zero();
    dg(1, 2, 2) = Dual<double>(2.0 * r, 0.0);
    dg(1, 3, 3) = Dual<double>(2.0 * r * sin_theta * sin_theta, 0.0);
    dg(2, 3, 3) = Dual<double>(2.0 * r * r * sin_theta * std::cos(theta), 0.0);
    
    ChristoffelSymbols gamma = TensorOps::christoffel(g, dg);
    
    // Γᶿ_rθ = 1/r (θ = index 2, r = index 1)
    EXPECT_NEAR(gamma.gamma(2, 1, 2).real, 1.0 / r, kEpsilon)
        << "Γᶿ_rθ should equal 1/r";
}

} // namespace sirius::test
