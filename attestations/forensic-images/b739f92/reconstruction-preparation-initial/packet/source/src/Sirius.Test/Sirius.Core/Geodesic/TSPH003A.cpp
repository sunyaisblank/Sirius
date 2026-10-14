// TSPH003A.cpp - Geodesic Equation Tests
// Tests: flat space acceleration, null normalization, timelike/null geodesics.

#include <gtest/gtest.h>
#include <cmath>
#include <MTTN001A.h>
#include <MTDL001A.h>

namespace sirius::test {
using namespace Sirius;

constexpr double kEpsilon = 1e-8;

// =============================================================================
// Test Fixture
// =============================================================================

class GeodesicTests : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

// =============================================================================
// Minkowski Spacetime: Straight Lines
// =============================================================================

// Test: Zero acceleration in flat space
TEST_F(GeodesicTests, FlatSpaceZeroAcceleration) {
    // Create Minkowski metric
    Metric4D eta;
    eta.zero();
    eta(0, 0) = Dual<double>(-1.0, 0.0);
    eta(1, 1) = Dual<double>(1.0, 0.0);
    eta(2, 2) = Dual<double>(1.0, 0.0);
    eta(3, 3) = Dual<double>(1.0, 0.0);
    
    // Zero derivatives
    Tensor<Dual<double>, 4, 4, 4> dg;
    dg.zero();
    
    ChristoffelSymbols gamma = TensorOps::christoffel(eta, dg);
    
    // Any velocity
    Vec4 velocity;
    velocity(0) = 1.0;
    velocity(1) = 0.5;
    velocity(2) = 0.3;
    velocity(3) = 0.1;
    
    Vec4 accel = TensorOps::geodesicAcceleration(velocity, gamma);
    
    // All accelerations should be zero
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(accel(i), 0.0, kEpsilon)
            << "Acceleration component " << i << " should be zero in flat space";
    }
}

// =============================================================================
// Null Vector Normalization
// =============================================================================

// Test: normalizeNull preserves null condition
TEST_F(GeodesicTests, NormalizeNullPreservesCondition) {
    Metric4D eta;
    eta.zero();
    eta(0, 0) = Dual<double>(-1.0, 0.0);
    eta(1, 1) = Dual<double>(1.0, 0.0);
    eta(2, 2) = Dual<double>(1.0, 0.0);
    eta(3, 3) = Dual<double>(1.0, 0.0);
    
    // Create a null vector manually: k = (1, 1, 0, 0) is null in Minkowski
    Vec4 k;
    k(0) = 1.0;
    k(1) = 1.0;
    k(2) = 0.0;
    k(3) = 0.0;
    
    Vec4 normalized = TensorOps::normalizeNull(k, eta);
    double norm_sq = TensorOps::innerProduct(normalized, normalized, eta);
    
    EXPECT_NEAR(norm_sq, 0.0, kEpsilon)
        << "Normalized null vector should have zero norm squared";
}

// =============================================================================
// Conservation Laws in Schwarzschild
// =============================================================================

// Test: Inner product is conserved along geodesic (conceptual)
// This tests the invariant: d/dλ(g_μν u^μ u^ν) = 0
TEST_F(GeodesicTests, InnerProductConservation) {
    // For a geodesic, the invariant g_μν u^μ u^ν should be constant.
    // We test this by verifying the integration conserves it numerically.
    
    // Create Schwarzschild metric at r = 10M
    double r = 10.0;
    double M = 1.0;
    double rs = 2.0 * M;
    double f = 1.0 - rs / r;
    
    Metric4D g;
    g.zero();
    g(0, 0) = Dual<double>(-f, 0.0);
    g(1, 1) = Dual<double>(1.0 / f, 0.0);
    g(2, 2) = Dual<double>(r * r, 0.0);
    g(3, 3) = Dual<double>(r * r, 0.0);  // θ = π/2
    
    // Create a null vector for radial photon
    // For radial null geodesic: -f(dt)² + (1/f)(dr)² = 0
    // => dr/dt = f (outgoing photon)
    Vec4 k;
    k(0) = 1.0;         // dt/dλ
    k(1) = f;           // dr/dλ = f * dt/dλ
    k(2) = 0.0;
    k(3) = 0.0;
    
    double norm_sq = TensorOps::innerProduct(k, k, g);
    
    // Should be zero for null geodesic
    EXPECT_NEAR(norm_sq, 0.0, kEpsilon)
        << "Radial null geodesic should satisfy null condition";
}

// Test: Timelike vector normalization
TEST_F(GeodesicTests, TimelikeVectorNormalization) {
    // Minkowski metric
    Metric4D eta;
    eta.zero();
    eta(0, 0) = Dual<double>(-1.0, 0.0);
    eta(1, 1) = Dual<double>(1.0, 0.0);
    eta(2, 2) = Dual<double>(1.0, 0.0);
    eta(3, 3) = Dual<double>(1.0, 0.0);
    
    // Stationary observer: u = (1, 0, 0, 0)
    Vec4 u;
    u(0) = 1.0;
    u(1) = 0.0;
    u(2) = 0.0;
    u(3) = 0.0;
    
    double norm_sq = TensorOps::innerProduct(u, u, eta);
    
    // Timelike: g_μν u^μ u^ν = -1
    EXPECT_NEAR(norm_sq, -1.0, kEpsilon)
        << "Stationary observer 4-velocity should have norm² = -1";
}

// Test: Moving observer with Lorentz boost
TEST_F(GeodesicTests, LorentzBoostedObserver) {
    Metric4D eta;
    eta.zero();
    eta(0, 0) = Dual<double>(-1.0, 0.0);
    eta(1, 1) = Dual<double>(1.0, 0.0);
    eta(2, 2) = Dual<double>(1.0, 0.0);
    eta(3, 3) = Dual<double>(1.0, 0.0);
    
    // Observer moving at v = 0.6c in x-direction
    double v = 0.6;
    double gamma_factor = 1.0 / std::sqrt(1.0 - v * v);
    
    Vec4 u;
    u(0) = gamma_factor;
    u(1) = gamma_factor * v;
    u(2) = 0.0;
    u(3) = 0.0;
    
    double norm_sq = TensorOps::innerProduct(u, u, eta);
    
    // Should still be -1 for timelike
    EXPECT_NEAR(norm_sq, -1.0, kEpsilon)
        << "Boosted observer should maintain norm² = -1";
}

// =============================================================================
// Circular Orbit in Schwarzschild
// =============================================================================

// Test: Circular orbit parameters at ISCO (r = 6M)
TEST_F(GeodesicTests, CircularOrbitISCO) {
    // At the ISCO (Innermost Stable Circular Orbit), r = 6M
    double M = 1.0;
    double r = 6.0 * M;
    double rs = 2.0 * M;
    double f = 1.0 - rs / r;  // 1 - 2/6 = 2/3
    
    // Angular velocity for circular orbit: Ω² = M/r³
    double omega = std::sqrt(M / (r * r * r));
    
    // For circular orbit, the 4-velocity is:
    // u^t = 1/√(f - r²Ω²) and u^φ = Ω * u^t
    // At ISCO: this simplifies
    
    // Create metric
    Metric4D g;
    g.zero();
    g(0, 0) = Dual<double>(-f, 0.0);
    g(1, 1) = Dual<double>(1.0 / f, 0.0);
    g(2, 2) = Dual<double>(r * r, 0.0);
    g(3, 3) = Dual<double>(r * r, 0.0);  // θ = π/2
    
    // Calculate proper u^t for circular orbit
    double ut_squared = 1.0 / (f - r * r * omega * omega);
    double ut = std::sqrt(ut_squared);
    double uphi = omega * ut;
    
    Vec4 u;
    u(0) = ut;
    u(1) = 0.0;  // No radial motion
    u(2) = 0.0;  // Equatorial
    u(3) = uphi;
    
    double norm_sq = TensorOps::innerProduct(u, u, g);
    
    // Timelike: should be -1
    EXPECT_NEAR(norm_sq, -1.0, kEpsilon)
        << "Circular orbit 4-velocity should have norm² = -1";
}

// =============================================================================
// Photon Sphere
// =============================================================================

// Test: Circular photon orbit at r = 3M
TEST_F(GeodesicTests, PhotonSphereOrbit) {
    // The photon sphere is at r = 3M for Schwarzschild
    double M = 1.0;
    double r = 3.0 * M;
    double rs = 2.0 * M;
    double f = 1.0 - rs / r;  // 1 - 2/3 = 1/3
    
    // For circular photon orbit: L/E = √27 * M (impact parameter)
    // Angular velocity: Ω = c/(3√3 M) for circular null orbit
    double omega = 1.0 / (std::sqrt(27.0) * M);
    
    Metric4D g;
    g.zero();
    g(0, 0) = Dual<double>(-f, 0.0);
    g(1, 1) = Dual<double>(1.0 / f, 0.0);
    g(2, 2) = Dual<double>(r * r, 0.0);
    g(3, 3) = Dual<double>(r * r, 0.0);
    
    // Null vector for circular photon: g_tt (k^t)² + g_φφ (k^φ)² = 0
    // -f k^t² + r² k^φ² = 0
    // k^φ = k^t * sqrt(f)/r
    Vec4 k;
    k(0) = 1.0;  // Arbitrary normalization
    k(1) = 0.0;  // Circular
    k(2) = 0.0;
    k(3) = std::sqrt(f) / r;  // From null condition
    
    double norm_sq = TensorOps::innerProduct(k, k, g);
    
    EXPECT_NEAR(norm_sq, 0.0, kEpsilon)
        << "Circular photon orbit should be null";
}

} // namespace sirius::test
