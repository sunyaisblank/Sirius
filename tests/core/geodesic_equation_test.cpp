// Geodesic Equation Tests
// Tests: flat space acceleration, null normalization, timelike/null geodesics.
// Ported from TSPH003A.cpp; assertions and tolerances unchanged.

#include "sirius/core/dual_number.h"
#include "sirius/core/tensor.h"

#include <gtest/gtest.h>

#include <cmath>

namespace sirius::test {
using namespace sirius::core;

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
    Metric4d eta;
    eta.Zero();
    eta(0, 0) = Dual<double>(-1.0, 0.0);
    eta(1, 1) = Dual<double>(1.0, 0.0);
    eta(2, 2) = Dual<double>(1.0, 0.0);
    eta(3, 3) = Dual<double>(1.0, 0.0);

    // Zero derivatives
    Tensor<Dual<double>, 4, 4, 4> dg;
    dg.Zero();

    ChristoffelSymbols gamma = TensorOps::Christoffel(eta, dg);

    // Any velocity
    Vec4 velocity;
    velocity(0) = 1.0;
    velocity(1) = 0.5;
    velocity(2) = 0.3;
    velocity(3) = 0.1;

    Vec4 accel = TensorOps::GeodesicAcceleration(velocity, gamma);

    // All accelerations should be zero
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(accel(i), 0.0, kEpsilon)
            << "Acceleration component " << i << " should be zero in flat space";
    }
}

// =============================================================================
// Null Vector Normalization
// =============================================================================

// Test: NormalizeNull preserves null condition
TEST_F(GeodesicTests, NormalizeNullPreservesCondition) {
    Metric4d eta;
    eta.Zero();
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

    Vec4 normalized = TensorOps::NormalizeNull(k, eta);
    double norm_sq = TensorOps::InnerProduct(normalized, normalized, eta);

    EXPECT_NEAR(norm_sq, 0.0, kEpsilon) << "Normalized null vector should have zero norm squared";

    Vec4 past = k;
    past(0) = -0.8;  // identifies the past-directed null family.
    const Vec4 past_normalized = TensorOps::NormalizeNull(past, eta);
    EXPECT_LT(past_normalized(0), 0.0);
    EXPECT_NEAR(TensorOps::InnerProduct(past_normalized, past_normalized, eta), 0.0, kEpsilon);
}

// =============================================================================
// Schwarzschild null-vector algebra
// =============================================================================

// Test: the TensorOps contraction recognises the exact radial null direction
// in a Schwarzschild Boyer-Lindquist metric sample. Integration conservation
// is owned by the live-path conservation suites.
TEST_F(GeodesicTests, RadialNullVectorSatisfiesSchwarzschildForm) {
    // Create Schwarzschild metric at r = 10M
    double r = 10.0;
    double M = 1.0;
    double rs = 2.0 * M;
    double f = 1.0 - rs / r;

    Metric4d g;
    g.Zero();
    g(0, 0) = Dual<double>(-f, 0.0);
    g(1, 1) = Dual<double>(1.0 / f, 0.0);
    g(2, 2) = Dual<double>(r * r, 0.0);
    g(3, 3) = Dual<double>(r * r, 0.0);  // θ = π/2

    // Create a null vector for radial photon
    // For radial null geodesic: -f(dt)² + (1/f)(dr)² = 0
    // => dr/dt = f (outgoing photon)
    Vec4 k;
    k(0) = 1.0;  // dt/dλ
    k(1) = f;    // dr/dλ = f * dt/dλ
    k(2) = 0.0;
    k(3) = 0.0;

    double norm_sq = TensorOps::InnerProduct(k, k, g);

    // Should be zero for null geodesic
    EXPECT_NEAR(norm_sq, 0.0, kEpsilon) << "Radial null geodesic should satisfy null condition";
}

// Test: Timelike vector normalization
TEST_F(GeodesicTests, TimelikeVectorNormalization) {
    // Minkowski metric
    Metric4d eta;
    eta.Zero();
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

    double norm_sq = TensorOps::InnerProduct(u, u, eta);

    // Timelike: g_μν u^μ u^ν = -1
    EXPECT_NEAR(norm_sq, -1.0, kEpsilon) << "Stationary observer 4-velocity should have norm² = -1";
}

// Test: Moving observer with Lorentz boost
TEST_F(GeodesicTests, LorentzBoostedObserver) {
    Metric4d eta;
    eta.Zero();
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

    double norm_sq = TensorOps::InnerProduct(u, u, eta);

    // Should still be -1 for timelike
    EXPECT_NEAR(norm_sq, -1.0, kEpsilon) << "Boosted observer should maintain norm² = -1";
}

}  // namespace sirius::test
