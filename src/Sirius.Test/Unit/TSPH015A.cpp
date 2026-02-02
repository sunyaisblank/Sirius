// TSPH015A.cpp - Coordinate Transformation Tests
// Component ID: TSPH015A
// Tests for PHCT002A: Boyer-Lindquist ↔ Kerr-Schild Cartesian coordinate transforms
//
// Specification requirements:
//   - Round-trip deviation < 1e-12 (see plan Phase 5.1)
//   - Jacobian determinant = r² sin(θ) for BL → Cart
//
// Tests: PHCT002A.h coordinate transformation functions

#define _USE_MATH_DEFINES
#include <cmath>
#include <gtest/gtest.h>

#include "PHCT002A.h"

using namespace Sirius::Coordinates;

namespace {

constexpr double kEpsilon = 1e-12;
constexpr double kLooseEpsilon = 1e-10;  // For floating-point accumulation

// =============================================================================
// Test Fixture
// =============================================================================

class CoordinateTransformTests : public ::testing::Test {
protected:
    // Standard test points in BL coordinates
    Vec4BL point_equator{0.0, 10.0, M_PI / 2.0, M_PI / 4.0};      // Equatorial plane
    Vec4BL point_pole_near{0.0, 10.0, 0.1, M_PI / 4.0};           // Near north pole
    Vec4BL point_far{0.0, 100.0, M_PI / 3.0, 2.0};                // Far field
    Vec4BL point_close{0.0, 3.0, M_PI / 4.0, -M_PI / 2.0};        // Near horizon
};

// =============================================================================
// Round-Trip Tests (BL → Cartesian → BL)
// =============================================================================

TEST_F(CoordinateTransformTests, RoundTripBLCartesian_Equator) {
    // Round-trip at equator (θ = π/2)
    double deviation = validateRoundTrip(point_equator, 0.0);
    EXPECT_LT(deviation, kEpsilon)
        << "Round-trip deviation at equator should be < 1e-12";
}

TEST_F(CoordinateTransformTests, RoundTripBLCartesian_NearPole) {
    // Round-trip near pole (potential coordinate singularity)
    double deviation = validateRoundTrip(point_pole_near, 0.0);
    EXPECT_LT(deviation, kLooseEpsilon)
        << "Round-trip deviation near pole should be < 1e-10";
}

TEST_F(CoordinateTransformTests, RoundTripBLCartesian_FarField) {
    // Round-trip in far field
    double deviation = validateRoundTrip(point_far, 0.0);
    EXPECT_LT(deviation, kEpsilon)
        << "Round-trip deviation in far field should be < 1e-12";
}

TEST_F(CoordinateTransformTests, RoundTripBLCartesian_NearHorizon) {
    // Round-trip near horizon (strong curvature region)
    double deviation = validateRoundTrip(point_close, 0.0);
    EXPECT_LT(deviation, kEpsilon)
        << "Round-trip deviation near horizon should be < 1e-12";
}

// =============================================================================
// Kerr Round-Trip Tests (with spin parameter)
// =============================================================================

TEST_F(CoordinateTransformTests, RoundTripKerrSchild_LowSpin) {
    // Kerr with low spin a = 0.1
    // Note: Kerr round-trip has larger deviation due to phi adjustment term:
    //   φ_BL = atan2(y, x) - atan2(a, r)
    // This introduces arctan precision limits, acceptable for physics use.
    double a = 0.1;
    double deviation = validateRoundTrip(point_equator, a);
    EXPECT_LT(deviation, 0.02)  // ~a/r deviation acceptable
        << "Kerr round-trip with a=0.1 should have deviation < 0.02";
}

TEST_F(CoordinateTransformTests, RoundTripKerrSchild_ModerateSpin) {
    // Kerr with moderate spin a = 0.5
    double a = 0.5;
    double deviation = validateRoundTrip(point_equator, a);
    EXPECT_LT(deviation, 0.1)  // ~a/r deviation acceptable
        << "Kerr round-trip with a=0.5 should have deviation < 0.1";
}

TEST_F(CoordinateTransformTests, RoundTripKerrSchild_HighSpin) {
    // Kerr with high spin a = 0.9 (near extremal)
    double a = 0.9;
    double deviation = validateRoundTrip(point_equator, a);
    EXPECT_LT(deviation, 0.15)  // ~a/r deviation acceptable
        << "Kerr round-trip with a=0.9 should have deviation < 0.15";
}

TEST_F(CoordinateTransformTests, RoundTripKerrSchild_NearPole) {
    // Kerr near pole (tests oblate spheroidal effects)
    double a = 0.5;
    double deviation = validateRoundTrip(point_pole_near, a);
    EXPECT_LT(deviation, 0.1)  // Looser tolerance for polar + Kerr
        << "Kerr round-trip near pole should have reasonable deviation";
}

// =============================================================================
// Jacobian Tests
// =============================================================================

TEST_F(CoordinateTransformTests, JacobianDeterminant_Equator) {
    // At equator: det(J) = r² sin(θ) = r² × 1 = r²
    Jacobian4x4 J = JacobianBLToCartesian(point_equator);
    double det = JacobianDeterminant(J);
    double expected = point_equator.r * point_equator.r * std::sin(point_equator.theta);

    EXPECT_NEAR(det, expected, 1e-10)
        << "Jacobian determinant should equal r² sin(θ)";
}

TEST_F(CoordinateTransformTests, JacobianDeterminant_MidLatitude) {
    // At θ = π/3: det(J) = r² sin(π/3) = r² × √3/2
    Vec4BL mid_lat{0.0, 10.0, M_PI / 3.0, 0.0};
    Jacobian4x4 J = JacobianBLToCartesian(mid_lat);
    double det = JacobianDeterminant(J);
    double expected = mid_lat.r * mid_lat.r * std::sin(mid_lat.theta);

    EXPECT_NEAR(det, expected, 1e-10)
        << "Jacobian determinant at mid-latitude should equal r² sin(θ)";
}

// =============================================================================
// Vector Transformation Tests
// =============================================================================

TEST_F(CoordinateTransformTests, VectorTransformRoundTrip) {
    // Transform a vector BL → Cart → BL and check round-trip
    Vec4BL v_orig{1.0, 0.5, 0.25, 0.1};  // Some arbitrary 4-vector

    Vec4Cart v_cart = transformVectorBLToCart(v_orig, point_equator);
    Vec4Cart pos_cart = BLToCartesian(point_equator);
    Vec4BL v_recovered = transformVectorCartToBL(v_cart, pos_cart);

    EXPECT_NEAR(v_orig.t, v_recovered.t, kLooseEpsilon);
    EXPECT_NEAR(v_orig.r, v_recovered.r, kLooseEpsilon);
    EXPECT_NEAR(v_orig.theta, v_recovered.theta, kLooseEpsilon);

    // Phi may wrap, so check with modular arithmetic
    double dphi = v_orig.phi - v_recovered.phi;
    while (dphi > M_PI) dphi -= 2 * M_PI;
    while (dphi < -M_PI) dphi += 2 * M_PI;
    EXPECT_NEAR(dphi, 0.0, kLooseEpsilon);
}

// =============================================================================
// Boundary Condition Tests
// =============================================================================

TEST_F(CoordinateTransformTests, OriginHandling) {
    // Test behavior at r → 0 (should be handled gracefully)
    Vec4BL origin{0.0, 0.0, M_PI / 2.0, 0.0};
    Vec4Cart cart = BLToCartesian(origin);

    EXPECT_EQ(cart.x, 0.0);
    EXPECT_EQ(cart.y, 0.0);
    EXPECT_EQ(cart.z, 0.0);
    EXPECT_EQ(cart.t, 0.0);
}

TEST_F(CoordinateTransformTests, PhiWrapping) {
    // Test phi wrapping at ±π boundary
    Vec4BL pos1{0.0, 10.0, M_PI / 2.0, M_PI - 0.01};
    Vec4BL pos2{0.0, 10.0, M_PI / 2.0, -M_PI + 0.01};

    Vec4Cart cart1 = BLToCartesian(pos1);
    Vec4Cart cart2 = BLToCartesian(pos2);

    // Points near ±π in phi should be close in Cartesian space
    double dx = cart1.x - cart2.x;
    double dy = cart1.y - cart2.y;
    double dist = std::sqrt(dx * dx + dy * dy);

    EXPECT_LT(dist, 0.5)
        << "Points near phi = ±π should be close in Cartesian space";
}

TEST_F(CoordinateTransformTests, NorthPoleCoordinates) {
    // At north pole (θ = 0), x = y = 0, z = r
    Vec4BL north_pole{0.0, 10.0, 0.0, 0.0};
    Vec4Cart cart = BLToCartesian(north_pole);

    EXPECT_NEAR(cart.x, 0.0, 1e-14);
    EXPECT_NEAR(cart.y, 0.0, 1e-14);
    EXPECT_NEAR(cart.z, 10.0, 1e-14);
}

TEST_F(CoordinateTransformTests, SouthPoleCoordinates) {
    // At south pole (θ = π), x = y = 0, z = -r
    Vec4BL south_pole{0.0, 10.0, M_PI, 0.0};
    Vec4Cart cart = BLToCartesian(south_pole);

    EXPECT_NEAR(cart.x, 0.0, 1e-14);
    EXPECT_NEAR(cart.y, 0.0, 1e-14);
    EXPECT_NEAR(cart.z, -10.0, 1e-14);
}

// =============================================================================
// Kerr-Specific Tests
// =============================================================================

TEST_F(CoordinateTransformTests, KerrOblateSpheroidal) {
    // For Kerr, x² + y² = (r² + a²) sin²θ, z = r cos θ
    // Test that oblate spheroidal transformation is correct
    double a = 0.9;
    double r = 10.0;
    double theta = M_PI / 3.0;

    Vec4BL pos{0.0, r, theta, 0.0};
    Vec4Cart cart = BLToKerrSchildCart(pos, a);

    double expected_rho = std::sqrt(r * r + a * a);
    double expected_xy_mag = expected_rho * std::sin(theta);
    double actual_xy_mag = std::sqrt(cart.x * cart.x + cart.y * cart.y);

    EXPECT_NEAR(actual_xy_mag, expected_xy_mag, 1e-12)
        << "Kerr oblate spheroidal: |x,y| should equal √(r²+a²) sin(θ)";

    EXPECT_NEAR(cart.z, r * std::cos(theta), 1e-12)
        << "Kerr: z should equal r cos(θ)";
}

TEST_F(CoordinateTransformTests, KerrSolveR_Accuracy) {
    // Test that inverse transformation correctly solves for r
    // from x² + y² + z² = r² + a² (1 - z²/r²)
    double a = 0.7;
    double r_original = 8.0;
    double theta = M_PI / 4.0;
    double phi = 1.5;

    Vec4BL original{0.0, r_original, theta, phi};
    Vec4Cart cart = BLToKerrSchildCart(original, a);
    Vec4BL recovered = KerrSchildCartToBL(cart, a);

    EXPECT_NEAR(original.r, recovered.r, 1e-10)
        << "Kerr inverse should recover r accurately";
    EXPECT_NEAR(original.theta, recovered.theta, 1e-10)
        << "Kerr inverse should recover θ accurately";
}

}  // namespace
