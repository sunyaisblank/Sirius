// TSPH015A.cpp - Coordinate Transformation Tests
// Component ID: TSPH015A
// Tests for PHCT002A: Boyer-Lindquist ↔ Kerr-Schild Cartesian coordinate transforms
//
// Specification requirements:
//   - Round-trip deviation < 1e-12 (see plan Phase 5.1)
//   - Jacobian determinant = r² sin(θ) for BL → Cart
//
// Tests: PHCT002A.h coordinate transformation functions

#include "sirius/core/coordinates.h"

#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

using namespace sirius::core::coordinates;

namespace {

constexpr double kEpsilon = 1e-12;
constexpr double kLooseEpsilon = 1e-10;  // For floating-point accumulation

// =============================================================================
// Test Fixture
// =============================================================================

class CoordinateTransformTests : public ::testing::Test {
  protected:
    // Standard test points in BL coordinates
    Vec4Bl point_equator{0.0, 10.0, std::numbers::pi / 2.0,
                         std::numbers::pi / 4.0};                    // Equatorial plane
    Vec4Bl point_pole_near{0.0, 10.0, 0.1, std::numbers::pi / 4.0};  // Near north pole
    Vec4Bl point_far{0.0, 100.0, std::numbers::pi / 3.0, 2.0};       // Far field
    Vec4Bl point_close{0.0, 3.0, std::numbers::pi / 4.0, -std::numbers::pi / 2.0};  // Near horizon
};

// =============================================================================
// Round-Trip Tests (BL → Cartesian → BL)
// =============================================================================

TEST_F(CoordinateTransformTests, RoundTripBLCartesian_Equator) {
    // Round-trip at equator (θ = π/2)
    double deviation = ValidateRoundTrip(point_equator, 0.0);
    EXPECT_LT(deviation, kEpsilon) << "Round-trip deviation at equator should be < 1e-12";
}

TEST_F(CoordinateTransformTests, RoundTripBLCartesian_NearPole) {
    // Round-trip near pole (potential coordinate singularity)
    double deviation = ValidateRoundTrip(point_pole_near, 0.0);
    EXPECT_LT(deviation, kLooseEpsilon) << "Round-trip deviation near pole should be < 1e-10";
}

TEST_F(CoordinateTransformTests, RoundTripBLCartesian_FarField) {
    // Round-trip in far field
    double deviation = ValidateRoundTrip(point_far, 0.0);
    EXPECT_LT(deviation, kEpsilon) << "Round-trip deviation in far field should be < 1e-12";
}

TEST_F(CoordinateTransformTests, RoundTripBLCartesian_NearHorizon) {
    // Round-trip near horizon (strong curvature region)
    double deviation = ValidateRoundTrip(point_close, 0.0);
    EXPECT_LT(deviation, kEpsilon) << "Round-trip deviation near horizon should be < 1e-12";
}

// =============================================================================
// Kerr Round-Trip Tests (with spin parameter)
// =============================================================================

TEST_F(CoordinateTransformTests, RoundTripKerrSchild_LowSpin) {
    double a = 0.1;
    double deviation = ValidateRoundTrip(point_equator, a);
    EXPECT_LT(deviation, kEpsilon);
}

TEST_F(CoordinateTransformTests, RoundTripKerrSchild_ModerateSpin) {
    // Kerr with moderate spin a = 0.5
    double a = 0.5;
    double deviation = ValidateRoundTrip(point_equator, a);
    EXPECT_LT(deviation, kEpsilon);
}

TEST_F(CoordinateTransformTests, RoundTripKerrSchild_HighSpin) {
    // Kerr with high spin a = 0.9 (near extremal)
    double a = 0.9;
    double deviation = ValidateRoundTrip(point_equator, a);
    EXPECT_LT(deviation, kEpsilon);
}

TEST_F(CoordinateTransformTests, KerrDiskRadiusIsSpheroidalNotCylindrical) {
    constexpr double spin = 0.9;
    const Vec4Cart cart = BlToKerrSchildCart(point_equator, spin);
    const double cylindrical_radius = std::hypot(cart.x, cart.y);

    EXPECT_NEAR(KerrSchildRadius(cart, spin), point_equator.r, kEpsilon);
    EXPECT_NEAR(cylindrical_radius, std::sqrt(point_equator.r * point_equator.r + spin * spin),
                kEpsilon);
    EXPECT_GT(cylindrical_radius - KerrSchildRadius(cart, spin), 0.0);
}

TEST_F(CoordinateTransformTests, KerrRadiusPreservesExactZeroAndScaleCovariance) {
    constexpr double spin = 0.9;
    constexpr double scale = 1.0e-140;
    const Vec4Cart disk{0.0, 0.5 * spin, 0.0, 0.0};
    const Vec4Cart unit{0.0, 2.0, 1.0, 0.5};
    const Vec4Cart scaled{0.0, 2.0 * scale, scale, 0.5 * scale};

    EXPECT_DOUBLE_EQ(KerrSchildRadius(disk, spin), 0.0);
    EXPECT_FALSE(TryKerrSchildRadiusDifferential(disk, spin).has_value());

    const auto unit_differential = TryKerrSchildRadiusDifferential(unit, 0.5);
    const auto scaled_differential = TryKerrSchildRadiusDifferential(scaled, 0.5 * scale);
    ASSERT_TRUE(unit_differential);
    ASSERT_TRUE(scaled_differential);
    EXPECT_NEAR(scaled_differential->radius / scale, unit_differential->radius, 3.0e-15);
    EXPECT_NEAR(scaled_differential->dx, unit_differential->dx, 3.0e-15);
    EXPECT_NEAR(scaled_differential->dy, unit_differential->dy, 3.0e-15);
    EXPECT_NEAR(scaled_differential->dz, unit_differential->dz, 3.0e-15);

    const Vec4Cart unit_vector{0.2, 0.3, -0.1, 0.05};
    const Vec4Cart scaled_vector{0.2 * scale, 0.3 * scale, -0.1 * scale, 0.05 * scale};
    const Vec4Bl unit_bl = TransformVectorKerrSchildCartToBl(unit_vector, unit, 1.0, 0.5);
    const Vec4Bl scaled_bl =
        TransformVectorKerrSchildCartToBl(scaled_vector, scaled, scale, 0.5 * scale);
    EXPECT_NEAR(scaled_bl.t / scale, unit_bl.t, 3.0e-14);
    EXPECT_NEAR(scaled_bl.r / scale, unit_bl.r, 3.0e-14);
    EXPECT_NEAR(scaled_bl.theta, unit_bl.theta, 3.0e-14);
    EXPECT_NEAR(scaled_bl.phi, unit_bl.phi, 3.0e-14);

    EXPECT_DEATH((void)KerrSchildCartToBl(disk, spin), "precondition.*enforced, terminating");
    const Vec4Cart axis{0.0, 0.0, 0.0, 2.0};
    EXPECT_DEATH((void)KerrSchildCartToBl(axis, 0.5), "precondition.*enforced, terminating");
    EXPECT_DEATH((void)TransformVectorKerrSchildCartToBl(unit_vector, axis, 1.0, 0.5),
                 "precondition.*enforced, terminating");
}

TEST_F(CoordinateTransformTests, RoundTripKerrSchild_NearPole) {
    // Kerr near pole (tests oblate spheroidal effects)
    double a = 0.5;
    double deviation = ValidateRoundTrip(point_pole_near, a);
    EXPECT_LT(deviation, kEpsilon);
}

// =============================================================================
// Jacobian Tests
// =============================================================================

TEST_F(CoordinateTransformTests, JacobianDeterminant_Equator) {
    // At equator: det(J) = r² sin(θ) = r² × 1 = r²
    Jacobian4x4 J = JacobianBlToCartesian(point_equator);
    double det = JacobianDeterminant(J);
    double expected = point_equator.r * point_equator.r * std::sin(point_equator.theta);

    EXPECT_NEAR(det, expected, 1e-10) << "Jacobian determinant should equal r² sin(θ)";
}

TEST_F(CoordinateTransformTests, JacobianDeterminant_MidLatitude) {
    // At θ = π/3: det(J) = r² sin(π/3) = r² × √3/2
    Vec4Bl mid_lat{0.0, 10.0, std::numbers::pi / 3.0, 0.0};
    Jacobian4x4 J = JacobianBlToCartesian(mid_lat);
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
    Vec4Bl v_orig{1.0, 0.5, 0.25, 0.1};  // Some arbitrary 4-vector

    Vec4Cart v_cart = TransformVectorBlToCart(v_orig, point_equator);
    Vec4Cart pos_cart = BlToCartesian(point_equator);
    Vec4Bl v_recovered = TransformVectorCartToBl(v_cart, pos_cart);

    EXPECT_NEAR(v_orig.t, v_recovered.t, kLooseEpsilon);
    EXPECT_NEAR(v_orig.r, v_recovered.r, kLooseEpsilon);
    EXPECT_NEAR(v_orig.theta, v_recovered.theta, kLooseEpsilon);

    // Phi may wrap, so check with modular arithmetic
    double dphi = v_orig.phi - v_recovered.phi;
    while (dphi > std::numbers::pi) dphi -= 2 * std::numbers::pi;
    while (dphi < -std::numbers::pi) dphi += 2 * std::numbers::pi;
    EXPECT_NEAR(dphi, 0.0, kLooseEpsilon);
}

// =============================================================================
// Boundary Condition Tests
// =============================================================================

TEST_F(CoordinateTransformTests, OriginHandling) {
    // Test behavior at r → 0 (should be handled gracefully)
    Vec4Bl origin{0.0, 0.0, std::numbers::pi / 2.0, 0.0};
    Vec4Cart cart = BlToCartesian(origin);

    EXPECT_EQ(cart.x, 0.0);
    EXPECT_EQ(cart.y, 0.0);
    EXPECT_EQ(cart.z, 0.0);
    EXPECT_EQ(cart.t, 0.0);

    Vec4Bl bl;
    const Vec4Bl& const_bl = bl;
    Vec4Cart ks;
    const Vec4Cart& const_ks = ks;
    EXPECT_DEATH(static_cast<void>(bl[-1]), "violated");
    EXPECT_DEATH(static_cast<void>(const_bl[4]), "violated");
    EXPECT_DEATH(static_cast<void>(ks[-1]), "violated");
    EXPECT_DEATH(static_cast<void>(const_ks[4]), "violated");
}

TEST_F(CoordinateTransformTests, PhiWrapping) {
    // Test phi wrapping at ±π boundary
    Vec4Bl pos1{0.0, 10.0, std::numbers::pi / 2.0, std::numbers::pi - 0.01};
    Vec4Bl pos2{0.0, 10.0, std::numbers::pi / 2.0, -std::numbers::pi + 0.01};

    Vec4Cart cart1 = BlToCartesian(pos1);
    Vec4Cart cart2 = BlToCartesian(pos2);

    // Points near ±π in phi should be close in Cartesian space
    double dx = cart1.x - cart2.x;
    double dy = cart1.y - cart2.y;
    double dist = std::sqrt(dx * dx + dy * dy);

    EXPECT_LT(dist, 0.5) << "Points near phi = ±π should be close in Cartesian space";
}

TEST_F(CoordinateTransformTests, NorthPoleCoordinates) {
    // At north pole (θ = 0), x = y = 0, z = r
    Vec4Bl north_pole{0.0, 10.0, 0.0, 0.0};
    Vec4Cart cart = BlToCartesian(north_pole);

    EXPECT_NEAR(cart.x, 0.0, 1e-14);
    EXPECT_NEAR(cart.y, 0.0, 1e-14);
    EXPECT_NEAR(cart.z, 10.0, 1e-14);
}

TEST_F(CoordinateTransformTests, SouthPoleCoordinates) {
    // At south pole (θ = π), x = y = 0, z = -r
    Vec4Bl south_pole{0.0, 10.0, std::numbers::pi, 0.0};
    Vec4Cart cart = BlToCartesian(south_pole);

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
    double theta = std::numbers::pi / 3.0;

    Vec4Bl pos{0.0, r, theta, 0.0};
    Vec4Cart cart = BlToKerrSchildCart(pos, a);

    double expected_rho = std::sqrt(r * r + a * a);
    double expected_xy_mag = expected_rho * std::sin(theta);
    double actual_xy_mag = std::sqrt(cart.x * cart.x + cart.y * cart.y);

    EXPECT_NEAR(actual_xy_mag, expected_xy_mag, 1e-12)
        << "Kerr oblate spheroidal: |x,y| should equal √(r²+a²) sin(θ)";

    EXPECT_NEAR(cart.z, r * std::cos(theta), 1e-12) << "Kerr: z should equal r cos(θ)";
}

TEST_F(CoordinateTransformTests, KerrSolveR_Accuracy) {
    // Test that inverse transformation correctly solves for r
    // from x² + y² + z² = r² + a² (1 - z²/r²)
    double a = 0.7;
    double r_original = 8.0;
    double theta = std::numbers::pi / 4.0;
    double phi = 1.5;

    Vec4Bl original{0.0, r_original, theta, phi};
    Vec4Cart cart = BlToKerrSchildCart(original, a);
    Vec4Bl recovered = KerrSchildCartToBl(cart, a);

    EXPECT_NEAR(original.r, recovered.r, 1e-10) << "Kerr inverse should recover r accurately";
    EXPECT_NEAR(original.theta, recovered.theta, 1e-10)
        << "Kerr inverse should recover θ accurately";
}

}  // namespace
