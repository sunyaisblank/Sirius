// TSCT001A.cpp - Coordinate Transform Interface Tests
// Component ID: TSCT001A (Test/Coordinate/Transforms)
// Tests: PHCT001A.h (Coordinates namespace)

#include <gtest/gtest.h>
#include <cmath>
#include <PHCT001A.h>

namespace sirius::test {
using namespace Sirius;

constexpr double kEps = 1e-10;
constexpr double kPi = 3.14159265358979323846;

// =============================================================================
// Cartesian ↔ Spherical Round-Trip
// =============================================================================

class CoordUtilityTests : public ::testing::Test {};

TEST_F(CoordUtilityTests, RoundTripCartesianSpherical) {
    double test_points[][3] = {
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0},
        {1.0, 1.0, 1.0}, {3.0, -2.0, 5.0}, {0.1, 0.2, 0.3}
    };

    for (auto& p : test_points) {
        double r, theta, phi;
        Sirius::Coordinates::cartesianToSpherical(p[0], p[1], p[2], r, theta, phi);

        double x2, y2, z2;
        Sirius::Coordinates::sphericalToCartesian(r, theta, phi, x2, y2, z2);

        EXPECT_NEAR(x2, p[0], kEps) << "x mismatch for (" << p[0] << "," << p[1] << "," << p[2] << ")";
        EXPECT_NEAR(y2, p[1], kEps) << "y mismatch";
        EXPECT_NEAR(z2, p[2], kEps) << "z mismatch";
    }
}

TEST_F(CoordUtilityTests, RoundTripSphericalCartesian) {
    double test_points[][3] = {
        {5.0, kPi / 4.0, kPi / 3.0},
        {1.0, kPi / 2.0, 0.0},
        {10.0, 0.1, 2.5},
        {0.5, kPi - 0.1, -1.0}
    };

    for (auto& p : test_points) {
        double x, y, z;
        Sirius::Coordinates::sphericalToCartesian(p[0], p[1], p[2], x, y, z);

        double r2, theta2, phi2;
        Sirius::Coordinates::cartesianToSpherical(x, y, z, r2, theta2, phi2);

        EXPECT_NEAR(r2, p[0], kEps);
        EXPECT_NEAR(theta2, p[1], kEps);
        EXPECT_NEAR(phi2, p[2], kEps);
    }
}

TEST_F(CoordUtilityTests, OriginHandling) {
    double r, theta, phi;
    Sirius::Coordinates::cartesianToSpherical(0.0, 0.0, 0.0, r, theta, phi);
    // r should be clamped to small positive (1e-10)
    EXPECT_GE(r, 0.0);
    EXPECT_FALSE(std::isnan(theta));
    EXPECT_FALSE(std::isnan(phi));
}

TEST_F(CoordUtilityTests, NorthPole) {
    // (0, 0, 1) → (1, 0, 0) i.e. r=1, theta=0, phi=0
    double r, theta, phi;
    Sirius::Coordinates::cartesianToSpherical(0.0, 0.0, 1.0, r, theta, phi);
    EXPECT_NEAR(r, 1.0, kEps);
    EXPECT_NEAR(theta, 0.0, kEps);
}

TEST_F(CoordUtilityTests, SouthPole) {
    // (0, 0, -1) → (1, pi, 0)
    double r, theta, phi;
    Sirius::Coordinates::cartesianToSpherical(0.0, 0.0, -1.0, r, theta, phi);
    EXPECT_NEAR(r, 1.0, kEps);
    EXPECT_NEAR(theta, kPi, kEps);
}

TEST_F(CoordUtilityTests, KnownValueXAxis) {
    // (1, 0, 0) → (1, pi/2, 0)
    double r, theta, phi;
    Sirius::Coordinates::cartesianToSpherical(1.0, 0.0, 0.0, r, theta, phi);
    EXPECT_NEAR(r, 1.0, kEps);
    EXPECT_NEAR(theta, kPi / 2.0, kEps);
    EXPECT_NEAR(phi, 0.0, kEps);
}

TEST_F(CoordUtilityTests, KnownValueZAxis) {
    // (0, 0, 1) → (1, 0, 0)
    double r, theta, phi;
    Sirius::Coordinates::cartesianToSpherical(0.0, 0.0, 1.0, r, theta, phi);
    EXPECT_NEAR(r, 1.0, kEps);
    EXPECT_NEAR(theta, 0.0, kEps);
}

// =============================================================================
// 4-Vector Versions
// =============================================================================

TEST_F(CoordUtilityTests, FourVectorTimePassthrough) {
    double cart[4] = {42.0, 1.0, 2.0, 3.0};
    double sph[4];
    Sirius::Coordinates::cartesianToSpherical4(cart, sph);
    EXPECT_DOUBLE_EQ(sph[0], 42.0);

    double cart2[4];
    Sirius::Coordinates::sphericalToCartesian4(sph, cart2);
    EXPECT_DOUBLE_EQ(cart2[0], 42.0);
    EXPECT_NEAR(cart2[1], cart[1], kEps);
    EXPECT_NEAR(cart2[2], cart[2], kEps);
    EXPECT_NEAR(cart2[3], cart[3], kEps);
}

// =============================================================================
// Velocity Transform Round-Trip
// =============================================================================

TEST_F(CoordUtilityTests, VelocityRoundTrip) {
    double x = 3.0, y = 4.0, z = 5.0;
    double ut = 1.0, ux = 0.5, uy = -0.3, uz = 0.2;

    // Cart → Sph
    double ut_s, ur, utheta, uphi;
    Sirius::Coordinates::velocityCartToSph(x, y, z, ut, ux, uy, uz,
                                           ut_s, ur, utheta, uphi);

    // Get spherical position
    double r, theta, phi;
    Sirius::Coordinates::cartesianToSpherical(x, y, z, r, theta, phi);

    // Sph → Cart
    double ut2, ux2, uy2, uz2;
    Sirius::Coordinates::velocitySphToCart(r, theta, phi, ut_s, ur, utheta, uphi,
                                           ut2, ux2, uy2, uz2);

    EXPECT_NEAR(ut2, ut, 1e-8);
    EXPECT_NEAR(ux2, ux, 1e-8);
    EXPECT_NEAR(uy2, uy, 1e-8);
    EXPECT_NEAR(uz2, uz, 1e-8);
}

TEST_F(CoordUtilityTests, VelocityTimeComponentUnchanged) {
    double ut_s, ur, utheta, uphi;
    Sirius::Coordinates::velocityCartToSph(1.0, 0.0, 0.0, 7.0, 1.0, 0.0, 0.0,
                                           ut_s, ur, utheta, uphi);
    EXPECT_DOUBLE_EQ(ut_s, 7.0);
}

// =============================================================================
// Kerr Radius
// =============================================================================

TEST_F(CoordUtilityTests, KerrRadiusZeroSpin) {
    // For a = 0, Kerr radius = Euclidean radius
    double r = Sirius::Coordinates::computeKerrRadius(3.0, 4.0, 0.0, 0.0);
    EXPECT_NEAR(r, 5.0, kEps);
}

TEST_F(CoordUtilityTests, KerrRadiusPositive) {
    double r = Sirius::Coordinates::computeKerrRadius(1.0, 1.0, 1.0, 0.9);
    EXPECT_GT(r, 0.0);
    EXPECT_FALSE(std::isnan(r));
}

// =============================================================================
// Equirectangular Mapping
// =============================================================================

TEST_F(CoordUtilityTests, EquirectangularBounds) {
    double directions[][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, -1}};
    for (auto& d : directions) {
        double u, v;
        Sirius::Coordinates::directionToEquirectangular(d[0], d[1], d[2], u, v);
        EXPECT_GE(u, 0.0) << "u below 0 for direction (" << d[0] << "," << d[1] << "," << d[2] << ")";
        EXPECT_LE(u, 1.0);
        EXPECT_GE(v, 0.0);
        EXPECT_LE(v, 1.0);
    }
}

TEST_F(CoordUtilityTests, EquirectangularPoles) {
    double u, v;
    // +z direction: theta = 0, v = 0
    Sirius::Coordinates::directionToEquirectangular(0.0, 0.0, 1.0, u, v);
    EXPECT_NEAR(v, 0.0, kEps);

    // -z direction: theta = pi, v = 1
    Sirius::Coordinates::directionToEquirectangular(0.0, 0.0, -1.0, u, v);
    EXPECT_NEAR(v, 1.0, kEps);
}

// =============================================================================
// Accretion Disk Helpers
// =============================================================================

TEST_F(CoordUtilityTests, InAccretionDisk) {
    // Inside disk: z near 0, r_cyl between inner and outer
    EXPECT_TRUE(Sirius::Coordinates::inAccretionDisk(5.0, 0.0, 0.01, 3.0, 10.0, 0.1));
    // Too high
    EXPECT_FALSE(Sirius::Coordinates::inAccretionDisk(5.0, 0.0, 1.0, 3.0, 10.0, 0.1));
    // Too close
    EXPECT_FALSE(Sirius::Coordinates::inAccretionDisk(1.0, 0.0, 0.01, 3.0, 10.0, 0.1));
    // Too far
    EXPECT_FALSE(Sirius::Coordinates::inAccretionDisk(15.0, 0.0, 0.01, 3.0, 10.0, 0.1));
}

TEST_F(CoordUtilityTests, CylindricalRadius) {
    EXPECT_NEAR(Sirius::Coordinates::cylindricalRadius(3.0, 4.0), 5.0, kEps);
}

TEST_F(CoordUtilityTests, AzimuthalAngle) {
    EXPECT_NEAR(Sirius::Coordinates::azimuthalAngle(1.0, 0.0), 0.0, kEps);
    EXPECT_NEAR(Sirius::Coordinates::azimuthalAngle(0.0, 1.0), kPi / 2.0, kEps);
}

// =============================================================================
// SIMD vs Scalar Parity (PHCT003A)
// =============================================================================

} // namespace sirius::test

#include <PHCT003A.h>

namespace sirius::test {
using namespace Sirius;

TEST_F(CoordUtilityTests, SIMDBatchRoundTrip) {
    // Forward transform then inverse should recover original coordinates
    const int N = 16;
    float r[N], theta[N], phi[N];
    float x[N], y[N], z[N];
    float r2[N], theta2[N], phi2[N];

    // Generate test points
    for (int i = 0; i < N; ++i) {
        r[i] = 1.0f + i * 2.0f;
        theta[i] = 0.1f + i * 0.18f;  // 0.1 to ~2.88 (avoid poles)
        phi[i] = i * 0.4f;
    }

    Sirius::Coordinates::SIMD::transformBatch(r, theta, phi, x, y, z, N);
    Sirius::Coordinates::SIMD::inverseTransformBatch(x, y, z, r2, theta2, phi2, N);

    for (int i = 0; i < N; ++i) {
        EXPECT_NEAR(r2[i], r[i], 1e-3f)
            << "r mismatch at index " << i;
        EXPECT_NEAR(theta2[i], theta[i], 1e-3f)
            << "theta mismatch at index " << i;
        // phi wraps to [0, 2π) in inverse
        float expected_phi = std::fmod(phi[i], 6.283185307f);
        if (expected_phi < 0) expected_phi += 6.283185307f;
        EXPECT_NEAR(phi2[i], expected_phi, 1e-3f)
            << "phi mismatch at index " << i;
    }
}

TEST_F(CoordUtilityTests, SIMDScalarParity) {
    // transformBatch (may use AVX2) should match transformBatchScalar
    const int N = 24;  // Not a multiple of 8 to test remainder handling
    float r[N], theta[N], phi[N];
    float x_simd[N], y_simd[N], z_simd[N];
    float x_scalar[N], y_scalar[N], z_scalar[N];

    for (int i = 0; i < N; ++i) {
        r[i] = 0.5f + i * 1.5f;
        theta[i] = 0.2f + i * 0.12f;
        phi[i] = i * 0.25f;
    }

    Sirius::Coordinates::SIMD::transformBatch(r, theta, phi, x_simd, y_simd, z_simd, N);
    Sirius::Coordinates::SIMD::transformBatchScalar(r, theta, phi, x_scalar, y_scalar, z_scalar, N);

    // AVX2 uses fast Taylor sin/cos (~1e-4 accuracy), so tolerance is relaxed
    for (int i = 0; i < N; ++i) {
        EXPECT_NEAR(x_simd[i], x_scalar[i], 1e-2f)
            << "x mismatch at index " << i;
        EXPECT_NEAR(y_simd[i], y_scalar[i], 1e-2f)
            << "y mismatch at index " << i;
        EXPECT_NEAR(z_simd[i], z_scalar[i], 1e-2f)
            << "z mismatch at index " << i;
    }
}

} // namespace sirius::test
