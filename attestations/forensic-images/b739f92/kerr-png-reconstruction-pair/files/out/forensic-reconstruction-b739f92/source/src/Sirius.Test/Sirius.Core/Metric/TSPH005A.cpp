// TSPH005A.cpp - Kerr-Schild Metric Tests
// Tests: Kerr radius, null vector, asymptotic flatness, Schwarzschild limit.

#define _USE_MATH_DEFINES
#include <cmath>
#include <gtest/gtest.h>

// Test parameters
constexpr double kEpsilon = 1e-8;
constexpr double kPI = 3.14159265358979323846;

// =============================================================================
// Test Fixture
// =============================================================================

class KerrSchildTests : public ::testing::Test {
protected:
    // Kerr radius from Cartesian coordinates
    // r^2 = (R^2 - a^2 + sqrt((R^2-a^2)^2 + 4*a^2*z^2)) / 2
    double kerr_radius(double x, double y, double z, double a) {
        double R2 = x*x + y*y + z*z;
        double a2 = a * a;
        double Rm2 = x*x + y*y - z*z;
        double disc = a2*a2 - 2.0*a2*Rm2 + R2*R2;
        double r2 = (-a2 + std::sqrt(disc) + R2) / 2.0;
        return std::sqrt(std::max(r2, 1e-10));
    }
    
    // Null vector l^mu
    void null_vector(double x, double y, double z, double a, double l[4]) {
        double r = kerr_radius(x, y, z, a);
        double r2 = r * r;
        double a2 = a * a;
        double denom = r2 + a2;
        
        l[0] = 1.0;
        l[1] = (r*x + a*y) / denom;
        l[2] = (r*y - a*x) / denom;
        l[3] = z / r;
    }
    
    // Kerr-Schild function f
    double ks_f(double x, double y, double z, double M, double a) {
        double r = kerr_radius(x, y, z, a);
        double r2 = r * r;
        double rs = 2.0 * M;
        return rs * r2 * r / (r2*r2 + a*a*z*z);
    }
};

// =============================================================================
// Kerr Radius Tests
// =============================================================================

TEST_F(KerrSchildTests, KerrRadiusEquatorialPlane) {
    // In equatorial plane (z=0), r = sqrt(x^2 + y^2) for a=0
    double x = 3.0, y = 4.0, z = 0.0;
    double a = 0.0;  // Schwarzschild limit
    
    double r_expected = std::sqrt(x*x + y*y);  // = 5.0
    double r_computed = kerr_radius(x, y, z, a);
    
    EXPECT_NEAR(r_computed, r_expected, kEpsilon)
        << "Kerr radius should equal Euclidean distance in equatorial plane for a=0";
}

TEST_F(KerrSchildTests, KerrRadiusOnAxis) {
    // On z-axis (x=y=0)
    double x = 0.0, y = 0.0, z = 5.0;
    double a = 0.5;
    
    double r_computed = kerr_radius(x, y, z, a);
    
    // Should satisfy r^2 + a^2 = z^2 + a^2 in limit where x=y=0
    EXPECT_GT(r_computed, 0.0) << "Kerr radius should be positive on axis";
    EXPECT_NEAR(r_computed, std::abs(z), 0.5) 
        << "On axis, Kerr radius approximately equals |z| for small a";
}

TEST_F(KerrSchildTests, KerrRadiusSymmetry) {
    // Axisymmetry: r(x,y,z) = r(-x,-y,z) = r(x,y,-z)
    double a = 0.7;
    
    double r1 = kerr_radius(3.0, 4.0, 2.0, a);
    double r2 = kerr_radius(-3.0, -4.0, 2.0, a);
    double r3 = kerr_radius(3.0, 4.0, -2.0, a);
    
    EXPECT_NEAR(r1, r2, kEpsilon) << "Should be symmetric under (x,y)->(-x,-y)";
    EXPECT_NEAR(r1, r3, kEpsilon) << "Should be symmetric under z->-z";
}

// =============================================================================
// Null Vector Tests
// =============================================================================

TEST_F(KerrSchildTests, NullVectorIsNull) {
    // l^u should satisfy g_uv l^u l^v = 0 with respect to Minkowski metric
    // Since g = eta + f*l*l, and l is null w.r.t. eta, l*l gives correct zero
    double x = 5.0, y = 0.0, z = 0.0;
    double a = 0.5;
    
    double l[4];
    null_vector(x, y, z, a, l);
    
    // Check l is null in Minkowski: eta_uv l^u l^v = -l0^2 + l1^2 + l2^2 + l3^2 = 0
    double minkowski_norm = -l[0]*l[0] + l[1]*l[1] + l[2]*l[2] + l[3]*l[3];
    
    EXPECT_NEAR(minkowski_norm, 0.0, kEpsilon)
        << "Null vector should be null with respect to Minkowski metric";
}

TEST_F(KerrSchildTests, NullVectorTimeComponent) {
    // l^0 = 1 for all points
    double l[4];
    
    null_vector(1.0, 2.0, 3.0, 0.5, l);
    EXPECT_EQ(l[0], 1.0) << "Temporal component of null vector is always 1";
    
    null_vector(10.0, 0.0, 0.0, 0.9, l);
    EXPECT_EQ(l[0], 1.0) << "Temporal component of null vector is always 1";
}

// =============================================================================
// Asymptotic Flatness Tests
// =============================================================================

TEST_F(KerrSchildTests, AsymptoticFlatness) {
    // At large r, f -> 0, so g -> eta (Minkowski)
    double M = 1.0, a = 0.5;
    double r_large = 1000.0;  // Far from black hole
    
    double f = ks_f(r_large, 0.0, 0.0, M, a);
    
    EXPECT_LT(std::abs(f), 1e-2)
        << "Kerr-Schild function f -> 0 at large r";
}

TEST_F(KerrSchildTests, SchwarzschildLimit) {
    // For a=0, should reduce to Schwarzschild
    double M = 1.0, a = 0.0;
    double x = 4.0, y = 0.0, z = 0.0;  // r = 4
    
    double f = ks_f(x, y, z, M, a);
    
    // For Schwarzschild, f = rs/r = 2M/r = 0.5
    double r = kerr_radius(x, y, z, a);
    double f_expected = 2.0 * M / r;
    
    EXPECT_NEAR(f, f_expected, kEpsilon)
        << "At a=0, Kerr-Schild function should match Schwarzschild rs/r";
}

// =============================================================================
// Metric Component Tests
// =============================================================================

TEST_F(KerrSchildTests, MetricDiagonalPositiveDefiniteSpatial) {
    // g_ii = 1 + f*l_i^2 > 0 for spatial components
    double M = 1.0, a = 0.5;
    double x = 4.0, y = 3.0, z = 2.0;
    
    double f = ks_f(x, y, z, M, a);
    double l[4];
    null_vector(x, y, z, a, l);
    
    double g_xx = 1.0 + f * l[1] * l[1];
    double g_yy = 1.0 + f * l[2] * l[2];
    double g_zz = 1.0 + f * l[3] * l[3];
    
    EXPECT_GT(g_xx, 0.0) << "g_xx should be positive";
    EXPECT_GT(g_yy, 0.0) << "g_yy should be positive";
    EXPECT_GT(g_zz, 0.0) << "g_zz should be positive";
}

TEST_F(KerrSchildTests, MetricTimeComponentSign) {
    // g_tt = -1 + f*l_0^2 = -1 + f should be negative outside horizon
    double M = 1.0, a = 0.5;
    double r_horizon = M + std::sqrt(M*M - a*a);
    
    // Test outside horizon
    double r_test = r_horizon * 2.0;
    double x = r_test, y = 0.0, z = 0.0;
    
    double f = ks_f(x, y, z, M, a);
    double g_tt = -1.0 + f;  // l[0] = 1
    
    EXPECT_LT(g_tt, 0.0) << "g_tt should be negative outside horizon";
}

// =============================================================================
// Coordinate Transform Tests
// =============================================================================

TEST_F(KerrSchildTests, OblateSpheroidalRelation) {
    // The oblate spheroidal coordinate relation is:
    // x^2 + y^2 = (r^2 + a^2) * sin^2(theta)
    // z = r * cos(theta)
    
    double a = 0.5;
    double r_true = 5.0;
    double theta = kPI / 4.0;  // 45 degrees
    
    // Construct Cartesian from known r, theta
    double x = std::sqrt(r_true*r_true + a*a) * std::sin(theta);
    double y = 0.0;
    double z = r_true * std::cos(theta);
    
    double r_computed = kerr_radius(x, y, z, a);
    
    EXPECT_NEAR(r_computed, r_true, 1e-3)
        << "Kerr radius should match for oblate spheroidal construction";
}

// main() is provided by GTest::gtest_main
