// TSPH014A.cpp - Double-Precision Metric Tests
// Tests for PHMT000B/PHMT100B double-precision metric infrastructure
// Verifies: metric inverse identity, Christoffel symmetry, Hamiltonian conservation

#include <gtest/gtest.h>
#include <cmath>
#include <MTTP001A.h>
#include <PHMT000B.h>
#include <PHMT100B.h>

using namespace sirius::math;
using namespace sirius::physics;

namespace {

// Test fixture for double-precision Kerr metric
class KerrMetricDTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create metrics with various spin parameters
        schwarzschild = std::make_unique<KerrMetricD>(1.0, 0.0);
        kerr_moderate = std::make_unique<KerrMetricD>(1.0, 0.5);
        kerr_extreme = std::make_unique<KerrMetricD>(1.0, 0.999);
    }
    
    std::unique_ptr<KerrMetricD> schwarzschild;
    std::unique_ptr<KerrMetricD> kerr_moderate;
    std::unique_ptr<KerrMetricD> kerr_extreme;
};

//==============================================================================
// Test: Metric Inverse Identity
// Verifies: g^μα g_αν = δ^μ_ν to tolerance < 10^-14
//==============================================================================

TEST_F(KerrMetricDTest, MetricInverseIdentity_Schwarzschild) {
    const double tolerance = 1e-14;
    
    // Test at various positions
    std::vector<Vec4d> positions = {
        Vec4d(0, 10.0, M_PI/2, 0),      // Equatorial plane
        Vec4d(0, 5.0, M_PI/4, 0),       // Off-equator
        Vec4d(0, 3.0, M_PI/3, M_PI/2),  // Near photon sphere
        Vec4d(0, 100.0, M_PI/2, 0),     // Far field
    };
    
    for (const auto& x : positions) {
        double g[4][4], g_inv[4][4];
        schwarzschild->evaluate(x, g, g_inv);
        
        // Compute product g^μα g_αν
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                double sum = 0;
                for (int alpha = 0; alpha < 4; ++alpha) {
                    sum += g_inv[mu][alpha] * g[alpha][nu];
                }
                
                double expected = (mu == nu) ? 1.0 : 0.0;
                EXPECT_NEAR(sum, expected, tolerance)
                    << "Failed at r=" << x.r << ", theta=" << x.theta
                    << ", mu=" << mu << ", nu=" << nu;
            }
        }
    }
}

TEST_F(KerrMetricDTest, MetricInverseIdentity_Kerr) {
    const double tolerance = 1e-14;
    
    std::vector<Vec4d> positions = {
        Vec4d(0, 10.0, M_PI/2, 0),
        Vec4d(0, 5.0, M_PI/4, 0),
        Vec4d(0, 3.0, M_PI/3, M_PI/2),
    };
    
    for (const auto& x : positions) {
        double g[4][4], g_inv[4][4];
        kerr_moderate->evaluate(x, g, g_inv);
        
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                double sum = 0;
                for (int alpha = 0; alpha < 4; ++alpha) {
                    sum += g_inv[mu][alpha] * g[alpha][nu];
                }
                
                double expected = (mu == nu) ? 1.0 : 0.0;
                EXPECT_NEAR(sum, expected, tolerance)
                    << "Failed for a=0.5 at r=" << x.r << ", theta=" << x.theta;
            }
        }
    }
}

//==============================================================================
// Test: Christoffel Symmetry
// Verifies: Γ^μ_νρ = Γ^μ_ρν (symmetric in lower indices)
//==============================================================================

TEST_F(KerrMetricDTest, ChristoffelSymmetry) {
    const double tolerance = 1e-12;
    
    Vec4d x(0, 6.0, M_PI/3, 0);  // Test position
    double Gamma[4][4][4];
    kerr_moderate->christoffel(x, Gamma);
    
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            for (int rho = 0; rho < 4; ++rho) {
                EXPECT_NEAR(Gamma[mu][nu][rho], Gamma[mu][rho][nu], tolerance)
                    << "Christoffel asymmetry at mu=" << mu 
                    << ", nu=" << nu << ", rho=" << rho;
            }
        }
    }
}

//==============================================================================
// Test: Hamiltonian for Null Geodesic
// Verifies: H = (1/2) g^μν p_μ p_ν ≈ 0 for properly initialised null ray
//==============================================================================

TEST_F(KerrMetricDTest, HamiltonianNullGeodesic) {
    const double tolerance = 1e-10;
    
    // Set up a null ray at r=10M pointing inward
    Vec4d x(0, 10.0, M_PI/2, 0);
    
    // Get metric at position
    double g[4][4], g_inv[4][4];
    kerr_moderate->evaluate(x, g, g_inv);
    
    // Set up covariant momentum for ingoing null ray
    // For null: g^μν p_μ p_ν = 0
    // At equator in Kerr: we need g^tt p_t^2 + 2 g^tφ p_t p_φ + g^rr p_r^2 + g^φφ p_φ^2 = 0
    
    Vec4d p;
    p.t = -1.0;  // E = 1
    p.phi = 0.2; // Some angular momentum
    p.theta = 0;
    
    // Solve for p_r from null condition
    // g^tt p_t^2 + 2 g^tφ p_t p_φ + g^rr p_r^2 + g^φφ p_φ^2 = 0
    double A = g_inv[1][1];
    double C = g_inv[0][0]*p.t*p.t + 2*g_inv[0][3]*p.t*p.phi + g_inv[3][3]*p.phi*p.phi;
    
    if (A > 0 && C < 0) {
        p.r = std::sqrt(-C / A);  // Ingoing ray
    } else {
        p.r = 0;  // Fallback
    }
    
    double H = kerr_moderate->hamiltonian(x, p);
    EXPECT_NEAR(H, 0.0, tolerance) << "Hamiltonian for null ray should be ~0";
}

//==============================================================================
// Test: Horizon Radius
// Verifies: r_+ = M + √(M² - a²)
//==============================================================================

TEST_F(KerrMetricDTest, HorizonRadius) {
    const double tolerance = 1e-12;
    
    // Schwarzschild: r_+ = 2M
    EXPECT_NEAR(schwarzschild->horizonRadius(), 2.0, tolerance);
    
    // Kerr a=0.5: r_+ = 1 + √(1 - 0.25) = 1 + √0.75 ≈ 1.866
    double expected = 1.0 + std::sqrt(0.75);
    EXPECT_NEAR(kerr_moderate->horizonRadius(), expected, tolerance);
    
    // Kerr a=0.999: r_+ → 1 as a → 1
    double expected_extreme = 1.0 + std::sqrt(1 - 0.999*0.999);
    EXPECT_NEAR(kerr_extreme->horizonRadius(), expected_extreme, tolerance);
}

//==============================================================================
// Test: ISCO Radius
// Verifies: Schwarzschild ISCO = 6M, decreases with spin
//==============================================================================

TEST_F(KerrMetricDTest, ISCORadius) {
    const double tolerance = 1e-6;
    
    // Schwarzschild: ISCO = 6M
    EXPECT_NEAR(schwarzschild->iscoRadius(), 6.0, tolerance);
    
    // Kerr with spin: ISCO < 6M (prograde)
    double isco_moderate = kerr_moderate->iscoRadius();
    EXPECT_LT(isco_moderate, 6.0);
    EXPECT_GT(isco_moderate, 1.0);  // Must be outside horizon
    
    // Extreme Kerr: ISCO → M
    double isco_extreme = kerr_extreme->iscoRadius();
    EXPECT_LT(isco_extreme, 2.0);
    EXPECT_GT(isco_extreme, 1.0);
}

//==============================================================================
// Test: Photon Sphere Radius
// Verifies: Schwarzschild photon sphere = 3M
//==============================================================================

TEST_F(KerrMetricDTest, PhotonSphereRadius) {
    const double tolerance = 1e-6;
    
    // Schwarzschild: photon sphere at 3M
    EXPECT_NEAR(schwarzschild->photonSphereRadius(), 3.0, tolerance);
    
    // Kerr: photon sphere < 3M for prograde
    EXPECT_LT(kerr_moderate->photonSphereRadius(), 3.0);
    EXPECT_LT(kerr_extreme->photonSphereRadius(), 2.0);
}

//==============================================================================
// Test: Vec4d Arithmetic
//==============================================================================

TEST(Vec4dTest, Arithmetic) {
    Vec4d a(1, 2, 3, 4);
    Vec4d b(0.5, 1, 1.5, 2);
    
    Vec4d sum = a + b;
    EXPECT_DOUBLE_EQ(sum.t, 1.5);
    EXPECT_DOUBLE_EQ(sum.r, 3.0);
    EXPECT_DOUBLE_EQ(sum.theta, 4.5);
    EXPECT_DOUBLE_EQ(sum.phi, 6.0);
    
    Vec4d diff = a - b;
    EXPECT_DOUBLE_EQ(diff.t, 0.5);
    EXPECT_DOUBLE_EQ(diff.r, 1.0);
    
    Vec4d scaled = a * 2.0;
    EXPECT_DOUBLE_EQ(scaled.t, 2.0);
    EXPECT_DOUBLE_EQ(scaled.r, 4.0);
}

TEST(Vec4dTest, IndexedAccess) {
    Vec4d v(1, 2, 3, 4);
    EXPECT_DOUBLE_EQ(v[0], 1.0);  // t
    EXPECT_DOUBLE_EQ(v[1], 2.0);  // r
    EXPECT_DOUBLE_EQ(v[2], 3.0);  // theta
    EXPECT_DOUBLE_EQ(v[3], 4.0);  // phi
}

//==============================================================================
// Test: Mat4d Operations
//==============================================================================

TEST(Mat4dTest, Identity) {
    Mat4d I = Mat4d::identity();
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_DOUBLE_EQ(I(i, j), expected);
        }
    }
}

TEST(Mat4dTest, MatrixVectorMultiply) {
    Mat4d I = Mat4d::identity();
    Vec4d v(1, 2, 3, 4);
    Vec4d result = I * v;
    
    EXPECT_DOUBLE_EQ(result.t, 1.0);
    EXPECT_DOUBLE_EQ(result.r, 2.0);
    EXPECT_DOUBLE_EQ(result.theta, 3.0);
    EXPECT_DOUBLE_EQ(result.phi, 4.0);
}

TEST(Mat4dTest, Determinant) {
    Mat4d I = Mat4d::identity();
    EXPECT_DOUBLE_EQ(I.determinant(), 1.0);
    
    Mat4d scaled = I * 2.0;
    EXPECT_DOUBLE_EQ(scaled.determinant(), 16.0);  // 2^4
}

} // namespace
