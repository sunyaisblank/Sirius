// Double-precision Kerr metric oracle: inverse identity, Christoffel
// symmetry, Hamiltonian null condition, horizon/ISCO/photon-sphere radii.
// Ported from TSPH014A.cpp; assertions and tolerances unchanged.

#include "sirius/oracle/kerr_boyer_lindquist.h"
#include "sirius/oracle/metric_interface.h"
#include "sirius/oracle/transport_types.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace sirius::oracle;

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
        Vec4d(0, 10.0, M_PI / 2, 0),        // Equatorial plane
        Vec4d(0, 5.0, M_PI / 4, 0),         // Off-equator
        Vec4d(0, 3.0, M_PI / 3, M_PI / 2),  // Near photon sphere
        Vec4d(0, 100.0, M_PI / 2, 0),       // Far field
    };

    for (const auto& x : positions) {
        double g[4][4], g_inv[4][4];
        schwarzschild->Evaluate(x, g, g_inv);

        // Compute product g^μα g_αν
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                double sum = 0;
                for (int alpha = 0; alpha < 4; ++alpha) {
                    sum += g_inv[mu][alpha] * g[alpha][nu];
                }

                double expected = (mu == nu) ? 1.0 : 0.0;
                EXPECT_NEAR(sum, expected, tolerance)
                    << "Failed at r=" << x.r << ", theta=" << x.theta << ", mu=" << mu
                    << ", nu=" << nu;
            }
        }
    }
}

TEST_F(KerrMetricDTest, MetricInverseIdentity_Kerr) {
    const double tolerance = 1e-14;

    std::vector<Vec4d> positions = {
        Vec4d(0, 10.0, M_PI / 2, 0),
        Vec4d(0, 5.0, M_PI / 4, 0),
        Vec4d(0, 3.0, M_PI / 3, M_PI / 2),
    };

    for (const auto& x : positions) {
        double g[4][4], g_inv[4][4];
        kerr_moderate->Evaluate(x, g, g_inv);

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

    Vec4d x(0, 6.0, M_PI / 3, 0);  // Test position
    double Gamma[4][4][4];
    kerr_moderate->Christoffel(x, Gamma);

    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            for (int rho = 0; rho < 4; ++rho) {
                EXPECT_NEAR(Gamma[mu][nu][rho], Gamma[mu][rho][nu], tolerance)
                    << "Christoffel asymmetry at mu=" << mu << ", nu=" << nu << ", rho=" << rho;
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
    Vec4d x(0, 10.0, M_PI / 2, 0);

    // Get metric at position
    double g[4][4], g_inv[4][4];
    kerr_moderate->Evaluate(x, g, g_inv);

    // Set up covariant momentum for ingoing null ray
    // For null: g^μν p_μ p_ν = 0
    // At equator in Kerr: we need g^tt p_t^2 + 2 g^tφ p_t p_φ + g^rr p_r^2 + g^φφ p_φ^2 = 0

    Vec4d p;
    p.t = -1.0;   // E = 1
    p.phi = 0.2;  // Some angular momentum
    p.theta = 0;

    // Solve for p_r from null condition
    // g^tt p_t^2 + 2 g^tφ p_t p_φ + g^rr p_r^2 + g^φφ p_φ^2 = 0
    double A = g_inv[1][1];
    double C =
        g_inv[0][0] * p.t * p.t + 2 * g_inv[0][3] * p.t * p.phi + g_inv[3][3] * p.phi * p.phi;

    if (A > 0 && C < 0) {
        p.r = std::sqrt(-C / A);  // Ingoing ray
    } else {
        p.r = 0;  // Fallback
    }

    double H = kerr_moderate->Hamiltonian(x, p);
    EXPECT_NEAR(H, 0.0, tolerance) << "Hamiltonian for null ray should be ~0";
}

//==============================================================================
// Test: Horizon Radius
// Verifies: r_+ = M + √(M² - a²)
//==============================================================================

TEST_F(KerrMetricDTest, HorizonRadius) {
    const double tolerance = 1e-12;

    // Schwarzschild: r_+ = 2M
    EXPECT_NEAR(schwarzschild->HorizonRadius(), 2.0, tolerance);

    // Kerr a=0.5: r_+ = 1 + √(1 - 0.25) = 1 + √0.75 ≈ 1.866
    double expected = 1.0 + std::sqrt(0.75);
    EXPECT_NEAR(kerr_moderate->HorizonRadius(), expected, tolerance);

    // Kerr a=0.999: r_+ → 1 as a → 1
    double expected_extreme = 1.0 + std::sqrt(1 - 0.999 * 0.999);
    EXPECT_NEAR(kerr_extreme->HorizonRadius(), expected_extreme, tolerance);
}

//==============================================================================
// Test: ISCO Radius
// Verifies: Schwarzschild ISCO = 6M, decreases with spin
//==============================================================================

TEST_F(KerrMetricDTest, ISCORadius) {
    const double tolerance = 1e-6;

    // Schwarzschild: ISCO = 6M
    EXPECT_NEAR(schwarzschild->IscoRadius(), 6.0, tolerance);

    // Kerr with spin: ISCO < 6M (prograde)
    double isco_moderate = kerr_moderate->IscoRadius();
    EXPECT_LT(isco_moderate, 6.0);
    EXPECT_GT(isco_moderate, 1.0);  // Must be outside horizon

    // Extreme Kerr: ISCO → M
    double isco_extreme = kerr_extreme->IscoRadius();
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
    EXPECT_NEAR(schwarzschild->PhotonSphereRadius(), 3.0, tolerance);

    // Kerr: photon sphere < 3M for prograde
    EXPECT_LT(kerr_moderate->PhotonSphereRadius(), 3.0);
    EXPECT_LT(kerr_extreme->PhotonSphereRadius(), 2.0);
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
    Mat4d I = Mat4d::Identity();
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            EXPECT_DOUBLE_EQ(I(i, j), expected);
        }
    }
}

TEST(Mat4dTest, MatrixVectorMultiply) {
    Mat4d I = Mat4d::Identity();
    Vec4d v(1, 2, 3, 4);
    Vec4d result = I * v;

    EXPECT_DOUBLE_EQ(result.t, 1.0);
    EXPECT_DOUBLE_EQ(result.r, 2.0);
    EXPECT_DOUBLE_EQ(result.theta, 3.0);
    EXPECT_DOUBLE_EQ(result.phi, 4.0);
}

TEST(Mat4dTest, Determinant) {
    Mat4d I = Mat4d::Identity();
    EXPECT_DOUBLE_EQ(I.Determinant(), 1.0);

    Mat4d scaled = I * 2.0;
    EXPECT_DOUBLE_EQ(scaled.Determinant(), 16.0);  // 2^4
}

}  // namespace

// The connection must agree with finite differences of the metric it claims
// to derive from; the July-era analytic formulas failed this (Gamma^phi_
// theta_phi lost its leading cot(theta) term, an O(1) error at a = 0), which
// went unnoticed because the integration path uses dHdq. Central differences
// with h = 1e-6 resolve the true connection to ~1e-9; the 1e-6 bar leaves
// three orders of margin while failing the O(1) defect decisively.
TEST_F(KerrMetricDTest, ChristoffelMatchesFiniteDifferencesOfMetric) {
    KerrMetricD metric(1.0, 0.9);
    const Vec4d x(0.0, 5.0, 1.1, 0.3);

    double Gamma[4][4][4];
    metric.Christoffel(x, Gamma);

    const double h = 1e-6;
    double dg[4][4][4];
    for (int sigma = 0; sigma < 4; ++sigma) {
        Vec4d xp = x, xm = x;
        xp[sigma] += h;
        xm[sigma] -= h;
        double gp[4][4], gm[4][4], gi[4][4];
        metric.Evaluate(xp, gp, gi);
        metric.Evaluate(xm, gm, gi);
        for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu)
                dg[sigma][mu][nu] = (gp[mu][nu] - gm[mu][nu]) / (2.0 * h);
    }

    double g[4][4], g_inv[4][4];
    metric.Evaluate(x, g, g_inv);

    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            for (int rho = 0; rho < 4; ++rho) {
                double fd = 0.0;
                for (int sigma = 0; sigma < 4; ++sigma) {
                    fd += g_inv[mu][sigma] *
                          (dg[nu][sigma][rho] + dg[rho][sigma][nu] - dg[sigma][nu][rho]);
                }
                fd *= 0.5;
                EXPECT_NEAR(Gamma[mu][nu][rho], fd,
                            1e-6 * std::max(1.0, std::abs(fd)))
                    << "Gamma^" << mu << "_" << nu << rho
                    << " disagrees with finite differences of the metric";
            }
        }
    }
}
