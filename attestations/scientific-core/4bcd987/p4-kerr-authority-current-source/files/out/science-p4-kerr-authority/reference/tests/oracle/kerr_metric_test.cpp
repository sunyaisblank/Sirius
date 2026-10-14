// Double-precision Kerr metric oracle: inverse identity, Christoffel
// symmetry, Hamiltonian null condition, horizon/ISCO/photon-sphere radii.
// Ported from TSPH014A.cpp; assertions and tolerances unchanged.

#include "sirius/core/coordinates.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/relativistic_transfer.h"
#include "sirius/oracle/kerr_boyer_lindquist.h"
#include "sirius/oracle/metric_interface.h"
#include "sirius/oracle/transport_types.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <numbers>

using namespace sirius::oracle;

namespace {

TEST(KerrZamoTransfer, ArbitraryLatitudeMatchesIndependentBoyerLindquistContraction) {
    constexpr double kMass = 1.0;
    constexpr double kSpin = 0.7;
    constexpr double kRadius = 6.0;
    constexpr double kCosTheta = 0.35;
    constexpr double kEnergy = 1.0;
    constexpr double kAngularMomentum = 1.4;
    constexpr double kObserverFrequency = 0.93;
    const Vec4d event(0.0, kRadius, std::acos(kCosTheta), 0.4);
    KerrMetricD oracle(kMass, kSpin);
    double metric[4][4];
    double inverse[4][4];
    oracle.Evaluate(event, metric, inverse);

    // Derive the locally non-rotating observer only from the independent
    // oracle metric: u_phi=0 gives Omega=-g_tphi/g_phiphi, and u.u=-1 fixes u^t.
    const double omega = -metric[0][3] / metric[3][3];
    const double norm = metric[0][0] + 2.0 * omega * metric[0][3] + omega * omega * metric[3][3];
    ASSERT_LT(norm, 0.0);
    const double time_component = 1.0 / std::sqrt(-norm);
    const double expected_frequency = time_component * (kEnergy - omega * kAngularMomentum);

    const auto transfer = sirius::core::relativity::KerrZamoFrequencyTransfer(
        kObserverFrequency, kEnergy, kAngularMomentum, kMass, kSpin, kRadius, kCosTheta);
    ASSERT_TRUE(transfer.has_value());
    EXPECT_NEAR(transfer->frame.angular_velocity, omega, 2.0e-15);
    EXPECT_NEAR(transfer->frame.time_component, time_component, 2.0e-15);
    EXPECT_NEAR(transfer->frame_frequency, expected_frequency, 2.0e-15);
    EXPECT_NEAR(transfer->g, kObserverFrequency / expected_frequency, 2.0e-15);

    const auto equatorial = sirius::core::relativity::KerrZamoFrequencyTransfer(
        kObserverFrequency, kEnergy, kAngularMomentum, kMass, kSpin, kRadius, 0.0);
    const auto disk = sirius::core::relativity::KerrDiskTransfer(
        kObserverFrequency, kEnergy, kAngularMomentum, kMass, kSpin, kRadius);
    ASSERT_TRUE(equatorial.has_value());
    ASSERT_TRUE(disk.has_value());
    EXPECT_DOUBLE_EQ(equatorial->frame_frequency, disk->zamo_frequency);
    EXPECT_DOUBLE_EQ(equatorial->g, disk->zamo_g);
}

TEST(KerrZamoTransfer, KerrSchildSlicingNormalIsNotTheOffEquatorialZamo) {
    constexpr double kMass = 1.0;
    constexpr double kSpin = 0.7;
    constexpr double kRadius = 6.0;
    constexpr double kCosTheta = 0.35;
    const sirius::core::coordinates::Vec4Bl event{0.0, kRadius, std::acos(kCosTheta), 0.4};
    const auto cart = sirius::core::coordinates::BlToKerrSchildCart(event, kSpin);

    sirius::core::KerrSchildFamily live_metric(sirius::core::KerrSchildParams::Kerr(kMass, kSpin));
    sirius::core::Vec4 position;
    position(0) = cart.t;
    position(1) = cart.x;
    position(2) = cart.y;
    position(3) = cart.z;
    sirius::core::Metric4d inverse;
    ASSERT_TRUE(live_metric.InverseMetric(position, inverse));
    const double lapse = 1.0 / std::sqrt(-inverse(0, 0).real);
    sirius::core::coordinates::Vec4Cart slicing_normal;
    slicing_normal.t = -lapse * inverse(0, 0).real;
    slicing_normal.x = -lapse * inverse(1, 0).real;
    slicing_normal.y = -lapse * inverse(2, 0).real;
    slicing_normal.z = -lapse * inverse(3, 0).real;
    const auto slicing_normal_bl = sirius::core::coordinates::TransformVectorKerrSchildCartToBl(
        slicing_normal, cart, kMass, kSpin);

    const auto zamo = sirius::core::relativity::KerrZamoFrequencyTransfer(
        1.0, 1.0, 0.0, kMass, kSpin, kRadius, kCosTheta);
    ASSERT_TRUE(zamo.has_value());
    EXPECT_GT(std::abs(slicing_normal_bl.r), 1.0e-3)
        << "the former volume diagnostic imported radial slicing motion";
    EXPECT_NE(slicing_normal_bl.phi / slicing_normal_bl.t, zamo->frame.angular_velocity);
}

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
        Vec4d(0, 10.0, std::numbers::pi / 2, 0),                    // Equatorial plane
        Vec4d(0, 5.0, std::numbers::pi / 4, 0),                     // Off-equator
        Vec4d(0, 3.0, std::numbers::pi / 3, std::numbers::pi / 2),  // Near photon sphere
        Vec4d(0, 100.0, std::numbers::pi / 2, 0),                   // Far field
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
        Vec4d(0, 10.0, std::numbers::pi / 2, 0),
        Vec4d(0, 5.0, std::numbers::pi / 4, 0),
        Vec4d(0, 3.0, std::numbers::pi / 3, std::numbers::pi / 2),
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

    Vec4d x(0, 6.0, std::numbers::pi / 3, 0);  // Test position
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
    Vec4d x(0, 10.0, std::numbers::pi / 2, 0);

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

TEST_F(KerrMetricDTest, ExactExtremalParametersAreNotSilentlyRewritten) {
    KerrMetricD extremal(1.0, 1.0);
    EXPECT_DOUBLE_EQ(extremal.mass(), 1.0);
    EXPECT_DOUBLE_EQ(extremal.spin(), 1.0);
    EXPECT_DOUBLE_EQ(extremal.charge(), 0.0);
    EXPECT_DOUBLE_EQ(extremal.HorizonRadius(), 1.0);
    EXPECT_DOUBLE_EQ(extremal.InnerHorizonRadius(), 1.0);
    EXPECT_DOUBLE_EQ(extremal.IscoRadius(), 1.0);
    EXPECT_DOUBLE_EQ(extremal.PhotonSphereRadius(), 1.0);

    KerrMetricD retrograde_extremal(1.0, -1.0);
    EXPECT_DOUBLE_EQ(retrograde_extremal.HorizonRadius(), 1.0);
    EXPECT_DOUBLE_EQ(retrograde_extremal.IscoRadius(), 9.0);
    EXPECT_DOUBLE_EQ(retrograde_extremal.PhotonSphereRadius(), 4.0);

    KerrMetricD flat_limit(0.0, 0.0);
    EXPECT_DOUBLE_EQ(flat_limit.mass(), 0.0);
    EXPECT_DOUBLE_EQ(flat_limit.spin(), 0.0);
    EXPECT_DOUBLE_EQ(flat_limit.HorizonRadius(), 0.0);
    EXPECT_TRUE(flat_limit.IsValid(Vec4d(0.0, 10.0, std::numbers::pi / 2.0, 0.0)));
    EXPECT_FALSE(flat_limit.IsValid(
        Vec4d(std::numeric_limits<double>::quiet_NaN(), 10.0, std::numbers::pi / 2.0, 0.0)));
    EXPECT_DEATH((void)flat_limit.ErgosphereRadius(std::numbers::pi / 2),
                 "precondition.*enforced, terminating");
    EXPECT_DEATH((void)flat_limit.IscoRadius(), "precondition.*enforced, terminating");
    EXPECT_DEATH((void)flat_limit.PhotonSphereRadius(), "precondition.*enforced, terminating");

    EXPECT_DEATH((void)KerrMetricD(1.0, 1.0001), "precondition.*enforced, terminating");
    EXPECT_DEATH((void)KerrMetricD(1.0, 0.5, 0.1), "precondition.*enforced, terminating");
    EXPECT_DEATH((void)KerrMetricD(-1.0, 0.0), "precondition.*enforced, terminating");
    EXPECT_DEATH((void)KerrMetricD(1.0, std::numeric_limits<double>::quiet_NaN()),
                 "precondition.*enforced, terminating");
}

TEST_F(KerrMetricDTest, BoyerLindquistAxisAndNonFinitePhaseSpaceDeclineWithoutClamping) {
    const Vec4d valid(0.0, 10.0, 2.0 * kBoyerLindquistPoleMargin, 0.0);
    const Vec4d north_axis(0.0, 10.0, 0.0, 0.0);
    const Vec4d inside_pole_margin(0.0, 10.0, 0.5 * kBoyerLindquistPoleMargin, 0.0);
    const Vec4d non_finite_event(std::numeric_limits<double>::quiet_NaN(), 10.0,
                                 std::numbers::pi / 2.0, 0.0);

    EXPECT_TRUE(IsRepresentedKerrBoyerLindquistEvent(1.0, 0.0, valid));
    EXPECT_TRUE(schwarzschild->IsValid(valid));
    EXPECT_FALSE(IsRepresentedKerrBoyerLindquistEvent(1.0, 0.0, north_axis));
    EXPECT_FALSE(IsRepresentedKerrBoyerLindquistEvent(1.0, 0.0, inside_pole_margin));
    EXPECT_FALSE(IsRepresentedKerrBoyerLindquistEvent(1.0, 0.0, non_finite_event));
    EXPECT_FALSE(schwarzschild->IsValid(north_axis));
    EXPECT_FALSE(schwarzschild->IsValid(non_finite_event));

    EXPECT_DEATH(schwarzschild->Evaluate(north_axis, nullptr, nullptr),
                 "precondition.*enforced, terminating");
    EXPECT_DEATH(KerrMetricDerivatives(1.0, 0.0, north_axis, nullptr),
                 "precondition.*enforced, terminating");
    EXPECT_DEATH((void)schwarzschild->Kretschmann(inside_pole_margin),
                 "precondition.*enforced, terminating");

    Vec4d non_finite_momentum(-1.0, 0.0, 0.0, 0.0);
    non_finite_momentum.r = std::numeric_limits<double>::infinity();
    EXPECT_DEATH((void)schwarzschild->Hamiltonian(valid, non_finite_momentum),
                 "precondition.*enforced, terminating");
    EXPECT_DEATH((void)schwarzschild->dHdq(valid, non_finite_momentum),
                 "precondition.*enforced, terminating");
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
                EXPECT_NEAR(Gamma[mu][nu][rho], fd, 1e-6 * std::max(1.0, std::abs(fd)))
                    << "Gamma^" << mu << "_" << nu << rho
                    << " disagrees with finite differences of the metric";
            }
        }
    }
}

// The hand-derived second metric derivatives must agree with central
// differences of the analytic first derivatives; they feed the Riemann
// assembly, so a slipped term here becomes a curvature error downstream.
// Both derivative orders are checked so the (s1, s2) symmetry of the storage
// is exercised, not assumed.
TEST_F(KerrMetricDTest, SecondDerivativesMatchFiniteDifference) {
    struct Case {
        double a, r, theta;
    };
    const Case cases[] = {{0.0, 5.0, 1.1}, {0.9, 5.0, 1.1}, {0.5, 3.2, 0.6}, {0.9, 12.0, 2.4}};

    const double h = 1e-6;
    for (const auto& c : cases) {
        double ddg[4][4][4][4];
        Vec4d x(0.0, c.r, c.theta, 0.3);
        sirius::oracle::KerrMetricSecondDerivatives(1.0, c.a, x, ddg);

        for (int s = 1; s <= 2; ++s) {
            Vec4d xp = x, xm = x;
            xp[s] += h;
            xm[s] -= h;
            double dgp[4][4][4], dgm[4][4][4];
            sirius::oracle::KerrMetricDerivatives(1.0, c.a, xp, dgp);
            sirius::oracle::KerrMetricDerivatives(1.0, c.a, xm, dgm);

            for (int s2 = 0; s2 < 4; ++s2) {
                for (int mu = 0; mu < 4; ++mu) {
                    for (int nu = 0; nu < 4; ++nu) {
                        double fd = (dgp[s2][mu][nu] - dgm[s2][mu][nu]) / (2.0 * h);
                        EXPECT_NEAR(ddg[s][s2][mu][nu], fd, 1e-6 * std::max(1.0, std::abs(fd)))
                            << "dd g[" << s << "][" << s2 << "][" << mu << "][" << nu
                            << "] disagrees with finite differences (a=" << c.a << ", r=" << c.r
                            << ", theta=" << c.theta << ")";
                    }
                }
            }
        }
    }
}
