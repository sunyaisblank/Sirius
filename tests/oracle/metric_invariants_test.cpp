// =============================================================================
// TSDG007A.cpp - Comprehensive Metric Tensor Validation
// Component ID: TSDG007A (Test/Diagnostic/Metric Validation)
// =============================================================================
//
// PURPOSE:
// Validates all metric tensor implementations against mathematical invariants
// defined in PHCN001A.h and tolerance constants in PHCN001A.h.
//
// TESTS:
// 1. Metric symmetry: max|g_μν - g_νμ| < SYMMETRY_TOL (1e-15)
// 2. Inverse accuracy: max|g^μα g_αν - δ^μ_ν| < INVERSE_TOL (1e-14)
// 3. Christoffel symmetry: max|Γ^λ_μν - Γ^λ_νμ| < CHRISTOFFEL_SYMMETRY_TOL (1e-15)
// 4. Lorentzian signature: exactly one negative eigenvalue
// 5. NaN/Inf detection at boundary conditions
// 6. Determinant non-degeneracy: |det(g)| > DETERMINANT_TOL (1e-30)
//
// METRICS TESTED:
// - KerrSchildFamily (PHMT100A): Minkowski, Schwarzschild, Kerr, Kerr-Newman
// - KerrMetricD (PHMT100B): Boyer-Lindquist Kerr
// - MorrisThorneFamily (PHMT101A): Ellis, ZeroTidal wormholes
// =============================================================================

#include "sirius/core/constants.h"
#include "sirius/core/dual_number.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/metrics/morris_thorne_family.h"
#include "sirius/core/tensor.h"
#include "sirius/oracle/kerr_boyer_lindquist.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace sirius::test {
using namespace sirius::core;
using sirius::oracle::KerrMetricD;

// =============================================================================
// Tolerance Constants (from PHCN001A.h)
// =============================================================================

using namespace sirius::core::constants;

constexpr double SYMMETRY_TOL = metric::kSymmetryTol;                         // 1e-15
constexpr double INVERSE_TOL = metric::kInverseTol;                           // 1e-14
constexpr double CHRISTOFFEL_SYMMETRY_TOL = metric::kChristoffelSymmetryTol;  // 1e-15
constexpr double DETERMINANT_TOL = metric::kDeterminantTol;                   // 1e-30
constexpr double SIGNATURE_TOL = metric::kSignatureTol;                       // 1e-10

// =============================================================================
// Test Fixture
// =============================================================================

class MetricValidationTests : public ::testing::Test {
  protected:
    void SetUp() override {}
    void TearDown() override {}

    // Sample points for testing (avoiding coordinate singularities)
    struct TestPoint {
        double t, x, y, z;     // Cartesian
        double r, theta, phi;  // Spherical (Boyer-Lindquist)
    };

    std::vector<TestPoint> getSamplePoints(double M = 1.0, double a = 0.0) {
        // Outer horizon for Kerr: r+ = M + sqrt(M² - a²)
        double r_plus = M + std::sqrt(std::max(M * M - a * a, 0.0));
        double r_safe = std::max(r_plus * 1.5, 3.0 * M);  // Safe distance from horizon

        return {
            // General positions (Cartesian: t, x, y, z | Spherical: r, θ, φ)
            {0, r_safe, 0, 0, r_safe, math::kHalfPi, 0},
            {0, 0, r_safe, 0, r_safe, math::kHalfPi, math::kHalfPi},
            {0, 0, 0, r_safe, r_safe, 0.1, 0},  // Near pole but not at pole
            {0, r_safe * 0.7, r_safe * 0.7, 0, r_safe, math::kHalfPi, math::kPi / 4},
            // Far field
            {0, 100 * M, 0, 0, 100 * M, math::kHalfPi, 0},
            {0, 50 * M, 50 * M, 0, std::sqrt(2) * 50 * M, math::kHalfPi, math::kPi / 4},
            // Intermediate distances
            {0, 10 * M, 0, 0, 10 * M, math::kHalfPi, 0},
            {0, 5 * M, 5 * M, 5 * M, std::sqrt(75) * M, std::acos(5 * M / (std::sqrt(75) * M)),
             math::kPi / 4},
        };
    }

    // Boundary test points (stress testing)
    std::vector<TestPoint> getBoundaryPoints(double M = 1.0, double a = 0.0) {
        double r_plus = M + std::sqrt(std::max(M * M - a * a, 0.0));

        return {
            // Very close to horizon (1.001 buffer per PHCN001A)
            {0, r_plus * 1.002, 0, 0, r_plus * 1.002, math::kHalfPi, 0},
            // Near pole (using POLE_EPSILON from PHCN001A)
            {0, 10 * M, 0, 0.01, 10 * M, 0.001, 0},
            // Large radius
            {0, 1e5 * M, 0, 0, 1e5 * M, math::kHalfPi, 0},
        };
    }
};

// =============================================================================
// Kerr-Schild Family Tests (PHMT100A)
// =============================================================================

// Test: Metric symmetry g_μν = g_νμ
TEST_F(MetricValidationTests, KerrSchild_MetricSymmetry) {
    std::vector<sirius::core::KerrSchildParams> configs = {
        sirius::core::KerrSchildParams::Minkowski(),
        sirius::core::KerrSchildParams::Schwarzschild(1.0),
        sirius::core::KerrSchildParams::Kerr(1.0, 0.5),
        sirius::core::KerrSchildParams::Kerr(1.0, 0.9),
        sirius::core::KerrSchildParams::ReissnerNordstrom(1.0, 0.5),
        sirius::core::KerrSchildParams::KerrNewman(1.0, 0.5, 0.3),
    };

    for (const auto& params : configs) {
        sirius::core::KerrSchildFamily metric(params);
        auto points = getSamplePoints(params.M, params.a);

        for (const auto& pt : points) {
            Tensor<double, 4> pos;
            pos(0) = pt.t;
            pos(1) = pt.x;
            pos(2) = pt.y;
            pos(3) = pt.z;

            Metric4d g;
            Tensor<Dual<double>, 4, 4, 4> dg;
            metric.Evaluate(pos, g, dg);

            double max_asym = 0;
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = mu + 1; nu < 4; ++nu) {
                    double asym = std::abs(g(mu, nu).real - g(nu, mu).real);
                    max_asym = std::max(max_asym, asym);
                }
            }

            EXPECT_LT(max_asym, SYMMETRY_TOL) << "Metric asymmetry for " << metric.GetName()
                                              << " at r=" << pt.r << ": " << max_asym;
        }
    }
}

// Test: Lorentzian signature (one negative eigenvalue)
TEST_F(MetricValidationTests, KerrSchild_LorentzianSignature) {
    std::vector<sirius::core::KerrSchildParams> configs = {
        sirius::core::KerrSchildParams::Schwarzschild(1.0),
        sirius::core::KerrSchildParams::Kerr(1.0, 0.5),
        sirius::core::KerrSchildParams::Kerr(1.0, 0.9),
    };

    for (const auto& params : configs) {
        sirius::core::KerrSchildFamily metric(params);
        auto points = getSamplePoints(params.M, params.a);

        for (const auto& pt : points) {
            Tensor<double, 4> pos;
            pos(0) = pt.t;
            pos(1) = pt.x;
            pos(2) = pt.y;
            pos(3) = pt.z;

            Metric4d g;
            Tensor<Dual<double>, 4, 4, 4> dg;
            metric.Evaluate(pos, g, dg);

            // For Kerr-Schild form: g = η + H·l⊗l where η = diag(-1,1,1,1)
            // The signature should remain (-,+,+,+)
            //
            // Simple check: g_tt should be negative (for static observer outside ergoregion)
            // and diagonal spatial components should be positive

            EXPECT_LT(g(0, 0).real, SIGNATURE_TOL)
                << "g_tt should be negative for " << metric.GetName() << " at r=" << pt.r;

            EXPECT_GT(g(1, 1).real, -SIGNATURE_TOL)
                << "g_xx should be positive for " << metric.GetName();
            EXPECT_GT(g(2, 2).real, -SIGNATURE_TOL)
                << "g_yy should be positive for " << metric.GetName();
            EXPECT_GT(g(3, 3).real, -SIGNATURE_TOL)
                << "g_zz should be positive for " << metric.GetName();
        }
    }
}

// Test: NaN/Inf detection at safe positions
TEST_F(MetricValidationTests, KerrSchild_NoNaNInf) {
    std::vector<sirius::core::KerrSchildParams> configs = {
        sirius::core::KerrSchildParams::Schwarzschild(1.0),
        sirius::core::KerrSchildParams::Kerr(1.0, 0.9),
    };

    for (const auto& params : configs) {
        sirius::core::KerrSchildFamily metric(params);
        auto points = getSamplePoints(params.M, params.a);

        for (const auto& pt : points) {
            Tensor<double, 4> pos;
            pos(0) = pt.t;
            pos(1) = pt.x;
            pos(2) = pt.y;
            pos(3) = pt.z;

            Metric4d g;
            Tensor<Dual<double>, 4, 4, 4> dg;
            metric.Evaluate(pos, g, dg);

            // Check all metric components
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    EXPECT_FALSE(std::isnan(g(mu, nu).real))
                        << "NaN in g_" << mu << nu << " for " << metric.GetName() << " at (" << pt.x
                        << "," << pt.y << "," << pt.z << ")";
                    EXPECT_FALSE(std::isinf(g(mu, nu).real))
                        << "Inf in g_" << mu << nu << " for " << metric.GetName();
                }
            }

            // Check all derivative components
            for (int lam = 0; lam < 4; ++lam) {
                for (int mu = 0; mu < 4; ++mu) {
                    for (int nu = 0; nu < 4; ++nu) {
                        EXPECT_FALSE(std::isnan(dg(lam, mu, nu).real))
                            << "NaN in dg_" << lam << mu << nu << " for " << metric.GetName();
                        EXPECT_FALSE(std::isinf(dg(lam, mu, nu).real))
                            << "Inf in dg_" << lam << mu << nu << " for " << metric.GetName();
                    }
                }
            }
        }
    }
}

// Test: Minkowski limit (M → 0)
TEST_F(MetricValidationTests, KerrSchild_MinkowskiLimit) {
    sirius::core::KerrSchildFamily metric(sirius::core::KerrSchildParams::Minkowski());

    // At any point, Minkowski metric should be η = diag(-1,1,1,1)
    std::vector<std::array<double, 3>> positions = {{10, 0, 0}, {0, 10, 0}, {5, 5, 5}, {100, 0, 0}};

    for (const auto& xyz : positions) {
        Tensor<double, 4> pos;
        pos(0) = 0;
        pos(1) = xyz[0];
        pos(2) = xyz[1];
        pos(3) = xyz[2];

        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.Evaluate(pos, g, dg);

        // Check diagonal components
        EXPECT_NEAR(g(0, 0).real, -1.0, 1e-14) << "g_tt != -1 for Minkowski";
        EXPECT_NEAR(g(1, 1).real, 1.0, 1e-14) << "g_xx != 1 for Minkowski";
        EXPECT_NEAR(g(2, 2).real, 1.0, 1e-14) << "g_yy != 1 for Minkowski";
        EXPECT_NEAR(g(3, 3).real, 1.0, 1e-14) << "g_zz != 1 for Minkowski";

        // Check off-diagonal components are zero
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                if (mu != nu) {
                    EXPECT_NEAR(g(mu, nu).real, 0.0, 1e-14)
                        << "g_" << mu << nu << " != 0 for Minkowski";
                }
            }
        }

        // All derivatives should be zero
        for (int lam = 0; lam < 4; ++lam) {
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    EXPECT_NEAR(dg(lam, mu, nu).real, 0.0, 1e-14)
                        << "dg_" << lam << mu << nu << " != 0 for Minkowski";
                }
            }
        }
    }
}

// Test: Schwarzschild weak field limit (r >> 2M)
TEST_F(MetricValidationTests, KerrSchild_SchwarzschildWeakField) {
    double M = 1.0;
    sirius::core::KerrSchildFamily metric(sirius::core::KerrSchildParams::Schwarzschild(M));

    // At r >> 2M, g_tt ≈ -(1 - 2M/r), g_rr ≈ 1 + 2M/r
    std::vector<double> radii = {100, 500, 1000};

    for (double r : radii) {
        Tensor<double, 4> pos;
        pos(0) = 0;
        pos(1) = r;
        pos(2) = 0;
        pos(3) = 0;

        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.Evaluate(pos, g, dg);

        double expected_g_tt = -(1.0 - 2.0 * M / r);
        double tol = 4.0 * M * M / (r * r);  // O(1/r²) corrections

        EXPECT_NEAR(g(0, 0).real, expected_g_tt, tol) << "Weak field g_tt error at r=" << r;

        // In Kerr-Schild Cartesian, g_rr ≈ 1 + 2M/r on the x-axis
        // But full tensor form includes off-diagonal corrections
        // Check spatial components are approximately δ_ij
        EXPECT_NEAR(g(1, 1).real, 1.0 + 2.0 * M / r, tol) << "Weak field g_xx error at r=" << r;
    }
}

// =============================================================================
// Boyer-Lindquist Kerr Tests (PHMT100B)
// =============================================================================

// Test: Metric symmetry for KerrMetricD
TEST_F(MetricValidationTests, KerrMetricD_MetricSymmetry) {
    std::vector<double> spins = {0.0, 0.3, 0.6, 0.9};

    for (double a : spins) {
        sirius::oracle::KerrMetricD metric(1.0, a);
        auto points = getSamplePoints(1.0, a);

        for (const auto& pt : points) {
            sirius::oracle::Vec4d x;
            x.t = pt.t;
            x.r = pt.r;
            x.theta = pt.theta;
            x.phi = pt.phi;

            double g[4][4], g_inv[4][4];
            metric.Evaluate(x, g, g_inv);

            double max_asym = 0;
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = mu + 1; nu < 4; ++nu) {
                    double asym = std::abs(g[mu][nu] - g[nu][mu]);
                    max_asym = std::max(max_asym, asym);
                }
            }

            EXPECT_LT(max_asym, SYMMETRY_TOL)
                << "Metric asymmetry for Kerr (a=" << a << ") at r=" << pt.r;
        }
    }
}

// Test: Inverse metric accuracy g^μα g_αν = δ^μ_ν
TEST_F(MetricValidationTests, KerrMetricD_InverseAccuracy) {
    std::vector<double> spins = {0.0, 0.5, 0.9};

    for (double a : spins) {
        sirius::oracle::KerrMetricD metric(1.0, a);
        auto points = getSamplePoints(1.0, a);

        for (const auto& pt : points) {
            sirius::oracle::Vec4d x;
            x.t = pt.t;
            x.r = pt.r;
            x.theta = pt.theta;
            x.phi = pt.phi;

            double g[4][4], g_inv[4][4];
            metric.Evaluate(x, g, g_inv);

            // Compute g^μα g_αν
            double max_error = 0;
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    double product = 0;
                    for (int alpha = 0; alpha < 4; ++alpha) {
                        product += g_inv[mu][alpha] * g[alpha][nu];
                    }

                    double expected = (mu == nu) ? 1.0 : 0.0;
                    double error = std::abs(product - expected);
                    max_error = std::max(max_error, error);
                }
            }

            EXPECT_LT(max_error, INVERSE_TOL) << "Inverse metric error for Kerr (a=" << a
                                              << ") at r=" << pt.r << ": " << max_error;
        }
    }
}

// Test: Christoffel symbol symmetry Γ^λ_μν = Γ^λ_νμ
TEST_F(MetricValidationTests, KerrMetricD_ChristoffelSymmetry) {
    std::vector<double> spins = {0.0, 0.5, 0.9};

    for (double a : spins) {
        sirius::oracle::KerrMetricD metric(1.0, a);
        auto points = getSamplePoints(1.0, a);

        for (const auto& pt : points) {
            sirius::oracle::Vec4d x;
            x.t = pt.t;
            x.r = pt.r;
            x.theta = pt.theta;
            x.phi = pt.phi;

            double Gamma[4][4][4];
            metric.Christoffel(x, Gamma);

            double max_asym = 0;
            for (int lam = 0; lam < 4; ++lam) {
                for (int mu = 0; mu < 4; ++mu) {
                    for (int nu = mu + 1; nu < 4; ++nu) {
                        double asym = std::abs(Gamma[lam][mu][nu] - Gamma[lam][nu][mu]);
                        max_asym = std::max(max_asym, asym);
                    }
                }
            }

            EXPECT_LT(max_asym, CHRISTOFFEL_SYMMETRY_TOL)
                << "Christoffel asymmetry for Kerr (a=" << a << ") at r=" << pt.r << ": "
                << max_asym;
        }
    }
}

// Test: Riemann tensor antisymmetry R^α_βγδ = -R^α_βδγ
TEST_F(MetricValidationTests, KerrMetricD_RiemannAntisymmetry) {
    std::vector<double> spins = {0.0, 0.5, 0.9};

    for (double a : spins) {
        sirius::oracle::KerrMetricD metric(1.0, a);
        auto points = getSamplePoints(1.0, a);

        for (const auto& pt : points) {
            sirius::oracle::Vec4d x;
            x.t = pt.t;
            x.r = pt.r;
            x.theta = pt.theta;
            x.phi = pt.phi;

            double violation = metric.VerifyRiemannSymmetries(x);

            EXPECT_LT(violation, 1e-12) << "Riemann antisymmetry violation for Kerr (a=" << a
                                        << ") at r=" << pt.r << ": " << violation;
        }
    }
}

// Test: Kretschmann scalar (analytic validation)
// Schwarzschild: K = 48M²/r⁶
// Kerr: K = 48M²(r² - a²cos²θ)[(r² + a²cos²θ)² - 16a²r²cos²θ] / Σ⁶
// Reference: Henry, R.C. (2000), Astrophys. J. 535:350-353, eq. 18 at Q = 0
TEST_F(MetricValidationTests, KerrMetricD_KretschmannScalar) {
    // Test 1: Schwarzschild exact formula K = 48M²/r⁶
    sirius::oracle::KerrMetricD schw(1.0, 0.0);

    std::vector<double> radii = {3.0, 5.0, 10.0, 20.0};

    for (double r : radii) {
        sirius::oracle::Vec4d x;
        x.t = 0;
        x.r = r;
        x.theta = math::kHalfPi;
        x.phi = 0;

        double K = schw.Kretschmann(x);
        double expected = 48.0 / std::pow(r, 6);  // M=1

        double rel_error = std::abs(K - expected) / expected;

        EXPECT_LT(rel_error, 1e-10) << "Schwarzschild Kretschmann error at r=" << r
                                    << ": computed=" << K << ", expected=" << expected;
    }

    // Test 2: Kerr at equator (θ = π/2, so cos²θ = 0)
    // At equator: K = 48M² × r² × r⁴ / (r²)⁶ = 48M² × r⁶ / r¹² = 48M²/r⁶
    // Same as Schwarzschild!
    sirius::oracle::KerrMetricD kerr(1.0, 0.9);

    for (double r : radii) {
        sirius::oracle::Vec4d x;
        x.t = 0;
        x.r = r;
        x.theta = math::kHalfPi;
        x.phi = 0;

        double K = kerr.Kretschmann(x);
        double expected = 48.0 / std::pow(r, 6);  // Same as Schwarzschild at equator

        double rel_error = std::abs(K - expected) / expected;

        EXPECT_LT(rel_error, 1e-10) << "Kerr equatorial Kretschmann error at r=" << r
                                    << ": computed=" << K << ", expected=" << expected;
    }

    // Test 3: Kerr off-equator (θ = π/4)
    // Full formula applies with non-zero cos²θ
    double theta = math::kPi / 4;
    double cos2th = 0.5;  // cos²(π/4) = 0.5
    double a = 0.5;
    double a2 = a * a;

    sirius::oracle::KerrMetricD kerr05(1.0, a);

    for (double r : radii) {
        sirius::oracle::Vec4d x;
        x.t = 0;
        x.r = r;
        x.theta = theta;
        x.phi = 0;

        double r2 = r * r;
        double Sigma = r2 + a2 * cos2th;
        double Sigma6 = std::pow(Sigma, 6);
        double r_term = r2 - a2 * cos2th;
        // Henry (2000) eq. 18 at Q = 0: the bracket is Σ² - 16a²r²cos²θ.
        // (A prior version squared r_term here, matching a defect in the
        // implementation — a circular pin. The Riemann-contraction gate now
        // anchors this expression independently.)
        double bracket = Sigma * Sigma - 16.0 * a2 * r2 * cos2th;
        double expected = 48.0 * r_term * bracket / Sigma6;  // M=1

        double K = kerr05.Kretschmann(x);

        double rel_error = std::abs(K - expected) / std::max(std::abs(expected), 1e-30);

        EXPECT_LT(rel_error, 1e-10) << "Kerr off-axis Kretschmann error at r=" << r
                                    << ": computed=" << K << ", expected=" << expected;
    }

    // Test 4: Kretschmann decreases with radius (monotonicity)
    double K_prev = std::numeric_limits<double>::max();
    for (double r : radii) {
        sirius::oracle::Vec4d x;
        x.t = 0;
        x.r = r;
        x.theta = math::kHalfPi;
        x.phi = 0;

        double K = schw.Kretschmann(x);
        EXPECT_LT(K, K_prev) << "Kretschmann should decrease with r";
        K_prev = K;
    }

    // Test 5: At large radius, Kretschmann approaches zero (asymptotic flatness)
    sirius::oracle::Vec4d x_far;
    x_far.t = 0;
    x_far.r = 1000.0;
    x_far.theta = math::kHalfPi;
    x_far.phi = 0;
    double K_far = schw.Kretschmann(x_far);
    double expected_far = 48.0 / std::pow(1000.0, 6);  // ~4.8e-17
    EXPECT_NEAR(K_far, expected_far, 1e-20) << "Kretschmann should match 48/r⁶ at large radius";
}

namespace {

// Reference Riemann tensor R^rho_sigma_mu_nu assembled from central differences
// of the analytic Christoffel connection plus the quadratic Gamma*Gamma terms:
//   R^rho_smn = d_m Gamma^rho_ns - d_n Gamma^rho_ms
//             + Gamma^rho_ml Gamma^l_ns - Gamma^rho_nl Gamma^l_ms
// The connection is exact, so the only error is the O(h^2) truncation of the
// two coordinate derivatives (r and theta; the metric is stationary and
// axisymmetric). This is the completeness gate the component-wise
// implementation historically lacked: it exercises every one of the 256 slots.
void RiemannFromChristoffelFD(const sirius::oracle::KerrMetricD& metric,
                              const sirius::oracle::Vec4d& x, double h, double R[4][4][4][4]) {
    double Gamma[4][4][4];
    metric.Christoffel(x, Gamma);

    // dGamma[s][m][n][r] = d_s Gamma^m_nr; only s = r, theta are non-zero.
    double dGamma[4][4][4][4] = {};
    for (int s = 1; s <= 2; ++s) {
        sirius::oracle::Vec4d xp = x;
        sirius::oracle::Vec4d xm = x;
        if (s == 1) {
            xp.r += h;
            xm.r -= h;
        } else {
            xp.theta += h;
            xm.theta -= h;
        }
        double Gp[4][4][4], Gm[4][4][4];
        metric.Christoffel(xp, Gp);
        metric.Christoffel(xm, Gm);
        for (int m = 0; m < 4; ++m)
            for (int n = 0; n < 4; ++n)
                for (int r = 0; r < 4; ++r)
                    dGamma[s][m][n][r] = (Gp[m][n][r] - Gm[m][n][r]) / (2.0 * h);
    }

    for (int rho = 0; rho < 4; ++rho) {
        for (int sig = 0; sig < 4; ++sig) {
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    double val = dGamma[mu][rho][nu][sig] - dGamma[nu][rho][mu][sig];
                    for (int lam = 0; lam < 4; ++lam) {
                        val += Gamma[rho][mu][lam] * Gamma[lam][nu][sig] -
                               Gamma[rho][nu][lam] * Gamma[lam][mu][sig];
                    }
                    R[rho][sig][mu][nu] = val;
                }
            }
        }
    }
}

double MaxAbsComponent(const double R[4][4][4][4]) {
    double m = 0.0;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
                for (int d = 0; d < 4; ++d) m = std::max(m, std::abs(R[a][b][c][d]));
    return m;
}

struct RiemannTestCase {
    double a, r, theta;
};

const std::vector<RiemannTestCase>& RiemannTestCases() {
    static const std::vector<RiemannTestCase> cases = {
        {0.0, 6.0, math::kHalfPi},    // Schwarzschild, equator
        {0.0, 3.1, 1.0},              // Schwarzschild, near photon sphere
        {0.5, 6.0, math::kPi / 3.0},  // Moderate spin, off-equator
        {0.9, 4.0, 1.1},              // High spin, strong field
        {0.9, 10.0, 0.4},             // High spin, toward the axis
    };
    return cases;
}

}  // namespace

// Test: every Riemann component matches the finite-difference assembly from
// the analytic connection. This is the gate that fails on any missing or
// mis-derived component, on or off the diagonal blocks.
TEST_F(MetricValidationTests, KerrMetricD_RiemannMatchesFiniteDifferenceChristoffel) {
    const double h = 1e-5;

    for (const auto& c : RiemannTestCases()) {
        sirius::oracle::KerrMetricD metric(1.0, c.a);
        sirius::oracle::Vec4d x;
        x.t = 0;
        x.r = c.r;
        x.theta = c.theta;
        x.phi = 0;

        double Rfd[4][4][4][4], Ran[4][4][4][4];
        RiemannFromChristoffelFD(metric, x, h, Rfd);
        metric.Riemann(x, Ran);

        // FD truncation is O(h^2) on connection third derivatives; at these
        // radii the curvature scale is O(1e-2), so 1e-7 relative headroom
        // dominates the truncation and rounding floor comfortably.
        double scale = MaxAbsComponent(Rfd);
        double tol = std::max(1e-10, 1e-7 * scale);

        double worst = 0.0;
        int wi[4] = {0, 0, 0, 0};
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int p = 0; p < 4; ++p)
                    for (int q = 0; q < 4; ++q) {
                        double err = std::abs(Ran[a][b][p][q] - Rfd[a][b][p][q]);
                        if (err > worst) {
                            worst = err;
                            wi[0] = a;
                            wi[1] = b;
                            wi[2] = p;
                            wi[3] = q;
                        }
                    }

        EXPECT_LT(worst, tol) << "Riemann mismatch vs finite differences for a=" << c.a
                              << " r=" << c.r << " theta=" << c.theta << " at R^" << wi[0] << "_"
                              << wi[1] << wi[2] << wi[3]
                              << ": analytic=" << Ran[wi[0]][wi[1]][wi[2]][wi[3]]
                              << ", fd=" << Rfd[wi[0]][wi[1]][wi[2]][wi[3]];
    }
}

// Test: Kerr is a vacuum solution, so the Ricci contraction of the Riemann
// tensor must vanish. A partial component table generically fails this.
TEST_F(MetricValidationTests, KerrMetricD_RiemannVacuumRicci) {
    for (const auto& c : RiemannTestCases()) {
        sirius::oracle::KerrMetricD metric(1.0, c.a);
        sirius::oracle::Vec4d x;
        x.t = 0;
        x.r = c.r;
        x.theta = c.theta;
        x.phi = 0;

        double R[4][4][4][4];
        metric.Riemann(x, R);
        double scale = MaxAbsComponent(R);

        double worst = 0.0;
        for (int sig = 0; sig < 4; ++sig) {
            for (int nu = 0; nu < 4; ++nu) {
                double ricci = 0.0;
                for (int rho = 0; rho < 4; ++rho) ricci += R[rho][sig][rho][nu];
                worst = std::max(worst, std::abs(ricci));
            }
        }

        EXPECT_LT(worst, std::max(1e-12, 1e-9 * scale))
            << "Non-zero Ricci for vacuum Kerr a=" << c.a << " r=" << c.r << " theta=" << c.theta
            << ": " << worst;
    }
}

// Test: with the first index lowered, the full symmetry group holds:
// antisymmetry in each index pair, pair-interchange symmetry, and the first
// Bianchi identity.
TEST_F(MetricValidationTests, KerrMetricD_RiemannLoweredSymmetries) {
    for (const auto& c : RiemannTestCases()) {
        sirius::oracle::KerrMetricD metric(1.0, c.a);
        sirius::oracle::Vec4d x;
        x.t = 0;
        x.r = c.r;
        x.theta = c.theta;
        x.phi = 0;

        double R[4][4][4][4], g[4][4], g_inv[4][4];
        metric.Riemann(x, R);
        metric.Evaluate(x, g, g_inv);

        double Rl[4][4][4][4];
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int p = 0; p < 4; ++p)
                    for (int q = 0; q < 4; ++q) {
                        double v = 0.0;
                        for (int m = 0; m < 4; ++m) v += g[a][m] * R[m][b][p][q];
                        Rl[a][b][p][q] = v;
                    }

        double scale = MaxAbsComponent(Rl);
        double tol = std::max(1e-12, 1e-9 * scale);

        double worst_first = 0.0, worst_pair = 0.0, worst_bianchi = 0.0;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int p = 0; p < 4; ++p)
                    for (int q = 0; q < 4; ++q) {
                        worst_first =
                            std::max(worst_first, std::abs(Rl[a][b][p][q] + Rl[b][a][p][q]));
                        worst_pair =
                            std::max(worst_pair, std::abs(Rl[a][b][p][q] - Rl[p][q][a][b]));
                        worst_bianchi =
                            std::max(worst_bianchi,
                                     std::abs(Rl[a][b][p][q] + Rl[a][p][q][b] + Rl[a][q][b][p]));
                    }

        EXPECT_LT(worst_first, tol) << "R_abcd != -R_bacd for a=" << c.a << " r=" << c.r;
        EXPECT_LT(worst_pair, tol) << "R_abcd != R_cdab for a=" << c.a << " r=" << c.r;
        EXPECT_LT(worst_bianchi, tol) << "First Bianchi violated for a=" << c.a << " r=" << c.r;
    }
}

// Test: contracting the tensor against itself reproduces the closed-form
// Kretschmann scalar. This ties the component table to the validated
// invariant: K = R_abcd R^abcd.
TEST_F(MetricValidationTests, KerrMetricD_KretschmannMatchesRiemannContraction) {
    for (const auto& c : RiemannTestCases()) {
        sirius::oracle::KerrMetricD metric(1.0, c.a);
        sirius::oracle::Vec4d x;
        x.t = 0;
        x.r = c.r;
        x.theta = c.theta;
        x.phi = 0;

        double R[4][4][4][4], g[4][4], g_inv[4][4];
        metric.Riemann(x, R);
        metric.Evaluate(x, g, g_inv);

        // Lower the first index, raise the last three, contract.
        double Rl[4][4][4][4];
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int p = 0; p < 4; ++p)
                    for (int q = 0; q < 4; ++q) {
                        double v = 0.0;
                        for (int m = 0; m < 4; ++m) v += g[a][m] * R[m][b][p][q];
                        Rl[a][b][p][q] = v;
                    }

        double K = 0.0;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b)
                for (int p = 0; p < 4; ++p)
                    for (int q = 0; q < 4; ++q) {
                        // R^abpq = g^aA g^bB g^pP g^qQ R_ABPQ, assembled per slot.
                        double up = 0.0;
                        for (int A = 0; A < 4; ++A)
                            for (int B = 0; B < 4; ++B)
                                for (int P = 0; P < 4; ++P)
                                    for (int Q = 0; Q < 4; ++Q)
                                        up += g_inv[a][A] * g_inv[b][B] * g_inv[p][P] *
                                              g_inv[q][Q] * Rl[A][B][P][Q];
                        K += Rl[a][b][p][q] * up;
                    }

        double expected = metric.Kretschmann(x);
        double rel = std::abs(K - expected) / std::max(std::abs(expected), 1e-30);
        EXPECT_LT(rel, 1e-9) << "Kretschmann contraction mismatch for a=" << c.a << " r=" << c.r
                             << " theta=" << c.theta << ": contracted=" << K
                             << ", analytic=" << expected;
    }
}

// Test: Horizon radius calculation
TEST_F(MetricValidationTests, KerrMetricD_HorizonRadius) {
    // Schwarzschild: r+ = 2M
    sirius::oracle::KerrMetricD schw(1.0, 0.0);
    EXPECT_NEAR(schw.HorizonRadius(), 2.0, 1e-14);

    // Kerr with a = 0.5: r+ = 1 + sqrt(1 - 0.25) = 1 + sqrt(0.75)
    sirius::oracle::KerrMetricD kerr05(1.0, 0.5);
    double expected_05 = 1.0 + std::sqrt(0.75);
    EXPECT_NEAR(kerr05.HorizonRadius(), expected_05, 1e-14);

    // Near-extremal Kerr a ≈ M: r+ ≈ M + sqrt(1 - a²)
    // Note: KerrMetricD constructor clamps spin to 0.9999*M for numerical stability
    sirius::oracle::KerrMetricD kerr_ext(1.0, 0.9999);
    double expected_ext = 1.0 + std::sqrt(1.0 - 0.9999 * 0.9999);  // ~1.0141
    EXPECT_NEAR(kerr_ext.HorizonRadius(), expected_ext, 1e-10);
}

// Test: ISCO radius calculation
TEST_F(MetricValidationTests, KerrMetricD_ISCORadius) {
    // Schwarzschild: r_ISCO = 6M
    sirius::oracle::KerrMetricD schw(1.0, 0.0);
    EXPECT_NEAR(schw.IscoRadius(), 6.0, 1e-10);

    // Prograde extremal Kerr: r_ISCO = M
    sirius::oracle::KerrMetricD kerr099(1.0, 0.999);
    EXPECT_LT(kerr099.IscoRadius(), 2.0);  // Should be close to M
}

// Test: No NaN/Inf in valid domain
TEST_F(MetricValidationTests, KerrMetricD_NoNaNInf) {
    std::vector<double> spins = {0.0, 0.5, 0.9};

    for (double a : spins) {
        sirius::oracle::KerrMetricD metric(1.0, a);
        auto points = getSamplePoints(1.0, a);

        for (const auto& pt : points) {
            sirius::oracle::Vec4d x;
            x.t = pt.t;
            x.r = pt.r;
            x.theta = pt.theta;
            x.phi = pt.phi;

            if (!metric.IsValid(x)) continue;  // Skip invalid points

            double g[4][4], g_inv[4][4];
            metric.Evaluate(x, g, g_inv);

            double Gamma[4][4][4];
            metric.Christoffel(x, Gamma);

            // Check metric
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    EXPECT_FALSE(std::isnan(g[mu][nu]))
                        << "NaN in g[" << mu << "][" << nu << "] for a=" << a;
                    EXPECT_FALSE(std::isinf(g[mu][nu]))
                        << "Inf in g[" << mu << "][" << nu << "] for a=" << a;
                }
            }

            // Check Christoffel
            for (int lam = 0; lam < 4; ++lam) {
                for (int mu = 0; mu < 4; ++mu) {
                    for (int nu = 0; nu < 4; ++nu) {
                        EXPECT_FALSE(std::isnan(Gamma[lam][mu][nu]))
                            << "NaN in Gamma[" << lam << "][" << mu << "][" << nu << "]";
                        EXPECT_FALSE(std::isinf(Gamma[lam][mu][nu]))
                            << "Inf in Gamma[" << lam << "][" << mu << "][" << nu << "]";
                    }
                }
            }
        }
    }
}

// =============================================================================
// Morris-Thorne Wormhole Tests (PHMT101A)
// =============================================================================

// Test: Metric symmetry for wormholes
TEST_F(MetricValidationTests, MorrisThorne_MetricSymmetry) {
    std::vector<sirius::core::MorrisThorneParams> configs = {
        sirius::core::MorrisThorneParams::Ellis(1.0),
        sirius::core::MorrisThorneParams::ZeroTidal(1.0),
    };

    for (const auto& params : configs) {
        sirius::core::MorrisThorneFamily metric(params);

        // Test points outside throat (r > b0)
        std::vector<double> radii = {1.5, 2.0, 5.0, 10.0, 100.0};

        for (double r : radii) {
            Tensor<double, 4> pos;
            pos(0) = 0;              // t
            pos(1) = r;              // r
            pos(2) = math::kHalfPi;  // θ = π/2
            pos(3) = 0;              // φ

            Metric4d g;
            Tensor<Dual<double>, 4, 4, 4> dg;
            metric.Evaluate(pos, g, dg);

            // Wormhole metric is diagonal, so symmetry is trivial
            // Just verify it's diagonal
            double max_offdiag = 0;
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    if (mu != nu) {
                        max_offdiag = std::max(max_offdiag, std::abs(g(mu, nu).real));
                    }
                }
            }

            EXPECT_LT(max_offdiag, SYMMETRY_TOL)
                << "Non-diagonal component for " << metric.GetName() << " at r=" << r;
        }
    }
}

// Test: Flare-out condition (traversability)
TEST_F(MetricValidationTests, MorrisThorne_FlareOutCondition) {
    // Ellis wormhole should satisfy flare-out
    sirius::core::MorrisThorneFamily ellis(sirius::core::MorrisThorneParams::Ellis(1.0));
    EXPECT_TRUE(ellis.ValidateFlareOutCondition())
        << "Ellis wormhole should satisfy flare-out condition";

    // Zero-tidal wormhole should satisfy flare-out
    sirius::core::MorrisThorneFamily zero_tidal(sirius::core::MorrisThorneParams::ZeroTidal(1.0));
    EXPECT_TRUE(zero_tidal.ValidateFlareOutCondition())
        << "Zero-tidal wormhole should satisfy flare-out condition";
}

// Test: Shape function at throat b(b0) = b0
TEST_F(MetricValidationTests, MorrisThorne_ThroatCondition) {
    double b0 = 2.0;

    // Ellis: b(r) = b0²/r → b(b0) = b0²/b0 = b0
    sirius::core::MorrisThorneFamily ellis(sirius::core::MorrisThorneParams::Ellis(b0));
    EXPECT_NEAR(ellis.ShapeFunction(b0), b0, 1e-14) << "Ellis throat condition violated";

    // Zero-tidal: b(r) = b0 → b(b0) = b0
    sirius::core::MorrisThorneFamily zero_tidal(sirius::core::MorrisThorneParams::ZeroTidal(b0));
    EXPECT_NEAR(zero_tidal.ShapeFunction(b0), b0, 1e-14) << "Zero-tidal throat condition violated";
}

// Test: No NaN/Inf outside throat
TEST_F(MetricValidationTests, MorrisThorne_NoNaNInf) {
    sirius::core::MorrisThorneFamily metric(sirius::core::MorrisThorneParams::Ellis(1.0));

    std::vector<double> radii = {1.01, 1.1, 2.0, 10.0, 100.0};
    std::vector<double> thetas = {0.1, math::kHalfPi, math::kPi - 0.1};

    for (double r : radii) {
        for (double theta : thetas) {
            Tensor<double, 4> pos;
            pos(0) = 0;
            pos(1) = r;
            pos(2) = theta;
            pos(3) = 0;

            Metric4d g;
            Tensor<Dual<double>, 4, 4, 4> dg;
            metric.Evaluate(pos, g, dg);

            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    EXPECT_FALSE(std::isnan(g(mu, nu).real))
                        << "NaN in wormhole metric at r=" << r << ", θ=" << theta;
                    EXPECT_FALSE(std::isinf(g(mu, nu).real))
                        << "Inf in wormhole metric at r=" << r << ", θ=" << theta;
                }
            }
        }
    }
}

// Test: Lorentzian signature for wormhole
TEST_F(MetricValidationTests, MorrisThorne_LorentzianSignature) {
    sirius::core::MorrisThorneFamily metric(sirius::core::MorrisThorneParams::Ellis(1.0));

    std::vector<double> radii = {1.1, 2.0, 10.0};

    for (double r : radii) {
        Tensor<double, 4> pos;
        pos(0) = 0;
        pos(1) = r;
        pos(2) = math::kHalfPi;
        pos(3) = 0;

        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.Evaluate(pos, g, dg);

        // Morris-Thorne: g_tt = -exp(2Φ) < 0
        EXPECT_LT(g(0, 0).real, 0) << "g_tt should be negative at r=" << r;

        // g_rr = 1/(1 - b/r) > 0 (for r > b)
        EXPECT_GT(g(1, 1).real, 0) << "g_rr should be positive at r=" << r;

        // g_θθ = r² > 0
        EXPECT_GT(g(2, 2).real, 0) << "g_θθ should be positive at r=" << r;

        // g_φφ = r² sin²θ ≥ 0
        EXPECT_GE(g(3, 3).real, 0) << "g_φφ should be non-negative at r=" << r;
    }
}

// =============================================================================
// Cross-Implementation Consistency Tests
// =============================================================================

// Test: Schwarzschild metric consistency between implementations
// KerrSchildFamily (Cartesian) vs KerrMetricD (Boyer-Lindquist)
TEST_F(MetricValidationTests, SchwarzschildConsistency) {
    double M = 1.0;
    sirius::core::KerrSchildFamily ks_metric(sirius::core::KerrSchildParams::Schwarzschild(M));
    sirius::oracle::KerrMetricD bl_metric(M, 0.0);

    // Test at equatorial positions
    std::vector<double> radii = {3.0, 6.0, 10.0, 20.0};

    for (double r : radii) {
        // Kerr-Schild evaluation (Cartesian: x = r, y = z = 0)
        Tensor<double, 4> pos_ks;
        pos_ks(0) = 0;
        pos_ks(1) = r;
        pos_ks(2) = 0;
        pos_ks(3) = 0;

        Metric4d g_ks;
        Tensor<Dual<double>, 4, 4, 4> dg_ks;
        ks_metric.Evaluate(pos_ks, g_ks, dg_ks);

        // Boyer-Lindquist evaluation
        sirius::oracle::Vec4d pos_bl;
        pos_bl.t = 0;
        pos_bl.r = r;
        pos_bl.theta = math::kHalfPi;
        pos_bl.phi = 0;

        double g_bl[4][4], g_inv_bl[4][4];
        bl_metric.Evaluate(pos_bl, g_bl, g_inv_bl);

        // For Schwarzschild on equator with φ=0:
        // g_tt should match
        double g_tt_ks = g_ks(0, 0).real;
        double g_tt_bl = g_bl[0][0];

        // Both should equal -(1 - 2M/r)
        double expected_g_tt = -(1.0 - 2.0 * M / r);

        EXPECT_NEAR(g_tt_ks, expected_g_tt, 1e-10) << "KS g_tt error at r=" << r;
        EXPECT_NEAR(g_tt_bl, expected_g_tt, 1e-10) << "BL g_tt error at r=" << r;
    }
}

// =============================================================================
// Determinant Non-Degeneracy Tests
// =============================================================================

// Test: Metric determinant is non-zero in valid domain
TEST_F(MetricValidationTests, KerrMetricD_DeterminantNonZero) {
    sirius::oracle::KerrMetricD metric(1.0, 0.5);
    auto points = getSamplePoints(1.0, 0.5);

    for (const auto& pt : points) {
        sirius::oracle::Vec4d x;
        x.t = pt.t;
        x.r = pt.r;
        x.theta = pt.theta;
        x.phi = pt.phi;

        if (!metric.IsValid(x)) continue;

        double g[4][4], g_inv[4][4];
        metric.Evaluate(x, g, g_inv);

        // Compute determinant of 4x4 metric
        // For block-diagonal Kerr metric: det(g) = det(t-φ block) × g_rr × g_θθ
        // det(t-φ block) = g_tt × g_φφ - g_tφ²
        double det_tphi = g[0][0] * g[3][3] - g[0][3] * g[0][3];
        double det_g = det_tphi * g[1][1] * g[2][2];

        EXPECT_GT(std::abs(det_g), DETERMINANT_TOL)
            << "Degenerate metric at r=" << pt.r << ", θ=" << pt.theta << ", det=" << det_g;

        // For Lorentzian signature, determinant should be negative
        EXPECT_LT(det_g, 0) << "Metric determinant should be negative for (-,+,+,+) signature";
    }
}

// =============================================================================
// Boundary Stress Tests
// =============================================================================

// Test: Behavior near horizon (should not produce NaN, may produce large values)
TEST_F(MetricValidationTests, KerrMetricD_NearHorizonStability) {
    sirius::oracle::KerrMetricD metric(1.0, 0.5);
    double r_plus = metric.HorizonRadius();

    // Points approaching horizon
    std::vector<double> buffer_factors = {1.01, 1.005, 1.002, 1.001};

    for (double factor : buffer_factors) {
        sirius::oracle::Vec4d x;
        x.t = 0;
        x.r = r_plus * factor;
        x.theta = math::kHalfPi;
        x.phi = 0;

        double g[4][4], g_inv[4][4];
        metric.Evaluate(x, g, g_inv);

        // Should not be NaN or Inf even close to horizon
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                EXPECT_FALSE(std::isnan(g[mu][nu])) << "NaN near horizon at r/r+ = " << factor;
                EXPECT_FALSE(std::isinf(g[mu][nu])) << "Inf near horizon at r/r+ = " << factor;
            }
        }

        // g_rr should be large but finite
        EXPECT_GT(g[1][1], 1.0) << "g_rr should be > 1 near horizon";
        EXPECT_LT(g[1][1], 1e10) << "g_rr should be bounded near horizon buffer";
    }
}

// Test: Behavior near poles (theta → 0 or π)
TEST_F(MetricValidationTests, KerrMetricD_NearPoleStability) {
    sirius::oracle::KerrMetricD metric(1.0, 0.5);

    // Points approaching poles
    std::vector<double> thetas = {0.01, 0.001, math::kPi - 0.01, math::kPi - 0.001};

    for (double theta : thetas) {
        sirius::oracle::Vec4d x;
        x.t = 0;
        x.r = 10.0;
        x.theta = theta;
        x.phi = 0;

        double g[4][4], g_inv[4][4];
        metric.Evaluate(x, g, g_inv);

        // Should not be NaN
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                EXPECT_FALSE(std::isnan(g[mu][nu])) << "NaN near pole at θ = " << theta;
            }
        }

        // g_φφ should approach 0 at poles (sin²θ → 0)
        if (theta < 0.1 || theta > math::kPi - 0.1) {
            EXPECT_LT(std::abs(g[3][3]), 1.0) << "g_φφ should be small near poles";
        }
    }
}

}  // namespace sirius::test
