// =============================================================================
// TSDG007A.cpp - Comprehensive Metric Tensor Validation
// Component ID: TSDG007A (Test/Diagnostic/Metric Validation)
// =============================================================================
//
// PURPOSE:
// Validates all metric tensor implementations against mathematical invariants
// defined in docs/specification.md and tolerance constants in PHCN001A.h.
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

#include <gtest/gtest.h>
#include <cmath>
#include <array>
#include <vector>
#include <limits>

#include <PHMT100A.h>
#include <PHMT100B.h>
#include <PHMT101A.h>
#include <PHCN001A.h>
#include <MTTN001A.h>
#include <MTDL001A.h>

namespace sirius::test {

// =============================================================================
// Tolerance Constants (from PHCN001A.h)
// =============================================================================

using namespace Sirius::Constants;

constexpr double SYMMETRY_TOL = Metric::SYMMETRY_TOL;              // 1e-15
constexpr double INVERSE_TOL = Metric::INVERSE_TOL;                // 1e-14
constexpr double CHRISTOFFEL_SYMMETRY_TOL = Metric::CHRISTOFFEL_SYMMETRY_TOL;  // 1e-15
constexpr double DETERMINANT_TOL = Metric::DETERMINANT_TOL;        // 1e-30
constexpr double SIGNATURE_TOL = Metric::SIGNATURE_TOL;            // 1e-10

// =============================================================================
// Test Fixture
// =============================================================================

class MetricValidationTests : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}

    // Sample points for testing (avoiding coordinate singularities)
    struct TestPoint {
        double t, x, y, z;      // Cartesian
        double r, theta, phi;   // Spherical (Boyer-Lindquist)
    };

    std::vector<TestPoint> getSamplePoints(double M = 1.0, double a = 0.0) {
        // Outer horizon for Kerr: r+ = M + sqrt(M² - a²)
        double r_plus = M + std::sqrt(std::max(M*M - a*a, 0.0));
        double r_safe = std::max(r_plus * 1.5, 3.0 * M);  // Safe distance from horizon

        return {
            // General positions (Cartesian: t, x, y, z | Spherical: r, θ, φ)
            {0, r_safe, 0, 0, r_safe, Math::HALF_PI, 0},
            {0, 0, r_safe, 0, r_safe, Math::HALF_PI, Math::HALF_PI},
            {0, 0, 0, r_safe, r_safe, 0.1, 0},  // Near pole but not at pole
            {0, r_safe*0.7, r_safe*0.7, 0, r_safe, Math::HALF_PI, Math::PI/4},
            // Far field
            {0, 100*M, 0, 0, 100*M, Math::HALF_PI, 0},
            {0, 50*M, 50*M, 0, std::sqrt(2)*50*M, Math::HALF_PI, Math::PI/4},
            // Intermediate distances
            {0, 10*M, 0, 0, 10*M, Math::HALF_PI, 0},
            {0, 5*M, 5*M, 5*M, std::sqrt(75)*M, std::acos(5*M/(std::sqrt(75)*M)), Math::PI/4},
        };
    }

    // Boundary test points (stress testing)
    std::vector<TestPoint> getBoundaryPoints(double M = 1.0, double a = 0.0) {
        double r_plus = M + std::sqrt(std::max(M*M - a*a, 0.0));

        return {
            // Very close to horizon (1.001 buffer per PHCN001A)
            {0, r_plus * 1.002, 0, 0, r_plus * 1.002, Math::HALF_PI, 0},
            // Near pole (using POLE_EPSILON from PHCN001A)
            {0, 10*M, 0, 0.01, 10*M, 0.001, 0},
            // Large radius
            {0, 1e5*M, 0, 0, 1e5*M, Math::HALF_PI, 0},
        };
    }
};

// =============================================================================
// Kerr-Schild Family Tests (PHMT100A)
// =============================================================================

// Test: Metric symmetry g_μν = g_νμ
TEST_F(MetricValidationTests, KerrSchild_MetricSymmetry) {
    std::vector<Sirius::KerrSchildParams> configs = {
        Sirius::KerrSchildParams::Minkowski(),
        Sirius::KerrSchildParams::Schwarzschild(1.0),
        Sirius::KerrSchildParams::Kerr(1.0, 0.5),
        Sirius::KerrSchildParams::Kerr(1.0, 0.9),
        Sirius::KerrSchildParams::ReissnerNordstrom(1.0, 0.5),
        Sirius::KerrSchildParams::KerrNewman(1.0, 0.5, 0.3),
    };

    for (const auto& params : configs) {
        Sirius::KerrSchildFamily metric(params);
        auto points = getSamplePoints(params.M, params.a);

        for (const auto& pt : points) {
            Tensor<double, 4> pos;
            pos(0) = pt.t; pos(1) = pt.x; pos(2) = pt.y; pos(3) = pt.z;

            Metric4D g;
            Tensor<Dual<double>, 4, 4, 4> dg;
            metric.evaluate(pos, g, dg);

            double max_asym = 0;
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = mu + 1; nu < 4; ++nu) {
                    double asym = std::abs(g(mu, nu).real - g(nu, mu).real);
                    max_asym = std::max(max_asym, asym);
                }
            }

            EXPECT_LT(max_asym, SYMMETRY_TOL)
                << "Metric asymmetry for " << metric.getName()
                << " at r=" << pt.r << ": " << max_asym;
        }
    }
}

// Test: Lorentzian signature (one negative eigenvalue)
TEST_F(MetricValidationTests, KerrSchild_LorentzianSignature) {
    std::vector<Sirius::KerrSchildParams> configs = {
        Sirius::KerrSchildParams::Schwarzschild(1.0),
        Sirius::KerrSchildParams::Kerr(1.0, 0.5),
        Sirius::KerrSchildParams::Kerr(1.0, 0.9),
    };

    for (const auto& params : configs) {
        Sirius::KerrSchildFamily metric(params);
        auto points = getSamplePoints(params.M, params.a);

        for (const auto& pt : points) {
            Tensor<double, 4> pos;
            pos(0) = pt.t; pos(1) = pt.x; pos(2) = pt.y; pos(3) = pt.z;

            Metric4D g;
            Tensor<Dual<double>, 4, 4, 4> dg;
            metric.evaluate(pos, g, dg);

            // For Kerr-Schild form: g = η + H·l⊗l where η = diag(-1,1,1,1)
            // The signature should remain (-,+,+,+)
            //
            // Simple check: g_tt should be negative (for static observer outside ergoregion)
            // and diagonal spatial components should be positive

            EXPECT_LT(g(0, 0).real, SIGNATURE_TOL)
                << "g_tt should be negative for " << metric.getName()
                << " at r=" << pt.r;

            EXPECT_GT(g(1, 1).real, -SIGNATURE_TOL)
                << "g_xx should be positive for " << metric.getName();
            EXPECT_GT(g(2, 2).real, -SIGNATURE_TOL)
                << "g_yy should be positive for " << metric.getName();
            EXPECT_GT(g(3, 3).real, -SIGNATURE_TOL)
                << "g_zz should be positive for " << metric.getName();
        }
    }
}

// Test: NaN/Inf detection at safe positions
TEST_F(MetricValidationTests, KerrSchild_NoNaNInf) {
    std::vector<Sirius::KerrSchildParams> configs = {
        Sirius::KerrSchildParams::Schwarzschild(1.0),
        Sirius::KerrSchildParams::Kerr(1.0, 0.9),
    };

    for (const auto& params : configs) {
        Sirius::KerrSchildFamily metric(params);
        auto points = getSamplePoints(params.M, params.a);

        for (const auto& pt : points) {
            Tensor<double, 4> pos;
            pos(0) = pt.t; pos(1) = pt.x; pos(2) = pt.y; pos(3) = pt.z;

            Metric4D g;
            Tensor<Dual<double>, 4, 4, 4> dg;
            metric.evaluate(pos, g, dg);

            // Check all metric components
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = 0; nu < 4; ++nu) {
                    EXPECT_FALSE(std::isnan(g(mu, nu).real))
                        << "NaN in g_" << mu << nu << " for " << metric.getName()
                        << " at (" << pt.x << "," << pt.y << "," << pt.z << ")";
                    EXPECT_FALSE(std::isinf(g(mu, nu).real))
                        << "Inf in g_" << mu << nu << " for " << metric.getName();
                }
            }

            // Check all derivative components
            for (int lam = 0; lam < 4; ++lam) {
                for (int mu = 0; mu < 4; ++mu) {
                    for (int nu = 0; nu < 4; ++nu) {
                        EXPECT_FALSE(std::isnan(dg(lam, mu, nu).real))
                            << "NaN in dg_" << lam << mu << nu << " for " << metric.getName();
                        EXPECT_FALSE(std::isinf(dg(lam, mu, nu).real))
                            << "Inf in dg_" << lam << mu << nu << " for " << metric.getName();
                    }
                }
            }
        }
    }
}

// Test: Minkowski limit (M → 0)
TEST_F(MetricValidationTests, KerrSchild_MinkowskiLimit) {
    Sirius::KerrSchildFamily metric(Sirius::KerrSchildParams::Minkowski());

    // At any point, Minkowski metric should be η = diag(-1,1,1,1)
    std::vector<std::array<double, 3>> positions = {
        {10, 0, 0}, {0, 10, 0}, {5, 5, 5}, {100, 0, 0}
    };

    for (const auto& xyz : positions) {
        Tensor<double, 4> pos;
        pos(0) = 0; pos(1) = xyz[0]; pos(2) = xyz[1]; pos(3) = xyz[2];

        Metric4D g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.evaluate(pos, g, dg);

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
    Sirius::KerrSchildFamily metric(Sirius::KerrSchildParams::Schwarzschild(M));

    // At r >> 2M, g_tt ≈ -(1 - 2M/r), g_rr ≈ 1 + 2M/r
    std::vector<double> radii = {100, 500, 1000};

    for (double r : radii) {
        Tensor<double, 4> pos;
        pos(0) = 0; pos(1) = r; pos(2) = 0; pos(3) = 0;

        Metric4D g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.evaluate(pos, g, dg);

        double expected_g_tt = -(1.0 - 2.0 * M / r);
        double tol = 4.0 * M * M / (r * r);  // O(1/r²) corrections

        EXPECT_NEAR(g(0, 0).real, expected_g_tt, tol)
            << "Weak field g_tt error at r=" << r;

        // In Kerr-Schild Cartesian, g_rr ≈ 1 + 2M/r on the x-axis
        // But full tensor form includes off-diagonal corrections
        // Check spatial components are approximately δ_ij
        EXPECT_NEAR(g(1, 1).real, 1.0 + 2.0 * M / r, tol)
            << "Weak field g_xx error at r=" << r;
    }
}

// =============================================================================
// Boyer-Lindquist Kerr Tests (PHMT100B)
// =============================================================================

// Test: Metric symmetry for KerrMetricD
TEST_F(MetricValidationTests, KerrMetricD_MetricSymmetry) {
    std::vector<double> spins = {0.0, 0.3, 0.6, 0.9};

    for (double a : spins) {
        sirius::physics::KerrMetricD metric(1.0, a);
        auto points = getSamplePoints(1.0, a);

        for (const auto& pt : points) {
            sirius::physics::Vec4d x;
            x.t = pt.t; x.r = pt.r; x.theta = pt.theta; x.phi = pt.phi;

            double g[4][4], g_inv[4][4];
            metric.evaluate(x, g, g_inv);

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
        sirius::physics::KerrMetricD metric(1.0, a);
        auto points = getSamplePoints(1.0, a);

        for (const auto& pt : points) {
            sirius::physics::Vec4d x;
            x.t = pt.t; x.r = pt.r; x.theta = pt.theta; x.phi = pt.phi;

            double g[4][4], g_inv[4][4];
            metric.evaluate(x, g, g_inv);

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

            EXPECT_LT(max_error, INVERSE_TOL)
                << "Inverse metric error for Kerr (a=" << a << ") at r=" << pt.r
                << ": " << max_error;
        }
    }
}

// Test: Christoffel symbol symmetry Γ^λ_μν = Γ^λ_νμ
TEST_F(MetricValidationTests, KerrMetricD_ChristoffelSymmetry) {
    std::vector<double> spins = {0.0, 0.5, 0.9};

    for (double a : spins) {
        sirius::physics::KerrMetricD metric(1.0, a);
        auto points = getSamplePoints(1.0, a);

        for (const auto& pt : points) {
            sirius::physics::Vec4d x;
            x.t = pt.t; x.r = pt.r; x.theta = pt.theta; x.phi = pt.phi;

            double Gamma[4][4][4];
            metric.christoffel(x, Gamma);

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
                << "Christoffel asymmetry for Kerr (a=" << a << ") at r=" << pt.r
                << ": " << max_asym;
        }
    }
}

// Test: Riemann tensor antisymmetry R^α_βγδ = -R^α_βδγ
TEST_F(MetricValidationTests, KerrMetricD_RiemannAntisymmetry) {
    std::vector<double> spins = {0.0, 0.5, 0.9};

    for (double a : spins) {
        sirius::physics::KerrMetricD metric(1.0, a);
        auto points = getSamplePoints(1.0, a);

        for (const auto& pt : points) {
            sirius::physics::Vec4d x;
            x.t = pt.t; x.r = pt.r; x.theta = pt.theta; x.phi = pt.phi;

            double violation = metric.verifyRiemannSymmetries(x);

            EXPECT_LT(violation, 1e-12)
                << "Riemann antisymmetry violation for Kerr (a=" << a
                << ") at r=" << pt.r << ": " << violation;
        }
    }
}

// Test: Kretschmann scalar (analytic validation)
// Schwarzschild: K = 48M²/r⁶
// Kerr: K = 48M²(r² - a²cos²θ)[(r² - a²cos²θ)² - 16a²r²cos²θ] / Σ⁶
// Reference: Henry, R.C. (2000), Astrophys. J. 535:350-353
TEST_F(MetricValidationTests, KerrMetricD_KretschmannScalar) {
    // Test 1: Schwarzschild exact formula K = 48M²/r⁶
    sirius::physics::KerrMetricD schw(1.0, 0.0);

    std::vector<double> radii = {3.0, 5.0, 10.0, 20.0};

    for (double r : radii) {
        sirius::physics::Vec4d x;
        x.t = 0; x.r = r; x.theta = Math::HALF_PI; x.phi = 0;

        double K = schw.kretschmann(x);
        double expected = 48.0 / std::pow(r, 6);  // M=1

        double rel_error = std::abs(K - expected) / expected;

        EXPECT_LT(rel_error, 1e-10)
            << "Schwarzschild Kretschmann error at r=" << r
            << ": computed=" << K << ", expected=" << expected;
    }

    // Test 2: Kerr at equator (θ = π/2, so cos²θ = 0)
    // At equator: K = 48M² × r² × r⁴ / (r²)⁶ = 48M² × r⁶ / r¹² = 48M²/r⁶
    // Same as Schwarzschild!
    sirius::physics::KerrMetricD kerr(1.0, 0.9);

    for (double r : radii) {
        sirius::physics::Vec4d x;
        x.t = 0; x.r = r; x.theta = Math::HALF_PI; x.phi = 0;

        double K = kerr.kretschmann(x);
        double expected = 48.0 / std::pow(r, 6);  // Same as Schwarzschild at equator

        double rel_error = std::abs(K - expected) / expected;

        EXPECT_LT(rel_error, 1e-10)
            << "Kerr equatorial Kretschmann error at r=" << r
            << ": computed=" << K << ", expected=" << expected;
    }

    // Test 3: Kerr off-equator (θ = π/4)
    // Full formula applies with non-zero cos²θ
    double theta = Math::PI / 4;
    double cos2th = 0.5;  // cos²(π/4) = 0.5
    double a = 0.5;
    double a2 = a * a;

    sirius::physics::KerrMetricD kerr05(1.0, a);

    for (double r : radii) {
        sirius::physics::Vec4d x;
        x.t = 0; x.r = r; x.theta = theta; x.phi = 0;

        double r2 = r * r;
        double Sigma = r2 + a2 * cos2th;
        double Sigma6 = std::pow(Sigma, 6);
        double r_term = r2 - a2 * cos2th;
        double bracket = r_term * r_term - 16.0 * a2 * r2 * cos2th;
        double expected = 48.0 * r_term * bracket / Sigma6;  // M=1

        double K = kerr05.kretschmann(x);

        double rel_error = std::abs(K - expected) / std::max(std::abs(expected), 1e-30);

        EXPECT_LT(rel_error, 1e-10)
            << "Kerr off-axis Kretschmann error at r=" << r
            << ": computed=" << K << ", expected=" << expected;
    }

    // Test 4: Kretschmann decreases with radius (monotonicity)
    double K_prev = std::numeric_limits<double>::max();
    for (double r : radii) {
        sirius::physics::Vec4d x;
        x.t = 0; x.r = r; x.theta = Math::HALF_PI; x.phi = 0;

        double K = schw.kretschmann(x);
        EXPECT_LT(K, K_prev) << "Kretschmann should decrease with r";
        K_prev = K;
    }

    // Test 5: At large radius, Kretschmann approaches zero (asymptotic flatness)
    sirius::physics::Vec4d x_far;
    x_far.t = 0; x_far.r = 1000.0; x_far.theta = Math::HALF_PI; x_far.phi = 0;
    double K_far = schw.kretschmann(x_far);
    double expected_far = 48.0 / std::pow(1000.0, 6);  // ~4.8e-17
    EXPECT_NEAR(K_far, expected_far, 1e-20)
        << "Kretschmann should match 48/r⁶ at large radius";
}

// Test: Horizon radius calculation
TEST_F(MetricValidationTests, KerrMetricD_HorizonRadius) {
    // Schwarzschild: r+ = 2M
    sirius::physics::KerrMetricD schw(1.0, 0.0);
    EXPECT_NEAR(schw.horizonRadius(), 2.0, 1e-14);

    // Kerr with a = 0.5: r+ = 1 + sqrt(1 - 0.25) = 1 + sqrt(0.75)
    sirius::physics::KerrMetricD kerr05(1.0, 0.5);
    double expected_05 = 1.0 + std::sqrt(0.75);
    EXPECT_NEAR(kerr05.horizonRadius(), expected_05, 1e-14);

    // Near-extremal Kerr a ≈ M: r+ ≈ M + sqrt(1 - a²)
    // Note: KerrMetricD constructor clamps spin to 0.9999*M for numerical stability
    sirius::physics::KerrMetricD kerr_ext(1.0, 0.9999);
    double expected_ext = 1.0 + std::sqrt(1.0 - 0.9999*0.9999);  // ~1.0141
    EXPECT_NEAR(kerr_ext.horizonRadius(), expected_ext, 1e-10);
}

// Test: ISCO radius calculation
TEST_F(MetricValidationTests, KerrMetricD_ISCORadius) {
    // Schwarzschild: r_ISCO = 6M
    sirius::physics::KerrMetricD schw(1.0, 0.0);
    EXPECT_NEAR(schw.iscoRadius(), 6.0, 1e-10);

    // Prograde extremal Kerr: r_ISCO = M
    sirius::physics::KerrMetricD kerr099(1.0, 0.999);
    EXPECT_LT(kerr099.iscoRadius(), 2.0);  // Should be close to M
}

// Test: No NaN/Inf in valid domain
TEST_F(MetricValidationTests, KerrMetricD_NoNaNInf) {
    std::vector<double> spins = {0.0, 0.5, 0.9};

    for (double a : spins) {
        sirius::physics::KerrMetricD metric(1.0, a);
        auto points = getSamplePoints(1.0, a);

        for (const auto& pt : points) {
            sirius::physics::Vec4d x;
            x.t = pt.t; x.r = pt.r; x.theta = pt.theta; x.phi = pt.phi;

            if (!metric.isValid(x)) continue;  // Skip invalid points

            double g[4][4], g_inv[4][4];
            metric.evaluate(x, g, g_inv);

            double Gamma[4][4][4];
            metric.christoffel(x, Gamma);

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
    std::vector<Sirius::MorrisThorneParams> configs = {
        Sirius::MorrisThorneParams::Ellis(1.0),
        Sirius::MorrisThorneParams::ZeroTidal(1.0),
    };

    for (const auto& params : configs) {
        Sirius::MorrisThorneFamily metric(params);

        // Test points outside throat (r > b0)
        std::vector<double> radii = {1.5, 2.0, 5.0, 10.0, 100.0};

        for (double r : radii) {
            Tensor<double, 4> pos;
            pos(0) = 0;        // t
            pos(1) = r;        // r
            pos(2) = Math::HALF_PI;  // θ = π/2
            pos(3) = 0;        // φ

            Metric4D g;
            Tensor<Dual<double>, 4, 4, 4> dg;
            metric.evaluate(pos, g, dg);

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
                << "Non-diagonal component for " << metric.getName()
                << " at r=" << r;
        }
    }
}

// Test: Flare-out condition (traversability)
TEST_F(MetricValidationTests, MorrisThorne_FlareOutCondition) {
    // Ellis wormhole should satisfy flare-out
    Sirius::MorrisThorneFamily ellis(Sirius::MorrisThorneParams::Ellis(1.0));
    EXPECT_TRUE(ellis.validateFlareOutCondition())
        << "Ellis wormhole should satisfy flare-out condition";

    // Zero-tidal wormhole should satisfy flare-out
    Sirius::MorrisThorneFamily zero_tidal(Sirius::MorrisThorneParams::ZeroTidal(1.0));
    EXPECT_TRUE(zero_tidal.validateFlareOutCondition())
        << "Zero-tidal wormhole should satisfy flare-out condition";
}

// Test: Shape function at throat b(b0) = b0
TEST_F(MetricValidationTests, MorrisThorne_ThroatCondition) {
    double b0 = 2.0;

    // Ellis: b(r) = b0²/r → b(b0) = b0²/b0 = b0
    Sirius::MorrisThorneFamily ellis(Sirius::MorrisThorneParams::Ellis(b0));
    EXPECT_NEAR(ellis.shapeFunction(b0), b0, 1e-14)
        << "Ellis throat condition violated";

    // Zero-tidal: b(r) = b0 → b(b0) = b0
    Sirius::MorrisThorneFamily zero_tidal(Sirius::MorrisThorneParams::ZeroTidal(b0));
    EXPECT_NEAR(zero_tidal.shapeFunction(b0), b0, 1e-14)
        << "Zero-tidal throat condition violated";
}

// Test: No NaN/Inf outside throat
TEST_F(MetricValidationTests, MorrisThorne_NoNaNInf) {
    Sirius::MorrisThorneFamily metric(Sirius::MorrisThorneParams::Ellis(1.0));

    std::vector<double> radii = {1.01, 1.1, 2.0, 10.0, 100.0};
    std::vector<double> thetas = {0.1, Math::HALF_PI, Math::PI - 0.1};

    for (double r : radii) {
        for (double theta : thetas) {
            Tensor<double, 4> pos;
            pos(0) = 0; pos(1) = r; pos(2) = theta; pos(3) = 0;

            Metric4D g;
            Tensor<Dual<double>, 4, 4, 4> dg;
            metric.evaluate(pos, g, dg);

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
    Sirius::MorrisThorneFamily metric(Sirius::MorrisThorneParams::Ellis(1.0));

    std::vector<double> radii = {1.1, 2.0, 10.0};

    for (double r : radii) {
        Tensor<double, 4> pos;
        pos(0) = 0; pos(1) = r; pos(2) = Math::HALF_PI; pos(3) = 0;

        Metric4D g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.evaluate(pos, g, dg);

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
    Sirius::KerrSchildFamily ks_metric(Sirius::KerrSchildParams::Schwarzschild(M));
    sirius::physics::KerrMetricD bl_metric(M, 0.0);

    // Test at equatorial positions
    std::vector<double> radii = {3.0, 6.0, 10.0, 20.0};

    for (double r : radii) {
        // Kerr-Schild evaluation (Cartesian: x = r, y = z = 0)
        Tensor<double, 4> pos_ks;
        pos_ks(0) = 0; pos_ks(1) = r; pos_ks(2) = 0; pos_ks(3) = 0;

        Metric4D g_ks;
        Tensor<Dual<double>, 4, 4, 4> dg_ks;
        ks_metric.evaluate(pos_ks, g_ks, dg_ks);

        // Boyer-Lindquist evaluation
        sirius::physics::Vec4d pos_bl;
        pos_bl.t = 0; pos_bl.r = r; pos_bl.theta = Math::HALF_PI; pos_bl.phi = 0;

        double g_bl[4][4], g_inv_bl[4][4];
        bl_metric.evaluate(pos_bl, g_bl, g_inv_bl);

        // For Schwarzschild on equator with φ=0:
        // g_tt should match
        double g_tt_ks = g_ks(0, 0).real;
        double g_tt_bl = g_bl[0][0];

        // Both should equal -(1 - 2M/r)
        double expected_g_tt = -(1.0 - 2.0 * M / r);

        EXPECT_NEAR(g_tt_ks, expected_g_tt, 1e-10)
            << "KS g_tt error at r=" << r;
        EXPECT_NEAR(g_tt_bl, expected_g_tt, 1e-10)
            << "BL g_tt error at r=" << r;
    }
}

// =============================================================================
// Determinant Non-Degeneracy Tests
// =============================================================================

// Test: Metric determinant is non-zero in valid domain
TEST_F(MetricValidationTests, KerrMetricD_DeterminantNonZero) {
    sirius::physics::KerrMetricD metric(1.0, 0.5);
    auto points = getSamplePoints(1.0, 0.5);

    for (const auto& pt : points) {
        sirius::physics::Vec4d x;
        x.t = pt.t; x.r = pt.r; x.theta = pt.theta; x.phi = pt.phi;

        if (!metric.isValid(x)) continue;

        double g[4][4], g_inv[4][4];
        metric.evaluate(x, g, g_inv);

        // Compute determinant of 4x4 metric
        // For block-diagonal Kerr metric: det(g) = det(t-φ block) × g_rr × g_θθ
        // det(t-φ block) = g_tt × g_φφ - g_tφ²
        double det_tphi = g[0][0] * g[3][3] - g[0][3] * g[0][3];
        double det_g = det_tphi * g[1][1] * g[2][2];

        EXPECT_GT(std::abs(det_g), DETERMINANT_TOL)
            << "Degenerate metric at r=" << pt.r << ", θ=" << pt.theta
            << ", det=" << det_g;

        // For Lorentzian signature, determinant should be negative
        EXPECT_LT(det_g, 0)
            << "Metric determinant should be negative for (-,+,+,+) signature";
    }
}

// =============================================================================
// Boundary Stress Tests
// =============================================================================

// Test: Behavior near horizon (should not produce NaN, may produce large values)
TEST_F(MetricValidationTests, KerrMetricD_NearHorizonStability) {
    sirius::physics::KerrMetricD metric(1.0, 0.5);
    double r_plus = metric.horizonRadius();

    // Points approaching horizon
    std::vector<double> buffer_factors = {1.01, 1.005, 1.002, 1.001};

    for (double factor : buffer_factors) {
        sirius::physics::Vec4d x;
        x.t = 0;
        x.r = r_plus * factor;
        x.theta = Math::HALF_PI;
        x.phi = 0;

        double g[4][4], g_inv[4][4];
        metric.evaluate(x, g, g_inv);

        // Should not be NaN or Inf even close to horizon
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                EXPECT_FALSE(std::isnan(g[mu][nu]))
                    << "NaN near horizon at r/r+ = " << factor;
                EXPECT_FALSE(std::isinf(g[mu][nu]))
                    << "Inf near horizon at r/r+ = " << factor;
            }
        }

        // g_rr should be large but finite
        EXPECT_GT(g[1][1], 1.0)
            << "g_rr should be > 1 near horizon";
        EXPECT_LT(g[1][1], 1e10)
            << "g_rr should be bounded near horizon buffer";
    }
}

// Test: Behavior near poles (theta → 0 or π)
TEST_F(MetricValidationTests, KerrMetricD_NearPoleStability) {
    sirius::physics::KerrMetricD metric(1.0, 0.5);

    // Points approaching poles
    std::vector<double> thetas = {0.01, 0.001, Math::PI - 0.01, Math::PI - 0.001};

    for (double theta : thetas) {
        sirius::physics::Vec4d x;
        x.t = 0;
        x.r = 10.0;
        x.theta = theta;
        x.phi = 0;

        double g[4][4], g_inv[4][4];
        metric.evaluate(x, g, g_inv);

        // Should not be NaN
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                EXPECT_FALSE(std::isnan(g[mu][nu]))
                    << "NaN near pole at θ = " << theta;
            }
        }

        // g_φφ should approach 0 at poles (sin²θ → 0)
        if (theta < 0.1 || theta > Math::PI - 0.1) {
            EXPECT_LT(std::abs(g[3][3]), 1.0)
                << "g_φφ should be small near poles";
        }
    }
}

} // namespace sirius::test
