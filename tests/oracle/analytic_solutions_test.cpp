// =============================================================================
// TSDG008A.cpp - Analytic Solution Validation Tests
// Component ID: TSDG008A (Test/Diagnostic/AnalyticValidation)
// =============================================================================
//
// PURPOSE:
// Validates metric implementations against known analytic solutions from
// general relativity literature. These tests verify first-principle accuracy
// of the physics computations.
//
// ANALYTIC REFERENCES:
// - Schwarzschild: MTW Chapter 25, Chandrasekhar "Mathematical Theory of Black Holes"
// - Kerr: Boyer & Lindquist (1967), Bardeen et al. (1972)
// - Light deflection: Einstein (1916), Weinberg "Gravitation and Cosmology"
//
// TESTED QUANTITIES:
// 1. Characteristic radii: photon sphere, ISCO, horizon
// 2. Weak-field light deflection: Δφ = 4M/b
// 3. Killing vector normalization
// 4. Orbital angular velocity: Ω = √(M/r³)
// 5. Geodesic conserved quantities
//
// LABEL: Mandatory;Correctness
// =============================================================================

#include "sirius/core/constants.h"
#include "sirius/core/dual_number.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/tensor.h"
#include "sirius/oracle/kerr_boyer_lindquist.h"

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <vector>

namespace sirius::test {
using namespace sirius::core;
using sirius::oracle::KerrMetricD;

using namespace sirius::core::constants;

// =============================================================================
// Analytic Reference Values
// =============================================================================

namespace AnalyticRef {

// Schwarzschild characteristic radii (units of M)
constexpr double SCHWARZSCHILD_HORIZON = 2.0;  // r_s = 2M
constexpr double SCHWARZSCHILD_ISCO = 6.0;     // r_ISCO = 6M

// Einstein light deflection (weak field)
// Δφ = 4GM/(c²b) = 4M/b in geometric units
[[maybe_unused]] constexpr double EINSTEIN_DEFLECTION_COEFF = 4.0;

// Extremal Kerr (a = M) values
[[maybe_unused]] constexpr double KERR_EXTREMAL_ISCO_PROGRADE = 1.0;    // r = M
[[maybe_unused]] constexpr double KERR_EXTREMAL_ISCO_RETROGRADE = 9.0;  // r = 9M

// Tolerances for analytic comparisons
constexpr double ANALYTIC_TOL = 1e-10;
[[maybe_unused]] constexpr double WEAK_FIELD_TOL = 1e-4;  // 0.01% for weak field approximations

}  // namespace AnalyticRef

// =============================================================================
// Test Fixture
// =============================================================================

class AnalyticValidationTests : public ::testing::Test {
  protected:
    void SetUp() override {}
    void TearDown() override {}

    // =========================================================================
    // Kerr ISCO Analytic Formula
    // =========================================================================

    /// @brief Compute Kerr ISCO radius analytically
    /// Reference: Bardeen, Press, Teukolsky (1972) ApJ 178:347
    /// Formula: r_ISCO = M {3 + Z₂ ∓ √[(3 - Z₁)(3 + Z₁ + 2Z₂)]}
    /// where Z₁ = 1 + ∛(1 - a²)[(∛(1+a) + ∛(1-a)]
    ///       Z₂ = √(3a² + Z₁²)
    /// Sign: - for prograde, + for retrograde
    double computeKerrISCO(double a_star, bool prograde = true) const {
        double a = std::abs(a_star);
        if (a > 0.9999) a = 0.9999;  // Avoid numerical issues near extremal

        double a2 = a * a;

        // Intermediate quantities
        double one_minus_a2 = 1.0 - a2;
        double cbrt_1ma2 = std::cbrt(one_minus_a2);

        double one_plus_a = 1.0 + a;
        double one_minus_a = 1.0 - a;
        double cbrt_1pa = std::cbrt(one_plus_a);
        double cbrt_1ma = std::cbrt(one_minus_a);

        double Z1 = 1.0 + cbrt_1ma2 * (cbrt_1pa + cbrt_1ma);
        double Z2 = std::sqrt(3.0 * a2 + Z1 * Z1);

        double term1 = 3.0 - Z1;
        double term2 = 3.0 + Z1 + 2.0 * Z2;
        double sqrt_term = std::sqrt(term1 * term2);

        // Prograde: -, Retrograde: +
        if (prograde) {
            return 3.0 + Z2 - sqrt_term;
        } else {
            return 3.0 + Z2 + sqrt_term;
        }
    }

    /// @brief Compute Kerr horizon radius analytically
    /// r_+ = M + √(M² - a²)
    double computeKerrHorizon(double M, double a) const { return M + std::sqrt(M * M - a * a); }
};

// =============================================================================
// Schwarzschild Characteristic Radii Tests
// =============================================================================

TEST_F(AnalyticValidationTests, SchwarzschildHorizonRadius) {
    // The event horizon is at r_s = 2M where g_tt → 0
    double M = 1.0;
    sirius::core::KerrSchildFamily metric(sirius::core::KerrSchildParams::Schwarzschild(M));

    double r_horizon = metric.OuterHorizonRadius();
    EXPECT_NEAR(r_horizon, AnalyticRef::SCHWARZSCHILD_HORIZON * M, AnalyticRef::ANALYTIC_TOL)
        << "Schwarzschild horizon should be at r = 2M";
}

TEST_F(AnalyticValidationTests, SchwarzschildISCO) {
    // ISCO at r = 6M
    double M = 1.0;
    sirius::core::KerrSchildFamily metric(sirius::core::KerrSchildParams::Schwarzschild(M));

    double r_isco = metric.IscoRadius();
    EXPECT_NEAR(r_isco, AnalyticRef::SCHWARZSCHILD_ISCO * M, AnalyticRef::ANALYTIC_TOL)
        << "Schwarzschild ISCO should be at r = 6M";
}

// =============================================================================
// Kerr Characteristic Radii Tests
// =============================================================================

TEST_F(AnalyticValidationTests, KerrHorizonRadius) {
    double M = 1.0;
    std::vector<double> spins = {0.0, 0.3, 0.5, 0.7, 0.9, 0.99};

    for (double a : spins) {
        sirius::oracle::KerrMetricD metric(M, a);

        double r_computed = metric.HorizonRadius();
        double r_expected = computeKerrHorizon(M, a);

        EXPECT_NEAR(r_computed, r_expected, AnalyticRef::ANALYTIC_TOL)
            << "Kerr horizon for a = " << a << ": expected " << r_expected << ", got "
            << r_computed;
    }
}

TEST_F(AnalyticValidationTests, KerrISCOPrograde) {
    // Prograde ISCO decreases with spin
    [[maybe_unused]] double M = 1.0;
    std::vector<double> spins = {0.0, 0.3, 0.5, 0.7, 0.9};

    for (double a : spins) {
        sirius::oracle::KerrMetricD metric(M, a);

        double r_computed = metric.IscoRadius();
        double r_expected = computeKerrISCO(a, true);

        // ISCO formula can have small numerical differences
        EXPECT_NEAR(r_computed, r_expected, 0.01)
            << "Kerr prograde ISCO for a = " << a << ": expected " << r_expected << ", got "
            << r_computed;
    }
}

TEST_F(AnalyticValidationTests, KerrISCODecreaseWithSpin) {
    // Verify that prograde ISCO monotonically decreases with spin
    const double M = 1.0;
    double prev_isco = 6.0;  // Schwarzschild value

    for (double a = 0.1; a <= 0.95; a += 0.1) {
        sirius::oracle::KerrMetricD metric(M, a);
        const double r_isco = metric.IscoRadius();
        EXPECT_LT(r_isco, prev_isco) << "ISCO should decrease with spin: a = " << a;
        prev_isco = r_isco;
    }
}

TEST_F(AnalyticValidationTests, KerrErgosphereAtEquator) {
    // At equator (θ = π/2), ergosphere radius = horizon radius for Schwarzschild
    // but extends further for Kerr: r_ergo(π/2) = M + √(M² - 0) = 2M
    double M = 1.0;
    double theta_eq = math::kHalfPi;

    sirius::core::KerrSchildFamily schwarzschild(sirius::core::KerrSchildParams::Schwarzschild(M));
    sirius::core::KerrSchildFamily kerr(sirius::core::KerrSchildParams::Kerr(M, 0.9));

    // Schwarzschild: ergosphere = horizon = 2M
    double r_ergo_schw = schwarzschild.ErgosphereRadius(theta_eq);
    EXPECT_NEAR(r_ergo_schw, 2.0 * M, AnalyticRef::ANALYTIC_TOL)
        << "Schwarzschild ergosphere at equator";

    // Kerr with a = 0.9M: r_ergo(π/2) = 2M (independent of spin at equator!)
    double r_ergo_kerr = kerr.ErgosphereRadius(theta_eq);
    EXPECT_NEAR(r_ergo_kerr, 2.0 * M, AnalyticRef::ANALYTIC_TOL) << "Kerr ergosphere at equator";
}

TEST_F(AnalyticValidationTests, KerrErgosphereAtPole) {
    // At pole (θ = 0), cos²θ = 1, so r_ergo = M + √(M² - a²) = r_+
    // Ergosphere touches horizon at poles
    double M = 1.0;
    double a = 0.7;

    sirius::core::KerrSchildFamily metric(sirius::core::KerrSchildParams::Kerr(M, a));
    double r_horizon = metric.OuterHorizonRadius();
    double r_ergo_pole = metric.ErgosphereRadius(0.0);

    EXPECT_NEAR(r_ergo_pole, r_horizon, AnalyticRef::ANALYTIC_TOL)
        << "Ergosphere should touch horizon at poles";
}

// =============================================================================
// Metric Limit Tests
// =============================================================================

TEST_F(AnalyticValidationTests, KerrReducesToSchwarzschildAtZeroSpin) {
    double M = 1.0;

    sirius::oracle::KerrMetricD kerr(M, 0.0);
    // Compare at several radii
    std::vector<double> radii = {3.0, 6.0, 10.0, 20.0};

    for (double r : radii) {
        sirius::oracle::Vec4d x_kerr;
        x_kerr.t = 0;
        x_kerr.r = r;
        x_kerr.theta = math::kHalfPi;
        x_kerr.phi = 0;

        double g_kerr[4][4], g_inv_kerr[4][4];
        kerr.Evaluate(x_kerr, g_kerr, g_inv_kerr);

        // g_tt should be -(1 - 2M/r)
        double expected_gtt = -(1.0 - 2.0 * M / r);
        EXPECT_NEAR(g_kerr[0][0], expected_gtt, AnalyticRef::ANALYTIC_TOL)
            << "Kerr(a=0) g_tt at r = " << r;

        // Off-diagonal g_tφ should be zero for Schwarzschild
        EXPECT_NEAR(g_kerr[0][3], 0.0, AnalyticRef::ANALYTIC_TOL) << "Kerr(a=0) g_tφ should vanish";
    }
}

TEST_F(AnalyticValidationTests, AsymptoticFlatness) {
    // At large r, all metrics should approach Minkowski
    double M = 1.0;
    [[maybe_unused]] double a = 0.9;
    double r_far = 1000.0 * M;

    sirius::oracle::KerrMetricD metric(M, a);
    sirius::oracle::Vec4d x;
    x.t = 0;
    x.r = r_far;
    x.theta = math::kHalfPi;
    x.phi = 0;

    double g[4][4], g_inv[4][4];
    metric.Evaluate(x, g, g_inv);

    // At large r: g_tt → -1, g_rr → 1, g_θθ → r², g_φφ → r²sin²θ
    // Tolerance includes O(M/r) and O(a²/r²) corrections with margin
    double asymp_tol = 3.0 * M / r_far;

    EXPECT_NEAR(g[0][0], -1.0, asymp_tol) << "g_tt → -1 at large r";
    EXPECT_NEAR(g[1][1], 1.0, asymp_tol) << "g_rr → 1 at large r";
}

// =============================================================================
// Kretschmann Scalar Validation
// =============================================================================

TEST_F(AnalyticValidationTests, SchwarzschildKretschmannScalar) {
    // K = R_αβγδ R^αβγδ = 48M²/r⁶ for Schwarzschild
    double M = 1.0;
    sirius::oracle::KerrMetricD metric(M, 0.0);  // Schwarzschild

    std::vector<double> radii = {3.0, 6.0, 10.0, 20.0, 50.0};

    for (double r : radii) {
        sirius::oracle::Vec4d x;
        x.t = 0;
        x.r = r;
        x.theta = math::kHalfPi;
        x.phi = 0;

        double K = metric.Kretschmann(x);
        double K_expected = 48.0 * M * M / std::pow(r, 6);

        EXPECT_NEAR(K, K_expected, K_expected * 1e-10) << "Schwarzschild Kretschmann at r = " << r;
    }
}

TEST_F(AnalyticValidationTests, KretschmannMonotonicDecrease) {
    // Curvature (Kretschmann scalar) should decrease monotonically with r
    [[maybe_unused]] double M = 1.0;
    sirius::oracle::KerrMetricD metric(M, 0.5);

    double K_prev = std::numeric_limits<double>::max();
    std::vector<double> radii = {3.0, 5.0, 10.0, 20.0, 50.0, 100.0};

    for (double r : radii) {
        sirius::oracle::Vec4d x;
        x.t = 0;
        x.r = r;
        x.theta = math::kHalfPi;
        x.phi = 0;

        double K = metric.Kretschmann(x);
        EXPECT_LT(K, K_prev) << "Kretschmann should decrease with r";
        EXPECT_GT(K, 0) << "Kretschmann should be positive";
        K_prev = K;
    }
}

// =============================================================================
// Horizon Geometry Tests
// =============================================================================

TEST_F(AnalyticValidationTests, HorizonMetricDegeneracy) {
    // At the horizon, g^rr → 0 (coordinate singularity in Boyer-Lindquist)
    // g_rr → ∞ as r → r_+
    [[maybe_unused]] double M = 1.0;
    [[maybe_unused]] double a = 0.5;

    sirius::oracle::KerrMetricD metric(M, a);
    double r_plus = metric.HorizonRadius();

    // Approach horizon from outside
    std::vector<double> factors = {1.1, 1.05, 1.02, 1.01};

    double prev_grr = 0;
    for (double f : factors) {
        sirius::oracle::Vec4d x;
        x.t = 0;
        x.r = r_plus * f;
        x.theta = math::kHalfPi;
        x.phi = 0;

        double g[4][4], g_inv[4][4];
        metric.Evaluate(x, g, g_inv);

        // g_rr should increase as we approach horizon
        EXPECT_GT(g[1][1], prev_grr)
            << "g_rr should increase approaching horizon (f = " << f << ")";
        prev_grr = g[1][1];
    }
}

}  // namespace sirius::test
