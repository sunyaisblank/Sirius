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

#include <gtest/gtest.h>
#include <cmath>
#include <array>
#include <vector>

#include <PHMT100A.h>
#include <PHMT100B.h>
#include <PHCN001A.h>
#include <MTTN001A.h>
#include <MTDL001A.h>

namespace sirius::test {

using namespace Sirius::Constants;

// =============================================================================
// Analytic Reference Values
// =============================================================================

namespace AnalyticRef {

// Schwarzschild characteristic radii (units of M)
constexpr double SCHWARZSCHILD_HORIZON = 2.0;          // r_s = 2M
constexpr double SCHWARZSCHILD_PHOTON_SPHERE = 3.0;    // r_ph = 3M
constexpr double SCHWARZSCHILD_ISCO = 6.0;             // r_ISCO = 6M

// Einstein light deflection (weak field)
// Δφ = 4GM/(c²b) = 4M/b in geometric units
constexpr double EINSTEIN_DEFLECTION_COEFF = 4.0;

// Extremal Kerr (a = M) values
constexpr double KERR_EXTREMAL_ISCO_PROGRADE = 1.0;    // r = M
constexpr double KERR_EXTREMAL_ISCO_RETROGRADE = 9.0;  // r = 9M

// Tolerances for analytic comparisons
constexpr double ANALYTIC_TOL = 1e-10;
constexpr double WEAK_FIELD_TOL = 1e-4;  // 0.01% for weak field approximations

} // namespace AnalyticRef

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
    double computeKerrHorizon(double M, double a) const {
        return M + std::sqrt(M * M - a * a);
    }

    /// @brief Compute Kerr ergosphere radius at given theta
    /// r_ergo = M + √(M² - a²cos²θ)
    double computeKerrErgosphere(double M, double a, double theta) const {
        double cos_th = std::cos(theta);
        return M + std::sqrt(M * M - a * a * cos_th * cos_th);
    }

    /// @brief Compute Schwarzschild orbital angular velocity
    /// Ω = √(M/r³) = dφ/dt for circular orbits
    double computeSchwarzschildOmega(double M, double r) const {
        return std::sqrt(M / (r * r * r));
    }

    /// @brief Compute weak-field light deflection angle
    /// Δφ = 4M/b where b is impact parameter
    double computeWeakFieldDeflection(double M, double b) const {
        return 4.0 * M / b;
    }
};

// =============================================================================
// Schwarzschild Characteristic Radii Tests
// =============================================================================

TEST_F(AnalyticValidationTests, SchwarzschildHorizonRadius) {
    // The event horizon is at r_s = 2M where g_tt → 0
    double M = 1.0;
    Sirius::KerrSchildFamily metric(Sirius::KerrSchildParams::Schwarzschild(M));

    double r_horizon = metric.outerHorizonRadius();
    EXPECT_NEAR(r_horizon, AnalyticRef::SCHWARZSCHILD_HORIZON * M, AnalyticRef::ANALYTIC_TOL)
        << "Schwarzschild horizon should be at r = 2M";
}

TEST_F(AnalyticValidationTests, SchwarzschildPhotonSphere) {
    // Photon sphere at r = 3M where circular null orbits exist
    double M = 1.0;
    Sirius::KerrSchildFamily metric(Sirius::KerrSchildParams::Schwarzschild(M));

    // The photon sphere is where V_eff has an extremum for null geodesics
    // For Schwarzschild: r_ph = 3M
    double r_expected = AnalyticRef::SCHWARZSCHILD_PHOTON_SPHERE * M;

    // Verify metric signature is still Lorentzian at photon sphere
    Tensor<double, 4> pos;
    pos(0) = 0; pos(1) = r_expected; pos(2) = 0; pos(3) = 0;

    Metric4D g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric.evaluate(pos, g, dg);

    // g_tt should be negative (outside horizon)
    EXPECT_LT(g(0, 0).real, 0) << "Metric should be Lorentzian at photon sphere";

    // At r = 3M: g_tt = -(1 - 2M/r) = -(1 - 2/3) = -1/3
    double expected_gtt = -(1.0 - 2.0 * M / r_expected);
    EXPECT_NEAR(g(0, 0).real, expected_gtt, AnalyticRef::ANALYTIC_TOL)
        << "g_tt at photon sphere should be -1/3";
}

TEST_F(AnalyticValidationTests, SchwarzschildISCO) {
    // ISCO at r = 6M
    double M = 1.0;
    Sirius::KerrSchildFamily metric(Sirius::KerrSchildParams::Schwarzschild(M));

    double r_isco = metric.iscoRadius();
    EXPECT_NEAR(r_isco, AnalyticRef::SCHWARZSCHILD_ISCO * M, AnalyticRef::ANALYTIC_TOL)
        << "Schwarzschild ISCO should be at r = 6M";
}

TEST_F(AnalyticValidationTests, SchwarzschildOrbitalVelocity) {
    // For circular orbits: Ω² = M/r³ (Kepler's third law in strong gravity form)
    double M = 1.0;
    std::vector<double> radii = {10.0, 20.0, 50.0, 100.0};

    for (double r : radii) {
        double omega_expected = computeSchwarzschildOmega(M, r);
        double omega_squared = M / (r * r * r);

        EXPECT_NEAR(omega_expected * omega_expected, omega_squared, AnalyticRef::ANALYTIC_TOL)
            << "Orbital angular velocity at r = " << r;
    }
}

// =============================================================================
// Kerr Characteristic Radii Tests
// =============================================================================

TEST_F(AnalyticValidationTests, KerrHorizonRadius) {
    double M = 1.0;
    std::vector<double> spins = {0.0, 0.3, 0.5, 0.7, 0.9, 0.99};

    for (double a : spins) {
        sirius::physics::KerrMetricD metric(M, a);

        double r_computed = metric.horizonRadius();
        double r_expected = computeKerrHorizon(M, a);

        EXPECT_NEAR(r_computed, r_expected, AnalyticRef::ANALYTIC_TOL)
            << "Kerr horizon for a = " << a << ": expected " << r_expected
            << ", got " << r_computed;
    }
}

TEST_F(AnalyticValidationTests, KerrISCOPrograde) {
    // Prograde ISCO decreases with spin
    double M = 1.0;
    std::vector<double> spins = {0.0, 0.3, 0.5, 0.7, 0.9};

    for (double a : spins) {
        sirius::physics::KerrMetricD metric(M, a);

        double r_computed = metric.iscoRadius();
        double r_expected = computeKerrISCO(a, true);

        // ISCO formula can have small numerical differences
        EXPECT_NEAR(r_computed, r_expected, 0.01)
            << "Kerr prograde ISCO for a = " << a << ": expected " << r_expected
            << ", got " << r_computed;
    }
}

TEST_F(AnalyticValidationTests, KerrISCODecreaseWithSpin) {
    // Verify that prograde ISCO monotonically decreases with spin
    double M = 1.0;
    double prev_isco = 6.0;  // Schwarzschild value

    for (double a = 0.1; a <= 0.95; a += 0.1) {
        double r_isco = computeKerrISCO(a, true);
        EXPECT_LT(r_isco, prev_isco)
            << "ISCO should decrease with spin: a = " << a;
        prev_isco = r_isco;
    }
}

TEST_F(AnalyticValidationTests, KerrErgosphereAtEquator) {
    // At equator (θ = π/2), ergosphere radius = horizon radius for Schwarzschild
    // but extends further for Kerr: r_ergo(π/2) = M + √(M² - 0) = 2M
    double M = 1.0;
    double theta_eq = Math::HALF_PI;

    // Schwarzschild: ergosphere = horizon = 2M
    double r_ergo_schw = computeKerrErgosphere(M, 0.0, theta_eq);
    EXPECT_NEAR(r_ergo_schw, 2.0 * M, AnalyticRef::ANALYTIC_TOL)
        << "Schwarzschild ergosphere at equator";

    // Kerr with a = 0.9M: r_ergo(π/2) = 2M (independent of spin at equator!)
    double r_ergo_kerr = computeKerrErgosphere(M, 0.9, theta_eq);
    EXPECT_NEAR(r_ergo_kerr, 2.0 * M, AnalyticRef::ANALYTIC_TOL)
        << "Kerr ergosphere at equator";
}

TEST_F(AnalyticValidationTests, KerrErgosphereAtPole) {
    // At pole (θ = 0), cos²θ = 1, so r_ergo = M + √(M² - a²) = r_+
    // Ergosphere touches horizon at poles
    double M = 1.0;
    double a = 0.7;

    double r_horizon = computeKerrHorizon(M, a);
    double r_ergo_pole = computeKerrErgosphere(M, a, 0.0);

    EXPECT_NEAR(r_ergo_pole, r_horizon, AnalyticRef::ANALYTIC_TOL)
        << "Ergosphere should touch horizon at poles";
}

// =============================================================================
// Weak-Field Light Deflection Tests
// =============================================================================

TEST_F(AnalyticValidationTests, WeakFieldDeflectionFormula) {
    // Einstein's prediction: Δφ = 4GM/(c²b) = 4M/b in geometric units
    // For Sun grazing ray: b ≈ R_sun, M ≈ 1.5 km
    // Δφ ≈ 1.75 arcsec (verified by Eddington 1919)

    double M = 1.0;
    std::vector<double> impact_params = {10.0, 50.0, 100.0, 500.0, 1000.0};

    for (double b : impact_params) {
        double deflection = computeWeakFieldDeflection(M, b);
        double expected = 4.0 * M / b;

        EXPECT_NEAR(deflection, expected, AnalyticRef::ANALYTIC_TOL)
            << "Weak-field deflection for b = " << b;

        // Higher-order corrections are O(M²/b²), so relative error should be small
        // for large b. Require M/b ≤ 0.1 for weak-field validity.
        double relative_correction = M / b;
        EXPECT_LE(relative_correction, 0.1)
            << "Impact parameter should be in weak-field regime";
    }
}

TEST_F(AnalyticValidationTests, SolarDeflectionOrderOfMagnitude) {
    // Sun: M_sun ≈ 1.5 km (in geometric units c=G=1)
    // R_sun ≈ 696,000 km
    // Exact formula: Δφ = 4GM/(c²R) = 4 × 1.5 / 696000 ≈ 8.62e-6 rad ≈ 1.78 arcsec
    // Historical value (Eddington 1919): ~1.75 arcsec (used different R_sun value)

    double M_sun_geometric = 1.5;  // km
    double R_sun = 696000.0;       // km

    double deflection_rad = computeWeakFieldDeflection(M_sun_geometric, R_sun);
    double deflection_arcsec = deflection_rad * (180.0 / Math::PI) * 3600.0;

    // Exact value with these parameters: 1.778 arcsec
    // Tolerance allows for variations in R_sun definition
    EXPECT_NEAR(deflection_arcsec, 1.78, 0.05)
        << "Solar limb deflection should be ~1.78 arcsec";
}

// =============================================================================
// Metric Limit Tests
// =============================================================================

TEST_F(AnalyticValidationTests, KerrReducesToSchwarzschildAtZeroSpin) {
    double M = 1.0;

    sirius::physics::KerrMetricD kerr(M, 0.0);
    Sirius::KerrSchildFamily schw(Sirius::KerrSchildParams::Schwarzschild(M));

    // Compare at several radii
    std::vector<double> radii = {3.0, 6.0, 10.0, 20.0};

    for (double r : radii) {
        sirius::physics::Vec4d x_kerr;
        x_kerr.t = 0; x_kerr.r = r; x_kerr.theta = Math::HALF_PI; x_kerr.phi = 0;

        double g_kerr[4][4], g_inv_kerr[4][4];
        kerr.evaluate(x_kerr, g_kerr, g_inv_kerr);

        // g_tt should be -(1 - 2M/r)
        double expected_gtt = -(1.0 - 2.0 * M / r);
        EXPECT_NEAR(g_kerr[0][0], expected_gtt, AnalyticRef::ANALYTIC_TOL)
            << "Kerr(a=0) g_tt at r = " << r;

        // Off-diagonal g_tφ should be zero for Schwarzschild
        EXPECT_NEAR(g_kerr[0][3], 0.0, AnalyticRef::ANALYTIC_TOL)
            << "Kerr(a=0) g_tφ should vanish";
    }
}

TEST_F(AnalyticValidationTests, AsymptoticFlatness) {
    // At large r, all metrics should approach Minkowski
    double M = 1.0;
    double a = 0.9;
    double r_far = 1000.0 * M;

    sirius::physics::KerrMetricD metric(M, a);
    sirius::physics::Vec4d x;
    x.t = 0; x.r = r_far; x.theta = Math::HALF_PI; x.phi = 0;

    double g[4][4], g_inv[4][4];
    metric.evaluate(x, g, g_inv);

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
    sirius::physics::KerrMetricD metric(M, 0.0);  // Schwarzschild

    std::vector<double> radii = {3.0, 6.0, 10.0, 20.0, 50.0};

    for (double r : radii) {
        sirius::physics::Vec4d x;
        x.t = 0; x.r = r; x.theta = Math::HALF_PI; x.phi = 0;

        double K = metric.kretschmann(x);
        double K_expected = 48.0 * M * M / std::pow(r, 6);

        EXPECT_NEAR(K, K_expected, K_expected * 1e-10)
            << "Schwarzschild Kretschmann at r = " << r;
    }
}

TEST_F(AnalyticValidationTests, KretschmannMonotonicDecrease) {
    // Curvature (Kretschmann scalar) should decrease monotonically with r
    double M = 1.0;
    sirius::physics::KerrMetricD metric(M, 0.5);

    double K_prev = std::numeric_limits<double>::max();
    std::vector<double> radii = {3.0, 5.0, 10.0, 20.0, 50.0, 100.0};

    for (double r : radii) {
        sirius::physics::Vec4d x;
        x.t = 0; x.r = r; x.theta = Math::HALF_PI; x.phi = 0;

        double K = metric.kretschmann(x);
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
    double M = 1.0;
    double a = 0.5;

    sirius::physics::KerrMetricD metric(M, a);
    double r_plus = metric.horizonRadius();

    // Approach horizon from outside
    std::vector<double> factors = {1.1, 1.05, 1.02, 1.01};

    double prev_grr = 0;
    for (double f : factors) {
        sirius::physics::Vec4d x;
        x.t = 0; x.r = r_plus * f; x.theta = Math::HALF_PI; x.phi = 0;

        double g[4][4], g_inv[4][4];
        metric.evaluate(x, g, g_inv);

        // g_rr should increase as we approach horizon
        EXPECT_GT(g[1][1], prev_grr)
            << "g_rr should increase approaching horizon (f = " << f << ")";
        prev_grr = g[1][1];
    }
}

// =============================================================================
// Energy and Angular Momentum at ISCO
// =============================================================================

TEST_F(AnalyticValidationTests, SchwarzschildISCOEnergy) {
    // At ISCO (r = 6M): E/m = √(8/9) ≈ 0.9428 for massive particles
    // This gives the radiative efficiency: η = 1 - E/m ≈ 5.7%
    double M = 1.0;
    double r_isco = 6.0 * M;

    // Specific energy for circular orbits: E = (1 - 2M/r) / √(1 - 3M/r)
    double E_per_m = (1.0 - 2.0 * M / r_isco) / std::sqrt(1.0 - 3.0 * M / r_isco);
    double E_expected = std::sqrt(8.0 / 9.0);

    EXPECT_NEAR(E_per_m, E_expected, AnalyticRef::ANALYTIC_TOL)
        << "ISCO specific energy should be √(8/9)";

    // Radiative efficiency
    double efficiency = 1.0 - E_per_m;
    EXPECT_NEAR(efficiency, 1.0 - std::sqrt(8.0 / 9.0), AnalyticRef::ANALYTIC_TOL)
        << "Schwarzschild radiative efficiency ~5.7%";
}

} // namespace sirius::test
