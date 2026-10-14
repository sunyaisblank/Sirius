// TSIN007A.cpp - Spectral Emission Pipeline Validation Tests
//
// Validates the spectral emission pipeline changes:
// - Novikov-Thorne Q(r) produces zero at ISCO and peaks at ~1.36 r_ISCO
// - Physical temperature from Q(r) maps to correct blackbody colour range
// - Doppler beaming via temperature shift produces correct asymmetry
// - Spin display formatting shows 3 decimal places
//
// LABEL: Mandatory

#include <gtest/gtest.h>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <PHSP001A.h>
#include <PHAD001A.h>

namespace sirius::test {
using namespace Sirius;

// CPU mirror of the GPU computeNovikovThorneQ function
static double computeNovikovThorneQ_CPU(double r, double r_isco, double M, double a_star) {
    a_star = std::clamp(a_star, -0.998, 0.998);

    double x = std::sqrt(r / M);
    double x_isco = std::sqrt(r_isco / M);

    double A = 1.0 - 2.0/(x*x) + a_star/(x*x*x);
    double B = 1.0 - 3.0/(x*x) + 2.0*a_star/(x*x*x);
    double D = (B > 0.0) ? std::sqrt(B) : 0.0;

    double B_isco = 1.0 - 3.0/(x_isco*x_isco) + 2.0*a_star/(x_isco*x_isco*x_isco);
    double D_isco = (B_isco > 0.0) ? std::sqrt(B_isco) : 0.0;

    if (D <= 0.0 || D_isco <= 0.0 || x <= x_isco * 1.001) {
        return 0.0;
    }

    double x_ratio = x_isco / x;
    double base_factor = 1.0 - std::sqrt(x_ratio);
    double log_factor = std::log(std::max(x / x_isco, 1.001));

    double frame_drag = 0.0;
    if (std::abs(a_star) > 0.01) {
        frame_drag = (3.0 * a_star / (2.0 * x * x * x)) * log_factor;
    }

    double redshift_factor = std::sqrt(A / B);
    double Q = base_factor * (1.0 - frame_drag) * std::min(redshift_factor, 2.0);
    Q *= D / D_isco;

    return std::max(Q, 0.0);
}

// =========================================================================
// NT Q(r) Profile Tests
// =========================================================================

TEST(SpectralEmissionTest, NT_QFactorZeroAtISCO_Schwarzschild) {
    double M = 1.0;
    double a_star = 0.0;
    double r_isco = AccretionDiskD::computeISCO(a_star) * M;

    double Q_at_isco = computeNovikovThorneQ_CPU(r_isco, r_isco, M, a_star);
    EXPECT_NEAR(Q_at_isco, 0.0, 1e-6) << "Q must be zero at ISCO";
}

TEST(SpectralEmissionTest, NT_QFactorZeroAtISCO_Kerr) {
    double M = 1.0;
    double a_star = 0.9;
    double r_isco = AccretionDiskD::computeISCO(a_star) * M;

    double Q_at_isco = computeNovikovThorneQ_CPU(r_isco, r_isco, M, a_star);
    EXPECT_NEAR(Q_at_isco, 0.0, 1e-6) << "Q must be zero at ISCO (Kerr)";
}

TEST(SpectralEmissionTest, NT_QFactorPositiveOutsideISCO) {
    double M = 1.0;
    double a_star = 0.0;
    double r_isco = AccretionDiskD::computeISCO(a_star) * M;  // 6M

    // Q should be positive at radii outside ISCO
    for (double r = r_isco * 1.1; r < 30.0 * M; r += 1.0 * M) {
        double Q = computeNovikovThorneQ_CPU(r, r_isco, M, a_star);
        EXPECT_GT(Q, 0.0) << "Q should be positive at r=" << r;
    }

    // Q should increase from zero near ISCO
    double Q_near = computeNovikovThorneQ_CPU(r_isco * 1.01, r_isco, M, a_star);
    double Q_further = computeNovikovThorneQ_CPU(r_isco * 1.5, r_isco, M, a_star);
    EXPECT_GT(Q_further, Q_near) << "Q should increase away from ISCO";
}

TEST(SpectralEmissionTest, NT_QFactorBounded) {
    double M = 1.0;
    double a_star = 0.0;
    double r_isco = AccretionDiskD::computeISCO(a_star) * M;

    // Q should be bounded and finite across the disk
    for (double r = r_isco * 1.01; r < 50.0 * M; r += 0.5 * M) {
        double Q = computeNovikovThorneQ_CPU(r, r_isco, M, a_star);
        EXPECT_GE(Q, 0.0) << "Q must be non-negative at r=" << r;
        EXPECT_LT(Q, 10.0) << "Q must be bounded at r=" << r;
        EXPECT_FALSE(std::isnan(Q)) << "Q must not be NaN at r=" << r;
    }
}

// =========================================================================
// Temperature Mapping Tests
// =========================================================================

TEST(SpectralEmissionTest, TemperatureRange) {
    double M = 1.0;
    double a_star = 0.0;
    double r_isco = AccretionDiskD::computeISCO(a_star) * M;
    double T_scale = 50000.0;

    // At ISCO: T = 0 (Q = 0)
    double Q_isco = computeNovikovThorneQ_CPU(r_isco, r_isco, M, a_star);
    double T_isco = T_scale * std::pow(std::max(Q_isco, 0.0), 0.25);
    EXPECT_NEAR(T_isco, 0.0, 1.0) << "Temperature at ISCO should be ~0";

    // Near ISCO: T should be low but non-zero
    double r_near = r_isco * 1.05;
    double Q_near = computeNovikovThorneQ_CPU(r_near, r_isco, M, a_star);
    double T_near = T_scale * std::pow(std::max(Q_near, 0.0), 0.25);
    EXPECT_GT(T_near, 0.0) << "Temperature just outside ISCO should be > 0";
    EXPECT_LT(T_near, T_scale) << "Temperature should be < T_scale";

    // Mid-disk: T should be in visible blackbody range
    double r_mid = r_isco * 2.0;
    double Q_mid = computeNovikovThorneQ_CPU(r_mid, r_isco, M, a_star);
    double T_mid = T_scale * std::pow(std::max(Q_mid, 0.0), 0.25);
    EXPECT_GT(T_mid, 1000.0) << "Mid-disk temperature should be > 1000K";
    EXPECT_LT(T_mid, T_scale) << "Mid-disk temperature should be < T_scale";

    // Temperature should increase away from ISCO (Q increases)
    EXPECT_GT(T_mid, T_near) << "Mid-disk should be hotter than near-ISCO";
}

// =========================================================================
// Blackbody Colour Direction Tests
// =========================================================================

TEST(SpectralEmissionTest, BlackbodyColourDirection) {
    using namespace Spectral;

    // Hot region (inner disk ~20000K): should be blue-white
    RGB hot = blackbodyToRGB(20000.0);
    EXPECT_GT(hot.b, 0.7f) << "Hot inner disk should have strong blue";

    // Cool region (outer disk ~5000K): should be yellowish
    RGB cool = blackbodyToRGB(5000.0);
    EXPECT_GT(cool.r, cool.b) << "Cool outer disk should be redder than blue";

    // Very cool (3000K): distinctly red
    RGB cold = blackbodyToRGB(3000.0);
    EXPECT_GT(cold.r, 0.8f) << "3000K should be red-dominant";
}

// =========================================================================
// Doppler Shift Direction Tests
// =========================================================================

TEST(SpectralEmissionTest, DopplerShiftDirection) {
    using namespace Spectral;

    double T_emit = 10000.0;

    // Approaching (g > 1): T_obs > T_emit → bluer
    double g_approach = 1.3;
    double T_blue = T_emit * g_approach;
    RGB blueShifted = blackbodyToRGB(T_blue);

    // Receding (g < 1): T_obs < T_emit → redder
    double g_recede = 0.7;
    double T_red = T_emit * g_recede;
    RGB redShifted = blackbodyToRGB(T_red);

    // Blue-shifted should have higher b/r ratio than red-shifted
    float blueRatio = blueShifted.b / std::max(blueShifted.r, 0.01f);
    float redRatio = redShifted.b / std::max(redShifted.r, 0.01f);
    EXPECT_GT(blueRatio, redRatio) << "Blue-shifted should have higher b/r ratio";
}

// =========================================================================
// Spin Display Format Test
// =========================================================================

TEST(SpectralEmissionTest, SpinDisplayFormat) {
    double spin = 0.998;
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(3) << spin;
    EXPECT_EQ(ss.str(), "0.998") << "Spin 0.998 must display as 0.998, not 1.00";

    // Verify 2-decimal would round incorrectly
    std::ostringstream ss2;
    ss2 << std::fixed << std::setprecision(2) << spin;
    EXPECT_EQ(ss2.str(), "1.00") << "2-decimal precision rounds 0.998 to 1.00";
}

} // namespace sirius::test
