// Spectral emission pipeline validation. Ported from TSIN007A.cpp.
//
// Validates the disk/blackbody emission pipeline the render path depends on:
// Novikov-Thorne Q(r) zero at the ISCO and rising outward, physical temperature
// mapping to a plausible blackbody colour, and Doppler-shift colour direction.
// Includes are core-only (blackbody, disk); assertions/tolerances unchanged.

#include "sirius/core/disk/novikov_thorne_disk.h"
#include "sirius/core/spectral/blackbody.h"
#include "sirius/core/spectral/colour_modes.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace {

using sirius::core::AccretionDiskD;

// CPU mirror of the GPU computeNovikovThorneQ function.
double computeNovikovThorneQ_CPU(double r, double r_isco, double M, double a_star) {
    a_star = std::clamp(a_star, -0.998, 0.998);

    double x = std::sqrt(r / M);
    double x_isco = std::sqrt(r_isco / M);

    double A = 1.0 - 2.0 / (x * x) + a_star / (x * x * x);
    double B = 1.0 - 3.0 / (x * x) + 2.0 * a_star / (x * x * x);
    double D = (B > 0.0) ? std::sqrt(B) : 0.0;

    double B_isco = 1.0 - 3.0 / (x_isco * x_isco) + 2.0 * a_star / (x_isco * x_isco * x_isco);
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

}  // namespace

// =========================================================================
// NT Q(r) profile tests.
// =========================================================================

TEST(SpectralEmissionTest, NT_QFactorZeroAtISCO_Schwarzschild) {
    double M = 1.0;
    double a_star = 0.0;
    double r_isco = AccretionDiskD::ComputeIsco(a_star) * M;

    double Q_at_isco = computeNovikovThorneQ_CPU(r_isco, r_isco, M, a_star);
    EXPECT_NEAR(Q_at_isco, 0.0, 1e-6) << "Q must be zero at ISCO";
}

TEST(SpectralEmissionTest, NT_QFactorZeroAtISCO_Kerr) {
    double M = 1.0;
    double a_star = 0.9;
    double r_isco = AccretionDiskD::ComputeIsco(a_star) * M;

    double Q_at_isco = computeNovikovThorneQ_CPU(r_isco, r_isco, M, a_star);
    EXPECT_NEAR(Q_at_isco, 0.0, 1e-6) << "Q must be zero at ISCO (Kerr)";
}

TEST(SpectralEmissionTest, NT_QFactorPositiveOutsideISCO) {
    double M = 1.0;
    double a_star = 0.0;
    double r_isco = AccretionDiskD::ComputeIsco(a_star) * M;  // 6M.

    for (double r = r_isco * 1.1; r < 30.0 * M; r += 1.0 * M) {
        double Q = computeNovikovThorneQ_CPU(r, r_isco, M, a_star);
        EXPECT_GT(Q, 0.0) << "Q should be positive at r=" << r;
    }

    double Q_near = computeNovikovThorneQ_CPU(r_isco * 1.01, r_isco, M, a_star);
    double Q_further = computeNovikovThorneQ_CPU(r_isco * 1.5, r_isco, M, a_star);
    EXPECT_GT(Q_further, Q_near) << "Q should increase away from ISCO";
}

TEST(SpectralEmissionTest, NT_QFactorBounded) {
    double M = 1.0;
    double a_star = 0.0;
    double r_isco = AccretionDiskD::ComputeIsco(a_star) * M;

    for (double r = r_isco * 1.01; r < 50.0 * M; r += 0.5 * M) {
        double Q = computeNovikovThorneQ_CPU(r, r_isco, M, a_star);
        EXPECT_GE(Q, 0.0) << "Q must be non-negative at r=" << r;
        EXPECT_LT(Q, 10.0) << "Q must be bounded at r=" << r;
        EXPECT_FALSE(std::isnan(Q)) << "Q must not be NaN at r=" << r;
    }
}

// =========================================================================
// Temperature mapping tests.
// =========================================================================

TEST(SpectralEmissionTest, TemperatureRange) {
    double M = 1.0;
    double a_star = 0.0;
    double r_isco = AccretionDiskD::ComputeIsco(a_star) * M;
    double T_scale = 50000.0;

    double Q_isco = computeNovikovThorneQ_CPU(r_isco, r_isco, M, a_star);
    double T_isco = T_scale * std::pow(std::max(Q_isco, 0.0), 0.25);
    EXPECT_NEAR(T_isco, 0.0, 1.0) << "Temperature at ISCO should be ~0";

    double r_near = r_isco * 1.05;
    double Q_near = computeNovikovThorneQ_CPU(r_near, r_isco, M, a_star);
    double T_near = T_scale * std::pow(std::max(Q_near, 0.0), 0.25);
    EXPECT_GT(T_near, 0.0) << "Temperature just outside ISCO should be > 0";
    EXPECT_LT(T_near, T_scale) << "Temperature should be < T_scale";

    double r_mid = r_isco * 2.0;
    double Q_mid = computeNovikovThorneQ_CPU(r_mid, r_isco, M, a_star);
    double T_mid = T_scale * std::pow(std::max(Q_mid, 0.0), 0.25);
    EXPECT_GT(T_mid, 1000.0) << "Mid-disk temperature should be > 1000K";
    EXPECT_LT(T_mid, T_scale) << "Mid-disk temperature should be < T_scale";

    EXPECT_GT(T_mid, T_near) << "Mid-disk should be hotter than near-ISCO";
}

// =========================================================================
// Blackbody colour direction tests.
// =========================================================================

TEST(SpectralEmissionTest, BlackbodyColourDirection) {
    using namespace sirius::core::spectral;

    // Hot region (inner disk ~20000K): blue-white.
    Rgb hot = BlackbodyToRgb(20000.0);
    EXPECT_GT(hot.b, 0.7f) << "Hot inner disk should have strong blue";

    // Cool region (outer disk ~5000K): yellowish.
    Rgb cool = BlackbodyToRgb(5000.0);
    EXPECT_GT(cool.r, cool.b) << "Cool outer disk should be redder than blue";

    // Very cool (3000K): distinctly red.
    Rgb cold = BlackbodyToRgb(3000.0);
    EXPECT_GT(cold.r, 0.8f) << "3000K should be red-dominant";

    using sirius::core::color_modes::ApplyColorMode;
    using sirius::core::color_modes::Mode;
    EXPECT_DEATH(static_cast<void>(ApplyColorMode(Mode::Polarisation, 1.0f, 1.0f, 1.0f)),
                 "violated");
    EXPECT_DEATH(static_cast<void>(ApplyColorMode(static_cast<Mode>(255), 1.0f, 1.0f, 1.0f)),
                 "violated");
}

// =========================================================================
// Doppler shift direction tests.
// =========================================================================

TEST(SpectralEmissionTest, DopplerShiftDirection) {
    using namespace sirius::core::spectral;

    double T_emit = 10000.0;

    double g_approach = 1.3;
    double T_blue = T_emit * g_approach;
    Rgb blueShifted = BlackbodyToRgb(T_blue);

    double g_recede = 0.7;
    double T_red = T_emit * g_recede;
    Rgb redShifted = BlackbodyToRgb(T_red);

    float blueRatio = blueShifted.b / std::max(blueShifted.r, 0.01f);
    float redRatio = redShifted.b / std::max(redShifted.r, 0.01f);
    EXPECT_GT(blueRatio, redRatio) << "Blue-shifted should have higher b/r ratio";
}

TEST(SpectralEmissionTest, TrueColorAppliesExactlyOneGFourthIntensityFactor) {
    using sirius::core::color_modes::ApplyColorMode;
    using sirius::core::color_modes::Mode;
    using sirius::core::spectral::BlackbodyToRgb;

    constexpr float emitted_temperature = 0.4f;
    constexpr float redshift = 1.25f;
    constexpr float emitted_intensity = 0.37f;
    constexpr float temperature_scale = 50000.0f;
    const float g4 = std::pow(redshift, 4.0f);
    const auto blackbody = BlackbodyToRgb(emitted_temperature * redshift * temperature_scale);
    const auto actual = ApplyColorMode(Mode::TrueColor, emitted_temperature, redshift,
                                       emitted_intensity, nullptr, temperature_scale);

    EXPECT_FLOAT_EQ(actual.r, blackbody.r * emitted_intensity * g4);
    EXPECT_FLOAT_EQ(actual.g, blackbody.g * emitted_intensity * g4);
    EXPECT_FLOAT_EQ(actual.b, blackbody.b * emitted_intensity * g4);
}

TEST(SpectralEmissionTest, MotionBlurAveragesNonlinearTemporalRadiance) {
    using sirius::core::color_modes::ApplyColorMode;
    using sirius::core::color_modes::AverageTemporalColorMode;
    using sirius::core::color_modes::Mode;

    constexpr float emitted_temperature = 0.4f;
    constexpr float emitted_intensity = 0.37f;
    constexpr float temperature_scale = 50000.0f;
    constexpr std::array<float, 2> redshifts = {0.5f, 1.5f};

    const auto actual = AverageTemporalColorMode(Mode::TrueColor, emitted_temperature, redshifts,
                                                 emitted_intensity, temperature_scale);
    const auto first = ApplyColorMode(Mode::TrueColor, emitted_temperature, redshifts[0],
                                      emitted_intensity, nullptr, temperature_scale);
    const auto second = ApplyColorMode(Mode::TrueColor, emitted_temperature, redshifts[1],
                                       emitted_intensity, nullptr, temperature_scale);
    EXPECT_FLOAT_EQ(actual.r, (first.r + second.r) * 0.5f);
    EXPECT_FLOAT_EQ(actual.g, (first.g + second.g) * 0.5f);
    EXPECT_FLOAT_EQ(actual.b, (first.b + second.b) * 0.5f);

    constexpr float averaged_redshift = (redshifts[0] + redshifts[1]) * 0.5f;
    const auto biased = ApplyColorMode(Mode::TrueColor, emitted_temperature, averaged_redshift,
                                       emitted_intensity, nullptr, temperature_scale);
    EXPECT_GT(std::abs(actual.r - biased.r) + std::abs(actual.g - biased.g) +
                  std::abs(actual.b - biased.b),
              0.01f)
        << "averaging g before the nonlinear radiance transform was not detected";
}

// =========================================================================
// Spin display format test.
// =========================================================================

TEST(SpectralEmissionTest, SpinDisplayFormat) {
    double spin = 0.998;
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(3) << spin;
    EXPECT_EQ(ss.str(), "0.998") << "Spin 0.998 must display as 0.998, not 1.00";

    std::ostringstream ss2;
    ss2 << std::fixed << std::setprecision(2) << spin;
    EXPECT_EQ(ss2.str(), "1.00") << "2-decimal precision rounds 0.998 to 1.00";
}
