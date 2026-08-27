// Spectral emission pipeline validation. Ported from TSIN007A.cpp.
//
// Validates the blackbody and relativistic-emission transforms used downstream
// of the disk model. The independent Page-Thorne flux oracle lives in the
// oracle suite rather than in this render-linked binary.

#include "sirius/core/spectral/blackbody.h"
#include "sirius/core/spectral/colour_modes.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>

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

TEST(SpectralEmissionTest, BolometricDiskAuthorityAppliesExactlyOneGFourthFactor) {
    using sirius::core::color_modes::ObservedBolometricIntensity;
    constexpr float emitted = 0.37f;
    constexpr float g = 1.25f;
    EXPECT_FLOAT_EQ(ObservedBolometricIntensity(emitted, g), emitted * g * g * g * g);
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
