#include "sirius/core/constants.h"
#include "sirius/core/spectral/blackbody.h"
#include "sirius/core/srgb_transfer.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

namespace sirius::core::test {
namespace {

double IndependentIecSrgbEncode(double linear) {
    if (linear <= 0.0031308) return 12.92 * linear;
    return 1.055 * std::pow(linear, 1.0 / 2.4) - 0.055;
}

TEST(SrgbTransferAuthority, CurveMatchesIndependentIecOracleAcrossBothBranches) {
    constexpr double kBreakpoint = 0.0031308;
    const std::array inputs{
        -0.25,
        0.0,
        std::nextafter(kBreakpoint, 0.0),
        kBreakpoint,
        std::nextafter(kBreakpoint, 1.0),
        0.018,
        0.18,
        0.5,
        1.0,
        1.5,
    };

    for (const double input : inputs) {
        const double expected = IndependentIecSrgbEncode(input);
        const double tolerance = 16.0 * constants::kEpsilonD * std::max(1.0, std::abs(expected));
        EXPECT_NEAR(colour::EncodeSrgbChannel(input), expected, tolerance) << input;
    }
}

TEST(SrgbTransferAuthority, EightBitQuantisationClipsAndDeclinesNonfiniteInputs) {
    struct Case {
        double linear;
        std::uint8_t expected;
    };
    constexpr std::array cases{
        Case{-1.0, 0},  Case{0.0, 0},   Case{0.0031308, 10}, Case{0.18, 118},
        Case{0.5, 188}, Case{1.0, 255}, Case{2.0, 255},
    };
    for (const Case& sample : cases) {
        const auto encoded = colour::TryEncodeSrgb8(sample.linear);
        ASSERT_TRUE(encoded.has_value()) << sample.linear;
        EXPECT_EQ(*encoded, sample.expected) << sample.linear;
    }

    EXPECT_FALSE(colour::TryEncodeSrgb8(std::numeric_limits<double>::quiet_NaN()).has_value());
    EXPECT_FALSE(colour::TryEncodeSrgb8(std::numeric_limits<double>::infinity()).has_value());
    EXPECT_FALSE(colour::TryEncodeSrgb8(-std::numeric_limits<double>::infinity()).has_value());
}

TEST(SrgbTransferAuthority, SpectralFacadeDelegatesWithoutChangingClippingSemantics) {
    constexpr std::array inputs{-0.25f, 0.0f, 0.0031308f, 0.18f, 0.5f, 1.0f, 1.5f};
    for (const float input : inputs) {
        EXPECT_FLOAT_EQ(spectral::SrgbGamma(input), colour::EncodeSrgbChannel(input));
    }

    const spectral::Rgb input{-0.25f, 0.18f, 1.5f};
    const spectral::Rgb encoded = spectral::LinearToSrgb(input);
    EXPECT_FLOAT_EQ(encoded.r, colour::EncodeClippedSrgbChannel(input.r));
    EXPECT_FLOAT_EQ(encoded.g, colour::EncodeClippedSrgbChannel(input.g));
    EXPECT_FLOAT_EQ(encoded.b, colour::EncodeClippedSrgbChannel(input.b));
    EXPECT_TRUE(std::isnan(spectral::SrgbGamma(std::numeric_limits<float>::quiet_NaN())));
}

}  // namespace
}  // namespace sirius::core::test
