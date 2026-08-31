// Non-render checks for the final live-view display boundary. The session
// supplies display-linear sRGB; the file-backed OpenGL shader may only apply
// the IEC transfer encode, and the optional film finish must preserve the
// finite unit-range input contract of that boundary.

#include "sirius/base/resource_locator.h"
#include "sirius/core/srgb_transfer.h"
#include "sirius/render/film_pipeline.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iterator>
#include <optional>
#include <regex>
#include <string>
#include <string_view>
#include <vector>

namespace sirius::app::test {
namespace {

std::optional<double> ShaderConstant(const std::string& source, std::string_view name) {
    const std::regex declaration(
        "const\\s+float\\s+" + std::string(name) +
        R"(\s*=\s*([0-9]+(?:\.[0-9]+)?)(?:\s*/\s*([0-9]+(?:\.[0-9]+)?))?\s*;)");
    std::smatch match;
    if (!std::regex_search(source, match, declaration)) return std::nullopt;
    const double numerator = std::stod(match[1].str());
    if (!match[2].matched) return numerator;
    return numerator / std::stod(match[2].str());
}

double MirrorEncode(double linear, double breakpoint, double slope, double scale, double offset,
                    double exponent) {
    linear = std::clamp(linear, 0.0, 1.0);
    if (linear <= breakpoint) return slope * linear;
    return scale * std::pow(linear, exponent) - offset;
}

TEST(ViewerDisplayContract, LiveShaderAppliesOnlyExactSrgbTransfer) {
    const auto path = base::ResolveResource("shaders/RDSD003A.frag");
    ASSERT_TRUE(path.has_value());
    std::ifstream input(*path, std::ios::binary);
    ASSERT_TRUE(input.good());
    const std::string source{std::istreambuf_iterator<char>(input),
                             std::istreambuf_iterator<char>()};

    const auto breakpoint = ShaderConstant(source, "kSrgbLinearBreakpoint");
    const auto slope = ShaderConstant(source, "kSrgbLinearSlope");
    const auto scale = ShaderConstant(source, "kSrgbPowerScale");
    const auto offset = ShaderConstant(source, "kSrgbPowerOffset");
    const auto exponent = ShaderConstant(source, "kSrgbPowerExponent");
    ASSERT_TRUE(breakpoint.has_value());
    ASSERT_TRUE(slope.has_value());
    ASSERT_TRUE(scale.has_value());
    ASSERT_TRUE(offset.has_value());
    ASSERT_TRUE(exponent.has_value());

    EXPECT_EQ(source.find("color / (color +"), std::string::npos);
    EXPECT_EQ(source.find("1.0 / 2.2"), std::string::npos);
    EXPECT_NE(source.find("if (linear <= kSrgbLinearBreakpoint)"), std::string::npos);
    EXPECT_NE(source.find("return kSrgbLinearSlope * linear;"), std::string::npos);
    EXPECT_NE(source.find("kSrgbPowerScale * pow(linear, kSrgbPowerExponent) - "
                          "kSrgbPowerOffset"),
              std::string::npos);

    for (int sample = 0; sample <= 4096; ++sample) {
        const double linear = -0.25 + 1.5 * static_cast<double>(sample) / 4096.0;
        const double mirror = MirrorEncode(linear, *breakpoint, *slope, *scale, *offset, *exponent);
        const double host = core::colour::EncodeClippedSrgbChannel(linear);
        EXPECT_NEAR(mirror, host, 2.0e-15) << linear;
    }
}

TEST(ViewerDisplayContract, FilmFinishClosesDisplayLinearRange) {
    render::FilmConfig config = render::FilmConfig::Interstellar();
    config.exposure = 2.0f;
    config.bloom_intensity = 5.0f;
    config.bloom_threshold = 0.0f;
    render::FilmPipeline pipeline(config);

    constexpr int kWidth = 4;
    constexpr int kHeight = 4;
    std::vector<float> pixels(static_cast<std::size_t>(kWidth * kHeight * 4));
    for (int pixel = 0; pixel < kWidth * kHeight; ++pixel) {
        const std::size_t index = static_cast<std::size_t>(pixel * 4);
        pixels[index + 0] = 4.0f;
        pixels[index + 1] = 2.0f;
        pixels[index + 2] = 0.5f;
        pixels[index + 3] = 0.25f;
    }

    pipeline.Apply(pixels.data(), kWidth, kHeight, 0);
    for (std::size_t index = 0; index < pixels.size(); index += 4) {
        for (std::size_t channel = 0; channel < 3; ++channel) {
            EXPECT_TRUE(std::isfinite(pixels[index + channel]));
            EXPECT_GE(pixels[index + channel], 0.0f);
            EXPECT_LE(pixels[index + channel], 1.0f);
        }
        EXPECT_FLOAT_EQ(pixels[index + 3], 0.25f);
    }
}

}  // namespace
}  // namespace sirius::app::test
