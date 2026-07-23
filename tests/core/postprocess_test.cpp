// Post-processing operator tests: tonemapping curves and dispatch, bloom,
// colour grading. Ported from TSPP001A.cpp; assertions and tolerances unchanged.

#include "sirius/core/postprocess.h"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

namespace sirius::test {
using namespace sirius::core;

constexpr float kEps = 1e-5f;

// =============================================================================
// Tonemapping Function Tests
// =============================================================================

class TonemapTests : public ::testing::Test {};

TEST_F(TonemapTests, ACESAtZero) { EXPECT_NEAR(sirius::core::tonemap::Aces(0.0f), 0.0f, kEps); }

TEST_F(TonemapTests, ACESMonotone) {
    float prev = sirius::core::tonemap::Aces(0.0f);
    for (float x = 0.1f; x <= 10.0f; x += 0.1f) {
        float val = sirius::core::tonemap::Aces(x);
        EXPECT_GE(val, prev - kEps) << "ACES not monotone at x = " << x;
        prev = val;
    }
}

TEST_F(TonemapTests, ACESBounded) {
    for (float x = 0.0f; x <= 100.0f; x += 0.5f) {
        float val = sirius::core::tonemap::Aces(x);
        EXPECT_GE(val, 0.0f);
        EXPECT_LE(val, 1.0f);
    }
}

TEST_F(TonemapTests, ReinhardAnalytic) {
    // Reinhard(x) = x / (1 + x)
    EXPECT_NEAR(sirius::core::tonemap::Reinhard(0.0f), 0.0f, kEps);
    EXPECT_NEAR(sirius::core::tonemap::Reinhard(1.0f), 0.5f, kEps);
    EXPECT_NEAR(sirius::core::tonemap::Reinhard(0.5f), 1.0f / 3.0f, kEps);
    EXPECT_NEAR(sirius::core::tonemap::Reinhard(10.0f), 10.0f / 11.0f, kEps);
    EXPECT_NEAR(sirius::core::tonemap::Reinhard(100.0f), 100.0f / 101.0f, kEps);
}

TEST_F(TonemapTests, FilmicAtZero) {
    EXPECT_NEAR(sirius::core::tonemap::Filmic(0.0f), 0.0f, 0.01f);
}

TEST_F(TonemapTests, FilmicMonotone) {
    float prev = sirius::core::tonemap::Filmic(0.0f);
    for (float x = 0.1f; x <= 10.0f; x += 0.1f) {
        float val = sirius::core::tonemap::Filmic(x);
        EXPECT_GE(val, prev - kEps) << "Filmic not monotone at x = " << x;
        prev = val;
    }
}

TEST_F(TonemapTests, FilmicBounded) {
    // Hable/Uncharted 2 Filmic overshoots 1.0 (asymptotes to ~1.29)
    // because curve(inf)/curve(whitePoint) > 1. Verify non-negative
    // and bounded by the analytic asymptotic limit.
    float asymptotic = sirius::core::tonemap::Filmic(1e6f);
    for (float x = 0.0f; x <= 100.0f; x += 1.0f) {
        float val = sirius::core::tonemap::Filmic(x);
        EXPECT_GE(val, -kEps);
        EXPECT_LE(val, asymptotic + kEps);
    }
}

// =============================================================================
// Tonemap::apply Dispatch Tests
// =============================================================================

TEST_F(TonemapTests, ApplyExposureScaling) {
    float r = 0.5f, g = 0.5f, b = 0.5f;
    float exposure = 2.0f;
    // Using None: just clamp after exposure scaling
    sirius::core::tonemap::Apply(r, g, b, sirius::core::TonemapType::None, exposure);
    EXPECT_NEAR(r, 1.0f, kEps);  // 0.5 * 2.0 = 1.0
}

TEST_F(TonemapTests, ApplyACESDispatch) {
    float r = 1.0f, g = 1.0f, b = 1.0f;
    sirius::core::tonemap::Apply(r, g, b, sirius::core::TonemapType::Aces, 1.0f);
    // All channels should be tonemapped identically for equal input
    EXPECT_NEAR(r, g, kEps);
    EXPECT_NEAR(g, b, kEps);
    EXPECT_GT(r, 0.0f);
    EXPECT_LT(r, 1.0f);
}

TEST_F(TonemapTests, ApplyReinhardDispatch) {
    float r = 1.0f, g = 0.0f, b = 0.0f;
    sirius::core::tonemap::Apply(r, g, b, sirius::core::TonemapType::Reinhard, 1.0f);
    EXPECT_NEAR(r, 0.5f, kEps);
    EXPECT_NEAR(g, 0.0f, kEps);
}

// =============================================================================
// BloomFilter Tests
// =============================================================================

class BloomFilterTests : public ::testing::Test {};

TEST_F(BloomFilterTests, DisabledBloomLeavesBufferUnchanged) {
    int w = 4, h = 4;
    std::vector<float> buffer(w * h * 4, 0.5f);
    std::vector<float> original = buffer;

    sirius::core::PostProcessConfig config;
    config.enable_bloom = false;

    sirius::core::BloomFilter::Apply(buffer, w, h, config);
    EXPECT_EQ(buffer, original);
}

TEST_F(BloomFilterTests, ZeroIntensityBloomLeavesUnchanged) {
    int w = 4, h = 4;
    std::vector<float> buffer(w * h * 4, 0.5f);
    std::vector<float> original = buffer;

    sirius::core::PostProcessConfig config;
    config.enable_bloom = true;
    config.bloom_intensity = 0.0f;

    sirius::core::BloomFilter::Apply(buffer, w, h, config);
    EXPECT_EQ(buffer, original);
}

TEST_F(BloomFilterTests, BrightPixelsCauseBloom) {
    int w = 8, h = 8;
    std::vector<float> buffer(w * h * 4, 0.0f);
    // Place a bright pixel at centre
    int cx = 4, cy = 4;
    int idx = (cy * w + cx) * 4;
    buffer[idx] = 5.0f;
    buffer[idx + 1] = 5.0f;
    buffer[idx + 2] = 5.0f;
    buffer[idx + 3] = 1.0f;

    sirius::core::PostProcessConfig config;
    config.enable_bloom = true;
    config.bloom_intensity = 1.0f;
    config.bloom_threshold = 0.5f;
    config.bloom_radius = 2;

    sirius::core::BloomFilter::Apply(buffer, w, h, config);

    // Neighbours should now have non-zero values from bloom spread
    int neighbour = ((cy)*w + (cx + 1)) * 4;
    EXPECT_GT(buffer[neighbour], 0.0f) << "Bloom should spread to neighbours";
}

// =============================================================================
// ColourGrading Tests
// =============================================================================

class ColourGradingTests : public ::testing::Test {};

TEST_F(ColourGradingTests, IdentityParametersPreserveValues) {
    sirius::core::PostProcessConfig config;
    config.saturation = 1.0f;
    config.contrast = 1.0f;
    config.lift = 0.0f;
    config.gain = 1.0f;

    float r = 0.5f, g = 0.3f, b = 0.7f;
    float r0 = r, g0 = g, b0 = b;
    sirius::core::ColourGrading::Apply(r, g, b, config);
    EXPECT_NEAR(r, r0, kEps);
    EXPECT_NEAR(g, g0, kEps);
    EXPECT_NEAR(b, b0, kEps);
}

TEST_F(ColourGradingTests, ZeroSaturationProducesGrey) {
    sirius::core::PostProcessConfig config;
    config.saturation = 0.0f;
    config.contrast = 1.0f;
    config.lift = 0.0f;
    config.gain = 1.0f;

    float r = 0.8f, g = 0.2f, b = 0.5f;
    sirius::core::ColourGrading::Apply(r, g, b, config);
    // All channels should be the luminance value
    EXPECT_NEAR(r, g, kEps);
    EXPECT_NEAR(g, b, kEps);
}

TEST_F(ColourGradingTests, OutputClamped) {
    sirius::core::PostProcessConfig config;
    config.gain = 10.0f;  // extreme gain
    config.lift = 0.0f;
    config.saturation = 1.0f;
    config.contrast = 1.0f;

    float r = 0.8f, g = 0.8f, b = 0.8f;
    sirius::core::ColourGrading::Apply(r, g, b, config);
    EXPECT_LE(r, 1.0f);
    EXPECT_GE(r, 0.0f);
}

}  // namespace sirius::test
