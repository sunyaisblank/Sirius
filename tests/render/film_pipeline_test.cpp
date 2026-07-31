// Film simulation unit tests. Ported from TSFL001A.cpp.
//
// Tests IMAX 70mm film simulation (format, stock, grain, halation, vignette,
// colour grade, presets, full pipeline). Film is intentionally host-side for
// both render backends, after linear readback and before display encoding.

#include "sirius/render/film_pipeline.h"

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

using namespace sirius::render;

namespace {

class FilmSimulationTest : public ::testing::Test {
  protected:
    FilmConfig config;
    std::unique_ptr<FilmPipeline> pipeline;

    void SetUp() override {
        config = FilmConfig::Interstellar();
        pipeline = std::make_unique<FilmPipeline>(config);
    }

    std::vector<float> createTestBuffer(int width, int height, float r, float g, float b,
                                        float a = 1.0f) {
        std::vector<float> buffer(width * height * 4);
        for (int i = 0; i < width * height; ++i) {
            buffer[i * 4 + 0] = r;
            buffer[i * 4 + 1] = g;
            buffer[i * 4 + 2] = b;
            buffer[i * 4 + 3] = a;
        }
        return buffer;
    }
};

TEST_F(FilmSimulationTest, IMAX_AspectRatio_143) {
    FilmConfig cfg;
    cfg.ApplyFormat(FilmFormat::IMAX70mm_15perf);
    EXPECT_NEAR(cfg.aspect_ratio, 1.43f, 0.01f);
    EXPECT_DEATH(cfg.ApplyFormat(static_cast<FilmFormat>(255)), "violated");
}

TEST_F(FilmSimulationTest, IMAX5perf_AspectRatio_220) {
    FilmConfig cfg;
    cfg.ApplyFormat(FilmFormat::IMAX70mm_5perf);
    EXPECT_NEAR(cfg.aspect_ratio, 2.20f, 0.01f);
}

TEST_F(FilmSimulationTest, ComputeHeight_Correct) {
    FilmConfig cfg;
    cfg.width = 4096;
    cfg.aspect_ratio = 1.43f;
    cfg.ComputeHeight();

    int expected_height = static_cast<int>(std::round(4096.0f / 1.43f));
    expected_height = (expected_height / 2) * 2;  // Even.
    EXPECT_EQ(cfg.height, static_cast<uint32_t>(expected_height));
}

TEST_F(FilmSimulationTest, ComputeHeight_EvenDimension) {
    FilmConfig cfg;
    cfg.width = 1920;
    cfg.aspect_ratio = 1.43f;
    cfg.ComputeHeight();

    EXPECT_EQ(cfg.height % 2, 0u);
}

TEST_F(FilmSimulationTest, KodakVision3_500T_Settings) {
    FilmConfig cfg;
    cfg.ApplyStock(FilmStock::KodakVision3_500T);

    EXPECT_FLOAT_EQ(cfg.iso, 500.0f);
    EXPECT_GT(cfg.grain_intensity, 0.0f);
    EXPECT_NEAR(cfg.color_temperature_K, 3200.0f, 100.0f);  // Tungsten.
    EXPECT_DEATH(cfg.ApplyStock(static_cast<FilmStock>(255)), "violated");
}

TEST_F(FilmSimulationTest, KodakVision3_50D_LowerGrain) {
    FilmConfig cfg_500T, cfg_50D;
    cfg_500T.ApplyStock(FilmStock::KodakVision3_500T);
    cfg_50D.ApplyStock(FilmStock::KodakVision3_50D);

    EXPECT_LT(cfg_50D.grain_intensity, cfg_500T.grain_intensity);
}

TEST_F(FilmSimulationTest, Grain_AddsNoise) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 0.5f, 0.5f, 0.5f);
    std::vector<float> original = buffer;

    config.grain_intensity = 0.1f;
    config.grain_enabled = true;
    pipeline->SetConfig(config);
    pipeline->ApplyGrain(buffer.data(), width, height, 12345);

    bool changed = false;
    for (size_t i = 0; i < buffer.size(); ++i) {
        if (std::abs(buffer[i] - original[i]) > 1e-6f) {
            changed = true;
            break;
        }
    }
    EXPECT_TRUE(changed);
}

TEST_F(FilmSimulationTest, Grain_DifferentFrames_DifferentNoise) {
    int width = 64, height = 64;
    auto buffer1 = createTestBuffer(width, height, 0.5f, 0.5f, 0.5f);
    auto buffer2 = createTestBuffer(width, height, 0.5f, 0.5f, 0.5f);

    config.grain_intensity = 0.1f;
    pipeline->SetConfig(config);

    pipeline->ApplyGrain(buffer1.data(), width, height, 1);
    pipeline->ApplyGrain(buffer2.data(), width, height, 2);

    bool different = false;
    for (size_t i = 0; i < buffer1.size(); ++i) {
        if (std::abs(buffer1[i] - buffer2[i]) > 1e-6f) {
            different = true;
            break;
        }
    }
    EXPECT_TRUE(different);
}

TEST_F(FilmSimulationTest, Grain_NonNegativeResult) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 0.1f, 0.1f, 0.1f);  // Dark image.

    config.grain_intensity = 0.5f;  // Strong grain.
    pipeline->SetConfig(config);
    pipeline->ApplyGrain(buffer.data(), width, height, 42);

    for (int i = 0; i < width * height * 4; i += 4) {
        EXPECT_GE(buffer[i + 0], 0.0f);
        EXPECT_GE(buffer[i + 1], 0.0f);
        EXPECT_GE(buffer[i + 2], 0.0f);
    }
}

TEST_F(FilmSimulationTest, Halation_AffectsBrightAreas) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 0.0f, 0.0f, 0.0f);

    int cx = width / 2, cy = height / 2;
    for (int dy = -2; dy <= 2; ++dy) {
        for (int dx = -2; dx <= 2; ++dx) {
            int idx = ((cy + dy) * width + (cx + dx)) * 4;
            buffer[idx + 0] = 2.0f;  // Over threshold.
            buffer[idx + 1] = 2.0f;
            buffer[idx + 2] = 2.0f;
        }
    }

    std::vector<float> original = buffer;

    config.halation_enabled = true;
    config.halation_threshold = 0.8f;
    config.halation_strength = 0.5f;
    config.halation_radius = 5.0f;
    pipeline->SetConfig(config);
    pipeline->ApplyHalation(buffer.data(), width, height);

    int sample_idx = ((cy + 5) * width + cx) * 4;
    EXPECT_GT(buffer[sample_idx], original[sample_idx]);
}

TEST_F(FilmSimulationTest, Halation_RedBias) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 0.0f, 0.0f, 0.0f);

    buffer[(32 * width + 32) * 4 + 0] = 2.0f;
    buffer[(32 * width + 32) * 4 + 1] = 2.0f;
    buffer[(32 * width + 32) * 4 + 2] = 2.0f;

    config.halation_color_r = 1.0f;
    config.halation_color_g = 0.5f;
    config.halation_color_b = 0.2f;
    pipeline->SetConfig(config);
    pipeline->ApplyHalation(buffer.data(), width, height);

    int sample_idx = (32 * width + 36) * 4;
    if (buffer[sample_idx] > 0.0f) {
        EXPECT_GE(buffer[sample_idx + 0], buffer[sample_idx + 2]);  // R >= B.
    }
}

TEST_F(FilmSimulationTest, Vignette_DarkensCorners) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 1.0f, 1.0f, 1.0f);

    config.vignette_enabled = true;
    config.vignette_strength = 0.5f;
    pipeline->SetConfig(config);
    pipeline->ApplyVignette(buffer.data(), width, height);

    float center = buffer[(32 * width + 32) * 4];
    float corner = buffer[0];

    EXPECT_GT(center, corner);
}

TEST_F(FilmSimulationTest, Vignette_CenterUnchanged) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 0.5f, 0.5f, 0.5f);

    config.vignette_enabled = true;
    config.vignette_strength = 0.3f;
    config.vignette_radius = 1.2f;
    pipeline->SetConfig(config);
    pipeline->ApplyVignette(buffer.data(), width, height);

    int cx = width / 2, cy = height / 2;
    float center_r = buffer[(cy * width + cx) * 4];
    EXPECT_NEAR(center_r, 0.5f, 0.1f);
}

TEST_F(FilmSimulationTest, Exposure_Positive_Brightens) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 0.2f, 0.2f, 0.2f);

    config.exposure = 1.0f;  // +1 stop.
    pipeline->SetConfig(config);
    pipeline->ApplyColorGrade(buffer.data(), width, height);

    float result = buffer[0];
    EXPECT_GT(result, 0.2f);
}

TEST_F(FilmSimulationTest, Exposure_Negative_Darkens) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 0.5f, 0.5f, 0.5f);

    config.exposure = -1.0f;  // -1 stop.
    pipeline->SetConfig(config);
    pipeline->ApplyColorGrade(buffer.data(), width, height);

    float result = buffer[0];
    EXPECT_LT(result, 0.5f);
}

TEST_F(FilmSimulationTest, Saturation_Zero_Grayscale) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 1.0f, 0.5f, 0.0f);  // Orange.

    config.saturation = 0.0f;
    config.exposure = 0.0f;
    pipeline->SetConfig(config);
    pipeline->ApplyColorGrade(buffer.data(), width, height);

    float r = buffer[0], g = buffer[1], b = buffer[2];
    EXPECT_NEAR(r, g, 0.1f);
    EXPECT_NEAR(g, b, 0.1f);
}

TEST_F(FilmSimulationTest, OutputClamped_0_1) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 10.0f, 10.0f, 10.0f);  // Overexposed.

    pipeline->ApplyColorGrade(buffer.data(), width, height);

    for (int i = 0; i < width * height * 4; i += 4) {
        EXPECT_LE(buffer[i + 0], 1.0f);
        EXPECT_LE(buffer[i + 1], 1.0f);
        EXPECT_LE(buffer[i + 2], 1.0f);
        EXPECT_GE(buffer[i + 0], 0.0f);
        EXPECT_GE(buffer[i + 1], 0.0f);
        EXPECT_GE(buffer[i + 2], 0.0f);
    }
}

TEST_F(FilmSimulationTest, Interstellar_Preset) {
    FilmConfig cfg = FilmConfig::Interstellar();

    EXPECT_EQ(cfg.format, FilmFormat::IMAX70mm_15perf);
    EXPECT_NEAR(cfg.aspect_ratio, 1.43f, 0.01f);
    EXPECT_TRUE(cfg.bloom_enabled);
    EXPECT_TRUE(cfg.vignette_enabled);
}

TEST_F(FilmSimulationTest, DigitalClean_NoEffects) {
    FilmConfig cfg = FilmConfig::DigitalClean();

    EXPECT_FALSE(cfg.enabled);
    EXPECT_FALSE(cfg.grain_enabled);
    EXPECT_FALSE(cfg.halation_enabled);
    EXPECT_FALSE(cfg.vignette_enabled);
}

TEST_F(FilmSimulationTest, FullPipeline_DoesNotCrash) {
    int width = 256, height = 179;  // IMAX-ish aspect.
    auto buffer = createTestBuffer(width, height, 0.5f, 0.3f, 0.1f);

    config.grain_enabled = true;
    config.halation_enabled = true;
    config.vignette_enabled = true;
    config.bloom_enabled = true;
    config.enabled = true;
    pipeline->SetConfig(config);

    EXPECT_NO_THROW(pipeline->Apply(buffer.data(), width, height, 0));

    for (size_t i = 0; i < buffer.size(); ++i) {
        EXPECT_FALSE(std::isnan(buffer[i]));
        EXPECT_FALSE(std::isinf(buffer[i]));
    }
}

}  // namespace
