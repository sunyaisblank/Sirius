// TSFL001A.cpp - Film Simulation Unit Tests
// Component ID: TSFL001A (Test/Film/Invariants)
//
// Tests for IMAX 70mm film simulation including grain, halation,
// and color grading.

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "../../Sirius.Render/Output/RDFL001A.h"

using namespace Sirius;

class FilmSimulationTest : public ::testing::Test {
protected:
    FilmConfig config;
    std::unique_ptr<FilmPipeline> pipeline;

    void SetUp() override {
        config = FilmConfig::Interstellar();
        pipeline = std::make_unique<FilmPipeline>(config);
    }

    // Helper to create test framebuffer
    std::vector<float> createTestBuffer(int width, int height, float r, float g, float b, float a = 1.0f) {
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

//==============================================================================
// Film Format Tests
//==============================================================================

TEST_F(FilmSimulationTest, IMAX_AspectRatio_143) {
    FilmConfig cfg;
    cfg.applyFormat(FilmFormat::IMAX70mm_15perf);
    EXPECT_NEAR(cfg.aspect_ratio, 1.43f, 0.01f);
}

TEST_F(FilmSimulationTest, IMAX5perf_AspectRatio_220) {
    FilmConfig cfg;
    cfg.applyFormat(FilmFormat::IMAX70mm_5perf);
    EXPECT_NEAR(cfg.aspect_ratio, 2.20f, 0.01f);
}

TEST_F(FilmSimulationTest, ComputeHeight_Correct) {
    FilmConfig cfg;
    cfg.width = 4096;
    cfg.aspect_ratio = 1.43f;
    cfg.computeHeight();

    // Height = width / aspect_ratio
    int expected_height = static_cast<int>(std::round(4096.0f / 1.43f));
    expected_height = (expected_height / 2) * 2;  // Even
    EXPECT_EQ(cfg.height, expected_height);
}

TEST_F(FilmSimulationTest, ComputeHeight_EvenDimension) {
    FilmConfig cfg;
    cfg.width = 1920;
    cfg.aspect_ratio = 1.43f;
    cfg.computeHeight();

    // Height should be even for video encoding
    EXPECT_EQ(cfg.height % 2, 0u);
}

//==============================================================================
// Film Stock Tests
//==============================================================================

TEST_F(FilmSimulationTest, KodakVision3_500T_Settings) {
    FilmConfig cfg;
    cfg.applyStock(FilmStock::KodakVision3_500T);

    EXPECT_FLOAT_EQ(cfg.iso, 500.0f);
    EXPECT_GT(cfg.grain_intensity, 0.0f);
    EXPECT_NEAR(cfg.color_temperature_K, 3200.0f, 100.0f);  // Tungsten
}

TEST_F(FilmSimulationTest, KodakVision3_50D_LowerGrain) {
    FilmConfig cfg_500T, cfg_50D;
    cfg_500T.applyStock(FilmStock::KodakVision3_500T);
    cfg_50D.applyStock(FilmStock::KodakVision3_50D);

    // Lower ISO = finer grain
    EXPECT_LT(cfg_50D.grain_intensity, cfg_500T.grain_intensity);
}

//==============================================================================
// Grain Tests
//==============================================================================

TEST_F(FilmSimulationTest, Grain_AddsNoise) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 0.5f, 0.5f, 0.5f);
    std::vector<float> original = buffer;

    config.grain_intensity = 0.1f;
    config.grain_enabled = true;
    pipeline->setConfig(config);
    pipeline->applyGrain(buffer.data(), width, height, 12345);

    // Buffer should have changed
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
    pipeline->setConfig(config);

    pipeline->applyGrain(buffer1.data(), width, height, 1);
    pipeline->applyGrain(buffer2.data(), width, height, 2);

    // Different seeds should produce different noise
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
    auto buffer = createTestBuffer(width, height, 0.1f, 0.1f, 0.1f);  // Dark image

    config.grain_intensity = 0.5f;  // Strong grain
    pipeline->setConfig(config);
    pipeline->applyGrain(buffer.data(), width, height, 42);

    // Result should be non-negative
    for (int i = 0; i < width * height * 4; i += 4) {
        EXPECT_GE(buffer[i + 0], 0.0f);
        EXPECT_GE(buffer[i + 1], 0.0f);
        EXPECT_GE(buffer[i + 2], 0.0f);
    }
}

//==============================================================================
// Halation Tests
//==============================================================================

TEST_F(FilmSimulationTest, Halation_AffectsBrightAreas) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 0.0f, 0.0f, 0.0f);

    // Create bright spot in center
    int cx = width / 2, cy = height / 2;
    for (int dy = -2; dy <= 2; ++dy) {
        for (int dx = -2; dx <= 2; ++dx) {
            int idx = ((cy + dy) * width + (cx + dx)) * 4;
            buffer[idx + 0] = 2.0f;  // Over threshold
            buffer[idx + 1] = 2.0f;
            buffer[idx + 2] = 2.0f;
        }
    }

    std::vector<float> original = buffer;

    config.halation_enabled = true;
    config.halation_threshold = 0.8f;
    config.halation_strength = 0.5f;
    config.halation_radius = 5.0f;
    pipeline->setConfig(config);
    pipeline->applyHalation(buffer.data(), width, height);

    // Areas around bright spot should have increased
    int sample_idx = ((cy + 5) * width + cx) * 4;
    EXPECT_GT(buffer[sample_idx], original[sample_idx]);
}

TEST_F(FilmSimulationTest, Halation_RedBias) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 0.0f, 0.0f, 0.0f);

    // Bright center
    buffer[(32 * width + 32) * 4 + 0] = 2.0f;
    buffer[(32 * width + 32) * 4 + 1] = 2.0f;
    buffer[(32 * width + 32) * 4 + 2] = 2.0f;

    config.halation_color_r = 1.0f;
    config.halation_color_g = 0.5f;
    config.halation_color_b = 0.2f;
    pipeline->setConfig(config);
    pipeline->applyHalation(buffer.data(), width, height);

    // Sample nearby pixel - should have more red than blue
    int sample_idx = (32 * width + 36) * 4;
    if (buffer[sample_idx] > 0.0f) {
        EXPECT_GE(buffer[sample_idx + 0], buffer[sample_idx + 2]);  // R >= B
    }
}

//==============================================================================
// Vignette Tests
//==============================================================================

TEST_F(FilmSimulationTest, Vignette_DarkensCorners) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 1.0f, 1.0f, 1.0f);

    config.vignette_enabled = true;
    config.vignette_strength = 0.5f;
    pipeline->setConfig(config);
    pipeline->applyVignette(buffer.data(), width, height);

    // Center should be brighter than corners
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
    pipeline->setConfig(config);
    pipeline->applyVignette(buffer.data(), width, height);

    // Center should be relatively unchanged
    int cx = width / 2, cy = height / 2;
    float center_r = buffer[(cy * width + cx) * 4];
    EXPECT_NEAR(center_r, 0.5f, 0.1f);
}

//==============================================================================
// Color Grading Tests
//==============================================================================

TEST_F(FilmSimulationTest, Exposure_Positive_Brightens) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 0.2f, 0.2f, 0.2f);

    config.exposure = 1.0f;  // +1 stop
    pipeline->setConfig(config);
    pipeline->applyColorGrade(buffer.data(), width, height);

    // Should be brighter (approximately 2x for +1 stop before tone mapping)
    float result = buffer[0];
    EXPECT_GT(result, 0.2f);
}

TEST_F(FilmSimulationTest, Exposure_Negative_Darkens) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 0.5f, 0.5f, 0.5f);

    config.exposure = -1.0f;  // -1 stop
    pipeline->setConfig(config);
    pipeline->applyColorGrade(buffer.data(), width, height);

    // Should be darker
    float result = buffer[0];
    EXPECT_LT(result, 0.5f);
}

TEST_F(FilmSimulationTest, Saturation_Zero_Grayscale) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 1.0f, 0.5f, 0.0f);  // Orange

    config.saturation = 0.0f;
    config.exposure = 0.0f;
    pipeline->setConfig(config);
    pipeline->applyColorGrade(buffer.data(), width, height);

    // After saturation=0, RGB should be equal (grayscale)
    float r = buffer[0], g = buffer[1], b = buffer[2];
    EXPECT_NEAR(r, g, 0.1f);
    EXPECT_NEAR(g, b, 0.1f);
}

TEST_F(FilmSimulationTest, OutputClamped_0_1) {
    int width = 64, height = 64;
    auto buffer = createTestBuffer(width, height, 10.0f, 10.0f, 10.0f);  // Overexposed

    pipeline->applyColorGrade(buffer.data(), width, height);

    // Output should be clamped to [0, 1]
    for (int i = 0; i < width * height * 4; i += 4) {
        EXPECT_LE(buffer[i + 0], 1.0f);
        EXPECT_LE(buffer[i + 1], 1.0f);
        EXPECT_LE(buffer[i + 2], 1.0f);
        EXPECT_GE(buffer[i + 0], 0.0f);
        EXPECT_GE(buffer[i + 1], 0.0f);
        EXPECT_GE(buffer[i + 2], 0.0f);
    }
}

//==============================================================================
// Preset Tests
//==============================================================================

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

//==============================================================================
// GPU Parameter Conversion Tests
//==============================================================================

// Feature flags for FilmParamsGPU
static constexpr uint32_t FILM_FEATURE_GRAIN    = 1 << 0;
static constexpr uint32_t FILM_FEATURE_HALATION = 1 << 1;
static constexpr uint32_t FILM_FEATURE_VIGNETTE = 1 << 2;
static constexpr uint32_t FILM_FEATURE_ENABLED  = 1 << 3;

// Helper to create GPU params from CPU config
static FilmParamsGPU createFilmParamsGPU(const FilmConfig& cpu, uint32_t frame_seed = 0) {
    FilmParamsGPU gpu;
    gpu.grain_intensity = cpu.grain_intensity;
    gpu.grain_size = cpu.grain_size;
    gpu.grain_uniformity = cpu.grain_uniformity;
    gpu.grain_seed = frame_seed;
    gpu.halation_radius = cpu.halation_radius;
    gpu.halation_strength = cpu.halation_strength;
    gpu.halation_threshold = cpu.halation_threshold;
    gpu.halation_color_r = cpu.halation_color_r;
    gpu.halation_color_g = cpu.halation_color_g;
    gpu.halation_color_b = cpu.halation_color_b;
    gpu.padding1 = gpu.padding2 = 0.0f;
    gpu.saturation = cpu.saturation;
    gpu.contrast = cpu.contrast;
    gpu.exposure = cpu.exposure;
    gpu.toe_strength = cpu.toe_strength;
    gpu.shoulder_strength = cpu.shoulder_strength;
    gpu.midtone_point = cpu.midtone_point;
    gpu.padding3 = gpu.padding4 = 0.0f;
    gpu.vignette_strength = cpu.vignette_strength;
    gpu.vignette_radius = cpu.vignette_radius;
    gpu.vignette_softness = cpu.vignette_softness;
    gpu.padding5 = 0.0f;
    gpu.features = 0;
    if (cpu.grain_enabled) gpu.features |= FILM_FEATURE_GRAIN;
    if (cpu.halation_enabled) gpu.features |= FILM_FEATURE_HALATION;
    if (cpu.vignette_enabled) gpu.features |= FILM_FEATURE_VIGNETTE;
    if (cpu.enabled) gpu.features |= FILM_FEATURE_ENABLED;
    return gpu;
}

TEST_F(FilmSimulationTest, GPUConversion_PreservesValues) {
    FilmParamsGPU gpu = createFilmParamsGPU(config, 42);

    EXPECT_FLOAT_EQ(gpu.grain_intensity, config.grain_intensity);
    EXPECT_FLOAT_EQ(gpu.halation_radius, config.halation_radius);
    EXPECT_FLOAT_EQ(gpu.saturation, config.saturation);
    EXPECT_FLOAT_EQ(gpu.vignette_strength, config.vignette_strength);
    EXPECT_EQ(gpu.grain_seed, 42u);
}

TEST_F(FilmSimulationTest, GPUConversion_FeatureFlags) {
    config.grain_enabled = true;
    config.halation_enabled = false;
    config.vignette_enabled = true;
    config.enabled = true;

    FilmParamsGPU gpu = createFilmParamsGPU(config, 0);

    EXPECT_TRUE(gpu.features & FILM_FEATURE_GRAIN);
    EXPECT_FALSE(gpu.features & FILM_FEATURE_HALATION);
    EXPECT_TRUE(gpu.features & FILM_FEATURE_VIGNETTE);
    EXPECT_TRUE(gpu.features & FILM_FEATURE_ENABLED);
}

//==============================================================================
// Full Pipeline Test
//==============================================================================

TEST_F(FilmSimulationTest, FullPipeline_DoesNotCrash) {
    int width = 256, height = 179;  // IMAX-ish aspect
    auto buffer = createTestBuffer(width, height, 0.5f, 0.3f, 0.1f);

    // Enable everything
    config.grain_enabled = true;
    config.halation_enabled = true;
    config.vignette_enabled = true;
    config.bloom_enabled = true;
    config.enabled = true;
    pipeline->setConfig(config);

    // Should not crash
    EXPECT_NO_THROW(pipeline->apply(buffer.data(), width, height, 0));

    // Output should be valid
    for (size_t i = 0; i < buffer.size(); ++i) {
        EXPECT_FALSE(std::isnan(buffer[i]));
        EXPECT_FALSE(std::isinf(buffer[i]));
    }
}
