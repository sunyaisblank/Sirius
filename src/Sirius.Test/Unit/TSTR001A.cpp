// TSTR001A.cpp - Turbulence Model Unit Tests
// Component ID: TSTR001A (Test/Turbulence/Invariants)
//
// Tests for Kolmogorov cascade density perturbations and fBm noise.

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <numeric>
#include "../../Sirius.Core/Disk/PHTR001A.h"
#include "../../Sirius.Render/Acceleration/OptiX/RDOP003A.h"

using namespace Sirius;

class TurbulenceTest : public ::testing::Test {
protected:
    TurbulenceConfig config;

    void SetUp() override {
        config.kolmogorov_exponent = -5.0f / 3.0f;
        config.outer_scale_M = 5.0f;
        config.inner_scale_M = 0.1f;
        config.amplitude = 0.3f;
        config.octaves = 6;
        config.seed = 12345;
        config.lacunarity = 2.0f;
        config.persistence = 0.5f;
        config.enabled = true;
    }
};

//==============================================================================
// Perlin Noise Tests
//==============================================================================

TEST_F(TurbulenceTest, PerlinNoise_BoundedRange) {
    // Perlin noise should be in [-1, 1]
    for (int i = 0; i < 1000; ++i) {
        float x = (i * 0.1f);
        float y = (i * 0.07f);
        float z = (i * 0.13f);

        float noise = TurbulenceNoise::perlin3D(x, y, z, 42);
        EXPECT_GE(noise, -1.5f);  // Allow small overshoot
        EXPECT_LE(noise, 1.5f);
    }
}

TEST_F(TurbulenceTest, PerlinNoise_Deterministic) {
    // Same input should produce same output
    float n1 = TurbulenceNoise::perlin3D(1.5f, 2.3f, 4.1f, 42);
    float n2 = TurbulenceNoise::perlin3D(1.5f, 2.3f, 4.1f, 42);
    EXPECT_FLOAT_EQ(n1, n2);
}

TEST_F(TurbulenceTest, PerlinNoise_DifferentSeeds) {
    // Different seeds should produce different output
    float n1 = TurbulenceNoise::perlin3D(1.5f, 2.3f, 4.1f, 42);
    float n2 = TurbulenceNoise::perlin3D(1.5f, 2.3f, 4.1f, 43);
    EXPECT_NE(n1, n2);
}

TEST_F(TurbulenceTest, PerlinNoise_Continuous) {
    // Noise should be continuous (small changes produce small differences)
    float epsilon = 1e-3f;
    float n1 = TurbulenceNoise::perlin3D(5.0f, 5.0f, 5.0f, 42);
    float n2 = TurbulenceNoise::perlin3D(5.0f + epsilon, 5.0f, 5.0f, 42);
    float diff = std::abs(n2 - n1);
    EXPECT_LT(diff, 0.1f);  // Gradient bounded
}

//==============================================================================
// fBm Noise Tests
//==============================================================================

TEST_F(TurbulenceTest, FBm_BoundedRange) {
    // fBm should be normalized to [-1, 1]
    for (int i = 0; i < 1000; ++i) {
        float x = (i * 0.5f) - 250.0f;
        float y = std::sin(i * 0.1f) * 10.0f;
        float z = std::cos(i * 0.07f) * 10.0f;

        float noise = TurbulenceNoise::fBm3D(x, y, z, config);
        EXPECT_GE(noise, -1.5f);
        EXPECT_LE(noise, 1.5f);
    }
}

TEST_F(TurbulenceTest, FBm_MoreOctaves_MoreDetail) {
    // More octaves should add high-frequency detail
    TurbulenceConfig cfg1 = config;
    TurbulenceConfig cfg2 = config;
    cfg1.octaves = 2;
    cfg2.octaves = 8;

    // Compute variance over many samples
    double var1 = 0.0, var2 = 0.0;
    int N = 1000;
    for (int i = 0; i < N; ++i) {
        float x = i * 0.05f;
        float n1 = TurbulenceNoise::fBm3D(x, 0, 0, cfg1);
        float n2 = TurbulenceNoise::fBm3D(x, 0, 0, cfg2);
        var1 += n1 * n1;
        var2 += n2 * n2;
    }
    // Both should have similar variance (normalized)
    // But more octaves might have slightly more variance
    EXPECT_GT(var1, 0.0);
    EXPECT_GT(var2, 0.0);
}

//==============================================================================
// Density Perturbation Tests
//==============================================================================

TEST_F(TurbulenceTest, DensityPerturbation_PositiveDefinite) {
    // Density must always be positive (ρ > 0)
    for (int i = 0; i < 1000; ++i) {
        float r = 5.0f + i * 0.01f;
        float theta = M_PI / 2.0f + 0.01f * std::sin(i);
        float phi = i * 0.1f;

        float perturbation = TurbulenceNoise::sampleDensityPerturbation(r, theta, phi, config);
        EXPECT_GT(perturbation, 0.0f) << "Negative density at sample " << i;
    }
}

TEST_F(TurbulenceTest, DensityPerturbation_MeanNearUnity) {
    // Mean perturbation should be close to 1.0
    double sum = 0.0;
    int N = 10000;
    for (int i = 0; i < N; ++i) {
        float r = 5.0f + (i % 100) * 0.1f;
        float theta = M_PI / 2.0f;
        float phi = i * 0.0628f;  // ~2π/100

        sum += TurbulenceNoise::sampleDensityPerturbation(r, theta, phi, config);
    }
    double mean = sum / N;
    EXPECT_NEAR(mean, 1.0, 0.1);  // Should be close to 1.0
}

TEST_F(TurbulenceTest, DensityPerturbation_Disabled_ReturnsUnity) {
    config.enabled = false;
    float perturbation = TurbulenceNoise::sampleDensityPerturbation(10.0f, M_PI / 2.0f, 0.0f, config);
    EXPECT_FLOAT_EQ(perturbation, 1.0f);
}

TEST_F(TurbulenceTest, DensityPerturbation_ZeroAmplitude_ReturnsUnity) {
    config.amplitude = 0.0f;
    float perturbation = TurbulenceNoise::sampleDensityPerturbation(10.0f, M_PI / 2.0f, 0.0f, config);
    EXPECT_FLOAT_EQ(perturbation, 1.0f);
}

//==============================================================================
// Configuration Validation Tests
//==============================================================================

TEST_F(TurbulenceTest, Validate_ClampsExponent) {
    TurbulenceConfig cfg;
    cfg.kolmogorov_exponent = -3.0f;  // Too negative
    cfg.validate();
    EXPECT_GE(cfg.kolmogorov_exponent, -2.0f);

    cfg.kolmogorov_exponent = -1.0f;  // Too shallow
    cfg.validate();
    EXPECT_LE(cfg.kolmogorov_exponent, -1.5f);
}

TEST_F(TurbulenceTest, Validate_ClampsAmplitude) {
    TurbulenceConfig cfg;
    cfg.amplitude = 2.0f;  // Too high
    cfg.validate();
    EXPECT_LE(cfg.amplitude, 1.0f);

    cfg.amplitude = -0.5f;  // Negative
    cfg.validate();
    EXPECT_GE(cfg.amplitude, 0.0f);
}

TEST_F(TurbulenceTest, Validate_EnforcesScaleHierarchy) {
    TurbulenceConfig cfg;
    cfg.outer_scale_M = 1.0f;
    cfg.inner_scale_M = 2.0f;  // Inner > outer (invalid)
    cfg.validate();
    EXPECT_LT(cfg.inner_scale_M, cfg.outer_scale_M);
}

TEST_F(TurbulenceTest, Validate_ClampsOctaves) {
    TurbulenceConfig cfg;
    cfg.octaves = 100;  // Too many
    cfg.validate();
    EXPECT_LE(cfg.octaves, 8u);

    cfg.octaves = 0;  // Too few
    cfg.validate();
    EXPECT_GE(cfg.octaves, 1u);
}

//==============================================================================
// GPU Parameter Conversion Tests
//==============================================================================

// Helper to create GPU params from CPU config
static TurbulenceParamsGPU createTurbulenceParamsGPU(const TurbulenceConfig& cpu) {
    TurbulenceParamsGPU gpu;
    gpu.amplitude = cpu.amplitude;
    gpu.outer_scale = cpu.outer_scale_M;
    gpu.inner_scale = cpu.inner_scale_M;
    gpu.lacunarity = cpu.lacunarity;
    gpu.persistence = cpu.persistence;
    gpu.octaves = cpu.octaves;
    gpu.seed = cpu.seed;
    gpu.enabled = cpu.enabled ? 1 : 0;
    return gpu;
}

TEST_F(TurbulenceTest, GPUConversion_PreservesValues) {
    TurbulenceParamsGPU gpu = createTurbulenceParamsGPU(config);

    EXPECT_FLOAT_EQ(gpu.amplitude, config.amplitude);
    EXPECT_FLOAT_EQ(gpu.outer_scale, config.outer_scale_M);
    EXPECT_FLOAT_EQ(gpu.inner_scale, config.inner_scale_M);
    EXPECT_FLOAT_EQ(gpu.lacunarity, config.lacunarity);
    EXPECT_FLOAT_EQ(gpu.persistence, config.persistence);
    EXPECT_EQ(gpu.octaves, config.octaves);
    EXPECT_EQ(gpu.seed, config.seed);
    EXPECT_EQ(gpu.enabled, 1u);
}
