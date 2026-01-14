// TSMB001A.cpp - Motion Blur Tests
// Tests for Phase 7: Temporal sampling, camera interpolation, disk rotation

#include <gtest/gtest.h>
#include <cmath>
#include "../../Sirius.Kernel/KNTA001A.h"

using namespace sirius::kernel;
using namespace sirius::spectral;

namespace {

//==============================================================================
// Shutter Configuration Tests
//==============================================================================

TEST(ShutterConfigTest, ShutterTime180Degrees) {
    ShutterConfig config;
    config.shutterAngle = 180.0;  // Half rotation
    config.frameRate = 24.0;
    
    // 180° = 0.5 of frame time = 0.5 / 24 ≈ 0.0208s
    double expected = 0.5 / 24.0;
    EXPECT_NEAR(config.shutterTime(), expected, 1e-6)
        << "180° shutter at 24fps should be ~20.8ms";
}

TEST(ShutterConfigTest, ShutterTime360Degrees) {
    ShutterConfig config;
    config.shutterAngle = 360.0;  // Full rotation
    config.frameRate = 24.0;
    
    // 360° = full frame = 1/24 ≈ 0.0417s
    double expected = 1.0 / 24.0;
    EXPECT_NEAR(config.shutterTime(), expected, 1e-6)
        << "360° shutter should equal frame duration";
}

//==============================================================================
// Temporal Sampler Tests
//==============================================================================

TEST(TemporalSamplerTest, GeneratesCorrectNumberOfSamples) {
    ShutterConfig config;
    config.temporalSamples = 16;
    
    TemporalSampler sampler(config);
    auto samples = sampler.generateSamples();
    
    EXPECT_EQ(samples.size(), 16)
        << "Should generate requested number of samples";
}

TEST(TemporalSamplerTest, SamplesWithinShutterTime) {
    ShutterConfig config;
    config.shutterAngle = 180.0;
    config.frameRate = 24.0;
    config.temporalSamples = 8;
    
    TemporalSampler sampler(config);
    auto samples = sampler.generateSamples();
    
    double shutterTime = config.shutterTime();
    
    for (const auto& s : samples) {
        EXPECT_GE(s.t, 0) << "Sample time should be >= 0";
        EXPECT_LE(s.t, shutterTime) << "Sample time should be <= shutter time";
    }
}

TEST(TemporalSamplerTest, StratifiedDistribution) {
    ShutterConfig config;
    config.temporalSamples = 4;
    
    TemporalSampler sampler(config);
    auto samples = sampler.generateUniformSamples();
    
    // Uniform samples should be evenly spaced
    double shutterTime = config.shutterTime();
    double expectedSpacing = shutterTime / 4;
    
    for (int i = 0; i < 4; ++i) {
        EXPECT_EQ(samples[i].stratumIndex, i)
            << "Stratum index should match";
        EXPECT_NEAR(samples[i].t, (i + 0.5) * expectedSpacing, 1e-10)
            << "Uniform samples should be centered in strata";
    }
}

TEST(TemporalSamplerTest, WeightsSumToOne) {
    ShutterConfig config;
    config.temporalSamples = 8;
    
    TemporalSampler sampler(config);
    auto samples = sampler.generateSamples();
    
    double totalWeight = 0;
    for (const auto& s : samples) {
        totalWeight += s.weight;
    }
    
    EXPECT_NEAR(totalWeight, 1.0, 1e-10)
        << "Sample weights should sum to 1";
}

//==============================================================================
// Temporal Accumulator Tests
//==============================================================================

TEST(TemporalAccumulatorTest, AccumulatesSamples) {
    TemporalAccumulator accumulator;
    
    // Add uniform samples
    SpectralRadiance sample = SpectralRadiance::blackbody(5500);
    
    accumulator.addSample(sample, 0.25);
    accumulator.addSample(sample, 0.25);
    accumulator.addSample(sample, 0.25);
    accumulator.addSample(sample, 0.25);
    
    EXPECT_EQ(accumulator.sampleCount(), 4);
    EXPECT_NEAR(accumulator.totalWeight(), 1.0, 1e-10);
    
    auto result = accumulator.getResult();
    EXPECT_NEAR(result.totalEnergy(), sample.totalEnergy(), 1e-6)
        << "Accumulated result should match original";
}

TEST(TemporalAccumulatorTest, Reset) {
    TemporalAccumulator accumulator;
    
    accumulator.addSample(SpectralRadiance::blackbody(6000));
    EXPECT_EQ(accumulator.sampleCount(), 1);
    
    accumulator.reset();
    EXPECT_EQ(accumulator.sampleCount(), 0);
    EXPECT_DOUBLE_EQ(accumulator.totalWeight(), 0);
}

//==============================================================================
// Camera Worldline Tests
//==============================================================================

TEST(CameraWorldlineTest, SingleKeyframe) {
    CameraWorldline worldline;
    
    CameraState state;
    state.r = 100;
    state.theta = M_PI / 4;
    state.phi = 0;
    
    worldline.addKeyframe(0, state);
    
    // Any time should return the single keyframe
    auto result = worldline.interpolate(0.5);
    EXPECT_DOUBLE_EQ(result.r, 100);
}

TEST(CameraWorldlineTest, LinearInterpolation) {
    CameraWorldline worldline;
    
    CameraState s0, s1;
    s0.r = 100; s0.theta = 0; s0.phi = 0;
    s1.r = 200; s1.theta = M_PI; s1.phi = M_PI;
    
    worldline.addKeyframe(0, s0);
    worldline.addKeyframe(1, s1);
    
    // Interpolate at t = 0.5
    auto result = worldline.interpolate(0.5);
    
    EXPECT_NEAR(result.r, 150, 1e-10)
        << "r should interpolate linearly";
    EXPECT_NEAR(result.theta, M_PI/2, 1e-10)
        << "theta should interpolate linearly";
    EXPECT_NEAR(result.phi, M_PI/2, 1e-10)
        << "phi should interpolate linearly";
}

TEST(CameraWorldlineTest, MultipleKeyframes) {
    CameraWorldline worldline;
    
    CameraState s0, s1, s2;
    s0.r = 100; s1.r = 200; s2.r = 150;
    s0.theta = s1.theta = s2.theta = M_PI/4;
    s0.phi = s1.phi = s2.phi = 0;
    
    worldline.addKeyframe(0, s0);
    worldline.addKeyframe(1, s1);
    worldline.addKeyframe(2, s2);
    
    EXPECT_EQ(worldline.keyframeCount(), 3);
    
    // t = 0.5 should be between s0 and s1
    auto r1 = worldline.interpolate(0.5);
    EXPECT_NEAR(r1.r, 150, 1e-10);
    
    // t = 1.5 should be between s1 and s2
    auto r2 = worldline.interpolate(1.5);
    EXPECT_NEAR(r2.r, 175, 1e-10);
}

//==============================================================================
// Disk Rotation Tests
//==============================================================================

TEST(DiskRotationTest, SchwarzschildKeplerianOmega) {
    DiskRotation rotation(1.0, 0.0);  // Schwarzschild
    
    // Ω = √(M/r³) for Schwarzschild
    // At r = 10M: Ω = √(1/1000) ≈ 0.0316
    double omega = rotation.omega(10.0);
    double expected = std::sqrt(1.0 / 1000.0);
    
    EXPECT_NEAR(omega, expected, 1e-10)
        << "Keplerian omega should be √(M/r³)";
}

TEST(DiskRotationTest, OrbitalPeriod) {
    DiskRotation rotation(1.0, 0.0);
    
    // Period T = 2π/Ω = 2π × √(r³/M)
    // At r = 10M: T = 2π × √1000 ≈ 198.7M
    double T = rotation.period(10.0);
    double expected = 2 * M_PI * std::sqrt(1000.0);
    
    EXPECT_NEAR(T, expected, 1e-6)
        << "Orbital period should be 2π√(r³/M)";
}

TEST(DiskRotationTest, DiskPhiEvolution) {
    DiskRotation rotation;
    
    double r = 10.0;
    double omega = rotation.omega(r);
    
    // After time T, disk should have rotated by Ω×T
    double t = 10.0;
    double phi = rotation.diskPhi(r, t);
    
    EXPECT_NEAR(phi, omega * t, 1e-10)
        << "Disk phi should increase as Ω×t";
}

TEST(DiskRotationTest, KerrOmega) {
    DiskRotation rotation(1.0, 0.9);  // Kerr a = 0.9
    
    // Kerr: Ω = √M / (r^(3/2) + a√M)
    double r = 10.0;
    double expected = std::sqrt(1.0) / (std::pow(10.0, 1.5) + 0.9 * std::sqrt(1.0));
    
    double omega = rotation.omega(r);
    EXPECT_NEAR(omega, expected, 1e-10)
        << "Kerr omega should use correct formula";
}

} // namespace
