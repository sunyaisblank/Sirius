// TSDK002A.cpp - Volumetric Disk Physics Tests
// Component ID: TSDK002A
// Tests for: PHAD002A.h (VolumetricDiskD)
//
// Tests verify:
// 1. Scale height computation
// 2. Density profile (Gaussian vertical)
// 3. Temperature profile (vertical gradient)
// 4. Optical depth integration
// 5. Photosphere height detection

#include <gtest/gtest.h>
#include <cmath>

#include "PHAD002A.h"

using namespace sirius::physics;

// =============================================================================
// Test Fixture
// =============================================================================

class VolumetricDiskTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Default configuration for tests
        VolumetricDiskD::Config config;
        config.M = 1.0;
        config.a_star = 0.0;
        config.Mdot = 1e-8;
        config.r_inner = 6.0;  // ISCO for Schwarzschild
        config.r_outer = 100.0;
        config.H_over_r = 0.1;
        config.r_ref = 10.0;
        config.H_power = 0.25;
        config.tau_midplane = 10.0;
        config.tau_power = -1.5;
        config.T_atm_ratio = 0.8;
        config.z_truncation = 3.0;

        disk = std::make_unique<VolumetricDiskD>(config);
    }

    std::unique_ptr<VolumetricDiskD> disk;
};

// =============================================================================
// Scale Height Tests
// =============================================================================

TEST_F(VolumetricDiskTest, ScaleHeightAtReferenceRadius) {
    // At r_ref = 10, H should be H_over_r × r_ref = 0.1 × 10 = 1.0
    double H = disk->scaleHeight(10.0);
    EXPECT_NEAR(H, 1.0, 1e-6);
}

TEST_F(VolumetricDiskTest, ScaleHeightFlaring) {
    // H/r ∝ r^H_power with H_power = 0.25
    // At r = 16 (16/10 = 1.6), H/r should be 0.1 × 1.6^0.25 ≈ 0.1124
    // So H ≈ 0.1124 × 16 ≈ 1.80

    double H_at_16 = disk->scaleHeight(16.0);
    double expected_H_over_r = 0.1 * std::pow(16.0 / 10.0, 0.25);
    double expected_H = expected_H_over_r * 16.0;

    EXPECT_NEAR(H_at_16, expected_H, 0.01);
}

TEST_F(VolumetricDiskTest, ScaleHeightMonotonicallyIncreases) {
    // For H_power > 0, H should increase with r
    double H_6 = disk->scaleHeight(6.0);
    double H_10 = disk->scaleHeight(10.0);
    double H_20 = disk->scaleHeight(20.0);
    double H_50 = disk->scaleHeight(50.0);

    EXPECT_LT(H_6, H_10);
    EXPECT_LT(H_10, H_20);
    EXPECT_LT(H_20, H_50);
}

// =============================================================================
// Density Profile Tests
// =============================================================================

TEST_F(VolumetricDiskTest, DensityMaxAtMidplane) {
    // Density should be maximum at z = 0
    double rho_midplane = disk->density(10.0, 0.0);
    double rho_above = disk->density(10.0, 0.5);
    double rho_below = disk->density(10.0, -0.5);

    EXPECT_GT(rho_midplane, rho_above);
    EXPECT_GT(rho_midplane, rho_below);
}

TEST_F(VolumetricDiskTest, DensityGaussianSymmetry) {
    // Density should be symmetric about z = 0
    double rho_plus = disk->density(10.0, 0.5);
    double rho_minus = disk->density(10.0, -0.5);

    EXPECT_NEAR(rho_plus, rho_minus, 1e-10);
}

TEST_F(VolumetricDiskTest, DensityGaussianProfile) {
    // Density at z = H should be exp(-0.5) ≈ 0.606 × density at z = 0
    double H = disk->scaleHeight(10.0);
    double rho_0 = disk->density(10.0, 0.0);
    double rho_H = disk->density(10.0, H);

    double expected_ratio = std::exp(-0.5);
    double actual_ratio = rho_H / rho_0;

    EXPECT_NEAR(actual_ratio, expected_ratio, 0.01);
}

TEST_F(VolumetricDiskTest, DensityZeroOutsideDisk) {
    // Density should be zero outside disk bounds
    double rho_inner = disk->density(3.0, 0.0);   // Inside ISCO
    double rho_outer = disk->density(150.0, 0.0); // Beyond outer radius

    EXPECT_EQ(rho_inner, 0.0);
    EXPECT_EQ(rho_outer, 0.0);
}

TEST_F(VolumetricDiskTest, DensityZeroBeyondTruncation) {
    // Density should be zero beyond 3σ truncation
    double H = disk->scaleHeight(10.0);
    double z_max = 3.0 * H;

    double rho_at_3H = disk->density(10.0, z_max);
    double rho_beyond = disk->density(10.0, z_max + 0.1);

    EXPECT_GT(rho_at_3H, 0.0);
    EXPECT_EQ(rho_beyond, 0.0);
}

// =============================================================================
// Temperature Profile Tests
// =============================================================================

TEST_F(VolumetricDiskTest, TemperatureMaxAtMidplane) {
    // Temperature should be maximum at midplane
    double T_mid = disk->temperature(10.0, 0.0);
    double T_above = disk->temperature(10.0, 0.5);

    EXPECT_GT(T_mid, T_above);
}

TEST_F(VolumetricDiskTest, TemperatureVerticalGradient) {
    // Temperature at surface (z = 3H) should be T_atm_ratio × T_mid
    double H = disk->scaleHeight(10.0);
    double T_mid = disk->temperature(10.0, 0.0);
    double T_surface = disk->temperature(10.0, 3.0 * H);

    // T_atm_ratio = 0.8
    double expected_ratio = 0.8;
    double actual_ratio = T_surface / T_mid;

    EXPECT_NEAR(actual_ratio, expected_ratio, 0.05);
}

TEST_F(VolumetricDiskTest, TemperatureDecreasesWithRadius) {
    // Temperature should decrease outward (Novikov-Thorne profile)
    // Note: T → 0 at ISCO (r=6), so we test at radii > ISCO
    // Peak temperature is around r ~ 1.36 × ISCO ≈ 8.2M for Schwarzschild
    double T_8 = disk->temperature(8.0, 0.0);
    double T_10 = disk->temperature(10.0, 0.0);
    double T_20 = disk->temperature(20.0, 0.0);

    // For Novikov-Thorne, peak is around r ~ 9.5M, so T_8 might be < T_10
    // At larger radii, temperature always decreases
    EXPECT_GT(T_10, T_20);
    EXPECT_GT(T_20, 0.0);
}

// =============================================================================
// Optical Depth Tests
// =============================================================================

TEST_F(VolumetricDiskTest, VerticalOpticalDepthAtMidplane) {
    // τ(r_ref, 0) should be tau_midplane / 2 (one side only)
    double tau = disk->verticalOpticalDepth(10.0, 0.0);

    EXPECT_NEAR(tau, 5.0, 0.5);  // Half of tau_midplane = 10
}

TEST_F(VolumetricDiskTest, VerticalOpticalDepthDecreaseswithZ) {
    // τ should decrease as we move away from midplane
    double tau_0 = disk->verticalOpticalDepth(10.0, 0.0);
    double tau_H = disk->verticalOpticalDepth(10.0, 1.0);  // z = H
    double tau_2H = disk->verticalOpticalDepth(10.0, 2.0); // z = 2H

    EXPECT_GT(tau_0, tau_H);
    EXPECT_GT(tau_H, tau_2H);
}

TEST_F(VolumetricDiskTest, VerticalOpticalDepthRadialScaling) {
    // τ_midplane ∝ r^(-1.5)
    double tau_10 = disk->verticalOpticalDepth(10.0, 0.0);
    double tau_20 = disk->verticalOpticalDepth(20.0, 0.0);

    // At 20M vs 10M, expect (20/10)^(-1.5) = 2^(-1.5) ≈ 0.354
    double expected_ratio = std::pow(2.0, -1.5);
    double actual_ratio = tau_20 / tau_10;

    EXPECT_NEAR(actual_ratio, expected_ratio, 0.1);
}

// =============================================================================
// Photosphere Height Tests
// =============================================================================

TEST_F(VolumetricDiskTest, PhotosphereDefined) {
    // For τ_midplane = 10 >> 2/3, photosphere should be well-defined
    double z_phot = disk->photosphereHeight(10.0);

    EXPECT_GT(z_phot, 0.0);
}

TEST_F(VolumetricDiskTest, PhotosphereBelowSurface) {
    // Photosphere should be below disk surface (3H)
    double H = disk->scaleHeight(10.0);
    double z_phot = disk->photosphereHeight(10.0);

    EXPECT_LT(z_phot, 3.0 * H);
}

// =============================================================================
// Boundary Tests
// =============================================================================

TEST_F(VolumetricDiskTest, IsInDiskVolume) {
    // Test various positions
    EXPECT_TRUE(disk->isInDiskVolume(10.0, 0.0));   // Midplane
    EXPECT_TRUE(disk->isInDiskVolume(10.0, 0.5));   // Above midplane

    double H = disk->scaleHeight(10.0);
    EXPECT_TRUE(disk->isInDiskVolume(10.0, 2.9 * H));  // Near surface
    EXPECT_FALSE(disk->isInDiskVolume(10.0, 3.1 * H)); // Beyond surface

    EXPECT_FALSE(disk->isInDiskVolume(3.0, 0.0));   // Inside ISCO
    EXPECT_FALSE(disk->isInDiskVolume(150.0, 0.0)); // Beyond outer radius
}

TEST_F(VolumetricDiskTest, SurfaceHeightConsistent) {
    // Surface height should be z_truncation × H
    double H = disk->scaleHeight(10.0);
    double z_surface = disk->surfaceHeight(10.0);

    EXPECT_NEAR(z_surface, 3.0 * H, 1e-6);
}

// =============================================================================
// Source Function Tests
// =============================================================================

TEST_F(VolumetricDiskTest, SourceFunctionPositive) {
    // Source function should be positive inside disk
    double S = disk->sourceFunction(10.0, 0.0);
    EXPECT_GT(S, 0.0);
}

TEST_F(VolumetricDiskTest, SourceFunctionDecreasesWithZ) {
    // Source function ∝ T^4, should decrease with |z|
    double S_0 = disk->sourceFunction(10.0, 0.0);
    double S_H = disk->sourceFunction(10.0, 1.0);

    EXPECT_GT(S_0, S_H);
}

// =============================================================================
// Edge Cases
// =============================================================================

TEST_F(VolumetricDiskTest, ZeroRadiusHandled) {
    // Should not crash at r = 0
    double H = disk->scaleHeight(0.0);
    EXPECT_EQ(H, 0.0);
}

TEST_F(VolumetricDiskTest, NegativeZHandled) {
    // Should work for z < 0 (below midplane)
    double rho = disk->density(10.0, -1.0);
    EXPECT_GT(rho, 0.0);

    double T = disk->temperature(10.0, -1.0);
    EXPECT_GT(T, 0.0);
}
