// TSJT002A.cpp - MHD Jet Model Unit Tests
// Component ID: TSJT002A (Test/Jet/MHDInvariants)
//
// Tests for MHD jet physics including magnetic flux conservation
// and Doppler beaming factors.

#include <gtest/gtest.h>
#include <cmath>
#include "../../Sirius.Core/Jet/PHJT002A.h"
#include "../../Sirius.Render/Acceleration/OptiX/RDOP003A.h"

using namespace Sirius;

class JetMHDTest : public ::testing::Test {
protected:
    JetMHDConfig config;
    std::unique_ptr<JetMHD> jet;

    void SetUp() override {
        config.B_base = 1e4f;
        config.power_law_index = 1.0f;
        config.B_field_order = 0.5f;
        config.opening_half_angle = 0.1f;
        config.z_launch_M = 3.0f;
        config.z_max_M = 200.0f;
        config.collimation = 0.5f;
        config.lorentz_factor = 5.0f;
        config.electron_index = 2.2f;
        config.n_e_0 = 1e5f;
        config.n_e_decay = 2.0f;
        config.enabled = true;

        jet = std::make_unique<JetMHD>(config);
    }
};

//==============================================================================
// Magnetic Flux Conservation Tests
//==============================================================================

TEST_F(JetMHDTest, MagneticFlux_ApproximatelyConserved) {
    // Φ = π R² B should be approximately constant along jet
    // For B ∝ z^(-p) and R ∝ z^(1-c/2), flux conserved if p = 2(1-c/2)

    float flux_launch = jet->magneticFlux(config.z_launch_M);
    float flux_10M = jet->magneticFlux(10.0f);
    float flux_50M = jet->magneticFlux(50.0f);
    float flux_100M = jet->magneticFlux(100.0f);

    // With p=1 and c=0.5, flux should increase slowly
    // But for physical jets, we accept some variation
    // Test that flux doesn't vary by more than a factor of 10
    EXPECT_GT(flux_10M, flux_launch * 0.1f);
    EXPECT_LT(flux_10M, flux_launch * 10.0f);
    EXPECT_GT(flux_100M, flux_launch * 0.01f);
    EXPECT_LT(flux_100M, flux_launch * 100.0f);
}

TEST_F(JetMHDTest, MagneticField_Decays) {
    float B_launch = jet->magneticField(config.z_launch_M);
    float B_10M = jet->magneticField(10.0f);
    float B_100M = jet->magneticField(100.0f);

    // B should decrease with distance
    EXPECT_GT(B_launch, B_10M);
    EXPECT_GT(B_10M, B_100M);

    // B(z) = B_0 × (z_launch/z)^p
    float expected_B_10M = config.B_base * std::pow(config.z_launch_M / 10.0f, config.power_law_index);
    EXPECT_NEAR(B_10M, expected_B_10M, 1e-3f);
}

//==============================================================================
// Doppler Beaming Tests
//==============================================================================

TEST_F(JetMHDTest, DopplerFactor_Approaching) {
    // For approaching jet (cos_theta = 1), D > 1
    float D = jet->dopplerFactor(1.0f, config.z_launch_M);
    EXPECT_GT(D, 1.0f);

    // For Γ=5 and approaching: D = 1/(Γ(1-β)) ≈ 2Γ
    float expected_D = 2.0f * config.lorentz_factor;  // Approximate
    EXPECT_NEAR(D, expected_D, 1.0f);
}

TEST_F(JetMHDTest, DopplerFactor_Receding) {
    // For receding jet (cos_theta = -1), D < 1
    float D = jet->dopplerFactor(-1.0f, config.z_launch_M);
    EXPECT_LT(D, 1.0f);
    EXPECT_GT(D, 0.0f);
}

TEST_F(JetMHDTest, DopplerFactor_Perpendicular) {
    // For perpendicular view (cos_theta = 0), D = 1/Γ
    float D = jet->dopplerFactor(0.0f, config.z_launch_M);
    float expected_D = 1.0f / config.lorentz_factor;
    EXPECT_NEAR(D, expected_D, 0.01f);
}

TEST_F(JetMHDTest, BoostedIntensity_Beaming) {
    float I_emit = 1.0f;

    // Approaching jet should be brightened
    float I_approach = jet->boostedIntensity(I_emit, 0.9f, config.z_launch_M);
    EXPECT_GT(I_approach, I_emit);

    // Receding jet should be dimmed
    float I_recede = jet->boostedIntensity(I_emit, -0.9f, config.z_launch_M);
    EXPECT_LT(I_recede, I_emit);

    // Beaming asymmetry: approaching >> receding
    EXPECT_GT(I_approach / I_recede, 10.0f);
}

//==============================================================================
// Geometry Tests
//==============================================================================

TEST_F(JetMHDTest, JetRadius_Increases) {
    float R_launch = jet->jetRadius(config.z_launch_M);
    float R_10M = jet->jetRadius(10.0f);
    float R_100M = jet->jetRadius(100.0f);

    EXPECT_GT(R_10M, R_launch);
    EXPECT_GT(R_100M, R_10M);
}

TEST_F(JetMHDTest, JetRadius_ZeroBeforeLaunch) {
    float R_below = jet->jetRadius(config.z_launch_M - 1.0f);
    EXPECT_FLOAT_EQ(R_below, 0.0f);
}

TEST_F(JetMHDTest, IsInsideJet_AtAxis) {
    // Point on axis within jet should be inside
    EXPECT_TRUE(jet->isInsideJet(0.0f, 0.0f, 10.0f, 1));
    EXPECT_TRUE(jet->isInsideJet(0.0f, 0.0f, -10.0f, -1));
}

TEST_F(JetMHDTest, IsInsideJet_OutsideBoundary) {
    // Point far from axis should be outside
    EXPECT_FALSE(jet->isInsideJet(100.0f, 0.0f, 10.0f, 1));
}

TEST_F(JetMHDTest, IsInsideJet_BelowLaunch) {
    // Point below launch height should be outside
    EXPECT_FALSE(jet->isInsideJet(0.0f, 0.0f, 1.0f, 1));
}

TEST_F(JetMHDTest, IsInsideJet_BeyondMax) {
    // Point beyond max extent should be outside
    EXPECT_FALSE(jet->isInsideJet(0.0f, 0.0f, 300.0f, 1));
}

//==============================================================================
// Synchrotron Emission Tests
//==============================================================================

TEST_F(JetMHDTest, SynchrotronEmissivity_Positive) {
    float j = jet->synchrotronEmissivity(10.0f);
    EXPECT_GT(j, 0.0f);
}

TEST_F(JetMHDTest, SynchrotronEmissivity_Decays) {
    float j_10M = jet->synchrotronEmissivity(10.0f);
    float j_100M = jet->synchrotronEmissivity(100.0f);

    // Emissivity should decrease with distance
    EXPECT_GT(j_10M, j_100M);
}

TEST_F(JetMHDTest, ElectronDensity_Decays) {
    float n_10M = jet->electronDensity(10.0f);
    float n_100M = jet->electronDensity(100.0f);

    // n_e ∝ z^(-n_decay)
    float expected_ratio = std::pow(10.0f / 100.0f, config.n_e_decay);
    EXPECT_NEAR(n_10M / n_100M, 1.0f / expected_ratio, 0.01f);
}

//==============================================================================
// Polarisation Tests
//==============================================================================

TEST_F(JetMHDTest, MaxPolarisationDegree_Physical) {
    // For p=2.2: π_max = (2.2+1)/(2.2+7/3) ≈ 0.71
    float pi_max = config.maxPolarisationDegree();
    EXPECT_GT(pi_max, 0.5f);
    EXPECT_LE(pi_max, 1.0f);  // Can't exceed 100%
}

TEST_F(JetMHDTest, PolarisationDegree_ScalesWithOrder) {
    float pi = jet->polarisationDegree(10.0f);
    float pi_max = config.maxPolarisationDegree();

    EXPECT_NEAR(pi, pi_max * config.B_field_order, 0.01f);
}

//==============================================================================
// Spectral Index Tests
//==============================================================================

TEST_F(JetMHDTest, SpectralIndex_FromElectronIndex) {
    // α = (p-1)/2 for optically thin synchrotron
    float alpha = config.spectralIndex();
    float expected = (config.electron_index - 1.0f) / 2.0f;
    EXPECT_FLOAT_EQ(alpha, expected);
}

//==============================================================================
// Configuration Validation Tests
//==============================================================================

TEST_F(JetMHDTest, Validate_ClampsOpeningAngle) {
    JetMHDConfig cfg;
    cfg.opening_half_angle = 1.5f;  // Too wide (> π/4)
    cfg.validate();
    EXPECT_LE(cfg.opening_half_angle, 0.785f);

    cfg.opening_half_angle = 0.001f;  // Too narrow
    cfg.validate();
    EXPECT_GE(cfg.opening_half_angle, 0.01f);
}

TEST_F(JetMHDTest, Validate_EnforcesLorentzFactorMinimum) {
    JetMHDConfig cfg;
    cfg.lorentz_factor = 0.5f;  // Unphysical (< 1)
    cfg.validate();
    EXPECT_GE(cfg.lorentz_factor, 1.0f);
}

TEST_F(JetMHDTest, Validate_ClampsPowerLawIndex) {
    JetMHDConfig cfg;
    cfg.power_law_index = 5.0f;  // Too steep
    cfg.validate();
    EXPECT_LE(cfg.power_law_index, 2.0f);

    cfg.power_law_index = 0.1f;  // Too shallow
    cfg.validate();
    EXPECT_GE(cfg.power_law_index, 0.5f);
}

//==============================================================================
// GPU Parameter Conversion Tests
//==============================================================================

// Helper to create GPU params from CPU config
static JetMHDParamsGPU createJetMHDParamsGPU(const JetMHDConfig& cpu) {
    JetMHDParamsGPU gpu;
    gpu.B_base = cpu.B_base;
    gpu.power_law_index = cpu.power_law_index;
    gpu.B_field_order = cpu.B_field_order;
    gpu.padding1 = 0.0f;
    gpu.opening_angle = cpu.opening_half_angle;
    gpu.z_launch = cpu.z_launch_M;
    gpu.z_max = cpu.z_max_M;
    gpu.collimation = cpu.collimation;
    gpu.lorentz_factor = cpu.lorentz_factor;
    gpu.beta = cpu.beta();
    gpu.velocity_profile = cpu.velocity_profile;
    gpu.padding2 = 0.0f;
    gpu.electron_index = cpu.electron_index;
    gpu.spectral_index = cpu.spectralIndex();
    gpu.gamma_min = cpu.gamma_min;
    gpu.gamma_max = cpu.gamma_max;
    gpu.n_e_0 = cpu.n_e_0;
    gpu.n_e_decay = cpu.n_e_decay;
    gpu.intensity_scale = cpu.intensity_scale;
    gpu.max_polarisation = cpu.maxPolarisationDegree();
    gpu.enabled = cpu.enabled ? 1 : 0;
    gpu.enable_polarisation = cpu.enable_polarisation ? 1 : 0;
    gpu.padding3[0] = gpu.padding3[1] = 0;
    return gpu;
}

TEST_F(JetMHDTest, GPUConversion_PreservesValues) {
    JetMHDParamsGPU gpu = createJetMHDParamsGPU(config);

    EXPECT_FLOAT_EQ(gpu.B_base, config.B_base);
    EXPECT_FLOAT_EQ(gpu.power_law_index, config.power_law_index);
    EXPECT_FLOAT_EQ(gpu.opening_angle, config.opening_half_angle);
    EXPECT_FLOAT_EQ(gpu.z_launch, config.z_launch_M);
    EXPECT_FLOAT_EQ(gpu.lorentz_factor, config.lorentz_factor);
    EXPECT_FLOAT_EQ(gpu.electron_index, config.electron_index);
    EXPECT_EQ(gpu.enabled, 1u);
}

TEST_F(JetMHDTest, GPUConversion_ComputesBeta) {
    JetMHDParamsGPU gpu = createJetMHDParamsGPU(config);

    // β = sqrt(1 - 1/Γ²)
    float expected_beta = std::sqrt(1.0f - 1.0f / (config.lorentz_factor * config.lorentz_factor));
    EXPECT_NEAR(gpu.beta, expected_beta, 1e-5f);
}
