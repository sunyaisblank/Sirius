// TSSM001A.cpp - SMBH Parameters Unit Tests
// Component ID: TSSM001A (Test/SMBH/Invariants)
//
// Tests for supermassive black hole astrophysical scaling and invariants.

#include <gtest/gtest.h>
#include <cmath>
#include "../../Sirius.Core/Metric/PHSM001A.h"
#include "../../Sirius.Render/Acceleration/OptiX/RDOP003A.h"

using namespace Sirius;

class SMBHParamsTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

//==============================================================================
// ISCO Radius Tests
//==============================================================================

TEST_F(SMBHParamsTest, ISCO_Schwarzschild) {
    // For a = 0 (Schwarzschild), ISCO should be exactly 6M
    double r_isco = SMBHParams::computeISCO(0.0f);
    EXPECT_NEAR(r_isco, 6.0, 1e-10);
}

TEST_F(SMBHParamsTest, ISCO_Extremal_Prograde) {
    // For a = 1 (extremal prograde), ISCO approaches 1M
    double r_isco = SMBHParams::computeISCO(0.998f);
    EXPECT_GT(r_isco, 1.0);
    EXPECT_LT(r_isco, 2.0);
}

TEST_F(SMBHParamsTest, ISCO_Extremal_Retrograde) {
    // For a = -1 (extremal retrograde), ISCO approaches 9M
    double r_isco = SMBHParams::computeISCO(-0.998f);
    EXPECT_GT(r_isco, 8.0);
    EXPECT_LT(r_isco, 10.0);
}

TEST_F(SMBHParamsTest, ISCO_Monotonic_In_Spin) {
    // ISCO should decrease monotonically with increasing prograde spin
    double prev_isco = SMBHParams::computeISCO(-0.9f);
    for (float a = -0.8f; a <= 0.9f; a += 0.1f) {
        double isco = SMBHParams::computeISCO(a);
        EXPECT_LE(isco, prev_isco + 1e-6);  // Allow small numerical tolerance
        prev_isco = isco;
    }
}

//==============================================================================
// Horizon Radius Tests
//==============================================================================

TEST_F(SMBHParamsTest, Horizon_Schwarzschild) {
    SMBHParams params;
    params.mass_solar = 1.0e6;
    params.spin_parameter = 0.0f;
    params.computeDerived();

    // For Schwarzschild, r_horizon = 2M
    EXPECT_NEAR(params.r_horizon_M, 2.0, 1e-10);
}

TEST_F(SMBHParamsTest, Horizon_Extremal) {
    SMBHParams params;
    params.mass_solar = 1.0e6;
    params.spin_parameter = 0.998f;
    params.computeDerived();

    // For extremal Kerr, r_horizon approaches M
    EXPECT_NEAR(params.r_horizon_M, 1.0, 0.1);
}

TEST_F(SMBHParamsTest, InnerHorizon_Schwarzschild) {
    // For Schwarzschild, inner horizon = 0 (coincides with singularity)
    double r_inner = SMBHParams::computeInnerHorizon(0.0f);
    EXPECT_NEAR(r_inner, 0.0, 1e-10);
}

//==============================================================================
// Mass Scaling Tests
//==============================================================================

TEST_F(SMBHParamsTest, GravitationalRadius_Scaling) {
    SMBHParams params1, params2;
    params1.mass_solar = 1.0e6;
    params2.mass_solar = 1.0e8;
    params1.computeDerived();
    params2.computeDerived();

    // r_g should scale linearly with mass
    EXPECT_NEAR(params2.r_g / params1.r_g, 100.0, 1e-6);
}

TEST_F(SMBHParamsTest, EddingtonLuminosity_Scaling) {
    SMBHParams params1, params2;
    params1.mass_solar = 1.0e6;
    params2.mass_solar = 1.0e8;
    params1.computeDerived();
    params2.computeDerived();

    // L_Edd should scale linearly with mass
    EXPECT_NEAR(params2.L_Edd / params1.L_Edd, 100.0, 1e-6);
}

//==============================================================================
// Thorne Limit Tests
//==============================================================================

TEST_F(SMBHParamsTest, SpinClamped_To_ThorneLimitPositive) {
    SMBHParams params;
    params.spin_parameter = 1.5f;  // Exceeds Thorne limit
    params.computeDerived();

    EXPECT_LE(params.spin_parameter, 0.998f);
}

TEST_F(SMBHParamsTest, SpinClamped_To_ThorneLimitNegative) {
    SMBHParams params;
    params.spin_parameter = -1.5f;  // Exceeds Thorne limit
    params.computeDerived();

    EXPECT_GE(params.spin_parameter, -0.998f);
}

//==============================================================================
// Preset Configuration Tests
//==============================================================================

TEST_F(SMBHParamsTest, M87Star_Configuration) {
    SMBHParams m87 = SMBHParams::M87Star();

    EXPECT_NEAR(m87.mass_solar, 6.5e9, 1e7);
    EXPECT_NEAR(m87.spin_parameter, 0.90f, 0.01f);
    EXPECT_NEAR(m87.inclination_deg, 17.0f, 0.1f);

    // Derived quantities should be computed
    EXPECT_GT(m87.r_g, 0.0);
    EXPECT_GT(m87.r_isco, 0.0);
    EXPECT_LT(m87.r_isco_M, 6.0);  // Spinning BH has smaller ISCO
}

TEST_F(SMBHParamsTest, SgrAStar_Configuration) {
    SMBHParams sgra = SMBHParams::SgrAStar();

    EXPECT_NEAR(sgra.mass_solar, 4.0e6, 1e4);
    EXPECT_NEAR(sgra.spin_parameter, 0.50f, 0.01f);
}

TEST_F(SMBHParamsTest, Gargantua_NearExtremal) {
    SMBHParams garg = SMBHParams::Gargantua();

    // Near-extremal spin
    EXPECT_GT(garg.spin_parameter, 0.99f);
    EXPECT_LE(garg.spin_parameter, 0.998f);

    // ISCO should be very close to horizon
    EXPECT_LT(garg.r_isco_M, 2.0);
}

//==============================================================================
// GPU Parameter Conversion Tests
//==============================================================================

// Helper to create GPU params from CPU params
static SMBHParamsGPU createSMBHParamsGPU(const SMBHParams& cpu) {
    SMBHParamsGPU gpu;
    gpu.mass_M = 1.0f;  // Normalized
    gpu.spin = cpu.spin_parameter;
    gpu.inclination_rad = cpu.inclination_deg * 3.14159265f / 180.0f;
    gpu.r_isco = static_cast<float>(cpu.r_isco_M);
    gpu.r_horizon = static_cast<float>(cpu.r_horizon_M);
    gpu.r_inner_horizon = static_cast<float>(SMBHParams::computeInnerHorizon(cpu.spin_parameter));
    gpu.angular_size_rg = static_cast<float>(cpu.angularSizeOfRg());
    gpu.padding = 0.0f;
    return gpu;
}

TEST_F(SMBHParamsTest, GPUConversion_PreservesValues) {
    SMBHParams cpu;
    cpu.mass_solar = 1.0e8;
    cpu.spin_parameter = 0.7f;
    cpu.inclination_deg = 45.0f;
    cpu.distance_Mpc = 10.0f;
    cpu.computeDerived();

    SMBHParamsGPU gpu = createSMBHParamsGPU(cpu);

    EXPECT_FLOAT_EQ(gpu.mass_M, 1.0f);
    EXPECT_FLOAT_EQ(gpu.spin, 0.7f);
    EXPECT_NEAR(gpu.inclination_rad, 0.785398f, 0.001f);  // π/4
    EXPECT_NEAR(gpu.r_isco, cpu.r_isco_M, 1e-4f);
    EXPECT_NEAR(gpu.r_horizon, cpu.r_horizon_M, 1e-4f);
}

//==============================================================================
// Angular Size Tests
//==============================================================================

TEST_F(SMBHParamsTest, AngularSize_InverseSqaureDistance) {
    SMBHParams params1, params2;
    params1.mass_solar = 1.0e8;
    params1.distance_Mpc = 10.0f;
    params2.mass_solar = 1.0e8;
    params2.distance_Mpc = 20.0f;
    params1.computeDerived();
    params2.computeDerived();

    // Angular size should scale as 1/distance
    EXPECT_NEAR(params1.angularSizeOfRg() / params2.angularSizeOfRg(), 2.0, 1e-6);
}
