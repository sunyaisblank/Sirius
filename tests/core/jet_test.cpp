// Relativistic jet emission model tests: Doppler factor, geometry, emission,
// ray marching. Ported from TSJT001A.cpp; assertions and tolerances unchanged.

#include "sirius/core/jet.h"

#include <gtest/gtest.h>

#include <cmath>

namespace sirius::test {
using namespace sirius::core;

constexpr float kEps = 1e-4f;

// =============================================================================
// Doppler Factor Tests
// =============================================================================

class JetDopplerTests : public ::testing::Test {
protected:
    sirius::core::RelativisticJet jet;
    float beta;
    float gamma;

    void SetUp() override {
        sirius::core::JetConfig config;
        config.lorentz_factor = 5.0f;
        jet.SetConfig(config);
        gamma = 5.0f;
        beta = std::sqrt(1.0f - 1.0f / (gamma * gamma));
    }
};

TEST_F(JetDopplerTests, HeadOnApproachBoosted) {
    // cos_theta = 1 (looking directly at approaching jet)
    float D = jet.DopplerFactor(1.0f);
    // D = 1 / [Gamma * (1 - beta)] > 1
    float expected = 1.0f / (gamma * (1.0f - beta));
    EXPECT_NEAR(D, expected, kEps);
    EXPECT_GT(D, 1.0f);
}

TEST_F(JetDopplerTests, RecedingDeBoosted) {
    // cos_theta = -1 (looking at receding jet)
    float D = jet.DopplerFactor(-1.0f);
    float expected = 1.0f / (gamma * (1.0f + beta));
    EXPECT_NEAR(D, expected, kEps);
    EXPECT_LT(D, 1.0f);
}

TEST_F(JetDopplerTests, TransverseDirection) {
    // cos_theta = 0 (perpendicular)
    float D = jet.DopplerFactor(0.0f);
    float expected = 1.0f / gamma;
    EXPECT_NEAR(D, expected, kEps);
}

TEST_F(JetDopplerTests, AnalyticFormula) {
    // Test at arbitrary angle
    float cos_theta = 0.6f;
    float D = jet.DopplerFactor(cos_theta);
    float expected = 1.0f / (gamma * (1.0f - beta * cos_theta));
    EXPECT_NEAR(D, expected, kEps);
}

// =============================================================================
// Geometry Tests
// =============================================================================

class JetGeometryTests : public ::testing::Test {
protected:
    sirius::core::RelativisticJet jet;
};

TEST_F(JetGeometryTests, InsideJetOnAxis) {
    // Point on z-axis at height above r_launch
    float h = jet.GetConfig().r_launch + 5.0f;
    EXPECT_TRUE(jet.IsInsideJet(0.0f, 0.0f, h, 1));
}

TEST_F(JetGeometryTests, OutsideJetBelowLaunch) {
    EXPECT_FALSE(jet.IsInsideJet(0.0f, 0.0f, 1.0f, 1));
}

TEST_F(JetGeometryTests, OutsideJetFarOffAxis) {
    float h = jet.GetConfig().r_launch + 5.0f;
    EXPECT_FALSE(jet.IsInsideJet(100.0f, 0.0f, h, 1));
}

TEST_F(JetGeometryTests, SouthernJet) {
    float h = jet.GetConfig().r_launch + 5.0f;
    EXPECT_TRUE(jet.IsInsideJet(0.0f, 0.0f, -h, -1));
    EXPECT_FALSE(jet.IsInsideJet(0.0f, 0.0f, -h, 1)); // wrong sign
}

TEST_F(JetGeometryTests, JetRadiusMonotone) {
    float r_launch = jet.GetConfig().r_launch;
    float prev = jet.JetRadius(r_launch);
    for (float h = r_launch + 1.0f; h <= 50.0f; h += 1.0f) {
        float r = jet.JetRadius(h);
        EXPECT_GE(r, prev - kEps) << "Jet radius not monotone at h = " << h;
        prev = r;
    }
}

TEST_F(JetGeometryTests, JetRadiusZeroBelowLaunch) {
    EXPECT_FLOAT_EQ(jet.JetRadius(0.0f), 0.0f);
    EXPECT_FLOAT_EQ(jet.JetRadius(jet.GetConfig().r_launch), 0.0f);
}

// =============================================================================
// Emission Tests
// =============================================================================

class JetEmissionTests : public ::testing::Test {
protected:
    sirius::core::RelativisticJet jet;
};

TEST_F(JetEmissionTests, MagneticFieldAtLaunch) {
    EXPECT_FLOAT_EQ(jet.MagneticField(jet.GetConfig().r_launch),
                    jet.GetConfig().B_field_0);
}

TEST_F(JetEmissionTests, MagneticFieldDecays) {
    float B1 = jet.MagneticField(10.0f);
    float B2 = jet.MagneticField(50.0f);
    EXPECT_GT(B1, B2);
}

TEST_F(JetEmissionTests, ElectronDensityAtLaunch) {
    EXPECT_FLOAT_EQ(jet.ElectronDensity(jet.GetConfig().r_launch),
                    jet.GetConfig().n_e_0);
}

TEST_F(JetEmissionTests, ElectronDensityDecays) {
    float n1 = jet.ElectronDensity(10.0f);
    float n2 = jet.ElectronDensity(50.0f);
    EXPECT_GT(n1, n2);
}

TEST_F(JetEmissionTests, SynchrotronEmissivityPositive) {
    float j = jet.SynchrotronEmissivity(10.0f);
    EXPECT_GT(j, 0.0f);
}

TEST_F(JetEmissionTests, BeamingApproachingBrighterThanReceding) {
    float I_emit = 1.0f;
    float I_approaching = jet.BoostedIntensity(I_emit, 1.0f, true);
    float I_receding = jet.BoostedIntensity(I_emit, -1.0f, true);
    EXPECT_GT(I_approaching, I_receding);
}

TEST_F(JetEmissionTests, PolarisationDegreeFormula) {
    float p = jet.GetConfig().spectral_index;
    float B_order = jet.GetConfig().B_field_order;
    float expected = (p + 1.0f) / (p + 7.0f / 3.0f) * B_order;
    EXPECT_NEAR(jet.PolarisationDegree(), expected, kEps);
}

TEST_F(JetEmissionTests, VelocityDirection) {
    float vx, vy, vz;
    jet.GetVelocity(vx, vy, vz, 1);
    EXPECT_FLOAT_EQ(vx, 0.0f);
    EXPECT_FLOAT_EQ(vy, 0.0f);
    EXPECT_GT(vz, 0.0f);

    jet.GetVelocity(vx, vy, vz, -1);
    EXPECT_LT(vz, 0.0f);
}

// =============================================================================
// Ray Marching Tests
// =============================================================================

class JetRayMarchTests : public ::testing::Test {};

TEST_F(JetRayMarchTests, EmissionOutsideJetIsZero) {
    sirius::core::RelativisticJet jet;
    // Ray entirely outside jet volume
    float emission = sirius::core::jet_ray_marching::IntegrateJetEmission(
        jet,
        100.0f, 100.0f, 0.0f,   // start far off-axis
        100.0f, 100.0f, 100.0f, // end far off-axis
        0.0f, 0.0f, 1000.0f,    // observer
        16
    );
    EXPECT_FLOAT_EQ(emission, 0.0f);
}

} // namespace sirius::test
