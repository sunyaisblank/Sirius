// TSJT001A.cpp - Relativistic Jet Emission Model Tests
// Component ID: TSJT001A (Test/Jet/JetEmission)
// Tests: PHJT001A.h (RelativisticJet, JetRayMarching)

#include <gtest/gtest.h>
#include <cmath>
#include <PHJT001A.h>

namespace sirius::test {
using namespace Sirius;

constexpr float kEps = 1e-4f;

// =============================================================================
// Doppler Factor Tests
// =============================================================================

class JetDopplerTests : public ::testing::Test {
protected:
    Sirius::RelativisticJet jet;
    float beta;
    float gamma;

    void SetUp() override {
        Sirius::JetConfig config;
        config.lorentz_factor = 5.0f;
        jet.setConfig(config);
        gamma = 5.0f;
        beta = std::sqrt(1.0f - 1.0f / (gamma * gamma));
    }
};

TEST_F(JetDopplerTests, HeadOnApproachBoosted) {
    // cos_theta = 1 (looking directly at approaching jet)
    float D = jet.dopplerFactor(1.0f);
    // D = 1 / [Gamma * (1 - beta)] > 1
    float expected = 1.0f / (gamma * (1.0f - beta));
    EXPECT_NEAR(D, expected, kEps);
    EXPECT_GT(D, 1.0f);
}

TEST_F(JetDopplerTests, RecedingDeBoosted) {
    // cos_theta = -1 (looking at receding jet)
    float D = jet.dopplerFactor(-1.0f);
    float expected = 1.0f / (gamma * (1.0f + beta));
    EXPECT_NEAR(D, expected, kEps);
    EXPECT_LT(D, 1.0f);
}

TEST_F(JetDopplerTests, TransverseDirection) {
    // cos_theta = 0 (perpendicular)
    float D = jet.dopplerFactor(0.0f);
    float expected = 1.0f / gamma;
    EXPECT_NEAR(D, expected, kEps);
}

TEST_F(JetDopplerTests, AnalyticFormula) {
    // Test at arbitrary angle
    float cos_theta = 0.6f;
    float D = jet.dopplerFactor(cos_theta);
    float expected = 1.0f / (gamma * (1.0f - beta * cos_theta));
    EXPECT_NEAR(D, expected, kEps);
}

// =============================================================================
// Geometry Tests
// =============================================================================

class JetGeometryTests : public ::testing::Test {
protected:
    Sirius::RelativisticJet jet;
};

TEST_F(JetGeometryTests, InsideJetOnAxis) {
    // Point on z-axis at height above r_launch
    float h = jet.getConfig().r_launch + 5.0f;
    EXPECT_TRUE(jet.isInsideJet(0.0f, 0.0f, h, 1));
}

TEST_F(JetGeometryTests, OutsideJetBelowLaunch) {
    EXPECT_FALSE(jet.isInsideJet(0.0f, 0.0f, 1.0f, 1));
}

TEST_F(JetGeometryTests, OutsideJetFarOffAxis) {
    float h = jet.getConfig().r_launch + 5.0f;
    EXPECT_FALSE(jet.isInsideJet(100.0f, 0.0f, h, 1));
}

TEST_F(JetGeometryTests, SouthernJet) {
    float h = jet.getConfig().r_launch + 5.0f;
    EXPECT_TRUE(jet.isInsideJet(0.0f, 0.0f, -h, -1));
    EXPECT_FALSE(jet.isInsideJet(0.0f, 0.0f, -h, 1)); // wrong sign
}

TEST_F(JetGeometryTests, JetRadiusMonotone) {
    float r_launch = jet.getConfig().r_launch;
    float prev = jet.jetRadius(r_launch);
    for (float h = r_launch + 1.0f; h <= 50.0f; h += 1.0f) {
        float r = jet.jetRadius(h);
        EXPECT_GE(r, prev - kEps) << "Jet radius not monotone at h = " << h;
        prev = r;
    }
}

TEST_F(JetGeometryTests, JetRadiusZeroBelowLaunch) {
    EXPECT_FLOAT_EQ(jet.jetRadius(0.0f), 0.0f);
    EXPECT_FLOAT_EQ(jet.jetRadius(jet.getConfig().r_launch), 0.0f);
}

// =============================================================================
// Emission Tests
// =============================================================================

class JetEmissionTests : public ::testing::Test {
protected:
    Sirius::RelativisticJet jet;
};

TEST_F(JetEmissionTests, MagneticFieldAtLaunch) {
    EXPECT_FLOAT_EQ(jet.magneticField(jet.getConfig().r_launch),
                    jet.getConfig().B_field_0);
}

TEST_F(JetEmissionTests, MagneticFieldDecays) {
    float B1 = jet.magneticField(10.0f);
    float B2 = jet.magneticField(50.0f);
    EXPECT_GT(B1, B2);
}

TEST_F(JetEmissionTests, ElectronDensityAtLaunch) {
    EXPECT_FLOAT_EQ(jet.electronDensity(jet.getConfig().r_launch),
                    jet.getConfig().n_e_0);
}

TEST_F(JetEmissionTests, ElectronDensityDecays) {
    float n1 = jet.electronDensity(10.0f);
    float n2 = jet.electronDensity(50.0f);
    EXPECT_GT(n1, n2);
}

TEST_F(JetEmissionTests, SynchrotronEmissivityPositive) {
    float j = jet.synchrotronEmissivity(10.0f);
    EXPECT_GT(j, 0.0f);
}

TEST_F(JetEmissionTests, BeamingApproachingBrighterThanReceding) {
    float I_emit = 1.0f;
    float I_approaching = jet.boostedIntensity(I_emit, 1.0f, true);
    float I_receding = jet.boostedIntensity(I_emit, -1.0f, true);
    EXPECT_GT(I_approaching, I_receding);
}

TEST_F(JetEmissionTests, PolarisationDegreeFormula) {
    float p = jet.getConfig().spectral_index;
    float B_order = jet.getConfig().B_field_order;
    float expected = (p + 1.0f) / (p + 7.0f / 3.0f) * B_order;
    EXPECT_NEAR(jet.polarisationDegree(), expected, kEps);
}

TEST_F(JetEmissionTests, VelocityDirection) {
    float vx, vy, vz;
    jet.getVelocity(vx, vy, vz, 1);
    EXPECT_FLOAT_EQ(vx, 0.0f);
    EXPECT_FLOAT_EQ(vy, 0.0f);
    EXPECT_GT(vz, 0.0f);

    jet.getVelocity(vx, vy, vz, -1);
    EXPECT_LT(vz, 0.0f);
}

// =============================================================================
// Ray Marching Tests
// =============================================================================

class JetRayMarchTests : public ::testing::Test {};

TEST_F(JetRayMarchTests, EmissionOutsideJetIsZero) {
    Sirius::RelativisticJet jet;
    // Ray entirely outside jet volume
    float emission = Sirius::JetRayMarching::integrateJetEmission(
        jet,
        100.0f, 100.0f, 0.0f,   // start far off-axis
        100.0f, 100.0f, 100.0f, // end far off-axis
        0.0f, 0.0f, 1000.0f,    // observer
        16
    );
    EXPECT_FLOAT_EQ(emission, 0.0f);
}

} // namespace sirius::test
