// Corona model tests. Ported from TSCR001A.cpp; assertions and tolerances
// unchanged. Tests corona.h (CoronaConfig, corona_physics).

#include "sirius/core/disk/corona.h"

#include <gtest/gtest.h>

#include <cmath>

namespace sirius::test {
using namespace sirius::core;

constexpr float kEps = 1e-4f;
constexpr float kPi = 3.14159265358979f;

// =============================================================================
// CoronaConfig Tests
// =============================================================================

class CoronaConfigTests : public ::testing::Test {};

TEST_F(CoronaConfigTests, ValidateClampsTemperature) {
    sirius::core::CoronaConfig c;
    c.temperature_keV = 1.0f; // below minimum 10
    c.Validate();
    EXPECT_GE(c.temperature_keV, 10.0f);

    c.temperature_keV = 1000.0f; // above maximum 500
    c.Validate();
    EXPECT_LE(c.temperature_keV, 500.0f);
}

TEST_F(CoronaConfigTests, ValidateClampsOpticalDepth) {
    sirius::core::CoronaConfig c;
    c.optical_depth = 0.01f;
    c.Validate();
    EXPECT_GE(c.optical_depth, 0.1f);

    c.optical_depth = 10.0f;
    c.Validate();
    EXPECT_LE(c.optical_depth, 5.0f);
}

TEST_F(CoronaConfigTests, ValidateEnsuresOuterGreaterThanInner) {
    sirius::core::CoronaConfig c;
    c.inner_radius_M = 10.0f;
    c.outer_radius_M = 5.0f; // invalid: outer < inner
    c.Validate();
    EXPECT_GT(c.outer_radius_M, c.inner_radius_M);
}

TEST_F(CoronaConfigTests, ComptonizationParameter) {
    sirius::core::CoronaConfig c;
    c.temperature_keV = 100.0f;
    c.optical_depth = 0.5f;

    // theta_e = 100/511 ≈ 0.1957
    // tau_eff = max(0.5, 0.25) = 0.5
    // y = 4 * 0.1957 * 0.5 = 0.3914
    float y = c.ComptonizationParameter();
    float theta_e = 100.0f / 511.0f;
    float expected = 4.0f * theta_e * 0.5f;
    EXPECT_NEAR(y, expected, kEps);
}

TEST_F(CoronaConfigTests, ComptonizationHighOpticalDepth) {
    sirius::core::CoronaConfig c;
    c.temperature_keV = 50.0f;
    c.optical_depth = 2.0f;

    // tau_eff = max(2, 4) = 4
    float y = c.ComptonizationParameter();
    float theta_e = 50.0f / 511.0f;
    float expected = 4.0f * theta_e * 4.0f;
    EXPECT_NEAR(y, expected, kEps);
}

TEST_F(CoronaConfigTests, SpectralIndexFinite) {
    sirius::core::CoronaConfig c;
    c.temperature_keV = 100.0f;
    c.optical_depth = 0.5f;

    float gamma = c.SpectralIndex();
    EXPECT_FALSE(std::isnan(gamma));
    EXPECT_FALSE(std::isinf(gamma));
    EXPECT_GT(gamma, 0.0f);
}

// =============================================================================
// corona_physics Geometry Tests
// =============================================================================

class CoronaGeometryTests : public ::testing::Test {
protected:
    sirius::core::CoronaConfig config;
    float isco = 6.0f;

    void SetUp() override {
        config.enabled = true;
        config.inner_radius_M = 0.0f; // use ISCO
        config.outer_radius_M = 20.0f;
        config.scale_height_M = 5.0f;
    }
};

TEST_F(CoronaGeometryTests, DisabledReturnsNoContainment) {
    config.enabled = false;
    EXPECT_FALSE(sirius::core::corona_physics::IsInsideCorona(10.0f, kPi / 2.0f, 0.0f, config, isco));
}

TEST_F(CoronaGeometryTests, OutsideRadialBoundsRejected) {
    // Too close
    EXPECT_FALSE(sirius::core::corona_physics::IsInsideCorona(3.0f, kPi / 2.0f, 0.0f, config, isco));
    // Too far
    EXPECT_FALSE(sirius::core::corona_physics::IsInsideCorona(30.0f, kPi / 2.0f, 0.0f, config, isco));
}

TEST_F(CoronaGeometryTests, SlabGeometryEquatorialInside) {
    config.geometry = sirius::core::CoronaGeometry::Slab;
    // Equatorial plane (theta = pi/2, z = 0)
    EXPECT_TRUE(sirius::core::corona_physics::IsInsideCorona(10.0f, kPi / 2.0f, 0.0f, config, isco));
}

TEST_F(CoronaGeometryTests, SlabGeometryPolarOutside) {
    config.geometry = sirius::core::CoronaGeometry::Slab;
    // Near pole (theta ≈ 0, z ≈ r >> scale_height)
    EXPECT_FALSE(sirius::core::corona_physics::IsInsideCorona(10.0f, 0.1f, 0.0f, config, isco));
}

TEST_F(CoronaGeometryTests, LamppostNearSourceInside) {
    config.geometry = sirius::core::CoronaGeometry::Lamppost;
    config.lamppost_height_M = 10.0f;
    config.scale_height_M = 3.0f;
    // Point near lamppost (r=10, theta≈0 → z≈10, rho≈0)
    EXPECT_TRUE(sirius::core::corona_physics::IsInsideCorona(10.0f, 0.05f, 0.0f, config, isco));
}

TEST_F(CoronaGeometryTests, SphereContainment) {
    config.geometry = sirius::core::CoronaGeometry::Sphere;
    // Inside sphere
    EXPECT_TRUE(sirius::core::corona_physics::IsInsideCorona(10.0f, kPi / 4.0f, 0.0f, config, isco));
    // Outside sphere
    EXPECT_FALSE(sirius::core::corona_physics::IsInsideCorona(25.0f, kPi / 4.0f, 0.0f, config, isco));
}

TEST_F(CoronaGeometryTests, ExtendedScalesWithRadius) {
    config.geometry = sirius::core::CoronaGeometry::Extended;
    // Near equator at moderate radius: should be inside
    EXPECT_TRUE(sirius::core::corona_physics::IsInsideCorona(10.0f, kPi / 2.0f, 0.0f, config, isco));
}

// =============================================================================
// Emissivity and Optical Depth Tests
// =============================================================================

class CoronaEmissivityTests : public ::testing::Test {
protected:
    sirius::core::CoronaConfig config;
    float isco = 6.0f;

    void SetUp() override {
        config.enabled = true;
        config.inner_radius_M = 0.0f;
        config.outer_radius_M = 20.0f;
        config.scale_height_M = 5.0f;
        config.emissivity_index = 3.0f;
        config.intensity_scale = 1.0f;
        config.geometry = sirius::core::CoronaGeometry::Extended;
    }
};

TEST_F(CoronaEmissivityTests, ZeroOutsideBounds) {
    EXPECT_FLOAT_EQ(sirius::core::corona_physics::Emissivity(3.0f, kPi / 2.0f, config, isco), 0.0f);
    EXPECT_FLOAT_EQ(sirius::core::corona_physics::Emissivity(30.0f, kPi / 2.0f, config, isco), 0.0f);
}

TEST_F(CoronaEmissivityTests, DisabledReturnsZero) {
    config.enabled = false;
    EXPECT_FLOAT_EQ(sirius::core::corona_physics::Emissivity(10.0f, kPi / 2.0f, config, isco), 0.0f);
}

TEST_F(CoronaEmissivityTests, DecreasesWithRadius) {
    float e1 = sirius::core::corona_physics::Emissivity(7.0f, kPi / 2.0f, config, isco);
    float e2 = sirius::core::corona_physics::Emissivity(15.0f, kPi / 2.0f, config, isco);
    EXPECT_GT(e1, e2) << "Emissivity should decrease with radius (power-law)";
}

TEST_F(CoronaEmissivityTests, GaussianVerticalFalloff) {
    // Equatorial (z=0) should have higher emissivity than off-plane
    float equatorial = sirius::core::corona_physics::Emissivity(10.0f, kPi / 2.0f, config, isco);
    float offplane = sirius::core::corona_physics::Emissivity(10.0f, kPi / 3.0f, config, isco);
    EXPECT_GE(equatorial, offplane);
}

TEST_F(CoronaEmissivityTests, OpticalDepthDisabledIsZero) {
    config.enabled = false;
    float tau = sirius::core::corona_physics::OpticalDepthAlongRay(
        7.0f, kPi / 2.0f, 0.0f, 15.0f, kPi / 2.0f, 0.0f, config, isco);
    EXPECT_FLOAT_EQ(tau, 0.0f);
}

TEST_F(CoronaEmissivityTests, OpticalDepthPositiveInsideCorona) {
    float tau = sirius::core::corona_physics::OpticalDepthAlongRay(
        7.0f, kPi / 2.0f, 0.0f, 15.0f, kPi / 2.0f, 0.0f, config, isco);
    EXPECT_GT(tau, 0.0f);
}

TEST_F(CoronaEmissivityTests, ScatteredIntensityZeroForZeroTau) {
    float I = sirius::core::corona_physics::ScatteredIntensity(1.0f, 0.0f, config);
    EXPECT_FLOAT_EQ(I, 0.0f);
}

TEST_F(CoronaEmissivityTests, ScatteredIntensityIncreasesWithTau) {
    float I1 = sirius::core::corona_physics::ScatteredIntensity(1.0f, 0.1f, config);
    float I2 = sirius::core::corona_physics::ScatteredIntensity(1.0f, 1.0f, config);
    EXPECT_GT(I2, I1);
}

} // namespace sirius::test
