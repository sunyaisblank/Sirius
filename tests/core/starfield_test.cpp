// Starfield generator tests: config validation, star colour/intensity, catalog
// generation. Ported from TSSF001A.cpp; assertions and tolerances unchanged.

#include "sirius/core/starfield.h"

#include <gtest/gtest.h>

#include <cmath>

namespace sirius::test {
using namespace sirius::core;

constexpr float kEps = 1e-5f;

// =============================================================================
// StarfieldConfig Tests
// =============================================================================

class StarfieldConfigTests : public ::testing::Test {};

TEST_F(StarfieldConfigTests, ValidateClampsStarCount) {
    sirius::core::StarfieldConfig c;
    c.star_count = 20000000u; // above 10^7 limit
    c.Validate();
    EXPECT_LE(c.star_count, 10000000u);
}

TEST_F(StarfieldConfigTests, ValidateClampsMinDistance) {
    sirius::core::StarfieldConfig c;
    c.min_distance_pc = -5.0f;
    c.Validate();
    EXPECT_GT(c.min_distance_pc, 0.0f);
}

TEST_F(StarfieldConfigTests, ValidateEnsuresMaxGreaterThanMin) {
    sirius::core::StarfieldConfig c;
    c.min_distance_pc = 100.0f;
    c.max_distance_pc = 50.0f;
    c.Validate();
    EXPECT_GT(c.max_distance_pc, c.min_distance_pc);
}

TEST_F(StarfieldConfigTests, ValidateClampsMagnitudeLimit) {
    sirius::core::StarfieldConfig c;
    c.magnitude_limit = 25.0f;
    c.Validate();
    EXPECT_LE(c.magnitude_limit, 20.0f);
}

TEST_F(StarfieldConfigTests, ValidateClampsAperture) {
    sirius::core::StarfieldConfig c;
    c.aperture_mm = 5000.0f;
    c.Validate();
    EXPECT_LE(c.aperture_mm, 1000.0f);
}

// =============================================================================
// StarEntry Tests
// =============================================================================

class StarEntryTests : public ::testing::Test {};

TEST_F(StarEntryTests, ComputeColorProducesValidRGB) {
    sirius::core::StarEntry star{};
    star.temperature_K = 5778.0f; // Solar temperature
    float r, g, b;
    star.ComputeColor(r, g, b);

    EXPECT_GE(r, 0.0f);
    EXPECT_LE(r, 1.0f);
    EXPECT_GE(g, 0.0f);
    EXPECT_LE(g, 1.0f);
    EXPECT_GE(b, 0.0f);
    EXPECT_LE(b, 1.0f);
    EXPECT_FALSE(std::isnan(r));
    EXPECT_FALSE(std::isnan(g));
    EXPECT_FALSE(std::isnan(b));
}

TEST_F(StarEntryTests, HotStarIsBluer) {
    sirius::core::StarEntry hot{};
    hot.temperature_K = 30000.0f;
    float rh, gh, bh;
    hot.ComputeColor(rh, gh, bh);

    sirius::core::StarEntry cool{};
    cool.temperature_K = 3000.0f;
    float rc, gc, bc;
    cool.ComputeColor(rc, gc, bc);

    // Hot star: blue channel dominant (b = 1.0)
    EXPECT_FLOAT_EQ(bh, 1.0f);
    // Cool star: red channel dominant (r = 1.0)
    EXPECT_FLOAT_EQ(rc, 1.0f);
}

TEST_F(StarEntryTests, ZeroTemperatureDefaultsToSolar) {
    sirius::core::StarEntry star{};
    star.temperature_K = 0.0f;
    float r, g, b;
    star.ComputeColor(r, g, b);
    EXPECT_FALSE(std::isnan(r));
}

TEST_F(StarEntryTests, IntensityFromMagnitude) {
    sirius::core::StarEntry star{};
    star.magnitude = 0.0f;
    EXPECT_NEAR(star.Intensity(), 1.0f, kEps);

    // Brighter star has higher intensity (lower magnitude)
    star.magnitude = -1.0f;
    float bright = star.Intensity();
    star.magnitude = 1.0f;
    float dim = star.Intensity();
    EXPECT_GT(bright, dim);
}

TEST_F(StarEntryTests, IntensityMagnitudeRelation) {
    // m2 - m1 = 5 → flux ratio = 100
    sirius::core::StarEntry s1{};
    s1.magnitude = 0.0f;
    sirius::core::StarEntry s2{};
    s2.magnitude = 5.0f;
    float ratio = s1.Intensity() / s2.Intensity();
    EXPECT_NEAR(ratio, 100.0f, 0.1f);
}

// =============================================================================
// StarfieldGenerator Tests
// =============================================================================

class StarfieldGeneratorTests : public ::testing::Test {
protected:
    sirius::core::StarfieldConfig config;

    void SetUp() override {
        config.star_count = 1000;
        config.magnitude_limit = 20.0f; // accept all
        config.seed = 42;
    }
};

TEST_F(StarfieldGeneratorTests, GeneratesNonEmptyCatalog) {
    sirius::core::StarfieldGenerator gen(config);
    auto stars = gen.Generate();
    EXPECT_GT(stars.size(), 0u);
}

TEST_F(StarfieldGeneratorTests, CatalogSizeBounded) {
    sirius::core::StarfieldGenerator gen(config);
    auto stars = gen.Generate();
    EXPECT_LE(stars.size(), static_cast<size_t>(config.star_count));
}

TEST_F(StarfieldGeneratorTests, DirectionVectorsNormalised) {
    sirius::core::StarfieldGenerator gen(config);
    auto stars = gen.Generate();
    for (size_t i = 0; i < std::min(stars.size(), size_t(100)); ++i) {
        float mag = std::sqrt(
            stars[i].direction_x * stars[i].direction_x +
            stars[i].direction_y * stars[i].direction_y +
            stars[i].direction_z * stars[i].direction_z
        );
        EXPECT_NEAR(mag, 1.0f, 1e-4f) << "Star " << i << " direction not unit";
    }
}

TEST_F(StarfieldGeneratorTests, AllTemperaturesPositive) {
    sirius::core::StarfieldGenerator gen(config);
    auto stars = gen.Generate();
    for (const auto& s : stars) {
        EXPECT_GT(s.temperature_K, 0.0f);
    }
}

TEST_F(StarfieldGeneratorTests, AllDistancesPositive) {
    sirius::core::StarfieldGenerator gen(config);
    auto stars = gen.Generate();
    for (const auto& s : stars) {
        EXPECT_GT(s.distance_pc, 0.0f);
    }
}

TEST_F(StarfieldGeneratorTests, DeterministicWithSameSeed) {
    sirius::core::StarfieldGenerator gen1(config);
    auto stars1 = gen1.Generate();

    sirius::core::StarfieldGenerator gen2(config);
    auto stars2 = gen2.Generate();

    ASSERT_EQ(stars1.size(), stars2.size());
    for (size_t i = 0; i < stars1.size(); ++i) {
        EXPECT_FLOAT_EQ(stars1[i].direction_x, stars2[i].direction_x);
        EXPECT_FLOAT_EQ(stars1[i].magnitude, stars2[i].magnitude);
    }
}

TEST_F(StarfieldGeneratorTests, DifferentSeedsDifferentCatalogs) {
    config.seed = 42;
    sirius::core::StarfieldGenerator gen1(config);
    auto stars1 = gen1.Generate();

    config.seed = 99;
    sirius::core::StarfieldGenerator gen2(config);
    auto stars2 = gen2.Generate();

    // At least some stars should differ
    bool any_different = false;
    size_t n = std::min(stars1.size(), stars2.size());
    for (size_t i = 0; i < n; ++i) {
        if (stars1[i].direction_x != stars2[i].direction_x) {
            any_different = true;
            break;
        }
    }
    EXPECT_TRUE(any_different);
}

TEST_F(StarfieldGeneratorTests, NoNaNInCatalog) {
    sirius::core::StarfieldGenerator gen(config);
    auto stars = gen.Generate();
    for (const auto& s : stars) {
        EXPECT_FALSE(std::isnan(s.direction_x));
        EXPECT_FALSE(std::isnan(s.direction_y));
        EXPECT_FALSE(std::isnan(s.direction_z));
        EXPECT_FALSE(std::isnan(s.distance_pc));
        EXPECT_FALSE(std::isnan(s.magnitude));
        EXPECT_FALSE(std::isnan(s.temperature_K));
    }
}

} // namespace sirius::test
