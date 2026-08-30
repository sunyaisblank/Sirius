// Starfield generator tests: config validation, star colour/intensity, catalog
// generation. Ported from TSSF001A.cpp; assertions and tolerances unchanged.

#include "sirius/core/starfield.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

namespace sirius::test {
using namespace sirius::core;

constexpr float kEps = 1e-5f;

// =============================================================================
// StarfieldConfig Tests
// =============================================================================

class StarfieldConfigTests : public ::testing::Test {};

TEST_F(StarfieldConfigTests, PointCatalogueProjectionHasNoUnconsumedSamplingControls) {
    PointStarfieldConfig point;
    point.star_count = 17;
    point.min_distance_pc = 2.0f;
    point.max_distance_pc = 200.0f;
    point.brightness_scale = 75.0f;
    point.seed = 9;
    ASSERT_TRUE(IsRepresentedPointStarfieldConfig(point));

    const StarfieldConfig expanded = ExpandPointStarfieldConfig(point);
    EXPECT_EQ(expanded.star_count, point.star_count);
    EXPECT_FLOAT_EQ(expanded.min_distance_pc, point.min_distance_pc);
    EXPECT_FLOAT_EQ(expanded.max_distance_pc, point.max_distance_pc);
    EXPECT_FLOAT_EQ(expanded.brightness_scale, point.brightness_scale);
    EXPECT_EQ(expanded.seed, point.seed);

    point.min_distance_pc = std::numeric_limits<float>::quiet_NaN();
    EXPECT_FALSE(IsRepresentedPointStarfieldConfig(point));
    EXPECT_DEATH((void)ExpandPointStarfieldConfig(point), "precondition.*enforced, terminating");
}

TEST_F(StarfieldConfigTests, OutOfRangeStarCountFailsClosed) {
    sirius::core::StarfieldConfig c;
    c.star_count = 20000000u;  // above 10^7 limit
    EXPECT_FALSE(IsRepresentedStarfieldConfig(c));
    EXPECT_DEATH((void)StarfieldGenerator(c), "precondition.*enforced, terminating");
}

TEST_F(StarfieldConfigTests, NonPositiveMinimumDistanceFailsClosed) {
    sirius::core::StarfieldConfig c;
    c.min_distance_pc = -5.0f;
    EXPECT_FALSE(IsRepresentedStarfieldConfig(c));
    EXPECT_DEATH((void)StarfieldGenerator(c), "precondition.*enforced, terminating");
}

TEST_F(StarfieldConfigTests, UnorderedDistanceDomainFailsClosed) {
    sirius::core::StarfieldConfig c;
    c.min_distance_pc = 100.0f;
    c.max_distance_pc = 50.0f;
    EXPECT_FALSE(IsRepresentedStarfieldConfig(c));
    EXPECT_DEATH((void)StarfieldGenerator(c), "precondition.*enforced, terminating");
}

TEST_F(StarfieldConfigTests, NarrowOrderedDistanceDomainIsPreservedExactly) {
    PointStarfieldConfig point;
    point.star_count = 1;
    point.min_distance_pc = 1.0f;
    point.max_distance_pc = 1.5f;
    ASSERT_TRUE(IsRepresentedPointStarfieldConfig(point));

    const StarfieldConfig expanded = ExpandPointStarfieldConfig(point);
    EXPECT_FLOAT_EQ(expanded.min_distance_pc, point.min_distance_pc);
    EXPECT_FLOAT_EQ(expanded.max_distance_pc, point.max_distance_pc);
    EXPECT_TRUE(IsRepresentedStarfieldConfig(expanded));
}

TEST_F(StarfieldConfigTests, OutOfRangeMagnitudeLimitFailsClosed) {
    sirius::core::StarfieldConfig c;
    c.magnitude_limit = 25.0f;
    EXPECT_FALSE(IsRepresentedStarfieldConfig(c));
    EXPECT_DEATH((void)StarfieldGenerator(c), "precondition.*enforced, terminating");
}

TEST_F(StarfieldConfigTests, OutOfRangeApertureFailsClosed) {
    sirius::core::StarfieldConfig c;
    c.aperture_mm = 5000.0f;
    EXPECT_FALSE(IsRepresentedStarfieldConfig(c));
    EXPECT_DEATH((void)StarfieldGenerator(c), "precondition.*enforced, terminating");
}

TEST_F(StarfieldConfigTests, NonFiniteScalarFailsClosedWithoutRewritingTheRequest) {
    sirius::core::StarfieldConfig c;
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const float infinity = std::numeric_limits<float>::infinity();
    c.min_distance_pc = nan;
    c.max_distance_pc = infinity;
    c.magnitude_limit = nan;
    c.aperture_mm = infinity;
    c.focus_distance_pc = nan;
    c.brightness_scale = infinity;
    EXPECT_FALSE(IsRepresentedStarfieldConfig(c));
    EXPECT_DEATH((void)StarfieldGenerator(c), "precondition.*enforced, terminating");
    EXPECT_TRUE(std::isnan(c.min_distance_pc));
    EXPECT_TRUE(std::isinf(c.max_distance_pc));
}

// =============================================================================
// StarEntry Tests
// =============================================================================

class StarEntryTests : public ::testing::Test {};

TEST_F(StarEntryTests, ComputeColorProducesValidRGB) {
    sirius::core::StarEntry star{};
    star.temperature_K = 5778.0f;  // Solar temperature
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

TEST_F(StarEntryTests, ComputeColorUsesTheIntegratedBlackbodyAuthority) {
    sirius::core::StarEntry star{};
    star.temperature_K = 5778.0f;
    float r = 0.0f;
    float g = 0.0f;
    float b = 0.0f;
    star.ComputeColor(r, g, b);
    const auto expected = sirius::core::spectral::BlackbodyToRgb(star.temperature_K);
    EXPECT_FLOAT_EQ(r, expected.r);
    EXPECT_FLOAT_EQ(g, expected.g);
    EXPECT_FLOAT_EQ(b, expected.b);
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

TEST_F(StarEntryTests, InvalidTemperatureFailsClosedInsteadOfDefaultingToSolar) {
    sirius::core::StarEntry invalid{};
    invalid.temperature_K = 0.0f;
    float invalid_r, invalid_g, invalid_b;
    EXPECT_DEATH(invalid.ComputeColor(invalid_r, invalid_g, invalid_b),
                 "precondition.*enforced, terminating");

    invalid.temperature_K = std::numeric_limits<float>::quiet_NaN();
    EXPECT_DEATH(invalid.ComputeColor(invalid_r, invalid_g, invalid_b),
                 "precondition.*enforced, terminating");
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

    s1.magnitude = -100.0f;
    EXPECT_DEATH((void)s1.Intensity(), "precondition.*enforced, terminating");
    s1.magnitude = std::numeric_limits<float>::quiet_NaN();
    EXPECT_DEATH((void)s1.Intensity(), "precondition.*enforced, terminating");
}

// =============================================================================
// StarfieldGenerator Tests
// =============================================================================

class StarfieldGeneratorTests : public ::testing::Test {
  protected:
    sirius::core::StarfieldConfig config;

    void SetUp() override {
        config.star_count = 1000;
        config.magnitude_limit = 20.0f;  // accept all
        config.seed = 42;
    }
};

TEST_F(StarfieldGeneratorTests, GeneratesNonEmptyCatalog) {
    sirius::core::StarfieldGenerator gen(config);
    auto stars = gen.Generate();
    EXPECT_GT(stars.size(), 0u);
}

TEST_F(StarfieldGeneratorTests, SpatialIndexOwnsValidatedCatalogueSnapshot) {
    config.star_count = 5;
    sirius::core::StarfieldGenerator gen(config);
    auto stars = gen.GenerateCatalogue();
    ASSERT_EQ(stars.size(), 5u);
    ASSERT_TRUE(std::all_of(stars.begin(), stars.end(), IsRepresentedStarEntry));
    const StarfieldSpatialIndex index(stars);
    ASSERT_EQ(index.Size(), stars.size());

    float exhaustive_r = 0.0f;
    float exhaustive_g = 0.0f;
    float exhaustive_b = 0.0f;
    float indexed_r = 0.0f;
    float indexed_g = 0.0f;
    float indexed_b = 0.0f;
    gen.AccumulateThroughBeam(1.0f, 0.0f, 0.0f, 0.5f, stars, exhaustive_r, exhaustive_g,
                              exhaustive_b);
    gen.AccumulateThroughBeam(1.0f, 0.0f, 0.0f, 0.5f, index, indexed_r, indexed_g, indexed_b);
    const float scale =
        std::max({std::abs(exhaustive_r), std::abs(exhaustive_g), std::abs(exhaustive_b), 1.0f});
    EXPECT_NEAR(indexed_r, exhaustive_r, 1.0e-6f * scale);
    EXPECT_NEAR(indexed_g, exhaustive_g, 1.0e-6f * scale);
    EXPECT_NEAR(indexed_b, exhaustive_b, 1.0e-6f * scale);

    const float indexed_direction_x = index.Stars()[0].direction_x;
    stars[0].direction_x = 0.0f;
    EXPECT_FLOAT_EQ(index.Stars()[0].direction_x, indexed_direction_x);

    stars[0].direction_x = std::numeric_limits<float>::quiet_NaN();
    EXPECT_FALSE(IsRepresentedStarEntry(stars[0]));
    EXPECT_DEATH((void)StarfieldSpatialIndex(stars), "precondition.*enforced, terminating");
}

TEST_F(StarfieldGeneratorTests, EllipticalFilterUsesTheBeamSachsBasis) {
    // This direction's least-aligned coordinate axis is y. The former point
    // filter projected z here, so it interpreted a beam angle in a different
    // tangent basis and rotated an anisotropic footprint.
    std::array<float, 3> direction{0.8f, 0.1f, 0.59f};
    float norm = std::sqrt(direction[0] * direction[0] + direction[1] * direction[1] +
                           direction[2] * direction[2]);
    for (float& component : direction) component /= norm;
    const auto basis = relativity::MakeCelestialTangentBasis(direction);
    ASSERT_TRUE(basis.has_value());
    EXPECT_GT(std::abs(basis->first[1]), 0.9f)
        << "fixture no longer distinguishes the old z-projection basis";

    constexpr float angular_offset = 0.02f;
    StarEntry star{};
    for (std::size_t component = 0; component < direction.size(); ++component) {
        const float value = std::cos(angular_offset) * direction[component] +
                            std::sin(angular_offset) * basis->first[component];
        if (component == 0) star.direction_x = value;
        if (component == 1) star.direction_y = value;
        if (component == 2) star.direction_z = value;
    }
    star.distance_pc = 10.0f;
    star.magnitude = 0.0f;
    star.color_bv = 0.65f;
    star.temperature_K = 5778.0f;
    ASSERT_TRUE(IsRepresentedStarEntry(star));

    config.star_count = 1;
    config.brightness_scale = 1.0f;
    StarfieldGenerator generator(config);
    const StarfieldSpatialIndex index(std::vector<StarEntry>{star});
    constexpr float sigma_major = 0.05f;
    constexpr float sigma_minor = 0.005f;

    float color_r = 0.0f;
    float color_g = 0.0f;
    float color_b = 0.0f;
    star.ComputeColor(color_r, color_g, color_b);
    const float expected_major_weight =
        std::exp(-0.5f * angular_offset * angular_offset / (sigma_major * sigma_major));
    const float expected_minor_weight =
        std::exp(-0.5f * angular_offset * angular_offset / (sigma_minor * sigma_minor));

    float major_r = 0.0f;
    float major_g = 0.0f;
    float major_b = 0.0f;
    generator.AccumulateThroughBeam(direction[0], direction[1], direction[2], sigma_major,
                                    sigma_minor, 0.0f, index, major_r, major_g, major_b);
    EXPECT_NEAR(major_r, color_r * expected_major_weight, 2.0e-5f);
    EXPECT_NEAR(major_g, color_g * expected_major_weight, 2.0e-5f);
    EXPECT_NEAR(major_b, color_b * expected_major_weight, 2.0e-5f);

    float minor_r = 0.0f;
    float minor_g = 0.0f;
    float minor_b = 0.0f;
    generator.AccumulateThroughBeam(direction[0], direction[1], direction[2], sigma_major,
                                    sigma_minor, static_cast<float>(std::numbers::pi / 2.0), index,
                                    minor_r, minor_g, minor_b);
    EXPECT_NEAR(minor_r, color_r * expected_minor_weight, 2.0e-5f);
    EXPECT_NEAR(minor_g, color_g * expected_minor_weight, 2.0e-5f);
    EXPECT_NEAR(minor_b, color_b * expected_minor_weight, 2.0e-5f);
    EXPECT_GT(major_r + major_g + major_b, 100.0f * (minor_r + minor_g + minor_b));
}

TEST_F(StarfieldGeneratorTests, CatalogSizeBounded) {
    sirius::core::StarfieldGenerator gen(config);
    auto stars = gen.Generate();
    ASSERT_FALSE(stars.empty());
    EXPECT_LE(stars.size(), static_cast<std::size_t>(config.star_count));
}

TEST_F(StarfieldGeneratorTests, DirectionVectorsNormalised) {
    sirius::core::StarfieldGenerator gen(config);
    auto stars = gen.Generate();
    ASSERT_FALSE(stars.empty());
    for (std::size_t i = 0; i < std::min(stars.size(), std::size_t(100)); ++i) {
        float mag = std::sqrt(stars[i].direction_x * stars[i].direction_x +
                              stars[i].direction_y * stars[i].direction_y +
                              stars[i].direction_z * stars[i].direction_z);
        EXPECT_NEAR(mag, 1.0f, 1e-4f) << "Star " << i << " direction not unit";
    }
}

TEST_F(StarfieldGeneratorTests, AllTemperaturesPositive) {
    sirius::core::StarfieldGenerator gen(config);
    auto stars = gen.Generate();
    ASSERT_FALSE(stars.empty());
    for (const auto& s : stars) {
        EXPECT_GT(s.temperature_K, 0.0f);
    }
}

TEST_F(StarfieldGeneratorTests, AllDistancesPositive) {
    sirius::core::StarfieldGenerator gen(config);
    auto stars = gen.Generate();
    ASSERT_FALSE(stars.empty());
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
    ASSERT_FALSE(stars1.empty());
    for (std::size_t i = 0; i < stars1.size(); ++i) {
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
    std::size_t n = std::min(stars1.size(), stars2.size());
    for (std::size_t i = 0; i < n; ++i) {
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
    ASSERT_FALSE(stars.empty());
    for (const auto& s : stars) {
        EXPECT_FALSE(std::isnan(s.direction_x));
        EXPECT_FALSE(std::isnan(s.direction_y));
        EXPECT_FALSE(std::isnan(s.direction_z));
        EXPECT_FALSE(std::isnan(s.distance_pc));
        EXPECT_FALSE(std::isnan(s.magnitude));
        EXPECT_FALSE(std::isnan(s.temperature_K));
    }
}

}  // namespace sirius::test
