// Configuration validation tests: the ConfigLoader::Validate contract accepts a
// default configuration and rejects each out-of-range or unknown field with a
// message. These exercise the CLI/config objects linked into this target (July
// remediation ledger item 8).

#include "sirius/app/config/config_loader.h"
#include "sirius/app/config/config_schema.h"

#include <gtest/gtest.h>

namespace sirius::app::test {

TEST(ConfigValidation, DefaultConfigurationIsValid) {
    SiriusConfig config = SiriusConfig::defaults();
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, WidthBelowMinimumRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.render.width = 64;  // Below the 128 floor.
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, NonPowerOfTwoTileSizeRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.render.tileSize = 48;  // In range but not a power of two.
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, UnknownMetricNameRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.metric.name = "Gargantua";
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, KnownMetricAliasAccepted) {
    SiriusConfig config = SiriusConfig::defaults();
    config.metric.name = "Wormhole";  // Morris-Thorne alias in the registry.
    // The name is valid; validation does not gate on backend support (that is
    // the session's decline), so no metric-name error is produced.
    auto errors = ConfigLoader::Validate(config);
    for (const auto& e : errors) {
        EXPECT_EQ(e.find("is not a known metric"), std::string::npos) << e;
    }
}

TEST(ConfigValidation, SpinAboveNearExtremalRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.metric.name = "Kerr";
    config.metric.spin = 1.5;
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, RotatingLambdaRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.metric.spin = 0.5;
    config.metric.lambda = 0.01;  // Lambda requires spin = 0.
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, ObserverInsidePoleBufferRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.observer.inclination = 0.0;  // At the pole.
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, UnknownTonemapperRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.postprocess.tonemapper = "Cineon";
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, VulkanBackendNameAccepted) {
    // Go-live (2026-07-18): vulkan joined the accepted set; device absence is
    // the render's decline, not a config-name error.
    SiriusConfig config = SiriusConfig::defaults();
    config.backend.preferred = "vulkan";
    for (const auto& e : ConfigLoader::Validate(config)) {
        EXPECT_EQ(e.find("backend.preferred"), std::string::npos) << e;
    }
}

TEST(ConfigValidation, UnknownBackendNameRejected) {
    SiriusConfig config = SiriusConfig::defaults();
    config.backend.preferred = "optix";  // Retired; never silently remapped.
    EXPECT_FALSE(ConfigLoader::Validate(config).empty());
}

TEST(ConfigValidation, CpuBackendAccepted) {
    SiriusConfig config = SiriusConfig::defaults();
    config.backend.preferred = "cpu";
    EXPECT_TRUE(ConfigLoader::Validate(config).empty());
}

}  // namespace sirius::app::test
