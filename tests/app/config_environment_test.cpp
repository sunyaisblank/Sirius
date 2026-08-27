// Environment-override tests: ApplyEnvironmentOverrides reads SIRIUS_* variables
// and writes them onto the configuration, leaving unset fields at their prior
// value. Each test scopes its own variables so the estate stays order-independent.

#include "sirius/app/config/config_loader.h"
#include "sirius/app/config/config_schema.h"

#include <gtest/gtest.h>

#include "support/scoped_environment.h"

namespace sirius::app::test {

using sirius::test::ScopedEnvironmentVariable;

TEST(ConfigEnvironment, IntegerOverridesApplied) {
    ScopedEnvironmentVariable width("SIRIUS_WIDTH", "1280");
    ScopedEnvironmentVariable samples("SIRIUS_SAMPLES", "32");

    SiriusConfig config = SiriusConfig::Defaults();
    ASSERT_TRUE(ConfigLoader::ApplyEnvironmentOverrides(config).has_value());

    EXPECT_EQ(config.render.width, 1280);
    EXPECT_EQ(config.render.samples_per_pixel, 32);
    // Untouched fields keep their defaults.
    EXPECT_EQ(config.render.height, 1080);
}

TEST(ConfigEnvironment, MetricNameAndSpinOverridesApplied) {
    ScopedEnvironmentVariable metric("SIRIUS_METRIC", "Kerr");
    ScopedEnvironmentVariable spin("SIRIUS_SPIN", "0.85");

    SiriusConfig config = SiriusConfig::Defaults();
    ASSERT_TRUE(ConfigLoader::ApplyEnvironmentOverrides(config).has_value());

    EXPECT_EQ(config.metric.name, "Kerr");
    EXPECT_DOUBLE_EQ(config.metric.spin, 0.85);
}

TEST(ConfigEnvironment, BooleanOverrideParsed) {
    ScopedEnvironmentVariable bloom("SIRIUS_BLOOM", "off");

    SiriusConfig config = SiriusConfig::Defaults();
    config.postprocess.enable_bloom = true;
    ASSERT_TRUE(ConfigLoader::ApplyEnvironmentOverrides(config).has_value());

    EXPECT_FALSE(config.postprocess.enable_bloom);
}

TEST(ConfigEnvironment, ColorModeOverrideApplied) {
    ScopedEnvironmentVariable color_mode("SIRIUS_COLOR_MODE", "Polarisation");
    SiriusConfig config = SiriusConfig::Defaults();
    ASSERT_TRUE(ConfigLoader::ApplyEnvironmentOverrides(config).has_value());
    EXPECT_EQ(config.color_mode, "Polarisation");
}

TEST(ConfigEnvironment, MalformedIntegerDeclines) {
    ScopedEnvironmentVariable width("SIRIUS_WIDTH", "not-a-number");

    SiriusConfig config = SiriusConfig::Defaults();
    EXPECT_FALSE(ConfigLoader::ApplyEnvironmentOverrides(config).has_value());
}

TEST(ConfigEnvironment, TrailingNumericGarbageDeclines) {
    ScopedEnvironmentVariable spin("SIRIUS_SPIN", "0.5junk");
    SiriusConfig config = SiriusConfig::Defaults();
    EXPECT_FALSE(ConfigLoader::ApplyEnvironmentOverrides(config).has_value());
}

TEST(ConfigEnvironment, NonFiniteNumberDeclines) {
    ScopedEnvironmentVariable spin("SIRIUS_SPIN", "nan");
    SiriusConfig config = SiriusConfig::Defaults();
    EXPECT_FALSE(ConfigLoader::ApplyEnvironmentOverrides(config).has_value());
}

TEST(ConfigEnvironment, MalformedBooleanDeclines) {
    ScopedEnvironmentVariable bloom("SIRIUS_BLOOM", "sometimes");
    SiriusConfig config = SiriusConfig::Defaults();
    EXPECT_FALSE(ConfigLoader::ApplyEnvironmentOverrides(config).has_value());
}

TEST(ConfigEnvironment, FailedOverrideLeavesConfigurationUnchanged) {
    ScopedEnvironmentVariable width("SIRIUS_WIDTH", "1280");
    ScopedEnvironmentVariable bloom("SIRIUS_BLOOM", "sometimes");
    SiriusConfig config = SiriusConfig::Defaults();
    const SiriusConfig original = config;

    EXPECT_FALSE(ConfigLoader::ApplyEnvironmentOverrides(config).has_value());
    EXPECT_EQ(config.render.width, original.render.width);
    EXPECT_EQ(config.postprocess.enable_bloom, original.postprocess.enable_bloom);
}

}  // namespace sirius::app::test
