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

    SiriusConfig config = SiriusConfig::defaults();
    ConfigLoader::ApplyEnvironmentOverrides(config);

    EXPECT_EQ(config.render.width, 1280);
    EXPECT_EQ(config.render.samplesPerPixel, 32);
    // Untouched fields keep their defaults.
    EXPECT_EQ(config.render.height, 1080);
}

TEST(ConfigEnvironment, MetricNameAndSpinOverridesApplied) {
    ScopedEnvironmentVariable metric("SIRIUS_METRIC", "Kerr");
    ScopedEnvironmentVariable spin("SIRIUS_SPIN", "0.85");

    SiriusConfig config = SiriusConfig::defaults();
    ConfigLoader::ApplyEnvironmentOverrides(config);

    EXPECT_EQ(config.metric.name, "Kerr");
    EXPECT_DOUBLE_EQ(config.metric.spin, 0.85);
}

TEST(ConfigEnvironment, BooleanOverrideParsed) {
    ScopedEnvironmentVariable bloom("SIRIUS_BLOOM", "off");

    SiriusConfig config = SiriusConfig::defaults();
    config.postprocess.enableBloom = true;
    ConfigLoader::ApplyEnvironmentOverrides(config);

    EXPECT_FALSE(config.postprocess.enableBloom);
}

TEST(ConfigEnvironment, ColorModeOverrideApplied) {
    ScopedEnvironmentVariable color_mode("SIRIUS_COLOR_MODE", "Polarisation");
    SiriusConfig config = SiriusConfig::defaults();
    ConfigLoader::ApplyEnvironmentOverrides(config);
    EXPECT_EQ(config.colorMode, "Polarisation");
}

TEST(ConfigEnvironment, MalformedIntegerDeclines) {
    ScopedEnvironmentVariable width("SIRIUS_WIDTH", "not-a-number");

    SiriusConfig config = SiriusConfig::defaults();
    EXPECT_THROW(ConfigLoader::ApplyEnvironmentOverrides(config), std::invalid_argument);
}

TEST(ConfigEnvironment, TrailingNumericGarbageDeclines) {
    ScopedEnvironmentVariable spin("SIRIUS_SPIN", "0.5junk");
    SiriusConfig config = SiriusConfig::defaults();
    EXPECT_THROW(ConfigLoader::ApplyEnvironmentOverrides(config), std::invalid_argument);
}

TEST(ConfigEnvironment, NonFiniteNumberDeclines) {
    ScopedEnvironmentVariable spin("SIRIUS_SPIN", "nan");
    SiriusConfig config = SiriusConfig::defaults();
    EXPECT_THROW(ConfigLoader::ApplyEnvironmentOverrides(config), std::invalid_argument);
}

TEST(ConfigEnvironment, MalformedBooleanDeclines) {
    ScopedEnvironmentVariable bloom("SIRIUS_BLOOM", "sometimes");
    SiriusConfig config = SiriusConfig::defaults();
    EXPECT_THROW(ConfigLoader::ApplyEnvironmentOverrides(config), std::invalid_argument);
}

}  // namespace sirius::app::test
