// Environment-override tests: ApplyEnvironmentOverrides reads SIRIUS_* variables
// and writes them onto the configuration, leaving unset fields at their prior
// value. Each test scopes its own variables so the estate stays order-independent.

#include "sirius/app/config/config_loader.h"
#include "sirius/app/config/config_schema.h"

#include <gtest/gtest.h>

#include <cstdlib>

namespace sirius::app::test {

// RAII guard that sets an environment variable and restores it on scope exit.
class ScopedEnv {
  public:
    ScopedEnv(const char* name, const char* value) : name_(name) { ::setenv(name_, value, 1); }
    ~ScopedEnv() { ::unsetenv(name_); }
    ScopedEnv(const ScopedEnv&) = delete;
    ScopedEnv& operator=(const ScopedEnv&) = delete;

  private:
    const char* name_;
};

TEST(ConfigEnvironment, IntegerOverridesApplied) {
    ScopedEnv width("SIRIUS_WIDTH", "1280");
    ScopedEnv samples("SIRIUS_SAMPLES", "32");

    SiriusConfig config = SiriusConfig::defaults();
    ConfigLoader::ApplyEnvironmentOverrides(config);

    EXPECT_EQ(config.render.width, 1280);
    EXPECT_EQ(config.render.samplesPerPixel, 32);
    // Untouched fields keep their defaults.
    EXPECT_EQ(config.render.height, 1080);
}

TEST(ConfigEnvironment, MetricNameAndSpinOverridesApplied) {
    ScopedEnv metric("SIRIUS_METRIC", "Kerr");
    ScopedEnv spin("SIRIUS_SPIN", "0.85");

    SiriusConfig config = SiriusConfig::defaults();
    ConfigLoader::ApplyEnvironmentOverrides(config);

    EXPECT_EQ(config.metric.name, "Kerr");
    EXPECT_DOUBLE_EQ(config.metric.spin, 0.85);
}

TEST(ConfigEnvironment, BooleanOverrideParsed) {
    ScopedEnv bloom("SIRIUS_BLOOM", "off");

    SiriusConfig config = SiriusConfig::defaults();
    config.postprocess.enableBloom = true;
    ConfigLoader::ApplyEnvironmentOverrides(config);

    EXPECT_FALSE(config.postprocess.enableBloom);
}

TEST(ConfigEnvironment, MalformedIntegerLeavesDefault) {
    ScopedEnv width("SIRIUS_WIDTH", "not-a-number");

    SiriusConfig config = SiriusConfig::defaults();
    ConfigLoader::ApplyEnvironmentOverrides(config);

    // A parse failure declines the override rather than corrupting the field.
    EXPECT_EQ(config.render.width, 1920);
}

}  // namespace sirius::app::test
