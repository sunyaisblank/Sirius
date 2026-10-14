// Render-command argument parsing tests. These are the CLI parsing tests the
// July ledger (item 8) recorded as never linked to the CLI objects; they are
// linked here. The parse path is driven through the public Execute contract: a
// trailing unrecognised option stops parsing before any render runs, leaving the
// config populated for inspection, and the GPU-request and unknown-metric paths
// return non-zero without rendering.

#include <gtest/gtest.h>

#include "sirius/app/cli/render_command.h"
#include "sirius/app/config/config_schema.h"

namespace sirius::app::test {

// Drive Execute with a trailing sentinel option so parsing populates the config
// then declines before validation or rendering.
constexpr const char* kStopSentinel = "--stop-before-render";

TEST(RenderCommandParse, BasicFlagsMapToConfig) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::defaults();

    std::vector<std::string> args = {"-w",   "640", "-h",   "360",  "-s", "16",
                                     "-m",   "Kerr", "-a",  "0.9",  "-d", "40",
                                     "--fov", "75",  "--cpu", kStopSentinel};

    int rc = cmd.Execute(args, globals, config);

    EXPECT_EQ(rc, 1);  // Sentinel stops the run.
    EXPECT_EQ(config.render.width, 640);
    EXPECT_EQ(config.render.height, 360);
    EXPECT_EQ(config.render.samplesPerPixel, 16);
    EXPECT_EQ(config.metric.name, "Kerr");
    EXPECT_DOUBLE_EQ(config.metric.spin, 0.9);
    EXPECT_DOUBLE_EQ(config.observer.distance, 40.0);
    EXPECT_DOUBLE_EQ(config.observer.fov, 75.0);
    EXPECT_EQ(config.backend.preferred, "cpu");
}

TEST(RenderCommandParse, VolumetricAndFilmFlagsSetEnables) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::defaults();

    std::vector<std::string> args = {"--volumetric", "--turbulence", "--film",
                                     "--film-preset", "Interstellar", kStopSentinel};

    int rc = cmd.Execute(args, globals, config);

    EXPECT_EQ(rc, 1);
    EXPECT_TRUE(config.volumetric.enabled);
    EXPECT_TRUE(config.volumetric.enableTurbulence);
    EXPECT_TRUE(config.film.enabled);
    EXPECT_EQ(config.film.preset, "Interstellar");
}

TEST(RenderCommandParse, ExplicitGpuRequestDeclinesNonZero) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::defaults();

    // --gpu selects a GPU backend; the Vulkan path is not yet wired, so the
    // command declines cleanly with a non-zero exit rather than rendering on CPU.
    int rc = cmd.Execute({"--gpu"}, globals, config);
    EXPECT_EQ(rc, 1);
}

TEST(RenderCommandParse, BackendVulkanNameDeclinesNonZero) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::defaults();

    int rc = cmd.Execute({"--backend", "vulkan"}, globals, config);
    EXPECT_EQ(rc, 1);
}

TEST(RenderCommandParse, UnknownMetricFailsValidation) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::defaults();

    // Parses cleanly, then validation rejects the unknown metric name.
    int rc = cmd.Execute({"--cpu", "-m", "NotAMetric"}, globals, config);
    EXPECT_EQ(rc, 1);
    EXPECT_EQ(config.metric.name, "NotAMetric");
}

TEST(RenderCommandParse, UnknownOptionRejected) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::defaults();

    int rc = cmd.Execute({"--not-a-flag"}, globals, config);
    EXPECT_EQ(rc, 1);
}

}  // namespace sirius::app::test
