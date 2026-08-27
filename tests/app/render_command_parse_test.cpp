// Render-command argument parsing tests. These are the CLI parsing tests the
// July ledger (item 8) recorded as never linked to the CLI objects; they are
// linked here. The parse path is driven through the public Execute contract: a
// trailing unrecognised option stops parsing before any render runs, leaving the
// config populated for inspection, and the GPU-request and unknown-metric paths
// return non-zero without rendering.

#include "sirius/app/cli/render_command.h"
#include "sirius/app/cli/view_command.h"
#include "sirius/app/config/config_schema.h"
#include "sirius/app/viewer/interactive_viewer.h"

#include <gtest/gtest.h>

#include "support/scoped_environment.h"

#ifdef SIRIUS_HAS_VULKAN_BACKEND
#include "sirius/backend/device.h"
#endif

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <string>
#include <thread>

namespace sirius::app::test {

// Drive Execute with a trailing sentinel option so parsing populates the config
// then declines before validation or rendering.
constexpr const char* kStopSentinel = "--stop-before-render";

TEST(RenderCommandParse, BasicFlagsMapToConfig) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();

    std::vector<std::string> args = {
        "-w",  "640", "-h", "360",   "-s", "16",    "-m",           "Kerr",         "-a",
        "0.9", "-d",  "40", "--fov", "75", "--cpu", "--color-mode", "Polarisation", kStopSentinel};

    int rc = cmd.Execute(args, globals, config);

    EXPECT_EQ(rc, 1);  // Sentinel stops the run.
    EXPECT_EQ(config.render.width, 640);
    EXPECT_EQ(config.render.height, 360);
    EXPECT_EQ(config.render.samples_per_pixel, 16);
    EXPECT_EQ(config.metric.name, "Kerr");
    EXPECT_DOUBLE_EQ(config.metric.spin, 0.9);
    EXPECT_DOUBLE_EQ(config.observer.distance, 40.0);
    EXPECT_DOUBLE_EQ(config.observer.fov, 75.0);
    EXPECT_EQ(config.backend.preferred, "cpu");
    EXPECT_EQ(config.color_mode, "Polarisation");
}

TEST(RenderCommandParse, RepresentedVolumetricAndFilmFlagsSetEnables) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();

    std::vector<std::string> args = {"--volumetric",  "--turbulence", "--corona",  "--film",
                                     "--film-preset", "Interstellar", "--no-disk", kStopSentinel};

    int rc = cmd.Execute(args, globals, config);

    EXPECT_EQ(rc, 1);
    EXPECT_TRUE(config.volumetric.enabled);
    EXPECT_TRUE(config.volumetric.enable_turbulence);
    EXPECT_TRUE(config.volumetric.enable_corona);
    EXPECT_TRUE(config.film.enabled);
    EXPECT_EQ(config.film.preset, "Interstellar");
    EXPECT_FALSE(config.disk_enabled);
}

TEST(RenderCommandParse, MotionBlurAndWormholeTopologyReachTheValidatedSchema) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();

    EXPECT_EQ(cmd.Execute({"--motion-blur", "--shutter-time", "0.25", "--motion-samples", "7",
                           "--wormhole-topology", "two-sheet", kStopSentinel},
                          globals, config),
              1);
    EXPECT_TRUE(config.motion_blur.enabled);
    EXPECT_FLOAT_EQ(config.motion_blur.shutter_time, 0.25f);
    EXPECT_EQ(config.motion_blur.samples, 7);
    EXPECT_EQ(config.metric.wormhole_topology, "TwoSheet");
}

TEST(RenderCommandParse, ExplicitGpuRequestRunsVulkanWhenDevicePresent) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();

    // --gpu is wired to the Vulkan render path: with a device present it renders
    // a small Kerr scene (exit 0); with none it declines cleanly (exit 1). It
    // never silently falls back to CPU.
    bool device_present = false;
#ifdef SIRIUS_HAS_VULKAN_BACKEND
    if (auto devices = backend::EnumerateVulkanDevices();
        devices.has_value() && !devices->empty()) {
        device_present = true;
    }
#endif
    const std::string out = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp") +
                            "/sirius_gpu_parse.ppm";
    int rc = cmd.Execute({"--gpu", "-m", "Kerr", "-w", "128", "-h", "128", "-s", "1", "-o", out},
                         globals, config);
    EXPECT_EQ(rc, device_present ? 0 : 1);
    EXPECT_EQ(config.backend.preferred, "vulkan");
}

TEST(RenderCommandParse, BackendVulkanDeclinesMetricOffTheRenderPath) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();

    // --backend vulkan routes to the Vulkan render path, which carries the
    // registry gpu_supported render set. A charge metric (Reissner-Nordstrom,
    // gpu_supported = false) declines cleanly with a non-zero exit whether or not
    // a device is present: no device declines at the CLI, a device declines in
    // the render path before any dispatch. It never falls back to CPU silently.
    int rc = cmd.Execute(
        {"--backend", "vulkan", "-m", "Reissner-Nordstrom", "--no-disk", "-w", "128", "-h", "128"},
        globals, config);
    EXPECT_EQ(rc, 1);
    EXPECT_EQ(config.backend.preferred, "vulkan");
}

TEST(RenderCommandParse, ReusedCommandDoesNotRetainAnEarlierGpuRequest) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig first = SiriusConfig::Defaults();
    SiriusConfig second = SiriusConfig::Defaults();

    EXPECT_EQ(cmd.Execute({"--gpu", kStopSentinel}, globals, first), 1);
    sirius::test::ScopedEnvironmentVariable icd("VK_ICD_FILENAMES",
                                                "/sirius/intentional/missing-vulkan-icd.json");
    sirius::test::ScopedEnvironmentVariable drivers(
        "VK_DRIVER_FILES", "/sirius/intentional/missing-vulkan-driver.json");
    sirius::test::ScopedEnvironmentVariable additional_drivers("VK_ADD_DRIVER_FILES", "");
    const auto output = std::filesystem::temp_directory_path() / "sirius_reused_render_command.ppm";
    std::filesystem::remove(output);
    EXPECT_EQ(cmd.Execute({"--metric", "Minkowski", "--no-disk", "--no-bloom", "--width", "128",
                           "--height", "128", "--samples", "1", "--output", output.string()},
                          globals, second),
              0);
    EXPECT_TRUE(std::filesystem::exists(output));
    std::filesystem::remove(output);
}

TEST(RenderCommandParse, CliCpuOverridesLowerLayerVulkanBackend) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();
    config.backend.preferred = "vulkan";

    sirius::test::ScopedEnvironmentVariable icd("VK_ICD_FILENAMES",
                                                "/sirius/intentional/missing-vulkan-icd.json");
    sirius::test::ScopedEnvironmentVariable drivers(
        "VK_DRIVER_FILES", "/sirius/intentional/missing-vulkan-driver.json");
    sirius::test::ScopedEnvironmentVariable additional_drivers("VK_ADD_DRIVER_FILES", "");

    const auto output = std::filesystem::temp_directory_path() / "sirius_cli_cpu_precedence.ppm";
    std::filesystem::remove(output);
    EXPECT_EQ(cmd.Execute({"--cpu", "--metric", "Minkowski", "--no-disk", "--no-bloom", "--width",
                           "128", "--height", "128", "--samples", "1", "--output", output.string()},
                          globals, config),
              0);
    EXPECT_EQ(config.backend.preferred, "cpu");
    EXPECT_TRUE(std::filesystem::exists(output));
    std::filesystem::remove(output);
}

TEST(RenderCommandParse, UnknownMetricFailsValidation) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();

    // Parses cleanly, then validation rejects the unknown metric name.
    int rc = cmd.Execute({"--cpu", "-m", "NotAMetric"}, globals, config);
    EXPECT_EQ(rc, 1);
    EXPECT_EQ(config.metric.name, "NotAMetric");
}

TEST(RenderCommandParse, UnknownOptionRejected) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();

    int rc = cmd.Execute({"--not-a-flag"}, globals, config);
    EXPECT_EQ(rc, 1);
}

TEST(RenderCommandParse, TrailingNumericGarbageRejected) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();
    EXPECT_EQ(cmd.Execute({"--width", "512junk"}, globals, config), 1);
}

TEST(RenderCommandParse, NonFiniteNumericValueRejected) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();
    EXPECT_EQ(cmd.Execute({"--spin", "nan"}, globals, config), 1);
}

TEST(RenderCommandParse, UnexpectedPositionalArgumentRejected) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();
    EXPECT_EQ(cmd.Execute({"unclaimed-scene-name"}, globals, config), 1);
}

TEST(RenderCommandParse, ExplicitMassOnMasslessMetricIsNotSilentlyDiscarded) {
    RenderCommand cmd;
    SiriusConfig config = SiriusConfig::Defaults();
    GlobalOptions globals;
    EXPECT_NE(cmd.Execute({"--cpu", "--metric", "Minkowski", "--mass", "1", "--no-disk", "--width",
                           "128", "--height", "128", "--samples", "1"},
                          globals, config),
              0);
}

TEST(ViewCommandOperational, StrictParsingAndSessionProjection) {
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();
    ViewCommand view;
    EXPECT_EQ(view.Execute({"--width", "512junk"}, globals, config), 1);
    EXPECT_EQ(view.Execute({"--unknown-view-option"}, globals, config), 1);

    ViewerConfig viewer_config;
    viewer_config.backend = render::RenderBackend::Vulkan;
    viewer_config.metric_id = core::MetricId::Kerr;
    viewer_config.enable_disk = false;
    viewer_config.enable_volumetric = true;
    viewer_config.enable_jets = true;
    InteractiveViewer viewer;
    ASSERT_TRUE(viewer.Initialise(viewer_config));
    viewer.SetCameraPosition(40.0, 1.0, 0.25);
    const render::SessionConfig projected = viewer.CreateSessionConfig(320, 180, 3);
    EXPECT_EQ(projected.backend, render::RenderBackend::Vulkan);
    EXPECT_EQ(projected.metric_id, core::MetricId::Kerr);
    EXPECT_DOUBLE_EQ(projected.observer_azimuth, 0.25);
    EXPECT_FALSE(projected.enable_disk);
    EXPECT_FALSE(projected.enable_volumetric_disk);
    EXPECT_TRUE(projected.enable_jets);
    EXPECT_FALSE(projected.write_output);

    ViewerConfig invalid = viewer_config;
    invalid.refinement_levels = 0;
    InteractiveViewer invalid_viewer;
    EXPECT_FALSE(invalid_viewer.Initialise(invalid));
}

TEST(ViewCommandOperational, HeadlessRefinementProducesASynchronisedFrame) {
    ViewerConfig config;
    config.preview_width = 64;
    config.preview_height = 64;
    config.final_width = 64;
    config.final_height = 64;
    config.refinement_levels = 1;
    config.samples_per_level = 1;
    config.backend = render::RenderBackend::Cpu;
    config.metric_id = core::MetricId::Schwarzschild;
    config.black_hole_spin = 0.0;

    InteractiveViewer viewer;
    ASSERT_TRUE(viewer.Initialise(config));
    ASSERT_TRUE(viewer.Start());

    // This is a liveness bound, not a performance assertion. Keep a real
    // render and complete frame publication while avoiding coupling the
    // correctness gate to hosted-runner throughput under three sanitizers.
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(60);
    while (!viewer.GetRefinementState().complete && viewer.GetLastError().empty() &&
           std::chrono::steady_clock::now() < deadline) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    const auto refinement = viewer.GetRefinementState();
    viewer.Stop();

    EXPECT_TRUE(viewer.GetLastError().empty()) << viewer.GetLastError();
    EXPECT_TRUE(refinement.complete);
    EXPECT_EQ(refinement.current_width, 64);
    EXPECT_EQ(refinement.current_height, 64);
    EXPECT_EQ(refinement.current_samples_per_pixel, 1);
    EXPECT_EQ(viewer.GetFrameBufferSnapshot().size(), 64u * 64u * 4u);
}

TEST(ViewCommandOperational, VulkanRefinementPublishesProgressiveFrames) {
#ifndef SIRIUS_HAS_VULKAN_BACKEND
    GTEST_SKIP() << "Vulkan backend was not compiled";
#else
    const auto devices = backend::EnumerateVulkanDevices();
    if (!devices || devices->empty()) {
        GTEST_SKIP() << "no Vulkan device present";
    }

    ViewerConfig config;
    config.preview_width = 64;
    config.preview_height = 64;
    config.final_width = 96;
    config.final_height = 64;
    config.refinement_levels = 2;
    config.samples_per_level = 1;
    config.backend = render::RenderBackend::Vulkan;
    config.metric_id = core::MetricId::Schwarzschild;
    config.black_hole_spin = 0.0;
    config.enable_disk = false;

    std::atomic<int> frame_count = 0;
    std::atomic<int> final_width = 0;
    std::atomic<int> final_height = 0;
    InteractiveViewer viewer;
    ASSERT_TRUE(viewer.Initialise(config));
    viewer.SetFrameCallback([&](const float* data, int width, int height) {
        if (data != nullptr) {
            final_width = width;
            final_height = height;
            ++frame_count;
        }
    });
    ASSERT_TRUE(viewer.Start());

    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(60);
    while (!viewer.GetRefinementState().complete && viewer.GetLastError().empty() &&
           std::chrono::steady_clock::now() < deadline) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    const auto refinement = viewer.GetRefinementState();
    viewer.Stop();

    EXPECT_TRUE(viewer.GetLastError().empty()) << viewer.GetLastError();
    EXPECT_TRUE(refinement.complete);
    EXPECT_EQ(refinement.current_width, 96);
    EXPECT_EQ(refinement.current_height, 64);
    EXPECT_GE(frame_count.load(), 2);
    EXPECT_EQ(final_width.load(), 96);
    EXPECT_EQ(final_height.load(), 64);
    EXPECT_EQ(viewer.GetFrameBufferSnapshot().size(), 96u * 64u * 4u);
#endif
}

TEST(ViewCommandOperational, InputStateHandlesPressRepeatReleaseMouseAndScroll) {
    ViewerConfig config;
    config.move_speed = 2.0f;
    config.mouse_sensitivity = 0.01f;
    InteractiveViewer viewer;
    ASSERT_TRUE(viewer.Initialise(config));

    const auto initial = viewer.GetCameraState();
    viewer.ProcessKey(87, 1);  // W press.
    viewer.UpdateCamera(0.5f);
    const auto pressed = viewer.GetCameraState();
    EXPECT_LT(pressed.r, initial.r);
    EXPECT_TRUE(viewer.GetRefinementState().needs_restart)
        << "camera input did not request progressive-refinement restart";

    viewer.ProcessKey(87, 2);  // W repeat must remain held.
    viewer.UpdateCamera(0.5f);
    const auto repeated = viewer.GetCameraState();
    EXPECT_LT(repeated.r, pressed.r);

    viewer.ProcessKey(87, 0);  // W release.
    viewer.UpdateCamera(0.5f);
    EXPECT_DOUBLE_EQ(viewer.GetCameraState().r, repeated.r);
    viewer.UpdateCamera(std::numeric_limits<float>::quiet_NaN());
    viewer.UpdateCamera(-1.0f);
    EXPECT_DOUBLE_EQ(viewer.GetCameraState().r, repeated.r);

    viewer.ProcessMouseMove(100.0, 100.0, false);
    viewer.ProcessMouseMove(110.0, 105.0, true);
    const auto mouse = viewer.GetCameraState();
    EXPECT_GT(mouse.phi, repeated.phi);
    EXPECT_GT(mouse.theta, repeated.theta);

    viewer.ProcessMouseMove(std::numeric_limits<double>::quiet_NaN(), 0.0, true);
    EXPECT_DOUBLE_EQ(viewer.GetCameraState().phi, mouse.phi);

    viewer.ProcessScroll(1000.0);
    EXPECT_FLOAT_EQ(viewer.GetCameraState().fov, 10.0f);
    viewer.ProcessScroll(-1000.0);
    EXPECT_FLOAT_EQ(viewer.GetCameraState().fov, 120.0f);
    viewer.ProcessScroll(std::numeric_limits<double>::infinity());
    EXPECT_FLOAT_EQ(viewer.GetCameraState().fov, 120.0f);

    viewer.SetCameraPosition(std::numeric_limits<double>::infinity(), 1.0, 0.0);
    EXPECT_DOUBLE_EQ(viewer.GetCameraState().r, repeated.r);
    viewer.SetCameraPosition(1.0, 0.0, 0.5);
    EXPECT_DOUBLE_EQ(viewer.GetCameraState().r, 5.0);
    EXPECT_DOUBLE_EQ(viewer.GetCameraState().theta, 0.1);

    ViewerConfig scaled = config;
    scaled.black_hole_mass = 10.0;
    scaled.observer_distance = 50.0;
    InteractiveViewer scaled_viewer;
    ASSERT_TRUE(scaled_viewer.Initialise(scaled));
    scaled_viewer.SetCameraPosition(1.0, 1.0, 0.0);
    EXPECT_DOUBLE_EQ(scaled_viewer.GetCameraState().r, 50.0);
    scaled_viewer.SetCameraPosition(20000.0, 1.0, 0.0);
    EXPECT_DOUBLE_EQ(scaled_viewer.GetCameraState().r, 10000.0);
}

TEST(RenderCommandParse, RetiredBackendNamesAreNotRemapped) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();
    EXPECT_EQ(cmd.Execute({"--backend", "optix"}, globals, config), 1);
    EXPECT_EQ(config.backend.preferred, "optix");
}

TEST(RenderCommandParse, MalformedCameraBetaRejected) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();
    EXPECT_EQ(cmd.Execute({"--camera-beta", "0.1,,0.2"}, globals, config), 1);
    EXPECT_EQ(cmd.Execute({"--camera-beta", "0.1,0.2,0.3,0.4"}, globals, config), 1);
}

}  // namespace sirius::app::test
