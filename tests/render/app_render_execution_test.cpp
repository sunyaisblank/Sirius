// Render-capable CLI and progressive-viewer execution tests. These cases live
// in the rendering-only executable by construction: the application suite may
// validate parsing, projection, and fail-closed declines, but it cannot dispatch
// CPU/Vulkan work or publish a rendered frame.

#include "sirius/app/cli/render_command.h"
#include "sirius/app/config/session_config_adapter.h"
#include "sirius/app/viewer/interactive_viewer.h"

#include <gtest/gtest.h>

#include "support/scoped_environment.h"

#ifdef SIRIUS_HAS_VULKAN_BACKEND
#include "sirius/backend/device.h"
#endif

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <future>
#include <thread>
#include <vector>

namespace sirius::app::test {

constexpr const char* kStopSentinel = "--stop-before-render";

TEST(RenderCommandParse, ExplicitGpuRequestRunsVulkanWhenDevicePresent) {
    sirius::test::ScopedEnvironmentVariable precision("SIRIUS_PRECISION", "fp32");
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();

    // --gpu is wired to Vulkan: the selected device must admit retained Kerr
    // arithmetic before it can render. Enumeration alone cannot turn a factory
    // refusal into success, and an explicit request never falls back to CPU.
    bool admitted = false;
#ifdef SIRIUS_HAS_VULKAN_BACKEND
    if (auto devices = backend::EnumerateVulkanDevices();
        devices.has_value() && !devices->empty()) {
        const auto selected = backend::ResolveVulkanDeviceIndex(*devices);
        ASSERT_TRUE(selected) << selected.error().Description();
#ifdef SIRIUS_HAS_RETAINED_COMPUTE
        admitted = true;  // Native controls or the embedded integer arithmetic.
#endif
    }
#endif
#ifdef SIRIUS_TEST_REQUIRE_VULKAN_RUNTIME
    ASSERT_TRUE(admitted) << "the required runtime profile cannot admit retained Kerr rendering";
#endif
    const auto output = std::filesystem::temp_directory_path() / "sirius_gpu_parse.ppm";
    std::filesystem::remove(output);
    const int rc = cmd.Execute(
        {"--gpu", "-m", "Kerr", "-w", "128", "-h", "128", "-s", "1", "-o", output.string()},
        globals, config);
    EXPECT_EQ(rc, admitted ? 0 : 1);
    EXPECT_EQ(config.backend.preferred, "vulkan");
    if (admitted) {
        ASSERT_TRUE(std::filesystem::exists(output));
        EXPECT_GT(std::filesystem::file_size(output), 1024u);
    } else {
        EXPECT_FALSE(std::filesystem::exists(output));
    }
    std::filesystem::remove(output);
}

TEST(RenderSessionProbe, BackendAutoPublishesCpuWithNoIcd) {
    sirius::test::ScopedEnvironmentVariable icd("VK_ICD_FILENAMES", "/sirius/missing-icd.json");
    sirius::test::ScopedEnvironmentVariable drivers("VK_DRIVER_FILES", "/sirius/missing-icd.json");
    sirius::test::ScopedEnvironmentVariable additional("VK_ADD_DRIVER_FILES", "");
    sirius::test::ScopedEnvironmentVariable selector("SIRIUS_VULKAN_DEVICE", nullptr);
    sirius::test::ScopedEnvironmentVariable precision("SIRIUS_PRECISION", "fp32");
    auto request = SiriusConfig::Defaults();
    request.metric.name = "Minkowski";
    request.metric.mass = 0;
    request.disk_enabled = false;
    request.render.samples_per_pixel = 1;
    const auto automatic = MakeSessionConfig(request);
    ASSERT_TRUE(automatic) << automatic.error().Description();
    ASSERT_EQ(automatic->backend, render::RenderBackend::Cpu);
    // Operator dimensions retain their schema bounds. This typed in-memory
    // session narrows only the publication witness to four straight flat rays.
    auto config = *automatic;
    config.width = config.height = 2;
    config.write_output = false;
    config.enable_bloom = false;
    render::RenderSession session;
    ASSERT_TRUE(session.Configure(config));
    unsigned publications = 0;
    session.SetProgressCallback([&](float, int complete, int, double) {
        if (complete > 0) ++publications;
    });
    ASSERT_EQ(session.Execute(), render::SessionState::Complete);
    EXPECT_GT(publications, 0U);
    const auto pixels = session.GetDisplayBuffer().SnapshotFloatData();
    ASSERT_EQ(pixels.size(), 16U);
    float total_radiance = 0;
    for (unsigned pixel = 0; pixel < 4; ++pixel) {
        for (unsigned channel = 0; channel < 3; ++channel) {
            EXPECT_TRUE(std::isfinite(pixels[4 * pixel + channel]));
            EXPECT_GE(pixels[4 * pixel + channel], 0);
            total_radiance += pixels[4 * pixel + channel];
        }
        EXPECT_EQ(pixels[4 * pixel + 3], 1);
    }
    EXPECT_GT(total_radiance, 0);
}

TEST(RenderCommandParse, BackendVulkanDeclinesMetricOffTheRenderPath) {
    RenderCommand cmd;
    GlobalOptions globals;
    SiriusConfig config = SiriusConfig::Defaults();

    // --backend vulkan enters the render session, which carries the registry
    // gpu_supported set. A charge metric (Reissner-Nordstrom,
    // gpu_supported=false) declines before device submission whether or not a
    // device is present, and it never falls back to CPU silently. Entering the
    // render path makes this a rendering-suite concern even without an image.
    const int rc = cmd.Execute(
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
    // Keep enough independent tiles to use the CPU workers at this resolution.
    config.session_template.tile_size = 8;
    config.session_template.thread_count = 2;
    config.session_template.enable_film_finish = true;
    // Bloom and the configured film finish must match the final display pipeline.
    config.session_template.enable_bloom = true;

    InteractiveViewer viewer;
    ASSERT_TRUE(viewer.Initialise(config));
    std::atomic<unsigned> previews{0}, completed_frames{0};
    std::atomic<bool> partial_preview{false};
    std::vector<float> last_preview;
    viewer.SetPreviewCallback(
        [&](const float* data, int width, int height, std::uint64_t generation) {
            EXPECT_EQ(generation, viewer.GetPreviewGeneration());
            EXPECT_EQ(width, 64);
            EXPECT_EQ(height, 64);
            EXPECT_EQ(completed_frames.load(), 0U);
            EXPECT_FALSE(viewer.GetRefinementState().complete);
            EXPECT_TRUE(viewer.GetFrameBufferSnapshot().empty())
                << "provisional radiance must not commit the final frame owner";
            last_preview.assign(data, data + static_cast<std::size_t>(width) * height * 4);
            EXPECT_TRUE(std::all_of(last_preview.begin(), last_preview.end(),
                                    [](float value) { return std::isfinite(value); }));
            bool completed_alpha = false, unfinished_alpha = false;
            for (std::size_t alpha = 3; alpha < last_preview.size(); alpha += 4) {
                EXPECT_TRUE(last_preview[alpha] == 0 || last_preview[alpha] == 1);
                completed_alpha |= last_preview[alpha] == 1;
                unfinished_alpha |= last_preview[alpha] == 0;
            }
            if (completed_alpha && unfinished_alpha) partial_preview = true;
            ++previews;
        });
    viewer.SetFrameCallback([&](const float*, int, int) { ++completed_frames; });
    ASSERT_TRUE(viewer.Start());

    // Joint retained transport does substantially more work than the former
    // scalar preview. This bound permits a complete physical frame, including
    // sanitizer overhead; dedicated cancellation tests bound interruption.
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::minutes(30);
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
    EXPECT_GT(previews.load(), 0U);
    EXPECT_TRUE(partial_preview.load())
        << "preview must expose completed tiles while other physical tiles remain unfinished";
    EXPECT_EQ(completed_frames.load(), 1U);
    // The last complete radiance preview runs the exact configured tone/bloom/
    // grade pipeline, while its private owner never feeds back into the final.
    EXPECT_EQ(last_preview, viewer.GetFrameBufferSnapshot());
    RecordProperty("preview_scope",
                   "in-flight preview has both completed/unfinished alpha ledger; separate "
                   "owner uses configured bloom/film pipeline before final frame callback; "
                   "this is not a frame-speed gate");

    ASSERT_TRUE(viewer.Initialise(config));
    completed_frames = 0;
    std::vector<std::uint64_t> generations;
    std::promise<void> stopped;
    auto stopped_future = stopped.get_future();
    viewer.SetPreviewCallback([&](const float*, int, int, std::uint64_t generation) {
        EXPECT_EQ(generation, viewer.GetPreviewGeneration());
        EXPECT_TRUE(viewer.GetFrameBufferSnapshot().empty());
        generations.push_back(generation);
        if (generations.size() == 1) {
            viewer.Restart();
            EXPECT_NE(generation, viewer.GetPreviewGeneration());
        } else {
            EXPECT_GT(generation, generations.front());
            viewer.SetPreviewCallback({});
            viewer.Stop();  // Reentrant owner callback requests stop without self-join.
            EXPECT_NE(generation, viewer.GetPreviewGeneration());
            stopped.set_value();
        }
    });
    ASSERT_TRUE(viewer.Start());
    const bool callback_stopped =
        stopped_future.wait_for(std::chrono::minutes(30)) == std::future_status::ready;
    viewer.Stop();  // External owner joins; captures remain alive until this returns.
    EXPECT_TRUE(callback_stopped);
    ASSERT_EQ(generations.size(), 2U);
    EXPECT_EQ(completed_frames.load(), 0U);
    EXPECT_TRUE(viewer.GetFrameBufferSnapshot().empty());
    EXPECT_TRUE(viewer.GetLastError().empty()) << viewer.GetLastError();
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

    // Both frames use joint retained transport. On Dozen the interval work
    // dominates the old scalar preview's shader-translation cost. This bound
    // permits complete physical frames; cancellation has separate interval-
    // and session-level checks, and throughput needs native qualification.
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::hours(6);
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

}  // namespace sirius::app::test
