#include "sirius/render/render_evidence.h"

#include "sirius/render/session/render_session.h"
#include "sirius/render/vulkan_renderer.h"

#include <gtest/gtest.h>

#include <nlohmann/json.hpp>

#include <iostream>

using namespace sirius::render;

// Synthetic protocol controls only: this case neither renders nor opens a GPU.
// The operational verifier consumes these records from the production writers,
// preventing its valid fixture from drifting away from the native wire format.
TEST(RenderEvidence, RetainedWireRecordsFeedAttestationControls) {
    for (const bool imax : {false, true}) {
        SessionConfig config;
        config.backend = RenderBackend::Vulkan;
        config.metric_id = sirius::core::MetricId::Kerr;
        config.black_hole_spin = .9;
        config.width = imax ? 5616 : 1920;
        config.height = imax ? 4096 : 1080;
        config.samples_per_pixel = 4;
        config.camera_fov = 60;
        config.enable_disk = false;
        config.point_starfield = true;
        config.ray_bundles = true;
        config.camera_beta_forward = .1;
        config.camera_beta_up = .02;
        config.camera_beta_right = -.01;
        config.lens_type = sirius::core::LensType::ThinLens;
        config.camera_focus_distance = 30;
        VulkanRenderStats stats;
        stats.device_name = "AMD Radeon 780M";
        stats.metric_name = "Kerr";
        stats.retained_intervals = true;
        stats.tile_plan.budget_bytes = 2147483648;
        stats.tile_plan.usable_bytes = 1073741824;
        stats.explicit_buffer_allocation_bytes = 67108864;
        stats.work_tile_edge = 32;
        stats.tiles_rendered = imax ? 22528 : 2040;
        stats.continuation_capacity = 64;
        stats.retained_timing.projection_capacity = 128;
        stats.maximum_dispatch_rays = 64;
        stats.band_dispatches = imax ? 32736 : 3000;
        stats.retained_stage_dispatches =
            imax ? std::array<std::int64_t, 7>{0, 5000, 5000, 5000, 5000, 12736, 0}
                 : std::array<std::int64_t, 7>{0, 500, 500, 500, 500, 1000, 0};
        auto& paired = stats.endpoint_dense_timing;
        paired.submissions = imax ? 3000 : 300;
        paired.submit_wait_ms = imax ? 30000 : 3000;
        paired.maximum_submit_wait_ms = 25;
        paired.pipeline_setup_ms = 3;
        paired.command_setup_ms = 4;
        paired.cleanup_ms = 5;
        paired.dispatch_total_ms = paired.submit_wait_ms + 12;
        paired.pipeline_creations = 2;
        paired.target_overshoots = 0;
        stats.queue_submissions =
            static_cast<std::uint64_t>(stats.band_dispatches) - paired.submissions;
        stats.camera_batches = imax ? 12736 : 1000;
        stats.accepted_intervals = imax ? 10000 : 1000;
        stats.dispatch_seconds = imax ? 90.0 : 15.0;
        stats.maximum_dispatch_ms = imax ? 75.0 : 50.0;
        stats.dispatch_target_overshoots = imax ? 3 : 2;
        stats.dispatch_subdivisions = imax ? 2 : 1;
        stats.initialization_seconds = imax ? .3 : 0;
        stats.seconds = imax ? 100.5 : 20.5;

        const auto scene = SessionSceneEvidenceJson(config, 100000);
        const auto completion = VulkanRenderEvidenceJson(config, stats);
        const auto decoded = nlohmann::json::parse(completion);
        EXPECT_EQ(decoded["work_items"], imax ? 22528 : 2040);
        EXPECT_EQ(decoded["work_tile_edge"], 32);
        EXPECT_EQ(decoded["maximum_dispatch_rays"], 64);
        EXPECT_EQ(decoded["ray_capacity"], 64);
        EXPECT_EQ(decoded["retained_timing"]["projection_capacity"], 128);
        EXPECT_EQ(decoded["source_owner"], "host");
        EXPECT_EQ(decoded["route"], "retained");
        EXPECT_TRUE(decoded["dispatches"].is_number_integer());
        EXPECT_EQ(decoded["dispatches"], stats.band_dispatches);
        EXPECT_EQ(decoded["queue_submissions"], stats.queue_submissions);
        EXPECT_EQ(decoded["shared_endpoint_dense"]["submissions"], paired.submissions);
        EXPECT_EQ(decoded["shared_endpoint_dense"]["submit_wait_ms"], paired.submit_wait_ms);
        EXPECT_EQ(decoded["shared_endpoint_dense"]["dispatch_total_ms"], paired.dispatch_total_ms);
        std::cout << kSceneEvidencePrefix << scene << '\n'
                  << kSourceSceneEvidencePrefix << scene << '\n'
                  << kVulkanEvidencePrefix << completion << '\n';
    }
}

TEST(RenderEvidence, DeviceIdentityEscapesJsonWithoutChangingItsValue) {
    SessionConfig config;
    VulkanRenderStats stats;
    stats.device_name = "adapter \"A\"\\driver\n\t";
    stats.device_index = 7;
    stats.precision = PrecisionRung::Fp64;
    const auto encoded = VulkanRenderEvidenceJson(config, stats);
    EXPECT_EQ(encoded.find('\n'), std::string::npos);
    const auto decoded = nlohmann::json::parse(encoded);
    EXPECT_EQ(decoded["device_name"], stats.device_name);
    EXPECT_EQ(decoded["device_index"], 7);
    EXPECT_EQ(decoded["precision"], "fp64");
    EXPECT_EQ(decoded["route"], "legacy");
    EXPECT_EQ(decoded["source_owner"], "device");
    EXPECT_FALSE(decoded.contains("retained_preparation"));

    // Initialization is visible without inventing governed rays or folding
    // its submission time into physical dispatch measurements.
    stats.retained_intervals = true;
    stats.initialization_dispatches = 5;
    stats.initialization_seconds = 2;
    stats.initialization_submit_wait_ms = 1500;
    stats.retained_preparation.wall_ms = 2000;
    auto& ray_camera = stats.retained_preparation.stages[static_cast<std::size_t>(
        sirius::backend::RetainedCompute::KernelStage::kRayCamera)];
    ray_camera.attempts = 1;
    ray_camera.dispatch_attempts = 1;
    ray_camera.completed_dispatches = 1;
    ray_camera.completed = 1;
    ray_camera.header_restored = true;
    ray_camera.timing.pipeline_setup_ms = 25;
    ray_camera.timing.command_setup_ms = 2;
    ray_camera.timing.submit_wait_ms = 1100;
    ray_camera.timing.cleanup_ms = 3;
    ray_camera.timing.total_ms = 1130;
    ray_camera.timing.pipeline_created = true;
    ray_camera.write_buffer_calls = 2;
    ray_camera.write_buffer_ms = 4;
    ray_camera.write_buffer_bytes = 8;
    const auto prepared = nlohmann::json::parse(VulkanRenderEvidenceJson(config, stats));
    EXPECT_EQ(prepared["initialization_dispatches"], 5);
    EXPECT_EQ(prepared["retained_preparation"]["wall_ms"], 2000);
    const auto& observed = prepared["retained_preparation"]["stages"][static_cast<std::size_t>(
        sirius::backend::RetainedCompute::KernelStage::kRayCamera)];
    EXPECT_EQ(observed["stage"], "ray_camera");
    EXPECT_EQ(observed["completed_dispatches"], 1);
    EXPECT_EQ(observed["header_restored"], true);
    EXPECT_EQ(observed["pipeline_setup_ms"], 25);
    EXPECT_EQ(observed["command_setup_ms"], 2);
    EXPECT_EQ(observed["submit_wait_ms"], 1100);
    EXPECT_EQ(observed["cleanup_ms"], 3);
    EXPECT_EQ(observed["dispatch_total_ms"], 1130);
    EXPECT_EQ(observed["pipeline_created"], true);
    EXPECT_EQ(observed["write_buffer_calls"], 2);
    EXPECT_EQ(observed["write_buffer_ms"], 4);
    EXPECT_EQ(observed["write_buffer_bytes"], 8);
    EXPECT_EQ(prepared["dispatches"], 0);
    EXPECT_EQ(prepared["dispatch_seconds"], 0);
    EXPECT_EQ(prepared["retained_timing"]["submit_wait_ms"], 0);
}
