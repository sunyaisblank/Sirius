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
        stats.maximum_dispatch_rays = 64;
        stats.band_dispatches = imax ? 32736 : 3000;
        stats.retained_stage_dispatches =
            imax ? std::array<std::int64_t, 6>{0, 5000, 5000, 5000, 5000, 12736}
                 : std::array<std::int64_t, 6>{0, 500, 500, 500, 500, 1000};
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
        EXPECT_EQ(decoded["source_owner"], "host");
        EXPECT_EQ(decoded["route"], "retained");
        EXPECT_TRUE(decoded["dispatches"].is_number_integer());
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
}
