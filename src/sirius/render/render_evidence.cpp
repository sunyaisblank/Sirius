#include "sirius/render/render_evidence.h"

#include "sirius/render/session/render_session.h"
#include "sirius/render/vulkan_renderer.h"

#include <nlohmann/json.hpp>

#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>

namespace sirius::render {

std::string SessionSceneEvidenceJson(const SessionConfig& config, std::size_t point_star_count) {
    const char* backend = nullptr;
    switch (config.backend) {
        case RenderBackend::Cpu:
            backend = "Cpu";
            break;
        case RenderBackend::Vulkan:
            backend = "Vulkan";
            break;
        default:
            SIRIUS_ASSERT(false);
            backend = "Invalid";
            break;
    }

    const char* lens = nullptr;
    switch (config.lens_type) {
        case core::LensType::Pinhole:
            lens = "Pinhole";
            break;
        case core::LensType::ThinLens:
            lens = "ThinLens";
            break;
        case core::LensType::Fisheye:
            lens = "Fisheye";
            break;
        default:
            SIRIUS_ASSERT(false);
            lens = "Invalid";
            break;
    }

    std::ostringstream evidence;
    evidence.imbue(std::locale::classic());
    evidence << std::setprecision(std::numeric_limits<double>::max_digits10)
             << "{\"schema\":\"sirius-render-scene-v1\"" << ",\"backend\":\"" << backend << "\""
             << ",\"metric\":\"" << core::MetricInfoFor(config.metric_id).canonical_name << "\""
             << ",\"spin\":" << config.black_hole_spin << ",\"width\":" << config.width
             << ",\"height\":" << config.height
             << ",\"samples_per_pixel\":" << config.samples_per_pixel
             << ",\"field_of_view\":" << static_cast<double>(config.camera_fov)
             << ",\"disk_enabled\":" << (config.enable_disk ? "true" : "false")
             << ",\"ray_bundles\":" << (config.ray_bundles ? "true" : "false")
             << ",\"point_starfield\":" << (config.point_starfield ? "true" : "false")
             << ",\"point_star_count\":" << point_star_count << ",\"point_brightness_scale\":"
             << static_cast<double>(config.point_starfield_config.brightness_scale)
             << ",\"point_seed\":" << config.point_starfield_config.seed
             << ",\"point_min_distance_pc\":"
             << static_cast<double>(config.point_starfield_config.min_distance_pc)
             << ",\"point_max_distance_pc\":"
             << static_cast<double>(config.point_starfield_config.max_distance_pc)
             << ",\"camera_beta\":[" << config.camera_beta_forward << ',' << config.camera_beta_up
             << ',' << config.camera_beta_right << "]" << ",\"lens\":\"" << lens << "\""
             << ",\"focal_length\":" << static_cast<double>(config.camera_focal_length)
             << ",\"aperture\":" << static_cast<double>(config.camera_aperture)
             << ",\"focus_distance\":" << static_cast<double>(config.camera_focus_distance) << '}';
    return evidence.str();
}

std::string VulkanRenderEvidenceJson(const SessionConfig& config, const VulkanRenderStats& stats) {
    const char* precision = "invalid";
    switch (stats.precision) {
        case PrecisionRung::Fp32:
            precision = "fp32";
            break;
        case PrecisionRung::Fp32Comp:
            precision = "fp32-comp";
            break;
        case PrecisionRung::Fp64:
            precision = "fp64";
            break;
    }
    // Human progress is deliberately separate from this versioned wire record.
    // In particular, retained ray rows are neither pixels nor residency tiles.
    const nlohmann::ordered_json evidence = {
        {"schema", "sirius-vulkan-render-v1"},
        {"route", stats.retained_intervals ? "retained" : "legacy"},
        {"source_owner", stats.retained_intervals ? "host" : "device"},
        {"device_name", stats.device_name},
        {"device_index", stats.device_index},
        {"metric", stats.metric_name},
        {"width", config.width},
        {"height", config.height},
        {"precision", precision},
        {"budget_bytes", stats.tile_plan.budget_bytes},
        {"usable_bytes", stats.tile_plan.usable_bytes},
        {"allocated_bytes", stats.explicit_buffer_allocation_bytes},
        {"work_items", stats.tiles_rendered},
        {"work_tile_edge", stats.work_tile_edge},
        {"ray_capacity", stats.continuation_capacity},
        {"maximum_dispatch_rays", stats.maximum_dispatch_rays},
        {"dispatches", stats.band_dispatches},
        {"retained_stage_dispatches", stats.retained_stage_dispatches},
        {"camera_batches", stats.camera_batches},
        {"accepted_intervals", stats.accepted_intervals},
        {"wall_seconds", stats.seconds},
        {"dispatch_seconds", stats.dispatch_seconds},
        {"maximum_dispatch_ms", stats.maximum_dispatch_ms},
        {"target_overshoots", stats.dispatch_target_overshoots},
        {"batch_subdivisions", stats.dispatch_subdivisions},
        {"safety_fallbacks", stats.dispatch_fallbacks},
        {"initialization_dispatches", stats.initialization_dispatches},
        {"initialization_seconds", stats.initialization_seconds},
        {"initialization_submit_wait_ms", stats.initialization_submit_wait_ms}};
    return evidence.dump();
}

}  // namespace sirius::render
