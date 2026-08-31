#include "sirius/app/config/session_config_adapter.h"

#include "sirius/app/config/config_loader.h"
#include "sirius/core/constants.h"
#include "sirius/core/metrics/registry.h"
#include "sirius/core/postprocess.h"
#include "sirius/render/vulkan_renderer.h"

#ifdef SIRIUS_HAS_VULKAN_BACKEND
#include "sirius/backend/device.h"
#endif

#include <iostream>
#include <optional>

namespace sirius::app {

namespace {

std::optional<render::DiskTemperatureModel> ParseTemperatureModel(const std::string& value) {
    if (value == "NovikovThorne" || value == "NT") {
        return render::DiskTemperatureModel::NovikovThorne;
    }
    if (value == "ShakuraSunyaev" || value == "SS") {
        return render::DiskTemperatureModel::ShakuraSunyaev;
    }
    return std::nullopt;
}

std::optional<core::color_modes::Mode> ParseColorMode(const std::string& value) {
    if (value == "TrueColor") return core::color_modes::Mode::TrueColor;
    if (value == "TemperatureMap") return core::color_modes::Mode::TemperatureMap;
    if (value == "RedshiftMap") return core::color_modes::Mode::RedshiftMap;
    if (value == "Polarisation") return core::color_modes::Mode::Polarisation;
    return std::nullopt;
}

std::optional<core::WormholeTopology> ParseWormholeTopology(const std::string& value) {
    if (value == "OneSheetCapture") return core::WormholeTopology::OneSheetCapture;
    if (value == "TwoSheet") return core::WormholeTopology::TwoSheet;
    return std::nullopt;
}

}  // namespace

base::Expected<render::FilmConfig> ProjectFilmFinishConfig(const FilmFinishConfig& config) {
    render::FilmConfig film;
    if (config.preset == "Interstellar") {
        film = render::FilmConfig::Interstellar();
    } else if (config.preset == "SpaceOdyssey2001") {
        film = render::FilmConfig::SpaceOdyssey2001();
    } else {
        return base::Fail(base::ErrorDomain::kConfiguration, "project film finish",
                          "unknown film finish preset '" + config.preset + "'");
    }
    if (config.grain_intensity.has_value()) film.grain_intensity = *config.grain_intensity;
    if (config.halation_strength.has_value()) film.halation_strength = *config.halation_strength;
    if (config.vignette_strength.has_value()) film.vignette_strength = *config.vignette_strength;
    return film;
}

base::Expected<render::SessionConfig> MakeSessionConfig(const SiriusConfig& config) {
    if (const auto errors = ConfigLoader::Validate(config); !errors.empty()) {
        return base::Fail(base::ErrorDomain::kConfiguration, "create render session",
                          errors.front());
    }

    render::SessionConfig session;
    session.width = config.render.width;
    session.height = config.render.height;
    session.samples_per_pixel = config.render.samples_per_pixel;
    session.tile_size = config.render.tile_size;
    session.thread_count = config.render.thread_count;
    session.output_path = config.render.output_path;
    session.enable_disk = config.disk_enabled;

    const auto metric_id = core::ParseMetricName(config.metric.name);
    if (!metric_id.has_value()) {
        return base::Fail(base::ErrorDomain::kConfiguration, "create render session",
                          "unknown metric name '" + config.metric.name + "'");
    }
    session.metric_id = *metric_id;
    session.black_hole_mass = config.metric.mass;
    session.black_hole_spin = config.metric.spin;
    session.black_hole_charge = config.metric.charge;
    session.cosmological_constant = config.metric.cosmological_constant;
    const auto temperature_model = ParseTemperatureModel(config.metric.temperature_model);
    if (!temperature_model) {
        return base::Fail(
            base::ErrorDomain::kConfiguration, "create render session",
            "unknown disk temperature model '" + config.metric.temperature_model + "'");
    }
    session.temperature_model = *temperature_model;
    session.disk_temperature_scale = config.metric.disk_temperature;
    session.throat_radius = config.metric.throat_radius;
    const auto wormhole_topology = ParseWormholeTopology(config.metric.wormhole_topology);
    if (!wormhole_topology) {
        return base::Fail(base::ErrorDomain::kConfiguration, "create render session",
                          "unknown wormhole topology '" + config.metric.wormhole_topology + "'");
    }
    session.wormhole_topology = *wormhole_topology;
    session.warp_velocity = config.metric.warp_velocity;
    session.bubble_radius = config.metric.bubble_radius;
    session.bubble_sigma = config.metric.bubble_sigma;

    session.observer_distance = config.observer.distance;
    session.observer_inclination = config.inclination_radians();
    session.observer_azimuth = config.observer.azimuth * core::constants::math::kPi / 180.0;
    session.camera_fov = static_cast<float>(config.observer.fov);
    session.camera_beta_forward = config.observer.camera_beta_forward;
    session.camera_beta_up = config.observer.camera_beta_up;
    session.camera_beta_right = config.observer.camera_beta_right;
    const auto lens_type = core::ParseLensType(config.observer.lens_model);
    if (!lens_type) {
        return base::Fail(base::ErrorDomain::kConfiguration, "create render session",
                          "unknown lens model '" + config.observer.lens_model + "'");
    }
    session.lens_type = *lens_type;
    session.camera_focal_length = config.observer.focal_length;
    session.camera_aperture = config.observer.aperture;
    session.camera_focus_distance = config.observer.focus_distance;

    const auto tonemapper = core::ParseTonemapType(config.postprocess.tonemapper);
    if (!tonemapper) {
        return base::Fail(base::ErrorDomain::kConfiguration, "create render session",
                          "unknown tonemapper '" + config.postprocess.tonemapper + "'");
    }
    session.tonemapper = *tonemapper;
    session.enable_bloom = config.postprocess.enable_bloom;
    session.bloom_intensity = config.postprocess.bloom_intensity;
    session.bloom_threshold = config.postprocess.bloom_threshold;
    session.exposure = config.postprocess.exposure;
    session.contrast = config.postprocess.contrast;
    session.saturation = config.postprocess.saturation;

    session.enable_volumetric_disk = config.volumetric.enabled;
    session.volumetric_h_over_r = config.volumetric.h_over_r;
    session.volumetric_h_power = config.volumetric.h_power;
    session.volumetric_tau_midplane = config.volumetric.tau_midplane;
    session.volumetric_samples = config.volumetric.samples;
    session.enable_turbulence = config.volumetric.enable_turbulence;
    session.enable_corona = config.volumetric.enable_corona;
    session.enable_motion_blur = config.motion_blur.enabled;
    session.shutter_time = config.motion_blur.shutter_time;
    session.motion_blur_samples = config.motion_blur.samples;

    session.enable_film_finish = config.film.enabled;
    const auto film_config = ProjectFilmFinishConfig(config.film);
    if (!film_config) {
        return std::unexpected(film_config.error());
    }
    session.film_config = *film_config;
    session.doppler_beaming = config.doppler_beaming;
    session.point_starfield = config.point_starfield;
    session.ray_bundles = config.ray_bundles;
    const auto color_mode = ParseColorMode(config.color_mode);
    if (!color_mode) {
        return base::Fail(base::ErrorDomain::kConfiguration, "create render session",
                          "unknown color mode '" + config.color_mode + "'");
    }
    session.color_mode = *color_mode;
    session.enable_polarisation = session.color_mode == core::color_modes::Mode::Polarisation;

    session.backend = render::RenderBackend::Cpu;
    if (config.backend.preferred == "vulkan") {
        session.backend = render::RenderBackend::Vulkan;
    } else if (config.backend.preferred == "auto") {
#ifdef SIRIUS_HAS_VULKAN_BACKEND
        if (!core::MetricInfoFor(session.metric_id).gpu_supported) {
            std::cout << "[Session] backend auto: metric '"
                      << core::MetricInfoFor(session.metric_id).canonical_name
                      << "' is CPU-only (registry); using the CPU path\n";
        } else if (auto compatible = render::ValidateVulkanRenderConfig(session); !compatible) {
            std::cout << "[Session] backend auto: " << compatible.error().detail()
                      << "; using the CPU path\n";
        } else if (auto devices = backend::EnumerateVulkanDevices();
                   devices.has_value() && !devices->empty()) {
            session.backend = render::RenderBackend::Vulkan;
        } else {
            std::cout << "[Session] backend auto: no Vulkan device visible; using the CPU path\n";
        }
#endif
    }
    return session;
}

}  // namespace sirius::app
