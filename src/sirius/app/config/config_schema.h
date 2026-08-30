#pragma once

// Operator-facing configuration belongs to the app layer. Its C++ names follow
// the Sirius style guide; the explicit JSON codec preserves the established
// camelCase file format without leaking those legacy spellings downstream.

#include "sirius/core/camera.h"
#include "sirius/core/disk/disk_defaults.h"
#include "sirius/core/feature_defaults.h"
#include "sirius/core/metrics/registry.h"

#include <nlohmann/json.hpp>

#include <numbers>
#include <optional>
#include <string>

namespace sirius::app {

struct RenderConfig {
    int width = 1920;
    int height = 1080;
    int samples_per_pixel = 64;
    int tile_size = 64;
    int thread_count = 0;
    std::string output_path = "render.ppm";
};

struct MetricConfig {
    std::string name = "Schwarzschild";
    double mass = 1.0;  // Geometric coordinate length; zero when M is not a parameter.
    double spin = 0.0;
    double charge = 0.0;
    double cosmological_constant = 0.0;
    std::string temperature_model = "NovikovThorne";
    float disk_temperature = core::kDefaultDiskTemperatureKelvin;
    double throat_radius = core::kDefaultMorrisThorneThroatRadius;
    std::string wormhole_topology = "OneSheetCapture";
    double warp_velocity = core::kDefaultAlcubierreWarpVelocity;
    double bubble_radius = core::kDefaultAlcubierreBubbleRadius;
    double bubble_sigma = core::kDefaultAlcubierreBubbleSigma;  // Inverse coordinate length.
};

struct ObserverConfig {
    double distance = 50.0;  // Coordinate radius r, not the dimensionless ratio r/M.
    double inclination = 90.0;
    double azimuth = 0.0;
    double fov = 60.0;
    double camera_beta_forward = 0.0;
    double camera_beta_up = 0.0;
    double camera_beta_right = 0.0;
    std::string lens_model = "Pinhole";
    float focal_length = core::kDefaultCameraFocalLength;
    float aperture = core::kDefaultCameraAperture;
    float focus_distance =
        core::kDefaultCameraFocusDistance;  // Geometric coordinate length from the camera.
};

struct PostProcessConfig {
    bool enable_bloom = true;
    float bloom_intensity = core::kDefaultBloomIntensity;
    float bloom_threshold = core::kDefaultBloomThreshold;
    float exposure = 1.0f;
    float contrast = 1.0f;
    float saturation = 1.0f;
    std::string tonemapper = "ACES";
};

struct VolumetricConfig {
    bool enabled = false;
    float h_over_r = core::kDefaultVolumetricHOverR;  // H/r at the disk inner edge.
    float h_power = core::kDefaultVolumetricHPower;   // Radial exponent before saturation.
    float tau_midplane =
        core::kDefaultVolumetricTauMidplane;  // Full vertical depth at the inner edge.
    int samples = core::kDefaultVolumetricSamples;
    bool enable_turbulence = false;
    bool enable_corona = false;
};

struct MotionBlurConfig {
    bool enabled = false;
    float shutter_time = core::kDefaultMotionBlurShutterTime;
    int samples = core::kDefaultMotionBlurSamples;
};

struct FilmFinishConfig {
    bool enabled = false;
    std::string preset = "Interstellar";
    // Absent means inherit the selected preset. This distinction is required:
    // a concrete default here would overwrite every non-default preset even
    // when the operator supplied no override.
    std::optional<float> grain_intensity;
    std::optional<float> halation_strength;
    std::optional<float> vignette_strength;
};

struct BackendConfig {
    std::string preferred = "auto";
    bool enable_denoiser = false;
    int cuda_device = 0;
};

struct SiriusConfig {
    RenderConfig render;
    MetricConfig metric;
    ObserverConfig observer;
    PostProcessConfig postprocess;
    BackendConfig backend;
    VolumetricConfig volumetric;
    MotionBlurConfig motion_blur;
    FilmFinishConfig film;
    bool disk_enabled = true;
    bool doppler_beaming = true;
    bool point_starfield = false;
    bool ray_bundles = false;
    std::string color_mode = "TrueColor";

    [[nodiscard]] static SiriusConfig Defaults() { return {}; }
    [[nodiscard]] double inclination_radians() const {
        return observer.inclination * std::numbers::pi_v<double> / 180.0;
    }
    [[nodiscard]] double fov_radians() const {
        return observer.fov * std::numbers::pi_v<double> / 180.0;
    }
};

namespace detail {

template <class T>
concept JsonDeserializable = requires(const nlohmann::json& json, T& value) { json.get_to(value); };

template <JsonDeserializable T>
void Read(const nlohmann::json& json, const char* key, T& value) {
    if (const auto item = json.find(key); item != json.end()) {
        item->get_to(value);
    }
}

}  // namespace detail

inline void to_json(nlohmann::json& json, const RenderConfig& config) {
    json = {{"width", config.width},
            {"height", config.height},
            {"samplesPerPixel", config.samples_per_pixel},
            {"tileSize", config.tile_size},
            {"threadCount", config.thread_count},
            {"outputPath", config.output_path}};
}

inline void from_json(const nlohmann::json& json, RenderConfig& config) {
    config = {};
    detail::Read(json, "width", config.width);
    detail::Read(json, "height", config.height);
    detail::Read(json, "samplesPerPixel", config.samples_per_pixel);
    detail::Read(json, "tileSize", config.tile_size);
    detail::Read(json, "threadCount", config.thread_count);
    detail::Read(json, "outputPath", config.output_path);
}

inline void to_json(nlohmann::json& json, const MetricConfig& config) {
    json = {{"name", config.name},
            {"mass", config.mass},
            {"spin", config.spin},
            {"charge", config.charge},
            {"lambda", config.cosmological_constant},
            {"temperatureModel", config.temperature_model},
            {"diskTemperature", config.disk_temperature},
            {"throatRadius", config.throat_radius},
            {"wormholeTopology", config.wormhole_topology},
            {"warpVelocity", config.warp_velocity},
            {"bubbleRadius", config.bubble_radius},
            {"bubbleSigma", config.bubble_sigma}};
}

inline void from_json(const nlohmann::json& json, MetricConfig& config) {
    config = {};
    detail::Read(json, "name", config.name);
    detail::Read(json, "mass", config.mass);
    detail::Read(json, "spin", config.spin);
    detail::Read(json, "charge", config.charge);
    detail::Read(json, "lambda", config.cosmological_constant);
    detail::Read(json, "temperatureModel", config.temperature_model);
    detail::Read(json, "diskTemperature", config.disk_temperature);
    detail::Read(json, "throatRadius", config.throat_radius);
    detail::Read(json, "wormholeTopology", config.wormhole_topology);
    detail::Read(json, "warpVelocity", config.warp_velocity);
    detail::Read(json, "bubbleRadius", config.bubble_radius);
    detail::Read(json, "bubbleSigma", config.bubble_sigma);
}

inline void to_json(nlohmann::json& json, const ObserverConfig& config) {
    json = {{"distance", config.distance},
            {"inclination", config.inclination},
            {"azimuth", config.azimuth},
            {"fov", config.fov},
            {"cameraBetaForward", config.camera_beta_forward},
            {"cameraBetaUp", config.camera_beta_up},
            {"cameraBetaRight", config.camera_beta_right},
            {"lensModel", config.lens_model},
            {"focalLength", config.focal_length},
            {"aperture", config.aperture},
            {"focusDistance", config.focus_distance}};
}

inline void from_json(const nlohmann::json& json, ObserverConfig& config) {
    config = {};
    detail::Read(json, "distance", config.distance);
    detail::Read(json, "inclination", config.inclination);
    detail::Read(json, "azimuth", config.azimuth);
    detail::Read(json, "fov", config.fov);
    detail::Read(json, "cameraBetaForward", config.camera_beta_forward);
    detail::Read(json, "cameraBetaUp", config.camera_beta_up);
    detail::Read(json, "cameraBetaRight", config.camera_beta_right);
    detail::Read(json, "lensModel", config.lens_model);
    detail::Read(json, "focalLength", config.focal_length);
    detail::Read(json, "aperture", config.aperture);
    detail::Read(json, "focusDistance", config.focus_distance);
}

inline void to_json(nlohmann::json& json, const PostProcessConfig& config) {
    json = {{"enableBloom", config.enable_bloom},
            {"bloomIntensity", config.bloom_intensity},
            {"bloomThreshold", config.bloom_threshold},
            {"exposure", config.exposure},
            {"contrast", config.contrast},
            {"saturation", config.saturation},
            {"tonemapper", config.tonemapper}};
}

inline void from_json(const nlohmann::json& json, PostProcessConfig& config) {
    config = {};
    detail::Read(json, "enableBloom", config.enable_bloom);
    detail::Read(json, "bloomIntensity", config.bloom_intensity);
    detail::Read(json, "bloomThreshold", config.bloom_threshold);
    detail::Read(json, "exposure", config.exposure);
    detail::Read(json, "contrast", config.contrast);
    detail::Read(json, "saturation", config.saturation);
    detail::Read(json, "tonemapper", config.tonemapper);
}

inline void to_json(nlohmann::json& json, const VolumetricConfig& config) {
    json = {{"enabled", config.enabled},
            {"hOverR", config.h_over_r},
            {"hPower", config.h_power},
            {"tauMidplane", config.tau_midplane},
            {"samples", config.samples},
            {"enableTurbulence", config.enable_turbulence},
            {"enableCorona", config.enable_corona}};
}

inline void from_json(const nlohmann::json& json, VolumetricConfig& config) {
    config = {};
    detail::Read(json, "enabled", config.enabled);
    detail::Read(json, "hOverR", config.h_over_r);
    detail::Read(json, "hPower", config.h_power);
    detail::Read(json, "tauMidplane", config.tau_midplane);
    detail::Read(json, "samples", config.samples);
    detail::Read(json, "enableTurbulence", config.enable_turbulence);
    detail::Read(json, "enableCorona", config.enable_corona);
}

inline void to_json(nlohmann::json& json, const MotionBlurConfig& config) {
    json = {{"enabled", config.enabled},
            {"shutterTime", config.shutter_time},
            {"samples", config.samples}};
}

inline void from_json(const nlohmann::json& json, MotionBlurConfig& config) {
    config = {};
    detail::Read(json, "enabled", config.enabled);
    detail::Read(json, "shutterTime", config.shutter_time);
    detail::Read(json, "samples", config.samples);
}

inline void to_json(nlohmann::json& json, const FilmFinishConfig& config) {
    json = {{"enabled", config.enabled}, {"preset", config.preset}};
    if (config.grain_intensity.has_value()) json["grainIntensity"] = *config.grain_intensity;
    if (config.halation_strength.has_value()) json["halationStrength"] = *config.halation_strength;
    if (config.vignette_strength.has_value()) json["vignetteStrength"] = *config.vignette_strength;
}

inline void from_json(const nlohmann::json& json, FilmFinishConfig& config) {
    config = {};
    detail::Read(json, "enabled", config.enabled);
    detail::Read(json, "preset", config.preset);
    if (const auto value = json.find("grainIntensity"); value != json.end())
        config.grain_intensity = value->get<float>();
    if (const auto value = json.find("halationStrength"); value != json.end())
        config.halation_strength = value->get<float>();
    if (const auto value = json.find("vignetteStrength"); value != json.end())
        config.vignette_strength = value->get<float>();
}

inline void to_json(nlohmann::json& json, const BackendConfig& config) {
    json = {{"preferred", config.preferred},
            {"enableDenoiser", config.enable_denoiser},
            {"cudaDevice", config.cuda_device}};
}

inline void from_json(const nlohmann::json& json, BackendConfig& config) {
    config = {};
    detail::Read(json, "preferred", config.preferred);
    detail::Read(json, "enableDenoiser", config.enable_denoiser);
    detail::Read(json, "cudaDevice", config.cuda_device);
}

inline void to_json(nlohmann::json& json, const SiriusConfig& config) {
    json = {{"render", config.render},
            {"metric", config.metric},
            {"observer", config.observer},
            {"postprocess", config.postprocess},
            {"backend", config.backend},
            {"volumetric", config.volumetric},
            {"motionBlur", config.motion_blur},
            {"film", config.film},
            {"diskEnabled", config.disk_enabled},
            {"dopplerBeaming", config.doppler_beaming},
            {"pointStarfield", config.point_starfield},
            {"rayBundles", config.ray_bundles},
            {"colorMode", config.color_mode}};
}

inline void from_json(const nlohmann::json& json, SiriusConfig& config) {
    config = {};
    detail::Read(json, "render", config.render);
    detail::Read(json, "metric", config.metric);
    detail::Read(json, "observer", config.observer);
    detail::Read(json, "postprocess", config.postprocess);
    detail::Read(json, "backend", config.backend);
    detail::Read(json, "volumetric", config.volumetric);
    detail::Read(json, "motionBlur", config.motion_blur);
    detail::Read(json, "film", config.film);
    detail::Read(json, "diskEnabled", config.disk_enabled);
    detail::Read(json, "dopplerBeaming", config.doppler_beaming);
    detail::Read(json, "pointStarfield", config.point_starfield);
    detail::Read(json, "rayBundles", config.ray_bundles);
    detail::Read(json, "colorMode", config.color_mode);
}

struct GlobalOptions {
    bool verbose = false;
    bool json_output = false;
    bool no_color = false;
    std::string config_path;
    bool show_help = false;
    bool show_version = false;
};

}  // namespace sirius::app
