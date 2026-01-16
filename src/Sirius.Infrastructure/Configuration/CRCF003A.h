// CRCF003A.h - Configuration Schema
// Component ID: CRCF003A (Configuration/Schema)
//
// Unified configuration structure with JSON serialization support.
// Replaces fragmented SessionConfig and RenderConfig with single schema.

#pragma once

#include <nlohmann/json.hpp>
#include <string>

namespace Sirius::Configuration {

/// @brief Render settings
struct RenderConfig {
    int width = 1920;
    int height = 1080;
    int samplesPerPixel = 64;
    int tileSize = 64;
    int threadCount = 0;  ///< 0 = auto-detect
    std::string outputPath = "render.ppm";

    NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(RenderConfig,
        width, height, samplesPerPixel, tileSize, threadCount, outputPath)
};

/// @brief Metric (spacetime) settings
struct MetricConfig {
    std::string name = "Schwarzschild";
    double mass = 1.0;
    double spin = 0.0;
    double charge = 0.0;

    NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(MetricConfig,
        name, mass, spin, charge)
};

/// @brief Observer/camera settings
struct ObserverConfig {
    double distance = 50.0;      ///< Distance from origin in M
    double inclination = 90.0;   ///< Inclination angle in degrees
    double azimuth = 0.0;        ///< Azimuthal angle in degrees
    double fov = 60.0;           ///< Field of view in degrees

    NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(ObserverConfig,
        distance, inclination, azimuth, fov)
};

/// @brief Post-processing settings
struct PostProcessConfig {
    bool enableBloom = true;
    float bloomIntensity = 0.3f;
    float exposure = 1.0f;
    std::string tonemapper = "ACES";  ///< ACES, Reinhard, Uncharted2

    NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(PostProcessConfig,
        enableBloom, bloomIntensity, exposure, tonemapper)
};

/// @brief Backend (renderer) settings
struct BackendConfig {
    std::string preferred = "auto";  ///< auto, optix, cpu
    bool enableDenoiser = false;
    int cudaDevice = 0;

    NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(BackendConfig,
        preferred, enableDenoiser, cudaDevice)
};

/// @brief Root configuration structure
struct SiriusConfig {
    RenderConfig render;
    MetricConfig metric;
    ObserverConfig observer;
    PostProcessConfig postprocess;
    BackendConfig backend;

    NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(SiriusConfig,
        render, metric, observer, postprocess, backend)

    /// @brief Get default configuration
    static SiriusConfig defaults() {
        return SiriusConfig{};
    }

    /// @brief Convert observer inclination from degrees to radians
    double inclinationRadians() const {
        return observer.inclination * 3.14159265358979323846 / 180.0;
    }

    /// @brief Convert FOV from degrees to radians
    double fovRadians() const {
        return observer.fov * 3.14159265358979323846 / 180.0;
    }
};

/// @brief Global CLI options (not persisted to config file)
struct GlobalOptions {
    bool verbose = false;
    bool jsonOutput = false;
    bool noColor = false;
    std::string configPath;  ///< Override config file path
    bool legacyMode = false; ///< Use legacy RenderJob
    bool showHelp = false;
    bool showVersion = false;
};

} // namespace Sirius::Configuration
