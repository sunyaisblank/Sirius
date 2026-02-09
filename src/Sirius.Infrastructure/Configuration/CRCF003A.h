// CRCF003A.h - Configuration Schema
// Component ID: CRCF003A (Configuration/Schema)
//
// Unified configuration structure with JSON serialization support.
// Replaces fragmented SessionConfig and RenderConfig with single schema.

#pragma once

#include <nlohmann/json.hpp>
#include <PHCN001A.h>
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
    std::string temperatureModel = "NovikovThorne";  ///< NovikovThorne or ShakuraSunyaev
    float diskTemperature = 50000.0f;                ///< T_scale (Kelvin)

    // Morris-Thorne wormhole parameters
    double throatRadius = 1.0;    ///< b0: throat radius in M

    // Alcubierre warp drive parameters
    double warpVelocity = 0.5;    ///< vs: warp bubble velocity (sub-luminal)
    double bubbleRadius = 1.0;    ///< R: warp bubble radius
    double bubbleSigma = 0.5;     ///< sigma: wall thickness

    NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(MetricConfig,
        name, mass, spin, charge, temperatureModel, diskTemperature,
        throatRadius, warpVelocity, bubbleRadius, bubbleSigma)
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
    float bloomThreshold = 0.3f;
    float exposure = 1.0f;
    float contrast = 1.0f;
    float saturation = 1.0f;
    std::string tonemapper = "ACES";  ///< ACES, Reinhard, Uncharted2, Filmic, AgX

    NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(PostProcessConfig,
        enableBloom, bloomIntensity, bloomThreshold, exposure, contrast, saturation, tonemapper)
};

/// @brief Volumetric disk settings
struct VolumetricConfig {
    bool enabled = false;
    float hOverR = 0.1f;           ///< Scale height ratio H/r
    float hPower = 0.25f;          ///< Flaring index
    float tauMidplane = 10.0f;     ///< Midplane optical depth
    int samples = 64;              ///< Ray marching samples
    bool enableTurbulence = false; ///< Enable turbulent density perturbations
    bool enableCorona = false;     ///< Enable corona emission

    NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(VolumetricConfig,
        enabled, hOverR, hPower, tauMidplane, samples, enableTurbulence, enableCorona)
};

/// @brief Film simulation settings
struct FilmSimConfig {
    bool enabled = false;
    std::string preset = "Interstellar";  ///< Interstellar, SpaceOdyssey2001, DigitalClean
    float grainIntensity = 0.15f;
    float halationStrength = 0.15f;
    float vignetteStrength = 0.3f;

    NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(FilmSimConfig,
        enabled, preset, grainIntensity, halationStrength, vignetteStrength)
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
    VolumetricConfig volumetric;
    FilmSimConfig film;

    NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT(SiriusConfig,
        render, metric, observer, postprocess, backend, volumetric, film)

    /// @brief Get default configuration
    static SiriusConfig defaults() {
        return SiriusConfig{};
    }

    /// @brief Convert observer inclination from degrees to radians
    double inclinationRadians() const {
        return observer.inclination * Sirius::Constants::Math::PI / 180.0;
    }

    /// @brief Convert FOV from degrees to radians
    double fovRadians() const {
        return observer.fov * Sirius::Constants::Math::PI / 180.0;
    }
};

/// @brief Global CLI options (not persisted to config file)
struct GlobalOptions {
    bool verbose = false;
    bool jsonOutput = false;
    bool noColor = false;
    std::string configPath;  ///< Override config file path
    bool showHelp = false;
    bool showVersion = false;
};

} // namespace Sirius::Configuration
