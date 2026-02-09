// CRCF002A.cpp - Configuration Loader
// Component ID: CRCF002A (Configuration/Loader)
//
// Validation aligns with docs/specification.md tolerance requirements.

#include "CRCF002A.h"
#include "CRPF001A.h"
#include <PHCN001A.h>

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <cctype>

namespace Sirius::Configuration {

// Static member initialization
std::optional<fs::path> ConfigLoader::s_loadedPath = std::nullopt;

SiriusConfig ConfigLoader::load(const std::optional<std::string>& overridePath) {
    SiriusConfig config = SiriusConfig::defaults();
    s_loadedPath = std::nullopt;

    // Determine which config file to load
    std::optional<fs::path> configPath;

    if (overridePath.has_value() && !overridePath->empty()) {
        // Use override path if provided
        fs::path path(*overridePath);
        if (fs::exists(path)) {
            configPath = path;
        }
    } else {
        // Search standard locations
        configPath = Platform::PathResolver::findConfigFile();
    }

    // Load from file if found
    if (configPath.has_value()) {
        try {
            std::ifstream file(configPath.value());
            if (file) {
                nlohmann::json j;
                file >> j;
                mergeConfig(config, j);
                s_loadedPath = configPath;
            }
        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed to load config from "
                      << configPath->string() << ": " << e.what() << std::endl;
        }
    }

    // Apply environment variable overrides
    applyEnvironmentOverrides(config);

    return config;
}

SiriusConfig ConfigLoader::loadFromFile(const fs::path& path) {
    SiriusConfig config = SiriusConfig::defaults();

    try {
        std::ifstream file(path);
        if (file) {
            nlohmann::json j;
            file >> j;
            config = j.get<SiriusConfig>();
        }
    } catch (const std::exception& e) {
        std::cerr << "Error loading config from " << path.string()
                  << ": " << e.what() << std::endl;
    }

    return config;
}

bool ConfigLoader::saveToFile(const SiriusConfig& config, const fs::path& path) {
    try {
        // Ensure parent directory exists
        if (path.has_parent_path()) {
            std::error_code ec;
            fs::create_directories(path.parent_path(), ec);
        }

        std::ofstream file(path);
        if (!file) {
            return false;
        }

        nlohmann::json j = config;
        file << j.dump(2);  // Pretty print with 2-space indent
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

void ConfigLoader::applyEnvironmentOverrides(SiriusConfig& config) {
    // Render settings
    if (auto val = getEnvInt("SIRIUS_WIDTH")) {
        config.render.width = *val;
    }
    if (auto val = getEnvInt("SIRIUS_HEIGHT")) {
        config.render.height = *val;
    }
    if (auto val = getEnvInt("SIRIUS_SAMPLES")) {
        config.render.samplesPerPixel = *val;
    }
    if (auto val = getEnvInt("SIRIUS_TILE_SIZE")) {
        config.render.tileSize = *val;
    }
    if (auto val = getEnvInt("SIRIUS_THREADS")) {
        config.render.threadCount = *val;
    }
    if (auto val = getEnv("SIRIUS_OUTPUT")) {
        config.render.outputPath = *val;
    }

    // Metric settings
    if (auto val = getEnv("SIRIUS_METRIC")) {
        config.metric.name = *val;
    }
    if (auto val = getEnvDouble("SIRIUS_MASS")) {
        config.metric.mass = *val;
    }
    if (auto val = getEnvDouble("SIRIUS_SPIN")) {
        config.metric.spin = *val;
    }
    if (auto val = getEnvDouble("SIRIUS_CHARGE")) {
        config.metric.charge = *val;
    }

    // Observer settings
    if (auto val = getEnvDouble("SIRIUS_DISTANCE")) {
        config.observer.distance = *val;
    }
    if (auto val = getEnvDouble("SIRIUS_INCLINATION")) {
        config.observer.inclination = *val;
    }
    if (auto val = getEnvDouble("SIRIUS_AZIMUTH")) {
        config.observer.azimuth = *val;
    }
    if (auto val = getEnvDouble("SIRIUS_FOV")) {
        config.observer.fov = *val;
    }

    // Post-process settings
    if (auto val = getEnvBool("SIRIUS_BLOOM")) {
        config.postprocess.enableBloom = *val;
    }
    if (auto val = getEnvDouble("SIRIUS_EXPOSURE")) {
        config.postprocess.exposure = static_cast<float>(*val);
    }

    // Backend settings
    if (auto val = getEnv("SIRIUS_BACKEND")) {
        config.backend.preferred = *val;
    }
    if (auto val = getEnvInt("SIRIUS_CUDA_DEVICE")) {
        config.backend.cudaDevice = *val;
    }
}

std::optional<fs::path> ConfigLoader::getLoadedConfigPath() {
    return s_loadedPath;
}

std::vector<std::string> ConfigLoader::validate(const SiriusConfig& config) {
    std::vector<std::string> errors;

    // =========================================================================
    // Render Validation (per docs/specification.md)
    // =========================================================================
    // Resolution constraints: [128, 8192] × [128, 8192]
    constexpr int MIN_RESOLUTION = 128;
    constexpr int MAX_RESOLUTION = 8192;

    if (config.render.width < MIN_RESOLUTION || config.render.width > MAX_RESOLUTION) {
        errors.push_back("render.width must be between " + std::to_string(MIN_RESOLUTION) +
                         " and " + std::to_string(MAX_RESOLUTION) + " (spec requirement)");
    }
    if (config.render.height < MIN_RESOLUTION || config.render.height > MAX_RESOLUTION) {
        errors.push_back("render.height must be between " + std::to_string(MIN_RESOLUTION) +
                         " and " + std::to_string(MAX_RESOLUTION) + " (spec requirement)");
    }

    // Samples per pixel: practical limits
    if (config.render.samplesPerPixel < 1 || config.render.samplesPerPixel > 4096) {
        errors.push_back("render.samplesPerPixel must be between 1 and 4096");
    }

    // Tile size: must be power of 2 for GPU efficiency
    if (config.render.tileSize < 8 || config.render.tileSize > 256) {
        errors.push_back("render.tileSize must be between 8 and 256");
    }
    if ((config.render.tileSize & (config.render.tileSize - 1)) != 0) {
        errors.push_back("render.tileSize should be a power of 2 for GPU efficiency");
    }

    // =========================================================================
    // Metric Validation (per docs/specification.md Parameter Ranges)
    // =========================================================================
    static const std::vector<std::string> validMetrics = {
        "Minkowski", "Schwarzschild", "Kerr", "Reissner-Nordstrom", "Kerr-Newman",
        "Morris-Thorne", "MorrisThorne", "Alcubierre", "Wormhole", "WarpDrive"
    };
    bool validMetric = std::find(validMetrics.begin(), validMetrics.end(),
                                  config.metric.name) != validMetrics.end();
    if (!validMetric) {
        errors.push_back("metric.name must be one of: Minkowski, Schwarzschild, Kerr, "
                         "Reissner-Nordstrom, Kerr-Newman, Morris-Thorne, Alcubierre");
    }

    // Mass M: [0.1, 100] in geometric units
    constexpr double MIN_MASS = 0.1;
    constexpr double MAX_MASS = 100.0;
    if (config.metric.mass < MIN_MASS || config.metric.mass > MAX_MASS) {
        errors.push_back("metric.mass must be between 0.1 and 100 (geometric units)");
    }

    // Spin parameter a/M: [0, 0.998] (near-extremal limit)
    // Extremal Kerr has coordinate singularity issues
    constexpr double MAX_SPIN = 0.998;
    if (config.metric.spin < 0 || config.metric.spin > MAX_SPIN) {
        errors.push_back("metric.spin must be between 0 and 0.998 (near-extremal limit)");
    }

    // Charge Q/M: [0, 0.999] (near-extremal limit)
    constexpr double MAX_CHARGE = 0.999;
    if (config.metric.charge < 0 || config.metric.charge > MAX_CHARGE) {
        errors.push_back("metric.charge must be between 0 and 0.999");
    }

    // Combined extremality: a² + Q² < M² for valid horizons
    if (config.metric.spin * config.metric.spin +
        config.metric.charge * config.metric.charge >= 0.999) {
        errors.push_back("Combined spin² + charge² must be < 0.999 (sub-extremal condition)");
    }

    // =========================================================================
    // Observer Validation (per docs/specification.md)
    // =========================================================================
    // Observer distance: [5M, 1000M] for valid rendering
    constexpr double MIN_DISTANCE_FACTOR = 5.0;    // In units of M
    constexpr double MAX_DISTANCE_FACTOR = 1000.0; // In units of M

    double minDistance = MIN_DISTANCE_FACTOR * config.metric.mass;
    double maxDistance = MAX_DISTANCE_FACTOR * config.metric.mass;

    if (config.observer.distance < minDistance || config.observer.distance > maxDistance) {
        errors.push_back("observer.distance must be between 5M and 1000M (currently " +
                         std::to_string(config.observer.distance) + ", M=" +
                         std::to_string(config.metric.mass) + ")");
    }

    // Inclination: (0, 180) degrees, avoiding poles
    constexpr double POLE_BUFFER = 0.1; // Degrees from pole
    if (config.observer.inclination <= POLE_BUFFER ||
        config.observer.inclination >= 180.0 - POLE_BUFFER) {
        errors.push_back("observer.inclination must be between 0.1 and 179.9 degrees "
                         "(avoiding coordinate singularity at poles)");
    }

    // FOV: reasonable range for rendering
    if (config.observer.fov < 1 || config.observer.fov > 170) {
        errors.push_back("observer.fov must be between 1 and 170 degrees");
    }

    // =========================================================================
    // Post-process Validation
    // =========================================================================
    if (config.postprocess.exposure <= 0 || config.postprocess.exposure > 100) {
        errors.push_back("postprocess.exposure must be between 0 and 100 stops");
    }
    if (config.postprocess.bloomIntensity < 0 || config.postprocess.bloomIntensity > 5) {
        errors.push_back("postprocess.bloomIntensity must be between 0 and 5");
    }

    // =========================================================================
    // Backend Validation
    // =========================================================================
    static const std::vector<std::string> validBackends = {"auto", "optix", "cpu"};
    bool validBackend = std::find(validBackends.begin(), validBackends.end(),
                                   config.backend.preferred) != validBackends.end();
    if (!validBackend) {
        errors.push_back("backend.preferred must be one of: auto, optix, cpu");
    }

    return errors;
}

std::string ConfigLoader::generateDefaultConfig() {
    SiriusConfig config = SiriusConfig::defaults();
    nlohmann::json j = config;
    return j.dump(2);
}

void ConfigLoader::mergeConfig(SiriusConfig& target, const nlohmann::json& source) {
    // Use nlohmann::json's update mechanism for partial merging
    nlohmann::json targetJson = target;

    // Recursively merge source into target
    for (auto& [key, value] : source.items()) {
        if (targetJson.contains(key) && targetJson[key].is_object() && value.is_object()) {
            // Recursively merge objects
            for (auto& [subKey, subValue] : value.items()) {
                targetJson[key][subKey] = subValue;
            }
        } else {
            targetJson[key] = value;
        }
    }

    target = targetJson.get<SiriusConfig>();
}

std::optional<std::string> ConfigLoader::getEnv(const std::string& name) {
    const char* value = std::getenv(name.c_str());
    if (value && value[0] != '\0') {
        return std::string(value);
    }
    return std::nullopt;
}

std::optional<int> ConfigLoader::getEnvInt(const std::string& name) {
    if (auto str = getEnv(name)) {
        try {
            return std::stoi(*str);
        } catch (...) {
            return std::nullopt;
        }
    }
    return std::nullopt;
}

std::optional<double> ConfigLoader::getEnvDouble(const std::string& name) {
    if (auto str = getEnv(name)) {
        try {
            return std::stod(*str);
        } catch (...) {
            return std::nullopt;
        }
    }
    return std::nullopt;
}

std::optional<bool> ConfigLoader::getEnvBool(const std::string& name) {
    if (auto str = getEnv(name)) {
        std::string lower = *str;
        std::transform(lower.begin(), lower.end(), lower.begin(),
                       [](unsigned char c) { return std::tolower(c); });

        if (lower == "true" || lower == "1" || lower == "yes" || lower == "on") {
            return true;
        }
        if (lower == "false" || lower == "0" || lower == "no" || lower == "off") {
            return false;
        }
    }
    return std::nullopt;
}

} // namespace Sirius::Configuration
