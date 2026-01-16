// CRCF002A.cpp - Configuration Loader
// Component ID: CRCF002A (Configuration/Loader)

#include "CRCF002A.h"
#include "CRPF001A.h"

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

    // Render validation
    if (config.render.width < 1 || config.render.width > 16384) {
        errors.push_back("render.width must be between 1 and 16384");
    }
    if (config.render.height < 1 || config.render.height > 16384) {
        errors.push_back("render.height must be between 1 and 16384");
    }
    if (config.render.samplesPerPixel < 1 || config.render.samplesPerPixel > 65536) {
        errors.push_back("render.samplesPerPixel must be between 1 and 65536");
    }
    if (config.render.tileSize < 8 || config.render.tileSize > 512) {
        errors.push_back("render.tileSize must be between 8 and 512");
    }

    // Metric validation
    static const std::vector<std::string> validMetrics = {
        "Minkowski", "Schwarzschild", "Kerr", "Reissner-Nordstrom", "Kerr-Newman"
    };
    bool validMetric = false;
    for (const auto& m : validMetrics) {
        if (config.metric.name == m) {
            validMetric = true;
            break;
        }
    }
    if (!validMetric) {
        errors.push_back("metric.name must be one of: Minkowski, Schwarzschild, Kerr, Reissner-Nordstrom, Kerr-Newman");
    }
    if (config.metric.mass <= 0) {
        errors.push_back("metric.mass must be positive");
    }
    if (config.metric.spin < 0 || config.metric.spin > 1) {
        errors.push_back("metric.spin must be between 0 and 1");
    }
    if (config.metric.charge < 0) {
        errors.push_back("metric.charge must be non-negative");
    }

    // Observer validation
    if (config.observer.distance <= 0) {
        errors.push_back("observer.distance must be positive");
    }
    if (config.observer.inclination < 0 || config.observer.inclination > 180) {
        errors.push_back("observer.inclination must be between 0 and 180 degrees");
    }
    if (config.observer.fov < 1 || config.observer.fov > 180) {
        errors.push_back("observer.fov must be between 1 and 180 degrees");
    }

    // Post-process validation
    if (config.postprocess.exposure <= 0 || config.postprocess.exposure > 100) {
        errors.push_back("postprocess.exposure must be between 0 and 100");
    }

    // Backend validation
    static const std::vector<std::string> validBackends = {"auto", "optix", "cpu"};
    bool validBackend = false;
    for (const auto& b : validBackends) {
        if (config.backend.preferred == b) {
            validBackend = true;
            break;
        }
    }
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
