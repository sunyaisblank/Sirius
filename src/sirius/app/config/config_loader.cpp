// Configuration loader and validator. Ported from CRCF002A.cpp.
// Validation tolerances trace to core/constants.h (mirrored here as named
// literals, exactly as the legacy loader documented them).

#include "sirius/app/config/config_loader.h"

#include "sirius/app/platform_paths.h"
#include "sirius/core/metrics/registry.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

namespace sirius::app {

std::optional<fs::path> ConfigLoader::loaded_path_ = std::nullopt;

SiriusConfig ConfigLoader::Load(const std::optional<std::string>& override_path) {
    SiriusConfig config = SiriusConfig::defaults();
    loaded_path_ = std::nullopt;

    std::optional<fs::path> config_path;

    if (override_path.has_value() && !override_path->empty()) {
        fs::path path(*override_path);
        if (fs::exists(path)) {
            config_path = path;
        }
    } else {
        config_path = PlatformPaths::FindConfigFile();
    }

    if (config_path.has_value()) {
        try {
            std::ifstream file(config_path.value());
            if (file) {
                nlohmann::json j;
                file >> j;
                MergeConfig(config, j);
                loaded_path_ = config_path;
            }
        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed to load config from " << config_path->string() << ": "
                      << e.what() << std::endl;
        }
    }

    ApplyEnvironmentOverrides(config);

    return config;
}

SiriusConfig ConfigLoader::LoadFromFile(const fs::path& path) {
    SiriusConfig config = SiriusConfig::defaults();

    try {
        std::ifstream file(path);
        if (file) {
            nlohmann::json j;
            file >> j;
            config = j.get<SiriusConfig>();
        }
    } catch (const std::exception& e) {
        std::cerr << "Error loading config from " << path.string() << ": " << e.what() << std::endl;
    }

    return config;
}

bool ConfigLoader::SaveToFile(const SiriusConfig& config, const fs::path& path) {
    try {
        if (path.has_parent_path()) {
            std::error_code ec;
            fs::create_directories(path.parent_path(), ec);
        }

        std::ofstream file(path);
        if (!file) {
            return false;
        }

        nlohmann::json j = config;
        file << j.dump(2);
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

void ConfigLoader::ApplyEnvironmentOverrides(SiriusConfig& config) {
    if (auto val = GetEnvInt("SIRIUS_WIDTH")) {
        config.render.width = *val;
    }
    if (auto val = GetEnvInt("SIRIUS_HEIGHT")) {
        config.render.height = *val;
    }
    if (auto val = GetEnvInt("SIRIUS_SAMPLES")) {
        config.render.samplesPerPixel = *val;
    }
    if (auto val = GetEnvInt("SIRIUS_TILE_SIZE")) {
        config.render.tileSize = *val;
    }
    if (auto val = GetEnvInt("SIRIUS_THREADS")) {
        config.render.threadCount = *val;
    }
    if (auto val = GetEnv("SIRIUS_OUTPUT")) {
        config.render.outputPath = *val;
    }

    if (auto val = GetEnv("SIRIUS_METRIC")) {
        config.metric.name = *val;
    }
    if (auto val = GetEnvDouble("SIRIUS_MASS")) {
        config.metric.mass = *val;
    }
    if (auto val = GetEnvDouble("SIRIUS_SPIN")) {
        config.metric.spin = *val;
    }
    if (auto val = GetEnvDouble("SIRIUS_CHARGE")) {
        config.metric.charge = *val;
    }

    if (auto val = GetEnvDouble("SIRIUS_DISTANCE")) {
        config.observer.distance = *val;
    }
    if (auto val = GetEnvDouble("SIRIUS_INCLINATION")) {
        config.observer.inclination = *val;
    }
    if (auto val = GetEnvDouble("SIRIUS_AZIMUTH")) {
        config.observer.azimuth = *val;
    }
    if (auto val = GetEnvDouble("SIRIUS_FOV")) {
        config.observer.fov = *val;
    }

    if (auto val = GetEnvBool("SIRIUS_BLOOM")) {
        config.postprocess.enableBloom = *val;
    }
    if (auto val = GetEnvDouble("SIRIUS_EXPOSURE")) {
        config.postprocess.exposure = static_cast<float>(*val);
    }

    if (auto val = GetEnv("SIRIUS_BACKEND")) {
        config.backend.preferred = *val;
    }
    if (auto val = GetEnvInt("SIRIUS_CUDA_DEVICE")) {
        config.backend.cudaDevice = *val;
    }
}

std::optional<fs::path> ConfigLoader::GetLoadedConfigPath() {
    return loaded_path_;
}

std::vector<std::string> ConfigLoader::Validate(const SiriusConfig& config) {
    std::vector<std::string> errors;

    // --- Render validation --------------------------------------------------
    constexpr int kMinResolution = 128;
    constexpr int kMaxResolution = 8192;

    if (config.render.width < kMinResolution || config.render.width > kMaxResolution) {
        errors.push_back("render.width must be between " + std::to_string(kMinResolution) + " and " +
                         std::to_string(kMaxResolution) + " (spec requirement)");
    }
    if (config.render.height < kMinResolution || config.render.height > kMaxResolution) {
        errors.push_back("render.height must be between " + std::to_string(kMinResolution) +
                         " and " + std::to_string(kMaxResolution) + " (spec requirement)");
    }

    if (config.render.samplesPerPixel < 1 || config.render.samplesPerPixel > 4096) {
        errors.push_back("render.samplesPerPixel must be between 1 and 4096");
    }

    if (config.render.tileSize < 8 || config.render.tileSize > 256) {
        errors.push_back("render.tileSize must be between 8 and 256");
    }
    if ((config.render.tileSize & (config.render.tileSize - 1)) != 0) {
        errors.push_back("render.tileSize should be a power of 2 for GPU efficiency");
    }

    // --- Metric validation (single authority: the identity registry) --------
    auto metric_id = core::ParseMetricName(config.metric.name);
    if (!metric_id.has_value()) {
        errors.push_back("metric.name '" + config.metric.name +
                         "' is not a known metric; accepted names: " + core::KnownMetricNames());
    }

    constexpr double kMinMass = 0.1;
    constexpr double kMaxMass = 100.0;
    const bool massless = metric_id.has_value() && (*metric_id == core::MetricId::Minkowski ||
                                                    *metric_id == core::MetricId::DeSitter);
    if (!massless && (config.metric.mass < kMinMass || config.metric.mass > kMaxMass)) {
        errors.push_back("metric.mass must be between 0.1 and 100 (geometric units)");
    }

    // Spin a/M capped at the near-extremal limit; extremal Kerr is singular.
    constexpr double kMaxSpin = 0.998;
    if (config.metric.spin < 0 || config.metric.spin > kMaxSpin) {
        errors.push_back("metric.spin must be between 0 and 0.998 (near-extremal limit)");
    }

    constexpr double kMaxCharge = 0.999;
    if (config.metric.charge < 0 || config.metric.charge > kMaxCharge) {
        errors.push_back("metric.charge must be between 0 and 0.999");
    }

    // Combined extremality a^2 + Q^2 < M^2 for a valid horizon.
    if (config.metric.spin * config.metric.spin + config.metric.charge * config.metric.charge >=
        0.999) {
        errors.push_back("Combined spin² + charge² must be < 0.999 (sub-extremal condition)");
    }

    static const std::vector<std::string> valid_temperature_models = {"NovikovThorne", "NT",
                                                                       "ShakuraSunyaev", "SS"};
    if (std::find(valid_temperature_models.begin(), valid_temperature_models.end(),
                  config.metric.temperatureModel) == valid_temperature_models.end()) {
        errors.push_back(
            "metric.temperatureModel must be one of: NovikovThorne, NT, ShakuraSunyaev, SS");
    }

    // Only the a = 0 Kerr-Schild form is exact, so lambda with spin is rejected.
    constexpr double kMaxLambda = 0.1;
    if (std::abs(config.metric.lambda) > kMaxLambda) {
        errors.push_back("metric.lambda must be between -0.1 and 0.1");
    }
    if (std::abs(config.metric.lambda) > 0 && config.metric.spin != 0) {
        errors.push_back(
            "metric.lambda requires metric.spin = 0 (rotating de Sitter forms are not represented)");
    }

    // --- Observer validation ------------------------------------------------
    constexpr double kMinDistanceFactor = 5.0;
    constexpr double kMaxDistanceFactor = 1000.0;

    double min_distance = kMinDistanceFactor * config.metric.mass;
    double max_distance = kMaxDistanceFactor * config.metric.mass;

    if (config.observer.distance < min_distance || config.observer.distance > max_distance) {
        errors.push_back("observer.distance must be between 5M and 1000M (currently " +
                         std::to_string(config.observer.distance) + ", M=" +
                         std::to_string(config.metric.mass) + ")");
    }

    constexpr double kPoleBuffer = 0.1;
    if (config.observer.inclination <= kPoleBuffer ||
        config.observer.inclination >= 180.0 - kPoleBuffer) {
        errors.push_back(
            "observer.inclination must be between 0.1 and 179.9 degrees "
            "(avoiding coordinate singularity at poles)");
    }

    if (config.observer.fov < 1 || config.observer.fov > 170) {
        errors.push_back("observer.fov must be between 1 and 170 degrees");
    }

    // --- Post-process validation --------------------------------------------
    if (config.postprocess.exposure <= 0 || config.postprocess.exposure > 100) {
        errors.push_back("postprocess.exposure must be between 0 and 100 stops");
    }
    if (config.postprocess.bloomIntensity < 0 || config.postprocess.bloomIntensity > 5) {
        errors.push_back("postprocess.bloomIntensity must be between 0 and 5");
    }
    static const std::vector<std::string> valid_tonemappers = {"ACES",      "Reinhard", "Filmic",
                                                               "Uncharted2", "None",     "Linear"};
    if (std::find(valid_tonemappers.begin(), valid_tonemappers.end(),
                  config.postprocess.tonemapper) == valid_tonemappers.end()) {
        errors.push_back("postprocess.tonemapper must be one of: ACES, Reinhard, Filmic, None");
    }

    // --- Backend validation -------------------------------------------------
    // OptiX is retired. Since the go-live flip (2026-07-18) the accepted set is
    // auto (resolves to Vulkan where a device and a gpu-dispatchable metric
    // align, CPU otherwise), cpu (pinned), and vulkan (explicit).
    static const std::vector<std::string> valid_backends = {"auto", "cpu", "vulkan"};
    bool valid_backend = std::find(valid_backends.begin(), valid_backends.end(),
                                   config.backend.preferred) != valid_backends.end();
    if (!valid_backend) {
        errors.push_back("backend.preferred must be one of: auto, cpu, vulkan");
    }

    return errors;
}

std::string ConfigLoader::GenerateDefaultConfig() {
    SiriusConfig config = SiriusConfig::defaults();
    nlohmann::json j = config;
    return j.dump(2);
}

void ConfigLoader::MergeConfig(SiriusConfig& target, const nlohmann::json& source) {
    nlohmann::json target_json = target;

    for (auto& [key, value] : source.items()) {
        if (target_json.contains(key) && target_json[key].is_object() && value.is_object()) {
            for (auto& [sub_key, sub_value] : value.items()) {
                target_json[key][sub_key] = sub_value;
            }
        } else {
            target_json[key] = value;
        }
    }

    target = target_json.get<SiriusConfig>();
}

std::optional<std::string> ConfigLoader::GetEnv(const std::string& name) {
    const char* value = std::getenv(name.c_str());
    if (value && value[0] != '\0') {
        return std::string(value);
    }
    return std::nullopt;
}

std::optional<int> ConfigLoader::GetEnvInt(const std::string& name) {
    if (auto str = GetEnv(name)) {
        try {
            return std::stoi(*str);
        } catch (...) {
            return std::nullopt;
        }
    }
    return std::nullopt;
}

std::optional<double> ConfigLoader::GetEnvDouble(const std::string& name) {
    if (auto str = GetEnv(name)) {
        try {
            return std::stod(*str);
        } catch (...) {
            return std::nullopt;
        }
    }
    return std::nullopt;
}

std::optional<bool> ConfigLoader::GetEnvBool(const std::string& name) {
    if (auto str = GetEnv(name)) {
        std::string lower = *str;
        std::transform(lower.begin(), lower.end(), lower.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

        if (lower == "true" || lower == "1" || lower == "yes" || lower == "on") {
            return true;
        }
        if (lower == "false" || lower == "0" || lower == "no" || lower == "off") {
            return false;
        }
    }
    return std::nullopt;
}

}  // namespace sirius::app
