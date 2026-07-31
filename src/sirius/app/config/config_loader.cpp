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
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_set>

namespace sirius::app {

namespace {

[[nodiscard]] std::string JoinErrors(const std::vector<std::string>& errors) {
    std::ostringstream joined;
    for (std::size_t i = 0; i < errors.size(); ++i) {
        if (i != 0) joined << "; ";
        joined << errors[i];
    }
    return joined.str();
}

[[nodiscard]] nlohmann::json ParseJsonStrict(std::istream& input) {
    std::vector<std::unordered_set<std::string>> object_keys;
    std::optional<std::string> duplicate;
    const nlohmann::json::parser_callback_t callback =
        [&object_keys, &duplicate](int depth, nlohmann::json::parse_event_t event,
                                   nlohmann::json& parsed) {
            const auto index = static_cast<std::size_t>(depth + 1);
            if (event == nlohmann::json::parse_event_t::object_start) {
                if (object_keys.size() <= index) object_keys.resize(index + 1);
                object_keys[index].clear();
            } else if (event == nlohmann::json::parse_event_t::key) {
                if (object_keys.size() <= static_cast<std::size_t>(depth)) {
                    object_keys.resize(static_cast<std::size_t>(depth) + 1);
                }
                const std::string key = parsed.get<std::string>();
                if (!object_keys[static_cast<std::size_t>(depth)].insert(key).second &&
                    !duplicate.has_value()) {
                    duplicate = key;
                }
            } else if (event == nlohmann::json::parse_event_t::object_end &&
                       object_keys.size() > index) {
                object_keys[index].clear();
            }
            return true;
        };
    nlohmann::json result = nlohmann::json::parse(input, callback);
    if (duplicate.has_value()) {
        throw std::invalid_argument("duplicate configuration field '" + *duplicate + "'");
    }
    return result;
}

void RequireKnownKeys(const nlohmann::json& object, std::string_view path,
                      std::initializer_list<std::string_view> allowed) {
    if (!object.is_object()) {
        throw std::invalid_argument(std::string(path) + " must be a JSON object");
    }
    for (const auto& [key, value] : object.items()) {
        (void)value;
        if (std::find(allowed.begin(), allowed.end(), key) == allowed.end()) {
            throw std::invalid_argument("unknown configuration field '" + std::string(path) + "." +
                                        key + "'");
        }
    }
}

void RequireKnownConfigShape(const nlohmann::json& source) {
    RequireKnownKeys(source, "config",
                     {"render", "metric", "observer", "postprocess", "backend", "volumetric",
                      "film", "motionBlur", "diskEnabled", "dopplerBeaming", "pointStarfield",
                      "rayBundles", "colorMode"});
    if (source.contains("render")) {
        RequireKnownKeys(
            source.at("render"), "render",
            {"width", "height", "samplesPerPixel", "tileSize", "threadCount", "outputPath"});
    }
    if (source.contains("metric")) {
        RequireKnownKeys(
            source.at("metric"), "metric",
            {"name", "mass", "spin", "charge", "lambda", "temperatureModel", "diskTemperature",
             "throatRadius", "wormholeTopology", "warpVelocity", "bubbleRadius", "bubbleSigma"});
    }
    if (source.contains("observer")) {
        RequireKnownKeys(
            source.at("observer"), "observer",
            {"distance", "inclination", "azimuth", "fov", "cameraBetaForward", "cameraBetaUp",
             "cameraBetaRight", "lensModel", "focalLength", "aperture", "focusDistance"});
    }
    if (source.contains("postprocess")) {
        RequireKnownKeys(source.at("postprocess"), "postprocess",
                         {"enableBloom", "bloomIntensity", "bloomThreshold", "exposure", "contrast",
                          "saturation", "tonemapper"});
    }
    if (source.contains("backend")) {
        RequireKnownKeys(source.at("backend"), "backend",
                         {"preferred", "enableDenoiser", "cudaDevice"});
    }
    if (source.contains("volumetric")) {
        RequireKnownKeys(source.at("volumetric"), "volumetric",
                         {"enabled", "hOverR", "hPower", "tauMidplane", "samples",
                          "enableTurbulence", "enableCorona"});
    }
    if (source.contains("motionBlur")) {
        RequireKnownKeys(source.at("motionBlur"), "motionBlur",
                         {"enabled", "shutterTime", "samples"});
    }
    if (source.contains("film")) {
        RequireKnownKeys(
            source.at("film"), "film",
            {"enabled", "preset", "grainIntensity", "halationStrength", "vignetteStrength"});
    }
}

}  // namespace

std::optional<fs::path> ConfigLoader::loaded_path_ = std::nullopt;

SiriusConfig ConfigLoader::Load(const std::optional<std::string>& override_path) {
    SiriusConfig config = SiriusConfig::defaults();
    loaded_path_ = std::nullopt;

    std::optional<fs::path> config_path;

    if (override_path.has_value() && !override_path->empty()) {
        fs::path path(*override_path);
        if (!fs::is_regular_file(path)) {
            throw std::runtime_error("explicit configuration file does not exist: " +
                                     path.string());
        }
        config_path = path;
    } else {
        config_path = PlatformPaths::FindConfigFile();
    }

    if (config_path.has_value()) {
        std::ifstream file(config_path.value());
        if (!file) {
            throw std::runtime_error("cannot open configuration file: " + config_path->string());
        }
        const nlohmann::json j = ParseJsonStrict(file);
        MergeConfig(config, j);
        loaded_path_ = config_path;
    }

    ApplyEnvironmentOverrides(config);
    const auto errors = Validate(config);
    if (!errors.empty()) {
        throw std::invalid_argument("invalid Sirius configuration: " + JoinErrors(errors));
    }

    return config;
}

SiriusConfig ConfigLoader::LoadFromFile(const fs::path& path) {
    if (!fs::is_regular_file(path)) {
        throw std::runtime_error("configuration file does not exist: " + path.string());
    }
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("cannot open configuration file: " + path.string());
    }
    const nlohmann::json j = ParseJsonStrict(file);
    SiriusConfig config = SiriusConfig::defaults();
    MergeConfig(config, j);
    const auto errors = Validate(config);
    if (!errors.empty()) {
        throw std::invalid_argument("invalid Sirius configuration: " + JoinErrors(errors));
    }
    return config;
}

bool ConfigLoader::SaveToFile(const SiriusConfig& config, const fs::path& path) {
    try {
        if (!Validate(config).empty()) return false;
        if (path.has_parent_path()) {
            std::error_code ec;
            fs::create_directories(path.parent_path(), ec);
            if (ec) return false;
        }

        std::ofstream file(path);
        if (!file) {
            return false;
        }

        nlohmann::json j = config;
        file << j.dump(2);
        file.flush();
        return file.good();
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
    if (auto val = GetEnv("SIRIUS_COLOR_MODE")) {
        config.colorMode = *val;
    }
    if (auto val = GetEnvInt("SIRIUS_CUDA_DEVICE")) {
        config.backend.cudaDevice = *val;
    }
}

std::optional<fs::path> ConfigLoader::GetLoadedConfigPath() { return loaded_path_; }

std::vector<std::string> ConfigLoader::Validate(const SiriusConfig& config) {
    std::vector<std::string> errors;
    const auto finite = [&errors](double value, std::string_view field) {
        if (!std::isfinite(value)) {
            errors.push_back(std::string(field) + " must be finite");
            return false;
        }
        return true;
    };

    // --- Render validation --------------------------------------------------
    constexpr int kMinResolution = 128;
    constexpr int kMaxResolution = 8192;

    if (config.render.width < kMinResolution || config.render.width > kMaxResolution) {
        errors.push_back("render.width must be between " + std::to_string(kMinResolution) +
                         " and " + std::to_string(kMaxResolution) + " (spec requirement)");
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
    if (config.render.threadCount < 0 || config.render.threadCount > 1024) {
        errors.push_back("render.threadCount must be between 0 (automatic) and 1024");
    }
    if (config.render.outputPath.empty() ||
        config.render.outputPath.find('\0') != std::string::npos) {
        errors.push_back("render.outputPath must not be empty");
    } else {
        std::string extension = fs::path(config.render.outputPath).extension().string();
        std::transform(extension.begin(), extension.end(), extension.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (extension != ".ppm" && extension != ".png" && extension != ".exr") {
            errors.push_back("render.outputPath extension must be .ppm, .png, or .exr");
        }
    }

    // --- Metric validation (single authority: the identity registry) --------
    auto metric_id = core::ParseMetricName(config.metric.name);
    if (!metric_id.has_value()) {
        errors.push_back("metric.name '" + config.metric.name +
                         "' is not a known metric; accepted names: " + core::KnownMetricNames());
    }
    if (metric_id.has_value()) {
        if (const auto issue = core::MetricParameterIssue(
                *metric_id, config.metric.spin, config.metric.charge, config.metric.lambda);
            issue.has_value()) {
            errors.emplace_back(*issue);
        }
    }

    constexpr double kMinMass = 0.1;
    constexpr double kMaxMass = 100.0;
    const bool massless = metric_id.has_value() && (*metric_id == core::MetricId::Minkowski ||
                                                    *metric_id == core::MetricId::DeSitter);
    if (finite(config.metric.mass, "metric.mass")) {
        if (massless && config.metric.mass != 0.0) {
            errors.push_back("metric.mass must be zero for Minkowski and de-Sitter");
        } else if (!massless && (config.metric.mass < kMinMass || config.metric.mass > kMaxMass)) {
            errors.push_back("metric.mass must be between 0.1 and 100 (geometric units)");
        }
    }

    // Spin a/M capped at the near-extremal limit; extremal Kerr is singular.
    constexpr double kMaxSpin = 0.998;
    if (finite(config.metric.spin, "metric.spin") &&
        (config.metric.spin < 0 || config.metric.spin > kMaxSpin)) {
        errors.push_back("metric.spin must be between 0 and 0.998 (near-extremal limit)");
    }

    constexpr double kMaxCharge = 0.999;
    if (finite(config.metric.charge, "metric.charge") &&
        (config.metric.charge < 0 || config.metric.charge > kMaxCharge)) {
        errors.push_back("metric.charge must be between 0 and 0.999");
    }

    // Combined extremality a^2 + Q^2 < M^2 for a valid horizon.
    if (std::isfinite(config.metric.spin) && std::isfinite(config.metric.charge) &&
        config.metric.spin * config.metric.spin + config.metric.charge * config.metric.charge >=
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
    if (finite(config.metric.lambda, "metric.lambda") &&
        std::abs(config.metric.lambda) > kMaxLambda) {
        errors.push_back("metric.lambda must be between -0.1 and 0.1");
    }
    if (std::abs(config.metric.lambda) > 0 && config.metric.spin != 0) {
        errors.push_back(
            "metric.lambda requires metric.spin = 0 (rotating de Sitter forms are not "
            "represented)");
    }
    if (finite(config.metric.diskTemperature, "metric.diskTemperature") &&
        (config.metric.diskTemperature < 100.0f || config.metric.diskTemperature > 1.0e8f)) {
        errors.push_back("metric.diskTemperature must be between 100 and 100000000 Kelvin");
    }
    if (config.diskEnabled && metric_id.has_value() &&
        core::DiskSupportFor(*metric_id) != core::DiskSupport::PageThorne) {
        errors.push_back(
            "diskEnabled must be false when the selected metric has no represented "
            "Page-Thorne accretion-disk emission model");
    }
    if (finite(config.metric.throatRadius, "metric.throatRadius") &&
        (config.metric.throatRadius <= 0.0 || config.metric.throatRadius > 1000.0)) {
        errors.push_back("metric.throatRadius must be greater than 0 and at most 1000");
    }
    if (config.metric.wormholeTopology != "OneSheetCapture" &&
        config.metric.wormholeTopology != "TwoSheet") {
        errors.push_back("metric.wormholeTopology must be one of: OneSheetCapture, TwoSheet");
    } else if (config.metric.wormholeTopology == "TwoSheet") {
        errors.push_back(
            "metric.wormholeTopology TwoSheet is not represented: Sirius currently renders "
            "one exterior sheet with the throat as a dark capture surface");
    }
    if (finite(config.metric.warpVelocity, "metric.warpVelocity") &&
        std::abs(config.metric.warpVelocity) > 10.0) {
        errors.push_back("metric.warpVelocity magnitude must be at most 10");
    }
    if (finite(config.metric.bubbleRadius, "metric.bubbleRadius") &&
        (config.metric.bubbleRadius <= 0.0 || config.metric.bubbleRadius > 1000.0)) {
        errors.push_back("metric.bubbleRadius must be greater than 0 and at most 1000");
    }
    if (finite(config.metric.bubbleSigma, "metric.bubbleSigma") &&
        (config.metric.bubbleSigma <= 0.0 || config.metric.bubbleSigma > 1000.0)) {
        errors.push_back("metric.bubbleSigma must be greater than 0 and at most 1000");
    }

    // --- Observer validation ------------------------------------------------
    constexpr double kMinDistanceFactor = 5.0;
    constexpr double kMaxDistanceFactor = 1000.0;

    // Flat/de-Sitter scenes intentionally project a zero mass into the typed
    // session. Their coordinates still need a finite operator scale, so use
    // one geometric unit instead of making the valid interval [0, 0].
    const double distance_scale = massless ? 1.0 : config.metric.mass;
    const double min_distance = kMinDistanceFactor * distance_scale;
    const double max_distance = kMaxDistanceFactor * distance_scale;

    if (finite(config.observer.distance, "observer.distance") &&
        (config.observer.distance < min_distance || config.observer.distance > max_distance)) {
        errors.push_back("observer.distance must be between 5M and 1000M (currently " +
                         std::to_string(config.observer.distance) +
                         ", M=" + std::to_string(config.metric.mass) + ")");
    }

    constexpr double kPoleBuffer = 0.1;
    if (finite(config.observer.inclination, "observer.inclination") &&
        (config.observer.inclination <= kPoleBuffer ||
         config.observer.inclination >= 180.0 - kPoleBuffer)) {
        errors.push_back(
            "observer.inclination must be between 0.1 and 179.9 degrees "
            "(avoiding coordinate singularity at poles)");
    }

    if (finite(config.observer.azimuth, "observer.azimuth") &&
        (config.observer.azimuth < -360.0 || config.observer.azimuth > 360.0)) {
        errors.push_back("observer.azimuth must be between -360 and 360 degrees");
    }
    if (finite(config.observer.fov, "observer.fov") &&
        (config.observer.fov < 1 || config.observer.fov > 170)) {
        errors.push_back("observer.fov must be between 1 and 170 degrees");
    }
    const bool beta_forward_finite =
        finite(config.observer.cameraBetaForward, "observer.cameraBetaForward");
    const bool beta_up_finite = finite(config.observer.cameraBetaUp, "observer.cameraBetaUp");
    const bool beta_right_finite =
        finite(config.observer.cameraBetaRight, "observer.cameraBetaRight");
    const double beta_squared =
        config.observer.cameraBetaForward * config.observer.cameraBetaForward +
        config.observer.cameraBetaUp * config.observer.cameraBetaUp +
        config.observer.cameraBetaRight * config.observer.cameraBetaRight;
    if (beta_forward_finite && beta_up_finite && beta_right_finite &&
        (!std::isfinite(beta_squared) || beta_squared >= 1.0)) {
        errors.push_back("observer camera beta magnitude must be less than 1");
    }
    static const std::vector<std::string> valid_lens_models = {"Pinhole", "ThinLens", "Fisheye"};
    if (std::find(valid_lens_models.begin(), valid_lens_models.end(), config.observer.lensModel) ==
        valid_lens_models.end()) {
        errors.push_back("observer.lensModel must be one of: Pinhole, ThinLens, Fisheye");
    }
    if (finite(config.observer.focalLength, "observer.focalLength") &&
        (config.observer.focalLength <= 0.0f || config.observer.focalLength > 10000.0f)) {
        errors.push_back("observer.focalLength must be greater than 0 and at most 10000");
    }
    if (finite(config.observer.aperture, "observer.aperture") &&
        (config.observer.aperture <= 0.0f || config.observer.aperture > 128.0f)) {
        errors.push_back("observer.aperture must be greater than 0 and at most 128");
    }
    if (finite(config.observer.focusDistance, "observer.focusDistance") &&
        (config.observer.focusDistance <= 0.0f || config.observer.focusDistance > 1.0e6f)) {
        errors.push_back("observer.focusDistance must be greater than 0 and at most 1000000");
    }

    // --- Post-process validation --------------------------------------------
    if (finite(config.postprocess.exposure, "postprocess.exposure") &&
        (config.postprocess.exposure <= 0 || config.postprocess.exposure > 100)) {
        errors.push_back("postprocess.exposure must be between 0 and 100 stops");
    }
    if (finite(config.postprocess.bloomIntensity, "postprocess.bloomIntensity") &&
        (config.postprocess.bloomIntensity < 0 || config.postprocess.bloomIntensity > 5)) {
        errors.push_back("postprocess.bloomIntensity must be between 0 and 5");
    }
    if (finite(config.postprocess.bloomThreshold, "postprocess.bloomThreshold") &&
        (config.postprocess.bloomThreshold < 0 || config.postprocess.bloomThreshold > 100)) {
        errors.push_back("postprocess.bloomThreshold must be between 0 and 100");
    }
    if (finite(config.postprocess.contrast, "postprocess.contrast") &&
        (config.postprocess.contrast < 0 || config.postprocess.contrast > 4)) {
        errors.push_back("postprocess.contrast must be between 0 and 4");
    }
    if (finite(config.postprocess.saturation, "postprocess.saturation") &&
        (config.postprocess.saturation < 0 || config.postprocess.saturation > 4)) {
        errors.push_back("postprocess.saturation must be between 0 and 4");
    }
    static const std::vector<std::string> valid_tonemappers = {"ACES",       "Reinhard", "Filmic",
                                                               "Uncharted2", "None",     "Linear"};
    if (std::find(valid_tonemappers.begin(), valid_tonemappers.end(),
                  config.postprocess.tonemapper) == valid_tonemappers.end()) {
        errors.push_back(
            "postprocess.tonemapper must be one of: ACES, Reinhard, Filmic, Uncharted2, None, "
            "Linear");
    }

    // --- Volumetric and film validation ------------------------------------
    if (finite(config.volumetric.hOverR, "volumetric.hOverR") &&
        (config.volumetric.hOverR <= 0 || config.volumetric.hOverR > 2)) {
        errors.push_back("volumetric.hOverR must be greater than 0 and at most 2");
    }
    if (finite(config.volumetric.hPower, "volumetric.hPower") &&
        (config.volumetric.hPower < -2 || config.volumetric.hPower > 4)) {
        errors.push_back("volumetric.hPower must be between -2 and 4");
    }
    if (finite(config.volumetric.tauMidplane, "volumetric.tauMidplane") &&
        (config.volumetric.tauMidplane < 0 || config.volumetric.tauMidplane > 1.0e6f)) {
        errors.push_back("volumetric.tauMidplane must be between 0 and 1000000");
    }
    if (config.volumetric.samples < 1 || config.volumetric.samples > 4096) {
        errors.push_back("volumetric.samples must be between 1 and 4096");
    }
    if ((config.volumetric.enableTurbulence || config.volumetric.enableCorona) &&
        !config.volumetric.enabled) {
        errors.push_back("volumetric.enabled must be true when turbulence or corona is enabled");
    }
    if (config.volumetric.enabled && !config.diskEnabled) {
        errors.push_back("diskEnabled must be true when volumetric disk rendering is enabled");
    }
    if (finite(config.motionBlur.shutterTime, "motionBlur.shutterTime") &&
        (config.motionBlur.shutterTime < 0.0f || config.motionBlur.shutterTime > 1000.0f)) {
        errors.push_back("motionBlur.shutterTime must be between 0 and 1000");
    }
    if (config.motionBlur.samples < 1 || config.motionBlur.samples > 4096) {
        errors.push_back("motionBlur.samples must be between 1 and 4096");
    }
    if (config.motionBlur.enabled && !config.diskEnabled) {
        errors.push_back("diskEnabled must be true when motion blur is enabled");
    }

    static const std::vector<std::string> valid_color_modes = {
        "TrueColor", "TemperatureMap", "RedshiftMap", "Narrowband", "Polarisation"};
    if (std::find(valid_color_modes.begin(), valid_color_modes.end(), config.colorMode) ==
        valid_color_modes.end()) {
        errors.push_back(
            "colorMode must be one of: TrueColor, TemperatureMap, RedshiftMap, Narrowband, "
            "Polarisation");
    }
    if (config.colorMode == "Polarisation") {
        if (!config.diskEnabled) {
            errors.push_back("colorMode Polarisation requires diskEnabled");
        }
        if (config.volumetric.enabled) {
            errors.push_back(
                "colorMode Polarisation requires the thin disk; volumetric.enabled must be false");
        }
        if (config.motionBlur.enabled) {
            errors.push_back(
                "colorMode Polarisation is not represented with temporal disk motion blur");
        }
        if (metric_id.has_value() && *metric_id != core::MetricId::Schwarzschild &&
            *metric_id != core::MetricId::Kerr) {
            errors.push_back(
                "colorMode Polarisation is represented only for Schwarzschild and Kerr");
        }
    }

    static const std::vector<std::string> valid_film_presets = {"Interstellar", "SpaceOdyssey2001",
                                                                "DigitalClean"};
    if (std::find(valid_film_presets.begin(), valid_film_presets.end(), config.film.preset) ==
        valid_film_presets.end()) {
        errors.push_back(
            "film.preset must be one of: Interstellar, SpaceOdyssey2001, DigitalClean");
    }
    if (finite(config.film.grainIntensity, "film.grainIntensity") &&
        (config.film.grainIntensity < 0 || config.film.grainIntensity > 1)) {
        errors.push_back("film.grainIntensity must be between 0 and 1");
    }
    if (finite(config.film.halationStrength, "film.halationStrength") &&
        (config.film.halationStrength < 0 || config.film.halationStrength > 5)) {
        errors.push_back("film.halationStrength must be between 0 and 5");
    }
    if (finite(config.film.vignetteStrength, "film.vignetteStrength") &&
        (config.film.vignetteStrength < 0 || config.film.vignetteStrength > 2)) {
        errors.push_back("film.vignetteStrength must be between 0 and 2");
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
    if (config.backend.enableDenoiser) {
        errors.push_back("backend.enableDenoiser is not implemented and must be false");
    }
    if (config.backend.cudaDevice != 0) {
        errors.push_back("backend.cudaDevice is retired with CUDA/OptiX and must be 0");
    }

    return errors;
}

std::string ConfigLoader::GenerateDefaultConfig() {
    SiriusConfig config = SiriusConfig::defaults();
    nlohmann::json j = config;
    return j.dump(2);
}

void ConfigLoader::MergeConfig(SiriusConfig& target, const nlohmann::json& source) {
    RequireKnownConfigShape(source);
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
            std::size_t consumed = 0;
            const int value = std::stoi(*str, &consumed);
            if (consumed != str->size()) throw std::invalid_argument("trailing characters");
            return value;
        } catch (...) {
            throw std::invalid_argument("invalid integer environment variable " + name + "='" +
                                        *str + "'");
        }
    }
    return std::nullopt;
}

std::optional<double> ConfigLoader::GetEnvDouble(const std::string& name) {
    if (auto str = GetEnv(name)) {
        try {
            std::size_t consumed = 0;
            const double value = std::stod(*str, &consumed);
            if (consumed != str->size() || !std::isfinite(value)) {
                throw std::invalid_argument("not a finite complete number");
            }
            return value;
        } catch (...) {
            throw std::invalid_argument("invalid numeric environment variable " + name + "='" +
                                        *str + "'");
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
        throw std::invalid_argument("invalid boolean environment variable " + name + "='" + *str +
                                    "' (accepted: true/false, 1/0, yes/no, on/off)");
    }
    return std::nullopt;
}

}  // namespace sirius::app
