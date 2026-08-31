// Configuration loader and validator. Ported from CRCF002A.cpp.
// Validation tolerances trace to core/constants.h (mirrored here as named
// literals, exactly as the legacy loader documented them).

#include "sirius/app/config/config_loader.h"

#include "sirius/app/platform_paths.h"
#include "sirius/core/camera.h"
#include "sirius/core/disk/disk_defaults.h"
#include "sirius/core/feature_defaults.h"
#include "sirius/core/metrics/registry.h"
#include "sirius/core/postprocess.h"

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
#include <utility>

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

void ApplyCompatibleMassDefault(MetricConfig& metric) {
    const auto id = core::ParseMetricName(metric.name);
    if (!id.has_value()) return;
    if (!core::MetricUsesMass(*id)) {
        metric.mass = 0.0;
    } else if (metric.mass == 0.0) {
        metric.mass = MetricConfig{}.mass;
    }
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

base::Expected<SiriusConfig> ConfigLoader::Load(const std::optional<std::string>& override_path) {
    SiriusConfig config = SiriusConfig::Defaults();
    loaded_path_ = std::nullopt;

    std::optional<fs::path> config_path;

    if (override_path.has_value() && !override_path->empty()) {
        fs::path path(*override_path);
        if (!fs::is_regular_file(path)) {
            return base::Fail(base::ErrorDomain::kIo, "load configuration",
                              "explicit file does not exist: " + path.string());
        }
        config_path = path;
    } else {
        config_path = PlatformPaths::FindConfigFile();
    }

    if (config_path.has_value()) {
        std::ifstream file(config_path.value());
        if (!file) {
            return base::Fail(base::ErrorDomain::kIo, "load configuration",
                              "cannot open file: " + config_path->string());
        }
        try {
            const nlohmann::json json = ParseJsonStrict(file);
            MergeConfig(config, json);
        } catch (const std::exception& error) {
            return base::Fail(base::ErrorDomain::kConfiguration, "parse configuration",
                              config_path->string() + ": " + error.what());
        }
    }

    if (auto environment = ApplyEnvironmentOverrides(config); !environment) {
        return std::unexpected(environment.error());
    }
    const auto errors = Validate(config);
    if (!errors.empty()) {
        return base::Fail(base::ErrorDomain::kConfiguration, "validate configuration",
                          JoinErrors(errors));
    }

    loaded_path_ = config_path;
    return config;
}

base::Expected<SiriusConfig> ConfigLoader::LoadFromFile(const fs::path& path) {
    if (!fs::is_regular_file(path)) {
        return base::Fail(base::ErrorDomain::kIo, "load configuration",
                          "file does not exist: " + path.string());
    }
    std::ifstream file(path);
    if (!file) {
        return base::Fail(base::ErrorDomain::kIo, "load configuration",
                          "cannot open file: " + path.string());
    }
    SiriusConfig config = SiriusConfig::Defaults();
    try {
        const nlohmann::json json = ParseJsonStrict(file);
        MergeConfig(config, json);
    } catch (const std::exception& error) {
        return base::Fail(base::ErrorDomain::kConfiguration, "parse configuration",
                          path.string() + ": " + error.what());
    }
    const auto errors = Validate(config);
    if (!errors.empty()) {
        return base::Fail(base::ErrorDomain::kConfiguration, "validate configuration",
                          path.string() + ": " + JoinErrors(errors));
    }
    return config;
}

base::Expected<void> ConfigLoader::SaveToFile(const SiriusConfig& config, const fs::path& path) {
    const auto errors = Validate(config);
    if (!errors.empty()) {
        return base::Fail(base::ErrorDomain::kConfiguration, "save configuration",
                          JoinErrors(errors));
    }
    try {
        if (path.has_parent_path()) {
            std::error_code ec;
            fs::create_directories(path.parent_path(), ec);
            if (ec) {
                return base::Fail(base::ErrorDomain::kIo, "create configuration directory",
                                  path.parent_path().string() + ": " + ec.message());
            }
        }

        std::ofstream file(path);
        if (!file) {
            return base::Fail(base::ErrorDomain::kIo, "open configuration for writing",
                              path.string());
        }

        nlohmann::json json = config;
        file << json.dump(2);
        file.flush();
        if (!file.good()) {
            return base::Fail(base::ErrorDomain::kIo, "write configuration", path.string());
        }
    } catch (const std::exception& error) {
        return base::Fail(base::ErrorDomain::kIo, "save configuration",
                          path.string() + ": " + error.what());
    }
    return {};
}

base::Expected<void> ConfigLoader::ApplyEnvironmentOverrides(SiriusConfig& config) {
    SiriusConfig updated = config;
    try {
        if (auto value = GetEnvInt("SIRIUS_WIDTH")) updated.render.width = *value;
        if (auto value = GetEnvInt("SIRIUS_HEIGHT")) updated.render.height = *value;
        if (auto value = GetEnvInt("SIRIUS_SAMPLES")) updated.render.samples_per_pixel = *value;
        if (auto value = GetEnvInt("SIRIUS_TILE_SIZE")) updated.render.tile_size = *value;
        if (auto value = GetEnvInt("SIRIUS_THREADS")) updated.render.thread_count = *value;
        if (auto value = GetEnv("SIRIUS_OUTPUT")) updated.render.output_path = *value;

        const auto metric_override = GetEnv("SIRIUS_METRIC");
        const auto mass_override = GetEnvDouble("SIRIUS_MASS");
        if (metric_override) updated.metric.name = *metric_override;
        if (mass_override) {
            updated.metric.mass = *mass_override;
        } else if (metric_override) {
            ApplyCompatibleMassDefault(updated.metric);
        }
        if (auto value = GetEnvDouble("SIRIUS_SPIN")) updated.metric.spin = *value;
        if (auto value = GetEnvDouble("SIRIUS_CHARGE")) updated.metric.charge = *value;

        if (auto value = GetEnvDouble("SIRIUS_DISTANCE")) updated.observer.distance = *value;
        if (auto value = GetEnvDouble("SIRIUS_INCLINATION")) {
            updated.observer.inclination = *value;
        }
        if (auto value = GetEnvDouble("SIRIUS_AZIMUTH")) updated.observer.azimuth = *value;
        if (auto value = GetEnvDouble("SIRIUS_FOV")) updated.observer.fov = *value;

        if (auto value = GetEnvBool("SIRIUS_BLOOM")) {
            updated.postprocess.enable_bloom = *value;
        }
        if (auto value = GetEnvDouble("SIRIUS_EXPOSURE")) {
            updated.postprocess.exposure = static_cast<float>(*value);
        }

        if (auto value = GetEnv("SIRIUS_BACKEND")) updated.backend.preferred = *value;
        if (auto value = GetEnv("SIRIUS_COLOR_MODE")) updated.color_mode = *value;
        if (auto value = GetEnvInt("SIRIUS_CUDA_DEVICE")) updated.backend.cuda_device = *value;
    } catch (const std::exception& error) {
        return base::Fail(base::ErrorDomain::kConfiguration, "apply environment overrides",
                          error.what());
    }
    config = std::move(updated);
    return {};
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

    if (config.render.samples_per_pixel < 1 || config.render.samples_per_pixel > 4096) {
        errors.push_back("render.samples_per_pixel must be between 1 and 4096");
    }

    if (config.render.tile_size < 8 || config.render.tile_size > 256) {
        errors.push_back("render.tile_size must be between 8 and 256");
    }
    if ((config.render.tile_size & (config.render.tile_size - 1)) != 0) {
        errors.push_back("render.tile_size should be a power of 2 for GPU efficiency");
    }
    if (config.render.thread_count < 0 || config.render.thread_count > 1024) {
        errors.push_back("render.thread_count must be between 0 (automatic) and 1024");
    }
    if (config.render.output_path.empty() ||
        config.render.output_path.find('\0') != std::string::npos) {
        errors.push_back("render.output_path must not be empty");
    } else {
        std::string extension = fs::path(config.render.output_path).extension().string();
        std::transform(extension.begin(), extension.end(), extension.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (extension != ".ppm" && extension != ".png" && extension != ".exr") {
            errors.push_back("render.output_path extension must be .ppm, .png, or .exr");
        }
    }

    // --- Metric validation (single authority: the identity registry) --------
    auto metric_id = core::ParseMetricName(config.metric.name);
    if (!metric_id.has_value()) {
        errors.push_back("metric.name '" + config.metric.name +
                         "' is not a known metric; accepted names: " + core::KnownMetricNames());
    }
    if (metric_id.has_value()) {
        if (const auto issue =
                core::MetricParameterIssue(*metric_id, config.metric.spin, config.metric.charge,
                                           config.metric.cosmological_constant);
            issue.has_value()) {
            errors.emplace_back(*issue);
        }
    }

    constexpr double kMinMass = 0.1;
    constexpr double kMaxMass = 100.0;
    const bool uses_mass = metric_id.has_value() && core::MetricUsesMass(*metric_id);
    if (finite(config.metric.mass, "metric.mass")) {
        if (metric_id.has_value() && !uses_mass && config.metric.mass != 0.0) {
            errors.push_back("metric.mass must be zero for metrics without a mass parameter");
        } else if ((!metric_id.has_value() || uses_mass) &&
                   (config.metric.mass < kMinMass || config.metric.mass > kMaxMass)) {
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
                  config.metric.temperature_model) == valid_temperature_models.end()) {
        errors.push_back(
            "metric.temperature_model must be one of: NovikovThorne, NT, ShakuraSunyaev, SS");
    }

    // Only the positive-Lambda spherical Kottler sector is publicly named.
    // Identity compatibility above rejects rotating/charged cosmological
    // requests; the horizon authority below rejects the Nariai boundary.
    constexpr double kMaxLambda = 0.1;
    if (finite(config.metric.cosmological_constant, "metric.lambda") &&
        (config.metric.cosmological_constant < 0.0 ||
         config.metric.cosmological_constant > kMaxLambda)) {
        errors.push_back("metric.lambda must be between 0 and 0.1");
    }
    if (metric_id.has_value() && std::isfinite(config.metric.mass) &&
        std::isfinite(config.metric.cosmological_constant)) {
        if (const auto issue = core::MetricHorizonIssue(*metric_id, config.metric.mass,
                                                        config.metric.cosmological_constant);
            issue.has_value()) {
            errors.emplace_back(*issue);
        }
    }
    if (finite(config.metric.disk_temperature, "metric.disk_temperature") &&
        (config.metric.disk_temperature < 100.0f || config.metric.disk_temperature > 1.0e8f)) {
        errors.push_back("metric.disk_temperature must be between 100 and 100000000 Kelvin");
    }
    if (config.disk_enabled && metric_id.has_value() &&
        core::DiskSupportFor(*metric_id) != core::DiskSupport::PageThorne) {
        errors.push_back(
            "diskEnabled must be false when the selected metric has no represented "
            "Page-Thorne accretion-disk emission model");
    }
    const bool default_temperature_model = config.metric.temperature_model == "NovikovThorne";
    if (!config.disk_enabled &&
        (!default_temperature_model ||
         config.metric.disk_temperature != core::kDefaultDiskTemperatureKelvin)) {
        errors.push_back("metric temperature model and disk temperature require diskEnabled=true");
    }
    if (!config.disk_enabled && !config.doppler_beaming) {
        errors.push_back("dopplerBeaming is a disk-emission control and requires diskEnabled=true");
    }
    if (config.ray_bundles && metric_id.has_value() && *metric_id != core::MetricId::Minkowski &&
        *metric_id != core::MetricId::Schwarzschild && *metric_id != core::MetricId::Kerr) {
        errors.push_back(
            "rayBundles require Minkowski, Schwarzschild, or Kerr: covariant curvature "
            "transport is not represented for the selected metric");
    }
    const bool throat_radius_finite = finite(config.metric.throat_radius, "metric.throat_radius");
    if (throat_radius_finite &&
        (config.metric.throat_radius < core::kMinMorrisThorneThroatRadius ||
         config.metric.throat_radius > core::kMaxMorrisThorneThroatRadius)) {
        errors.push_back("metric.throat_radius must be between 0.1 and 1000");
    }
    if (config.metric.wormhole_topology != "OneSheetCapture" &&
        config.metric.wormhole_topology != "TwoSheet") {
        errors.push_back("metric.wormhole_topology must be one of: OneSheetCapture, TwoSheet");
    } else if (config.metric.wormhole_topology == "TwoSheet" && metric_id.has_value() &&
               *metric_id == core::MetricId::MorrisThorne) {
        errors.push_back(
            "metric.wormhole_topology TwoSheet is not represented: Sirius currently renders "
            "one isotropic Ellis output sheet with the throat as a dark capture surface");
    }
    if (finite(config.metric.warp_velocity, "metric.warp_velocity") &&
        std::abs(config.metric.warp_velocity) > 10.0) {
        errors.push_back("metric.warp_velocity magnitude must be at most 10");
    }
    const bool bubble_radius_finite = finite(config.metric.bubble_radius, "metric.bubble_radius");
    if (bubble_radius_finite &&
        (config.metric.bubble_radius <= 0.0 || config.metric.bubble_radius > 1000.0)) {
        errors.push_back("metric.bubble_radius must be greater than 0 and at most 1000");
    }
    const bool bubble_sigma_finite = finite(config.metric.bubble_sigma, "metric.bubble_sigma");
    if (bubble_sigma_finite &&
        (config.metric.bubble_sigma <= 0.0 || config.metric.bubble_sigma > 1000.0)) {
        errors.push_back("metric.bubble_sigma must be greater than 0 and at most 1000");
    }
    if (metric_id.has_value()) {
        if (const auto issue = core::MetricSpecificParameterIssue(
                *metric_id, config.metric.throat_radius,
                config.metric.wormhole_topology == "OneSheetCapture", config.metric.warp_velocity,
                config.metric.bubble_radius, config.metric.bubble_sigma);
            issue.has_value()) {
            errors.emplace_back(*issue);
        }
        if (*metric_id == core::MetricId::Alcubierre && bubble_radius_finite &&
            bubble_sigma_finite && config.metric.bubble_radius > 0.0 &&
            config.metric.bubble_radius <= 1000.0 && config.metric.bubble_sigma > 0.0 &&
            config.metric.bubble_sigma <= 1000.0) {
            if (const auto issue = core::AlcubierreScaleIssue(config.metric.bubble_radius,
                                                              config.metric.bubble_sigma);
                issue.has_value()) {
                errors.emplace_back(*issue);
            }
        }
    }

    // --- Observer validation ------------------------------------------------
    // observer.distance is the coordinate radius r, not a dimensionless ratio.
    // The launch observer remains exterior to the governed central scene and
    // the finite trace budget scales with the same authority.
    const double distance_scale =
        metric_id.has_value() ? core::MetricSceneLengthScale(
                                    *metric_id, config.metric.mass, config.metric.throat_radius,
                                    config.metric.bubble_radius, config.metric.bubble_sigma)
                              : 1.0;
    const core::MetricObserverRadiusIssue observer_radius_issue =
        metric_id.has_value()
            ? core::MetricObserverRadiusIssueFor(
                  *metric_id, config.metric.mass, config.metric.cosmological_constant,
                  config.observer.distance, config.metric.throat_radius,
                  config.metric.bubble_radius, config.metric.bubble_sigma)
            : core::MetricObserverRadiusIssue::NaturalScale;

    if (finite(config.observer.distance, "observer.distance") && std::isfinite(distance_scale) &&
        distance_scale > 0.0 && observer_radius_issue != core::MetricObserverRadiusIssue::None) {
        if (observer_radius_issue == core::MetricObserverRadiusIssue::CosmologicalHorizon) {
            const double horizon = *core::MetricCosmologicalHorizonRadius(
                *metric_id, config.metric.mass, config.metric.cosmological_constant);
            errors.push_back(
                "positive-lambda observer.distance must remain at or below 0.99*r_c inside the "
                "cosmological horizon (r=" +
                std::to_string(config.observer.distance) + ", r_c=" + std::to_string(horizon) +
                ")");
        } else if (uses_mass) {
            errors.push_back(
                "observer.distance coordinate radius must satisfy 5*M <= r <= 1000*M "
                "(r=" +
                std::to_string(config.observer.distance) +
                ", M=" + std::to_string(config.metric.mass) + ")");
        } else if (metric_id == core::MetricId::MorrisThorne) {
            errors.push_back(
                "Morris-Thorne observer.distance must satisfy 5*b0 <= rho <= 1000*b0 "
                "for an exterior isotropic one-sheet observer (rho=" +
                std::to_string(config.observer.distance) +
                ", b0=" + std::to_string(config.metric.throat_radius) + ")");
        } else if (metric_id == core::MetricId::Alcubierre) {
            errors.push_back(
                "Alcubierre observer.distance must satisfy 5*L <= r <= 1000*L, "
                "L=max(R,1/sigma), for an exterior observer (r=" +
                std::to_string(config.observer.distance) + ", L=" + std::to_string(distance_scale) +
                ")");
        } else {
            errors.push_back(
                "observer.distance coordinate radius must be between 5 and 1000 geometric "
                "coordinate units (r=" +
                std::to_string(config.observer.distance) + ")");
        }
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
        finite(config.observer.camera_beta_forward, "observer.camera_beta_forward");
    const bool beta_up_finite = finite(config.observer.camera_beta_up, "observer.camera_beta_up");
    const bool beta_right_finite =
        finite(config.observer.camera_beta_right, "observer.camera_beta_right");
    const double beta_squared =
        config.observer.camera_beta_forward * config.observer.camera_beta_forward +
        config.observer.camera_beta_up * config.observer.camera_beta_up +
        config.observer.camera_beta_right * config.observer.camera_beta_right;
    if (beta_forward_finite && beta_up_finite && beta_right_finite &&
        (!std::isfinite(beta_squared) || beta_squared >= 1.0)) {
        errors.push_back("observer camera beta magnitude must be less than 1");
    }
    const auto lens_type = core::ParseLensType(config.observer.lens_model);
    if (!lens_type.has_value()) {
        errors.push_back("observer.lens_model must be one of: Pinhole, ThinLens, Fisheye");
    }
    if (finite(config.observer.focal_length, "observer.focal_length") &&
        (config.observer.focal_length <= 0.0f || config.observer.focal_length > 10000.0f)) {
        errors.push_back("observer.focal_length must be greater than 0 and at most 10000");
    }
    if (finite(config.observer.aperture, "observer.aperture") &&
        (config.observer.aperture <= 0.0f || config.observer.aperture > 128.0f)) {
        errors.push_back("observer.aperture must be greater than 0 and at most 128");
    }
    if (finite(config.observer.focus_distance, "observer.focus_distance") &&
        (config.observer.focus_distance <= 0.0f || config.observer.focus_distance > 1.0e6f)) {
        errors.push_back("observer.focus_distance must be greater than 0 and at most 1000000");
    }
    if (lens_type.has_value()) {
        if (const auto issue = core::LensSpecificParameterIssue(
                *lens_type, config.observer.focal_length, config.observer.aperture,
                config.observer.focus_distance);
            issue.has_value()) {
            errors.emplace_back(*issue);
        }
    }

    // --- Post-process validation --------------------------------------------
    if (finite(config.postprocess.exposure, "postprocess.exposure") &&
        (config.postprocess.exposure <= 0 || config.postprocess.exposure > 100)) {
        errors.push_back("postprocess.exposure must be between 0 and 100 stops");
    }
    if (finite(config.postprocess.bloom_intensity, "postprocess.bloom_intensity") &&
        (config.postprocess.bloom_intensity < 0 || config.postprocess.bloom_intensity > 5)) {
        errors.push_back("postprocess.bloom_intensity must be between 0 and 5");
    }
    if (finite(config.postprocess.bloom_threshold, "postprocess.bloom_threshold") &&
        (config.postprocess.bloom_threshold < 0 || config.postprocess.bloom_threshold > 100)) {
        errors.push_back("postprocess.bloom_threshold must be between 0 and 100");
    }
    if (!config.postprocess.enable_bloom &&
        (config.postprocess.bloom_intensity != core::kDefaultBloomIntensity ||
         config.postprocess.bloom_threshold != core::kDefaultBloomThreshold)) {
        errors.push_back("postprocess bloom intensity and threshold require enableBloom=true");
    }
    if (finite(config.postprocess.contrast, "postprocess.contrast") &&
        (config.postprocess.contrast < 0 || config.postprocess.contrast > 4)) {
        errors.push_back("postprocess.contrast must be between 0 and 4");
    }
    if (finite(config.postprocess.saturation, "postprocess.saturation") &&
        (config.postprocess.saturation < 0 || config.postprocess.saturation > 4)) {
        errors.push_back("postprocess.saturation must be between 0 and 4");
    }
    if (!core::ParseTonemapType(config.postprocess.tonemapper).has_value()) {
        errors.push_back("postprocess.tonemapper must be one of: " +
                         core::SupportedTonemapperNames());
    }

    // --- Volumetric and film validation ------------------------------------
    if (finite(config.volumetric.h_over_r, "volumetric.h_over_r") &&
        (config.volumetric.h_over_r < 0.01f || config.volumetric.h_over_r > 0.5f)) {
        errors.push_back("volumetric.h_over_r must be between 0.01 and 0.5");
    }
    if (finite(config.volumetric.h_power, "volumetric.h_power") &&
        (config.volumetric.h_power < -2 || config.volumetric.h_power > 4)) {
        errors.push_back("volumetric.h_power must be between -2 and 4");
    }
    if (finite(config.volumetric.tau_midplane, "volumetric.tau_midplane") &&
        (config.volumetric.tau_midplane < 0 || config.volumetric.tau_midplane > 1.0e6f)) {
        errors.push_back("volumetric.tau_midplane must be between 0 and 1000000");
    }
    if (config.volumetric.samples < 1 || config.volumetric.samples > 4096) {
        errors.push_back("volumetric.samples must be between 1 and 4096");
    }
    if (!config.volumetric.enabled &&
        (config.volumetric.h_over_r != core::kDefaultVolumetricHOverR ||
         config.volumetric.h_power != core::kDefaultVolumetricHPower ||
         config.volumetric.tau_midplane != core::kDefaultVolumetricTauMidplane ||
         config.volumetric.samples != core::kDefaultVolumetricSamples)) {
        errors.push_back("volumetric parameters require volumetric.enabled=true");
    }
    if (config.volumetric.enable_corona) {
        errors.push_back(
            "volumetric.enable_corona is not represented: frequency-dependent covariant "
            "Compton transfer is required");
    }
    if (config.volumetric.enable_turbulence && !config.volumetric.enabled) {
        errors.push_back("volumetric.enabled must be true when turbulence is enabled");
    }
    if (config.volumetric.enabled && !config.disk_enabled) {
        errors.push_back("diskEnabled must be true when volumetric disk rendering is enabled");
    }
    if (finite(config.motion_blur.shutter_time, "motionBlur.shutter_time") &&
        (config.motion_blur.shutter_time < 0.0f || config.motion_blur.shutter_time > 1000.0f)) {
        errors.push_back("motionBlur.shutter_time must be between 0 and 1000");
    }
    if (config.motion_blur.samples < 1 || config.motion_blur.samples > 4096) {
        errors.push_back("motionBlur.samples must be between 1 and 4096");
    }
    if (!config.motion_blur.enabled &&
        (config.motion_blur.shutter_time != core::kDefaultMotionBlurShutterTime ||
         config.motion_blur.samples != core::kDefaultMotionBlurSamples)) {
        errors.push_back("motionBlur parameters require motionBlur.enabled=true");
    }
    if (config.motion_blur.enabled && !config.disk_enabled) {
        errors.push_back("diskEnabled must be true when motion blur is enabled");
    }
    if (config.motion_blur.enabled) {
        errors.push_back(
            "motionBlur is not represented for the stationary axisymmetric disk: a covariant "
            "temporal emissivity model is required before shutter integration can be enabled");
    }

    static const std::vector<std::string> valid_color_modes = {
        "TrueColor", "TemperatureMap", "RedshiftMap", "Narrowband", "Polarisation"};
    if (std::find(valid_color_modes.begin(), valid_color_modes.end(), config.color_mode) ==
        valid_color_modes.end()) {
        errors.push_back(
            "colorMode must be one of: TrueColor, TemperatureMap, RedshiftMap, Narrowband, "
            "Polarisation");
    }
    if (config.color_mode == "Narrowband") {
        errors.push_back(
            "colorMode Narrowband is not represented: line emission requires ionisation, "
            "abundance, density, and frequency-dependent transfer rather than a temperature "
            "palette");
    }
    if (!config.disk_enabled && config.color_mode != "TrueColor") {
        errors.push_back("diagnostic color modes require diskEnabled=true");
    }
    if (config.color_mode == "Polarisation") {
        if (!config.disk_enabled) {
            errors.push_back("colorMode Polarisation requires diskEnabled");
        }
        if (config.volumetric.enabled) {
            errors.push_back(
                "colorMode Polarisation requires the thin disk; volumetric.enabled must be false");
        }
        if (config.motion_blur.enabled) {
            errors.push_back(
                "colorMode Polarisation is not represented with temporal disk motion blur");
        }
        if (metric_id.has_value() && *metric_id != core::MetricId::Schwarzschild &&
            *metric_id != core::MetricId::Kerr) {
            errors.push_back(
                "colorMode Polarisation is represented only for Schwarzschild and Kerr");
        }
    }

    static const std::vector<std::string> valid_film_presets = {"Interstellar", "SpaceOdyssey2001"};
    if (std::find(valid_film_presets.begin(), valid_film_presets.end(), config.film.preset) ==
        valid_film_presets.end()) {
        errors.push_back("film.preset must be one of: Interstellar, SpaceOdyssey2001");
    }
    if (config.film.grain_intensity.has_value() &&
        finite(*config.film.grain_intensity, "film.grain_intensity") &&
        (*config.film.grain_intensity < 0 || *config.film.grain_intensity > 1)) {
        errors.push_back("film.grain_intensity must be between 0 and 1");
    }
    if (config.film.halation_strength.has_value() &&
        finite(*config.film.halation_strength, "film.halation_strength") &&
        (*config.film.halation_strength < 0 || *config.film.halation_strength > 5)) {
        errors.push_back("film.halation_strength must be between 0 and 5");
    }
    if (config.film.vignette_strength.has_value() &&
        finite(*config.film.vignette_strength, "film.vignette_strength") &&
        (*config.film.vignette_strength < 0 || *config.film.vignette_strength > 2)) {
        errors.push_back("film.vignette_strength must be between 0 and 2");
    }
    const FilmFinishConfig default_film;
    if (!config.film.enabled &&
        (config.film.preset != default_film.preset || config.film.grain_intensity.has_value() ||
         config.film.halation_strength.has_value() || config.film.vignette_strength.has_value())) {
        errors.push_back("film-finish parameters require film.enabled=true");
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
    if (config.backend.enable_denoiser) {
        errors.push_back("backend.enable_denoiser is not implemented and must be false");
    }
    if (config.backend.cuda_device != 0) {
        errors.push_back("backend.cuda_device is retired with CUDA/OptiX and must be 0");
    }

    return errors;
}

std::string ConfigLoader::GenerateDefaultConfig() {
    SiriusConfig config = SiriusConfig::Defaults();
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
    if (const auto metric = source.find("metric");
        metric != source.end() && metric->contains("name") && !metric->contains("mass")) {
        ApplyCompatibleMassDefault(target.metric);
    }
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
