// CRCF002A.h - Configuration Loader
// Component ID: CRCF002A (Configuration/Loader)
//
// Loads and merges configuration from multiple sources:
// 1. Built-in defaults
// 2. System config file
// 3. User config file
// 4. Local config file
// 5. Environment variables (SIRIUS_* prefix)
// 6. Command-line arguments (highest priority)

#pragma once

#include "CRCF003A.h"
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace Sirius::Configuration {

namespace fs = std::filesystem;

/// @brief Configuration loading and merging
class ConfigLoader {
public:
    /// @brief Load configuration from all sources (merged)
    /// @param overridePath Optional path to override config file search
    /// @return Merged configuration
    static SiriusConfig load(const std::optional<std::string>& overridePath = std::nullopt);

    /// @brief Load configuration from a specific JSON file
    /// @param path Path to JSON configuration file
    /// @return Loaded configuration, or defaults on error
    static SiriusConfig loadFromFile(const fs::path& path);

    /// @brief Save configuration to a JSON file
    /// @param config Configuration to save
    /// @param path Destination path
    /// @return true if successful
    static bool saveToFile(const SiriusConfig& config, const fs::path& path);

    /// @brief Apply environment variable overrides
    /// @param config Configuration to modify in-place
    static void applyEnvironmentOverrides(SiriusConfig& config);

    /// @brief Get the path to the loaded config file (if any)
    /// @return Path to config file used, or nullopt if using defaults
    static std::optional<fs::path> getLoadedConfigPath();

    /// @brief Validate configuration values
    /// @param config Configuration to validate
    /// @return Vector of validation error messages (empty if valid)
    static std::vector<std::string> validate(const SiriusConfig& config);

    /// @brief Generate default configuration file content (pretty-printed JSON)
    /// @return JSON string for default configuration
    static std::string generateDefaultConfig();

private:
    /// @brief Merge source config into target (source takes precedence for non-default values)
    static void mergeConfig(SiriusConfig& target, const nlohmann::json& source);

    /// @brief Get environment variable as string
    static std::optional<std::string> getEnv(const std::string& name);

    /// @brief Get environment variable as integer
    static std::optional<int> getEnvInt(const std::string& name);

    /// @brief Get environment variable as double
    static std::optional<double> getEnvDouble(const std::string& name);

    /// @brief Get environment variable as boolean
    static std::optional<bool> getEnvBool(const std::string& name);

    /// @brief Last loaded config path
    static std::optional<fs::path> s_loadedPath;
};

} // namespace Sirius::Configuration
