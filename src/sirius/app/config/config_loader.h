#pragma once

// Loads and merges configuration from defaults, config files, environment
// variables, and (via the CLI) command-line flags. Ported from CRCF002A.h;
// methods take the new-tree PascalCase spelling.
//
// Layering, lowest priority first: built-in defaults, system/user/local config
// file, SIRIUS_* environment variables, command-line arguments.

#include "sirius/app/config/config_schema.h"

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace sirius::app {

namespace fs = std::filesystem;

// Configuration loading, merging, and validation.
class ConfigLoader {
  public:
    // Load merged configuration; override_path forces a specific config file.
    // Missing, malformed, unknown, or invalid values throw: configuration never
    // degrades to defaults after the operator supplied an input.
    static SiriusConfig Load(const std::optional<std::string>& override_path = std::nullopt);

    // Load and validate a specific JSON file; throws on every read/parse error.
    static SiriusConfig LoadFromFile(const fs::path& path);

    // Serialise a configuration to a JSON file. Returns false on IO failure.
    static bool SaveToFile(const SiriusConfig& config, const fs::path& path);

    // Apply SIRIUS_* environment overrides in place; malformed set variables
    // throw rather than being treated as unset.
    static void ApplyEnvironmentOverrides(SiriusConfig& config);

    // Path of the config file the last Load() used, or nullopt for defaults.
    static std::optional<fs::path> GetLoadedConfigPath();

    // Validate configuration; returns one message per violated constraint.
    static std::vector<std::string> Validate(const SiriusConfig& config);

    // Pretty-printed default configuration JSON.
    static std::string GenerateDefaultConfig();

  private:
    // Merge a partial JSON document over target (source wins per leaf key).
    static void MergeConfig(SiriusConfig& target, const nlohmann::json& source);

    static std::optional<std::string> GetEnv(const std::string& name);
    static std::optional<int> GetEnvInt(const std::string& name);
    static std::optional<double> GetEnvDouble(const std::string& name);
    static std::optional<bool> GetEnvBool(const std::string& name);

    static std::optional<fs::path> loaded_path_;
};

}  // namespace sirius::app
