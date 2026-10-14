// CRCL004A.cpp - Config Command
// Component ID: CRCL004A (Cli/ConfigCommand)

#include "CRCL004A.h"
#include "CRCL005A.h"
#include "CRCF002A.h"
#include "CRPF001A.h"

#include <nlohmann/json.hpp>

#include <iostream>
#include <fstream>

namespace Sirius::Cli {

std::string ConfigCommand::usage() const {
    return R"(Usage: sirius config [subcommand]

Manage Sirius configuration files.

Subcommands:
  show          Display effective configuration
  validate      Validate a configuration file
  init          Create default configuration file
  paths         Show configuration search paths

Examples:
  sirius config show
  sirius config validate ./sirius.json
  sirius config init
  sirius config init --output ~/.config/sirius/config.json
)";
}

int ConfigCommand::execute(
    const std::vector<std::string>& args,
    const Configuration::GlobalOptions& globals,
    Configuration::SiriusConfig& config)
{
    // Default to show if no subcommand
    std::string subcommand = "show";
    if (!args.empty() && args[0][0] != '-') {
        subcommand = args[0];
    }

    if (subcommand == "show") {
        return showConfig(globals, config);
    } else if (subcommand == "validate") {
        return validateConfig(args, globals);
    } else if (subcommand == "init") {
        return initConfig(args, globals);
    } else if (subcommand == "paths") {
        return showPaths(globals);
    } else {
        Output::error("Unknown subcommand: " + subcommand);
        std::cout << usage() << std::endl;
        return 1;
    }
}

int ConfigCommand::showConfig(
    const Configuration::GlobalOptions& globals,
    const Configuration::SiriusConfig& config)
{
    if (globals.jsonOutput) {
        nlohmann::json j = config;

        // Add metadata
        auto loadedPath = Configuration::ConfigLoader::getLoadedConfigPath();
        if (loadedPath.has_value()) {
            j["_meta"]["source"] = loadedPath->string();
        } else {
            j["_meta"]["source"] = "defaults";
        }
        j["_meta"]["search_paths"] = nlohmann::json::array();
        for (const auto& path : Platform::PathResolver::configSearchPaths()) {
            j["_meta"]["search_paths"].push_back(path.string());
        }

        Output::printJson(j.dump(2));
    } else {
        // Show source info
        auto loadedPath = Configuration::ConfigLoader::getLoadedConfigPath();
        if (loadedPath.has_value()) {
            Output::info("Loaded from: " + loadedPath->string());
        } else {
            Output::info("Using default configuration (no config file found)");
        }
        std::cout << std::endl;

        Output::printConfig(config);
    }

    return 0;
}

int ConfigCommand::validateConfig(
    const std::vector<std::string>& args,
    const Configuration::GlobalOptions& globals)
{
    // Get file path from args
    std::string filePath;
    for (size_t i = 1; i < args.size(); ++i) {
        if (args[i][0] != '-') {
            filePath = args[i];
            break;
        }
    }

    if (filePath.empty()) {
        // Validate default config path
        auto defaultPath = Platform::PathResolver::findConfigFile();
        if (defaultPath.has_value()) {
            filePath = defaultPath->string();
        } else {
            Output::error("No configuration file specified and none found in search paths");
            return 1;
        }
    }

    // Check if file exists
    if (!std::filesystem::exists(filePath)) {
        if (globals.jsonOutput) {
            nlohmann::json j;
            j["valid"] = false;
            j["file"] = filePath;
            j["errors"] = nlohmann::json::array();
            j["errors"].push_back("File not found");
            Output::printJson(j.dump(2));
        } else {
            Output::error("File not found: " + filePath);
        }
        return 1;
    }

    // Try to load and parse
    try {
        std::ifstream file(filePath);
        nlohmann::json j;
        file >> j;

        // Try to parse as SiriusConfig
        Configuration::SiriusConfig config = j.get<Configuration::SiriusConfig>();

        // Validate
        auto errors = Configuration::ConfigLoader::validate(config);

        if (globals.jsonOutput) {
            nlohmann::json result;
            result["valid"] = errors.empty();
            result["file"] = filePath;
            result["errors"] = errors;
            Output::printJson(result.dump(2));
        } else {
            if (errors.empty()) {
                Output::success("Configuration is valid: " + filePath);
            } else {
                Output::error("Configuration has errors:");
                for (const auto& err : errors) {
                    std::cout << "  - " << err << std::endl;
                }
            }
        }

        return errors.empty() ? 0 : 1;

    } catch (const nlohmann::json::exception& e) {
        if (globals.jsonOutput) {
            nlohmann::json result;
            result["valid"] = false;
            result["file"] = filePath;
            result["errors"] = nlohmann::json::array();
            result["errors"].push_back(std::string("JSON parse error: ") + e.what());
            Output::printJson(result.dump(2));
        } else {
            Output::error("JSON parse error: " + std::string(e.what()));
        }
        return 1;
    }
}

int ConfigCommand::initConfig(
    const std::vector<std::string>& args,
    const Configuration::GlobalOptions& globals)
{
    // Determine output path
    std::string outputPath = "sirius.json";

    for (size_t i = 1; i < args.size(); ++i) {
        if ((args[i] == "-o" || args[i] == "--output") && i + 1 < args.size()) {
            outputPath = args[++i];
        }
    }

    // Check if file already exists
    if (std::filesystem::exists(outputPath)) {
        if (globals.jsonOutput) {
            nlohmann::json result;
            result["success"] = false;
            result["file"] = outputPath;
            result["error"] = "File already exists";
            Output::printJson(result.dump(2));
        } else {
            Output::error("File already exists: " + outputPath);
            Output::info("Use a different path or delete the existing file");
        }
        return 1;
    }

    // Generate and write default config
    Configuration::SiriusConfig defaultConfig = Configuration::SiriusConfig::defaults();

    if (Configuration::ConfigLoader::saveToFile(defaultConfig, outputPath)) {
        if (globals.jsonOutput) {
            nlohmann::json result;
            result["success"] = true;
            result["file"] = outputPath;
            Output::printJson(result.dump(2));
        } else {
            Output::success("Created configuration file: " + outputPath);
        }
        return 0;
    } else {
        if (globals.jsonOutput) {
            nlohmann::json result;
            result["success"] = false;
            result["file"] = outputPath;
            result["error"] = "Failed to write file";
            Output::printJson(result.dump(2));
        } else {
            Output::error("Failed to create configuration file");
        }
        return 1;
    }
}

int ConfigCommand::showPaths(const Configuration::GlobalOptions& globals) {
    auto paths = Platform::PathResolver::configSearchPaths();
    auto foundPath = Platform::PathResolver::findConfigFile();

    if (globals.jsonOutput) {
        nlohmann::json j;
        j["search_paths"] = nlohmann::json::array();
        for (const auto& path : paths) {
            j["search_paths"].push_back({
                {"path", path.string()},
                {"exists", std::filesystem::exists(path)}
            });
        }
        if (foundPath.has_value()) {
            j["active_config"] = foundPath->string();
        } else {
            j["active_config"] = nullptr;
        }
        j["user_config_dir"] = Platform::PathResolver::userConfigDirectory().string();
        j["system_config_dir"] = Platform::PathResolver::systemConfigDirectory().string();

        Output::printJson(j.dump(2));
    } else {
        std::cout << "Configuration search paths (in priority order):\n\n";
        for (const auto& path : paths) {
            bool exists = std::filesystem::exists(path);
            if (exists) {
                Output::success(path.string() + " (found)");
            } else {
                std::cout << "  " << path.string() << " (not found)\n";
            }
        }
        std::cout << "\n";
        std::cout << "User config directory:   " << Platform::PathResolver::userConfigDirectory().string() << "\n";
        std::cout << "System config directory: " << Platform::PathResolver::systemConfigDirectory().string() << "\n";
    }

    return 0;
}

} // namespace Sirius::Cli
