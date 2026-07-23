// The config command. Ported from CRCL004A.cpp.

#include "sirius/app/cli/config_command.h"

#include "sirius/app/cli/cli_output.h"
#include "sirius/app/config/config_loader.h"
#include "sirius/app/platform_paths.h"

#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>

namespace sirius::app {

namespace cli = cli_output;
namespace fs = std::filesystem;

std::string ConfigCommand::Usage() const {
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

int ConfigCommand::Execute(const std::vector<std::string>& args, const GlobalOptions& globals,
                           SiriusConfig& config) {
    std::string subcommand = "show";
    if (!args.empty() && args[0][0] != '-') {
        subcommand = args[0];
    }

    if (subcommand == "show") {
        return ShowConfig(globals, config);
    } else if (subcommand == "validate") {
        return ValidateConfig(args, globals);
    } else if (subcommand == "init") {
        return InitConfig(args, globals);
    } else if (subcommand == "paths") {
        return ShowPaths(globals);
    } else {
        cli::Error("Unknown subcommand: " + subcommand);
        std::cout << Usage() << std::endl;
        return 1;
    }
}

int ConfigCommand::ShowConfig(const GlobalOptions& globals, const SiriusConfig& config) {
    if (globals.json_output) {
        nlohmann::json j = config;

        auto loaded_path = ConfigLoader::GetLoadedConfigPath();
        if (loaded_path.has_value()) {
            j["_meta"]["source"] = loaded_path->string();
        } else {
            j["_meta"]["source"] = "defaults";
        }
        j["_meta"]["search_paths"] = nlohmann::json::array();
        for (const auto& path : PlatformPaths::ConfigSearchPaths()) {
            j["_meta"]["search_paths"].push_back(path.string());
        }

        cli::PrintJson(j.dump(2));
    } else {
        auto loaded_path = ConfigLoader::GetLoadedConfigPath();
        if (loaded_path.has_value()) {
            cli::Info("Loaded from: " + loaded_path->string());
        } else {
            cli::Info("Using default configuration (no config file found)");
        }
        std::cout << std::endl;

        cli::PrintConfig(config);
    }

    return 0;
}

int ConfigCommand::ValidateConfig(const std::vector<std::string>& args,
                                  const GlobalOptions& globals) {
    std::string file_path;
    for (size_t i = 1; i < args.size(); ++i) {
        if (args[i][0] != '-') {
            file_path = args[i];
            break;
        }
    }

    if (file_path.empty()) {
        auto default_path = PlatformPaths::FindConfigFile();
        if (default_path.has_value()) {
            file_path = default_path->string();
        } else {
            cli::Error("No configuration file specified and none found in search paths");
            return 1;
        }
    }

    if (!fs::exists(file_path)) {
        if (globals.json_output) {
            nlohmann::json j;
            j["valid"] = false;
            j["file"] = file_path;
            j["errors"] = nlohmann::json::array();
            j["errors"].push_back("File not found");
            cli::PrintJson(j.dump(2));
        } else {
            cli::Error("File not found: " + file_path);
        }
        return 1;
    }

    try {
        std::ifstream file(file_path);
        nlohmann::json j;
        file >> j;

        SiriusConfig config = j.get<SiriusConfig>();

        auto errors = ConfigLoader::Validate(config);

        if (globals.json_output) {
            nlohmann::json result;
            result["valid"] = errors.empty();
            result["file"] = file_path;
            result["errors"] = errors;
            cli::PrintJson(result.dump(2));
        } else {
            if (errors.empty()) {
                cli::Success("Configuration is valid: " + file_path);
            } else {
                cli::Error("Configuration has errors:");
                for (const auto& err : errors) {
                    std::cout << "  - " << err << std::endl;
                }
            }
        }

        return errors.empty() ? 0 : 1;

    } catch (const nlohmann::json::exception& e) {
        if (globals.json_output) {
            nlohmann::json result;
            result["valid"] = false;
            result["file"] = file_path;
            result["errors"] = nlohmann::json::array();
            result["errors"].push_back(std::string("JSON parse error: ") + e.what());
            cli::PrintJson(result.dump(2));
        } else {
            cli::Error("JSON parse error: " + std::string(e.what()));
        }
        return 1;
    }
}

int ConfigCommand::InitConfig(const std::vector<std::string>& args, const GlobalOptions& globals) {
    std::string output_path = "sirius.json";

    for (size_t i = 1; i < args.size(); ++i) {
        if ((args[i] == "-o" || args[i] == "--output") && i + 1 < args.size()) {
            output_path = args[++i];
        }
    }

    if (fs::exists(output_path)) {
        if (globals.json_output) {
            nlohmann::json result;
            result["success"] = false;
            result["file"] = output_path;
            result["error"] = "File already exists";
            cli::PrintJson(result.dump(2));
        } else {
            cli::Error("File already exists: " + output_path);
            cli::Info("Use a different path or delete the existing file");
        }
        return 1;
    }

    SiriusConfig default_config = SiriusConfig::defaults();

    if (ConfigLoader::SaveToFile(default_config, output_path)) {
        if (globals.json_output) {
            nlohmann::json result;
            result["success"] = true;
            result["file"] = output_path;
            cli::PrintJson(result.dump(2));
        } else {
            cli::Success("Created configuration file: " + output_path);
        }
        return 0;
    } else {
        if (globals.json_output) {
            nlohmann::json result;
            result["success"] = false;
            result["file"] = output_path;
            result["error"] = "Failed to write file";
            cli::PrintJson(result.dump(2));
        } else {
            cli::Error("Failed to create configuration file");
        }
        return 1;
    }
}

int ConfigCommand::ShowPaths(const GlobalOptions& globals) {
    auto paths = PlatformPaths::ConfigSearchPaths();
    auto found_path = PlatformPaths::FindConfigFile();

    if (globals.json_output) {
        nlohmann::json j;
        j["search_paths"] = nlohmann::json::array();
        for (const auto& path : paths) {
            j["search_paths"].push_back({{"path", path.string()}, {"exists", fs::exists(path)}});
        }
        if (found_path.has_value()) {
            j["active_config"] = found_path->string();
        } else {
            j["active_config"] = nullptr;
        }
        j["user_config_dir"] = PlatformPaths::UserConfigDirectory().string();
        j["system_config_dir"] = PlatformPaths::SystemConfigDirectory().string();

        cli::PrintJson(j.dump(2));
    } else {
        std::cout << "Configuration search paths (in priority order):\n\n";
        for (const auto& path : paths) {
            bool exists = fs::exists(path);
            if (exists) {
                cli::Success(path.string() + " (found)");
            } else {
                std::cout << "  " << path.string() << " (not found)\n";
            }
        }
        std::cout << "\n";
        std::cout << "User config directory:   " << PlatformPaths::UserConfigDirectory().string()
                  << "\n";
        std::cout << "System config directory: " << PlatformPaths::SystemConfigDirectory().string()
                  << "\n";
    }

    return 0;
}

}  // namespace sirius::app
