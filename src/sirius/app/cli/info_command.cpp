// The info command. Ported from CRCL003A.cpp.
//
// OptiX and the CUDA runtime are retired: system and backend reporting now name
// the CPU tracer (always present) and, when the Vulkan backend is compiled in,
// enumerate Vulkan-visible devices through the backend seam. The metric table
// prints from the single-authority registry, unchanged.

#include "sirius/app/cli/info_command.h"

#include "sirius/app/cli/cli_output.h"
#include "sirius/app/config/config_loader.h"
#include "sirius/app/platform_paths.h"
#include "sirius/core/metrics/registry.h"

#ifdef SIRIUS_HAS_VULKAN_BACKEND
#include "sirius/backend/device.h"
#endif

#include <nlohmann/json.hpp>

#include <iostream>

namespace sirius::app {

namespace cli = cli_output;

std::string InfoCommand::Usage() const {
    return R"(Usage: sirius info [subcommand]

Display system and configuration information.

Subcommands:
  system        Display system capabilities (platform, backends)
  metrics       List available spacetime metrics
  config        Show current configuration
  backends      List available render backends

Examples:
  sirius info system
  sirius info metrics
  sirius info config
  sirius info system --json
)";
}

int InfoCommand::Execute(const std::vector<std::string>& args, const GlobalOptions& globals,
                         SiriusConfig& config) {
    std::string subcommand = "system";
    if (!args.empty() && args[0][0] != '-') {
        subcommand = args[0];
    }

    if (subcommand == "system") {
        return ShowSystem(globals);
    } else if (subcommand == "metrics") {
        return ShowMetrics(globals);
    } else if (subcommand == "config") {
        return ShowConfig(globals, config);
    } else if (subcommand == "backends") {
        return ShowBackends(globals);
    } else {
        cli::Error("Unknown subcommand: " + subcommand);
        std::cout << Usage() << std::endl;
        return 1;
    }
}

int InfoCommand::ShowSystem(const GlobalOptions& globals) {
    if (globals.json_output) {
        nlohmann::json j;

        j["platform"] = PlatformPaths::PlatformName();
        j["wsl2"] = PlatformPaths::IsWsl2();

        j["backends"]["cpu"] = {{"available", true}};

#ifdef SIRIUS_HAS_VULKAN_BACKEND
        auto devices = backend::EnumerateVulkanDevices();
        if (devices.has_value()) {
            j["backends"]["vulkan"]["compiled"] = true;
            j["backends"]["vulkan"]["device_count"] = devices->size();
            nlohmann::json device_array = nlohmann::json::array();
            for (const auto& d : *devices) {
                device_array.push_back({{"name", d.name},
                                        {"kind", backend::ToString(d.kind)},
                                        {"device_local_bytes", d.device_local_bytes},
                                        {"supports_fp64", d.supports_fp64}});
            }
            j["backends"]["vulkan"]["devices"] = device_array;
        } else {
            j["backends"]["vulkan"]["compiled"] = true;
            j["backends"]["vulkan"]["error"] = devices.error().Description();
        }
#else
        j["backends"]["vulkan"]["compiled"] = false;
#endif

        cli::PrintJson(j.dump(2));
    } else {
        std::vector<cli::TableRow> rows;

        rows.push_back({"Platform", PlatformPaths::PlatformName(), false});
        if (PlatformPaths::IsWsl2()) {
            rows.push_back({"Environment", "WSL2", false});
        }

        rows.push_back({"Available Backends", "", true});
        rows.push_back({"CPU", "Available (reference path)", false});

#ifdef SIRIUS_HAS_VULKAN_BACKEND
        auto devices = backend::EnumerateVulkanDevices();
        if (devices.has_value()) {
            if (devices->empty()) {
                rows.push_back({"Vulkan", "Compiled in (no devices found)", false});
            } else {
                for (const auto& d : *devices) {
                    rows.push_back(
                        {"Vulkan", d.name + " (" + backend::ToString(d.kind) + ")", false});
                }
            }
        } else {
            rows.push_back({"Vulkan", "Compiled in (enumeration declined)", false});
        }
#else
        rows.push_back({"Vulkan", "Not compiled", false});
#endif

        cli::PrintTable("System Information", rows);
    }

    return 0;
}

int InfoCommand::ShowMetrics(const GlobalOptions& globals) {
    // The identity registry is the single metric catalogue; this command has no
    // list of its own.
    const auto& registry = core::MetricRegistry();

    if (globals.json_output) {
        nlohmann::json j = nlohmann::json::array();
        for (const auto& m : registry) {
            nlohmann::json aliases = nlohmann::json::array();
            for (const char* alias : m.aliases) {
                if (alias != nullptr) aliases.push_back(alias);
            }
            j.push_back({{"name", m.canonical_name},
                         {"aliases", aliases},
                         {"parameters", m.parameters},
                         {"backends", {{"cpu", m.cpu_supported}, {"gpu", m.gpu_supported}}}});
        }
        cli::PrintJson(j.dump(2));
    } else {
        std::vector<cli::TableRow> rows;
        for (const auto& m : registry) {
            std::string backends =
                m.gpu_supported ? (m.cpu_supported ? "CPU, GPU" : "GPU") : "CPU only";
            rows.push_back(
                {m.canonical_name, std::string(m.parameters) + " [" + backends + "]", false});
        }
        cli::PrintTable("Available Metrics", rows);
    }

    return 0;
}

int InfoCommand::ShowConfig(const GlobalOptions& globals, const SiriusConfig& config) {
    if (globals.json_output) {
        nlohmann::json j = config;

        auto loaded_path = ConfigLoader::GetLoadedConfigPath();
        if (loaded_path.has_value()) {
            j["_source"] = loaded_path->string();
        } else {
            j["_source"] = "defaults";
        }

        cli::PrintJson(j.dump(2));
    } else {
        auto loaded_path = ConfigLoader::GetLoadedConfigPath();
        if (loaded_path.has_value()) {
            cli::Info("Configuration loaded from: " + loaded_path->string());
        } else {
            cli::Info("Using default configuration");
        }
        std::cout << std::endl;

        cli::PrintConfig(config);
    }

    return 0;
}

int InfoCommand::ShowBackends(const GlobalOptions& globals) {
    if (globals.json_output) {
        nlohmann::json j;

        j["cpu"] = {{"available", true},
                    {"description", "CPU software ray tracing"},
                    {"features", {"Multi-threaded", "Portable"}}};

#ifdef SIRIUS_HAS_VULKAN_BACKEND
        auto devices = backend::EnumerateVulkanDevices();
        j["vulkan"]["available"] = true;
        j["vulkan"]["description"] = "Vulkan compute backend (not yet wired to the render session)";
        j["vulkan"]["device_count"] = devices.has_value() ? devices->size() : 0;
#else
        j["vulkan"] = {{"available", false},
                       {"description", "Vulkan compute backend"},
                       {"reason", "Vulkan development files not found at build time"}};
#endif

        cli::PrintJson(j.dump(2));
    } else {
        std::vector<cli::TableRow> rows;

        rows.push_back({"CPU", "Available (multi-threaded reference path)", false});

#ifdef SIRIUS_HAS_VULKAN_BACKEND
        auto devices = backend::EnumerateVulkanDevices();
        std::size_t count = devices.has_value() ? devices->size() : 0;
        rows.push_back(
            {"Vulkan", "Compiled in, " + std::to_string(count) + " device(s) (not yet wired)",
             false});
#else
        rows.push_back({"Vulkan", "Not available (build without Vulkan)", false});
#endif

        cli::PrintTable("Render Backends", rows);
    }

    return 0;
}

}  // namespace sirius::app
