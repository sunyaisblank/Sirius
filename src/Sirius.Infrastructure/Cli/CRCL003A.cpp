// CRCL003A.cpp - Info Command
// Component ID: CRCL003A (Cli/InfoCommand)

#include "CRCL003A.h"
#include "CRCL005A.h"
#include "CRCF002A.h"
#include "CRPF001A.h"

#include <nlohmann/json.hpp>

#include <iostream>
#include <sstream>

#ifdef SIRIUS_HAS_OPTIX
#include <cuda_runtime.h>
#endif

namespace Sirius::Cli {

std::string InfoCommand::usage() const {
    return R"(Usage: sirius info [subcommand]

Display system and configuration information.

Subcommands:
  system        Display system capabilities (GPU, CUDA, OptiX)
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

int InfoCommand::execute(
    const std::vector<std::string>& args,
    const Configuration::GlobalOptions& globals,
    Configuration::SiriusConfig& config)
{
    // Default to system info if no subcommand
    std::string subcommand = "system";
    if (!args.empty() && args[0][0] != '-') {
        subcommand = args[0];
    }

    if (subcommand == "system") {
        return showSystem(globals);
    } else if (subcommand == "metrics") {
        return showMetrics(globals);
    } else if (subcommand == "config") {
        return showConfig(globals, config);
    } else if (subcommand == "backends") {
        return showBackends(globals);
    } else {
        Output::error("Unknown subcommand: " + subcommand);
        std::cout << usage() << std::endl;
        return 1;
    }
}

int InfoCommand::showSystem(const Configuration::GlobalOptions& globals) {
    if (globals.jsonOutput) {
        nlohmann::json j;

        j["platform"] = Platform::PathResolver::platformName();
        j["wsl2"] = Platform::PathResolver::isWSL2();

#ifdef SIRIUS_HAS_OPTIX
        int cudaVersion = 0;
        cudaRuntimeGetVersion(&cudaVersion);
        int major = cudaVersion / 1000;
        int minor = (cudaVersion % 1000) / 10;

        j["cuda"]["available"] = true;
        j["cuda"]["version"] = std::to_string(major) + "." + std::to_string(minor);

        int deviceCount = 0;
        cudaGetDeviceCount(&deviceCount);
        j["gpu"]["count"] = deviceCount;

        if (deviceCount > 0) {
            cudaDeviceProp prop;
            cudaGetDeviceProperties(&prop, 0);
            j["gpu"]["name"] = prop.name;
            j["gpu"]["vram_mb"] = prop.totalGlobalMem / (1024 * 1024);
            j["gpu"]["compute_capability"] = std::to_string(prop.major) + "." + std::to_string(prop.minor);
        }

        j["optix"]["available"] = true;
#else
        j["cuda"]["available"] = false;
        j["optix"]["available"] = false;
#endif

        j["backends"]["optix"] = {
#ifdef SIRIUS_HAS_OPTIX
            {"available", true},
            {"preferred", true}
#else
            {"available", false},
            {"preferred", false}
#endif
        };
        j["backends"]["cpu"] = {
            {"available", true},
#ifdef SIRIUS_HAS_OPTIX
            {"preferred", false}
#else
            {"preferred", true}
#endif
        };

        Output::printJson(j.dump(2));
    } else {
        Output::printSystemInfo();
    }

    return 0;
}

int InfoCommand::showMetrics(const Configuration::GlobalOptions& globals) {
    struct MetricInfo {
        std::string name;
        std::string description;
        bool hasSpinParam;
        bool hasChargeParam;
    };

    std::vector<MetricInfo> metrics = {
        {"Minkowski", "Flat spacetime (special relativity)", false, false},
        {"Schwarzschild", "Non-rotating black hole", false, false},
        {"Kerr", "Rotating black hole", true, false},
        {"Reissner-Nordstrom", "Charged non-rotating black hole", false, true},
        {"Kerr-Newman", "Rotating charged black hole", true, true}
    };

    if (globals.jsonOutput) {
        nlohmann::json j = nlohmann::json::array();
        for (const auto& m : metrics) {
            j.push_back({
                {"name", m.name},
                {"description", m.description},
                {"parameters", {
                    {"spin", m.hasSpinParam},
                    {"charge", m.hasChargeParam}
                }}
            });
        }
        Output::printJson(j.dump(2));
    } else {
        std::vector<Output::TableRow> rows;
        for (const auto& m : metrics) {
            std::string params;
            if (m.hasSpinParam && m.hasChargeParam) {
                params = "M, a, Q";
            } else if (m.hasSpinParam) {
                params = "M, a";
            } else if (m.hasChargeParam) {
                params = "M, Q";
            } else {
                params = "M";
            }
            rows.push_back({m.name, m.description + " [" + params + "]"});
        }
        Output::printTable("Available Metrics", rows);
    }

    return 0;
}

int InfoCommand::showConfig(
    const Configuration::GlobalOptions& globals,
    const Configuration::SiriusConfig& config)
{
    if (globals.jsonOutput) {
        nlohmann::json j = config;

        // Add config source info
        auto loadedPath = Configuration::ConfigLoader::getLoadedConfigPath();
        if (loadedPath.has_value()) {
            j["_source"] = loadedPath->string();
        } else {
            j["_source"] = "defaults";
        }

        Output::printJson(j.dump(2));
    } else {
        // Show config source
        auto loadedPath = Configuration::ConfigLoader::getLoadedConfigPath();
        if (loadedPath.has_value()) {
            Output::info("Configuration loaded from: " + loadedPath->string());
        } else {
            Output::info("Using default configuration");
        }
        std::cout << std::endl;

        Output::printConfig(config);
    }

    return 0;
}

int InfoCommand::showBackends(const Configuration::GlobalOptions& globals) {
    if (globals.jsonOutput) {
        nlohmann::json j;

#ifdef SIRIUS_HAS_OPTIX
        j["optix"] = {
            {"available", true},
            {"description", "NVIDIA OptiX RTX ray tracing"},
            {"features", {"RTX acceleration", "Hardware BVH", "AI denoising"}}
        };
#else
        j["optix"] = {
            {"available", false},
            {"description", "NVIDIA OptiX RTX ray tracing"},
            {"reason", "CUDA/OptiX not found at build time"}
        };
#endif

        j["cpu"] = {
            {"available", true},
            {"description", "CPU software ray tracing"},
            {"features", {"Multi-threaded", "Portable"}}
        };

        Output::printJson(j.dump(2));
    } else {
        std::vector<Output::TableRow> rows;

#ifdef SIRIUS_HAS_OPTIX
        rows.push_back({"OptiX", "Available (RTX acceleration, AI denoising)"});
#else
        rows.push_back({"OptiX", "Not available (build without CUDA)"});
#endif
        rows.push_back({"CPU", "Available (multi-threaded fallback)"});

        Output::printTable("Render Backends", rows);
    }

    return 0;
}

} // namespace Sirius::Cli
