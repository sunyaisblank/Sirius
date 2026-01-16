// CRCL002A.cpp - Render Command
// Component ID: CRCL002A (Cli/RenderCommand)

#include "CRCL002A.h"
#include "CRCL005A.h"
#include "CRCF002A.h"
#include "SNRS001A.h"
#include "IRRJ001A.h"

#include <nlohmann/json.hpp>

#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Sirius::Cli {

std::string RenderCommand::usage() const {
    return R"(Usage: sirius render [options]

Render a black hole visualization using geodesic ray tracing.

Options:
  -o, --output <path>       Output file path (default: render.ppm)
  -w, --width <n>           Image width (default: 1920)
  -h, --height <n>          Image height (default: 1080)
  -s, --samples <n>         Samples per pixel (default: 64)
  -t, --tile-size <n>       Tile size in pixels (default: 64)
  -m, --metric <name>       Metric type: Schwarzschild, Kerr (default: Schwarzschild)
  -d, --distance <r>        Observer distance in M (default: 50)
  -i, --inclination <deg>   Observer inclination (default: 90)
  -a, --spin <a>            Black hole spin 0-1 (default: 0)
  --fov <deg>               Camera field of view (default: 60)
  --exposure <e>            Exposure value (default: 1.0)
  --no-bloom                Disable bloom post-processing
  --legacy                  Use legacy RenderJob instead of FSM

Examples:
  sirius render -o test.ppm -s 32
  sirius render -m Kerr -a 0.9 -s 256 -o kerr.ppm
  sirius render -w 3840 -h 2160 -s 512 --fov 90
)";
}

int RenderCommand::execute(
    const std::vector<std::string>& args,
    const Configuration::GlobalOptions& globals,
    Configuration::SiriusConfig& config)
{
    // Parse command arguments
    if (!parseArgs(args, globals, config)) {
        return 1;
    }

    // Validate configuration
    auto errors = Configuration::ConfigLoader::validate(config);
    if (!errors.empty()) {
        for (const auto& err : errors) {
            Output::error(err);
        }
        return 1;
    }

    // Print configuration if verbose or not JSON output
    if (!globals.jsonOutput) {
        printConfig(config, globals.verbose);
    }

    // Execute render
    if (globals.legacyMode) {
        return executeLegacy(config, globals);
    } else {
        return executeSession(config, globals);
    }
}

bool RenderCommand::parseArgs(
    const std::vector<std::string>& args,
    const Configuration::GlobalOptions& globals,
    Configuration::SiriusConfig& config)
{
    for (size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];

        try {
            if ((arg == "-o" || arg == "--output") && i + 1 < args.size()) {
                config.render.outputPath = args[++i];
            }
            else if ((arg == "-w" || arg == "--width") && i + 1 < args.size()) {
                config.render.width = std::stoi(args[++i]);
            }
            else if ((arg == "-h" || arg == "--height") && i + 1 < args.size()) {
                config.render.height = std::stoi(args[++i]);
            }
            else if ((arg == "-s" || arg == "--samples") && i + 1 < args.size()) {
                config.render.samplesPerPixel = std::stoi(args[++i]);
            }
            else if ((arg == "-t" || arg == "--tile-size") && i + 1 < args.size()) {
                config.render.tileSize = std::stoi(args[++i]);
            }
            else if ((arg == "-m" || arg == "--metric") && i + 1 < args.size()) {
                config.metric.name = args[++i];
            }
            else if ((arg == "-d" || arg == "--distance") && i + 1 < args.size()) {
                config.observer.distance = std::stod(args[++i]);
            }
            else if ((arg == "-i" || arg == "--inclination") && i + 1 < args.size()) {
                config.observer.inclination = std::stod(args[++i]);
            }
            else if ((arg == "-a" || arg == "--spin") && i + 1 < args.size()) {
                config.metric.spin = std::stod(args[++i]);
            }
            else if (arg == "--fov" && i + 1 < args.size()) {
                config.observer.fov = std::stod(args[++i]);
            }
            else if (arg == "--exposure" && i + 1 < args.size()) {
                config.postprocess.exposure = std::stof(args[++i]);
            }
            else if (arg == "--no-bloom") {
                config.postprocess.enableBloom = false;
            }
            else if (arg == "--no-gpu" || arg == "--cpu") {
                config.backend.preferred = "cpu";
            }
            else if (arg == "--legacy") {
                // Handled at global level, but also accept here
            }
            else if (arg.substr(0, 1) == "-") {
                Output::error("Unknown option: " + arg);
                return false;
            }
        } catch (const std::exception& e) {
            Output::error("Invalid value for " + arg + ": " + e.what());
            return false;
        }
    }

    return true;
}

void RenderCommand::printConfig(const Configuration::SiriusConfig& config, bool verbose) {
    Output::rule("Sirius Render");
    std::cout << std::endl;

    std::cout << "Configuration:" << std::endl;
    std::cout << "  Resolution:  " << config.render.width << " x " << config.render.height << std::endl;
    std::cout << "  Samples:     " << config.render.samplesPerPixel << " spp" << std::endl;
    std::cout << "  Tile size:   " << config.render.tileSize << " px" << std::endl;
    std::cout << "  Metric:      " << config.metric.name;
    if (config.metric.spin > 0) {
        std::cout << " (a=" << std::fixed << std::setprecision(2) << config.metric.spin << ")";
    }
    std::cout << std::endl;
    std::cout << "  Observer:    r=" << std::fixed << std::setprecision(1)
              << config.observer.distance << "M, θ="
              << config.observer.inclination << "°" << std::endl;
    std::cout << "  FOV:         " << config.observer.fov << "°" << std::endl;
    std::cout << "  Output:      " << config.render.outputPath << std::endl;
    std::cout << std::endl;

    if (verbose) {
        std::cout << "Post-processing:" << std::endl;
        std::cout << "  Bloom:       " << (config.postprocess.enableBloom ? "enabled" : "disabled") << std::endl;
        std::cout << "  Exposure:    " << std::fixed << std::setprecision(2) << config.postprocess.exposure << std::endl;
        std::cout << std::endl;
    }
}

int RenderCommand::executeSession(
    const Configuration::SiriusConfig& config,
    const Configuration::GlobalOptions& globals)
{
    // Convert to SessionConfig
    Sirius::SessionConfig sessionConfig;
    sessionConfig.width = config.render.width;
    sessionConfig.height = config.render.height;
    sessionConfig.samplesPerPixel = config.render.samplesPerPixel;
    sessionConfig.tileSize = config.render.tileSize;
    sessionConfig.outputPath = config.render.outputPath;
    sessionConfig.metricName = config.metric.name;
    sessionConfig.blackHoleMass = config.metric.mass;
    sessionConfig.blackHoleSpin = config.metric.spin;
    sessionConfig.observerDistance = config.observer.distance;
    sessionConfig.observerInclination = config.observer.inclination * M_PI / 180.0;
    sessionConfig.cameraFOV = static_cast<float>(config.observer.fov);
    sessionConfig.enableBloom = config.postprocess.enableBloom;
    sessionConfig.bloomIntensity = config.postprocess.bloomIntensity;
    sessionConfig.exposure = config.postprocess.exposure;

    // Backend selection: "cpu" disables GPU, "auto" or "optix" enables it
    sessionConfig.useGPU = (config.backend.preferred != "cpu");

    // Create and configure session
    Sirius::RenderSession session;
    session.configure(sessionConfig);

    // Set up progress display
    auto startTime = std::chrono::steady_clock::now();

    if (!globals.jsonOutput) {
        session.setProgressCallback([&](float progress, int completed, int total, double eta) {
            auto now = std::chrono::steady_clock::now();
            double elapsed = std::chrono::duration<double>(now - startTime).count();

            Output::ProgressState state;
            state.progress = progress;
            state.tilesCompleted = completed;
            state.tilesTotal = total;
            state.elapsedSeconds = elapsed;
            state.etaSeconds = eta;

            Output::printProgress(state);
        });
    }

    // Set up completion callback
    bool success = false;
    std::string errorMessage;

    session.setCompletionCallback([&](Sirius::SessionState state, const std::string& message) {
        if (!globals.jsonOutput) {
            Output::clearProgress();
        }

        if (state == Sirius::SessionState::Complete) {
            success = true;
        } else {
            success = false;
            errorMessage = message;
        }
    });

    // Execute
    Sirius::SessionState result = session.execute();

    // Report results
    if (globals.jsonOutput) {
        nlohmann::json j;
        j["success"] = (result == Sirius::SessionState::Complete);
        j["output"] = config.render.outputPath;
        j["state"] = Sirius::stateName(result);
        if (!errorMessage.empty()) {
            j["error"] = errorMessage;
        }
        Output::printJson(j.dump(2));
    } else {
        std::cout << std::endl;
        Output::rule();
        if (result == Sirius::SessionState::Complete) {
            Output::success("Render complete: " + config.render.outputPath);
        } else if (result == Sirius::SessionState::Cancelled) {
            Output::warning("Render cancelled");
        } else {
            Output::error("Render failed: " + errorMessage);
        }
    }

    return (result == Sirius::SessionState::Complete) ? 0 : 1;
}

int RenderCommand::executeLegacy(
    const Configuration::SiriusConfig& config,
    const Configuration::GlobalOptions& globals)
{
    if (!globals.jsonOutput) {
        Output::warning("Using legacy RenderJob mode");
    }

    // Convert to RenderConfig
    Sirius::RenderConfig renderConfig;
    renderConfig.width = config.render.width;
    renderConfig.height = config.render.height;
    renderConfig.samplesPerPixel = config.render.samplesPerPixel;
    renderConfig.outputPath = config.render.outputPath;
    renderConfig.metricName = config.metric.name;
    renderConfig.a = config.metric.spin;
    renderConfig.observerPosition[1] = config.observer.distance;
    renderConfig.observerPosition[2] = config.observer.inclination * M_PI / 180.0;

    // Create and execute job
    Sirius::RenderJob job(renderConfig);

    bool success = false;
    std::string errorMessage;

    job.setCompletionCallback([&](Sirius::JobStatus status, const std::string& message) {
        if (status == Sirius::JobStatus::Completed) {
            success = true;
        } else {
            success = false;
            errorMessage = message;
        }
    });

    job.execute();

    // Report results
    if (globals.jsonOutput) {
        nlohmann::json j;
        j["success"] = success;
        j["output"] = config.render.outputPath;
        j["mode"] = "legacy";
        if (!errorMessage.empty()) {
            j["error"] = errorMessage;
        }
        Output::printJson(j.dump(2));
    } else {
        Output::rule();
        if (success) {
            Output::success("Render complete: " + config.render.outputPath);
        } else {
            Output::error("Render failed: " + errorMessage);
        }
    }

    return success ? 0 : 1;
}

} // namespace Sirius::Cli
