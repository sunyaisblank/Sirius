// CREP001A.cpp - Program Entry Point
// Component ID: CREP001A (Application/Entry Point)
//
// Supports batch rendering via FSM RenderSession (default) or legacy RenderJob.
// Use --interactive for deprecated real-time mode.

#include "CRAP001A.h"
#include "IRRJ001A.h"
#include "SNRS001A.h"
#include <iostream>
#include <exception>
#include <string>
#include <cstring>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void printUsage(const char* programName) {
    std::cout << "Sirius Static Renderer\n";
    std::cout << "======================\n\n";
    std::cout << "Usage:\n";
    std::cout << "  " << programName << " [options]\n\n";
    std::cout << "Options:\n";
    std::cout << "  --help, -h          Show this help message\n";
    std::cout << "  --output, -o PATH   Output file path (default: render.ppm)\n";
    std::cout << "  --width W           Image width (default: 1920)\n";
    std::cout << "  --height H          Image height (default: 1080)\n";
    std::cout << "  --samples N         Samples per pixel (default: 64)\n";
    std::cout << "  --tile-size T       Tile size in pixels (default: 64)\n";
    std::cout << "  --metric NAME       Metric to use (default: Schwarzschild)\n";
    std::cout << "  --distance R        Observer distance in M (default: 50)\n";
    std::cout << "  --inclination THETA Observer inclination in degrees (default: 90)\n";
    std::cout << "  --spin A            Black hole spin (0 to 1, default: 0)\n";
    std::cout << "  --fov FOV           Camera field of view in degrees (default: 60)\n";
    std::cout << "  --exposure E        Exposure value (default: 1.0)\n";
    std::cout << "  --no-bloom          Disable bloom post-processing\n";
    std::cout << "  --legacy            Use legacy RenderJob instead of FSM RenderSession\n";
    std::cout << "  --interactive       Run in interactive mode (deprecated)\n";
    std::cout << "\n";
}

int runSessionMode(const Sirius::SessionConfig& config) {
    std::cout << "================================================================================\n";
    std::cout << "                      SIRIUS FSM RENDER SESSION                                \n";
    std::cout << "================================================================================\n\n";
    
    std::cout << "Configuration:\n";
    std::cout << "  Resolution:  " << config.width << " x " << config.height << "\n";
    std::cout << "  Tile size:   " << config.tileSize << " px\n";
    std::cout << "  Samples:     " << config.samplesPerPixel << " spp\n";
    std::cout << "  Metric:      " << config.metricName << "\n";
    std::cout << "  Observer:    r=" << config.observerDistance << "M, "
              << "θ=" << (config.observerInclination * 180.0 / M_PI) << "°\n";
    std::cout << "  FOV:         " << config.cameraFOV << "°\n";
    std::cout << "  Output:      " << config.outputPath << "\n";
    std::cout << "\n";
    
    Sirius::RenderSession session;
    session.configure(config);
    
    session.setCompletionCallback([](Sirius::SessionState state, const std::string& message) {
        std::cout << "\n================================================================================\n";
        if (state == Sirius::SessionState::Complete) {
            std::cout << "                              RENDER COMPLETE                                   \n";
        } else if (state == Sirius::SessionState::Cancelled) {
            std::cout << "                              RENDER CANCELLED                                  \n";
        } else {
            std::cout << "                               RENDER FAILED                                    \n";
            std::cerr << "Error: " << message << "\n";
        }
        std::cout << "================================================================================\n";
    });
    
    Sirius::SessionState result = session.execute();
    return (result == Sirius::SessionState::Complete) ? 0 : 1;
}

int runLegacyMode(const Sirius::RenderConfig& config) {
    std::cout << "================================================================================\n";
    std::cout << "                      SIRIUS LEGACY RENDER (RenderJob)                         \n";
    std::cout << "================================================================================\n\n";
    
    std::cout << "Configuration:\n";
    std::cout << "  Resolution: " << config.width << " x " << config.height << "\n";
    std::cout << "  Samples:    " << config.samplesPerPixel << " spp\n";
    std::cout << "  Metric:     " << config.metricName << "\n";
    std::cout << "  Output:     " << config.outputPath << "\n";
    std::cout << "\n";
    
    Sirius::RenderJob job(config);
    
    job.setCompletionCallback([](Sirius::JobStatus status, const std::string& message) {
        std::cout << "\n================================================================================\n";
        if (status == Sirius::JobStatus::Completed) {
            std::cout << "                              RENDER COMPLETE                                   \n";
        } else {
            std::cerr << "Render failed: " << message << "\n";
        }
        std::cout << "================================================================================\n";
    });
    
    job.execute();
    return (job.getStatus() == Sirius::JobStatus::Completed) ? 0 : 1;
}

int runInteractiveMode() {
    std::cout << "================================================================================\n";
    std::cout << "                      SIRIUS INTERACTIVE MODE (DEPRECATED)                      \n";
    std::cout << "================================================================================\n";
    std::cout << "\nWARNING: Interactive mode is deprecated and will be removed.\n";
    std::cout << "         Please use batch rendering for production work.\n\n";
    
    Application app;
    app.run();
    return 0;
}

int main(int argc, char* argv[]) {
    try {
        Sirius::SessionConfig sessionConfig;
        Sirius::RenderConfig legacyConfig;
        bool interactiveMode = false;
        bool legacyMode = false;
        
        // Parse command line arguments
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            
            if (arg == "--help" || arg == "-h") {
                printUsage(argv[0]);
                return 0;
            }
            else if (arg == "--interactive") {
                interactiveMode = true;
            }
            else if (arg == "--legacy") {
                legacyMode = true;
            }
            else if ((arg == "--output" || arg == "-o") && i + 1 < argc) {
                std::string path = argv[++i];
                sessionConfig.outputPath = path;
                legacyConfig.outputPath = path;
            }
            else if (arg == "--width" && i + 1 < argc) {
                int w = std::stoi(argv[++i]);
                sessionConfig.width = w;
                legacyConfig.width = w;
            }
            else if (arg == "--height" && i + 1 < argc) {
                int h = std::stoi(argv[++i]);
                sessionConfig.height = h;
                legacyConfig.height = h;
            }
            else if (arg == "--samples" && i + 1 < argc) {
                int s = std::stoi(argv[++i]);
                sessionConfig.samplesPerPixel = s;
                legacyConfig.samplesPerPixel = s;
            }
            else if (arg == "--tile-size" && i + 1 < argc) {
                sessionConfig.tileSize = std::stoi(argv[++i]);
            }
            else if (arg == "--metric" && i + 1 < argc) {
                std::string m = argv[++i];
                sessionConfig.metricName = m;
                legacyConfig.metricName = m;
            }
            else if (arg == "--distance" && i + 1 < argc) {
                double d = std::stod(argv[++i]);
                sessionConfig.observerDistance = d;
                legacyConfig.observerPosition[1] = d;
            }
            else if (arg == "--inclination" && i + 1 < argc) {
                double degrees = std::stod(argv[++i]);
                double radians = degrees * M_PI / 180.0;
                sessionConfig.observerInclination = radians;
                legacyConfig.observerPosition[2] = radians;
            }
            else if (arg == "--spin" && i + 1 < argc) {
                double a = std::stod(argv[++i]);
                sessionConfig.blackHoleSpin = a;
                legacyConfig.a = a;
            }
            else if (arg == "--fov" && i + 1 < argc) {
                sessionConfig.cameraFOV = std::stof(argv[++i]);
            }
            else if (arg == "--exposure" && i + 1 < argc) {
                sessionConfig.exposure = std::stof(argv[++i]);
            }
            else if (arg == "--no-bloom") {
                sessionConfig.enableBloom = false;
            }
            else {
                std::cerr << "Unknown option: " << arg << "\n";
                printUsage(argv[0]);
                return 1;
            }
        }
        
        if (interactiveMode) {
            return runInteractiveMode();
        } else if (legacyMode) {
            return runLegacyMode(legacyConfig);
        } else {
            return runSessionMode(sessionConfig);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Unknown fatal error occurred." << std::endl;
        return 1;
    }
}