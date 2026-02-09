// SREP001A.cpp - GPU-Accelerated Offline Renderer Entry Point
// Component ID: SREP001A (Session/EntryPoint)
// 
// Main entry point for the SiriusRender executable.
// Uses the Acceleration BackendManager to support multiple backends (OptiX, etc.)

#include <iostream>
#include <string>
#include <chrono>
#include <thread>
#include <vector>
#include <cstring>
#include <iomanip>
#include <cmath>

#include <PHCN001A.h>
using Sirius::Constants::Math::PI;

// Architecture Includes
#include "../Acceleration/Backend/ACBM001A.h"
#include "../Buffer/BFIO001A.h"
#include "PPOP001A.h"
#include "../Acceleration/OptiX/RDOP003A.h" // For LaunchParams struct definition if needed
#include "../Output/OUIB001A.h"  // ImageBuffer
#include "../Output/OUEW001A.h"  // EXR Writer
#include "../Output/RDFL001A.h"  // Film Pipeline


// Dependencies
#include "stb_image.h"

using namespace Sirius::Acceleration;
using namespace Sirius::Buffer;
// using namespace Sirius::PostProcess; // Does not exist


//==============================================================================
// Command Line Parsing
//==============================================================================

struct RenderOptions {
    LaunchConfig config;
    std::string outputPath = "render_output";
    std::string backgroundPath = "Sirius.Render/Texture/Starfield.png";
    bool useBackground = true;
    bool verbose = true;
    bool outputEXR = false;
    
    Sirius::PostProcessConfig ppSettings;
};


void printUsage() {
    std::cout << "SiriusRender [options]\n"
              << "  --output <file>      Output filename (default: render.ppm)\n"
              << "  --width <n>          Width (default: 1920)\n"
              << "  --height <n>         Height (default: 1080)\n"
              << "  --spin <a>           Black hole spin (default: 0.999)\n"
              << "  --distance <r>       Observer distance (default: 200)\n"
              << "  --inclination <deg>  Observer inclination angle (0-90, default: 90)\n"
              << "  --azimuth <deg>      Camera azimuth offset (default: 0)\n"
              << "  --fov <deg>          Field of view (default: 60)\n"
              << "  --samples <n>        Samples per pixel (default: 64)\n"
              << "  --no-bloom           Disable bloom\n"
              << "  --backend <name>     Force backend (optix/cuda)\n"
              << "\nCinematic Features:\n"
              << "  --turbulence <val>   Enable disk turbulence (amplitude 0-1)\n"
              << "  --starfield          Enable enhanced procedural starfield\n"
              << "  --film               Enable IMAX film simulation\n"
              << "  --film-grain <val>   Film grain intensity (0-1, default: 0.025)\n"
              << "  --film-halation <v>  Film halation strength (0-1, default: 0.15)\n";
}

RenderOptions parseArgs(int argc, char* argv[]) {
    RenderOptions opts;
    
    for(int i=1; i<argc; ++i) {
        std::string arg = argv[i];
        if(arg == "--output" && i+1<argc) opts.outputPath = argv[++i];
        else if(arg == "--width" && i+1<argc) opts.config.width = std::atoi(argv[++i]);
        else if(arg == "--height" && i+1<argc) opts.config.height = std::atoi(argv[++i]);
        else if(arg == "--spin" && i+1<argc) opts.config.blackHoleSpin = std::atof(argv[++i]);
        else if(arg == "--distance" && i+1<argc) opts.config.observerDistance = std::atof(argv[++i]);
        else if(arg == "--fov" && i+1<argc) opts.config.cameraFOV = std::atof(argv[++i]);

        else if(arg == "--samples" && i+1<argc) opts.config.samplesPerPixel = std::atoi(argv[++i]);
        else if(arg == "--inclination" && i+1<argc) {
            float deg = std::atof(argv[++i]);
            opts.config.observerInclination = deg * PI / 180.0f;  // Convert degrees to radians
        }
        else if(arg == "--azimuth" && i+1<argc) {
            float deg = std::atof(argv[++i]);
            opts.config.observerAzimuth = deg * PI / 180.0f;  // Convert degrees to radians
        }
        else if(arg == "--no-bloom") opts.ppSettings.enableBloom = false;

        // Cinematic features
        else if(arg == "--turbulence" && i+1<argc) {
            opts.config.enableTurbulence = true;
            opts.config.turbulenceAmplitude = std::atof(argv[++i]);
        }
        else if(arg == "--starfield") opts.config.enableStarfield = true;
        else if(arg == "--film") opts.config.enableFilm = true;
        else if(arg == "--film-grain" && i+1<argc) {
            opts.config.enableFilm = true;
            opts.config.filmGrainIntensity = std::atof(argv[++i]);
        }
        else if(arg == "--film-halation" && i+1<argc) {
            opts.config.enableFilm = true;
            opts.config.filmHalationStrength = std::atof(argv[++i]);
        }
        else if(arg == "--metric" && i+1<argc) {
            std::string m = argv[++i];
            if (m == "Minkowski") opts.config.metricType = (int)Sirius::MetricType::Minkowski;
            else if (m == "Schwarzschild") opts.config.metricType = (int)Sirius::MetricType::Schwarzschild;
            else if (m == "Kerr") opts.config.metricType = (int)Sirius::MetricType::Kerr;
            else if (m == "KerrSchild") opts.config.metricType = (int)Sirius::MetricType::KerrSchild;
            else if (m == "ReissnerNordstrom") opts.config.metricType = (int)Sirius::MetricType::ReissnerNordstrom;
            else if (m == "Alcubierre") opts.config.metricType = (int)Sirius::MetricType::Alcubierre;
            else if (m == "DeSitter") opts.config.metricType = (int)Sirius::MetricType::DeSitter;
            else if (m == "EllisDrainhole") opts.config.metricType = (int)Sirius::MetricType::EllisDrainhole;
        }
        else if(arg == "--help") { printUsage(); exit(0); }

    }
    
    // Default Metric if not specified: Kerr (2)
    // Default Family: KerrSchild (0)
    if (opts.config.metricType == 0) {
        // Did user specify Minkowski explicitly? Likely not differentiating 0 vs unset.
        // Let's default to Kerr for a demo unless they said --metric Minkowski
        // But 0 is Minkowski. 
        // We can check if args contained --metric. 
        // For simplicity, let's just set the default in initialization or here if 0.
        // Actually, let's force a default in declaration.
    }
    
    if (opts.outputPath.find(".exr") != std::string::npos) {
        opts.outputEXR = true;
    }
    
    return opts;
}

//==============================================================================
// Main
//==============================================================================
int main(int argc, char* argv[]) {
    RenderOptions opts = parseArgs(argc, argv);
    
    if(opts.verbose) {
        std::cout << "SiriusRender Enterprise Edition\n"
                  << "Backend: " << backendName(getBestBackend()) << "\n"
                  << "Resolution: " << opts.config.width << "x" << opts.config.height << "\n"
                  << "Samples: " << opts.config.samplesPerPixel << "\n";
    }
    
    // Initialize Accelerator
    auto& manager = BackendManager::instance();
    IAccelerator* accelerator = manager.getAccelerator(BackendType::OptiX);
    
    if (!accelerator) {
        std::cerr << "Fatal: No generic accelerator backend available.\n";
        return 1;
    }
    
    if (!accelerator->initialise(opts.config.width, opts.config.height)) {
        std::cerr << "Fatal: Failed to initialise backend: " << accelerator->getLastError() << "\n";
        return 1;
    }
    
    // Load Background
    if (opts.useBackground) {
        int w, h, c;
        unsigned char* data = stbi_load(opts.backgroundPath.c_str(), &w, &h, &c, 4);
        if (data) {
            accelerator->uploadBackground(data, w, h);
            stbi_image_free(data);
        } else if (opts.verbose) {
           // Try relative path
           std::string relPath = "../" + opts.backgroundPath;
           data = stbi_load(relPath.c_str(), &w, &h, &c, 4);
           if(data) {
               accelerator->uploadBackground(data, w, h);
               stbi_image_free(data);
           } else {
               std::cerr << "Warning: Could not load starfield: " << opts.backgroundPath << "\n";
           }
        }
    }
    
    // Render Loop
    std::cout << "Rendering..." << std::flush;
    auto start = std::chrono::high_resolution_clock::now();

    // Reset accumulation before starting fresh render
    accelerator->resetAccumulation();

    // Progressive accumulation: each launch adds 1 sample, backend accumulates
    int totalSamples = opts.config.samplesPerPixel;

    for (int s = 0; s < totalSamples; ++s) {
        accelerator->launch(opts.config);

        if (opts.verbose && (s % 10 == 0 || s == totalSamples - 1)) {
            float progress = 100.0f * (s + 1) / totalSamples;
            std::cout << "\rRendering: " << (s + 1) << "/" << totalSamples
                      << " samples (" << std::fixed << std::setprecision(1) << progress << "%)" << std::flush;
        }
    }
    
    accelerator->synchronise();
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "\nDone in " << diff.count() << "s\n";
    
    // Retrieve Buffer
    float* gpuBuffer = accelerator->getFrameBuffer();
    if (!gpuBuffer) {
         std::cerr << "Error: Null framebuffer\n";
         return 1;
    }

    // =========================================================================
    // Film Pipeline Post-Processing (Cinematic Features Phase 8)
    // =========================================================================
    if (opts.config.enableFilm) {
        std::cout << "Applying film simulation..." << std::flush;

        Sirius::FilmConfig filmConfig = Sirius::FilmConfig::Interstellar();
        filmConfig.grain_intensity = opts.config.filmGrainIntensity;
        filmConfig.halation_strength = opts.config.filmHalationStrength;
        filmConfig.halation_radius = opts.config.filmHalationRadius;
        filmConfig.vignette_strength = opts.config.filmVignetteStrength;

        // Cinematic grading: darker, higher contrast for Interstellar look
        filmConfig.exposure = -0.7f;       // Darken overall exposure
        filmConfig.contrast = 1.4f;        // Increase contrast for deep blacks
        filmConfig.saturation = 0.85f;     // Slightly desaturate for film look
        filmConfig.toe_strength = 0.6f;    // Stronger shadow compression
        filmConfig.shoulder_strength = 0.4f; // Softer highlight rolloff
        filmConfig.enabled = true;

        Sirius::FilmPipeline pipeline(filmConfig);
        pipeline.apply(gpuBuffer, opts.config.width, opts.config.height, 0);

        std::cout << " done.\n";
    }

    // Buffer IO
    if (opts.outputEXR) {
        // Convert float buffer to ImageBufferRGBA and write EXR
        sirius::render::ImageBufferRGBA buffer;
        buffer.allocate(opts.config.width, opts.config.height);
        std::memcpy(buffer.pixels.data(), gpuBuffer,
                    static_cast<size_t>(opts.config.width) * opts.config.height * 4 * sizeof(float));

        sirius::render::EXRMetadata meta;
        meta.blackHoleSpin = opts.config.blackHoleSpin;
        meta.observerDistance = opts.config.observerDistance;
        meta.samplesPerPixel = opts.config.samplesPerPixel;
        meta.renderTimeSeconds = diff.count();

        if (sirius::render::EXRWriter::writeEXR(opts.outputPath, buffer, meta)) {
            std::cout << "Wrote EXR: " << opts.outputPath << "\n";
        } else {
            std::cerr << "Failed to write EXR, falling back to PPM\n";
            BufferWriter::writePPM(opts.outputPath + ".ppm", gpuBuffer, opts.config.width, opts.config.height);
        }
    } else {
        BufferWriter::writePPM(opts.outputPath, gpuBuffer, opts.config.width, opts.config.height);
    }
    
    manager.release();
    return 0;
}
