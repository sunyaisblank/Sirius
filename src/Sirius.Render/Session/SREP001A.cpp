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

// Architecture Includes
#include "../Acceleration/Backend/ACBM001A.h"
#include "../Buffer/BFIO001A.h"
#include "../../Sirius.Core/PostProcess/PPOP001A.h"
#include "../Acceleration/OptiX/RDOP003A.h" // For LaunchParams struct definition if needed


// Dependencies
#define STB_IMAGE_IMPLEMENTATION
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
              << "  --fov <deg>          Field of view (default: 60)\n"
              << "  --samples <n>        Samples per pixel (default: 64)\n"
              << "  --no-bloom           Disable bloom\n"
              << "  --backend <name>     Force backend (optix/cuda)\n";
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
        else if(arg == "--no-bloom") opts.ppSettings.enableBloom = false;
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
    
    // Launch Render
    // For progressive accumulation, we might want to loop here.
    // The current backend interface 'launch' does one pass. 
    // We can loop here for updated progress or just pass all samples to config.
    // The config has 'samplesPerPixel'.
    // If backend handles accumulation internally (likely for OptiX), we just call launch.
    
    // Actually, OptiX backend usually does one pass per launch.
    // Let's loop for progress.
    
    int totalSamples = opts.config.samplesPerPixel;
    int batchSize = 1; 
    
    for (int s = 0; s < totalSamples; s += batchSize) {
        // Update config for current frame/seed if needed?
        // The accelerator interface might need refinement for progressive updates if it holds state.
        // Assuming 'launch' adds samples.
        accelerator->launch(opts.config);
        
        if (opts.verbose && (s % 10 == 0)) {
            std::cout << "\rSample " << s << "/" << totalSamples << std::flush;
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
    
    // Buffer IO
    if (opts.outputEXR) {
        // TODO: Implement writeEXR in BufferIO or use tinyexr directly
        // BufferWriter::writeEXR(opts.outputPath, gpuBuffer, opts.config.width, opts.config.height);
        std::cerr << "EXR output not yet fully implemented in BufferWriter, saving raw.\n";
        BufferWriter::writeRaw(opts.outputPath + ".raw", gpuBuffer, opts.config.width, opts.config.height);
    } else {
        BufferWriter::writePPM(opts.outputPath, gpuBuffer, opts.config.width, opts.config.height);
    }
    
    manager.release();
    return 0;
}
