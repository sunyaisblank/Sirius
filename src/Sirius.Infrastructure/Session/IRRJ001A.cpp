// IRRJ001A.cpp - Render Job Implementation
//
// Executes a single render task with black hole visualisation.
// Integrates ICamera for ray generation and PostProcessor for output.

#include "IRRJ001A.h"
#include <CMBS001A.h>   // ICamera interface
#include <PPOP001A.h>   // Post-processing pipeline
// Note: EXR output deferred - using PPM for now

#include <chrono>
#include <thread>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <memory>

// Unified constants
#include <PHCN001A.h>
using Sirius::Constants::Math::PI;

namespace Sirius {

//==============================================================================
// Constructor / Destructor
//==============================================================================
RenderJob::RenderJob(const RenderConfig& config)
    : m_Config(config)
{
    m_ImageBuffer.resize(config.width * config.height * 4, 0.0f);

    // =================================================================
    // Priority 1 Fix: Initialize physics components
    // =================================================================

    // Create metric from config (Kerr-Schild family)
    KerrSchildParams params;
    params.M = config.M;
    params.a = config.a * config.M;  // spin as a/M
    params.Q = 0.0;  // No charge for now
    params.Lambda = 0.0;  // No cosmological constant
    m_Metric = std::make_unique<KerrSchildFamily>(params);

    // Create camera
    CameraConfig camConfig;
    camConfig.r = config.observerPosition[1];
    camConfig.theta = config.observerPosition[2];
    camConfig.phi = config.observerPosition[3];
    camConfig.fov = config.fov;
    camConfig.width = config.width;
    camConfig.height = config.height;
    m_Camera = std::make_unique<PinholeCamera>(camConfig);

    // Create geodesic tracer
    TracerConfig tracerConfig;
    tracerConfig.escape_radius = 100.0f;
    tracerConfig.horizon_factor = 1.05f;
    tracerConfig.max_steps = config.maxIntegrationSteps;
    tracerConfig.integrator.abs_tolerance = static_cast<float>(config.integrationTolerance);
    tracerConfig.integrator.rel_tolerance = static_cast<float>(config.integrationTolerance);

    // Compute ISCO for disk inner radius
    double a_dim = config.a;
    if (std::abs(a_dim) < 0.01) {
        tracerConfig.disk_inner = static_cast<float>(6.0 * config.M);  // Schwarzschild ISCO
    } else {
        // Kerr ISCO (prograde, simplified)
        tracerConfig.disk_inner = static_cast<float>(config.M * (3.0 + 3.0 * std::pow(1.0 - a_dim*a_dim, 1.0/3.0)));
    }
    tracerConfig.disk_outer = static_cast<float>(config.diskOuterRadius);
    tracerConfig.enable_disk = config.enableDisk;

    m_Tracer = std::make_unique<GeodesicTracer>(m_Metric.get(), tracerConfig);
}

RenderJob::~RenderJob() {
    if (isRunning()) {
        cancel();
    }
}

//==============================================================================
// Status Queries
//==============================================================================
bool RenderJob::isComplete() const {
    JobStatus status = m_Status.load();
    return status == JobStatus::Completed || 
           status == JobStatus::Failed || 
           status == JobStatus::Cancelled;
}

bool RenderJob::isRunning() const {
    return m_Status.load() == JobStatus::Running;
}

double RenderJob::getElapsedTime() const {
    if (m_StartTime == 0.0) return 0.0;
    
    auto now = std::chrono::high_resolution_clock::now();
    double currentTime = std::chrono::duration<double>(now.time_since_epoch()).count();
    
    if (isComplete()) {
        return m_EndTime - m_StartTime;
    }
    return currentTime - m_StartTime;
}

double RenderJob::getEstimatedTimeRemaining() const {
    float progress = m_Progress.load();
    if (progress <= 0.0f) return -1.0;
    
    double elapsed = getElapsedTime();
    return elapsed * (1.0 - progress) / progress;
}

//==============================================================================
// Callbacks
//==============================================================================
void RenderJob::setProgressCallback(ProgressCallback callback) {
    m_ProgressCallback = std::move(callback);
}

void RenderJob::setCompletionCallback(CompletionCallback callback) {
    m_CompletionCallback = std::move(callback);
}

//==============================================================================
// Execution Control
//==============================================================================
void RenderJob::pause() {
    if (m_Status.load() == JobStatus::Running) {
        m_PauseRequested.store(true);
    }
}

void RenderJob::resume() {
    if (m_Status.load() == JobStatus::Paused) {
        m_PauseRequested.store(false);
        m_Status.store(JobStatus::Running);
    }
}

void RenderJob::cancel() {
    m_Status.store(JobStatus::Cancelled);
}

//==============================================================================
// Main Execution
//==============================================================================
void RenderJob::execute() {
    if (m_Status.load() != JobStatus::Pending) {
        std::cerr << "[RenderJob] Cannot execute: job not in Pending state" << std::endl;
        return;
    }
    
    m_Status.store(JobStatus::Running);
    
    auto startClock = std::chrono::high_resolution_clock::now();
    m_StartTime = std::chrono::duration<double>(startClock.time_since_epoch()).count();
    
    std::cout << "================================================================================\n";
    std::cout << "                       SIRIUS BLACK HOLE RENDERER                               \n";
    std::cout << "================================================================================\n";
    std::cout << "[RenderJob] Resolution: " << m_Config.width << "x" << m_Config.height << std::endl;
    std::cout << "[RenderJob] Samples:    " << m_Config.samplesPerPixel << " spp" << std::endl;
    std::cout << "[RenderJob] Metric:     " << m_Config.metricName << " (M=" << m_Config.M << ", a=" << m_Config.a << ")" << std::endl;
    std::cout << "[RenderJob] Observer:   r=" << m_Config.observerPosition[1] << "M, "
              << "θ=" << (m_Config.observerPosition[2] * 180.0 / PI) << "°" << std::endl;
    std::cout << "[RenderJob] Output:     " << m_Config.outputPath << std::endl;
    std::cout << "--------------------------------------------------------------------------------\n";
    
    int totalScanlines = m_Config.height;
    int completedScanlines = 0;
    
    try {
        for (int y = 0; y < m_Config.height; ++y) {
            if (m_Status.load() == JobStatus::Cancelled) {
                std::cout << "\n[RenderJob] Cancelled at scanline " << y << std::endl;
                break;
            }
            
            while (m_PauseRequested.load()) {
                m_Status.store(JobStatus::Paused);
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
                if (m_Status.load() == JobStatus::Cancelled) break;
            }
            
            if (m_Status.load() == JobStatus::Cancelled) break;
            m_Status.store(JobStatus::Running);
            
            renderScanlineBlackHole(y);
            
            completedScanlines++;
            m_Progress.store(static_cast<float>(completedScanlines) / totalScanlines);
            m_SamplesComplete.store(completedScanlines * m_Config.width * m_Config.samplesPerPixel);
            
            if (completedScanlines % 20 == 0 || completedScanlines == totalScanlines) {
                int percent = static_cast<int>(m_Progress.load() * 100);
                double eta = getEstimatedTimeRemaining();
                std::cout << "\r[RenderJob] Progress: " << percent << "%";
                if (eta > 0) {
                    int mins = static_cast<int>(eta / 60);
                    int secs = static_cast<int>(eta) % 60;
                    std::cout << " | ETA: " << mins << "m " << secs << "s";
                }
                std::cout << std::flush;
            }
        }
        
        if (m_Status.load() != JobStatus::Cancelled) {
            std::cout << "\n";
            writeOutput();
            m_Status.store(JobStatus::Completed);
            
            auto endClock = std::chrono::high_resolution_clock::now();
            m_EndTime = std::chrono::duration<double>(endClock.time_since_epoch()).count();
            
            std::cout << "[RenderJob] Completed in " << getElapsedTime() << "s" << std::endl;
            
            if (m_CompletionCallback) {
                m_CompletionCallback(JobStatus::Completed, "Render completed successfully");
            }
        }
    }
    catch (const std::exception& e) {
        m_Status.store(JobStatus::Failed);
        std::cerr << "\n[RenderJob] Error: " << e.what() << std::endl;
        
        if (m_CompletionCallback) {
            m_CompletionCallback(JobStatus::Failed, e.what());
        }
    }
}

//==============================================================================
// Black Hole Scanline Rendering (Geodesic Integration - Priority 1 Fix)
//==============================================================================
void RenderJob::renderScanlineBlackHole(int y) {
    for (int x = 0; x < m_Config.width; ++x) {
        float r_acc = 0.0f, g_acc = 0.0f, b_acc = 0.0f;

        // Multi-sample anti-aliasing (stratified sampling)
        int spp = std::max(1, m_Config.samplesPerPixel);
        int grid_size = static_cast<int>(std::sqrt(static_cast<float>(spp)));
        if (grid_size < 1) grid_size = 1;

        for (int sy = 0; sy < grid_size; ++sy) {
            for (int sx = 0; sx < grid_size; ++sx) {
                // Sub-pixel offset
                float u = (sx + 0.5f) / grid_size;
                float v = (sy + 0.5f) / grid_size;

                // Generate camera ray for this pixel sample
                CameraRay camRay = m_Camera->generateRay(x, y, u, v);

                // Trace ray through curved spacetime
                TraceResult result = m_Tracer->trace(camRay);

                float sr = 0.0f, sg = 0.0f, sb = 0.0f;

                switch (result.outcome) {
                    case TraceResult::Outcome::HORIZON:
                        // Black hole shadow - pure black
                        sr = sg = sb = 0.0f;
                        break;

                    case TraceResult::Outcome::DISK_HIT: {
                        // Sample accretion disk emission
                        // Blackbody color approximation based on temperature
                        float T = result.disk_temperature;

                        // Apply redshift/blueshift correction
                        T *= result.redshift;

                        // Simple blackbody color mapping (Planckian locus approximation)
                        // Higher T = bluer, lower T = redder
                        sr = std::min(2.0f, T * 2.0f);
                        sg = std::min(1.5f, T * 1.5f);
                        sb = std::min(1.0f, T * 0.8f);

                        // Add slight phi-dependent variation (orbital motion)
                        float phi_factor = 1.0f + 0.2f * std::cos(result.disk_phi);
                        sr *= phi_factor;
                        sg *= phi_factor;
                        sb *= phi_factor;
                        break;
                    }

                    case TraceResult::Outcome::ESCAPED: {
                        // Sample background starfield
                        float brightness = sampleStarfield(result.final_direction);
                        sr = sg = sb = brightness;
                        break;
                    }

                    case TraceResult::Outcome::MAX_STEPS:
                    default:
                        // Integration limit reached - render as faint grey
                        sr = sg = sb = 0.01f;
                        break;
                }

                r_acc += sr;
                g_acc += sg;
                b_acc += sb;
            }
        }

        // Average samples
        int total_samples = grid_size * grid_size;
        float inv_samples = 1.0f / static_cast<float>(total_samples);

        int idx = (y * m_Config.width + x) * 4;
        m_ImageBuffer[idx + 0] = r_acc * inv_samples;
        m_ImageBuffer[idx + 1] = g_acc * inv_samples;
        m_ImageBuffer[idx + 2] = b_acc * inv_samples;
        m_ImageBuffer[idx + 3] = 1.0f;
    }
}

//==============================================================================
// Legacy Scanline Rendering (simple starfield)
//==============================================================================
void RenderJob::renderScanline(int y) {
    float aspectRatio = static_cast<float>(m_Config.width) / m_Config.height;
    float fovRad = m_Config.fov * static_cast<float>(PI) / 180.0f;
    float tanHalfFov = std::tan(fovRad / 2.0f);
    
    for (int x = 0; x < m_Config.width; ++x) {
        float u = (2.0f * x / m_Config.width - 1.0f) * aspectRatio * tanHalfFov;
        float v = (1.0f - 2.0f * y / m_Config.height) * tanHalfFov;
        
        auto hash = [](float a, float b) {
            int xi = static_cast<int>(std::floor(a * 50));
            int yi = static_cast<int>(std::floor(b * 50));
            return (xi * 73856093 ^ yi * 19349663) & 0xFFFFFF;
        };
        
        int h = hash(v, u);
        float star = (h % 1000 < 3) ? 1.0f : 0.0f;
        
        int idx = (y * m_Config.width + x) * 4;
        m_ImageBuffer[idx + 0] = star + 0.003f;
        m_ImageBuffer[idx + 1] = star * 0.9f + 0.002f;
        m_ImageBuffer[idx + 2] = star * 0.95f + 0.005f;
        m_ImageBuffer[idx + 3] = 1.0f;
    }
}

//==============================================================================
// Progress Reporting
//==============================================================================
void RenderJob::reportProgress() {
    float progress = m_Progress.load();
    int percent = static_cast<int>(progress * 100);
    double eta = getEstimatedTimeRemaining();
    
    std::cout << "[RenderJob] Progress: " << percent << "% | ETA: ";
    if (eta > 0) {
        int hours = static_cast<int>(eta / 3600);
        int mins = static_cast<int>((eta - hours * 3600) / 60);
        int secs = static_cast<int>(eta) % 60;
        if (hours > 0) std::cout << hours << "h ";
        if (mins > 0 || hours > 0) std::cout << mins << "m ";
        std::cout << secs << "s";
    } else {
        std::cout << "calculating...";
    }
    std::cout << std::endl;
    
    if (m_ProgressCallback) {
        m_ProgressCallback(progress, m_SamplesComplete.load(), 
                          m_Config.width * m_Config.height * m_Config.samplesPerPixel);
    }
}

//==============================================================================
// Output Writing with Post-Processing Pipeline
//==============================================================================
void RenderJob::writeOutput() {
    std::string path = m_Config.outputPath;
    
    // Apply post-processing pipeline
    applyPostProcessing();
    
    // Determine output format
    size_t dotPos = path.rfind('.');
    std::string ext = (dotPos != std::string::npos) ? path.substr(dotPos) : "";
    
    // Convert to lowercase for comparison
    std::string extLower = ext;
    for (char& c : extLower) c = static_cast<char>(std::tolower(c));
    
    if (extLower == ".exr") {
        writeEXR(path);
    } else {
        // Default to PPM
        if (ext != ".ppm") {
            path = (dotPos != std::string::npos) ? path.substr(0, dotPos) + ".ppm" : path + ".ppm";
        }
        writePPM(path);
    }
}

void RenderJob::applyPostProcessing() {
    // Create post-process configuration from RenderConfig
    PostProcessConfig ppConfig;
    
    // Tonemapping
    if (m_Config.enableTonemapping) {
        switch (m_Config.tonemapper) {
            case RenderConfig::TonemapOperator::ACES:
                ppConfig.tonemapper = TonemapType::ACES;
                break;
            case RenderConfig::TonemapOperator::Reinhard:
                ppConfig.tonemapper = TonemapType::Reinhard;
                break;
            case RenderConfig::TonemapOperator::Filmic:
                ppConfig.tonemapper = TonemapType::Filmic;
                break;
            default:
                ppConfig.tonemapper = TonemapType::ACES;
        }
    } else {
        ppConfig.tonemapper = TonemapType::None;
    }
    
    ppConfig.exposure = m_Config.exposure;
    ppConfig.gamma = 2.2f;
    
    // Bloom
    ppConfig.enableBloom = m_Config.enableBloom;
    ppConfig.bloomIntensity = m_Config.bloomIntensity;
    ppConfig.bloomThreshold = m_Config.bloomThreshold;
    ppConfig.bloomRadius = 8;
    
    // Apply full pipeline
    std::cout << "[RenderJob] Applying post-processing..." << std::endl;
    PostProcessor::process(m_ImageBuffer, m_Config.width, m_Config.height, ppConfig);
}

void RenderJob::writePPM(const std::string& path) {
    std::ofstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Failed to open output file: " + path);
    }
    
    file << "P6\n" << m_Config.width << " " << m_Config.height << "\n255\n";
    
    for (int y = 0; y < m_Config.height; ++y) {
        for (int x = 0; x < m_Config.width; ++x) {
            int idx = (y * m_Config.width + x) * 4;
            
            // Buffer already processed by PostProcessor (includes gamma)
            float r = std::clamp(m_ImageBuffer[idx + 0], 0.0f, 1.0f);
            float g = std::clamp(m_ImageBuffer[idx + 1], 0.0f, 1.0f);
            float b = std::clamp(m_ImageBuffer[idx + 2], 0.0f, 1.0f);
            
            unsigned char rgb[3] = {
                static_cast<unsigned char>(r * 255.0f),
                static_cast<unsigned char>(g * 255.0f),
                static_cast<unsigned char>(b * 255.0f)
            };
            file.write(reinterpret_cast<char*>(rgb), 3);
        }
    }
    
    std::cout << "[RenderJob] Wrote PPM output to: " << path << std::endl;
}

void RenderJob::writeEXR(const std::string& path) {
    // EXR support deferred - fall back to PPM
    // Note: SiriusRender executable supports EXR via tinyexr
    std::string ppmPath = path.substr(0, path.rfind('.')) + ".ppm";
    std::cout << "[RenderJob] EXR deferred - writing PPM instead." << std::endl;
    writePPM(ppmPath);
}

//==============================================================================
// Starfield Background Sampling
//==============================================================================
float RenderJob::sampleStarfield(const Vec4& direction) const {
    // Convert direction to spherical angles for starfield lookup
    // direction is in Cartesian (x, y, z)
    double x = direction(1);
    double y = direction(2);
    double z = direction(3);

    // Normalize
    double len = std::sqrt(x*x + y*y + z*z);
    if (len < 1e-10) return 0.0f;

    x /= len;
    y /= len;
    z /= len;

    // Convert to spherical (theta, phi)
    double theta = std::acos(std::clamp(z, -1.0, 1.0));
    double phi = std::atan2(y, x);
    if (phi < 0) phi += 2.0 * PI;

    // Procedural starfield based on direction hash
    // This creates a deterministic star pattern
    double hash_input = theta * 1000.0 + phi * 100.0;
    double hash = std::sin(hash_input * 12.9898) * 43758.5453;
    hash = hash - std::floor(hash);

    // Sparse stars
    if (hash > 0.997) {
        // Star brightness varies
        float brightness = static_cast<float>((hash - 0.997) / 0.003);

        // Add some color variation based on secondary hash
        double color_hash = std::sin(hash_input * 78.233 + 1.0) * 43758.5453;
        color_hash = color_hash - std::floor(color_hash);

        // Some stars are brighter
        if (color_hash > 0.9) {
            brightness *= 2.0f;
        }

        return std::min(brightness, 1.0f);
    }

    // Background glow (very faint)
    return 0.001f;
}

} // namespace Sirius
