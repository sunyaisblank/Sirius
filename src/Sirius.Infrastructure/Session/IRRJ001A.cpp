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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Sirius {

//==============================================================================
// Constructor / Destructor
//==============================================================================
RenderJob::RenderJob(const RenderConfig& config)
    : m_Config(config)
{
    m_ImageBuffer.resize(config.width * config.height * 4, 0.0f);
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
              << "θ=" << (m_Config.observerPosition[2] * 180.0 / M_PI) << "°" << std::endl;
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
// Black Hole Scanline Rendering (Geometric Approximation)
//==============================================================================
void RenderJob::renderScanlineBlackHole(int y) {
    float aspectRatio = static_cast<float>(m_Config.width) / m_Config.height;
    float fovRad = m_Config.fov * static_cast<float>(M_PI) / 180.0f;
    float tanHalfFov = std::tan(fovRad / 2.0f);
    
    // Observer parameters
    double r_obs = m_Config.observerPosition[1];  // Distance from black hole
    double theta_obs = m_Config.observerPosition[2];  // Inclination
    
    // Black hole parameters
    double M = m_Config.M;
    double a = m_Config.a;
    double r_s = 2.0 * M;  // Schwarzschild radius
    
    // Angular size of photon sphere as seen from observer
    double r_photon = 1.5 * r_s;  // Photon sphere radius (Schwarzschild)
    if (std::abs(a) > 0.01) {
        // Kerr correction (approximate)
        r_photon = 1.5 * r_s * (1.0 - a * a / (6.0 * M * M));
    }
    
    // Angular size of shadow (impact parameter)
    double b_crit = r_photon * std::sqrt(r_obs / (r_obs - r_s));
    double angular_size = std::atan2(b_crit, r_obs);
    
    // Disk parameters
    double r_isco = 6.0 * M;  // ISCO for Schwarzschild
    if (std::abs(a) > 0.01) {
        // Kerr ISCO (prograde approximation)
        r_isco = 6.0 * M * (1.0 - 0.5 * a / M);
    }
    double r_disk_inner = std::max(m_Config.diskInnerRadius, r_isco);
    double r_disk_outer = m_Config.diskOuterRadius;
    
    for (int x = 0; x < m_Config.width; ++x) {
        float r = 0.0f, g = 0.0f, b = 0.0f;
        
        for (int s = 0; s < m_Config.samplesPerPixel; ++s) {
            float jitterX = (s % 4) * 0.25f + 0.125f;
            float jitterY = (s / 4 % 4) * 0.25f + 0.125f;
            
            // Normalized screen coordinates
            float u = (2.0f * (x + jitterX) / m_Config.width - 1.0f) * aspectRatio * tanHalfFov;
            float v = (1.0f - 2.0f * (y + jitterY) / m_Config.height) * tanHalfFov;
            
            // Angular offset from centre
            double theta_ray = std::sqrt(u * u + v * v);
            double phi_ray = std::atan2(v, u);
            
            // Impact parameter
            double b = r_obs * std::tan(theta_ray);
            
            float sampleR = 0.0f, sampleG = 0.0f, sampleB = 0.0f;
            
            if (b < b_crit * 0.95) {
                // Inside shadow - black
                sampleR = sampleG = sampleB = 0.0f;
            } else if (b < b_crit * 1.05) {
                // Photon ring - bright ring at shadow edge
                float ring = 1.0f - std::abs(static_cast<float>(b / b_crit) - 1.0f) * 20.0f;
                ring = std::max(0.0f, ring);
                sampleR = ring * 1.5f;
                sampleG = ring * 1.3f;
                sampleB = ring * 1.0f;
            } else {
                // Outside shadow - check for disk and background
                
                // Disk intersection (thin disk at equatorial plane)
                double y_disk = b * std::sin(phi_ray);  // y coordinate in screen plane
                double disk_height = b * std::cos(theta_obs) * std::cos(phi_ray);
                
                // Check if ray passes through disk plane
                bool hitDisk = false;
                double disk_r = 0.0;
                
                if (std::abs(disk_height) < b * 0.1 + 1.0) {
                    // Approximate disk radius at this position
                    disk_r = std::sqrt(b * b + r_obs * r_obs - 2.0 * b * r_obs * std::cos(theta_ray));
                    
                    // Apply gravitational lensing (bends rays around black hole)
                    double lens_factor = 1.0 + r_s / (2.0 * disk_r);
                    disk_r *= lens_factor;
                    
                    if (disk_r >= r_disk_inner && disk_r <= r_disk_outer) {
                        hitDisk = true;
                    }
                }
                
                if (hitDisk) {
                    // Disk emission - temperature profile
                    double T = std::pow(disk_r / r_isco, -0.75);
                    T = std::max(0.0, T);
                    
                    // Doppler shift (approaching side brighter, receding dimmer)
                    double doppler = 1.0 + 0.3 * std::sin(phi_ray);
                    T *= doppler;
                    
                    // Colour based on temperature
                    float intensity = static_cast<float>(T * m_Config.diskMdot * 1e4);
                    if (T > 0.6) {
                        sampleR = intensity * 0.95f;
                        sampleG = intensity * 0.9f;
                        sampleB = intensity * 0.85f;
                    } else if (T > 0.3) {
                        sampleR = intensity;
                        sampleG = intensity * 0.6f;
                        sampleB = intensity * 0.2f;
                    } else {
                        sampleR = intensity;
                        sampleG = intensity * 0.3f;
                        sampleB = intensity * 0.1f;
                    }
                } else {
                    // Background starfield with gravitational lensing
                    double lensed_theta = theta_ray + r_s / (2.0 * b);
                    double lensed_phi = phi_ray;
                    
                    // Hash for star positions
                    auto hash = [](double x, double y) {
                        int xi = static_cast<int>(std::floor(x * 50));
                        int yi = static_cast<int>(std::floor(y * 50));
                        return (xi * 73856093 ^ yi * 19349663) & 0xFFFFFF;
                    };
                    
                    int h = hash(lensed_theta, lensed_phi);
                    float star = (h % 1000 < 3) ? 1.0f : 0.0f;
                    float brightness = 0.5f + 0.5f * ((h >> 8) % 256) / 255.0f;
                    
                    // Star colour
                    float temp = ((h >> 4) % 256) / 255.0f;
                    if (temp < 0.3f) {
                        sampleR = star * brightness * 0.7f;
                        sampleG = star * brightness * 0.85f;
                        sampleB = star * brightness;
                    } else if (temp < 0.7f) {
                        sampleR = star * brightness;
                        sampleG = star * brightness * 0.95f;
                        sampleB = star * brightness * 0.85f;
                    } else {
                        sampleR = star * brightness;
                        sampleG = star * brightness * 0.5f;
                        sampleB = star * brightness * 0.3f;
                    }
                    
                    // Faint background
                    sampleR += 0.003f;
                    sampleG += 0.002f;
                    sampleB += 0.005f;
                }
            }
            
            r += sampleR;
            g += sampleG;
            b += sampleB;
        }
        
        int idx = (y * m_Config.width + x) * 4;
        m_ImageBuffer[idx + 0] = r / m_Config.samplesPerPixel;
        m_ImageBuffer[idx + 1] = g / m_Config.samplesPerPixel;
        m_ImageBuffer[idx + 2] = b / m_Config.samplesPerPixel;
        m_ImageBuffer[idx + 3] = 1.0f;
    }
}

//==============================================================================
// Legacy Scanline Rendering (simple starfield)
//==============================================================================
void RenderJob::renderScanline(int y) {
    float aspectRatio = static_cast<float>(m_Config.width) / m_Config.height;
    float fovRad = m_Config.fov * static_cast<float>(M_PI) / 180.0f;
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

} // namespace Sirius
