// SNRS001A.cpp - Render Session Implementation
// Component ID: SNRS001A (Session/RenderSession)

#include "SNRS001A.h"
#include <PPOP001A.h>
#include <fstream>
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Sirius {

//==============================================================================
// Initialisation
//==============================================================================
void RenderSession::initialise() {
    try {
        std::cout << "[Session] Initialising render..." << std::endl;
        std::cout << "  Resolution: " << m_Config.width << " x " << m_Config.height << std::endl;
        std::cout << "  Tile size:  " << m_Config.tileSize << std::endl;
        std::cout << "  Samples:    " << m_Config.samplesPerPixel << std::endl;
        
        // Initialise tiles
        m_Tiles.initialise(m_Config.width, m_Config.height, m_Config.tileSize);
        std::cout << "  Tiles:      " << m_Tiles.getTileCount() << " (spiral order)" << std::endl;
        
        // Initialise display buffer
        m_Display.initialise(m_Config.width, m_Config.height);
        
        // Initialise progress tracker
        m_Progress.setTotals(m_Tiles.getTileCount(), m_Config.samplesPerPixel);
        m_Progress.start();
        
        // Transition to Scheduling
        m_FSM.process(SessionEvent::Ready);
    }
    catch (const std::exception& e) {
        m_ErrorMessage = e.what();
        m_FSM.process(SessionEvent::Error);
    }
}

//==============================================================================
// Tile Scheduling
//==============================================================================
void RenderSession::scheduleNextTile() {
    // Check for cancellation
    if (m_Progress.getCancellationToken().isCancelled()) {
        m_FSM.process(SessionEvent::Cancel);
        return;
    }
    
    // Get next tile
    Tile* tile = m_Tiles.getNextTile();
    
    if (tile == nullptr) {
        // All tiles complete
        m_FSM.process(SessionEvent::AllTilesComplete);
        return;
    }
    
    // Signal tile available and render
    m_FSM.process(SessionEvent::TileAvailable);
    
    // Render the tile synchronously (for now - can add threading later)
    renderTile(tile);
}

//==============================================================================
// Tile Rendering (Black Hole Approximation)
//==============================================================================
void RenderSession::renderTile(Tile* tile) {
    if (!tile) return;
    
    // Allocate tile buffer
    std::vector<float> tileBuffer(tile->width * tile->height * 4, 0.0f);
    
    // Black hole rendering parameters
    float aspectRatio = static_cast<float>(m_Config.width) / m_Config.height;
    float fovRad = m_Config.cameraFOV * static_cast<float>(M_PI) / 180.0f;
    float tanHalfFov = std::tan(fovRad / 2.0f);
    double r_obs = m_Config.observerDistance;
    double theta_obs = m_Config.observerInclination;
    double M = m_Config.blackHoleMass;
    double a = m_Config.blackHoleSpin;
    
    // Schwarzschild radius and shadow
    double r_s = 2.0 * M;
    double r_shadow = r_s * std::sqrt(27.0) / 2.0;
    if (std::abs(a) > 0.01) {
        r_shadow *= (1.0 + 0.1 * std::abs(a));
    }
    double angular_shadow = r_shadow / r_obs;
    
    // Disk parameters
    double r_isco = 6.0 * M;
    double r_disk_outer = 20.0 * M;
    
    // Render each pixel in tile
    for (int ty = 0; ty < tile->height; ++ty) {
        for (int tx = 0; tx < tile->width; ++tx) {
            int px = tile->x + tx;
            int py = tile->y + ty;
            
            // Normalised device coordinates
            float u = (2.0f * (px + 0.5f) / m_Config.width - 1.0f) * aspectRatio;
            float v = 1.0f - 2.0f * (py + 0.5f) / m_Config.height;
            
            // Ray angles
            double theta_ray = u * tanHalfFov;
            double phi_ray = v * tanHalfFov;
            double ray_angle = std::sqrt(theta_ray * theta_ray + phi_ray * phi_ray);
            
            float r = 0.0f, g = 0.0f, b = 0.0f;
            
            // Black hole shadow
            if (ray_angle < angular_shadow) {
                r = g = b = 0.0f;
            }
            // Photon ring
            else if (ray_angle < angular_shadow * 1.05) {
                float ring = 1.0f - (static_cast<float>(ray_angle / angular_shadow) - 1.0f) / 0.05f;
                r = 1.0f * ring;
                g = 0.8f * ring;
                b = 0.5f * ring;
            }
            // Accretion disk
            else if (std::abs(phi_ray) < 0.15 && ray_angle > angular_shadow) {
                double r_disk = r_obs * ray_angle;
                if (r_disk > r_isco && r_disk < r_disk_outer) {
                    // Temperature profile
                    double T = std::pow((r_isco / r_disk), 0.75);
                    
                    // Doppler shift
                    double v_orbital = std::sqrt(M / r_disk);
                    double doppler = 1.0 + v_orbital * theta_ray / std::abs(theta_ray + 0.001);
                    doppler = std::clamp(doppler, 0.3, 2.0);
                    
                    r = static_cast<float>(T * doppler * 2.0);
                    g = static_cast<float>(T * doppler * 1.5);
                    b = static_cast<float>(T * doppler * 0.8);
                }
            }
            // Starfield
            else {
                double hash = std::sin(px * 12.9898 + py * 78.233) * 43758.5453;
                hash = hash - std::floor(hash);
                if (hash > 0.997) {
                    float brightness = static_cast<float>((hash - 0.997) / 0.003);
                    r = g = b = brightness;
                }
            }
            
            // Store in tile buffer
            int idx = (ty * tile->width + tx) * 4;
            tileBuffer[idx + 0] = r;
            tileBuffer[idx + 1] = g;
            tileBuffer[idx + 2] = b;
            tileBuffer[idx + 3] = 1.0f;
        }
    }
    
    // Update display buffer
    m_Display.updateTile(tile->x, tile->y, tile->width, tile->height, tileBuffer.data());
    
    // Mark tile complete
    m_Tiles.completeTile(tile->id);
    m_Progress.completeTile(tile->pixelCount());
    
    // Print progress
    int percent = static_cast<int>(m_Progress.getProgress() * 100);
    std::cout << "\r[Session] Progress: " << percent << "% | ETA: " 
              << m_Progress.getETAString() << "    " << std::flush;
    
    // Transition back to scheduling
    m_FSM.process(SessionEvent::TileComplete);
}

//==============================================================================
// Output Writing
//==============================================================================
void RenderSession::writeOutput() {
    try {
        std::cout << "\n[Session] Writing output..." << std::endl;
        
        // Apply post-processing
        PostProcessConfig ppConfig;
        ppConfig.tonemapper = TonemapType::ACES;
        ppConfig.exposure = m_Config.exposure;
        ppConfig.enableBloom = m_Config.enableBloom;
        ppConfig.bloomIntensity = m_Config.bloomIntensity;
        ppConfig.bloomThreshold = 0.8f;
        
        PostProcessor::process(m_Display.getFloatBuffer(), m_Config.width, m_Config.height, ppConfig);
        
        // Write PPM
        std::ofstream file(m_Config.outputPath, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Failed to open output file: " + m_Config.outputPath);
        }
        
        file << "P6\n" << m_Config.width << " " << m_Config.height << "\n255\n";
        
        const float* data = m_Display.getFloatData();
        for (int y = 0; y < m_Config.height; ++y) {
            for (int x = 0; x < m_Config.width; ++x) {
                int idx = (y * m_Config.width + x) * 4;
                unsigned char rgb[3] = {
                    static_cast<unsigned char>(std::clamp(data[idx + 0], 0.0f, 1.0f) * 255.0f),
                    static_cast<unsigned char>(std::clamp(data[idx + 1], 0.0f, 1.0f) * 255.0f),
                    static_cast<unsigned char>(std::clamp(data[idx + 2], 0.0f, 1.0f) * 255.0f)
                };
                file.write(reinterpret_cast<char*>(rgb), 3);
            }
        }
        
        std::cout << "[Session] Wrote: " << m_Config.outputPath << std::endl;
        m_FSM.process(SessionEvent::OutputWritten);
    }
    catch (const std::exception& e) {
        m_ErrorMessage = e.what();
        m_FSM.process(SessionEvent::Error);
    }
}

//==============================================================================
// Session End
//==============================================================================
void RenderSession::onSessionEnd(SessionState state) {
    std::cout << "\n[Session] Finished with state: " << stateName(state) << std::endl;
    
    std::string message;
    switch (state) {
        case SessionState::Complete:
            message = "Render completed successfully";
            break;
        case SessionState::Failed:
            message = "Render failed: " + m_ErrorMessage;
            break;
        case SessionState::Cancelled:
            message = "Render cancelled by user";
            break;
        default:
            message = "Unknown end state";
    }
    
    if (m_CompletionCallback) {
        m_CompletionCallback(state, message);
    }
}

} // namespace Sirius
