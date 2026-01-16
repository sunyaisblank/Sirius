// RDOX001A.cpp - OptiX Accelerator Implementation
// Component ID: RDOX001A

// IMPORTANT: Include RDPTX001A.h first to define CUDA types before RDOP003A.h
// This ensures __CUDA_RUNTIME_H__ is defined when RDOP003A.h is processed
#include "RDPTX001A.h"  // PTX pre-compilation utilities (includes cuda_runtime.h)
#include "RDOX001A.h"
#include <cmath>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Sirius::Acceleration::OptiX {

OptiXAccelerator::OptiXAccelerator() {
    m_Handle = sirius_optix_create();
}

OptiXAccelerator::~OptiXAccelerator() {
    cleanup();
    if (m_Handle) {
        sirius_optix_destroy(m_Handle);
        m_Handle = nullptr;
    }
}

DeviceCapabilities OptiXAccelerator::getCapabilities() const {
    DeviceCapabilities caps;
    caps.deviceName = "NVIDIA OptiX Device";
    caps.supportsRTCores = true;
    caps.supportsFP64 = true;
    return caps;
}

bool OptiXAccelerator::initialise(int width, int height) {
    if (!m_Handle) return false;

    m_Width = width;
    m_Height = height;

    if (!sirius_optix_initialize(m_Handle, width, height)) {
        m_LastError = "Failed to initialise OptiX context";
        return false;
    }

    // Use PTX utilities for dynamic architecture selection
    std::string ptxDir;
#ifdef OPTIX_PTX_DIR
    ptxDir = OPTIX_PTX_DIR;
#endif

    // Find the best PTX file for the current device
    std::string ptxPath = sirius::render::findPTXFile("RDOP002A", ptxDir, 0);

    if (ptxPath.empty()) {
        // Fallback to legacy path resolution
        std::vector<std::string> fallbackPaths = {
            "Sirius.Render/ptx/RDOP002A.ptx",
            "./Sirius.Render/ptx/RDOP002A.ptx",
            "../Sirius.Render/ptx/RDOP002A.ptx"
        };
        for (const auto& path : fallbackPaths) {
            if (sirius::render::fileExists(path)) {
                ptxPath = path;
                break;
            }
        }
    }

    if (ptxPath.empty()) {
        m_LastError = "Could not find PTX file for device";
        return false;
    }

    std::cout << "[OptiX] Using PTX: " << ptxPath << std::endl;

    if (!sirius_optix_create_pipeline(m_Handle, ptxPath.c_str())) {
        m_LastError = "Failed to create pipeline with PTX: " + ptxPath;
        return false;
    }

    return true;
}

bool OptiXAccelerator::isInitialised() const {
    return m_Handle && sirius_optix_is_initialized(m_Handle);
}

void OptiXAccelerator::launch(const LaunchConfig& config) {
    if (!isInitialised()) return;

    Sirius::LaunchParams params = convertConfig(config);

    // Update accumulation state
    params.frameIndex = m_FrameCount;
    params.pathTracing.accumulatedFrames = m_FrameCount;

    // Select correct raygen program for metric type
    // 0=Minkowski, 1=Schwarzschild, 2=Kerr, etc.
    sirius_optix_set_metric_type(m_Handle, config.metricType);

    sirius_optix_launch(m_Handle, &params);

    // Increment frame counter after launch
    m_FrameCount++;
}

void OptiXAccelerator::resetAccumulation() {
    m_FrameCount = 0;
}

float* OptiXAccelerator::getFrameBuffer() {
    if (!isInitialised()) return nullptr;
    return sirius_optix_get_frame_buffer(m_Handle);
}

bool OptiXAccelerator::uploadBackground(const uint8_t* data, int width, int height) {
    if (!isInitialised()) return false;
    return sirius_optix_upload_background(m_Handle, data, width, height);
}

void OptiXAccelerator::synchronise() {
    // sirius_optix_launch is synchronous for now in the underlying implementation
}

void OptiXAccelerator::cleanup() {
    if (m_Handle) {
        sirius_optix_cleanup(m_Handle);
    }
}

Sirius::LaunchParams OptiXAccelerator::convertConfig(const LaunchConfig& config) {
    Sirius::LaunchParams params = Sirius::createDefaultLaunchParams();
    
    // Camera Setup
    float r = config.observerDistance;
    float theta = config.observerInclination;
    float phi = config.observerAzimuth;  // Use azimuth for asymmetric composition
    
    float sinTheta = std::sin(theta);
    float cosTheta = std::cos(theta);
    float sinPhi = std::sin(phi);
    float cosPhi = std::cos(phi);
    
    float camX = r * sinTheta * cosPhi;
    float camY = r * cosTheta;
    float camZ = r * sinTheta * sinPhi;
    
    params.camera.position = Sirius::make_float3(camX, camY, camZ);
    params.camera.fov = static_cast<float>(config.cameraFOV * M_PI / 180.0);
    params.camera.aspectRatio = static_cast<float>(m_Width) / m_Height;
    params.camera.observerTheta = theta;
    params.camera.observerPhi = phi;
    
    float camDist = std::sqrt(camX*camX + camY*camY + camZ*camZ);

    // Camera direction: points toward origin (black hole center)
    float dirX = -camX / camDist;
    float dirY = -camY / camDist;
    float dirZ = -camZ / camDist;
    params.camera.direction = Sirius::make_float3(dirX, dirY, dirZ);

    // Camera up: perpendicular to direction, in the "northward" direction along the meridian
    // For an observer at (r, θ, φ), the local "up" points toward decreasing θ (toward north pole)
    // In Cartesian: up = ∂/∂θ (normalized and made perpendicular to direction)
    // ∂r/∂θ = (r cosθ cosφ, -r sinθ, r cosθ sinφ) -> normalized: (cosθ cosφ, -sinθ, cosθ sinφ)
    // But we need to project out the component along direction (Gram-Schmidt)
    float upX_raw = cosTheta * cosPhi;  // Points "up" along the sphere
    float upY_raw = -sinTheta;
    float upZ_raw = cosTheta * sinPhi;

    // Project out component along direction: up = up_raw - (up_raw · dir) * dir
    float dotUpDir = upX_raw * dirX + upY_raw * dirY + upZ_raw * dirZ;
    float upX = upX_raw - dotUpDir * dirX;
    float upY = upY_raw - dotUpDir * dirY;
    float upZ = upZ_raw - dotUpDir * dirZ;

    // Normalize
    float upLen = std::sqrt(upX*upX + upY*upY + upZ*upZ);
    if (upLen > 0.0001f) {
        upX /= upLen;
        upY /= upLen;
        upZ /= upLen;
    } else {
        // Fallback for θ ≈ 0 or π (looking straight down/up)
        upX = 1.0f; upY = 0.0f; upZ = 0.0f;
    }
    params.camera.up = Sirius::make_float3(upX, upY, upZ);

    // Camera right: perpendicular to both direction and up (cross product)
    // right = direction × up
    float rightX = dirY * upZ - dirZ * upY;
    float rightY = dirZ * upX - dirX * upZ;
    float rightZ = dirX * upY - dirY * upX;
    params.camera.right = Sirius::make_float3(rightX, rightY, rightZ);
    
    // Metric
    params.metricType = static_cast<Sirius::MetricType>(config.metricType);
    params.metricParams.family = static_cast<Sirius::MetricFamily>(config.metricFamily);
    
    // Set family parameters
    params.metricParams.kerrSchild.M = config.blackHoleMass;
    params.metricParams.kerrSchild.a = config.blackHoleSpin;
    // Sync legacy
    params.metricParams.M = config.blackHoleMass;
    params.metricParams.a = config.blackHoleSpin;
    
    // Disk
    params.accretionDisk.innerRadius = config.diskInnerRadius;
    params.accretionDisk.outerRadius = config.diskOuterRadius;
    params.accretionDisk.innerTemperature = config.diskTemperature;
    params.accretionDisk.emissionCoefficient = config.diskEmission;
    
    // Integration
    params.integration.maxSteps = config.maxSteps;
    params.integration.maxStepSize = config.maxStepSize;
    
    // Background
    if (config.diskOuterRadius < 0) { // Hack to signal background usage?
         // Handled via separate calls
    }
    unsigned long long bgTex = sirius_optix_get_background_texture(m_Handle);
    if (bgTex) {
        params.useBackgroundTexture = true;
        params.backgroundTexture = bgTex;
    }
    
    // Path tracing accumulation
    params.pathTracing.samplesPerPixel = 1;           // 1 sample per launch
    params.pathTracing.enableAccumulation = true;     // Enable progressive accumulation
    params.pathTracing.useExponentialMA = false;      // Use simple averaging
    params.pathTracing.resetAccumulation = (m_FrameCount == 0);  // Reset on first frame
    params.pathTracing.seed = m_FrameCount;           // Vary seed per frame for Monte Carlo

    // =========================================================================
    // Cinematic Features (Phase 8)
    // =========================================================================

    // Turbulence parameters
    if (config.enableTurbulence) {
        params.volumetricDisk.turbulence.enabled = 1;
        params.volumetricDisk.turbulence.amplitude = config.turbulenceAmplitude;
        params.volumetricDisk.turbulence.outer_scale = config.turbulenceOuterScale;
        params.volumetricDisk.turbulence.inner_scale = config.turbulenceInnerScale;
        params.volumetricDisk.turbulence.octaves = config.turbulenceOctaves;
        params.volumetricDisk.turbulence.seed = config.turbulenceSeed;
        params.volumetricDisk.turbulence.lacunarity = 2.0f;
        params.volumetricDisk.turbulence.persistence = 0.5f;
    } else {
        params.volumetricDisk.turbulence.enabled = 0;
    }

    // Starfield parameters
    if (config.enableStarfield) {
        params.starfield.enabled = 1;
        params.starfield.brightness_scale = config.starfieldBrightness;
    } else {
        params.starfield.enabled = 0;
    }

    // Film simulation parameters
    if (config.enableFilm) {
        params.film.grain_intensity = config.filmGrainIntensity;
        params.film.halation_radius = config.filmHalationRadius;
        params.film.halation_strength = config.filmHalationStrength;
        params.film.vignette_strength = config.filmVignetteStrength;
        // Enable film features bitfield: bit 0=grain, bit 1=halation, bit 2=vignette, bit 3=enabled
        params.film.features = 0x0F;  // All features enabled
    } else {
        params.film.features = 0;
    }

    return params;
}

} // namespace Sirius::Acceleration::OptiX


// Factory functions moved to ACBF001A.cpp (Backend Factory)

