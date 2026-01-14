// RDOX001A.cpp - OptiX Accelerator Implementation
// Component ID: RDOX001A

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
    
    // Assume PTX is in standard location relative to executable
    // Note: In real implementation, resolve path properly
    const char* ptxPath = "Sirius.Render/ptx/RDOP002A.ptx"; 
    // Trying relative path first
    
    if (!sirius_optix_create_pipeline(m_Handle, ptxPath)) {
        // Fallback to absolute path or other locations if needed? 
        // For now rely on sirius_optix_create_pipeline logic
        m_LastError = "Failed to create pipeline";
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
    sirius_optix_launch(m_Handle, &params);
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
    float phi = 0.0f;
    
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
    params.camera.direction = Sirius::make_float3(-camX/camDist, -camY/camDist, -camZ/camDist);
    params.camera.up = Sirius::make_float3(0.0f, 1.0f, 0.0f);
    params.camera.right = Sirius::make_float3(cosTheta, 0.0f, -sinTheta);
    
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
    
    params.pathTracing.samplesPerPixel = 1; // Accumulation handled by caller loop
    
    return params;
}

} // namespace Sirius::Acceleration::OptiX


// Factory functions moved to ACBF001A.cpp (Backend Factory)

