// ACIB001A.h - Accelerator Backend Interface
// Component ID: ACIB001A (Acceleration/Interface/Backend)
//
// Abstract interface for hardware-accelerated ray tracing backends.
// Implements Strategy pattern for backend switching (OptiX, CUDA, future HIP/Metal).
//
// DESIGN RATIONALE:
// - Unified API regardless of underlying hardware
// - Runtime backend selection based on availability
// - Memory management abstracted behind interface

#pragma once

#include <cstdint>
#include <string>
#include <memory>
#include <vector>

namespace Sirius::Acceleration {

//==============================================================================
// Backend Type Enumeration
//==============================================================================
enum class BackendType : uint8_t {
    None,       ///< No backend selected
    OptiX,      ///< NVIDIA OptiX (RTX hardware)
    CUDA,       ///< Raw CUDA (compute-only)
    HIP,        ///< AMD HIP (future)
    OneAPI,     ///< Intel OneAPI (future)
    Metal,      ///< Apple Metal (future)
    CPU         ///< CPU fallback with SIMD (future)
};

constexpr const char* backendName(BackendType type) {
    switch (type) {
        case BackendType::None:   return "None";
        case BackendType::OptiX:  return "OptiX";
        case BackendType::CUDA:   return "CUDA";
        case BackendType::HIP:    return "HIP";
        case BackendType::OneAPI: return "OneAPI";
        case BackendType::Metal:  return "Metal";
        case BackendType::CPU:    return "CPU";
    }
    return "Unknown";
}

//==============================================================================
// Device Capabilities
//==============================================================================
struct DeviceCapabilities {
    std::string deviceName;
    size_t totalMemory = 0;         ///< Total VRAM in bytes
    size_t freeMemory = 0;          ///< Available VRAM
    int computeCapability = 0;      ///< CUDA compute capability (major*10 + minor)
    bool supportsRTCores = false;   ///< Hardware ray tracing support
    bool supportsFP64 = false;      ///< Double precision support
    int multiprocessorCount = 0;
    int maxThreadsPerBlock = 0;
};

//==============================================================================
// Launch Configuration
//==============================================================================
struct LaunchConfig {
    int width = 1920;
    int height = 1080;
    int samplesPerPixel = 64;
    int tileSize = 64;
    
    // Black hole parameters
    float blackHoleMass = 1.0f;
    float blackHoleSpin = 0.0f;
    float observerDistance = 100.0f;  // Far enough for cinematic wide shot
    float observerInclination = 1.5708f;  // 90°
    float observerAzimuth = 0.0f;         // Camera azimuth angle (radians)
    float cameraFOV = 45.0f;              // Narrower FOV for cinematic look

    // Accretion disk
    float diskInnerRadius = 6.0f;   // ISCO
    float diskOuterRadius = 30.0f;  // Extended for wider composition
    float diskTemperature = 6500.0f;  // Warm orange (Interstellar look)
    float diskEmission = 2.5f;        // Emission coefficient
    
    // Integration
    float maxStepSize = 0.5f;
    int maxSteps = 10000;

    // Metric Selection
    int metricType = 2; // Default to Kerr (2) instead of Minkowski (0)
    int metricFamily = 0; // 0=KerrSchild, 2=WarpDrive

    // Cinematic Features (Phase 8)
    // Turbulence
    bool enableTurbulence = false;
    float turbulenceAmplitude = 0.3f;
    float turbulenceOuterScale = 5.0f;
    float turbulenceInnerScale = 0.1f;
    uint32_t turbulenceOctaves = 6;
    uint32_t turbulenceSeed = 12345;

    // Starfield
    bool enableStarfield = false;
    float starfieldBrightness = 1.0f;

    // Film simulation
    bool enableFilm = false;
    float filmGrainIntensity = 0.025f;
    float filmHalationRadius = 8.0f;
    float filmHalationStrength = 0.15f;
    float filmVignetteStrength = 0.3f;
};

//==============================================================================
// IAccelerator - Abstract Backend Interface
//==============================================================================
class IAccelerator {
public:
    virtual ~IAccelerator() = default;
    
    /// @brief Get backend type
    virtual BackendType getType() const = 0;
    
    /// @brief Query device capabilities
    virtual DeviceCapabilities getCapabilities() const = 0;
    
    /// @brief Initialise backend with framebuffer dimensions
    /// @return true on success
    virtual bool initialise(int width, int height) = 0;
    
    /// @brief Check if backend is ready
    virtual bool isInitialised() const = 0;
    
    /// @brief Configure and launch render kernel
    virtual void launch(const LaunchConfig& config) = 0;
    
    /// @brief Get pointer to framebuffer (RGBA float)
    /// @note Caller does NOT own this memory
    virtual float* getFrameBuffer() = 0;
    
    /// @brief Get framebuffer size
    virtual size_t getFrameBufferSize() const = 0;
    
    /// @brief Upload background texture (starfield)
    virtual bool uploadBackground(const uint8_t* data, int width, int height) = 0;
    
    /// @brief Synchronise and wait for GPU completion
    virtual void synchronise() = 0;
    
    /// @brief Release all resources
    virtual void cleanup() = 0;

    /// @brief Get last error message
    virtual std::string getLastError() const = 0;

    /// @brief Reset accumulation buffer (call when camera changes or starting fresh render)
    virtual void resetAccumulation() {}
};

//==============================================================================
// Factory Function (implemented per backend)
//==============================================================================
std::unique_ptr<IAccelerator> createAccelerator(BackendType type);

/// @brief Query available backends on this system
std::vector<BackendType> getAvailableBackends();

/// @brief Get best available backend
BackendType getBestBackend();

} // namespace Sirius::Acceleration
