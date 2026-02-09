// RDRT001A.h - Renderer with OptiX Backend
//
// Relativistic ray tracing via null geodesic integration.
// Maps ray termination to colors: background (escaped), black (horizon),
// or emission (accretion disk).
//
// Backend: NVIDIA OptiX 7+ with CUDA acceleration

#pragma once

#include <string>
#include <memory>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <MTTN001A.h>

// Forward declarations
namespace Sirius { class IMetric; }
using namespace Sirius;

#ifdef SIRIUS_HAS_OPTIX
#include <cuda_runtime.h>
#include "RDOP003A.h"
#endif

class Renderer {
public:
    Renderer();
    ~Renderer();

    void init(int width, int height);
    void resize(int width, int height);
    void render(IMetric* metric, const Vec4& observerPos, const Vec4& observerVel, 
                float cameraYaw, float cameraPitch, float cameraFOV);
    void cleanup();
    
    // Background texture control
    bool loadBackgroundTexture(const std::string& path);
    void setUseBackgroundTexture(bool use) { m_UseBackgroundTexture = use; }
    bool getUseBackgroundTexture() const { return m_UseBackgroundTexture; }
    
    // Backend info
    bool isOptiXEnabled() const { return m_OptixEnabled; }
    
    // Lens flare control (Phase 6.5 - DNGR Cinematic Mode)
    void setLensFlareEnabled(bool enabled) { m_LensFlareEnabled = enabled; }
    bool isLensFlareEnabled() const { return m_LensFlareEnabled; }
    void setLensFlareIntensity(float intensity) { m_LensFlareIntensity = intensity; }
    void setLensFlareThreshold(float threshold) { m_LensFlareThreshold = threshold; }
    
    // Bloom control (Phase 7 - Cinematic Visual Quality)
    void setBloomEnabled(bool enabled) { m_BloomEnabled = enabled; }
    bool isBloomEnabled() const { return m_BloomEnabled; }
    void setBloomIntensity(float intensity) { m_BloomIntensity = intensity; }
    void setBloomThreshold(float threshold) { m_BloomThreshold = threshold; }
    
    // Disk mode control (Phase 7.5 - Planar vs Volumetric)
    void setDiskModePlanar(bool planar) { m_DiskModePlanar = planar; m_ResetAccumulation = true; }
    bool isDiskModePlanar() const { return m_DiskModePlanar; }
    
    // Ray bundle is always enabled (DNGR methodology - Phase 6.3)
    // Removed toggle: bundling is now default for all curved spacetimes

    // Explicit reset
    void resetAccumulation() { m_ResetAccumulation = true; }
    
    // Getters
    int getWidth() const { return m_Width; }
    int getHeight() const { return m_Height; }

private:
    void loadScreenShader(const std::string& vertPath, const std::string& fragPath);
    void loadLensFlareShader(const std::string& vertPath, const std::string& fragPath);
    void loadBloomShader(const std::string& vertPath, const std::string& fragPath);
    
#ifdef SIRIUS_HAS_OPTIX
    void renderOptiX(IMetric* metric, const Vec4& observerPos, const Vec4& observerVel,
                     float cameraYaw, float cameraPitch, float cameraFOV);
#endif

    // OpenGL objects (for display)
    GLuint m_ScreenProgram = 0;
    GLuint m_LensFlareProgram = 0;  // Phase 6.5 - Lens flare post-process
    GLuint m_BloomProgram = 0;      // Phase 7 - Bloom/glow post-process
    GLuint m_Texture = 0;
    GLuint m_Vao = 0;
    GLuint m_Vbo = 0;
    
    // Lens flare control (Phase 6.5)
    bool m_LensFlareEnabled = false;
    float m_LensFlareIntensity = 0.3f;
    float m_LensFlareThreshold = 0.8f;
    // Ray bundling always enabled - no toggle needed
    bool m_ResetAccumulation = false; // Trigger for accumulation reset
    
    // Bloom control (Phase 7)
    bool m_BloomEnabled = true;
    float m_BloomIntensity = 0.5f;
    float m_BloomThreshold = 0.7f;
    
    // Disk mode control (Phase 7.5)
    bool m_DiskModePlanar = true;  // true = Planar (default), false = Volumetric
    
    // Background texture
    GLuint m_BackgroundTexture = 0;
    int m_BackgroundWidth = 0;
    int m_BackgroundHeight = 0;
    bool m_UseBackgroundTexture = false;
    
    // Screen dimensions
    int m_Width = 0;
    int m_Height = 0;
    
    // OptiX backend
    void* m_OptixHandle = nullptr;
    bool m_OptixEnabled = false;
    int m_FrameCount = 0;           // Resets on camera movement (for accumulation)
    int m_AnimationFrame = 0;       // Never resets (for disk animation)
    int m_CurrentLoadedFrame = -1;
    
    // Camera state tracking for accumulation reset
    float m_LastCameraX = 0.0f;
    float m_LastCameraY = 0.0f;
    float m_LastCameraZ = 0.0f;
    float m_LastCameraYaw = 0.0f;
    float m_LastCameraPitch = 0.0f;
};