// IRRJ001A.h - Render Job
//
// Defines a single render task with configuration, execution, and progress tracking.
// Part of the static rendering architecture (replaces real-time loop).

#pragma once

#include <string>
#include <atomic>
#include <functional>
#include <vector>

namespace Sirius {

//==============================================================================
// Render Configuration
//==============================================================================
struct RenderConfig {
    // Scene parameters
    std::string metricName = "Schwarzschild";
    double M = 1.0;           // Black hole mass (geometric units)
    double a = 0.0;           // Spin parameter (|a| <= M)
    
    // Camera position (Boyer-Lindquist coordinates: t, r, θ, φ)
    double observerPosition[4] = {0.0, 50.0, 1.5708, 0.0};
    double observerVelocity[4] = {1.0, 0.0, 0.0, 0.0};
    float fov = 60.0f;        // Field of view (degrees)
    float yaw = 0.0f;         // Camera yaw (radians)
    float pitch = 0.0f;       // Camera pitch (radians)
    
    // Lens model
    enum class LensType { Pinhole, ThinLens, Fisheye, Equirectangular };
    LensType lens = LensType::Pinhole;
    float aperture = 0.0f;    // For ThinLens (0 = pinhole)
    float focusDistance = 10.0f;
    
    // Image dimensions
    int width = 1920;
    int height = 1080;
    
    // Quality parameters
    int samplesPerPixel = 256;
    int maxBounces = 1;       // 1-bounce default (per confirmed decision)
    double integrationTolerance = 1e-8;
    int maxIntegrationSteps = 10000;
    
    // Accretion disk
    bool enableDisk = true;
    double diskInnerRadius = 6.0;   // ISCO for Schwarzschild
    double diskOuterRadius = 20.0;
    double diskMdot = 1e-4;         // Accretion rate
    
    // Post-processing (applied during output)
    bool enableBloom = true;
    float bloomIntensity = 1.0f;
    float bloomThreshold = 1.0f;
    bool enableTonemapping = true;
    enum class TonemapOperator { ACES, Reinhard, Filmic };
    TonemapOperator tonemapper = TonemapOperator::ACES;
    float exposure = 1.0f;
    
    // Output
    std::string outputPath = "render.exr";
    bool enablePreview = false;
};

//==============================================================================
// Render Job Status
//==============================================================================
enum class JobStatus {
    Pending,
    Running,
    Paused,
    Completed,
    Failed,
    Cancelled
};

//==============================================================================
// Render Job
//==============================================================================
class RenderJob {
public:
    using ProgressCallback = std::function<void(float progress, int samplesComplete, int samplesTotal)>;
    using CompletionCallback = std::function<void(JobStatus status, const std::string& message)>;
    
    explicit RenderJob(const RenderConfig& config);
    ~RenderJob();
    
    // Non-copyable
    RenderJob(const RenderJob&) = delete;
    RenderJob& operator=(const RenderJob&) = delete;
    
    // Execution control
    void execute();
    void pause();
    void resume();
    void cancel();
    
    // Status queries
    JobStatus getStatus() const { return m_Status.load(); }
    float getProgress() const { return m_Progress.load(); }
    bool isComplete() const;
    bool isRunning() const;
    
    // Callbacks
    void setProgressCallback(ProgressCallback callback);
    void setCompletionCallback(CompletionCallback callback);
    
    // Configuration access
    const RenderConfig& getConfig() const { return m_Config; }
    
    // Statistics
    double getElapsedTime() const;
    double getEstimatedTimeRemaining() const;
    int getSamplesComplete() const { return m_SamplesComplete.load(); }
    
private:
    void renderScanline(int y);
    void renderScanlineBlackHole(int y);
    void applyPostProcessing();
    void writeOutput();
    void writePPM(const std::string& path);
    void writeEXR(const std::string& path);
    void reportProgress();
    
    RenderConfig m_Config;
    
    // Thread-safe status
    std::atomic<JobStatus> m_Status{JobStatus::Pending};
    std::atomic<float> m_Progress{0.0f};
    std::atomic<int> m_SamplesComplete{0};
    std::atomic<bool> m_PauseRequested{false};
    
    // Timing
    double m_StartTime = 0.0;
    double m_EndTime = 0.0;
    
    // Callbacks
    ProgressCallback m_ProgressCallback;
    CompletionCallback m_CompletionCallback;
    
    // Image buffer
    std::vector<float> m_ImageBuffer;  // RGBA float buffer
};

} // namespace Sirius
