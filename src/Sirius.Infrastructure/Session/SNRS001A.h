// SNRS001A.h - Render Session
// Component ID: SNRS001A (Session/RenderSession)
//
// FSM-based render session orchestrator. Manages tile scheduling,
// progress tracking, and display buffer updates.

#pragma once

#include "SMSM001A.h"
#include "SNST001A.h"
#include "SNEV001A.h"
#include "SNTL001A.h"
#include "SNPR001A.h"
#include "SNDP001A.h"
#include "IRRJ001A.h"

// Physics integration (new in Priority 1)
#include "PHMT100A.h"  // Kerr-Schild metric family
#include "CMBS001A.h"  // Camera interface
#include "GTRC001A.h"  // Geodesic tracer
#include "PHJT001A.h"  // Relativistic jet model
#include "PHPL001A.h"  // Stokes polarisation

// Cinematic Expansion (2026)
#include "../Sirius.Core/Metric/PHSM001A.h"      // SMBH parameters
#include "../Sirius.Core/Disk/PHTR001A.h"       // Turbulence model
#include "../Sirius.Core/Disk/PHCR001A.h"       // Corona model
#include "../Sirius.Core/Disk/PHAD003A.h"       // Volumetric disk
#include "../Sirius.Core/Jet/PHJT002A.h"        // MHD jet
#include "../Sirius.Core/Environment/PHSF001A.h" // Starfield
#include "../Sirius.Infrastructure/Configuration/CRFM001A.h" // Film config

// GPU acceleration (OptiX backend)
#include "../Sirius.Render/Acceleration/Backend/ACIB001A.h"

#include <thread>
#include <vector>
#include <functional>
#include <iostream>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <atomic>

namespace Sirius {

//==============================================================================
// Session Transition Table (Compile-time DFA)
//==============================================================================
constexpr std::array<Transition<SessionState, SessionEvent>, 14> SESSION_TRANSITIONS = {{
    // From Idle
    {SessionState::Idle, SessionEvent::Start, SessionState::Initialising},
    
    // From Initialising
    {SessionState::Initialising, SessionEvent::Ready, SessionState::Scheduling},
    {SessionState::Initialising, SessionEvent::Error, SessionState::Failed},
    
    // From Scheduling
    {SessionState::Scheduling, SessionEvent::TileAvailable, SessionState::Rendering},
    {SessionState::Scheduling, SessionEvent::AllTilesComplete, SessionState::Completing},
    
    // From Rendering
    {SessionState::Rendering, SessionEvent::TileComplete, SessionState::Scheduling},
    {SessionState::Rendering, SessionEvent::Pause, SessionState::Paused},
    {SessionState::Rendering, SessionEvent::Cancel, SessionState::Cancelled},
    {SessionState::Rendering, SessionEvent::Error, SessionState::Failed},
    
    // From Paused
    {SessionState::Paused, SessionEvent::Resume, SessionState::Rendering},
    {SessionState::Paused, SessionEvent::Cancel, SessionState::Cancelled},
    
    // From Completing
    {SessionState::Completing, SessionEvent::OutputWritten, SessionState::Complete},
    {SessionState::Completing, SessionEvent::Error, SessionState::Failed},
    
    // Cancel available from Scheduling
    {SessionState::Scheduling, SessionEvent::Cancel, SessionState::Cancelled}
}};

constexpr StateMachineConfig<SessionState, SessionEvent, 14> SESSION_CONFIG = {
    SessionState::Idle,
    SESSION_TRANSITIONS
};

//==============================================================================
// Render Session Configuration
//==============================================================================
struct SessionConfig {
    int width = 1920;
    int height = 1080;
    int tileSize = 64;
    int samplesPerPixel = 64;
    int threadCount = 0;    ///< 0 = auto-detect, 1 = single-threaded
    bool enableParallelRendering = true;  ///< Enable multi-threaded tile rendering

    // GPU Acceleration
    bool useGPU = true;     ///< Enable GPU acceleration (OptiX) when available
    std::string ptxPath;    ///< Path to PTX files (auto-detected if empty)
    std::string outputPath = "render.ppm";

    // Black hole parameters
    std::string metricName = "Schwarzschild";
    double blackHoleMass = 1.0;
    double blackHoleSpin = 0.0;
    double observerDistance = 50.0;
    double observerInclination = 1.5708;  // 90 degrees
    float cameraFOV = 60.0f;

    // Post-processing (cinematic defaults)
    bool enableBloom = true;
    float bloomIntensity = 0.5f;    ///< Bloom glow intensity
    float bloomThreshold = 0.3f;    ///< Lower threshold to capture disk glow
    float exposure = 3.0f;          ///< Higher exposure for visible disk
    float contrast = 1.1f;          ///< Slight contrast boost
    float saturation = 1.15f;       ///< Slight saturation boost for colors

    // Phase 3 P3: Motion blur (disk rotation)
    bool enableMotionBlur = false;      ///< Enable disk rotation motion blur
    float shutterTime = 0.1f;           ///< Shutter time in GM/c³ units (orbital period scale)
    int motionBlurSamples = 3;          ///< Number of temporal samples for motion blur

    // Phase 6: Volumetric disk
    bool enableVolumetricDisk = false;  ///< Enable 3D volumetric disk (vs thin disk)
    float volumetricHOverR = 0.1f;      ///< Scale height ratio H/r at reference radius
    float volumetricHPower = 0.25f;     ///< Flaring index: H/r ∝ r^H_power
    float volumetricTauMidplane = 10.0f; ///< Midplane optical depth at r_ref
    int volumetricSamples = 32;         ///< Number of ray marching samples

    // ==========================================================================
    // Relativistic Jets (Phase 7)
    // ==========================================================================
    bool enableJets = false;            ///< Enable bipolar relativistic jets
    float jetLorentzFactor = 5.0f;      ///< Bulk Lorentz factor Γ
    float jetOpeningAngle = 0.1f;       ///< Half-opening angle at launch [radians]
    float jetLaunchRadius = 3.0f;       ///< Jet launching radius [M]
    float jetMaxExtent = 200.0f;        ///< Maximum jet extent [M]
    float jetCollimation = 0.5f;        ///< Shape parameter (0=conical, 1=parabolic)
    float jetSpectralIndex = 2.2f;      ///< Electron power-law index p
    float jetIntensity = 1.0f;          ///< Jet emission intensity scale

    // ==========================================================================
    // Astronomical Coloring Modes (Phase 7)
    // ==========================================================================
    // Inspired by Hubble imaging pipeline:
    // - TrueColor: Physical blackbody → XYZ → sRGB (default)
    // - TemperatureMap: False color showing temperature distribution
    // - RedshiftMap: Visualize g-factor (Doppler + gravitational)
    // - Narrowband: Hubble palette mapping emission lines to RGB
    // - Polarisation: Show polarisation degree and EVPA
    // ==========================================================================
    enum class ColorMode {
        TrueColor,       ///< Physical blackbody colors (what a camera would see)
        TemperatureMap,  ///< False color temperature: cold→blue, hot→red
        RedshiftMap,     ///< g-factor visualization: redshift→red, blueshift→blue
        Narrowband,      ///< Hubble palette: map emission bands to RGB
        Polarisation     ///< Polarisation degree and angle visualization
    };
    ColorMode colorMode = ColorMode::TrueColor;

    // ==========================================================================
    // Polarisation Output (Phase 7)
    // ==========================================================================
    bool enablePolarisation = false;    ///< Track Stokes parameters during tracing
    bool outputPolarisationMap = false; ///< Write separate polarisation visualization

    // ==========================================================================
    // Cinematic Expansion (2026)
    // ==========================================================================

    // SMBH Astrophysical Parameters
    SMBHParams smbhParams;              ///< Supermassive black hole scaling

    // Enhanced Volumetric Disk
    bool enableVolumetricDiskAdvanced = false;  ///< Enable turbulent volumetric disk
    VolumetricDiskConfig volumetricDiskConfig;  ///< Turbulence, corona, full volumetric

    // MHD Relativistic Jet
    bool enableJetMHD = false;          ///< Enable MHD-based jet (replaces basic jet)
    JetMHDConfig jetMHDConfig;          ///< MHD jet parameters

    // Depth-Resolved Starfield
    bool enableStarfield = false;       ///< Enable parallax-correct starfield
    StarfieldConfig starfieldConfig;    ///< Starfield parameters

    // IMAX 70mm Film Simulation
    bool enableFilmSimulation = false;  ///< Enable film post-processing
    FilmConfig filmConfig;              ///< Film simulation parameters
};

//==============================================================================
// Render Session
//==============================================================================
class RenderSession {
public:
    using FSM = StateMachine<SessionState, SessionEvent, 14>;
    using CompletionCallback = std::function<void(SessionState finalState, const std::string& message)>;
    
    RenderSession() : m_FSM(SESSION_CONFIG) {
        setupActions();
    }
    
    /// @brief Configure session (must be in Idle state)
    bool configure(const SessionConfig& config) {
        if (m_FSM.getState() != SessionState::Idle) {
            return false;
        }
        m_Config = config;
        return true;
    }
    
    /// @brief Start the render session
    bool start() {
        return m_FSM.process(SessionEvent::Start);
    }
    
    /// @brief Pause rendering
    bool pause() {
        return m_FSM.process(SessionEvent::Pause);
    }
    
    /// @brief Resume from pause
    bool resume() {
        return m_FSM.process(SessionEvent::Resume);
    }
    
    /// @brief Cancel rendering
    bool cancel() {
        m_Progress.getCancellationToken().cancel();
        return m_FSM.process(SessionEvent::Cancel);
    }
    
    /// @brief Wait for session to complete (blocking)
    void waitForCompletion() {
        if (m_RenderThread.joinable()) {
            m_RenderThread.join();
        }
    }
    
    /// @brief Execute synchronously (blocking)
    SessionState execute() {
        if (!start()) {
            return SessionState::Failed;
        }
        waitForCompletion();
        return m_FSM.getState();
    }
    
    // State queries
    SessionState getState() const { return m_FSM.getState(); }
    bool isComplete() const { return isTerminal(m_FSM.getState()); }
    bool isRunning() const { return isActive(m_FSM.getState()); }
    
    // Progress access
    float getProgress() const { return m_Progress.getProgress(); }
    double getETA() const { return m_Progress.getETA(); }
    const TileScheduler& getTileScheduler() const { return m_Tiles; }
    DisplayBuffer& getDisplayBuffer() { return m_Display; }
    
    // Callbacks
    void setCompletionCallback(CompletionCallback cb) { m_CompletionCallback = cb; }
    void setProgressCallback(ProgressCallback cb) { m_Progress.setCallback(cb); }
    
private:
    void setupActions() {
        m_FSM.setEntryAction([this](SessionState state) {
            onEnterState(state);
        });
        
        m_FSM.setTransitionAction([this](SessionState from, SessionEvent event, SessionState to) {
            std::cout << "[Session] " << stateName(from) << " --" 
                      << eventName(event) << "--> " << stateName(to) << std::endl;
        });
    }
    
    void onEnterState(SessionState state) {
        switch (state) {
            case SessionState::Initialising:
                initialise();
                break;
            case SessionState::Scheduling:
                scheduleNextTile();
                break;
            case SessionState::Rendering:
                // Tile is being processed
                break;
            case SessionState::Completing:
                writeOutput();
                break;
            case SessionState::Complete:
            case SessionState::Failed:
            case SessionState::Cancelled:
                onSessionEnd(state);
                break;
            default:
                break;
        }
    }
    
    void initialise();
    void scheduleNextTile();
    void renderTile(Tile* tile);
    void writeOutput();
    void onSessionEnd(SessionState state);
    
    FSM m_FSM;
    SessionConfig m_Config;
    TileScheduler m_Tiles;
    ProgressTracker m_Progress;
    DisplayBuffer m_Display;
    std::thread m_RenderThread;
    std::string m_ErrorMessage;
    CompletionCallback m_CompletionCallback;

    // Physics integration components (Priority 1 fix)
    std::unique_ptr<KerrSchildFamily> m_Metric;
    std::unique_ptr<GeodesicTracer> m_Tracer;
    std::unique_ptr<PinholeCamera> m_Camera;

    // Relativistic jet model (Phase 7)
    std::unique_ptr<RelativisticJet> m_Jet;

    // Polarisation output buffer (Phase 7)
    std::vector<StokesVector> m_PolarisationBuffer;

    // GPU acceleration (OptiX backend)
    std::unique_ptr<Acceleration::IAccelerator> m_Accelerator;
    bool m_UseGPU = false;  ///< True if GPU rendering is active

    // Starfield texture (equirectangular projection)
    std::vector<unsigned char> m_StarfieldData;
    int m_StarfieldWidth = 0;
    int m_StarfieldHeight = 0;
    bool m_StarfieldLoaded = false;

    // Load starfield texture from file
    bool loadStarfieldTexture(const std::string& path);

    // Sample starfield texture at given direction (returns RGB)
    void sampleStarfield(const Vec4& direction, float& r, float& g, float& b) const;

    // =========================================================================
    // Multi-threaded Tile Rendering
    // =========================================================================
    // Thread pool for parallel tile rendering. Each worker has its own
    // tracer instance (sharing the metric) to avoid contention.
    // =========================================================================

    /// @brief Render tiles in parallel using a thread pool
    void renderTilesParallel();

    /// @brief Render using GPU (OptiX) backend
    void renderGPU();

    /// @brief Worker thread function for parallel rendering
    void workerThread(int threadId);

    /// @brief Render a single tile (thread-safe version)
    void renderTileThreaded(Tile* tile, int threadId);

    // Thread pool components
    std::vector<std::thread> m_WorkerThreads;
    std::vector<std::unique_ptr<GeodesicTracer>> m_ThreadTracers;  ///< Per-thread tracers
    std::mutex m_TileMutex;                 ///< Protects tile acquisition
    std::mutex m_DisplayMutex;              ///< Protects display buffer updates
    std::atomic<bool> m_StopWorkers{false}; ///< Signal workers to stop
    std::atomic<int> m_ActiveWorkers{0};    ///< Number of workers currently rendering
    int m_NumThreads = 1;                   ///< Actual thread count
};

} // namespace Sirius
