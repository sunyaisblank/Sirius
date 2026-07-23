#pragma once

// FSM-based CPU render session: owns tile scheduling, progress tracking, the
// display buffer, and the physics components (metric, tracer, camera, jets).
// Ported from SNRS001A.h.
//
// The GPU/OptiX backend that legacy drove from here is removed: OptiX is retired
// and the Vulkan compute path arrives later through sirius::backend::device. The
// removal sites are marked in the implementation. The CPU render path is
// preserved byte-for-byte, including the spiral tile order, the thread-pool
// structure, and the per-pixel shading sequence.

#include "sirius/render/render_config.h"
#include "sirius/render/session/display_buffer.h"
#include "sirius/render/session/progress_tracker.h"
#include "sirius/render/session/session_events.h"
#include "sirius/render/session/session_states.h"
#include "sirius/render/session/state_machine.h"
#include "sirius/render/session/tile_scheduler.h"

// Physics integration.
#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera.h"
#include "sirius/core/disk/novikov_thorne_disk.h"       // ISCO computation.
#include "sirius/core/jet.h"                            // Relativistic jet model.
#include "sirius/core/metrics/astrophysical_scaling.h"  // SMBH parameters.
#include "sirius/core/metrics/kerr_schild_family.h"     // Kerr-Schild metric family.
#include "sirius/core/metrics/registry.h"               // Metric identity registry.
#include "sirius/core/metrics/warp_drive_family.h"      // Alcubierre warp drive family.
#include "sirius/core/polarisation/stokes.h"            // Stokes polarisation.
#include "sirius/core/postprocess.h"                    // Tonemapping types.
#include "sirius/core/spectral/colour_modes.h"          // Colour modes.
#include "sirius/core/starfield.h"                      // Starfield configuration.

#include <atomic>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace sirius::render {

// Compile-time transition table for the session DFA.
inline constexpr std::array<Transition<SessionState, SessionEvent>, 15> kSessionTransitions = {{
    // From Idle.
    {SessionState::Idle, SessionEvent::Start, SessionState::Initialising},

    // From Initialising.
    {SessionState::Initialising, SessionEvent::Ready, SessionState::Scheduling},
    {SessionState::Initialising, SessionEvent::Error, SessionState::Failed},

    // From Scheduling.
    {SessionState::Scheduling, SessionEvent::TileAvailable, SessionState::Rendering},
    {SessionState::Scheduling, SessionEvent::AllTilesComplete, SessionState::Completing},
    // A scheduling-time decline (a metric neither backend can render, or the
    // Vulkan render path declining before dispatch) fails cleanly with its
    // message rather than stranding the session mid-schedule.
    {SessionState::Scheduling, SessionEvent::Error, SessionState::Failed},

    // From Rendering.
    {SessionState::Rendering, SessionEvent::TileComplete, SessionState::Scheduling},
    {SessionState::Rendering, SessionEvent::Pause, SessionState::Paused},
    {SessionState::Rendering, SessionEvent::Cancel, SessionState::Cancelled},
    {SessionState::Rendering, SessionEvent::Error, SessionState::Failed},

    // From Paused.
    {SessionState::Paused, SessionEvent::Resume, SessionState::Rendering},
    {SessionState::Paused, SessionEvent::Cancel, SessionState::Cancelled},

    // From Completing.
    {SessionState::Completing, SessionEvent::OutputWritten, SessionState::Complete},
    {SessionState::Completing, SessionEvent::Error, SessionState::Failed},

    // Cancel available from Scheduling.
    {SessionState::Scheduling, SessionEvent::Cancel, SessionState::Cancelled},
}};

inline constexpr StateMachineConfig<SessionState, SessionEvent, 15> kSessionConfig = {
    SessionState::Idle, kSessionTransitions};

// Which backend the session drives. A two-value enum, not a bool, so the
// selection reads at the call site (docs/STYLE.md section 1). The Vulkan path
// dispatches the Slang trace kernel through sirius::backend::ComputeDevice; the
// CPU path is the reference tracer and is unchanged.
enum class RenderBackend {
    Cpu,     // the reference geodesic tracer (docs/ARCHITECTURE.md section 4)
    Vulkan,  // the Slang compute kernel on a ComputeDevice
};

// Configuration consumed by the render session.
struct SessionConfig {
    int width = 1920;
    int height = 1080;
    int tileSize = 64;
    int samplesPerPixel = 64;
    int threadCount = 0;                  // 0 = auto-detect, 1 = single-threaded.
    bool enableParallelRendering = true;  // Multi-threaded tile rendering.

    // Backend selection. `backend` decides the render path; useGPU/ptxPath are
    // legacy config-surface fields retained for parity. Since the go-live flip
    // (owner decision, 2026-07-18) `auto` in the config resolves to Vulkan
    // when a device is present and the registry marks the metric
    // gpu-dispatchable, falling back to Cpu with the reason logged;
    // `vulkan` selects the Vulkan path unconditionally and `cpu` pins the CPU
    // path. SIRIUS_RENDER_BACKEND overrides with the same values.
    RenderBackend backend = RenderBackend::Cpu;
    bool useGPU = true;
    std::string ptxPath;
    std::string outputPath = "render.ppm";

    // Spacetime identity and parameters; the id comes from the core registry.
    core::MetricId metricId = core::MetricId::Schwarzschild;
    double blackHoleMass = 1.0;
    double blackHoleSpin = 0.0;
    double blackHoleCharge = 0.0;       // Q/M (Reissner-Nordstrom, Kerr-Newman).
    double cosmologicalConstant = 0.0;  // Lambda (de Sitter family, spin = 0 only).
    double observerDistance = 50.0;
    double observerInclination = 1.5708;  // 90 degrees.
    float cameraFOV = 60.0f;

    // Camera four-velocity (P5): spatial beta in the local ray-component frame
    // (forward, up, right) = direction indices (1, 2, 3). Zero is a static
    // camera; the pinned render never sets these.
    double cameraBetaForward = 0.0;
    double cameraBetaUp = 0.0;
    double cameraBetaRight = 0.0;

    // Disk temperature model.
    std::string temperatureModel = "NovikovThorne";  // NovikovThorne or ShakuraSunyaev.
    float diskTemperatureScale = 50000.0f;           // T_scale (Kelvin).

    // Doppler beaming toggle (P4). True (default) keeps the full disk physics and
    // the pinned render; false suppresses the approaching/receding asymmetry.
    bool dopplerBeaming = true;

    // Exotic metric parameters.
    double throatRadius = 1.0;  // Morris-Thorne b0.
    double warpVelocity = 0.5;  // Alcubierre vs.
    double bubbleRadius = 1.0;  // Alcubierre R.
    double bubbleSigma = 0.5;   // Alcubierre sigma.

    // Post-processing (cinematic defaults).
    core::TonemapType tonemapper = core::TonemapType::Aces;  // Display transform (PPM/PNG).
    bool enableBloom = true;
    float bloomIntensity = 0.5f;
    float bloomThreshold = 0.3f;
    float exposure = 3.0f;
    float contrast = 1.1f;
    float saturation = 1.15f;

    // Motion blur (disk rotation).
    bool enableMotionBlur = false;
    float shutterTime = 0.1f;
    int motionBlurSamples = 3;

    // Volumetric disk.
    bool enableVolumetricDisk = false;
    float volumetricHOverR = 0.1f;
    float volumetricHPower = 0.25f;
    float volumetricTauMidplane = 10.0f;
    int volumetricSamples = 32;

    // Relativistic jets.
    bool enableJets = false;
    float jetLorentzFactor = 5.0f;
    float jetOpeningAngle = 0.1f;
    float jetLaunchRadius = 3.0f;
    float jetMaxExtent = 200.0f;
    float jetCollimation = 0.5f;
    float jetSpectralIndex = 2.2f;
    float jetIntensity = 1.0f;

    // Astronomical colouring mode.
    using ColorMode = core::color_modes::Mode;
    ColorMode colorMode = core::color_modes::Mode::TrueColor;

    // Polarisation output.
    bool enablePolarisation = false;
    bool outputPolarisationMap = false;

    // SMBH astrophysical scaling.
    core::SmbhParams smbhParams;

    // Turbulence and corona (GPU volumetric features; inert on the CPU path).
    bool enableTurbulence = false;
    bool enableCorona = false;

    // Depth-resolved starfield (GPU feature; inert on the CPU path).
    bool enableStarfield = false;
    core::StarfieldConfig starfieldConfig;

    // Filtered point-source star field (P3): render catalogue stars through the
    // beam footprint instead of the equirectangular texture. Default false keeps
    // the texture path and the pinned render.
    bool pointStarfield = false;

    // Ray bundles (P2): propagate geodesic deviation on the live path and derive
    // per-pixel beam footprints. Default false keeps the point-sampled path.
    bool rayBundles = false;

    // IMAX 70mm film simulation.
    bool enableFilmSimulation = false;
    FilmConfig filmConfig;

    // Build a SessionConfig from the unified render-config schema.
    static SessionConfig FromSiriusConfig(const SiriusConfig& config);
};

// Orchestrates a CPU render from configuration to written output.
class RenderSession {
  public:
    using FSM = StateMachine<SessionState, SessionEvent, 15>;
    using CompletionCallback =
        std::function<void(SessionState finalState, const std::string& message)>;

    RenderSession() : fsm_(kSessionConfig) { SetupActions(); }

    // Configure the session (must be in the Idle state).
    bool Configure(const SessionConfig& config) {
        if (fsm_.GetState() != SessionState::Idle) {
            return false;
        }
        config_ = config;
        return true;
    }

    bool Start() { return fsm_.Process(SessionEvent::Start); }
    bool Pause() { return fsm_.Process(SessionEvent::Pause); }
    bool Resume() { return fsm_.Process(SessionEvent::Resume); }

    bool Cancel() {
        progress_.GetCancellationToken().Cancel();
        return fsm_.Process(SessionEvent::Cancel);
    }

    // Block until the render thread (if any) finishes.
    void WaitForCompletion() {
        if (render_thread_.joinable()) {
            render_thread_.join();
        }
    }

    // Run synchronously and return the final state.
    SessionState Execute() {
        if (!Start()) {
            return SessionState::Failed;
        }
        WaitForCompletion();
        return fsm_.GetState();
    }

    SessionState GetState() const { return fsm_.GetState(); }
    bool IsComplete() const { return IsTerminal(fsm_.GetState()); }
    bool IsRunning() const { return IsActive(fsm_.GetState()); }

    float GetProgress() const { return progress_.GetProgress(); }
    double GetEta() const { return progress_.GetEta(); }
    const TileScheduler& GetTileScheduler() const { return tiles_; }
    DisplayBuffer& GetDisplayBuffer() { return display_; }

    void SetCompletionCallback(CompletionCallback cb) { completion_callback_ = cb; }
    void SetProgressCallback(ProgressCallback cb) { progress_.SetCallback(cb); }

  private:
    void SetupActions();
    void OnEnterState(SessionState state);

    void Initialise();
    void ScheduleNextTile();
    void RenderTile(Tile* tile);
    // Vulkan render path: dispatches the trace kernel per governed tile and
    // fills the display buffer, then fires AllTilesComplete so WriteOutput
    // applies the host display pipeline. Declines loudly (Error) when the
    // backend is absent, the metric is not on the Vulkan render path, or the
    // device cannot be opened.
    void RenderVulkanPath();
    void WriteOutput();
    void OnSessionEnd(SessionState state);

    // Shaded result of a single ray sample.
    struct PixelResult {
        float r = 0.0f, g = 0.0f, b = 0.0f;
        core::StokesVector stokes;
    };

    PixelResult ShadePixel(int px, int py, backend::GeodesicTracer* tracer) const;
    PixelResult ShadeDiskHit(const backend::TraceResult& result) const;
    PixelResult ShadeEscaped(const backend::TraceResult& result) const;

    FSM fsm_;
    SessionConfig config_;
    TileScheduler tiles_;
    ProgressTracker progress_;
    DisplayBuffer display_;
    std::thread render_thread_;
    std::string error_message_;
    CompletionCallback completion_callback_;

    // Physics components. metric_ is null when the CPU path cannot represent the
    // requested spacetime; the scheduler then refuses rather than substituting.
    std::unique_ptr<core::IMetric> metric_;
    std::unique_ptr<backend::GeodesicTracer> tracer_;
    std::unique_ptr<core::PinholeCamera> camera_;

    // Relativistic jet model.
    std::unique_ptr<core::RelativisticJet> jet_;

    // Polarisation output buffer.
    std::vector<core::StokesVector> polarisation_buffer_;

    // Starfield background texture (equirectangular RGBA).
    std::vector<unsigned char> starfield_data_;
    int starfield_width_ = 0;
    int starfield_height_ = 0;
    bool starfield_loaded_ = false;

    bool LoadStarfieldTexture(const std::string& path);
    void SampleStarfield(const core::Vec4& direction, float& r, float& g, float& b) const;

    // Filtered point-source star field (P3): accumulate catalogue stars through
    // the per-ray beam footprint (ray bundles on) or a pinhole sigma (bundles
    // off, the flicker baseline).
    std::vector<core::StarEntry> star_catalogue_;
    std::unique_ptr<core::StarfieldGenerator> star_generator_;
    double pixel_angular_size_ = 0.0;  // fov / height, radians.
    void SampleStarfieldPoints(const backend::TraceResult& result, float& r, float& g,
                               float& b) const;

    // Multi-threaded tile rendering: each worker has its own tracer sharing the
    // metric, so no synchronisation is needed on the physics path.
    void RenderTilesParallel();
    void WorkerThread(int thread_id);
    void RenderTileThreaded(Tile* tile, int thread_id);

    std::vector<std::thread> worker_threads_;
    std::vector<std::unique_ptr<backend::GeodesicTracer>> thread_tracers_;  // Per-thread tracers.
    std::mutex tile_mutex_;                  // Protects tile acquisition.
    std::mutex display_mutex_;               // Protects display buffer updates.
    std::atomic<bool> stop_workers_{false};  // Signal workers to stop.
    std::atomic<int> active_workers_{0};     // Workers currently rendering.
    int num_threads_ = 1;                    // Actual thread count.
};

}  // namespace sirius::render
