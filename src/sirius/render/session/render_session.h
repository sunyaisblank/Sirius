#pragma once

// FSM-based CPU render session: owns tile scheduling, progress tracking, the
// display buffer, and the physics components (metric, tracer, camera, jets).
// Ported from SNRS001A.h.
//
// OptiX is retired. Vulkan enters through sirius::backend::device behind the
// same asynchronous, cancellable session lifecycle as the CPU reference path.

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
#include "sirius/core/disk/novikov_thorne_disk.h"    // ISCO computation.
#include "sirius/core/jet.h"                         // Relativistic jet model.
#include "sirius/core/metrics/kerr_schild_family.h"  // Kerr-Schild metric family.
#include "sirius/core/metrics/registry.h"            // Metric identity registry.
#include "sirius/core/metrics/warp_drive_family.h"   // Alcubierre warp drive family.
#include "sirius/core/postprocess.h"                 // Tonemapping types.
#include "sirius/core/spectral/colour_modes.h"       // Colour modes.
#include "sirius/core/starfield.h"                   // Starfield configuration.

#include <atomic>
#include <condition_variable>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace sirius::render {

// Compile-time transition table for the session DFA.
inline constexpr std::array<Transition<SessionState, SessionEvent>, 14> kSessionTransitions = {{
    // From Idle.
    {SessionState::Idle, SessionEvent::Start, SessionState::Initialising},

    // From Initialising.
    {SessionState::Initialising, SessionEvent::Ready, SessionState::Scheduling},
    {SessionState::Initialising, SessionEvent::Error, SessionState::Failed},
    {SessionState::Initialising, SessionEvent::Cancel, SessionState::Cancelled},

    // From Scheduling.
    {SessionState::Scheduling, SessionEvent::TileAvailable, SessionState::Rendering},
    {SessionState::Scheduling, SessionEvent::AllTilesComplete, SessionState::Completing},
    // A scheduling-time decline (a metric neither backend can render, or the
    // Vulkan render path declining before dispatch) fails cleanly with its
    // message rather than stranding the session mid-schedule.
    {SessionState::Scheduling, SessionEvent::Error, SessionState::Failed},

    // From Rendering.
    {SessionState::Rendering, SessionEvent::TileComplete, SessionState::Scheduling},
    {SessionState::Rendering, SessionEvent::Cancel, SessionState::Cancelled},
    {SessionState::Rendering, SessionEvent::Error, SessionState::Failed},

    // From Completing.
    {SessionState::Completing, SessionEvent::OutputWritten, SessionState::Complete},
    {SessionState::Completing, SessionEvent::Error, SessionState::Failed},
    {SessionState::Completing, SessionEvent::Cancel, SessionState::Cancelled},

    // Cancel available from Scheduling.
    {SessionState::Scheduling, SessionEvent::Cancel, SessionState::Cancelled},
}};

inline constexpr StateMachineConfig<SessionState, SessionEvent, 14> kSessionConfig = {
    SessionState::Idle, kSessionTransitions};

// Which backend the session drives. A two-value enum, not a bool, so the
// selection reads at the call site (docs/STYLE.md section 1). The Vulkan path
// dispatches the Slang trace kernel through sirius::backend::ComputeDevice; the
// CPU path is the reference tracer and is unchanged.
enum class RenderBackend {
    Cpu,     // the reference geodesic tracer (docs/ARCHITECTURE.md section 4)
    Vulkan,  // the Slang compute kernel on a ComputeDevice
};

enum class WormholeTopology {
    OneSheetCapture,
    TwoSheet,
};

// Configuration consumed by the render session.
struct SessionConfig {
    int width = 1920;
    int height = 1080;
    int tileSize = 64;
    int samplesPerPixel = 64;
    int threadCount = 0;                  // 0 = auto-detect, 1 = single-threaded.
    bool enableParallelRendering = true;  // Multi-threaded tile rendering.
    bool writeOutput = true;              // False for in-memory progressive previews.

    // Backend selection. Since the go-live flip (owner decision, 2026-07-18)
    // `auto` in the external config resolves to Vulkan when a device is present
    // and the complete scene is represented, falling back to Cpu with the
    // reason logged; `vulkan` selects Vulkan unconditionally and `cpu` pins CPU.
    RenderBackend backend = RenderBackend::Cpu;
    std::string outputPath = "render.ppm";

    // Spacetime identity and parameters; the id comes from the core registry.
    core::MetricId metricId = core::MetricId::Schwarzschild;
    double blackHoleMass = 1.0;
    double blackHoleSpin = 0.0;
    double blackHoleCharge = 0.0;       // Q/M (Reissner-Nordstrom, Kerr-Newman).
    double cosmologicalConstant = 0.0;  // Lambda (de Sitter family, spin = 0 only).
    double observerDistance = 50.0;
    double observerInclination = 1.5708;  // 90 degrees.
    double observerAzimuth = 0.0;
    float cameraFOV = 60.0f;

    // Camera four-velocity (P5): spatial beta in the local ray-component frame
    // (forward, up, right) = direction indices (1, 2, 3). Zero is a static
    // camera; the pinned render never sets these.
    double cameraBetaForward = 0.0;
    double cameraBetaUp = 0.0;
    double cameraBetaRight = 0.0;
    core::LensType lensType = core::LensType::Pinhole;
    float cameraFocalLength = 50.0f;
    float cameraAperture = 2.8f;
    float cameraFocusDistance = 50.0f;

    // Disk temperature model.
    DiskTemperatureModel temperatureModel = DiskTemperatureModel::NovikovThorne;
    float diskTemperatureScale = 50000.0f;  // T_scale (Kelvin).
    bool enableDisk = true;

    // Doppler beaming toggle (P4). True (default) keeps the full disk physics and
    // the pinned render; false suppresses the approaching/receding asymmetry.
    bool dopplerBeaming = true;

    // Exotic metric parameters.
    double throatRadius = 1.0;  // Morris-Thorne b0.
    WormholeTopology wormholeTopology = WormholeTopology::OneSheetCapture;
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

    // Polarisation transport is a strict companion to ColorMode::Polarisation.
    // Typed callers must set both coherently; the external schema derives this
    // bool from its single colorMode field.
    bool enablePolarisation = false;

    // Optional deterministic density turbulence and inverse-Compton corona
    // contributions on both live volumetric transfer paths.
    bool enableTurbulence = false;
    bool enableCorona = false;

    // Depth-resolved starfield catalogue parameters.
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

// Typed boundary shared by CPU initialisation and Vulkan capability selection.
// Small positive dimensions remain legal for probes, while production limits
// and every enum/feature dependency are enforced before allocation or dispatch.
[[nodiscard]] std::optional<std::string> SessionConfigIssue(const SessionConfig& config);

// Orchestrates a CPU render from configuration to written output.
class RenderSession {
  public:
    using FSM = StateMachine<SessionState, SessionEvent, 14>;
    using CompletionCallback =
        std::function<void(SessionState finalState, const std::string& message)>;

    RenderSession() : fsm_(kSessionConfig) { SetupActions(); }
    ~RenderSession();

    RenderSession(const RenderSession&) = delete;
    RenderSession& operator=(const RenderSession&) = delete;
    RenderSession(RenderSession&&) = delete;
    RenderSession& operator=(RenderSession&&) = delete;

    // Configure the session (must be in the Idle state).
    bool Configure(const SessionConfig& config);

    // Launch asynchronously. Execute() is the synchronous convenience wrapper.
    bool Start();
    bool Cancel();

    // Block until the render thread (if any) finishes.
    void WaitForCompletion();

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

    void SetCompletionCallback(CompletionCallback cb) {
        std::lock_guard<std::mutex> lock(callback_mutex_);
        completion_callback_ = std::move(cb);
    }
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
    mutable std::mutex callback_mutex_;

    // Physics components. metric_ is null when the CPU path cannot represent the
    // requested spacetime; the scheduler then refuses rather than substituting.
    std::unique_ptr<core::IMetric> metric_;
    std::unique_ptr<backend::GeodesicTracer> tracer_;
    std::unique_ptr<core::ICamera> camera_;

    // Relativistic jet model.
    std::unique_ptr<core::RelativisticJet> jet_;

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
    std::unique_ptr<core::StarfieldSpatialIndex> star_index_;
    double pixel_angular_size_ = 0.0;  // fov / height, radians.
    void SampleStarfieldPoints(const backend::TraceResult& result, float& r, float& g,
                               float& b) const;

    // Multi-threaded tile rendering: each worker has its own tracer sharing the
    // metric, so no synchronisation is needed on the physics path.
    void RenderTilesParallel();
    void WorkerThread(int thread_id);
    [[nodiscard]] bool RenderTileThreaded(Tile* tile, int thread_id);

    std::vector<std::thread> worker_threads_;
    std::vector<std::unique_ptr<backend::GeodesicTracer>> thread_tracers_;  // Per-thread tracers.
    std::mutex tile_mutex_;                  // Protects tile acquisition.
    std::mutex display_mutex_;               // Protects display buffer updates.
    std::atomic<bool> stop_workers_{false};  // Signal workers to stop.
    std::atomic<int> active_workers_{0};     // Workers currently rendering.
    int num_threads_ = 1;                    // Actual thread count.
    mutable std::mutex lifecycle_mutex_;
    std::condition_variable lifecycle_cv_;
    bool join_in_progress_ = false;
    std::thread::id render_thread_id_;
};

}  // namespace sirius::render
