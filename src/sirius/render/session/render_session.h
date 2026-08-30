#pragma once

#include "sirius/base/error.h"

// FSM-based CPU render session: owns tile scheduling, progress tracking, the
// display buffer, and the physics components (metric, tracer, camera).
// Ported from SNRS001A.h.
//
// OptiX is retired. Vulkan enters through sirius::backend::device behind the
// same asynchronous, cancellable session lifecycle as the CPU reference path.

#include "sirius/render/film_config.h"
#include "sirius/render/session/display_buffer.h"
#include "sirius/render/session/progress_tracker.h"
#include "sirius/render/session/session_events.h"
#include "sirius/render/session/session_states.h"
#include "sirius/render/session/state_machine.h"
#include "sirius/render/session/tile_scheduler.h"

// Physics integration.
#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/core/camera.h"
#include "sirius/core/disk/disk_defaults.h"
#include "sirius/core/disk/novikov_thorne_disk.h"  // ISCO computation.
#include "sirius/core/feature_defaults.h"
#include "sirius/core/metrics/kerr_schild_family.h"  // Kerr-Schild metric family.
#include "sirius/core/metrics/registry.h"            // Metric identity registry.
#include "sirius/core/metrics/warp_drive_family.h"   // Alcubierre warp drive family.
#include "sirius/core/postprocess.h"                 // Tonemapping types.
#include "sirius/core/spectral/colour_modes.h"       // Colour modes.
#include "sirius/core/starfield.h"                   // Starfield configuration.

#include <atomic>
#include <condition_variable>
#include <cstddef>
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

enum class DiskTemperatureModel : std::uint8_t {
    NovikovThorne,
    ShakuraSunyaev,
};

// Configuration consumed by the render session.
struct SessionConfig {
    int width = 1920;
    int height = 1080;
    int tile_size = 64;
    int samples_per_pixel = 64;
    int thread_count = 0;                   // 0 = auto-detect, 1 = single-threaded.
    bool enable_parallel_rendering = true;  // Multi-threaded tile rendering.
    bool write_output = true;               // False for in-memory progressive previews.

    // Backend selection. Since the go-live flip (owner decision, 2026-07-18)
    // `auto` in the external config resolves to Vulkan when a device is present
    // and the complete scene is represented, falling back to Cpu with the
    // reason logged; `vulkan` selects Vulkan unconditionally and `cpu` pins CPU.
    RenderBackend backend = RenderBackend::Cpu;
    std::string output_path = "render.ppm";

    // Spacetime identity and parameters; the id comes from the core registry.
    core::MetricId metric_id = core::MetricId::Schwarzschild;
    double black_hole_mass = 1.0;  // Geometric coordinate length; zero if M is not a parameter.
    double black_hole_spin = 0.0;
    double black_hole_charge = 0.0;        // Q/M (Reissner-Nordstrom, Kerr-Newman).
    double cosmological_constant = 0.0;    // Lambda (de Sitter family, spin = 0 only).
    double observer_distance = 50.0;       // Coordinate radius r, not the ratio r/M.
    double observer_inclination = 1.5708;  // 90 degrees.
    double observer_azimuth = 0.0;
    float camera_fov = 60.0f;

    // Camera four-velocity (P5): spatial beta in the local ray-component frame
    // (forward, up, right) = direction indices (1, 2, 3). Zero is a static
    // camera; the pinned render never sets these.
    double camera_beta_forward = 0.0;
    double camera_beta_up = 0.0;
    double camera_beta_right = 0.0;
    core::LensType lens_type = core::LensType::Pinhole;
    float camera_focal_length = core::kDefaultCameraFocalLength;
    float camera_aperture = core::kDefaultCameraAperture;
    float camera_focus_distance =
        core::kDefaultCameraFocusDistance;  // Geometric coordinate length from the camera.

    // Disk temperature model.
    DiskTemperatureModel temperature_model = DiskTemperatureModel::NovikovThorne;
    // Effective temperature at 1.5 times the disk inner edge [K].
    float disk_temperature_scale = core::kDefaultDiskTemperatureKelvin;
    bool enable_disk = true;

    // Doppler beaming toggle (P4). True (default) keeps the full disk physics and
    // the pinned render; false suppresses the approaching/receding asymmetry.
    bool doppler_beaming = true;

    // Exotic metric parameters.
    double throat_radius = core::kDefaultMorrisThorneThroatRadius;  // Morris-Thorne b0.
    WormholeTopology wormhole_topology = WormholeTopology::OneSheetCapture;
    double warp_velocity = core::kDefaultAlcubierreWarpVelocity;  // Alcubierre vs.
    double bubble_radius = core::kDefaultAlcubierreBubbleRadius;  // Alcubierre R.
    double bubble_sigma = core::kDefaultAlcubierreBubbleSigma;    // Inverse wall length sigma.

    // Post-processing (cinematic defaults).
    core::TonemapType tonemapper = core::TonemapType::Aces;  // Display transform (PPM/PNG).
    bool enable_bloom = true;
    float bloom_intensity = core::kDefaultBloomIntensity;
    float bloom_threshold = core::kDefaultBloomThreshold;
    float exposure = 3.0f;
    float contrast = 1.1f;
    float saturation = 1.15f;

    // Motion blur (disk rotation).
    bool enable_motion_blur = false;
    float shutter_time = core::kDefaultMotionBlurShutterTime;
    int motion_blur_samples = core::kDefaultMotionBlurSamples;

    // Volumetric disk.
    bool enable_volumetric_disk = false;
    float volumetric_h_over_r = core::kDefaultVolumetricHOverR;
    float volumetric_h_power = core::kDefaultVolumetricHPower;
    float volumetric_tau_midplane = core::kDefaultVolumetricTauMidplane;
    int volumetric_samples = core::kDefaultVolumetricSamples;

    // Reserved typed request. The session rejects it until a covariant
    // geodesic radiative-transfer model is represented.
    bool enable_jets = false;

    // Astronomical colouring mode.
    using ColorMode = core::color_modes::Mode;
    ColorMode color_mode = core::color_modes::Mode::TrueColor;

    // Polarisation transport is a strict companion to ColorMode::Polarisation.
    // Typed callers must set both coherently; the external schema derives this
    // bool from its single colorMode field.
    bool enable_polarisation = false;

    // Optional deterministic procedural density modulation.
    bool enable_turbulence = false;
    // Reserved request: rejected until frequency-dependent covariant Compton
    // transfer is represented.
    bool enable_corona = false;

    // Filtered point-source star field (P3): render catalogue stars through the
    // beam footprint instead of the equirectangular texture. Default false keeps
    // the texture path and the pinned render.
    bool point_starfield = false;
    core::PointStarfieldConfig point_starfield_config;

    // Ray bundles (P2): propagate geodesic deviation on the live path and derive
    // per-pixel beam footprints. Default false keeps the point-sampled path.
    bool ray_bundles = false;

    // Optional bounded grain/halation/grade/vignette/bloom display finish.
    bool enable_film_finish = false;
    FilmConfig film_config = FilmConfig::Interstellar();
};

// Canonical, machine-readable witness emitted from the typed configuration
// that the session actually consumes. External attestation compares this event
// with its claims instead of trusting runbook metadata alone.
[[nodiscard]] std::string SessionSceneEvidenceJson(const SessionConfig& config,
                                                   std::size_t point_star_count);

// Typed boundary shared by CPU initialisation and Vulkan capability selection.
// Small positive dimensions remain legal for probes, while production limits
// and every enum/feature dependency are enforced before allocation or dispatch.
[[nodiscard]] std::optional<std::string> SessionConfigIssue(const SessionConfig& config);

// Orchestrates a CPU render from configuration to written output.
class RenderSession {
  public:
    using FSM = StateMachine<SessionState, SessionEvent, 14>;
    using CompletionCallback =
        std::function<void(SessionState final_state, const std::string& message)>;

    RenderSession() : fsm_(kSessionConfig) { SetupActions(); }
    ~RenderSession();

    RenderSession(const RenderSession&) = delete;
    RenderSession& operator=(const RenderSession&) = delete;
    RenderSession(RenderSession&&) = delete;
    RenderSession& operator=(RenderSession&&) = delete;

    // Configure and validate the session (must be in the Idle state).
    [[nodiscard]] base::Expected<void> Configure(const SessionConfig& config);

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

    [[nodiscard]] base::Expected<void> Initialise();
    void ScheduleNextTile();
    void RenderTile(Tile* tile);
    // Vulkan render path: dispatches the trace kernel per governed tile and
    // fills the display buffer, then fires AllTilesComplete so WriteOutput
    // applies the host display pipeline. Declines loudly (Error) when the
    // backend is absent, the metric is not on the Vulkan render path, or the
    // device cannot be opened.
    void RenderVulkanPath();
    [[nodiscard]] base::Expected<void> WriteOutput();
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
