#pragma once

// Real-time black hole viewer with progressive refinement and orbital camera
// controls. The public surface owns only viewer state and lifecycle; the render
// loop is compiled in interactive_viewer.cpp.

#include "sirius/render/session/render_session.h"

#include <atomic>
#include <chrono>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace sirius::app {

// Viewer configuration: resolution ladder, refinement schedule, camera speeds,
// and the spacetime defaults copied into each SessionConfig.
struct ViewerConfig {
    render::SessionConfig session_template;

    int preview_width = 320;
    int preview_height = 180;
    int final_width = 1920;
    int final_height = 1080;

    int refinement_levels = 4;
    int samples_per_level = 4;

    float move_speed = 2.0f;           // Movement speed [M/s].
    float mouse_sensitivity = 0.002f;  // Mouse look sensitivity.

    double black_hole_mass = 1.0;
    double black_hole_spin = 0.9;
    core::MetricId metric_id = core::MetricId::Kerr;
    double observer_distance = 50.0;
    double observer_inclination = 1.308;  // 75 degrees.
    double observer_azimuth = 0.0;
    float observer_fov = 60.0f;

    bool enable_disk = true;
    bool enable_volumetric = false;
    bool enable_jets = false;
    render::RenderBackend backend = render::RenderBackend::Cpu;
};

// Orbital camera state in spherical coordinates plus per-axis input velocities.
struct CameraState {
    double r = 50.0;       // Radial distance [M].
    double theta = 1.308;  // Inclination [radians].
    double phi = 0.0;      // Azimuthal angle [radians].
    float fov = 60.0f;     // Field of view [degrees].

    bool move_forward = false;
    bool move_backward = false;
    bool move_left = false;
    bool move_right = false;
    bool move_up = false;
    bool move_down = false;
};

// Progressive-refinement bookkeeping.
struct RefinementState {
    int current_level = 0;
    int current_sample = 0;
    int current_samples_per_pixel = 0;
    bool complete = false;
    bool needs_restart = false;

    int current_width = 0;
    int current_height = 0;

    std::chrono::steady_clock::time_point start_time;
    double elapsed_ms = 0.0;
};

// Owns a render thread that refines the image progressively and pushes each
// completed frame to a callback.
class InteractiveViewer {
  public:
    using FrameCallback = std::function<void(const float* data, int width, int height)>;

    InteractiveViewer();
    ~InteractiveViewer();

    InteractiveViewer(const InteractiveViewer&) = delete;
    InteractiveViewer& operator=(const InteractiveViewer&) = delete;
    InteractiveViewer(InteractiveViewer&&) = delete;
    InteractiveViewer& operator=(InteractiveViewer&&) = delete;

    bool Initialise(const ViewerConfig& config);
    void SetFrameCallback(FrameCallback callback);

    bool Start();
    void Stop();
    [[nodiscard]] bool IsRunning() const noexcept { return running_; }
    void Restart();

    void UpdateCamera(float dt);
    void ProcessKey(int key, int action);
    void ProcessMouseMove(double xpos, double ypos, bool dragging);
    void ProcessScroll(double yoffset);

    [[nodiscard]] CameraState GetCameraState() const;
    void SetCameraPosition(double r, double theta, double phi);

    [[nodiscard]] RefinementState GetRefinementState() const;

    [[nodiscard]] std::vector<float> GetFrameBufferSnapshot() const;
    [[nodiscard]] int GetFrameWidth() const;
    [[nodiscard]] int GetFrameHeight() const;
    [[nodiscard]] std::string GetLastError() const;

    // Deterministic projection from viewer state to the live render boundary.
    // Public so the non-windowed test estate can prove that operator controls
    // are not discarded before a preview begins.
    [[nodiscard]] render::SessionConfig CreateSessionConfig(int width, int height, int spp) const;

  private:
    void RenderThread();
    [[nodiscard]] bool RenderStep();
    void GetResolutionForLevel(int level, int& width, int& height) const;
    void NotifyFrame();

    ViewerConfig config_;
    CameraState camera_;
    RefinementState refinement_;

    std::unique_ptr<render::RenderSession> session_;
    std::thread render_thread_;
    std::atomic<bool> running_{false};
    std::atomic<bool> stop_requested_{false};
    std::atomic<bool> restart_requested_{false};

    mutable std::mutex state_mutex_;
    mutable std::mutex thread_mutex_;
    mutable std::mutex session_mutex_;
    mutable std::mutex frame_mutex_;
    std::vector<float> frame_buffer_;
    std::string last_error_;

    FrameCallback frame_callback_;

    double last_mouse_x_ = 0.0;
    double last_mouse_y_ = 0.0;
    bool first_mouse_ = true;
};

}  // namespace sirius::app
