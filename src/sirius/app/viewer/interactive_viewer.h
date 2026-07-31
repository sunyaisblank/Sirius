#pragma once

// Real-time black hole viewer with progressive refinement and orbital camera
// controls. Drives RenderSession at increasing resolution, low to high,
// restarting on any camera or parameter change. Ported from UIVW001A.h; methods
// take the new-tree PascalCase spelling and data members the snake_case_ form.
//
// The refinement loop, camera integration, and per-level resolution schedule are
// preserved exactly; only names and the session API (render::RenderSession)
// change from the legacy source.

#include "sirius/core/constants.h"
#include "sirius/render/session/render_session.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstring>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

struct GLFWwindow;

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

    double blackHoleMass = 1.0;
    double blackHoleSpin = 0.9;
    core::MetricId metricId = core::MetricId::Kerr;
    double observerDistance = 50.0;
    double observerInclination = 1.308;  // 75 degrees.
    double observerAzimuth = 0.0;
    float observerFov = 60.0f;

    bool enableDisk = true;
    bool enableVolumetric = false;
    bool enableJets = false;
    render::RenderBackend backend = render::RenderBackend::Cpu;
};

// Orbital camera state in spherical coordinates plus per-axis input velocities.
struct CameraState {
    double r = 50.0;       // Radial distance [M].
    double theta = 1.308;  // Inclination [radians].
    double phi = 0.0;      // Azimuthal angle [radians].
    float fov = 60.0f;     // Field of view [degrees].

    double dr = 0.0;
    double dtheta = 0.0;
    double dphi = 0.0;

    bool moveForward = false;
    bool moveBackward = false;
    bool moveLeft = false;
    bool moveRight = false;
    bool moveUp = false;
    bool moveDown = false;
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

    bool Initialise(const ViewerConfig& config);
    void AttachWindow(GLFWwindow* window);
    void SetFrameCallback(FrameCallback callback);

    bool Start();
    void Stop();
    bool IsRunning() const { return running_; }
    void Restart();

    void UpdateCamera(float dt);
    void ProcessKey(int key, int action);
    void ProcessMouseMove(double xpos, double ypos, bool dragging);
    void ProcessScroll(double yoffset);

    CameraState GetCameraState() const;
    void SetCameraPosition(double r, double theta, double phi);

    RefinementState GetRefinementState() const;

    std::vector<float> GetFrameBufferSnapshot() const;
    int GetFrameWidth() const;
    int GetFrameHeight() const;
    std::string GetLastError() const;

    // Deterministic projection from viewer state to the live render boundary.
    // Public so the non-windowed test estate can prove that operator controls
    // are not discarded before a preview begins.
    render::SessionConfig CreateSessionConfig(int width, int height, int spp) const;

  private:
    void RenderThread();
    bool RenderStep();
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
    GLFWwindow* window_ = nullptr;

    double last_mouse_x_ = 0.0;
    double last_mouse_y_ = 0.0;
    bool first_mouse_ = true;
};

// =============================================================================
// Implementation (header-only, as the legacy viewer was).
// =============================================================================

inline InteractiveViewer::InteractiveViewer() = default;

inline InteractiveViewer::~InteractiveViewer() { Stop(); }

inline bool InteractiveViewer::Initialise(const ViewerConfig& config) {
    const bool metric_supported =
        config.metricId == core::MetricId::Schwarzschild || config.metricId == core::MetricId::Kerr;
    const bool spin_valid =
        std::isfinite(config.blackHoleSpin) && config.blackHoleSpin >= 0.0 &&
        config.blackHoleSpin <= 0.998 &&
        (config.metricId != core::MetricId::Schwarzschild || config.blackHoleSpin == 0.0);
    const bool mass_valid = std::isfinite(config.blackHoleMass) && config.blackHoleMass >= 0.1 &&
                            config.blackHoleMass <= 100.0;
    const double minimum_distance = 5.0 * config.blackHoleMass;
    const double maximum_distance = 1000.0 * config.blackHoleMass;
    if (running_ || !metric_supported || !mass_valid || !spin_valid || config.preview_width < 64 ||
        config.preview_height < 64 || config.final_width < config.preview_width ||
        config.final_height < config.preview_height || config.final_width > 8192 ||
        config.final_height > 8192 || config.refinement_levels < 1 ||
        config.refinement_levels > 32 || config.samples_per_level < 1 ||
        config.samples_per_level > 4096 || !std::isfinite(config.move_speed) ||
        config.move_speed <= 0.0f || !std::isfinite(config.mouse_sensitivity) ||
        config.mouse_sensitivity <= 0.0f || !std::isfinite(config.observerDistance) ||
        config.observerDistance < minimum_distance || config.observerDistance > maximum_distance ||
        !std::isfinite(config.observerInclination) || config.observerInclination < 0.1 ||
        config.observerInclination > core::constants::math::kPi - 0.1 ||
        !std::isfinite(config.observerAzimuth) || !std::isfinite(config.observerFov) ||
        config.observerFov < 1.0f || config.observerFov > 170.0f) {
        return false;
    }

    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        config_ = config;
        camera_.r = config.observerDistance;
        camera_.theta = config.observerInclination;
        camera_.phi = config.observerAzimuth;
        camera_.fov = config.observerFov;

        refinement_.current_level = 0;
        refinement_.current_sample = 0;
        refinement_.current_samples_per_pixel = 0;
        refinement_.complete = false;
        refinement_.needs_restart = false;
        first_mouse_ = true;
        last_mouse_x_ = 0.0;
        last_mouse_y_ = 0.0;
    }
    {
        // GetLastError, frame snapshots, and the render thread all use the
        // frame mutex for this state. Reinitialisation starts a fresh epoch.
        std::lock_guard<std::mutex> lock(frame_mutex_);
        last_error_.clear();
        frame_buffer_.clear();
    }

    return true;
}

inline void InteractiveViewer::AttachWindow(GLFWwindow* window) {
    std::lock_guard<std::mutex> lock(state_mutex_);
    window_ = window;
}

inline void InteractiveViewer::SetFrameCallback(FrameCallback callback) {
    std::lock_guard<std::mutex> lock(frame_mutex_);
    frame_callback_ = std::move(callback);
}

inline bool InteractiveViewer::Start() {
    std::lock_guard<std::mutex> lock(thread_mutex_);
    if (running_) return false;
    if (render_thread_.joinable()) {
        render_thread_.join();
    }
    running_ = true;
    stop_requested_ = false;
    restart_requested_ = false;

    try {
        render_thread_ = std::thread(&InteractiveViewer::RenderThread, this);
    } catch (...) {
        running_ = false;
        return false;
    }
    return true;
}

inline void InteractiveViewer::Stop() {
    stop_requested_ = true;
    {
        std::lock_guard<std::mutex> lock(session_mutex_);
        if (session_) {
            (void)session_->Cancel();
        }
    }
    std::lock_guard<std::mutex> lock(thread_mutex_);
    if (render_thread_.joinable()) {
        render_thread_.join();
    }
    running_ = false;
}

inline void InteractiveViewer::Restart() {
    restart_requested_ = true;
    std::lock_guard<std::mutex> lock(session_mutex_);
    if (session_) {
        (void)session_->Cancel();
    }
}

inline void InteractiveViewer::SetCameraPosition(double r, double theta, double phi) {
    if (!std::isfinite(r) || !std::isfinite(theta) || !std::isfinite(phi)) return;
    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        camera_.r = std::clamp(r, 5.0 * config_.blackHoleMass, 1000.0 * config_.blackHoleMass);
        camera_.theta = std::clamp(theta, 0.1, core::constants::math::kPi - 0.1);
        camera_.phi = phi;
    }
    Restart();
}

inline void InteractiveViewer::UpdateCamera(float dt) {
    if (!std::isfinite(dt) || dt <= 0.0f) return;
    bool moved = false;
    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        if (camera_.moveForward) {
            camera_.r -= config_.move_speed * dt;
            moved = true;
        }
        if (camera_.moveBackward) {
            camera_.r += config_.move_speed * dt;
            moved = true;
        }
        if (camera_.moveLeft) {
            camera_.phi -= config_.move_speed * dt / camera_.r;
            moved = true;
        }
        if (camera_.moveRight) {
            camera_.phi += config_.move_speed * dt / camera_.r;
            moved = true;
        }
        if (camera_.moveUp) {
            camera_.theta = std::max(0.1, camera_.theta - config_.move_speed * 0.02 * dt);
            moved = true;
        }
        if (camera_.moveDown) {
            camera_.theta = std::min(core::constants::math::kPi - 0.1,
                                     camera_.theta + config_.move_speed * 0.02 * dt);
            moved = true;
        }
        camera_.r =
            std::clamp(camera_.r, 5.0 * config_.blackHoleMass, 1000.0 * config_.blackHoleMass);
    }

    if (moved) {
        Restart();
    }
}

inline void InteractiveViewer::ProcessKey(int key, int action) {
    // GLFW_RELEASE=0, GLFW_PRESS=1, GLFW_REPEAT=2. A repeat event preserves a
    // held key; treating it as release makes movement stutter under real input.
    if (action < 0 || action > 2) return;
    const bool pressed = action != 0;

    std::lock_guard<std::mutex> lock(state_mutex_);
    // GLFW key codes: W=87, A=65, S=83, D=68, Q=81, E=69.
    switch (key) {
        case 87:
            camera_.moveForward = pressed;
            break;
        case 83:
            camera_.moveBackward = pressed;
            break;
        case 65:
            camera_.moveLeft = pressed;
            break;
        case 68:
            camera_.moveRight = pressed;
            break;
        case 81:
            camera_.moveUp = pressed;
            break;
        case 69:
            camera_.moveDown = pressed;
            break;
        default:
            break;
    }
}

inline void InteractiveViewer::ProcessMouseMove(double xpos, double ypos, bool dragging) {
    if (!std::isfinite(xpos) || !std::isfinite(ypos)) return;
    bool moved = false;
    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        if (first_mouse_) {
            last_mouse_x_ = xpos;
            last_mouse_y_ = ypos;
            first_mouse_ = false;
            return;
        }

        if (dragging) {
            const double dx = xpos - last_mouse_x_;
            const double dy = ypos - last_mouse_y_;

            camera_.phi += dx * config_.mouse_sensitivity;
            camera_.theta = std::clamp(camera_.theta + dy * config_.mouse_sensitivity, 0.1,
                                       core::constants::math::kPi - 0.1);
            moved = true;
        }

        last_mouse_x_ = xpos;
        last_mouse_y_ = ypos;
    }
    if (moved) Restart();
}

inline void InteractiveViewer::ProcessScroll(double yoffset) {
    if (!std::isfinite(yoffset)) return;
    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        camera_.fov = std::clamp(camera_.fov - static_cast<float>(yoffset) * 2.0f, 10.0f, 120.0f);
    }
    Restart();
}

inline CameraState InteractiveViewer::GetCameraState() const {
    std::lock_guard<std::mutex> lock(state_mutex_);
    return camera_;
}

inline RefinementState InteractiveViewer::GetRefinementState() const {
    std::lock_guard<std::mutex> lock(state_mutex_);
    return refinement_;
}

inline std::vector<float> InteractiveViewer::GetFrameBufferSnapshot() const {
    std::lock_guard<std::mutex> lock(frame_mutex_);
    return frame_buffer_;
}

inline int InteractiveViewer::GetFrameWidth() const {
    std::lock_guard<std::mutex> lock(state_mutex_);
    return refinement_.current_width;
}

inline int InteractiveViewer::GetFrameHeight() const {
    std::lock_guard<std::mutex> lock(state_mutex_);
    return refinement_.current_height;
}

inline std::string InteractiveViewer::GetLastError() const {
    std::lock_guard<std::mutex> lock(frame_mutex_);
    return last_error_;
}

inline void InteractiveViewer::GetResolutionForLevel(int level, int& width, int& height) const {
    ViewerConfig config;
    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        config = config_;
    }
    float t = static_cast<float>(level) / std::max(1, config.refinement_levels - 1);
    width =
        static_cast<int>(config.preview_width + t * (config.final_width - config.preview_width));
    height =
        static_cast<int>(config.preview_height + t * (config.final_height - config.preview_height));

    // Round to a multiple of 8 for tile alignment.
    width = (width / 8) * 8;
    height = (height / 8) * 8;

    width = std::max(64, width);
    height = std::max(64, height);
}

inline render::SessionConfig InteractiveViewer::CreateSessionConfig(int width, int height,
                                                                    int spp) const {
    ViewerConfig viewer;
    CameraState camera;
    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        viewer = config_;
        camera = camera_;
    }
    render::SessionConfig config = viewer.session_template;
    config.width = width;
    config.height = height;
    config.samplesPerPixel = spp;
    config.tileSize = 32;
    config.writeOutput = false;

    config.metricId = viewer.metricId;
    config.blackHoleMass = viewer.blackHoleMass;
    config.blackHoleSpin = viewer.blackHoleSpin;
    config.observerDistance = camera.r;
    config.observerInclination = camera.theta;
    config.observerAzimuth = camera.phi;
    config.cameraFOV = camera.fov;
    config.enableDisk = viewer.enableDisk;
    config.enableVolumetricDisk = viewer.enableDisk && viewer.enableVolumetric;
    config.enableJets = viewer.enableJets || config.enableJets;

    config.enableBloom = true;
    config.exposure = 3.0f;

    config.enableParallelRendering = true;
    config.backend = viewer.backend;

    return config;
}

inline void InteractiveViewer::RenderThread() {
    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        refinement_.start_time = std::chrono::steady_clock::now();
    }

    while (!stop_requested_) {
        if (restart_requested_.exchange(false)) {
            std::lock_guard<std::mutex> lock(state_mutex_);
            refinement_.current_level = 0;
            refinement_.current_sample = 0;
            refinement_.current_samples_per_pixel = 0;
            refinement_.complete = false;
            refinement_.start_time = std::chrono::steady_clock::now();
        }

        bool complete = false;
        {
            std::lock_guard<std::mutex> lock(state_mutex_);
            complete = refinement_.complete;
        }
        if (complete) {
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
            continue;
        }

        if (!RenderStep()) {
            running_ = false;
            break;
        }

        if (restart_requested_) {
            continue;
        }
        {
            std::lock_guard<std::mutex> lock(state_mutex_);
            refinement_.current_sample++;
            if (refinement_.current_sample >= config_.samples_per_level) {
                refinement_.current_sample = 0;
                refinement_.current_level++;

                if (refinement_.current_level >= config_.refinement_levels) {
                    refinement_.complete = true;
                }
            }

            const auto now = std::chrono::steady_clock::now();
            refinement_.elapsed_ms =
                std::chrono::duration<double, std::milli>(now - refinement_.start_time).count();
        }
    }
    running_ = false;
}

inline bool InteractiveViewer::RenderStep() {
    int level = 0;
    int samples_per_pixel = 1;
    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        level = refinement_.current_level;
        samples_per_pixel = refinement_.current_sample + 1;
    }
    int width = 0;
    int height = 0;
    GetResolutionForLevel(level, width, height);

    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        refinement_.current_width = width;
        refinement_.current_height = height;
        refinement_.current_samples_per_pixel = samples_per_pixel;
    }

    render::SessionConfig config = CreateSessionConfig(width, height, samples_per_pixel);

    render::RenderSession* session = nullptr;
    {
        std::lock_guard<std::mutex> lock(session_mutex_);
        session_ = std::make_unique<render::RenderSession>();
        session = session_.get();
    }
    std::string session_error;
    session->SetCompletionCallback([&](render::SessionState state, const std::string& message) {
        if (state != render::SessionState::Complete) session_error = message;
    });
    if (!session->Configure(config)) {
        std::lock_guard<std::mutex> lock(frame_mutex_);
        last_error_ = "viewer render session configuration failed";
        return false;
    }
    const render::SessionState final_state = session->Execute();
    if (final_state != render::SessionState::Complete) {
        if (final_state == render::SessionState::Cancelled &&
            (stop_requested_ || restart_requested_)) {
            return !stop_requested_;
        }
        std::lock_guard<std::mutex> lock(frame_mutex_);
        last_error_ = session_error.empty() ? "viewer render session failed" : session_error;
        return false;
    }

    const std::vector<float> snapshot = session->GetDisplayBuffer().SnapshotFloatData();
    const size_t size = static_cast<size_t>(width) * static_cast<size_t>(height) * 4;
    if (snapshot.size() != size) {
        std::lock_guard<std::mutex> lock(frame_mutex_);
        last_error_ = "viewer render completed without a correctly sized display buffer";
        return false;
    }
    {
        std::lock_guard<std::mutex> lock(frame_mutex_);
        frame_buffer_ = snapshot;
    }

    NotifyFrame();
    return true;
}

inline void InteractiveViewer::NotifyFrame() {
    FrameCallback callback;
    std::vector<float> frame;
    {
        std::lock_guard<std::mutex> lock(frame_mutex_);
        callback = frame_callback_;
        frame = frame_buffer_;
    }
    if (callback) {
        const RefinementState refinement = GetRefinementState();
        callback(frame.data(), refinement.current_width, refinement.current_height);
    }
}

}  // namespace sirius::app
