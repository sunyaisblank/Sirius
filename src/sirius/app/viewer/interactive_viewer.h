#pragma once

// Real-time black hole viewer with progressive refinement and orbital camera
// controls. Drives the CPU RenderSession at increasing resolution, low to high,
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
#include <cstring>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

struct GLFWwindow;

namespace sirius::app {

// Viewer configuration: resolution ladder, refinement schedule, camera speeds,
// and the spacetime defaults copied into each SessionConfig.
struct ViewerConfig {
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
    double observerDistance = 50.0;
    double observerInclination = 1.308;  // 75 degrees.

    bool enableDisk = true;
    bool enableVolumetric = false;
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
    void SetFrameCallback(FrameCallback callback) { frame_callback_ = callback; }

    void Start();
    void Stop();
    bool IsRunning() const { return running_; }
    void Restart();

    void UpdateCamera(float dt);
    void ProcessKey(int key, int action);
    void ProcessMouseMove(double xpos, double ypos, bool dragging);
    void ProcessScroll(double yoffset);

    const CameraState& GetCameraState() const { return camera_; }
    void SetCameraPosition(double r, double theta, double phi);

    ViewerConfig& GetConfig() { return config_; }
    const RefinementState& GetRefinementState() const { return refinement_; }

    const float* GetFrameBuffer() const;
    int GetFrameWidth() const { return refinement_.current_width; }
    int GetFrameHeight() const { return refinement_.current_height; }

  private:
    void RenderThread();
    void RenderStep();
    void GetResolutionForLevel(int level, int& width, int& height) const;
    render::SessionConfig CreateSessionConfig(int width, int height, int spp) const;
    void NotifyFrame();

    ViewerConfig config_;
    CameraState camera_;
    RefinementState refinement_;

    std::unique_ptr<render::RenderSession> session_;
    std::thread render_thread_;
    std::atomic<bool> running_{false};
    std::atomic<bool> stop_requested_{false};
    std::atomic<bool> restart_requested_{false};

    std::mutex frame_mutex_;
    std::vector<float> frame_buffer_;

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
    config_ = config;

    camera_.r = config.observerDistance;
    camera_.theta = config.observerInclination;
    camera_.phi = 0.0;
    camera_.fov = 60.0f;

    refinement_.current_level = 0;
    refinement_.current_sample = 0;
    refinement_.complete = false;
    refinement_.needs_restart = false;

    return true;
}

inline void InteractiveViewer::AttachWindow(GLFWwindow* window) { window_ = window; }

inline void InteractiveViewer::Start() {
    if (running_) return;

    running_ = true;
    stop_requested_ = false;
    restart_requested_ = false;

    render_thread_ = std::thread(&InteractiveViewer::RenderThread, this);
}

inline void InteractiveViewer::Stop() {
    stop_requested_ = true;
    if (render_thread_.joinable()) {
        render_thread_.join();
    }
    running_ = false;
}

inline void InteractiveViewer::Restart() { restart_requested_ = true; }

inline void InteractiveViewer::SetCameraPosition(double r, double theta, double phi) {
    camera_.r = r;
    camera_.theta = theta;
    camera_.phi = phi;
    Restart();
}

inline void InteractiveViewer::UpdateCamera(float dt) {
    bool moved = false;

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

    camera_.r = std::max(5.0, camera_.r);

    if (moved) {
        Restart();
    }
}

inline void InteractiveViewer::ProcessKey(int key, int action) {
    bool pressed = (action == 1);  // GLFW_PRESS.

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
    if (first_mouse_) {
        last_mouse_x_ = xpos;
        last_mouse_y_ = ypos;
        first_mouse_ = false;
        return;
    }

    if (dragging) {
        double dx = xpos - last_mouse_x_;
        double dy = ypos - last_mouse_y_;

        camera_.phi += dx * config_.mouse_sensitivity;
        camera_.theta = std::clamp(camera_.theta + dy * config_.mouse_sensitivity, 0.1,
                                   core::constants::math::kPi - 0.1);

        Restart();
    }

    last_mouse_x_ = xpos;
    last_mouse_y_ = ypos;
}

inline void InteractiveViewer::ProcessScroll(double yoffset) {
    camera_.fov = std::clamp(camera_.fov - static_cast<float>(yoffset) * 2.0f, 10.0f, 120.0f);
    Restart();
}

inline const float* InteractiveViewer::GetFrameBuffer() const {
    return frame_buffer_.empty() ? nullptr : frame_buffer_.data();
}

inline void InteractiveViewer::GetResolutionForLevel(int level, int& width, int& height) const {
    float t = static_cast<float>(level) / std::max(1, config_.refinement_levels - 1);
    width =
        static_cast<int>(config_.preview_width + t * (config_.final_width - config_.preview_width));
    height = static_cast<int>(config_.preview_height +
                              t * (config_.final_height - config_.preview_height));

    // Round to a multiple of 8 for tile alignment.
    width = (width / 8) * 8;
    height = (height / 8) * 8;

    width = std::max(64, width);
    height = std::max(64, height);
}

inline render::SessionConfig InteractiveViewer::CreateSessionConfig(int width, int height,
                                                                    int spp) const {
    render::SessionConfig config;
    config.width = width;
    config.height = height;
    config.samplesPerPixel = spp;
    config.tileSize = 32;

    config.blackHoleMass = config_.blackHoleMass;
    config.blackHoleSpin = config_.blackHoleSpin;
    config.observerDistance = camera_.r;
    config.observerInclination = camera_.theta;
    config.cameraFOV = camera_.fov;

    config.enableBloom = true;
    config.exposure = 3.0f;

    config.enableParallelRendering = true;
    // useGPU is inert until the Vulkan seam lands; the CPU path renders here.
    config.useGPU = false;

    return config;
}

inline void InteractiveViewer::RenderThread() {
    refinement_.start_time = std::chrono::steady_clock::now();

    while (!stop_requested_) {
        if (restart_requested_) {
            restart_requested_ = false;
            refinement_.current_level = 0;
            refinement_.current_sample = 0;
            refinement_.complete = false;
            refinement_.start_time = std::chrono::steady_clock::now();
        }

        if (refinement_.complete) {
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
            continue;
        }

        RenderStep();

        refinement_.current_sample++;
        if (refinement_.current_sample >= config_.samples_per_level) {
            refinement_.current_sample = 0;
            refinement_.current_level++;

            if (refinement_.current_level >= config_.refinement_levels) {
                refinement_.complete = true;
            }
        }

        auto now = std::chrono::steady_clock::now();
        refinement_.elapsed_ms =
            std::chrono::duration<double, std::milli>(now - refinement_.start_time).count();
    }
}

inline void InteractiveViewer::RenderStep() {
    int width, height;
    GetResolutionForLevel(refinement_.current_level, width, height);

    refinement_.current_width = width;
    refinement_.current_height = height;

    render::SessionConfig config = CreateSessionConfig(width, height, 1);

    session_ = std::make_unique<render::RenderSession>();
    session_->Configure(config);
    session_->Execute();

    {
        std::lock_guard<std::mutex> lock(frame_mutex_);
        const float* src_data = session_->GetDisplayBuffer().GetFloatData();
        size_t size = static_cast<size_t>(width) * static_cast<size_t>(height) * 4;
        frame_buffer_.assign(src_data, src_data + size);
    }

    NotifyFrame();
}

inline void InteractiveViewer::NotifyFrame() {
    if (frame_callback_) {
        std::lock_guard<std::mutex> lock(frame_mutex_);
        frame_callback_(frame_buffer_.data(), refinement_.current_width,
                        refinement_.current_height);
    }
}

}  // namespace sirius::app
