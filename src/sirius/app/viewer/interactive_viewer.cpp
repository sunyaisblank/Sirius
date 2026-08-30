#include "sirius/app/viewer/interactive_viewer.h"

#include "sirius/core/constants.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <thread>
#include <utility>

namespace sirius::app {

namespace {

double NormaliseAzimuth(double phi) {
    return std::remainder(phi, 2.0 * core::constants::math::kPi);
}

render::SessionConfig ProjectSessionConfig(const ViewerConfig& viewer, const CameraState& camera,
                                           int width, int height, int samples_per_pixel) {
    render::SessionConfig config = viewer.session_template;
    config.width = width;
    config.height = height;
    config.samples_per_pixel = samples_per_pixel;
    config.write_output = false;
    config.output_path = render::SessionConfig{}.output_path;

    config.metric_id = viewer.metric_id;
    config.black_hole_mass = viewer.black_hole_mass;
    config.black_hole_spin = viewer.black_hole_spin;
    config.observer_distance = camera.r;
    config.observer_inclination = camera.theta;
    config.observer_azimuth = camera.phi;
    config.camera_fov = camera.fov;
    config.enable_disk = viewer.enable_disk;
    config.enable_volumetric_disk = viewer.enable_volumetric;
    config.backend = viewer.backend;
    if (config.backend == render::RenderBackend::Vulkan) {
        const render::SessionConfig defaults;
        config.tile_size = defaults.tile_size;
        config.thread_count = defaults.thread_count;
        config.enable_parallel_rendering = defaults.enable_parallel_rendering;
    } else if (!config.enable_parallel_rendering) {
        config.thread_count = 0;
    }

    return config;
}

}  // namespace

InteractiveViewer::InteractiveViewer() = default;

InteractiveViewer::~InteractiveViewer() { Stop(); }

bool InteractiveViewer::Initialise(const ViewerConfig& config) {
    const bool metric_supported = config.metric_id == core::MetricId::Schwarzschild ||
                                  config.metric_id == core::MetricId::Kerr;
    const bool spin_valid =
        std::isfinite(config.black_hole_spin) && config.black_hole_spin >= 0.0 &&
        config.black_hole_spin <= 0.998 &&
        (config.metric_id != core::MetricId::Schwarzschild || config.black_hole_spin == 0.0);
    const bool mass_valid = std::isfinite(config.black_hole_mass) &&
                            config.black_hole_mass >= 0.1 && config.black_hole_mass <= 100.0;
    const double minimum_distance = 5.0 * config.black_hole_mass;
    const double maximum_distance = 1000.0 * config.black_hole_mass;
    if (running_ || !metric_supported || !mass_valid || !spin_valid || config.preview_width < 64 ||
        config.preview_height < 64 || config.final_width < config.preview_width ||
        config.final_height < config.preview_height || config.final_width > 8192 ||
        config.final_height > 8192 || config.refinement_levels < 1 ||
        config.refinement_levels > 32 || config.samples_per_level < 1 ||
        config.samples_per_level > 4096 || !std::isfinite(config.move_speed) ||
        config.move_speed <= 0.0f || !std::isfinite(config.mouse_sensitivity) ||
        config.mouse_sensitivity <= 0.0f || !std::isfinite(config.observer_distance) ||
        config.observer_distance < minimum_distance ||
        config.observer_distance > maximum_distance ||
        !std::isfinite(config.observer_inclination) || config.observer_inclination < 0.1 ||
        config.observer_inclination > core::constants::math::kPi - 0.1 ||
        !std::isfinite(config.observer_azimuth) || !std::isfinite(config.observer_fov) ||
        config.observer_fov < 1.0f || config.observer_fov > 170.0f) {
        return false;
    }
    const CameraState initial_camera{
        .r = config.observer_distance,
        .theta = config.observer_inclination,
        .phi = config.observer_azimuth,
        .fov = config.observer_fov,
    };
    if (render::SessionConfigIssue(ProjectSessionConfig(config, initial_camera,
                                                        config.preview_width, config.preview_height,
                                                        1))
            .has_value()) {
        return false;
    }

    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        config_ = config;
        camera_.r = config.observer_distance;
        camera_.theta = config.observer_inclination;
        camera_.phi = config.observer_azimuth;
        camera_.fov = config.observer_fov;

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

void InteractiveViewer::SetFrameCallback(FrameCallback callback) {
    std::lock_guard<std::mutex> lock(frame_mutex_);
    frame_callback_ = std::move(callback);
}

bool InteractiveViewer::Start() {
    std::lock_guard<std::mutex> lock(thread_mutex_);
    if (running_) return false;
    if (render_thread_.joinable()) {
        render_thread_.join();
    }
    running_ = true;
    stop_requested_ = false;
    restart_requested_ = false;
    {
        std::lock_guard<std::mutex> state_lock(state_mutex_);
        refinement_.needs_restart = false;
    }

    try {
        render_thread_ = std::thread(&InteractiveViewer::RenderThread, this);
    } catch (...) {
        running_ = false;
        return false;
    }
    return true;
}

void InteractiveViewer::Stop() {
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

void InteractiveViewer::Restart() {
    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        refinement_.needs_restart = true;
    }
    restart_requested_ = true;
    std::lock_guard<std::mutex> lock(session_mutex_);
    if (session_) {
        (void)session_->Cancel();
    }
}

void InteractiveViewer::SetCameraPosition(double r, double theta, double phi) {
    if (!std::isfinite(r) || !std::isfinite(theta) || !std::isfinite(phi)) return;
    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        camera_.r = std::clamp(r, 5.0 * config_.black_hole_mass, 1000.0 * config_.black_hole_mass);
        camera_.theta = std::clamp(theta, 0.1, core::constants::math::kPi - 0.1);
        camera_.phi = NormaliseAzimuth(phi);
    }
    Restart();
}

void InteractiveViewer::UpdateCamera(float dt) {
    if (!std::isfinite(dt) || dt <= 0.0f) return;
    bool moved = false;
    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        const double movement = static_cast<double>(config_.move_speed) * static_cast<double>(dt);
        const double angular_radius = camera_.r;
        if (camera_.move_forward) {
            camera_.r -= movement;
            moved = true;
        }
        if (camera_.move_backward) {
            camera_.r += movement;
            moved = true;
        }
        if (camera_.move_left) {
            camera_.phi -= movement / angular_radius;
            moved = true;
        }
        if (camera_.move_right) {
            camera_.phi += movement / angular_radius;
            moved = true;
        }
        if (camera_.move_up) {
            camera_.theta = std::max(0.1, camera_.theta - movement * 0.02);
            moved = true;
        }
        if (camera_.move_down) {
            camera_.theta =
                std::min(core::constants::math::kPi - 0.1, camera_.theta + movement * 0.02);
            moved = true;
        }
        camera_.r =
            std::clamp(camera_.r, 5.0 * config_.black_hole_mass, 1000.0 * config_.black_hole_mass);
        camera_.phi = NormaliseAzimuth(camera_.phi);
    }

    if (moved) {
        Restart();
    }
}

void InteractiveViewer::ProcessKey(int key, int action) {
    // GLFW_RELEASE=0, GLFW_PRESS=1, GLFW_REPEAT=2. A repeat event preserves a
    // held key; treating it as release makes movement stutter under real input.
    if (action < 0 || action > 2) return;
    const bool pressed = action != 0;

    std::lock_guard<std::mutex> lock(state_mutex_);
    // GLFW key codes: W=87, A=65, S=83, D=68, Q=81, E=69.
    switch (key) {
        case 87:
            camera_.move_forward = pressed;
            break;
        case 83:
            camera_.move_backward = pressed;
            break;
        case 65:
            camera_.move_left = pressed;
            break;
        case 68:
            camera_.move_right = pressed;
            break;
        case 81:
            camera_.move_up = pressed;
            break;
        case 69:
            camera_.move_down = pressed;
            break;
        default:
            break;
    }
}

void InteractiveViewer::ProcessMouseMove(double xpos, double ypos, bool dragging) {
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
            const double next_phi = camera_.phi + dx * config_.mouse_sensitivity;
            const double next_theta = camera_.theta + dy * config_.mouse_sensitivity;
            if (std::isfinite(next_phi) && std::isfinite(next_theta)) {
                camera_.phi = NormaliseAzimuth(next_phi);
                camera_.theta = std::clamp(next_theta, 0.1, core::constants::math::kPi - 0.1);
                moved = true;
            }
        }

        last_mouse_x_ = xpos;
        last_mouse_y_ = ypos;
    }
    if (moved) Restart();
}

void InteractiveViewer::ProcessScroll(double yoffset) {
    if (!std::isfinite(yoffset)) return;
    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        camera_.fov = std::clamp(camera_.fov - static_cast<float>(yoffset) * 2.0f, 10.0f, 120.0f);
    }
    Restart();
}

CameraState InteractiveViewer::GetCameraState() const {
    std::lock_guard<std::mutex> lock(state_mutex_);
    return camera_;
}

RefinementState InteractiveViewer::GetRefinementState() const {
    std::lock_guard<std::mutex> lock(state_mutex_);
    return refinement_;
}

std::vector<float> InteractiveViewer::GetFrameBufferSnapshot() const {
    std::lock_guard<std::mutex> lock(frame_mutex_);
    return frame_buffer_;
}

int InteractiveViewer::GetFrameWidth() const {
    std::lock_guard<std::mutex> lock(state_mutex_);
    return refinement_.current_width;
}

int InteractiveViewer::GetFrameHeight() const {
    std::lock_guard<std::mutex> lock(state_mutex_);
    return refinement_.current_height;
}

std::string InteractiveViewer::GetLastError() const {
    std::lock_guard<std::mutex> lock(frame_mutex_);
    return last_error_;
}

void InteractiveViewer::GetResolutionForLevel(int level, int& width, int& height) const {
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

render::SessionConfig InteractiveViewer::CreateSessionConfig(int width, int height, int spp) const {
    ViewerConfig viewer;
    CameraState camera;
    {
        std::lock_guard<std::mutex> lock(state_mutex_);
        viewer = config_;
        camera = camera_;
    }
    return ProjectSessionConfig(viewer, camera, width, height, spp);
}

void InteractiveViewer::RenderThread() {
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
            refinement_.needs_restart = false;
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

bool InteractiveViewer::RenderStep() {
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
    if (auto configured = session->Configure(config); !configured) {
        std::lock_guard<std::mutex> lock(frame_mutex_);
        last_error_ = configured.error().Description();
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
    const std::size_t size = static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 4;
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

void InteractiveViewer::NotifyFrame() {
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
