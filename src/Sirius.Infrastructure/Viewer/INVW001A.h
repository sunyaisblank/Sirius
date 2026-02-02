// INVW001A.h - Interactive Viewer with Progressive Refinement
// Component ID: INVW001A (Infrastructure/Viewer/InteractiveViewer)
//
// Real-time black hole visualization with camera controls and progressive
// image refinement. Starts at low resolution and progressively refines.
//
// FEATURES:
// =========
// 1. Progressive refinement: Low res -> High res over time
// 2. WASD camera movement with mouse look
// 3. Real-time parameter adjustment
// 4. Accumulation buffer for noise reduction
// 5. Immediate restart on camera/parameter change

#pragma once

#include "SNRS001A.h"
#include <PHCN001A.h>
#include <atomic>
#include <thread>
#include <mutex>
#include <chrono>
#include <functional>

struct GLFWwindow;

namespace Sirius {

//==============================================================================
// Viewer Configuration
//==============================================================================
struct ViewerConfig {
    // Resolution levels for progressive refinement
    int preview_width = 320;    ///< Initial preview resolution
    int preview_height = 180;
    int final_width = 1920;     ///< Final high-res
    int final_height = 1080;

    // Progressive refinement
    int refinement_levels = 4;  ///< Number of progressive levels
    int samples_per_level = 4;  ///< Samples per pixel at each level

    // Camera controls
    float move_speed = 2.0f;    ///< Movement speed [M/s]
    float mouse_sensitivity = 0.002f;  ///< Mouse look sensitivity

    // Black hole defaults (copied to SessionConfig)
    double blackHoleMass = 1.0;
    double blackHoleSpin = 0.9;
    double observerDistance = 50.0;
    double observerInclination = 1.308;  // 75 degrees

    // Disk parameters
    bool enableDisk = true;
    bool enableVolumetric = false;
};

//==============================================================================
// Camera State
//==============================================================================
struct CameraState {
    double r = 50.0;        ///< Radial distance [M]
    double theta = 1.308;   ///< Inclination [radians]
    double phi = 0.0;       ///< Azimuthal angle [radians]
    float fov = 60.0f;      ///< Field of view [degrees]

    // For smooth movement
    double dr = 0.0;        ///< Radial velocity
    double dtheta = 0.0;    ///< Inclination velocity
    double dphi = 0.0;      ///< Azimuthal velocity

    // Movement state
    bool moveForward = false;
    bool moveBackward = false;
    bool moveLeft = false;
    bool moveRight = false;
    bool moveUp = false;
    bool moveDown = false;
};

//==============================================================================
// Progressive Refinement State
//==============================================================================
struct RefinementState {
    int current_level = 0;      ///< Current refinement level (0 = lowest)
    int current_sample = 0;     ///< Current sample at this level
    bool complete = false;      ///< True if fully refined
    bool needs_restart = false; ///< True if camera moved (restart)

    int current_width = 0;
    int current_height = 0;

    std::chrono::steady_clock::time_point start_time;
    double elapsed_ms = 0.0;
};

//==============================================================================
// Interactive Viewer
//==============================================================================
class InteractiveViewer {
public:
    using FrameCallback = std::function<void(const float* data, int width, int height)>;
    using ParameterChangeCallback = std::function<void()>;

    InteractiveViewer();
    ~InteractiveViewer();

    //--------------------------------------------------------------------------
    // Initialization
    //--------------------------------------------------------------------------

    /// @brief Initialize viewer with configuration
    bool initialize(const ViewerConfig& config);

    /// @brief Attach to GLFW window for input handling
    void attachWindow(GLFWwindow* window);

    /// @brief Set callback for frame updates (called when new frame available)
    void setFrameCallback(FrameCallback callback) { m_FrameCallback = callback; }

    //--------------------------------------------------------------------------
    // Control
    //--------------------------------------------------------------------------

    /// @brief Start rendering thread
    void start();

    /// @brief Stop rendering thread
    void stop();

    /// @brief Check if rendering is active
    bool isRunning() const { return m_Running; }

    /// @brief Force restart rendering (e.g., after parameter change)
    void restart();

    //--------------------------------------------------------------------------
    // Camera Control
    //--------------------------------------------------------------------------

    /// @brief Update camera from input (call every frame)
    void updateCamera(float dt);

    /// @brief Process keyboard input
    void processKey(int key, int action);

    /// @brief Process mouse movement
    void processMouseMove(double xpos, double ypos, bool dragging);

    /// @brief Process scroll wheel (zoom)
    void processScroll(double yoffset);

    /// @brief Get current camera state
    const CameraState& getCameraState() const { return m_Camera; }

    /// @brief Set camera position directly
    void setCameraPosition(double r, double theta, double phi);

    //--------------------------------------------------------------------------
    // Parameter Access
    //--------------------------------------------------------------------------

    /// @brief Get mutable reference to configuration (restart after changes)
    ViewerConfig& getConfig() { return m_Config; }

    /// @brief Get current refinement state
    const RefinementState& getRefinementState() const { return m_Refinement; }

    /// @brief Get current frame buffer
    const float* getFrameBuffer() const;
    int getFrameWidth() const { return m_Refinement.current_width; }
    int getFrameHeight() const { return m_Refinement.current_height; }

private:
    //--------------------------------------------------------------------------
    // Rendering
    //--------------------------------------------------------------------------

    /// @brief Main render thread function
    void renderThread();

    /// @brief Render one refinement step
    void renderStep();

    /// @brief Compute resolution for refinement level
    void getResolutionForLevel(int level, int& width, int& height) const;

    /// @brief Create session config from current state
    SessionConfig createSessionConfig(int width, int height, int spp) const;

    /// @brief Notify frame update callback
    void notifyFrame();

    //--------------------------------------------------------------------------
    // State
    //--------------------------------------------------------------------------

    ViewerConfig m_Config;
    CameraState m_Camera;
    RefinementState m_Refinement;

    std::unique_ptr<RenderSession> m_Session;
    std::thread m_RenderThread;
    std::atomic<bool> m_Running{false};
    std::atomic<bool> m_StopRequested{false};
    std::atomic<bool> m_RestartRequested{false};

    std::mutex m_FrameMutex;
    std::vector<float> m_FrameBuffer;

    FrameCallback m_FrameCallback;
    GLFWwindow* m_Window = nullptr;

    // Mouse state
    double m_LastMouseX = 0.0;
    double m_LastMouseY = 0.0;
    bool m_FirstMouse = true;
};

//==============================================================================
// Implementation
//==============================================================================

inline InteractiveViewer::InteractiveViewer() = default;

inline InteractiveViewer::~InteractiveViewer() {
    stop();
}

inline bool InteractiveViewer::initialize(const ViewerConfig& config) {
    m_Config = config;

    // Initialize camera from config
    m_Camera.r = config.observerDistance;
    m_Camera.theta = config.observerInclination;
    m_Camera.phi = 0.0;
    m_Camera.fov = 60.0f;

    // Initialize refinement state
    m_Refinement.current_level = 0;
    m_Refinement.current_sample = 0;
    m_Refinement.complete = false;
    m_Refinement.needs_restart = false;

    return true;
}

inline void InteractiveViewer::attachWindow(GLFWwindow* window) {
    m_Window = window;
}

inline void InteractiveViewer::start() {
    if (m_Running) return;

    m_Running = true;
    m_StopRequested = false;
    m_RestartRequested = false;

    m_RenderThread = std::thread(&InteractiveViewer::renderThread, this);
}

inline void InteractiveViewer::stop() {
    m_StopRequested = true;
    if (m_RenderThread.joinable()) {
        m_RenderThread.join();
    }
    m_Running = false;
}

inline void InteractiveViewer::restart() {
    m_RestartRequested = true;
}

inline void InteractiveViewer::setCameraPosition(double r, double theta, double phi) {
    m_Camera.r = r;
    m_Camera.theta = theta;
    m_Camera.phi = phi;
    restart();
}

inline void InteractiveViewer::updateCamera(float dt) {
    bool moved = false;

    // Apply movement based on input state
    if (m_Camera.moveForward) { m_Camera.r -= m_Config.move_speed * dt; moved = true; }
    if (m_Camera.moveBackward) { m_Camera.r += m_Config.move_speed * dt; moved = true; }
    if (m_Camera.moveLeft) { m_Camera.phi -= m_Config.move_speed * dt / m_Camera.r; moved = true; }
    if (m_Camera.moveRight) { m_Camera.phi += m_Config.move_speed * dt / m_Camera.r; moved = true; }
    if (m_Camera.moveUp) { m_Camera.theta = std::max(0.1, m_Camera.theta - m_Config.move_speed * 0.02 * dt); moved = true; }
    if (m_Camera.moveDown) { m_Camera.theta = std::min(Sirius::Constants::Math::PI - 0.1, m_Camera.theta + m_Config.move_speed * 0.02 * dt); moved = true; }

    // Clamp values
    m_Camera.r = std::max(5.0, m_Camera.r);

    if (moved) {
        restart();
    }
}

inline void InteractiveViewer::processKey(int key, int action) {
    bool pressed = (action == 1); // GLFW_PRESS

    // WASD + QE movement (GLFW key codes: W=87, A=65, S=83, D=68, Q=81, E=69)
    switch (key) {
        case 87: m_Camera.moveForward = pressed; break;   // W
        case 83: m_Camera.moveBackward = pressed; break;  // S
        case 65: m_Camera.moveLeft = pressed; break;      // A
        case 68: m_Camera.moveRight = pressed; break;     // D
        case 81: m_Camera.moveUp = pressed; break;        // Q
        case 69: m_Camera.moveDown = pressed; break;      // E
        default: break;
    }
}

inline void InteractiveViewer::processMouseMove(double xpos, double ypos, bool dragging) {
    if (m_FirstMouse) {
        m_LastMouseX = xpos;
        m_LastMouseY = ypos;
        m_FirstMouse = false;
        return;
    }

    if (dragging) {
        double dx = xpos - m_LastMouseX;
        double dy = ypos - m_LastMouseY;

        m_Camera.phi += dx * m_Config.mouse_sensitivity;
        m_Camera.theta = std::clamp(m_Camera.theta + dy * m_Config.mouse_sensitivity, 0.1, Sirius::Constants::Math::PI - 0.1);

        restart();
    }

    m_LastMouseX = xpos;
    m_LastMouseY = ypos;
}

inline void InteractiveViewer::processScroll(double yoffset) {
    // Zoom by changing FOV
    m_Camera.fov = std::clamp(m_Camera.fov - static_cast<float>(yoffset) * 2.0f, 10.0f, 120.0f);
    restart();
}

inline const float* InteractiveViewer::getFrameBuffer() const {
    return m_FrameBuffer.empty() ? nullptr : m_FrameBuffer.data();
}

inline void InteractiveViewer::getResolutionForLevel(int level, int& width, int& height) const {
    // Level 0 = preview, Level N-1 = final
    float t = static_cast<float>(level) / std::max(1, m_Config.refinement_levels - 1);
    width = static_cast<int>(m_Config.preview_width + t * (m_Config.final_width - m_Config.preview_width));
    height = static_cast<int>(m_Config.preview_height + t * (m_Config.final_height - m_Config.preview_height));

    // Round to multiple of 8 for tile alignment
    width = (width / 8) * 8;
    height = (height / 8) * 8;

    width = std::max(64, width);
    height = std::max(64, height);
}

inline SessionConfig InteractiveViewer::createSessionConfig(int width, int height, int spp) const {
    SessionConfig config;
    config.width = width;
    config.height = height;
    config.samplesPerPixel = spp;
    config.tileSize = 32;

    config.blackHoleMass = m_Config.blackHoleMass;
    config.blackHoleSpin = m_Config.blackHoleSpin;
    config.observerDistance = m_Camera.r;
    config.observerInclination = m_Camera.theta;
    config.cameraFOV = m_Camera.fov;

    config.enableBloom = true;
    config.exposure = 3.0f;

    config.enableParallelRendering = true;
    config.useGPU = true;  // Use GPU if available

    return config;
}

inline void InteractiveViewer::renderThread() {
    m_Refinement.start_time = std::chrono::steady_clock::now();

    while (!m_StopRequested) {
        // Check for restart request
        if (m_RestartRequested) {
            m_RestartRequested = false;
            m_Refinement.current_level = 0;
            m_Refinement.current_sample = 0;
            m_Refinement.complete = false;
            m_Refinement.start_time = std::chrono::steady_clock::now();
        }

        // Skip if already complete
        if (m_Refinement.complete) {
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
            continue;
        }

        // Render current step
        renderStep();

        // Advance to next step
        m_Refinement.current_sample++;
        if (m_Refinement.current_sample >= m_Config.samples_per_level) {
            m_Refinement.current_sample = 0;
            m_Refinement.current_level++;

            if (m_Refinement.current_level >= m_Config.refinement_levels) {
                m_Refinement.complete = true;
            }
        }

        // Update elapsed time
        auto now = std::chrono::steady_clock::now();
        m_Refinement.elapsed_ms = std::chrono::duration<double, std::milli>(now - m_Refinement.start_time).count();
    }
}

inline void InteractiveViewer::renderStep() {
    // Get resolution for current level
    int width, height;
    getResolutionForLevel(m_Refinement.current_level, width, height);

    m_Refinement.current_width = width;
    m_Refinement.current_height = height;

    // Create session config
    SessionConfig config = createSessionConfig(width, height, 1);

    // Create and run session
    m_Session = std::make_unique<RenderSession>();
    m_Session->configure(config);
    m_Session->execute();

    // Copy result to frame buffer
    {
        std::lock_guard<std::mutex> lock(m_FrameMutex);
        const float* srcData = m_Session->getDisplayBuffer().getFloatData();
        size_t size = width * height * 4;
        m_FrameBuffer.assign(srcData, srcData + size);
    }

    // Notify callback
    notifyFrame();
}

inline void InteractiveViewer::notifyFrame() {
    if (m_FrameCallback) {
        std::lock_guard<std::mutex> lock(m_FrameMutex);
        m_FrameCallback(m_FrameBuffer.data(), m_Refinement.current_width, m_Refinement.current_height);
    }
}

} // namespace Sirius
