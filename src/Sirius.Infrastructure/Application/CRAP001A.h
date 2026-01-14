// CRAP001A.h - Application Class (Static Rendering)
//
// Manages renderer, metrics, and optional preview window.
// Real-time mode deprecated - use RenderJob for batch rendering.

#pragma once
#include <memory>
#include <string>
#include <MTTN001A.h>

class Window;
class UIManager;
class PluginManager;
class IMetric;
class Renderer;

namespace Sirius {
    class RenderJob;
}

/// @brief Application for static black hole rendering
/// @note Real-time mode is deprecated. Use RenderJob directly for rendering.
class Application {
public:
    Application();
    ~Application();

    /// @brief Run the interactive mode (DEPRECATED)
    /// @deprecated Use RenderJob for batch rendering instead
    [[deprecated("Use RenderJob for batch rendering")]]
    void run();

    // Getters for system components
    Renderer* getRenderer() const { return m_Renderer.get(); }
    PluginManager* getPluginManager() const { return m_PluginManager.get(); }
    
    // Metric management
    IMetric* getCurrentMetric() const { return m_CurrentMetric; }
    const std::string& getCurrentMetricName() const { return m_CurrentMetricName; }
    void setCurrentMetric(const std::string& name);
    
    // Observer state management
    void setObserverPosition(const Vec4& position);
    void setObserverVelocity(const Vec4& velocity);
    Vec4 getObserverPosition() const { return m_ObserverPosition; }
    Vec4 getObserverVelocity() const { return m_ObserverVelocity; }
    
    // Camera parameters
    void setCameraFOV(float fov);
    float getCameraFOV() const { return m_CameraFOV; }

    // Window event handling (for preview mode)
    void resize(int width, int height);
    static void framebufferSizeCallback(struct GLFWwindow* window, int width, int height);

private:
    // Initialization helpers
    void initializeMetrics();
    void initializeRenderer();
    void initializeUI();
    void initializeObserver();
    void initializeCallbacks();
    
    // Legacy input handling (deprecated)
    void handleInput();
    void handleMouse(double xpos, double ypos);
    void updateObserver(float deltaTime);
    void handleKeyToggles(struct GLFWwindow* window);
    void handleMovement(struct GLFWwindow* window, float speed);

    // System components
    std::unique_ptr<Window> m_Window;
    std::unique_ptr<UIManager> m_UIManager;
    std::unique_ptr<PluginManager> m_PluginManager;
    std::unique_ptr<Renderer> m_Renderer;

    // Current metric state
    IMetric* m_CurrentMetric = nullptr;
    std::string m_CurrentMetricName;
    
    // Observer state
    Vec4 m_ObserverPosition;
    Vec4 m_ObserverVelocity;
    float m_CameraFOV = 60.0f;
    float m_CameraYaw = 0.0f;
    float m_CameraPitch = 0.0f;
    
    // Legacy input state (deprecated)
    bool m_MouseLook = false;
    double m_LastMouseX = 0.0;
    double m_LastMouseY = 0.0;
    bool m_FirstMouse = true;
    
    friend class UIManager;
};
