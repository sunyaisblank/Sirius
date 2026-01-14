// UIMN001A.h - UI Manager
//
// ImGui-based interface for metric selection, parameter adjustment,
// observer state display, and camera controls.
#ifndef UIMN001A_H
#define UIMN001A_H

#include "CRPM001A.h"
#include <MTTN001A.h>
#include <functional>
#include <string>

struct GLFWwindow;

class UIManager {
public:
    void init(GLFWwindow* window, PluginManager* pluginManager);
    void beginFrame();
    void render();
    void endFrame();
    ~UIManager();

    // Metric management
    IMetric* getSelectedMetric() const { return m_SelectedMetric; }
    void setCurrentMetric(IMetric* metric);
    
    // Observer state display
    void setObserverPosition(const Vec4& position) { m_ObserverPosition = position; }
    void setObserverVelocity(const Vec4& velocity) { m_ObserverVelocity = velocity; }
    void setCameraFOV(float fov) { m_CameraFOV = fov; }
    
    // Geodesic camera state (Phase 8.7)
    void setCameraMode(int mode) { m_CameraMode = mode; }
    void setProperTime(float tau) { m_ProperTime = tau; }
    void setTimeDilation(float dilation) { m_TimeDilation = dilation; }

    // Callbacks
    void setMetricChangeCallback(std::function<void(const std::string&)> callback) { m_MetricChangeCallback = callback; }
    void setCameraModeCallback(std::function<void(int)> callback) { m_CameraModeCallback = callback; }

private:
    void renderMetricSelector();
    void renderMetricParameters();
    void renderNumericalMetricUI();
    void renderParameterSliders(Config& params);
    void renderObserverInfo();
    void renderCameraControls();
    void renderSystemInfo();

    PluginManager* m_PluginManager = nullptr;
    IMetric* m_SelectedMetric = nullptr;
    int m_SelectedMetricIndex = -1;
    
    // Observer state (for display only)
    Vec4 m_ObserverPosition;
    Vec4 m_ObserverVelocity;
    float m_CameraFOV = 60.0f;
    
    // Geodesic camera state (Phase 8.7)
    int m_CameraMode = 0;      // 0=FreeWASD, 1=Freefall, 2=Circular, 3=ISCO, 4=Plunging
    float m_ProperTime = 0.0f;
    float m_TimeDilation = 1.0f;
    
    // UI state
    bool m_ShowObserverInfo = true;
    bool m_ShowMetricParams = true;
    bool m_ShowCameraControls = true;
    bool m_ShowSystemInfo = false;
    
    std::function<void(const std::string&)> m_MetricChangeCallback;
    std::function<void(int)> m_CameraModeCallback;
};

#endif // UIMN001A_H