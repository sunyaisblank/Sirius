// Sirius.UI/UIMN001A.cpp - UI Manager Implementation
#include "UIMN001A.h"
#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"
#include "CRPM001A.h"
#include <glad/glad.h>
#include <iostream>
#include <iomanip>
#include <sstream>

#ifdef SIRIUS_HAS_HDF5
#include <PHMT009A.h>
#endif

void UIManager::init(GLFWwindow* window, PluginManager* pluginManager) {
    m_PluginManager = pluginManager;
    
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    
    // Configure ImGui style
    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 5.0f;
    style.FrameRounding = 3.0f;
    style.ScrollbarRounding = 3.0f;
    style.Alpha = 0.95f;
    
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    // Set initial selected metric
    if (!m_PluginManager->getMetricNames().empty()) {
        m_SelectedMetricIndex = 0;
        m_SelectedMetric = m_PluginManager->getMetric(m_PluginManager->getMetricNames()[0]);
    }
    
    std::cout << "UIManager initialized successfully" << std::endl;
}

void UIManager::beginFrame() {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

void UIManager::render() {
    // Main control panel
    ImGui::Begin("Sirius", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
    
    // Menu bar
    if (ImGui::BeginMenuBar()) {
        if (ImGui::BeginMenu("View")) {
            ImGui::MenuItem("Observer Info", nullptr, &m_ShowObserverInfo);
            ImGui::MenuItem("Metric Parameters", nullptr, &m_ShowMetricParams);
            ImGui::MenuItem("Camera Controls", nullptr, &m_ShowCameraControls);
            ImGui::MenuItem("System Info", nullptr, &m_ShowSystemInfo);
            ImGui::EndMenu();
        }
        ImGui::EndMenuBar();
    }
    
    renderMetricSelector();
    
    if (m_ShowMetricParams) {
        ImGui::Separator();
        renderMetricParameters();
    }
    
    if (m_ShowObserverInfo) {
        ImGui::Separator();
        renderObserverInfo();
    }
    
    if (m_ShowCameraControls) {
        ImGui::Separator();
        renderCameraControls();
    }
    
    if (m_ShowSystemInfo) {
        ImGui::Separator();
        renderSystemInfo();
    }
    
    ImGui::End();
}

void UIManager::renderMetricSelector() {
    ImGui::Text("Spacetime Metric");
    
    const auto& metricNames = m_PluginManager->getMetricNames();
    if (ImGui::BeginCombo("##Metric", m_SelectedMetric ? m_SelectedMetric->getName() : "None")) {
        for (int i = 0; i < static_cast<int>(metricNames.size()); ++i) {
            const auto& name = metricNames[i];
            bool is_selected = (m_SelectedMetricIndex == i);
            if (ImGui::Selectable(name.c_str(), is_selected)) {
                if (m_SelectedMetricIndex != i) {
                    m_SelectedMetricIndex = i;
                    m_SelectedMetric = m_PluginManager->getMetric(name);
                    std::cout << "Selected metric: " << name << std::endl;
                    
                    if (m_MetricChangeCallback) {
                        m_MetricChangeCallback(name);
                    }
                }
            }
            if (is_selected) {
                ImGui::SetItemDefaultFocus();
            }
        }
        ImGui::EndCombo();
    }
}

void UIManager::renderMetricParameters() {
    if (!m_SelectedMetric) {
        ImGui::TextColored(ImVec4(0.7f, 0.7f, 0.7f, 1.0f), "No metric selected");
        return;
    }
    
    ImGui::Text("Metric Parameters");
    
    Config& params = const_cast<Config&>(m_SelectedMetric->getParameters());
    if (params.empty()) {
        ImGui::TextColored(ImVec4(0.7f, 0.7f, 0.7f, 1.0f), "No parameters available");
        return;
    }
    
    renderNumericalMetricUI();
    renderParameterSliders(params);
}

// Helper: Render numerical metric HDF5 loading UI
void UIManager::renderNumericalMetricUI() {
#ifdef SIRIUS_HAS_HDF5
    if (std::string(m_SelectedMetric->getName()) == "Numerical") {
        ImGui::Separator();
        ImGui::Text("Numerical Data Source");
        
        static char dataPath[256] = "../../../Sirius.Numerics/Tools/test_data";
        ImGui::InputText("Data Directory", dataPath, 256);
        
        NumericalMetric* numMetric = static_cast<NumericalMetric*>(m_SelectedMetric);
        
        if (ImGui::Button("Load Einstein Toolkit Data")) {
            numMetric->loadData(dataPath);
        }
        
        if (numMetric->isDataLoaded()) {
            ImGui::TextColored(ImVec4(0.2f, 1.0f, 0.2f, 1.0f), "Data Loaded");
        } else {
            ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.0f, 1.0f), "No Data Loaded");
        }
        ImGui::Separator();
    }
#endif
}

// Helper: Render parameter sliders for metric configuration
void UIManager::renderParameterSliders(Config& params) {
    for (auto& pair : params) {
        const std::string& name = pair.first;
        MetricParameter& param = pair.second;
        
        std::string label = name;
        if (name == "mass") label = "Mass (M)";
        else if (name == "charge") label = "Charge (Q)";
        else if (name == "spin") label = "Angular momentum (a)";
        
        bool changed = ImGui::SliderScalar(
            label.c_str(), 
            ImGuiDataType_Double, 
            &param.value, 
            &param.min, 
            &param.max, 
            "%.3f"
        );
        
        if (changed) {
            m_SelectedMetric->setParameter(name, param.value);
        }
        
        if (ImGui::IsItemHovered()) {
            std::stringstream tooltip;
            tooltip << std::fixed << std::setprecision(6);
            tooltip << "Value: " << param.value << "\n";
            tooltip << "Range: [" << param.min << ", " << param.max << "]";
            ImGui::SetTooltip("%s", tooltip.str().c_str());
        }
    }
}

void UIManager::renderObserverInfo() {
    ImGui::Text("Observer State");
    
    // Position display
    ImGui::Text("Position (t, r, θ, φ):");
    ImGui::Indent();
    ImGui::Text("t: %.3f", m_ObserverPosition(0));
    ImGui::Text("r: %.3f", m_ObserverPosition(1));
    ImGui::Text("θ: %.3f (%.1f°)", m_ObserverPosition(2), m_ObserverPosition(2) * 180.0 / 3.14159);
    ImGui::Text("φ: %.3f (%.1f°)", m_ObserverPosition(3), m_ObserverPosition(3) * 180.0 / 3.14159);
    ImGui::Unindent();
    
    // Velocity display
    ImGui::Text("4-Velocity (dt/dτ, dr/dτ, dθ/dτ, dφ/dτ):");
    ImGui::Indent();
    ImGui::Text("u⁰: %.3f", m_ObserverVelocity(0));
    ImGui::Text("u¹: %.3f", m_ObserverVelocity(1));
    ImGui::Text("u²: %.3f", m_ObserverVelocity(2));
    ImGui::Text("u³: %.3f", m_ObserverVelocity(3));
    ImGui::Unindent();
    
    // Derived quantities
    if (m_SelectedMetric) {
        // Calculate coordinate speed
        float spatial_speed = sqrt(
            m_ObserverVelocity(1) * m_ObserverVelocity(1) +
            m_ObserverVelocity(2) * m_ObserverVelocity(2) +
            m_ObserverVelocity(3) * m_ObserverVelocity(3)
        );
        ImGui::Text("Coordinate speed: %.3f", spatial_speed);
        
        // For Schwarzschild metric, show distance from event horizon
        auto params = m_SelectedMetric->getParameters();
        if (params.find("mass") != params.end()) {
            double mass = params.at("mass").value;
            double schwarzschild_radius = 2.0 * mass;
            double distance_from_horizon = m_ObserverPosition(1) - schwarzschild_radius;
            
            if (distance_from_horizon > 0) {
                ImGui::Text("Distance from horizon: %.3f", distance_from_horizon);
            } else {
                ImGui::TextColored(ImVec4(1.0f, 0.3f, 0.3f, 1.0f), "Inside event horizon!");
            }
        }
    }
}

void UIManager::renderCameraControls() {
    ImGui::Text("Camera Settings");
    
    float fov_degrees = m_CameraFOV * 180.0f / 3.14159f;
    ImGui::Text("Field of View: %.1f°", fov_degrees);
    
    ImGui::Separator();
    
    // Geodesic Camera Mode Selector (Phase 8.7)
    ImGui::Text("Camera Motion Mode");
    const char* cameraModes[] = {
        "Free WASD (Manual)",
        "Geodesic Freefall",
        "Circular Orbit",
        "ISCO Orbit",
        "Plunging (Infall)"
    };
    
    int prevMode = m_CameraMode;
    if (ImGui::Combo("##CameraMode", &m_CameraMode, cameraModes, 5)) {
        if (m_CameraMode != prevMode && m_CameraModeCallback) {
            m_CameraModeCallback(m_CameraMode);
        }
    }
    
    // Show mode description
    if (m_CameraMode == 0) {
        ImGui::TextColored(ImVec4(0.7f, 0.7f, 0.7f, 1.0f), "Use WASD/QE/Mouse to control camera");
    } else if (m_CameraMode == 1) {
        ImGui::TextColored(ImVec4(0.5f, 0.8f, 1.0f, 1.0f), "Camera follows geodesic from current position");
    } else if (m_CameraMode == 2) {
        ImGui::TextColored(ImVec4(0.3f, 1.0f, 0.3f, 1.0f), "Stable circular orbit at current radius");
    } else if (m_CameraMode == 3) {
        ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.3f, 1.0f), "Innermost Stable Circular Orbit");
    } else if (m_CameraMode == 4) {
        ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.3f, 1.0f), "Falling toward the black hole");
    }
    
    // Show proper time and time dilation for geodesic modes
    if (m_CameraMode != 0) {
        ImGui::Separator();
        ImGui::Text("Proper Time (τ): %.3f s", m_ProperTime);
        ImGui::Text("Time Dilation (dτ/dt): %.4f", m_TimeDilation);
        
        // Visual indicator for time dilation
        float dilationBar = m_TimeDilation;
        ImGui::ProgressBar(dilationBar, ImVec2(-1, 0), "");
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Time dilation relative to distant observer\n1.0 = same time flow\n0.0 = time frozen (at horizon)");
        }
    }
    
    ImGui::Separator();
    ImGui::Text("Controls:");
    ImGui::BulletText("TAB - Toggle mouse look");
    ImGui::BulletText("WASD - Move in r/φ directions");
    ImGui::BulletText("QE - Move in θ direction");
    ImGui::BulletText("Shift - Fast movement");
    ImGui::BulletText("Ctrl - Slow movement");
    ImGui::BulletText("Mouse - Look around (when mouse look is on)");
}

void UIManager::renderSystemInfo() {
    ImGui::Text("System Information");
    
    ImGuiIO& io = ImGui::GetIO();
    ImGui::Text("Frame rate: %.1f FPS", io.Framerate);
    ImGui::Text("Frame time: %.3f ms", 1000.0f / io.Framerate);
    
    // GPU information would be nice here
    const char* renderer = reinterpret_cast<const char*>(glGetString(GL_RENDERER));
    const char* vendor = reinterpret_cast<const char*>(glGetString(GL_VENDOR));
    const char* version = reinterpret_cast<const char*>(glGetString(GL_VERSION));
    
    if (renderer) ImGui::Text("GPU: %s", renderer);
    if (vendor) ImGui::Text("Vendor: %s", vendor);
    if (version) ImGui::Text("OpenGL: %s", version);
}

void UIManager::setCurrentMetric(IMetric* metric) {
    m_SelectedMetric = metric;
    
    // Update the combo box index
    if (metric && m_PluginManager) {
        const auto& names = m_PluginManager->getMetricNames();
        for (int i = 0; i < static_cast<int>(names.size()); ++i) {
            IMetric* testMetric = m_PluginManager->getMetric(names[i]);
            if (testMetric && std::string(testMetric->getName()) == std::string(metric->getName())) {
                m_SelectedMetricIndex = i;
                break;
            }
        }
    }
}

void UIManager::endFrame() {
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

UIManager::~UIManager() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}