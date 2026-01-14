// CRAP001A.cpp - Application Implementation (Legacy Interactive Mode)
//
// NOTE: Real-time mode is deprecated. Use RenderJob for batch rendering.
// This file is kept for backwards compatibility with --interactive flag.

#include "CRAP001A.h"
#include "CRWN001A.h"
#include "CRPM001A.h"
#include <RDRT001A.h>
#include <UIMN001A.h>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <chrono>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Application::Application() {
    std::cout << "WARNING: Interactive mode is deprecated. Use batch rendering instead." << std::endl;
    
    m_Window = std::make_unique<Window>(1280, 720, "Sirius (Deprecated Interactive Mode)");
    
    initializeMetrics();
    initializeRenderer();
    initializeUI();
    initializeObserver();
    initializeCallbacks();
}

void Application::initializeMetrics() {
    m_PluginManager = std::make_unique<PluginManager>();
    m_PluginManager->loadPlugins();

    auto metricNames = m_PluginManager->getMetricNames();
    if (!metricNames.empty()) {
        setCurrentMetric(metricNames[0]);
    }
}

void Application::initializeRenderer() {
    m_Renderer = std::make_unique<Renderer>();
    m_Renderer->init(1280, 720);
    
    if (!m_Renderer->loadBackgroundTexture("Sirius.Render/Texture/Starfield.png")) {
        std::cout << "Using procedural starfield (background texture not found)" << std::endl;
    }
}

void Application::initializeUI() {
    m_UIManager = std::make_unique<UIManager>();
    m_UIManager->init(m_Window->getGLFWwindow(), m_PluginManager.get());
    
    m_UIManager->setMetricChangeCallback([this](const std::string& name) {
        this->setCurrentMetric(name);
    });
    
    // Camera mode callback removed - geodesic camera deprecated
}

void Application::initializeObserver() {
    m_ObserverPosition = Vec4(); 
    m_ObserverPosition(0) = 0.0;
    m_ObserverPosition(1) = 10.0;
    m_ObserverPosition(2) = static_cast<double>(M_PI/2.0);
    m_ObserverPosition(3) = 0.0;
    
    m_ObserverVelocity = Vec4(); 
    m_ObserverVelocity(0) = 1.0;
    m_ObserverVelocity(1) = 0.0;
    m_ObserverVelocity(2) = 0.0;
    m_ObserverVelocity(3) = 0.0;
    
    m_CameraFOV = 60.0f * static_cast<float>(M_PI) / 180.0f;
    m_MouseLook = false;
    m_LastMouseX = 640.0;
    m_LastMouseY = 360.0;
    m_FirstMouse = true;
}

void Application::initializeCallbacks() {
    glfwSetWindowUserPointer(m_Window->getGLFWwindow(), this);
    glfwSetFramebufferSizeCallback(m_Window->getGLFWwindow(), framebufferSizeCallback);
}

Application::~Application() = default;

void Application::setCurrentMetric(const std::string& name) {
    m_CurrentMetric = m_PluginManager->createMetric(name);
    m_CurrentMetricName = name;
    
    if (m_CurrentMetric) {
        auto params = m_CurrentMetric->getParameters();
        
        if (params.find("mass") != params.end()) {
            double mass = params.at("mass").value;
            double schwarzschild_radius = 2.0 * mass;
            m_ObserverPosition(1) = std::max(10.0, schwarzschild_radius * 2.0);
        }
        
        if (m_UIManager) {
            m_UIManager->setCurrentMetric(m_CurrentMetric);
        }
        
        std::cout << "Switched to metric: " << name << std::endl;
    } else {
        std::cerr << "Failed to create metric: " << name << std::endl;
    }
}

void Application::setObserverPosition(const Vec4& position) {
    m_ObserverPosition = position;
    
    if (m_CurrentMetric) {
        auto params = m_CurrentMetric->getParameters();
        if (params.find("mass") != params.end()) {
            double mass = params.at("mass").value;
            double schwarzschild_radius = 2.0 * mass;
            m_ObserverPosition(1) = std::max(m_ObserverPosition(1), schwarzschild_radius * 1.1);
        }
    }
    
    m_ObserverPosition(2) = std::clamp(m_ObserverPosition(2), 0.01, M_PI - 0.01);
}

void Application::setObserverVelocity(const Vec4& velocity) {
    m_ObserverVelocity = velocity;
    
    if (m_CurrentMetric) {
        if (m_ObserverVelocity(0) <= 0.0) {
            m_ObserverVelocity(0) = 1.0;
        }
    }
}

void Application::setCameraFOV(float fov) {
    m_CameraFOV = std::clamp(fov, 0.1f, 3.0f);
}

void Application::resize(int width, int height) {
    if (m_Renderer) {
        m_Renderer->resize(width, height);
    }
}

void Application::framebufferSizeCallback(GLFWwindow* window, int width, int height) {
    Application* app = static_cast<Application*>(glfwGetWindowUserPointer(window));
    if (app) {
        app->resize(width, height);
    }
}

void Application::handleInput() {
    GLFWwindow* window = m_Window->getGLFWwindow();
    
    ImGuiIO& io = ImGui::GetIO();
    if (io.WantCaptureKeyboard || io.WantCaptureMouse) {
        return;
    }
    
    handleKeyToggles(window);
    
    if (!m_MouseLook) return;
    
    float speed = 0.1f;
    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS) speed *= 5.0f;
    if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS) speed *= 0.2f;
    
    handleMovement(window, speed);
    
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);
    handleMouse(xpos, ypos);
}

void Application::handleKeyToggles(GLFWwindow* window) {
    static bool tab_was_pressed = false;
    if (glfwGetKey(window, GLFW_KEY_TAB) == GLFW_PRESS) {
        if (!tab_was_pressed) {
            m_MouseLook = !m_MouseLook;
            glfwSetInputMode(window, GLFW_CURSOR, m_MouseLook ? GLFW_CURSOR_DISABLED : GLFW_CURSOR_NORMAL);
            tab_was_pressed = true;
        }
    } else {
        tab_was_pressed = false;
    }
    
    static bool b_was_pressed = false;
    if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS) {
        if (!b_was_pressed && m_Renderer) {
            bool enabled = !m_Renderer->isBloomEnabled();
            m_Renderer->setBloomEnabled(enabled);
            std::cout << "Bloom: " << (enabled ? "ON" : "OFF") << std::endl;
            b_was_pressed = true;
        }
    } else {
        b_was_pressed = false;
    }
    
    static bool l_was_pressed = false;
    if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS) {
        if (!l_was_pressed && m_Renderer) {
            bool enabled = !m_Renderer->isLensFlareEnabled();
            m_Renderer->setLensFlareEnabled(enabled);
            std::cout << "Lens Flare: " << (enabled ? "ON" : "OFF") << std::endl;
            l_was_pressed = true;
        }
    } else {
        l_was_pressed = false;
    }
    
    static bool p_was_pressed = false;
    if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
        if (!p_was_pressed && m_Renderer) {
            bool planar = !m_Renderer->isDiskModePlanar();
            m_Renderer->setDiskModePlanar(planar);
            std::cout << "Disk Mode: " << (planar ? "PLANAR" : "VOLUMETRIC") << std::endl;
            p_was_pressed = true;
        }
    } else {
        p_was_pressed = false;
    }
}

void Application::handleMovement(GLFWwindow* window, float speed) {
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) m_ObserverPosition(1) -= speed;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) m_ObserverPosition(1) += speed;
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) m_ObserverPosition(3) -= speed * 0.1f;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) m_ObserverPosition(3) += speed * 0.1f;
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) m_ObserverPosition(2) -= speed * 0.1f;
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) m_ObserverPosition(2) += speed * 0.1f;
    
    setObserverPosition(m_ObserverPosition);
}

void Application::handleMouse(double xpos, double ypos) {
    if (!m_MouseLook) return;
    
    if (m_FirstMouse) {
        m_LastMouseX = xpos;
        m_LastMouseY = ypos;
        m_FirstMouse = false;
    }
    
    float xoffset = static_cast<float>(xpos - m_LastMouseX);
    float yoffset = static_cast<float>(m_LastMouseY - ypos);
    
    m_LastMouseX = xpos;
    m_LastMouseY = ypos;
    
    float sensitivity = 0.005f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;
    
    m_CameraYaw += xoffset;
    m_CameraPitch += yoffset;
    
    m_CameraPitch = std::clamp(m_CameraPitch, static_cast<float>(-M_PI/2) + 0.01f, static_cast<float>(M_PI/2) - 0.01f);
    
    if (m_CameraYaw > M_PI) m_CameraYaw -= 2.0f * static_cast<float>(M_PI);
    if (m_CameraYaw < -M_PI) m_CameraYaw += 2.0f * static_cast<float>(M_PI);
}

void Application::updateObserver(float deltaTime) {
    // Simple time update - geodesic camera removed
    m_ObserverPosition(0) += deltaTime;
}

void Application::run() {
    auto lastTime = std::chrono::high_resolution_clock::now();
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "  DEPRECATED: Interactive Mode" << std::endl;
    std::cout << "  Use batch rendering for production." << std::endl;
    std::cout << "========================================\n" << std::endl;
    std::cout << "Controls:" << std::endl;
    std::cout << "  TAB - Toggle mouse look" << std::endl;
    std::cout << "  WASD - Move in r/phi directions" << std::endl;
    std::cout << "  QE - Move in theta direction" << std::endl;
    std::cout << "  Shift - Fast movement" << std::endl;
    std::cout << "  Ctrl - Slow movement" << std::endl;
    std::cout << "  B - Toggle bloom" << std::endl;
    std::cout << "  L - Toggle lens flare" << std::endl;
    std::cout << "  P - Toggle disk mode" << std::endl;
    
    while (!m_Window->shouldClose()) {
        auto currentTime = std::chrono::high_resolution_clock::now();
        float deltaTime = std::chrono::duration<float>(currentTime - lastTime).count();
        lastTime = currentTime;
        
        glClearColor(0.05f, 0.05f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        m_Window->pollEvents();
        handleInput();
        updateObserver(deltaTime);

        if (m_CurrentMetric && m_Renderer) {
            try {
                m_Renderer->render(m_CurrentMetric, m_ObserverPosition, m_ObserverVelocity, m_CameraYaw, m_CameraPitch, m_CameraFOV);
            } catch (const std::exception& e) {
                std::cerr << "Rendering error: " << e.what() << std::endl;
            }
        }

        if (m_UIManager) {
            m_UIManager->beginFrame();
            m_UIManager->setObserverPosition(m_ObserverPosition);
            m_UIManager->setObserverVelocity(m_ObserverVelocity);
            m_UIManager->setCameraFOV(m_CameraFOV);
            m_UIManager->render();
            m_UIManager->endFrame();
        }

        m_Window->swapBuffers();
    }
    
    std::cout << "Shutting down Sirius" << std::endl;
}