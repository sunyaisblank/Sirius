// CRCL006A.cpp - View Command (Interactive Viewer)
// Component ID: CRCL006A (Cli/ViewCommand)
//
// Launches an OpenGL window with real-time progressive black hole rendering.
// Camera controls: WASD movement, mouse drag for look, scroll for zoom.

#include "CRCL006A.h"
#include "CRCL005A.h"  // Output utilities
#include "INVW001A.h"  // Interactive viewer

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <sstream>

namespace Sirius::Cli {

// Required by GLFW C callback API — closures cannot be used as function pointers
static InteractiveViewer* g_Viewer = nullptr;
static bool g_MouseDragging = false;
static GLuint g_Texture = 0;
static bool g_TextureNeedsUpdate = false;
static std::mutex g_FrameMutex;
static std::vector<float> g_FrameData;
static int g_FrameWidth = 0;
static int g_FrameHeight = 0;

//==============================================================================
// GLFW Callbacks
//==============================================================================
static void keyCallback(GLFWwindow* /*window*/, int key, int /*scancode*/, int action, int /*mods*/) {
    if (g_Viewer) {
        // ESC to quit
        if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
            g_Viewer->stop();
            return;
        }
        g_Viewer->processKey(key, action);
    }
}

static void cursorPosCallback(GLFWwindow* /*window*/, double xpos, double ypos) {
    if (g_Viewer) {
        g_Viewer->processMouseMove(xpos, ypos, g_MouseDragging);
    }
}

static void mouseButtonCallback(GLFWwindow* /*window*/, int button, int action, int /*mods*/) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        g_MouseDragging = (action == GLFW_PRESS);
    }
}

static void scrollCallback(GLFWwindow* /*window*/, double /*xoffset*/, double yoffset) {
    if (g_Viewer) {
        g_Viewer->processScroll(yoffset);
    }
}

static void framebufferSizeCallback(GLFWwindow* /*window*/, int /*width*/, int /*height*/) {
    // Handled in main loop
}

//==============================================================================
// Frame callback (called from render thread)
//==============================================================================
static void frameCallback(const float* data, int width, int height) {
    std::lock_guard<std::mutex> lock(g_FrameMutex);
    size_t size = width * height * 4;
    g_FrameData.resize(size);
    std::memcpy(g_FrameData.data(), data, size * sizeof(float));
    g_FrameWidth = width;
    g_FrameHeight = height;
    g_TextureNeedsUpdate = true;
}

//==============================================================================
// ViewCommand Implementation
//==============================================================================

std::string ViewCommand::usage() const {
    std::ostringstream ss;
    ss << "Usage: sirius view [options]\n\n"
       << "Launch interactive viewer with real-time progressive rendering.\n\n"
       << "Controls:\n"
       << "  W/S          Move closer/further from black hole\n"
       << "  A/D          Orbit left/right (azimuthal)\n"
       << "  Q/E          Tilt up/down (inclination)\n"
       << "  Mouse Drag   Look around\n"
       << "  Scroll       Zoom (adjust FOV)\n"
       << "  ESC          Exit viewer\n\n"
       << "Options:\n"
       << "  --width <w>         Window width (default: 1280)\n"
       << "  --height <h>        Window height (default: 720)\n"
       << "  --spin <a>          Black hole spin a/M (default: 0.9)\n"
       << "  --distance <r>      Observer distance (default: 50M)\n"
       << "  --inclination <deg> Observer inclination (default: 75)\n"
       << "  --fov <deg>         Field of view (default: 60)\n"
       << "  --no-disk           Disable accretion disk\n"
       << "  --jets              Enable relativistic jets\n";
    return ss.str();
}

bool ViewCommand::parseArgs(const std::vector<std::string>& args,
                           const Configuration::GlobalOptions& /*globals*/,
                           Configuration::SiriusConfig& config) {
    for (size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];

        if (arg == "--width" && i + 1 < args.size()) {
            config.render.width = std::stoi(args[++i]);
        } else if (arg == "--height" && i + 1 < args.size()) {
            config.render.height = std::stoi(args[++i]);
        } else if (arg == "--spin" && i + 1 < args.size()) {
            config.metric.spin = std::stod(args[++i]);
        } else if (arg == "--distance" && i + 1 < args.size()) {
            config.observer.distance = std::stod(args[++i]);
        } else if (arg == "--inclination" && i + 1 < args.size()) {
            config.observer.inclination = std::stod(args[++i]);
        } else if (arg == "--fov" && i + 1 < args.size()) {
            config.observer.fov = std::stod(args[++i]);
        } else if (arg == "--help" || arg == "-h") {
            std::cout << usage();
            return false;
        }
    }
    return true;
}

int ViewCommand::execute(const std::vector<std::string>& args,
                        const Configuration::GlobalOptions& globals,
                        Configuration::SiriusConfig& config) {
    // Set viewer defaults (override config defaults)
    config.render.width = 1280;
    config.render.height = 720;
    config.metric.spin = 0.9;
    config.observer.distance = 50.0;
    config.observer.inclination = 75.0;  // degrees
    config.observer.fov = 60.0;

    // Parse arguments (overrides defaults)
    if (!parseArgs(args, globals, config)) {
        return 0;  // Help was shown
    }

    // Initialize GLFW
    if (!glfwInit()) {
        Output::error("Failed to initialize GLFW");
        return 1;
    }

    // Request OpenGL 3.3 core
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // Create window
    std::string title = "Sirius - Interactive Black Hole Viewer (a=" +
                        std::to_string(config.metric.spin) + ")";
    GLFWwindow* window = glfwCreateWindow(config.render.width, config.render.height,
                                           title.c_str(), nullptr, nullptr);
    if (!window) {
        Output::error("Failed to create GLFW window");
        glfwTerminate();
        return 1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);  // VSync

    // Load OpenGL functions
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        Output::error("Failed to initialize GLAD");
        glfwTerminate();
        return 1;
    }

    Output::info("OpenGL " + std::string((const char*)glGetString(GL_VERSION)));

    // Set callbacks
    glfwSetKeyCallback(window, keyCallback);
    glfwSetCursorPosCallback(window, cursorPosCallback);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetScrollCallback(window, scrollCallback);
    glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);

    // Create texture for display
    glGenTextures(1, &g_Texture);
    glBindTexture(GL_TEXTURE_2D, g_Texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    // Create and configure interactive viewer
    InteractiveViewer viewer;
    g_Viewer = &viewer;

    ViewerConfig viewConfig;
    viewConfig.preview_width = config.render.width / 4;
    viewConfig.preview_height = config.render.height / 4;
    viewConfig.final_width = config.render.width;
    viewConfig.final_height = config.render.height;
    viewConfig.blackHoleSpin = config.metric.spin;
    viewConfig.observerDistance = config.observer.distance;
    viewConfig.observerInclination = config.inclinationRadians();
    viewConfig.enableDisk = true;

    viewer.initialize(viewConfig);
    viewer.attachWindow(window);
    viewer.setFrameCallback(frameCallback);

    Output::success("Interactive viewer started");
    Output::info("Controls: WASD to move, mouse drag to look, scroll to zoom, ESC to quit");

    // Start rendering
    viewer.start();

    // Create simple fullscreen quad VAO
    float quadVertices[] = {
        // pos      // tex
        -1.0f, -1.0f, 0.0f, 1.0f,
         1.0f, -1.0f, 1.0f, 1.0f,
         1.0f,  1.0f, 1.0f, 0.0f,
        -1.0f,  1.0f, 0.0f, 0.0f
    };
    unsigned int quadIndices[] = { 0, 1, 2, 2, 3, 0 };

    GLuint VAO, VBO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(quadIndices), quadIndices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // Simple shader
    const char* vertexShaderSrc = R"(
        #version 330 core
        layout(location = 0) in vec2 aPos;
        layout(location = 1) in vec2 aTexCoord;
        out vec2 TexCoord;
        void main() {
            gl_Position = vec4(aPos, 0.0, 1.0);
            TexCoord = aTexCoord;
        }
    )";

    const char* fragmentShaderSrc = R"(
        #version 330 core
        in vec2 TexCoord;
        out vec4 FragColor;
        uniform sampler2D screenTexture;
        void main() {
            vec3 color = texture(screenTexture, TexCoord).rgb;
            // Simple tonemap (Reinhard)
            color = color / (color + vec3(1.0));
            // Gamma correction
            color = pow(color, vec3(1.0/2.2));
            FragColor = vec4(color, 1.0);
        }
    )";

    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSrc, nullptr);
    glCompileShader(vertexShader);

    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSrc, nullptr);
    glCompileShader(fragmentShader);

    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    // Main loop
    auto lastTime = glfwGetTime();

    while (!glfwWindowShouldClose(window) && viewer.isRunning()) {
        // Calculate delta time
        auto currentTime = glfwGetTime();
        float dt = static_cast<float>(currentTime - lastTime);
        lastTime = currentTime;

        // Update camera
        viewer.updateCamera(dt);

        // Poll events
        glfwPollEvents();

        // Update texture if new frame available
        {
            std::lock_guard<std::mutex> lock(g_FrameMutex);
            if (g_TextureNeedsUpdate && !g_FrameData.empty()) {
                glBindTexture(GL_TEXTURE_2D, g_Texture);
                glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, g_FrameWidth, g_FrameHeight,
                            0, GL_RGBA, GL_FLOAT, g_FrameData.data());
                g_TextureNeedsUpdate = false;
            }
        }

        // Render
        int fbWidth, fbHeight;
        glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
        glViewport(0, 0, fbWidth, fbHeight);
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glUseProgram(shaderProgram);
        glBindTexture(GL_TEXTURE_2D, g_Texture);
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        // Swap buffers
        glfwSwapBuffers(window);

        // Display status in title
        const auto& refinement = viewer.getRefinementState();
        const auto& camera = viewer.getCameraState();
        std::ostringstream status;
        status << "Sirius - r=" << std::fixed << std::setprecision(1) << camera.r
               << "M, Level " << refinement.current_level + 1 << "/"
               << viewConfig.refinement_levels
               << " (" << refinement.current_width << "x" << refinement.current_height << ")";
        if (refinement.complete) {
            status << " [Complete]";
        }
        glfwSetWindowTitle(window, status.str().c_str());
    }

    // Cleanup
    viewer.stop();
    g_Viewer = nullptr;

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
    glDeleteProgram(shaderProgram);
    glDeleteTextures(1, &g_Texture);

    glfwDestroyWindow(window);
    glfwTerminate();

    Output::success("Viewer closed");
    return 0;
}

} // namespace Sirius::Cli
