// The view command: OpenGL window with real-time progressive black hole
// rendering and orbital camera controls. Ported from CRCL006A.cpp.

#include "sirius/app/cli/view_command.h"

#include "sirius/app/cli/cli_output.h"
#include "sirius/app/viewer/interactive_viewer.h"

// Prevent GLFW from including a platform OpenGL header before GLAD owns the
// function declarations. This keeps the contract correct even when the format
// authority alphabetizes third-party includes.
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include <glad/glad.h>

#include <cstring>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <vector>

namespace sirius::app {

namespace cli = cli_output;

namespace {

// GLFW's C callback API cannot take closures, so viewer state the callbacks need
// lives at file scope. One viewer runs at a time.
InteractiveViewer* g_viewer = nullptr;
bool g_mouse_dragging = false;
GLuint g_texture = 0;
bool g_texture_needs_update = false;
std::mutex g_frame_mutex;
std::vector<float> g_frame_data;
int g_frame_width = 0;
int g_frame_height = 0;

void KeyCallback(GLFWwindow* /*window*/, int key, int /*scancode*/, int action, int /*mods*/) {
    if (g_viewer) {
        if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
            g_viewer->Stop();
            return;
        }
        g_viewer->ProcessKey(key, action);
    }
}

void CursorPosCallback(GLFWwindow* /*window*/, double xpos, double ypos) {
    if (g_viewer) {
        g_viewer->ProcessMouseMove(xpos, ypos, g_mouse_dragging);
    }
}

void MouseButtonCallback(GLFWwindow* /*window*/, int button, int action, int /*mods*/) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        g_mouse_dragging = (action == GLFW_PRESS);
    }
}

void ScrollCallback(GLFWwindow* /*window*/, double /*xoffset*/, double yoffset) {
    if (g_viewer) {
        g_viewer->ProcessScroll(yoffset);
    }
}

void FramebufferSizeCallback(GLFWwindow* /*window*/, int /*width*/, int /*height*/) {
    // Handled in the main loop via glfwGetFramebufferSize.
}

// Frame callback, invoked from the render thread; snapshots the frame for the
// GL thread to upload.
void FrameCallback(const float* data, int width, int height) {
    std::lock_guard<std::mutex> lock(g_frame_mutex);
    size_t size = static_cast<size_t>(width) * static_cast<size_t>(height) * 4;
    g_frame_data.resize(size);
    std::memcpy(g_frame_data.data(), data, size * sizeof(float));
    g_frame_width = width;
    g_frame_height = height;
    g_texture_needs_update = true;
}

}  // namespace

std::string ViewCommand::Usage() const {
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

bool ViewCommand::ParseArgs(const std::vector<std::string>& args, const GlobalOptions& /*globals*/,
                            SiriusConfig& config) {
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
        } else if (arg == "--help") {
            std::cout << Usage();
            return false;
        }
    }
    return true;
}

int ViewCommand::Execute(const std::vector<std::string>& args, const GlobalOptions& globals,
                         SiriusConfig& config) {
    // Viewer defaults override the config defaults.
    config.render.width = 1280;
    config.render.height = 720;
    config.metric.spin = 0.9;
    config.observer.distance = 50.0;
    config.observer.inclination = 75.0;  // Degrees.
    config.observer.fov = 60.0;

    if (!ParseArgs(args, globals, config)) {
        return 0;  // Help was shown.
    }

    if (!glfwInit()) {
        cli::Error("Failed to initialize GLFW");
        return 1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    std::string title =
        "Sirius - Interactive Black Hole Viewer (a=" + std::to_string(config.metric.spin) + ")";
    GLFWwindow* window = glfwCreateWindow(config.render.width, config.render.height, title.c_str(),
                                          nullptr, nullptr);
    if (!window) {
        cli::Error("Failed to create GLFW window");
        glfwTerminate();
        return 1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);  // VSync.

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        cli::Error("Failed to initialize GLAD");
        glfwTerminate();
        return 1;
    }

    cli::Info("OpenGL " + std::string((const char*)glGetString(GL_VERSION)));

    glfwSetKeyCallback(window, KeyCallback);
    glfwSetCursorPosCallback(window, CursorPosCallback);
    glfwSetMouseButtonCallback(window, MouseButtonCallback);
    glfwSetScrollCallback(window, ScrollCallback);
    glfwSetFramebufferSizeCallback(window, FramebufferSizeCallback);

    glGenTextures(1, &g_texture);
    glBindTexture(GL_TEXTURE_2D, g_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    InteractiveViewer viewer;
    g_viewer = &viewer;

    ViewerConfig view_config;
    view_config.preview_width = config.render.width / 4;
    view_config.preview_height = config.render.height / 4;
    view_config.final_width = config.render.width;
    view_config.final_height = config.render.height;
    view_config.blackHoleSpin = config.metric.spin;
    view_config.observerDistance = config.observer.distance;
    view_config.observerInclination = config.inclinationRadians();
    view_config.enableDisk = true;

    viewer.Initialise(view_config);
    viewer.AttachWindow(window);
    viewer.SetFrameCallback(FrameCallback);

    cli::Success("Interactive viewer started");
    cli::Info("Controls: WASD to move, mouse drag to look, scroll to zoom, ESC to quit");

    viewer.Start();

    // Fullscreen quad.
    float quad_vertices[] = {
        // pos      // tex
        -1.0f, -1.0f, 0.0f, 1.0f, 1.0f,  -1.0f, 1.0f, 1.0f,
        1.0f,  1.0f,  1.0f, 0.0f, -1.0f, 1.0f,  0.0f, 0.0f,
    };
    unsigned int quad_indices[] = {0, 1, 2, 2, 3, 0};

    GLuint vao, vbo, ebo;
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glGenBuffers(1, &ebo);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quad_vertices), quad_vertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(quad_indices), quad_indices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // Passthrough blit + Reinhard tonemap (matches the legacy inline shader; the
    // RDSD00* files under viewer/shaders/ are the file-loadable equivalents).
    const char* vertex_shader_src = R"(
        #version 330 core
        layout(location = 0) in vec2 aPos;
        layout(location = 1) in vec2 aTexCoord;
        out vec2 TexCoord;
        void main() {
            gl_Position = vec4(aPos, 0.0, 1.0);
            TexCoord = aTexCoord;
        }
    )";

    const char* fragment_shader_src = R"(
        #version 330 core
        in vec2 TexCoord;
        out vec4 FragColor;
        uniform sampler2D screenTexture;
        void main() {
            vec3 color = texture(screenTexture, TexCoord).rgb;
            color = color / (color + vec3(1.0));
            color = pow(color, vec3(1.0/2.2));
            FragColor = vec4(color, 1.0);
        }
    )";

    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, &vertex_shader_src, nullptr);
    glCompileShader(vertex_shader);

    GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, &fragment_shader_src, nullptr);
    glCompileShader(fragment_shader);

    GLuint shader_program = glCreateProgram();
    glAttachShader(shader_program, vertex_shader);
    glAttachShader(shader_program, fragment_shader);
    glLinkProgram(shader_program);

    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);

    auto last_time = glfwGetTime();

    while (!glfwWindowShouldClose(window) && viewer.IsRunning()) {
        auto current_time = glfwGetTime();
        float dt = static_cast<float>(current_time - last_time);
        last_time = current_time;

        viewer.UpdateCamera(dt);

        glfwPollEvents();

        {
            std::lock_guard<std::mutex> lock(g_frame_mutex);
            if (g_texture_needs_update && !g_frame_data.empty()) {
                glBindTexture(GL_TEXTURE_2D, g_texture);
                glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, g_frame_width, g_frame_height, 0,
                             GL_RGBA, GL_FLOAT, g_frame_data.data());
                g_texture_needs_update = false;
            }
        }

        int fb_width, fb_height;
        glfwGetFramebufferSize(window, &fb_width, &fb_height);
        glViewport(0, 0, fb_width, fb_height);
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glUseProgram(shader_program);
        glBindTexture(GL_TEXTURE_2D, g_texture);
        glBindVertexArray(vao);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        glfwSwapBuffers(window);

        const auto& refinement = viewer.GetRefinementState();
        const auto& camera = viewer.GetCameraState();
        std::ostringstream status;
        status << "Sirius - r=" << std::fixed << std::setprecision(1) << camera.r << "M, Level "
               << refinement.current_level + 1 << "/" << view_config.refinement_levels << " ("
               << refinement.current_width << "x" << refinement.current_height << ")";
        if (refinement.complete) {
            status << " [Complete]";
        }
        glfwSetWindowTitle(window, status.str().c_str());
    }

    viewer.Stop();
    g_viewer = nullptr;

    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
    glDeleteBuffers(1, &ebo);
    glDeleteProgram(shader_program);
    glDeleteTextures(1, &g_texture);

    glfwDestroyWindow(window);
    glfwTerminate();

    cli::Success("Viewer closed");
    return 0;
}

}  // namespace sirius::app
