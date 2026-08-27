// The view command: OpenGL window with real-time progressive black hole
// rendering and orbital camera controls. Ported from CRCL006A.cpp.

#include "sirius/app/cli/view_command.h"

#include "sirius/app/cli/cli_output.h"
#include "sirius/app/config/config_loader.h"
#include "sirius/app/config/session_config_adapter.h"
#include "sirius/app/viewer/interactive_viewer.h"
#include "sirius/base/resource_locator.h"

// Prevent GLFW from including a platform OpenGL header before GLAD owns the
// function declarations. This keeps the contract correct even when the format
// authority alphabetizes third-party includes.
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include <glad/glad.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace sirius::app {

namespace cli = cli_output;

namespace {

struct DisplayFrame {
    std::mutex mutex;
    std::vector<float> pixels;
    int width = 0;
    int height = 0;
    bool needs_upload = false;
};

struct ViewerCallbackState {
    InteractiveViewer* viewer = nullptr;
    std::ostream* transcript = nullptr;
    std::mutex transcript_mutex;
    bool mouse_dragging = false;
};

ViewerCallbackState* GetCallbackState(GLFWwindow* window) {
    return static_cast<ViewerCallbackState*>(glfwGetWindowUserPointer(window));
}

std::string_view ActionName(int action) {
    switch (action) {
        case GLFW_RELEASE:
            return "release";
        case GLFW_PRESS:
            return "press";
        case GLFW_REPEAT:
            return "repeat";
        default:
            return "unknown";
    }
}

std::string LoadShaderSource(const std::string& relative_path) {
    const auto path = base::ResolveResource(relative_path);
    if (!path) {
        throw std::runtime_error("required viewer shader is missing: " + relative_path);
    }
    std::ifstream file(*path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("required viewer shader is unreadable: " + path->string());
    }
    return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

GLuint CompileShader(GLenum type, const std::string& source, const std::string& label) {
    const GLuint shader = glCreateShader(type);
    const char* text = source.c_str();
    glShaderSource(shader, 1, &text, nullptr);
    glCompileShader(shader);

    GLint compiled = GL_FALSE;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
    if (compiled == GL_TRUE) return shader;

    GLint length = 0;
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);
    std::vector<char> log(static_cast<std::size_t>(std::max(length, 1)), '\0');
    glGetShaderInfoLog(shader, length, nullptr, log.data());
    glDeleteShader(shader);
    throw std::runtime_error("viewer shader compile failed (" + label + "): " + log.data());
}

GLuint BuildViewerShaderProgram() {
    const std::string vertex_source = LoadShaderSource("shaders/RDSD003A.vert");
    const std::string fragment_source = LoadShaderSource("shaders/RDSD003A.frag");
    const GLuint vertex = CompileShader(GL_VERTEX_SHADER, vertex_source, "RDSD003A.vert");
    GLuint fragment = 0;
    try {
        fragment = CompileShader(GL_FRAGMENT_SHADER, fragment_source, "RDSD003A.frag");
    } catch (...) {
        glDeleteShader(vertex);
        throw;
    }

    const GLuint program = glCreateProgram();
    glAttachShader(program, vertex);
    glAttachShader(program, fragment);
    glLinkProgram(program);
    glDeleteShader(vertex);
    glDeleteShader(fragment);

    GLint linked = GL_FALSE;
    glGetProgramiv(program, GL_LINK_STATUS, &linked);
    if (linked == GL_TRUE) return program;

    GLint length = 0;
    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &length);
    std::vector<char> log(static_cast<std::size_t>(std::max(length, 1)), '\0');
    glGetProgramInfoLog(program, length, nullptr, log.data());
    glDeleteProgram(program);
    throw std::runtime_error("viewer shader link failed: " + std::string(log.data()));
}

const std::string& RequireValue(const std::vector<std::string>& args, std::size_t& index,
                                const std::string& option) {
    if (++index >= args.size()) {
        throw std::invalid_argument(option + " requires a value");
    }
    return args[index];
}

int ParseInteger(const std::string& text) {
    std::size_t consumed = 0;
    const int value = std::stoi(text, &consumed);
    if (consumed != text.size()) throw std::invalid_argument("trailing characters");
    return value;
}

double ParseDouble(const std::string& text) {
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != text.size() || !std::isfinite(value)) {
        throw std::invalid_argument("expected one finite number");
    }
    return value;
}

void KeyCallback(GLFWwindow* window, int key, int /*scancode*/, int action, int /*mods*/) {
    if (ViewerCallbackState* state = GetCallbackState(window); state != nullptr) {
        if (state->transcript != nullptr) {
            std::lock_guard<std::mutex> lock(state->transcript_mutex);
            *state->transcript << "keyboard-callback key=" << key
                               << " action=" << ActionName(action) << '\n'
                               << std::flush;
        }
        if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
            state->viewer->Stop();
            glfwSetWindowShouldClose(window, GLFW_TRUE);
            return;
        }
        state->viewer->ProcessKey(key, action);
    }
}

void CursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
    if (ViewerCallbackState* state = GetCallbackState(window); state != nullptr) {
        if (state->transcript != nullptr) {
            std::lock_guard<std::mutex> lock(state->transcript_mutex);
            *state->transcript << "pointer-callback kind=cursor x=" << xpos << " y=" << ypos
                               << " dragging=" << std::boolalpha << state->mouse_dragging << '\n'
                               << std::flush;
        }
        state->viewer->ProcessMouseMove(xpos, ypos, state->mouse_dragging);
    }
}

void MouseButtonCallback(GLFWwindow* window, int button, int action, int /*mods*/) {
    if (ViewerCallbackState* state = GetCallbackState(window); state != nullptr) {
        if (state->transcript != nullptr) {
            std::lock_guard<std::mutex> lock(state->transcript_mutex);
            *state->transcript << "pointer-callback kind=button button=" << button
                               << " action=" << ActionName(action) << '\n'
                               << std::flush;
        }
        if (button == GLFW_MOUSE_BUTTON_LEFT) {
            state->mouse_dragging = (action == GLFW_PRESS);
        }
    }
}

void ScrollCallback(GLFWwindow* window, double /*xoffset*/, double yoffset) {
    if (ViewerCallbackState* state = GetCallbackState(window); state != nullptr) {
        if (state->transcript != nullptr) {
            std::lock_guard<std::mutex> lock(state->transcript_mutex);
            *state->transcript << "pointer-callback kind=scroll yoffset=" << yoffset << '\n'
                               << std::flush;
        }
        state->viewer->ProcessScroll(yoffset);
    }
}

void FramebufferSizeCallback(GLFWwindow* /*window*/, int /*width*/, int /*height*/) {
    // Handled in the main loop via glfwGetFramebufferSize.
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
       << "  --width <w>         Window width (default: 1920)\n"
       << "  --height <h>        Window height (default: 1080)\n"
       << "  --spin <a>          Black hole spin a/M (default: 0)\n"
       << "  --distance <r>      Observer distance (default: 50M)\n"
       << "  --inclination <deg> Observer inclination (default: 90)\n"
       << "  --fov <deg>         Field of view (default: 60)\n"
       << "  --no-disk           Disable accretion disk\n"
       << "  --jets              Enable relativistic jets\n"
       << "  --cpu                Pin the CPU render path\n"
       << "  --gpu                Require the Vulkan render path\n"
       << "  --backend <name>     Select auto, cpu, or vulkan\n";
    return ss.str();
}

bool ViewCommand::ParseArgs(const std::vector<std::string>& args, const GlobalOptions& /*globals*/,
                            SiriusConfig& config) {
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];

        try {
            if (arg == "--width") {
                config.render.width = ParseInteger(RequireValue(args, i, arg));
            } else if (arg == "--height") {
                config.render.height = ParseInteger(RequireValue(args, i, arg));
            } else if (arg == "--spin") {
                config.metric.spin = ParseDouble(RequireValue(args, i, arg));
                config.metric.name = config.metric.spin == 0.0 ? "Schwarzschild" : "Kerr";
            } else if (arg == "--distance") {
                config.observer.distance = ParseDouble(RequireValue(args, i, arg));
            } else if (arg == "--inclination") {
                config.observer.inclination = ParseDouble(RequireValue(args, i, arg));
            } else if (arg == "--fov") {
                config.observer.fov = ParseDouble(RequireValue(args, i, arg));
            } else if (arg == "--no-disk") {
                config.disk_enabled = false;
            } else if (arg == "--jets") {
                jets_enabled_ = true;
            } else if (arg == "--cpu") {
                config.backend.preferred = "cpu";
            } else if (arg == "--gpu") {
                config.backend.preferred = "vulkan";
            } else if (arg == "--backend") {
                config.backend.preferred = RequireValue(args, i, arg);
            } else {
                cli::Error("Unknown view option: " + arg);
                return false;
            }
        } catch (const std::exception& e) {
            cli::Error("Invalid value for " + arg + ": " + e.what());
            return false;
        }
    }
    return true;
}

int ViewCommand::Execute(const std::vector<std::string>& args, const GlobalOptions& globals,
                         SiriusConfig& config) {
    jets_enabled_ = false;
    DisplayFrame display_frame;
    std::ofstream input_transcript;

    if (!ParseArgs(args, globals, config)) {
        return 1;
    }
    const auto errors = ConfigLoader::Validate(config);
    if (!errors.empty()) {
        for (const auto& error : errors) cli::Error(error);
        return 1;
    }

    if (const char* path = std::getenv("SIRIUS_VIEWER_INPUT_TRANSCRIPT"); path != nullptr) {
        if (*path == '\0') {
            cli::Error("SIRIUS_VIEWER_INPUT_TRANSCRIPT must name a file when set");
            return 1;
        }
        input_transcript.open(path, std::ios::out | std::ios::trunc);
        if (!input_transcript) {
            cli::Error("Could not open the requested viewer input transcript");
            return 1;
        }
    }

    auto adapted = MakeSessionConfig(config);
    if (!adapted) {
        cli::Error(adapted.error().Description());
        return 1;
    }
    render::SessionConfig resolved = std::move(*adapted);
    if (resolved.metric_id != core::MetricId::Schwarzschild &&
        resolved.metric_id != core::MetricId::Kerr) {
        cli::Error("The interactive viewer currently represents Schwarzschild and Kerr only");
        return 1;
    }
    if (jets_enabled_ && resolved.backend == render::RenderBackend::Vulkan) {
        if (config.backend.preferred == "vulkan") {
            cli::Error("The Vulkan viewer does not represent relativistic jets; use --cpu");
            return 1;
        }
        resolved.backend = render::RenderBackend::Cpu;
        cli::Info("Viewer jets require the CPU render path; backend auto selected CPU");
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
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    const GLubyte* gl_version = glGetString(GL_VERSION);
    if (gl_version == nullptr) {
        cli::Error("OpenGL context did not report a version string");
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }
    cli::Info("OpenGL " + std::string(reinterpret_cast<const char*>(gl_version)));

    GLuint texture = 0;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    InteractiveViewer viewer;
    ViewerCallbackState callback_state{
        .viewer = &viewer,
        .transcript = input_transcript ? &input_transcript : nullptr,
    };
    glfwSetWindowUserPointer(window, &callback_state);
    glfwSetKeyCallback(window, KeyCallback);
    glfwSetCursorPosCallback(window, CursorPosCallback);
    glfwSetMouseButtonCallback(window, MouseButtonCallback);
    glfwSetScrollCallback(window, ScrollCallback);
    glfwSetFramebufferSizeCallback(window, FramebufferSizeCallback);
    if (callback_state.transcript != nullptr) {
        std::lock_guard<std::mutex> lock(callback_state.transcript_mutex);
        *callback_state.transcript
            << "window-created opengl-version=" << reinterpret_cast<const char*>(gl_version) << '\n'
            << std::flush;
    }

    ViewerConfig view_config;
    view_config.session_template = resolved;
    view_config.preview_width =
        std::min(config.render.width, std::max(64, config.render.width / 4));
    view_config.preview_height =
        std::min(config.render.height, std::max(64, config.render.height / 4));
    view_config.final_width = config.render.width;
    view_config.final_height = config.render.height;
    view_config.black_hole_mass = resolved.black_hole_mass;
    view_config.black_hole_spin = resolved.black_hole_spin;
    view_config.metric_id = resolved.metric_id;
    view_config.observer_distance = config.observer.distance;
    view_config.observer_inclination = config.inclination_radians();
    view_config.observer_azimuth = config.observer.azimuth * core::constants::math::kPi / 180.0;
    view_config.observer_fov = static_cast<float>(config.observer.fov);
    view_config.enable_disk = config.disk_enabled;
    view_config.enable_volumetric = config.volumetric.enabled;
    view_config.enable_jets = jets_enabled_;
    view_config.backend = resolved.backend;

    if (!viewer.Initialise(view_config)) {
        cli::Error("Viewer configuration is outside the represented operating domain");
        glfwSetWindowUserPointer(window, nullptr);
        glDeleteTextures(1, &texture);
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }
    viewer.SetFrameCallback([&display_frame, &callback_state, backend = view_config.backend](
                                const float* data, int width, int height) {
        {
            std::lock_guard<std::mutex> lock(display_frame.mutex);
            const std::size_t size =
                static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 4;
            display_frame.pixels.resize(size);
            std::memcpy(display_frame.pixels.data(), data, size * sizeof(float));
            display_frame.width = width;
            display_frame.height = height;
            display_frame.needs_upload = true;
        }
        if (callback_state.transcript != nullptr) {
            std::lock_guard<std::mutex> lock(callback_state.transcript_mutex);
            *callback_state.transcript
                << "frame-published backend="
                << (backend == render::RenderBackend::Vulkan ? "Vulkan" : "Cpu")
                << " width=" << width << " height=" << height << '\n'
                << std::flush;
        }
    });

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

    GLuint shader_program = 0;
    try {
        shader_program = BuildViewerShaderProgram();
    } catch (const std::exception& error) {
        cli::Error(error.what());
        glfwSetWindowUserPointer(window, nullptr);
        glDeleteVertexArrays(1, &vao);
        glDeleteBuffers(1, &vbo);
        glDeleteBuffers(1, &ebo);
        glDeleteTextures(1, &texture);
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }

    if (!viewer.Start()) {
        cli::Error("Failed to start the viewer render thread");
        glfwSetWindowUserPointer(window, nullptr);
        glDeleteVertexArrays(1, &vao);
        glDeleteBuffers(1, &vbo);
        glDeleteBuffers(1, &ebo);
        glDeleteProgram(shader_program);
        glDeleteTextures(1, &texture);
        glfwDestroyWindow(window);
        glfwTerminate();
        return 1;
    }
    cli::Success("Interactive viewer started");
    cli::Info("Controls: WASD to move, mouse drag to look, scroll to zoom, ESC to quit");

    auto last_time = glfwGetTime();

    while (!glfwWindowShouldClose(window) && viewer.IsRunning()) {
        auto current_time = glfwGetTime();
        float dt = static_cast<float>(current_time - last_time);
        last_time = current_time;

        viewer.UpdateCamera(dt);

        glfwPollEvents();

        {
            std::lock_guard<std::mutex> lock(display_frame.mutex);
            if (display_frame.needs_upload && !display_frame.pixels.empty()) {
                glBindTexture(GL_TEXTURE_2D, texture);
                glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, display_frame.width,
                             display_frame.height, 0, GL_RGBA, GL_FLOAT,
                             display_frame.pixels.data());
                display_frame.needs_upload = false;
            }
        }

        int fb_width, fb_height;
        glfwGetFramebufferSize(window, &fb_width, &fb_height);
        glViewport(0, 0, fb_width, fb_height);
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glUseProgram(shader_program);
        glBindTexture(GL_TEXTURE_2D, texture);
        glBindVertexArray(vao);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        glfwSwapBuffers(window);

        const auto refinement = viewer.GetRefinementState();
        const auto camera = viewer.GetCameraState();
        std::ostringstream status;
        const int displayed_level =
            std::min(refinement.current_level + 1, view_config.refinement_levels);
        status << "Sirius - r=" << std::fixed << std::setprecision(1) << camera.r << "M, Level "
               << displayed_level << "/" << view_config.refinement_levels << " ("
               << refinement.current_width << "x" << refinement.current_height << ")";
        if (refinement.complete) {
            status << " [Complete]";
        }
        glfwSetWindowTitle(window, status.str().c_str());
    }

    viewer.Stop();
    const std::string viewer_error = viewer.GetLastError();
    glfwSetWindowUserPointer(window, nullptr);

    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
    glDeleteBuffers(1, &ebo);
    glDeleteProgram(shader_program);
    glDeleteTextures(1, &texture);

    glfwDestroyWindow(window);
    glfwTerminate();

    if (!viewer_error.empty()) {
        cli::Error(viewer_error);
        return 1;
    }
    cli::Success("Viewer closed");
    return 0;
}

}  // namespace sirius::app
