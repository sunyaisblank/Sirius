// =============================================================================
// Sirius.Render/RDRT001A.cpp - Main Renderer with OptiX Integration
// =============================================================================
// This renderer now uses NVIDIA OptiX for GPU ray tracing instead of GLSL 
// compute shaders. The OptiX backend provides RTX hardware acceleration.
// =============================================================================

#include "RDRT001A.h"

// Include CUDA runtime BEFORE RDOP003A.h to get proper float3/float4 types
#ifdef SIRIUS_HAS_OPTIX
#include <cuda_runtime.h>
#endif

#include "RDOP003A.h"




// stb_image for texture loading
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// Standard library includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <cstring>
#include <cmath>

// Filesystem for path resolution
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <limits.h>
#endif

// OpenGL includes (for display only)
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// GLM for matrix operations
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// Project includes
#include "PHMT000A.h"
#include <MTTN001A.h>

using namespace Sirius;

// =============================================================================
// OptiX C API (from RDOP001A.cu)
// =============================================================================
extern "C" {
    typedef void* SiriusOptixHandle;
    
    SiriusOptixHandle sirius_optix_create();
    void sirius_optix_destroy(SiriusOptixHandle handle);
    bool sirius_optix_initialize(SiriusOptixHandle handle, int width, int height);
    bool sirius_optix_create_pipeline(SiriusOptixHandle handle, const char* ptxPath);
    void sirius_optix_launch(SiriusOptixHandle handle, const Sirius::LaunchParams* params);
    void sirius_optix_update_display(SiriusOptixHandle handle);
    void sirius_optix_resize(SiriusOptixHandle handle, int width, int height);
    void sirius_optix_cleanup(SiriusOptixHandle handle);
    float* sirius_optix_get_frame_buffer(SiriusOptixHandle handle);
    void sirius_optix_register_gl_texture(SiriusOptixHandle handle, unsigned int glTexture);
    bool sirius_optix_is_initialized(SiriusOptixHandle handle);
    bool sirius_optix_upload_background(SiriusOptixHandle handle, const unsigned char* data, int width, int height);
    unsigned long long sirius_optix_get_background_texture(SiriusOptixHandle handle);
    void sirius_optix_upload_numerical_metric(SiriusOptixHandle handle, 
                                              const Sirius::NumericalMetricHostData* hostData,
                                              Sirius::NumericalMetricData* outDeviceData);
    void sirius_optix_set_metric_type(SiriusOptixHandle handle, int type);
    
    // Denoiser API
    bool sirius_optix_init_denoiser(SiriusOptixHandle handle);
    void sirius_optix_denoise(SiriusOptixHandle handle, float blendFactor);
    void sirius_optix_set_denoiser_enabled(SiriusOptixHandle handle, bool enabled);
}

// =============================================================================
// Renderer Implementation
// =============================================================================

Renderer::Renderer() 
    : m_ComputeProgram(0), m_ScreenProgram(0), m_Texture(0), m_Vao(0), m_MetricUBO(0), 
      m_BackgroundTexture(0), m_BackgroundWidth(0), m_BackgroundHeight(0), m_UseBackgroundTexture(false),
      m_Width(0), m_Height(0), m_OptixHandle(nullptr), m_OptixEnabled(false), m_FrameCount(0) {
}

Renderer::~Renderer() {
    cleanup();
}

// Get directory containing the executable
static std::string getExecutableDir() {
#ifdef _WIN32
    char result[MAX_PATH];
    DWORD count = GetModuleFileNameA(NULL, result, MAX_PATH);
    if (count > 0) {
        std::string path(result, count);
        size_t pos = path.find_last_of("\\/");
        if (pos != std::string::npos) {
            return path.substr(0, pos + 1);
        }
    }
#else
    char result[PATH_MAX];
    ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
    if (count != -1) {
        std::string path(result, count);
        size_t pos = path.find_last_of("/");
        if (pos != std::string::npos) {
            return path.substr(0, pos + 1);
        }
    }
#endif
    return "";
}

// Helper to resolve resource paths across different build layouts
static std::string resolvePath(const std::string& relativePath) {
    std::vector<std::string> pathsToTry;
    std::string baseDir = getExecutableDir();
    
    // Check if it's an absolute path or exists relative to CWD
    pathsToTry.push_back(relativePath);
    
    // Check from executable directory
    pathsToTry.push_back(baseDir + relativePath);
    
    // Check parent directory (VS builds: exe in Release/Debug, resources in parent)
    pathsToTry.push_back(baseDir + "../" + relativePath);
    
    // Check grandparent directory  
    pathsToTry.push_back(baseDir + "../../" + relativePath);
    
    // Check great-grandparent directory (for deep build structures)
    pathsToTry.push_back(baseDir + "../../../" + relativePath);
    
    // Check even deeper
    pathsToTry.push_back(baseDir + "../../../../" + relativePath);
    
    // For VS multi-config: exe in bin/Release, resources in bin/Sirius.Render
    // Need to go up one level from Release and into Sirius.Render
    if (relativePath.find("Sirius.Render") == 0) {
        pathsToTry.push_back(baseDir + "../" + relativePath);
    }
    
    for (const auto& path : pathsToTry) {
        if (std::ifstream(path).good()) {
            return path;
        }
    }
    
    // Debug: print all paths we tried
    std::cerr << "[Renderer] Failed to find: " << relativePath << std::endl;
    std::cerr << "[Renderer] Base directory: " << baseDir << std::endl;
    std::cerr << "[Renderer] Paths tried:" << std::endl;
    for (const auto& path : pathsToTry) {
        std::cerr << "  - " << path << std::endl;
    }
    
    // Return the most likely path for better error messages
    return baseDir + "../" + relativePath;
}

// =============================================================================
// Shader Helper Functions
// =============================================================================

// Helper: Load shader source from file or return fallback
static std::string loadShaderSource(const std::string& path, const std::string& fallback) {
    std::ifstream file(path);
    if (!file.is_open()) {
        return fallback;
    }
    return std::string((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
}

// Helper: Compile a shader and return its ID (0 on failure)
static GLuint compileShader(GLenum type, const std::string& source, const std::string& name) {
    GLuint shader = glCreateShader(type);
    const char* src = source.c_str();
    glShaderSource(shader, 1, &src, nullptr);
    glCompileShader(shader);
    
    GLint status;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
    if (status != GL_TRUE) {
        char buffer[1024];
        glGetShaderInfoLog(shader, 1024, NULL, buffer);
        glDeleteShader(shader);
        std::cerr << "[" << name << "] Shader compile error: " << buffer << std::endl;
        return 0;
    }
    return shader;
}

// Helper: Link vertex and fragment shaders into a program
static GLuint linkShaderProgram(GLuint vertShader, GLuint fragShader, const std::string& name) {
    GLuint program = glCreateProgram();
    glAttachShader(program, vertShader);
    glAttachShader(program, fragShader);
    glLinkProgram(program);
    
    GLint status;
    glGetProgramiv(program, GL_LINK_STATUS, &status);
    if (status != GL_TRUE) {
        char buffer[1024];
        glGetProgramInfoLog(program, 1024, NULL, buffer);
        glDeleteProgram(program);
        std::cerr << "[" << name << "] Shader link error: " << buffer << std::endl;
        return 0;
    }
    
    glDeleteShader(vertShader);
    glDeleteShader(fragShader);
    return program;
}

// Default passthrough shaders
static const char* PASSTHROUGH_VERTEX_SHADER = R"(
#version 330 core
layout (location = 0) in vec2 aPos;
layout (location = 1) in vec2 aTexCoord;
out vec2 TexCoord;
void main() {
    gl_Position = vec4(aPos, 0.0, 1.0);
    TexCoord = aTexCoord;
}
)";

static const char* PASSTHROUGH_FRAGMENT_SHADER = R"(
#version 330 core
in vec2 TexCoord;
out vec4 FragColor;
uniform sampler2D screenTexture;
void main() {
    FragColor = texture(screenTexture, TexCoord);
}
)";

void Renderer::init(int width, int height) {
    m_Width = width;
    m_Height = height;

    // Match OpenGL viewport to the initial window size
    glViewport(0, 0, m_Width, m_Height);
    
    // ==========================================================================
    // Initialize OptiX Backend
    // ==========================================================================
#ifdef SIRIUS_HAS_OPTIX
    std::cout << "Initializing OptiX backend..." << std::endl;
    
    m_OptixHandle = sirius_optix_create();
    if (m_OptixHandle && sirius_optix_initialize(m_OptixHandle, width, height)) {
        // Build PTX path
        std::string ptxPath = resolvePath("Sirius.Render/ptx/RDOP002A.ptx");
        std::cout << "Loading PTX from: " << ptxPath << std::endl;
        
        if (sirius_optix_create_pipeline(m_OptixHandle, ptxPath.c_str())) {
            m_OptixEnabled = true;
            std::cout << "OptiX backend initialized successfully!" << std::endl;
            
            // Initialize and enable AI denoiser for noise reduction
            sirius_optix_set_denoiser_enabled(m_OptixHandle, true);
            if (sirius_optix_init_denoiser(m_OptixHandle)) {
                std::cout << "OptiX AI Denoiser enabled!" << std::endl;
            }
        } else {
            std::cerr << "OptiX pipeline creation failed" << std::endl;
            sirius_optix_destroy(m_OptixHandle);
            m_OptixHandle = nullptr;
        }
    } else {
        std::cerr << "OptiX initialization failed" << std::endl;
        if (m_OptixHandle) {
            sirius_optix_destroy(m_OptixHandle);
            m_OptixHandle = nullptr;
        }
    }
#endif
    
    // ==========================================================================
    // Create OpenGL Resources (for display)
    // ==========================================================================
    
    // Create output texture for display
    glGenTextures(1, &m_Texture);
    glBindTexture(GL_TEXTURE_2D, m_Texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    
#ifdef SIRIUS_HAS_OPTIX
    if (m_OptixEnabled && m_OptixHandle) {
        sirius_optix_register_gl_texture(m_OptixHandle, m_Texture);
    }
#endif
    
    // Create VAO for screen quad
    glGenVertexArrays(1, &m_Vao);
    glBindVertexArray(m_Vao);
    
    // Screen quad vertices
    float quadVertices[] = {
        -1.0f, -1.0f, 0.0f, 0.0f,
         1.0f, -1.0f, 1.0f, 0.0f,
        -1.0f,  1.0f, 0.0f, 1.0f,
        
        -1.0f,  1.0f, 0.0f, 1.0f,
         1.0f, -1.0f, 1.0f, 0.0f,
         1.0f,  1.0f, 1.0f, 1.0f
    };
    
    GLuint quadVBO;
    glGenBuffers(1, &quadVBO);
    glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);
    
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    glEnableVertexAttribArray(1);
    
    glBindVertexArray(0);
    
    // Load display shader (simple fullscreen quad)
    loadScreenShader(resolvePath("Sirius.Render/RDSD003A.vert"), resolvePath("Sirius.Render/RDSD003A.frag"));
    
    // Load lens flare post-process shader (Phase 6.5 - DNGR Cinematic Mode)
    loadLensFlareShader(resolvePath("Sirius.Render/RDSD004A.vert"), resolvePath("Sirius.Render/RDSD004A.frag"));
    
    // Load bloom post-process shader (Phase 7 - Cinematic Visual Quality)
    loadBloomShader(resolvePath("Sirius.Render/RDSD005A.vert"), resolvePath("Sirius.Render/RDSD005A.frag"));
    
    std::cout << "Renderer initialized (OptiX: " << (m_OptixEnabled ? "ON" : "OFF") 
              << ") " << width << "x" << height << std::endl;
}

void Renderer::resize(int width, int height) {
    if (width <= 0 || height <= 0 || (width == m_Width && height == m_Height)) {
        return;
    }

    m_Width = width;
    m_Height = height;
    m_FrameCount = 0;  // Reset accumulation

    glViewport(0, 0, m_Width, m_Height);

    // Resize display texture
    if (m_Texture != 0) {
        glDeleteTextures(1, &m_Texture);
    }
    glGenTextures(1, &m_Texture);
    glBindTexture(GL_TEXTURE_2D, m_Texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

#ifdef SIRIUS_HAS_OPTIX
    if (m_OptixEnabled && m_OptixHandle) {
        sirius_optix_resize(m_OptixHandle, width, height);
        sirius_optix_register_gl_texture(m_OptixHandle, m_Texture);
    }
#endif

    std::cout << "Renderer resized to " << width << "x" << height << std::endl;
}

void Renderer::loadScreenShader(const std::string& vertPath, const std::string& fragPath) {
    std::string vertSource = loadShaderSource(vertPath, PASSTHROUGH_VERTEX_SHADER);
    std::string fragSource = loadShaderSource(fragPath, PASSTHROUGH_FRAGMENT_SHADER);
    
    GLuint vertShader = compileShader(GL_VERTEX_SHADER, vertSource, "Screen");
    GLuint fragShader = compileShader(GL_FRAGMENT_SHADER, fragSource, "Screen");
    
    if (vertShader == 0 || fragShader == 0) {
        throw std::runtime_error("Screen shader compilation failed");
    }
    
    m_ScreenProgram = linkShaderProgram(vertShader, fragShader, "Screen");
    if (m_ScreenProgram == 0) {
        throw std::runtime_error("Screen shader linking failed");
    }
}

// Load lens flare post-process shader
void Renderer::loadLensFlareShader(const std::string& vertPath, const std::string& fragPath) {
    std::string vertSource = loadShaderSource(vertPath, PASSTHROUGH_VERTEX_SHADER);
    std::string fragSource = loadShaderSource(fragPath, PASSTHROUGH_FRAGMENT_SHADER);
    
    GLuint vertShader = compileShader(GL_VERTEX_SHADER, vertSource, "LensFlare");
    GLuint fragShader = compileShader(GL_FRAGMENT_SHADER, fragSource, "LensFlare");
    
    if (vertShader == 0 || fragShader == 0) {
        std::cerr << "[LensFlare] Shader compilation failed, using fallback" << std::endl;
        return;
    }
    
    m_LensFlareProgram = linkShaderProgram(vertShader, fragShader, "LensFlare");
}

// Load bloom post-process shader
void Renderer::loadBloomShader(const std::string& vertPath, const std::string& fragPath) {
    std::string vertSource = loadShaderSource(vertPath, PASSTHROUGH_VERTEX_SHADER);
    std::string fragSource = loadShaderSource(fragPath, PASSTHROUGH_FRAGMENT_SHADER);
    
    GLuint vertShader = compileShader(GL_VERTEX_SHADER, vertSource, "Bloom");
    GLuint fragShader = compileShader(GL_FRAGMENT_SHADER, fragSource, "Bloom");
    
    if (vertShader == 0 || fragShader == 0) {
        std::cerr << "[Bloom] Shader compilation failed, using fallback" << std::endl;
        return;
    }
    
    m_BloomProgram = linkShaderProgram(vertShader, fragShader, "Bloom");
}

bool Renderer::loadBackgroundTexture(const std::string& path) {
    std::string fullPath = resolvePath(path);
    
    // Load image using stb_image
    stbi_set_flip_vertically_on_load(true);
    int width, height, channels;
    unsigned char* data = stbi_load(fullPath.c_str(), &width, &height, &channels, 4);
    
    if (!data) {
        std::cerr << "Failed to load background texture: " << fullPath << std::endl;
        std::cerr << "stb_image error: " << stbi_failure_reason() << std::endl;
        return false;
    }
    
    m_BackgroundWidth = width;
    m_BackgroundHeight = height;
    
    // Delete existing texture if present
    if (m_BackgroundTexture != 0) {
        glDeleteTextures(1, &m_BackgroundTexture);
    }
    
    // Create OpenGL texture
    glGenTextures(1, &m_BackgroundTexture);
    glBindTexture(GL_TEXTURE_2D, m_BackgroundTexture);
    
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
    glGenerateMipmap(GL_TEXTURE_2D);
    
    stbi_image_free(data);
    
    m_UseBackgroundTexture = true;
    std::cout << "Background texture loaded: " << path << " (" << width << "x" << height << ")" << std::endl;

#ifdef SIRIUS_HAS_OPTIX
    // Upload to CUDA for OptiX use
    if (m_OptixHandle) {
        // Reload the raw data for CUDA (OpenGL texture is separate)
        stbi_set_flip_vertically_on_load(false);  // CUDA doesn't need flip
        unsigned char* cudaData = stbi_load(fullPath.c_str(), &width, &height, &channels, 4);
        if (cudaData) {
            sirius_optix_upload_background(m_OptixHandle, cudaData, width, height);
            stbi_image_free(cudaData);
        }
    }
#endif

    return true;
}

void Renderer::render(IMetric* metric, const Vec4& observerPos, const Vec4& observerVel, 
                      float cameraYaw, float cameraPitch, float cameraFOV) {
    
#ifdef SIRIUS_HAS_OPTIX
    if (m_OptixEnabled && m_OptixHandle) {
        renderOptiX(metric, observerPos, observerVel, cameraYaw, cameraPitch, cameraFOV);
        return;
    }
#endif

    // Fallback error - OptiX required
    std::cerr << "Error: OptiX backend not available. Build with SIRIUS_HAS_OPTIX=1" << std::endl;
    
    // Clear to error color
    glClearColor(0.5f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
}

void Renderer::uploadNumericalMetric(IMetric* metric) {
    // Numerical metrics removed - system uses only analytic families
    (void)metric;
}

#ifdef SIRIUS_HAS_OPTIX

// Helper: Resolve metric name to MetricType enum
static Sirius::MetricType resolveMetricType(const std::string& name) {
    if (name == "Minkowski") return Sirius::MetricType::Minkowski;
    if (name == "Schwarzschild") return Sirius::MetricType::Schwarzschild;
    if (name == "Kerr") return Sirius::MetricType::Kerr;
    if (name == "Reissner-Nordstrom") return Sirius::MetricType::ReissnerNordstrom;
    if (name == "Gödel") return Sirius::MetricType::Godel;
    if (name == "Taub-NUT") return Sirius::MetricType::TaubNUT;
    if (name == "Kerr-Schild") return Sirius::MetricType::KerrSchild;
    if (name == "Ellis-Drainhole") return Sirius::MetricType::EllisDrainhole;
    if (name == "Alcubierre") return Sirius::MetricType::Alcubierre;
    if (name == "de-Sitter") return Sirius::MetricType::DeSitter;
    return Sirius::MetricType::Minkowski;
}

// Helper: Extract metric parameters (M, a, Q) from IMetric config
static void extractMetricParams(IMetric* metric, Sirius::MetricParams& params) {
    const auto& cfg = metric->getParameters();
    auto it = cfg.find("mass");
    params.M = (it != cfg.end()) ? static_cast<float>(it->second.value) : 1.0f;
    params.rs = 2.0f * params.M;
    
    it = cfg.find("spin");
    params.a = (it != cfg.end()) ? static_cast<float>(it->second.value) : 0.9f;
    
    it = cfg.find("charge");
    params.Q = (it != cfg.end()) ? static_cast<float>(it->second.value) : 0.0f;
}

// Helper: Calculate ISCO radius for accretion disk
static float calculateISCO(float M, float a) {
    if (fabsf(a) < 1e-6f) {
        return 6.0f; // Schwarzschild ISCO
    }
    float a_star = std::max(-0.999f, std::min(0.999f, a / M));
    float Z1 = 1.0f + powf(1.0f - a_star*a_star, 1.0f/3.0f) * 
               (powf(1.0f + a_star, 1.0f/3.0f) + powf(1.0f - a_star, 1.0f/3.0f));
    float Z2 = sqrtf(3.0f * a_star*a_star + Z1*Z1);
    return 3.0f + Z2 - sqrtf((3.0f - Z1) * (3.0f + Z1 + 2.0f*Z2));
}
void Renderer::renderOptiX(IMetric* metric, const Vec4& observerPos, const Vec4& observerVel, 
                           float cameraYaw, float cameraPitch, float cameraFOV) {
    // Build launch parameters using correct structure from RDOP003A.h
    Sirius::LaunchParams params = Sirius::createDefaultLaunchParams();
    
    // Camera state - convert from spherical (r, theta, phi) to Cartesian for CUDA
    // observerPos format: (t, r, theta, phi) where indices 1,2,3 are spatial
    float r = static_cast<float>(observerPos(1));
    float theta = static_cast<float>(observerPos(2));
    float phi = static_cast<float>(observerPos(3));
    
    // Convert spherical to Cartesian for the raytracer
    float sinTheta = sinf(theta);
    float cosTheta = cosf(theta);
    float sinPhi = sinf(phi);
    float cosPhi = cosf(phi);
    
    float camX = r * sinTheta * cosPhi;
    float camY = r * cosTheta;
    float camZ = r * sinTheta * sinPhi;
    
    // =========================================================================
    // CAMERA MOVEMENT DETECTION - Reset accumulation when camera moves
    // =========================================================================
    // Without this, progressive accumulation makes viewport appear frozen
    const float posThreshold = 0.001f;
    const float rotThreshold = 0.0001f;
    
    bool cameraMoved = 
        fabsf(camX - m_LastCameraX) > posThreshold ||
        fabsf(camY - m_LastCameraY) > posThreshold ||
        fabsf(camZ - m_LastCameraZ) > posThreshold ||
        fabsf(cameraYaw - m_LastCameraYaw) > rotThreshold ||
        fabsf(cameraPitch - m_LastCameraPitch) > rotThreshold;
    
    if (cameraMoved) {
        m_FrameCount = 0;  // Reset accumulation for fresh rendering
    }
    
    // Store current camera state for next frame comparison
    m_LastCameraX = camX;
    m_LastCameraY = camY;
    m_LastCameraZ = camZ;
    m_LastCameraYaw = cameraYaw;
    m_LastCameraPitch = cameraPitch;
    
    m_FrameCount++;
    m_AnimationFrame++;  // Always increment for smooth animation
    params.frameIndex = m_FrameCount;
    
    params.camera.position = Sirius::make_float3(camX, camY, camZ);
    params.camera.velocity = Sirius::make_float3(
        static_cast<float>(observerVel(1)),
        static_cast<float>(observerVel(2)),
        static_cast<float>(observerVel(3))
    );
    params.camera.yaw = cameraYaw;
    params.camera.pitch = cameraPitch;
    params.camera.fov = cameraFOV;
    params.camera.aspectRatio = static_cast<float>(m_Width) / static_cast<float>(m_Height);
    
    // Pass observer's spherical coordinates for world-frame background sampling
    // This fixes the banding artifact by ensuring all rays are transformed to
    // a consistent world frame rather than using their escape positions
    params.camera.observerTheta = theta;
    params.camera.observerPhi = phi;
    
    // Calculate look direction from yaw/pitch
    float cosP = cosf(cameraPitch);
    float sinP = sinf(cameraPitch);
    float cosY = cosf(cameraYaw);
    float sinY = sinf(cameraYaw);
    params.camera.direction = Sirius::make_float3(cosP * sinY, sinP, cosP * cosY);
    params.camera.up = Sirius::make_float3(0.0f, 1.0f, 0.0f);
    params.camera.right = Sirius::make_float3(cosY, 0.0f, -sinY);
    
    // Metric configuration
    params.metricType = resolveMetricType(metric->getName());
    extractMetricParams(metric, params.metricParams);
    
    // Integration settings - uses IntegrationParams
    // Balanced for quality vs performance
    params.integration.type = Sirius::IntegratorType::Symplectic;
    
    // Scale maxDistance and maxSteps with observer position
    // Rays from distant observers need more steps to converge to horizon or escape
    float observerRadius = sqrtf(camX*camX + camY*camY + camZ*camZ);
    float stepsMultiplier = std::max(1.0f, observerRadius / 10.0f);
    params.integration.maxSteps = static_cast<int>(2000 * stepsMultiplier);
    params.integration.maxDistance = std::max(200.0f, observerRadius * 3.0f);
    
    params.integration.initialStepSize = 0.02f;  // Balanced for speed vs accuracy
    params.integration.minStepSize = 1e-5f;
    params.integration.maxStepSize = 0.5f;       // Allow larger steps
    params.integration.tolerance = 1e-6f;
    params.integration.horizonBuffer = 0.01f;
    params.integration.detectHorizon = true;
    
    // =========================================================================
    // KERR-SPECIFIC PARAMETERS
    // Now using finite-difference Christoffels (slower but correct).
    // Keeping moderate tolerance for good accuracy with frame-dragging.
    // =========================================================================
    if (params.metricType == Sirius::MetricType::Kerr) {
        params.integration.tolerance = 1e-6f;     // Standard tolerance
        params.integration.maxStepSize = 0.1f;    // Allow larger steps
        params.integration.initialStepSize = 0.02f;
    }
    
    // Path tracing settings
    params.pathTracing.samplesPerPixel = 1;
    params.pathTracing.maxBounces = 1;
    params.pathTracing.seed = m_FrameCount;
    params.pathTracing.seed = m_FrameCount;
    params.pathTracing.enableAccumulation = false;  // DISABLED: causes after-images
    
    // Explicit reset trigger from toggles
    if (m_ResetAccumulation) {
        params.pathTracing.resetAccumulation = true;
        m_ResetAccumulation = false; // Clear flag after use
        m_FrameCount = 0; // Also reset frame count for RNG stability
    } else {
        params.pathTracing.resetAccumulation = false;
    }
    
    // Wire up runtime toggles
    // Note: rayBundle.enabled uses default (true) from createDefaultLaunchParams()
    params.lensFlare.enabled = m_LensFlareEnabled;
    params.lensFlare.intensity = m_LensFlareIntensity;
    params.lensFlare.threshold = m_LensFlareThreshold;
    
    // Wire up disk mode (Phase 7.5 - Planar vs Volumetric)
    params.accretionDisk.diskMode = m_DiskModePlanar 
        ? Sirius::DiskMode::Planar 
        : Sirius::DiskMode::Volumetric;
    
    // Animation time - use m_AnimationFrame for continuous animation (NEVER resets)
    // This drives disk ring animation independent of camera movement/accumulation resets
    params.time = static_cast<float>(m_AnimationFrame) * 0.033f;  // ~30fps = 0.033s per frame
    
    // Background
    params.useBackgroundTexture = m_UseBackgroundTexture;
    if (m_UseBackgroundTexture && m_OptixHandle) {
        params.backgroundTexture = sirius_optix_get_background_texture(m_OptixHandle);
    }
    // params.backgroundColor already set to (0.02, 0.02, 0.08) in defaults
    
    // NOTE: frameIndex already set above, don't increment again
    
    // Auto-calculate ISCO for accretion disk
    if (params.accretionDisk.innerRadius <= 0.001f) {
        params.accretionDisk.innerRadius = calculateISCO(params.metricParams.M, params.metricParams.a);
        params.accretionDisk.heightScale = std::max(params.accretionDisk.heightScale, 0.005f);
    }
    
    // Launch OptiX
    sirius_optix_set_metric_type(m_OptixHandle, static_cast<int>(params.metricType));
    sirius_optix_launch(m_OptixHandle, &params);
    
    // Apply AI denoiser for noise reduction (Phase 4 Optimization)
    // blendFactor: 0.0 = fully denoised, 1.0 = original noisy
    sirius_optix_denoise(m_OptixHandle, 0.0f);
    
    // Copy result to OpenGL texture (Hardware Interop)
    sirius_optix_update_display(m_OptixHandle);
    
    // Draw to screen with post-processing (Phase 6.5 Lens Flare, Phase 7 Bloom)
    glViewport(0, 0, m_Width, m_Height);
    glClear(GL_COLOR_BUFFER_BIT);
    
    // Select shader: bloom > lens flare > basic screen
    // Bloom is the primary cinematic effect (enabled by default)
    GLuint activeProgram;
    if (m_BloomEnabled && m_BloomProgram != 0) {
        activeProgram = m_BloomProgram;
    } else if (m_LensFlareEnabled && m_LensFlareProgram != 0) {
        activeProgram = m_LensFlareProgram;
    } else {
        activeProgram = m_ScreenProgram;
    }
    
    glUseProgram(activeProgram);
    
    // Set texture uniform
    GLint texLoc = glGetUniformLocation(activeProgram, "screenTexture");
    if (texLoc != -1) {
        glUniform1i(texLoc, 0);
    }
    
    // Set bloom uniforms if using bloom shader
    if (m_BloomEnabled && m_BloomProgram != 0 && activeProgram == m_BloomProgram) {
        GLint enableLoc = glGetUniformLocation(activeProgram, "bloomEnabled");
        if (enableLoc != -1) glUniform1i(enableLoc, 1);
        
        GLint intensityLoc = glGetUniformLocation(activeProgram, "bloomIntensity");
        if (intensityLoc != -1) glUniform1f(intensityLoc, m_BloomIntensity);
        
        GLint thresholdLoc = glGetUniformLocation(activeProgram, "bloomThreshold");
        if (thresholdLoc != -1) glUniform1f(thresholdLoc, m_BloomThreshold);
        
        GLint texelSizeLoc = glGetUniformLocation(activeProgram, "texelSize");
        if (texelSizeLoc != -1) {
            glUniform2f(texelSizeLoc, 1.0f / m_Width, 1.0f / m_Height);
        }
    }
    // Set lens flare uniforms if using lens flare shader
    else if (m_LensFlareEnabled && m_LensFlareProgram != 0 && activeProgram == m_LensFlareProgram) {
        GLint enableLoc = glGetUniformLocation(activeProgram, "enableLensFlare");
        if (enableLoc != -1) glUniform1i(enableLoc, 1);
        
        GLint intensityLoc = glGetUniformLocation(activeProgram, "flareIntensity");
        if (intensityLoc != -1) glUniform1f(intensityLoc, m_LensFlareIntensity);
        
        GLint thresholdLoc = glGetUniformLocation(activeProgram, "flareThreshold");
        if (thresholdLoc != -1) glUniform1f(thresholdLoc, m_LensFlareThreshold);
        
        GLint texelSizeLoc = glGetUniformLocation(activeProgram, "texelSize");
        if (texelSizeLoc != -1) {
            glUniform2f(texelSizeLoc, 1.0f / m_Width, 1.0f / m_Height);
        }
    }
    
    glBindVertexArray(m_Vao);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, m_Texture);
    glDrawArrays(GL_TRIANGLES, 0, 6);
}
#endif

void Renderer::setupMetricUniforms(IMetric* /*metric*/, const Vec4& /*position*/) {
    // Legacy - no longer used with OptiX backend
}

void Renderer::loadComputeShader(const std::string& /*path*/) {
    // Legacy - no longer used with OptiX backend
}

void Renderer::setMetricParameters(const std::map<std::string, double>& /*params*/) {
    // Reset frame counter to restart accumulation
    m_FrameCount = 0;
}

void Renderer::cleanup() {
#ifdef SIRIUS_HAS_OPTIX
    if (m_OptixHandle) {
        sirius_optix_cleanup(m_OptixHandle);
        sirius_optix_destroy(m_OptixHandle);
        m_OptixHandle = nullptr;
        m_OptixEnabled = false;
    }
    
    // Free Pinned Memory using CUDA allocator
    if (m_PinnedMemoryBuffer) {
        cudaFreeHost(m_PinnedMemoryBuffer);
        m_PinnedMemoryBuffer = nullptr;
        m_PinnedMemoryCapacity = 0;
    }
#else
    // Fallback cleanup if compiled without OptiX (but with this code enabled)
    if (m_PinnedMemoryBuffer) {
        delete[] m_PinnedMemoryBuffer;
        m_PinnedMemoryBuffer = nullptr;
        m_PinnedMemoryCapacity = 0;
    }
#endif
    
    if (m_MetricUBO != 0) {
        glDeleteBuffers(1, &m_MetricUBO);
        m_MetricUBO = 0;
    }
    
    if (m_ComputeProgram != 0) {
        glDeleteProgram(m_ComputeProgram);
        m_ComputeProgram = 0;
    }
    
    if (m_ScreenProgram != 0) {
        glDeleteProgram(m_ScreenProgram);
        m_ScreenProgram = 0;
    }
    
    // Cleanup lens flare shader (Phase 6.5)
    if (m_LensFlareProgram != 0) {
        glDeleteProgram(m_LensFlareProgram);
        m_LensFlareProgram = 0;
    }
    
    // Cleanup bloom shader (Phase 7)
    if (m_BloomProgram != 0) {
        glDeleteProgram(m_BloomProgram);
        m_BloomProgram = 0;
    }
    
    if (m_Texture != 0) {
        glDeleteTextures(1, &m_Texture);
        m_Texture = 0;
    }
    
    if (m_BackgroundTexture != 0) {
        glDeleteTextures(1, &m_BackgroundTexture);
        m_BackgroundTexture = 0;
    }
    
    if (m_Vao != 0) {
        glDeleteVertexArrays(1, &m_Vao);
        m_Vao = 0;
    }
}