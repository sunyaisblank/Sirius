// SRHL001A.h - Headless Rendering Mode
// Component ID: SRHL001A
// Environment detection and headless rendering utilities
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_SRHL001A_H
#define SIRIUS_RENDER_SRHL001A_H

#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace sirius::render {

//==============================================================================
// Headless Environment Detection
//==============================================================================

/// Check if running in a headless environment (no display)
inline bool isHeadlessEnvironment() {
    // 1. Check explicit HEADLESS environment variable
    const char* headless = std::getenv("HEADLESS");
    if (headless && (std::string(headless) == "1" ||
                     std::string(headless) == "true" ||
                     std::string(headless) == "yes")) {
        return true;
    }

#ifdef _WIN32
    // Windows: Assume display is always available unless HEADLESS is set
    return false;
#else
    // 2. Check for DISPLAY environment variable (Linux/Unix)
    const char* display = std::getenv("DISPLAY");
    if (!display || display[0] == '\0') {
        return true;
    }

    // 3. Check for container environments
    if (isRunningInContainer()) {
        // In containers, check if display is actually accessible
        // Even if DISPLAY is set, X server may not be available
        return !isDisplayAccessible();
    }

    return false;
#endif
}

/// Check if running inside a container (Docker, Singularity, etc.)
inline bool isRunningInContainer() {
#ifdef _WIN32
    return false;
#else
    // Check for Docker
    std::ifstream dockerenv("/.dockerenv");
    if (dockerenv.good()) return true;

    // Check for container in cgroups (works for Docker, Kubernetes, etc.)
    std::ifstream cgroup("/proc/1/cgroup");
    if (cgroup.good()) {
        std::string line;
        while (std::getline(cgroup, line)) {
            if (line.find("docker") != std::string::npos ||
                line.find("kubepods") != std::string::npos ||
                line.find("containerd") != std::string::npos) {
                return true;
            }
        }
    }

    // Check for Singularity
    if (std::getenv("SINGULARITY_CONTAINER") != nullptr) {
        return true;
    }

    return false;
#endif
}

/// Check if X display is actually accessible (Linux)
inline bool isDisplayAccessible() {
#ifdef _WIN32
    return true;
#else
    const char* display = std::getenv("DISPLAY");
    if (!display) return false;

    // Try to connect to X server via a simple test
    // This is a lightweight check without requiring X11 headers
    std::string displayStr(display);

    // Local display (:0, :1, etc.)
    if (displayStr[0] == ':') {
        // Check if X socket exists
        std::string socketPath = "/tmp/.X11-unix/X" + displayStr.substr(1, 1);
        std::ifstream socket(socketPath);
        return socket.good();
    }

    // Remote display (hostname:0) - assume accessible if DISPLAY is set
    return true;
#endif
}

//==============================================================================
// Headless Mode Configuration
//==============================================================================

struct HeadlessConfig {
    bool forceHeadless = false;      // Override automatic detection
    bool useEGL = true;              // Use EGL for offscreen context (Linux)
    bool useOSMesa = false;          // Use software rendering fallback
    int virtualDisplayWidth = 1920;
    int virtualDisplayHeight = 1080;
    std::string outputPath;          // Default output directory
    std::string outputFormat = "exr"; // Default: exr, ppm, png
};

//==============================================================================
// Headless Renderer Context
// Provides GPU rendering without display window
//==============================================================================

class HeadlessContext {
public:
    HeadlessContext() = default;
    ~HeadlessContext() { shutdown(); }

    /// Initialize headless rendering context
    /// @param config Headless configuration
    /// @return true if initialization succeeded
    bool initialize(const HeadlessConfig& config = HeadlessConfig()) {
        m_config = config;

        if (m_config.forceHeadless || isHeadlessEnvironment()) {
            m_isHeadless = true;
        }

        // For CUDA/OptiX, we don't need a display context
        // Just ensure CUDA is initialized
#ifdef SIRIUS_HAS_CUDA
        int deviceCount = 0;
        cudaError_t err = cudaGetDeviceCount(&deviceCount);
        if (err != cudaSuccess || deviceCount == 0) {
            m_lastError = "No CUDA devices available";
            return false;
        }

        // Select best device
        cudaDeviceProp bestProps;
        int bestDevice = 0;
        int bestScore = 0;

        for (int i = 0; i < deviceCount; ++i) {
            cudaDeviceProp props;
            cudaGetDeviceProperties(&props, i);
            int score = props.multiProcessorCount * props.clockRate;
            if (score > bestScore) {
                bestScore = score;
                bestDevice = i;
                bestProps = props;
            }
        }

        cudaSetDevice(bestDevice);
        m_deviceId = bestDevice;
        m_deviceName = bestProps.name;
#endif

        m_initialized = true;
        return true;
    }

    /// Shutdown context
    void shutdown() {
        m_initialized = false;
    }

    /// Check if running in headless mode
    bool isHeadless() const { return m_isHeadless; }

    /// Check if context is initialized
    bool isInitialized() const { return m_initialized; }

    /// Get device name
    const std::string& deviceName() const { return m_deviceName; }

    /// Get device ID
    int deviceId() const { return m_deviceId; }

    /// Get last error message
    const std::string& lastError() const { return m_lastError; }

    /// Get configuration
    const HeadlessConfig& config() const { return m_config; }

private:
    HeadlessConfig m_config;
    bool m_isHeadless = false;
    bool m_initialized = false;
    int m_deviceId = 0;
    std::string m_deviceName;
    std::string m_lastError;
};

//==============================================================================
// Framebuffer Capture Utilities
//==============================================================================

/// Write framebuffer to PPM file (simple, always available)
inline bool writeFramebufferPPM(const std::string& filename,
                                 const float* data,
                                 int width, int height) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) return false;

    // PPM header
    file << "P6\n" << width << " " << height << "\n255\n";

    // Convert float RGBA to byte RGB
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int idx = (y * width + x) * 4;  // Assuming RGBA float
            float r = data[idx + 0];
            float g = data[idx + 1];
            float b = data[idx + 2];

            // Clamp and convert to 8-bit
            auto toByte = [](float v) -> unsigned char {
                v = std::max(0.0f, std::min(1.0f, v));
                return static_cast<unsigned char>(v * 255.0f + 0.5f);
            };

            file.put(toByte(r));
            file.put(toByte(g));
            file.put(toByte(b));
        }
    }

    return file.good();
}

/// Write framebuffer to PFM file (Portable Float Map - HDR)
inline bool writeFramebufferPFM(const std::string& filename,
                                 const float* data,
                                 int width, int height) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) return false;

    // PFM header (color, little-endian)
    file << "PF\n" << width << " " << height << "\n-1.0\n";

    // Write RGB float data (PFM is bottom-up)
    for (int y = height - 1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            int idx = (y * width + x) * 4;  // Assuming RGBA float
            file.write(reinterpret_cast<const char*>(&data[idx + 0]), sizeof(float));
            file.write(reinterpret_cast<const char*>(&data[idx + 1]), sizeof(float));
            file.write(reinterpret_cast<const char*>(&data[idx + 2]), sizeof(float));
        }
    }

    return file.good();
}

//==============================================================================
// Render Output Configuration
//==============================================================================

struct RenderOutputConfig {
    std::string basePath = "render";
    std::string format = "exr";     // exr, pfm, ppm, png
    bool includeTimestamp = true;
    bool includeFrameNumber = true;
    int frameDigits = 4;            // Zero-padding for frame numbers
};

/// Generate output filename
inline std::string generateOutputFilename(const RenderOutputConfig& config,
                                           int frame = 0) {
    std::string filename = config.basePath;

    if (config.includeTimestamp) {
        time_t now = time(nullptr);
        char timebuf[32];
        strftime(timebuf, sizeof(timebuf), "_%Y%m%d_%H%M%S", localtime(&now));
        filename += timebuf;
    }

    if (config.includeFrameNumber) {
        char framebuf[16];
        snprintf(framebuf, sizeof(framebuf), "_%0*d", config.frameDigits, frame);
        filename += framebuf;
    }

    filename += "." + config.format;
    return filename;
}

//==============================================================================
// Environment Information
//==============================================================================

struct EnvironmentInfo {
    bool isHeadless = false;
    bool isContainer = false;
    bool hasDisplay = false;
    bool hasCUDA = false;
    int cudaDeviceCount = 0;
    std::string hostname;
    std::string displayVar;
};

/// Gather environment information
inline EnvironmentInfo getEnvironmentInfo() {
    EnvironmentInfo info;

    info.isHeadless = isHeadlessEnvironment();
    info.isContainer = isRunningInContainer();
    info.hasDisplay = isDisplayAccessible();

    const char* display = std::getenv("DISPLAY");
    info.displayVar = display ? display : "(not set)";

#ifdef _WIN32
    char hostname[256];
    DWORD size = sizeof(hostname);
    GetComputerNameA(hostname, &size);
    info.hostname = hostname;
#else
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    info.hostname = hostname;
#endif

#ifdef SIRIUS_HAS_CUDA
    int count = 0;
    if (cudaGetDeviceCount(&count) == cudaSuccess) {
        info.hasCUDA = true;
        info.cudaDeviceCount = count;
    }
#endif

    return info;
}

/// Print environment information to stdout
inline void printEnvironmentInfo() {
    auto info = getEnvironmentInfo();

    printf("=== Sirius Render Environment ===\n");
    printf("Hostname:       %s\n", info.hostname.c_str());
    printf("DISPLAY:        %s\n", info.displayVar.c_str());
    printf("Container:      %s\n", info.isContainer ? "yes" : "no");
    printf("Display OK:     %s\n", info.hasDisplay ? "yes" : "no");
    printf("Headless Mode:  %s\n", info.isHeadless ? "yes" : "no");
    printf("CUDA Available: %s\n", info.hasCUDA ? "yes" : "no");
    if (info.hasCUDA) {
        printf("CUDA Devices:   %d\n", info.cudaDeviceCount);
    }
    printf("=================================\n");
}

} // namespace sirius::render

#endif // SIRIUS_RENDER_SRHL001A_H
