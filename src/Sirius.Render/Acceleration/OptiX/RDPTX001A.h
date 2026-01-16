// RDPTX001A.h - PTX Pre-compilation Utilities
// Component ID: RDPTX001A
// Multi-architecture PTX selection and path resolution
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_RDPTX001A_H
#define SIRIUS_RENDER_RDPTX001A_H

#include <cuda_runtime.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <limits.h>
#endif

namespace sirius::render {

//==============================================================================
// PTX Architecture Info
//==============================================================================

struct PTXArchInfo {
    int smMajor;           // SM major version
    int smMinor;           // SM minor version
    std::string archName;  // e.g., "sm_75", "sm_86"
    std::string ptxSuffix; // e.g., "_sm75", "_sm86"
};

// Supported GPU architectures for PTX pre-compilation
inline std::vector<PTXArchInfo> getSupportedArchitectures() {
    return {
        {7, 5, "sm_75", "_sm75"},  // Turing (RTX 20xx)
        {8, 0, "sm_80", "_sm80"},  // Ampere GA100
        {8, 6, "sm_86", "_sm86"},  // Ampere (RTX 30xx)
        {8, 9, "sm_89", "_sm89"},  // Ada Lovelace (RTX 40xx)
    };
}

//==============================================================================
// PTX File Discovery
//==============================================================================

/// Get the directory containing the executable
inline std::string getExecutableDirectory() {
#ifdef _WIN32
    char path[MAX_PATH];
    GetModuleFileNameA(nullptr, path, MAX_PATH);
    std::string fullPath(path);
    size_t pos = fullPath.find_last_of("\\/");
    return pos != std::string::npos ? fullPath.substr(0, pos) : ".";
#else
    char path[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
    if (len != -1) {
        path[len] = '\0';
        std::string fullPath(path);
        size_t pos = fullPath.find_last_of('/');
        return pos != std::string::npos ? fullPath.substr(0, pos) : ".";
    }
    return ".";
#endif
}

/// Get list of search paths for PTX files
inline std::vector<std::string> getPTXSearchPaths(const std::string& ptxDir = "") {
    std::vector<std::string> paths;
    std::string exeDir = getExecutableDirectory();

    // 1. Explicit PTX directory (from CMake define)
    if (!ptxDir.empty()) {
        paths.push_back(ptxDir);
    }

    // 2. Relative to current working directory
    paths.push_back("Sirius.Render/ptx");
    paths.push_back("./Sirius.Render/ptx");
    paths.push_back("../Sirius.Render/ptx");

    // 3. Relative to executable
    paths.push_back(exeDir + "/Sirius.Render/ptx");
    paths.push_back(exeDir + "/ptx");
    paths.push_back(exeDir + "/../Sirius.Render/ptx");

    // 4. System installation paths
#ifndef _WIN32
    paths.push_back("/usr/local/lib/sirius/ptx");
    paths.push_back("/usr/share/sirius/ptx");
#endif

    return paths;
}

/// Check if a file exists
inline bool fileExists(const std::string& path) {
    std::ifstream f(path);
    return f.good();
}

//==============================================================================
// Architecture Selection
//==============================================================================

/// Get the compute capability of the current device
inline bool getDeviceComputeCapability(int deviceId, int& smMajor, int& smMinor) {
    cudaDeviceProp props;
    cudaError_t err = cudaGetDeviceProperties(&props, deviceId);
    if (err != cudaSuccess) {
        return false;
    }
    smMajor = props.major;
    smMinor = props.minor;
    return true;
}

/// Find the best matching architecture for a device
/// Returns the highest compatible architecture that doesn't exceed device capability
inline PTXArchInfo selectBestArchitecture(int smMajor, int smMinor) {
    auto archs = getSupportedArchitectures();
    int deviceSM = smMajor * 10 + smMinor;

    PTXArchInfo bestMatch = archs[0];  // Default to lowest supported
    int bestSM = 0;

    for (const auto& arch : archs) {
        int archSM = arch.smMajor * 10 + arch.smMinor;
        // Select the highest architecture that doesn't exceed device capability
        if (archSM <= deviceSM && archSM > bestSM) {
            bestSM = archSM;
            bestMatch = arch;
        }
    }

    return bestMatch;
}

//==============================================================================
// PTX File Resolution
//==============================================================================

/// Find the appropriate PTX file for the current device
/// @param baseName Base name of the PTX file (e.g., "RDOP002A")
/// @param ptxDir Optional directory to search first
/// @return Full path to PTX file, or empty string if not found
inline std::string findPTXFile(const std::string& baseName,
                                const std::string& ptxDir = "",
                                int deviceId = 0) {
    // Get device compute capability
    int smMajor = 7, smMinor = 5;  // Default fallback
    getDeviceComputeCapability(deviceId, smMajor, smMinor);

    // Select best matching architecture
    PTXArchInfo arch = selectBestArchitecture(smMajor, smMinor);

    // Generate candidate filenames in order of preference
    std::vector<std::string> candidates = {
        baseName + arch.ptxSuffix + ".ptx",  // Architecture-specific (preferred)
        baseName + ".ptx"                      // Generic fallback
    };

    // Search paths
    auto searchPaths = getPTXSearchPaths(ptxDir);

    // Try each candidate in each search path
    for (const auto& candidate : candidates) {
        for (const auto& path : searchPaths) {
            std::string fullPath = path + "/" + candidate;
            if (fileExists(fullPath)) {
                return fullPath;
            }
        }
    }

    // Last resort: check if generic PTX exists in any search path
    for (const auto& path : searchPaths) {
        std::string fullPath = path + "/" + baseName + ".ptx";
        if (fileExists(fullPath)) {
            return fullPath;
        }
    }

    return "";  // Not found
}

/// Load PTX file contents
inline std::string loadPTXContents(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.good()) {
        return "";
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

//==============================================================================
// PTX Version Metadata
//==============================================================================

struct PTXMetadata {
    std::string filename;
    std::string archName;
    int smMajor = 0;
    int smMinor = 0;
    bool isValid = false;
};

/// Parse PTX file header to extract version info
inline PTXMetadata parsePTXMetadata(const std::string& ptxContent) {
    PTXMetadata meta;

    // PTX files start with: .version X.Y
    // and contain: .target sm_XX
    size_t versionPos = ptxContent.find(".version");
    size_t targetPos = ptxContent.find(".target");

    if (versionPos != std::string::npos && targetPos != std::string::npos) {
        // Extract target architecture
        size_t smPos = ptxContent.find("sm_", targetPos);
        if (smPos != std::string::npos) {
            // Parse sm_XX
            size_t numStart = smPos + 3;
            size_t numEnd = ptxContent.find_first_not_of("0123456789", numStart);
            if (numEnd == std::string::npos) numEnd = ptxContent.length();

            std::string smStr = ptxContent.substr(numStart, numEnd - numStart);
            if (smStr.length() >= 2) {
                meta.smMajor = smStr[0] - '0';
                meta.smMinor = std::stoi(smStr.substr(1));
                meta.archName = "sm_" + smStr;
                meta.isValid = true;
            }
        }
    }

    return meta;
}

//==============================================================================
// PTX Compilation Status
//==============================================================================

/// Check which PTX architectures are available
inline std::vector<PTXArchInfo> getAvailablePTXArchitectures(
    const std::string& baseName,
    const std::string& ptxDir = "") {

    std::vector<PTXArchInfo> available;
    auto allArchs = getSupportedArchitectures();
    auto searchPaths = getPTXSearchPaths(ptxDir);

    for (const auto& arch : allArchs) {
        std::string filename = baseName + arch.ptxSuffix + ".ptx";
        for (const auto& path : searchPaths) {
            if (fileExists(path + "/" + filename)) {
                available.push_back(arch);
                break;
            }
        }
    }

    // Also check for generic PTX
    for (const auto& path : searchPaths) {
        if (fileExists(path + "/" + baseName + ".ptx")) {
            // Generic PTX available as fallback
            break;
        }
    }

    return available;
}

} // namespace sirius::render

#endif // SIRIUS_RENDER_RDPTX001A_H
