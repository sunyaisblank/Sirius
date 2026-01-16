// CRPF001A.cpp - Platform Path Resolution
// Component ID: CRPF001A (Platform/Paths)

#include "CRPF001A.h"

#include <cstdlib>
#include <fstream>

#if defined(_WIN32)
    #include <windows.h>
#elif defined(__linux__)
    #include <unistd.h>
    #include <linux/limits.h>
#elif defined(__APPLE__)
    #include <mach-o/dyld.h>
    #include <limits.h>
#endif

namespace Sirius::Platform {

// Static member initialization
fs::path PathResolver::s_executablePath;
bool PathResolver::s_executablePathInitialized = false;

fs::path PathResolver::executableDirectory() {
    if (!s_executablePathInitialized) {
        #if defined(_WIN32)
            wchar_t path[MAX_PATH];
            DWORD len = GetModuleFileNameW(NULL, path, MAX_PATH);
            if (len > 0 && len < MAX_PATH) {
                s_executablePath = fs::path(path).parent_path();
            }
        #elif defined(__linux__)
            char path[PATH_MAX];
            ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
            if (len != -1) {
                path[len] = '\0';
                s_executablePath = fs::path(path).parent_path();
            }
        #elif defined(__APPLE__)
            char path[PATH_MAX];
            uint32_t size = sizeof(path);
            if (_NSGetExecutablePath(path, &size) == 0) {
                s_executablePath = fs::canonical(fs::path(path)).parent_path();
            }
        #endif
        s_executablePathInitialized = true;
    }
    return s_executablePath;
}

fs::path PathResolver::userConfigDirectory() {
    fs::path configDir;

    #if defined(_WIN32)
        const char* appdata = std::getenv("APPDATA");
        if (appdata) {
            configDir = fs::path(appdata) / "Sirius";
        }
    #elif defined(__APPLE__)
        const char* home = std::getenv("HOME");
        if (home) {
            configDir = fs::path(home) / "Library" / "Application Support" / "Sirius";
        }
    #else // Linux and other Unix
        // Check XDG_CONFIG_HOME first
        const char* xdgConfig = std::getenv("XDG_CONFIG_HOME");
        if (xdgConfig && xdgConfig[0] != '\0') {
            configDir = fs::path(xdgConfig) / "sirius";
        } else {
            const char* home = std::getenv("HOME");
            if (home) {
                configDir = fs::path(home) / ".config" / "sirius";
            }
        }
    #endif

    // Create directory if it doesn't exist
    if (!configDir.empty() && !fs::exists(configDir)) {
        std::error_code ec;
        fs::create_directories(configDir, ec);
    }

    return configDir;
}

fs::path PathResolver::systemConfigDirectory() {
    #if defined(_WIN32)
        const char* programData = std::getenv("PROGRAMDATA");
        if (programData) {
            return fs::path(programData) / "Sirius";
        }
        return fs::path("C:/ProgramData/Sirius");
    #else
        return fs::path("/etc/sirius");
    #endif
}

std::optional<fs::path> PathResolver::resolveResource(const std::string& relativePath) {
    // Search paths in priority order
    std::vector<fs::path> searchPaths = {
        // 1. Current working directory
        fs::current_path() / relativePath,
        // 2. Executable directory
        executableDirectory() / relativePath,
        // 3. Sirius.Render subdirectory (for shaders/textures)
        executableDirectory() / "Sirius.Render" / relativePath,
        // 4. One level up (for VS multi-config builds)
        executableDirectory().parent_path() / relativePath,
        // 5. Two levels up
        executableDirectory().parent_path().parent_path() / relativePath,
        // 6. Source tree location (development)
        executableDirectory().parent_path().parent_path().parent_path() / "src" / relativePath
    };

    for (const auto& path : searchPaths) {
        if (fs::exists(path)) {
            return fs::canonical(path);
        }
    }

    return std::nullopt;
}

std::vector<fs::path> PathResolver::configSearchPaths() {
    std::vector<fs::path> paths;

    // 1. Current directory (highest priority)
    paths.push_back(fs::current_path() / "sirius.json");

    // 2. Executable directory
    paths.push_back(executableDirectory() / "sirius.json");

    // 3. User config directory
    fs::path userConfig = userConfigDirectory();
    if (!userConfig.empty()) {
        paths.push_back(userConfig / "config.json");
    }

    // 4. System config directory
    paths.push_back(systemConfigDirectory() / "config.json");

    return paths;
}

std::optional<fs::path> PathResolver::findConfigFile() {
    for (const auto& path : configSearchPaths()) {
        if (fs::exists(path)) {
            return path;
        }
    }
    return std::nullopt;
}

std::string PathResolver::platformName() {
    #if defined(_WIN32)
        return "Windows";
    #elif defined(__APPLE__)
        return "macOS";
    #else
        return "Linux";
    #endif
}

bool PathResolver::isWSL2() {
    #if defined(__linux__)
        // Check for WSL2 indicators
        if (fs::exists("/proc/sys/fs/binfmt_misc/WSLInterop") ||
            fs::exists("/run/WSL")) {
            return true;
        }

        // Check /proc/version for "microsoft" or "WSL"
        std::ifstream versionFile("/proc/version");
        if (versionFile) {
            std::string line;
            std::getline(versionFile, line);
            if (line.find("microsoft") != std::string::npos ||
                line.find("WSL") != std::string::npos) {
                return true;
            }
        }
    #endif
    return false;
}

} // namespace Sirius::Platform
