// Platform path resolution. Ported from CRPF001A.cpp.

#include "sirius/app/platform_paths.h"

#include <cstdlib>
#include <fstream>

#if defined(_WIN32)
#include <windows.h>
#elif defined(__linux__)
#include <linux/limits.h>
#include <unistd.h>
#elif defined(__APPLE__)
#include <limits.h>
#include <mach-o/dyld.h>
#endif

namespace sirius::app {

fs::path PlatformPaths::executable_path_;
bool PlatformPaths::executable_path_initialised_ = false;

fs::path PlatformPaths::ExecutableDirectory() {
    if (!executable_path_initialised_) {
#if defined(_WIN32)
        wchar_t path[MAX_PATH];
        DWORD len = GetModuleFileNameW(NULL, path, MAX_PATH);
        if (len > 0 && len < MAX_PATH) {
            executable_path_ = fs::path(path).parent_path();
        }
#elif defined(__linux__)
        char path[PATH_MAX];
        ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
        if (len != -1) {
            path[len] = '\0';
            executable_path_ = fs::path(path).parent_path();
        }
#elif defined(__APPLE__)
        char path[PATH_MAX];
        uint32_t size = sizeof(path);
        if (_NSGetExecutablePath(path, &size) == 0) {
            executable_path_ = fs::canonical(fs::path(path)).parent_path();
        }
#endif
        executable_path_initialised_ = true;
    }
    return executable_path_;
}

fs::path PlatformPaths::UserConfigDirectory() {
    fs::path config_dir;

#if defined(_WIN32)
    const char* appdata = std::getenv("APPDATA");
    if (appdata) {
        config_dir = fs::path(appdata) / "Sirius";
    }
#elif defined(__APPLE__)
    const char* home = std::getenv("HOME");
    if (home) {
        config_dir = fs::path(home) / "Library" / "Application Support" / "Sirius";
    }
#else  // Linux and other Unix.
    const char* xdg_config = std::getenv("XDG_CONFIG_HOME");
    if (xdg_config && xdg_config[0] != '\0') {
        config_dir = fs::path(xdg_config) / "sirius";
    } else {
        const char* home = std::getenv("HOME");
        if (home) {
            config_dir = fs::path(home) / ".config" / "sirius";
        }
    }
#endif

    if (!config_dir.empty() && !fs::exists(config_dir)) {
        std::error_code ec;
        fs::create_directories(config_dir, ec);
    }

    return config_dir;
}

fs::path PlatformPaths::SystemConfigDirectory() {
#if defined(_WIN32)
    const char* program_data = std::getenv("PROGRAMDATA");
    if (program_data) {
        return fs::path(program_data) / "Sirius";
    }
    return fs::path("C:/ProgramData/Sirius");
#else
    return fs::path("/etc/sirius");
#endif
}

std::optional<fs::path> PlatformPaths::ResolveResource(const std::string& relative_path) {
    // Search order, highest priority first. The shaders/ subdir is the new-tree
    // home the app CMakeLists copies viewer shaders into; the Sirius.Render entry
    // preserves the legacy runtime-load path.
    std::vector<fs::path> search_paths = {
        fs::current_path() / relative_path,
        ExecutableDirectory() / relative_path,
        ExecutableDirectory() / "shaders" / relative_path,
        ExecutableDirectory() / "Sirius.Render" / relative_path,
        ExecutableDirectory().parent_path() / relative_path,
        ExecutableDirectory().parent_path().parent_path() / relative_path,
        ExecutableDirectory().parent_path().parent_path().parent_path() / "src" / relative_path,
    };

    for (const auto& path : search_paths) {
        if (fs::exists(path)) {
            return fs::canonical(path);
        }
    }

    return std::nullopt;
}

std::vector<fs::path> PlatformPaths::ConfigSearchPaths() {
    std::vector<fs::path> paths;

    paths.push_back(fs::current_path() / "sirius.json");
    paths.push_back(ExecutableDirectory() / "sirius.json");

    fs::path user_config = UserConfigDirectory();
    if (!user_config.empty()) {
        paths.push_back(user_config / "config.json");
    }

    paths.push_back(SystemConfigDirectory() / "config.json");

    return paths;
}

std::optional<fs::path> PlatformPaths::FindConfigFile() {
    for (const auto& path : ConfigSearchPaths()) {
        if (fs::exists(path)) {
            return path;
        }
    }
    return std::nullopt;
}

std::string PlatformPaths::PlatformName() {
#if defined(_WIN32)
    return "Windows";
#elif defined(__APPLE__)
    return "macOS";
#else
    return "Linux";
#endif
}

bool PlatformPaths::IsWsl2() {
#if defined(__linux__)
    if (fs::exists("/proc/sys/fs/binfmt_misc/WSLInterop") || fs::exists("/run/WSL")) {
        return true;
    }

    std::ifstream version_file("/proc/version");
    if (version_file) {
        std::string line;
        std::getline(version_file, line);
        if (line.find("microsoft") != std::string::npos || line.find("WSL") != std::string::npos) {
            return true;
        }
    }
#endif
    return false;
}

}  // namespace sirius::app
