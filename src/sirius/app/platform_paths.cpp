// Platform path resolution. Ported from CRPF001A.cpp.

#include "sirius/app/platform_paths.h"

#include "sirius/base/resource_locator.h"

#include <cstdlib>
#include <fstream>

namespace sirius::app {

fs::path PlatformPaths::ExecutableDirectory() { return base::ExecutableDirectory(); }

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
    return base::ResolveResource(relative_path);
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
