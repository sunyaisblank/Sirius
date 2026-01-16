// CRPF001A.h - Platform Path Resolution
// Component ID: CRPF001A (Platform/Paths)
//
// Cross-platform path resolution using std::filesystem.
// Provides unified access to executable, config, and resource paths.

#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace Sirius::Platform {

namespace fs = std::filesystem;

/// @brief Cross-platform path resolution utilities
class PathResolver {
public:
    /// @brief Get the directory containing the executable
    /// @return Absolute path to executable directory
    static fs::path executableDirectory();

    /// @brief Get user configuration directory
    /// Linux: ~/.config/sirius/
    /// Windows: %APPDATA%/Sirius/
    /// macOS: ~/Library/Application Support/Sirius/
    /// @return Path to user config directory (creates if needed)
    static fs::path userConfigDirectory();

    /// @brief Get system configuration directory
    /// Linux: /etc/sirius/
    /// Windows: %PROGRAMDATA%/Sirius/
    /// @return Path to system config directory
    static fs::path systemConfigDirectory();

    /// @brief Resolve a resource path (shader, texture, PTX)
    /// Searches in order: executable dir, resource subdir, install prefix
    /// @param relativePath Resource path relative to search roots
    /// @return Resolved absolute path if found, nullopt otherwise
    static std::optional<fs::path> resolveResource(const std::string& relativePath);

    /// @brief Get all potential config file locations (in priority order)
    /// @return Vector of paths to search for config files
    static std::vector<fs::path> configSearchPaths();

    /// @brief Get the first existing config file from search paths
    /// @return Path to found config file, nullopt if none exists
    static std::optional<fs::path> findConfigFile();

    /// @brief Get platform name string
    /// @return "Linux", "Windows", or "macOS"
    static std::string platformName();

    /// @brief Check if running in WSL2 environment
    /// @return true if running under WSL2
    static bool isWSL2();

private:
    /// @brief Cached executable path (computed once)
    static fs::path s_executablePath;
    static bool s_executablePathInitialized;
};

} // namespace Sirius::Platform
