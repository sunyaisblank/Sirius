#pragma once

// Cross-platform path resolution over std::filesystem: executable, config, and
// resource locations. Ported from CRPF001A.h; the class and methods take the
// new-tree PascalCase spelling, behaviour preserved.

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace sirius::app {

namespace fs = std::filesystem;

// Resolves the platform-dependent locations Sirius reads and writes.
class PlatformPaths {
  public:
    // Directory containing the running executable (absolute).
    static fs::path ExecutableDirectory();

    // User configuration directory (created if absent):
    //   Linux:   $XDG_CONFIG_HOME/sirius or ~/.config/sirius
    //   Windows: %APPDATA%/Sirius
    //   macOS:   ~/Library/Application Support/Sirius
    static fs::path UserConfigDirectory();

    // System configuration directory (Linux /etc/sirius, Windows %PROGRAMDATA%).
    static fs::path SystemConfigDirectory();

    // Resolve a resource (shader, texture) by searching the executable dir, the
    // new-tree shaders/ subdir, the legacy Sirius.Render subdir, and the source
    // tree. Returns the canonical path if found.
    static std::optional<fs::path> ResolveResource(const std::string& relative_path);

    // Config file search locations, highest priority first.
    static std::vector<fs::path> ConfigSearchPaths();

    // First existing config file from the search paths, if any.
    static std::optional<fs::path> FindConfigFile();

    // "Linux", "Windows", or "macOS".
    static std::string PlatformName();

    // Whether the process runs under WSL2.
    static bool IsWsl2();

  private:
    static fs::path executable_path_;
    static bool executable_path_initialised_;
};

}  // namespace sirius::app
