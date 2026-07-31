#pragma once

// Relocatable runtime-resource discovery shared by the CLI and render layers.
// Installed resources live below <prefix>/share/sirius; build-tree resources
// are staged below <executable>/resources. Candidates are always anchored to
// the executable volume; the working directory is never searched. Setting
// SIRIUS_RESOURCE_DIR selects one explicit resource root. Absolute/traversing
// names and symlinks escaping that root fail closed.

#include <filesystem>
#include <optional>
#include <string_view>
#include <vector>

namespace sirius::base {

namespace fs = std::filesystem;

[[nodiscard]] fs::path ExecutableDirectory();

[[nodiscard]] std::vector<fs::path> ResourceCandidates(std::string_view relative_path);

[[nodiscard]] std::optional<fs::path> ResolveResource(std::string_view relative_path);

}  // namespace sirius::base
