#pragma once

// Relocatable runtime-resource discovery shared by the CLI and render layers.
// Installed resources live below <prefix>/share/sirius; build-tree resources
// are staged below <executable>/resources. Candidates are always anchored to
// the executable volume; the working directory is never searched. Setting
// SIRIUS_RESOURCE_DIR selects one explicit resource root in development builds.
// Release builds ignore the override and remain bound to the executable volume.
// Absolute/traversing names and symlinks escaping the selected root fail closed.

#include <filesystem>
#include <optional>
#include <string_view>
#include <vector>

#ifndef SIRIUS_RELEASE_RESOURCE_LOCKED
#error "Sirius resource lookup requires an explicit alignment-derived release policy"
#endif
#if SIRIUS_RELEASE_RESOURCE_LOCKED != 0 && SIRIUS_RELEASE_RESOURCE_LOCKED != 1
#error "SIRIUS_RELEASE_RESOURCE_LOCKED must be exactly 0 or 1"
#endif

namespace sirius::base {

namespace fs = std::filesystem;

[[nodiscard]] fs::path ExecutablePath();
[[nodiscard]] fs::path ExecutableDirectory();

[[nodiscard]] std::vector<fs::path> ResourceCandidates(std::string_view relative_path);

// Resolve one safe relative resource below an already selected volume root.
// This is the common containment authority used by runtime discovery and by
// explicit-volume integrity verification.
[[nodiscard]] std::optional<fs::path> ResolveResourceFromRoot(const fs::path& root,
                                                              std::string_view relative_path);

[[nodiscard]] std::optional<fs::path> ResolveResource(std::string_view relative_path);

}  // namespace sirius::base
