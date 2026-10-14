#pragma once

#include "sirius/base/resource_locator.h"

#include <stdexcept>
#include <string>
#include <string_view>

namespace sirius::test {

// Qualification inputs travel with the consuming executable. Never consult a
// source/build directory, the working directory, or a development override.
inline std::string ResourcePath(std::string_view relative) {
    const auto executable = base::ExecutableDirectory();
    if (executable.empty()) throw std::runtime_error("Cannot locate test executable");
    const auto root = executable / "resources";
    if (const auto path = base::ResolveResourceFromRoot(root, relative)) {
        return path->string();
    }
    throw std::runtime_error("Missing or escaping test resource: " + (root / relative).string());
}

}  // namespace sirius::test
