#pragma once

#include "sirius/backend/vulkan/vulkan_device.h"

#include <cstdlib>
#include <string_view>

namespace sirius::backend::detail {

[[nodiscard]] inline bool PackagedSoftwareRequested() {
    const char* value = std::getenv("SIRIUS_PACKAGED_SOFTWARE_VULKAN");
    return value != nullptr && std::string_view(value) == "1";
}

#ifdef SIRIUS_HAS_PACKAGED_SOFTWARE_VULKAN
base::Expected<VkInstance> CreateProbedInstance(VkInstanceCreateInfo);
std::span<const std::uint32_t> SelectSealedRetainedModule(const DeviceInfo&,
                                                          std::span<const std::uint32_t>);
bool UsesSealedRetainedModules(const DeviceInfo&);
#endif

}  // namespace sirius::backend::detail
