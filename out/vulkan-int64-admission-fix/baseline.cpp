#include "sirius/backend/vulkan/vulkan_device.h"
#include <array>
#include <iostream>

int main() {
    // Independent SPIR-V wire fixture: Shader followed by Int64. This calls
    // only the precision validator, before any Vulkan instance or driver work.
    constexpr std::array<std::uint32_t, 9> words{
        0x07230203u, 0x00010500u, 0u, 1u, 0u, 0x00020011u, 1u, 0x00020011u, 11u};
    const auto result = sirius::backend::ValidateVulkanKernelPrecision(words, false);
    std::cout << "{\"requested_Int64\":true,\"logical_Int64_enabled\":false,"
                 "\"admitted\":" << (result ? "true" : "false") << ",\"Vulkan_calls\":0}\n";
}
