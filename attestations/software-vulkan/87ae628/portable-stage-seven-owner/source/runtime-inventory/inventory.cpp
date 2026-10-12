#include "sirius/backend/device.h"
#include <iomanip>
#include <iostream>
int main() {
    auto devices = sirius::backend::EnumerateVulkanDevices();
    if (!devices) { std::cerr << devices.error().Description(); return 1; }
    std::cout << "{\"devices\":[";
    bool first = true;
    for (const auto& info : *devices) {
        if (!first) std::cout << ',';
        first = false;
        std::cout << "{\"name\":" << std::quoted(info.name)
                  << ",\"driver_name\":" << std::quoted(info.driver_name)
                  << ",\"kind\":" << std::quoted(sirius::backend::ToString(info.kind))
                  << ",\"preserves_fp32_denormals\":" << std::boolalpha << info.preserves_fp32_denormals
                  << ",\"rounds_fp32_to_nearest\":" << info.rounds_fp32_to_nearest
                  << ",\"portable_owner_eligible\":" << (sirius::backend::RetainedUsesPortableArithmetic(info) && info.rounds_fp32_to_nearest)
                  << ",\"render_memory_bytes\":" << info.render_memory_bytes << '}';
    }
    std::cout << "]}\n";
}
