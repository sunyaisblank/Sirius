#include <vulkan/vulkan_core.h>
#include <cstddef>
#include <iostream>
int main() {
    std::cout << "{\"sType\":" << VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FLOAT_CONTROLS_PROPERTIES
        << ",\"bytes\":" << sizeof(VkPhysicalDeviceFloatControlsProperties)
        << ",\"first_field_offset\":" << offsetof(VkPhysicalDeviceFloatControlsProperties,denormBehaviorIndependence)
        << ",\"last_field_offset\":" << offsetof(VkPhysicalDeviceFloatControlsProperties,shaderRoundingModeRTZFloat64)
        << "}\n";
}
