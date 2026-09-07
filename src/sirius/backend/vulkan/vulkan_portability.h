#pragma once

// Internal creation boundary. Injectable Vulkan entry points let the CPU tests
// verify the exact extension names and flags sent to the loader and driver.
#include "sirius/base/error.h"

#include <vulkan/vulkan.h>

namespace sirius::backend::detail {

[[nodiscard]] base::Expected<VkInstance> CreateInstanceWithPortability(
    VkInstanceCreateInfo create_info,
    PFN_vkEnumerateInstanceExtensionProperties enumerate = vkEnumerateInstanceExtensionProperties,
    PFN_vkCreateInstance create = vkCreateInstance);

[[nodiscard]] base::Expected<VkDevice> CreateDeviceWithPortability(
    VkPhysicalDevice physical, VkDeviceCreateInfo create_info,
    PFN_vkEnumerateDeviceExtensionProperties enumerate = vkEnumerateDeviceExtensionProperties,
    PFN_vkCreateDevice create = vkCreateDevice);

}  // namespace sirius::backend::detail
