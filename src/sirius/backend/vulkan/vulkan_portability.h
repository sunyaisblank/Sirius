#pragma once

// Internal creation boundary. Injectable Vulkan entry points let the CPU tests
// verify the exact extension names and flags sent to the loader and driver.
#include "sirius/backend/device.h"
#include "sirius/base/error.h"

#include <vulkan/vulkan.h>

namespace sirius::backend::detail {

#ifdef VK_KHR_shader_fma
using ShaderFmaFeatures = VkPhysicalDeviceShaderFmaFeaturesKHR;
inline constexpr auto kShaderFmaFeaturesType =
    VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SHADER_FMA_FEATURES_KHR;
#else
// Vulkan-Headers v1.4.357 ABI, retained for supported older SDKs:
// https://docs.vulkan.org/refpages/latest/refpages/source/VkPhysicalDeviceShaderFmaFeaturesKHR.html
struct ShaderFmaFeatures {
    VkStructureType sType;
    void* pNext;
    VkBool32 shaderFmaFloat16;
    VkBool32 shaderFmaFloat32;
    VkBool32 shaderFmaFloat64;
};
inline constexpr auto kShaderFmaFeaturesType = static_cast<VkStructureType>(1000579000);
#endif

[[nodiscard]] base::Expected<VkInstance> CreateInstanceWithPortability(
    VkInstanceCreateInfo create_info,
    PFN_vkEnumerateInstanceExtensionProperties enumerate = vkEnumerateInstanceExtensionProperties,
    PFN_vkCreateInstance create = vkCreateInstance);

[[nodiscard]] base::Expected<VkDevice> CreateDeviceWithPortability(
    VkPhysicalDevice physical, VkDeviceCreateInfo create_info,
    PFN_vkEnumerateDeviceExtensionProperties enumerate = vkEnumerateDeviceExtensionProperties,
    PFN_vkCreateDevice create = vkCreateDevice, DeviceInfo* device_info = nullptr,
    PFN_vkGetPhysicalDeviceFeatures2 get_features = vkGetPhysicalDeviceFeatures2);

}  // namespace sirius::backend::detail
