#include <vulkan/vulkan.h>
#include <cstdio>
#include <vector>

int main() {
    VkApplicationInfo app{};
    app.sType = VK_STRUCTURE_TYPE_APPLICATION_INFO;
    app.apiVersion = VK_API_VERSION_1_2;
    VkInstanceCreateInfo create{};
    create.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
    create.pApplicationInfo = &app;
    VkInstance instance{};
    VkResult result = vkCreateInstance(&create, nullptr, &instance);
    if (result != VK_SUCCESS) { std::fprintf(stderr, "instance result=%d\n", result); return 1; }
    uint32_t count = 0;
    result = vkEnumeratePhysicalDevices(instance, &count, nullptr);
    if (result != VK_SUCCESS || count == 0) { vkDestroyInstance(instance, nullptr); return 2; }
    std::vector<VkPhysicalDevice> devices(count);
    result = vkEnumeratePhysicalDevices(instance, &count, devices.data());
    if (result != VK_SUCCESS) { vkDestroyInstance(instance, nullptr); return 3; }
    for (uint32_t i = 0; i < count; ++i) {
        VkPhysicalDeviceDriverProperties driver{};
        driver.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_DRIVER_PROPERTIES;
        VkPhysicalDeviceFloatControlsProperties fp{};
        fp.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FLOAT_CONTROLS_PROPERTIES;
        fp.pNext = &driver;
        VkPhysicalDeviceProperties2 properties{};
        properties.sType = VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2;
        properties.pNext = &fp;
        vkGetPhysicalDeviceProperties2(devices[i], &properties);
        VkPhysicalDeviceFeatures features{};
        vkGetPhysicalDeviceFeatures(devices[i], &features);
        std::printf("device=%u name=%s vendor=%u device_id=%u api=%u driver_id=%d driver=%s info=%s fp32_denorm_preserve=%u fp32_rte=%u fp64=%u fp64_denorm_preserve=%u fp64_rte=%u\n",
                    i, properties.properties.deviceName, properties.properties.vendorID,
                    properties.properties.deviceID, properties.properties.apiVersion,
                    driver.driverID, driver.driverName, driver.driverInfo,
                    fp.shaderDenormPreserveFloat32, fp.shaderRoundingModeRTEFloat32,
                    features.shaderFloat64, fp.shaderDenormPreserveFloat64,
                    fp.shaderRoundingModeRTEFloat64);
    }
    vkDestroyInstance(instance, nullptr);
}
