// Compiler diagnostic only: Endpoint's unchanged ni<329 guard returns before
// any retained arithmetic or output writes for this one-word input.
#include "sirius/backend/device.h"
#include <array>
#include <dlfcn.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::runtime_error("usage: short_input_jit_probe module");
        std::ifstream file(argv[1], std::ios::binary | std::ios::ate);
        if (!file) throw std::runtime_error("module open failed");
        const auto bytes = file.tellg();
        if (bytes <= 0 || bytes > 8 * 1024 * 1024 || bytes % 4 != 0)
            throw std::runtime_error("module size invalid");
        std::vector<std::uint32_t> shader(static_cast<std::size_t>(bytes) / 4);
        file.seekg(0);
        if (!file.read(reinterpret_cast<char*>(shader.data()), bytes) || shader.front() != 0x07230203u)
            throw std::runtime_error("module read/magic invalid");
        using namespace sirius::backend;
        const auto inventory = EnumerateVulkanDevices();
        if (!inventory) throw std::runtime_error(inventory.error().Description());
        const auto selected = ResolveVulkanDeviceIndex(*inventory);
        if (!selected || *selected >= inventory->size()) throw std::runtime_error("device selection invalid");
        auto opened = CreateVulkanDevice(*selected);
        if (!opened) throw std::runtime_error(opened.error().Description());
        auto& device = **opened;
        const auto& info = device.Info();
        if (info.kind != DeviceKind::kSoftware || info.driver_id != 13 ||
            info.driver_name != "llvmpipe" || !info.rounds_fp32_to_nearest)
            throw std::runtime_error("actual software RTE32 provider mismatch");
        void* driver = dlopen("libvulkan_lvp.so", RTLD_NOW | RTLD_NOLOAD);
        Dl_info resident{};
        if (!driver || !dladdr(dlsym(driver, "vk_icdGetInstanceProcAddr"), &resident) || !resident.dli_fname)
            throw std::runtime_error("resident provider missing");
        const std::string provider = resident.dli_fname;
        dlclose(driver);
        auto status = device.SetBufferAllocationLimit(8 * 1024 * 1024);
        if (!status) throw std::runtime_error(status.error().Description());
        const auto kernel = device.LoadKernel(shader);
        if (!kernel) throw std::runtime_error(kernel.error().Description());
        const auto in = device.CreateBuffer(4, BufferUsage::kStorage);
        const auto out = device.CreateBuffer(4, BufferUsage::kStorage);
        if (!in || !out) throw std::runtime_error("diagnostic buffers failed");
        const std::array<std::uint32_t, 1> input{1};
        std::array<std::uint32_t, 1> output{0xdadbaabau};
        status = device.WriteBuffer(*in, std::as_bytes(std::span(input)));
        if (!status) throw std::runtime_error(status.error().Description());
        status = device.WriteBuffer(*out, std::as_bytes(std::span(output)));
        if (!status) throw std::runtime_error(status.error().Description());
        std::cout << std::unitbuf << "{\"scope\":\"short-input compiler diagnostic\",\"device\":"
                  << std::quoted(info.name) << ",\"provider\":" << std::quoted(provider)
                  << ",\"input_words\":1,\"rows_header\":1,\"output_words\":1,\"groups\":1}\n";
        const std::array buffers{*in, *out};
        status = device.Dispatch(*kernel, buffers, 1, 1, 1);
        if (!status) throw std::runtime_error(status.error().Description());
        status = device.ReadBuffer(*out, std::as_writable_bytes(std::span(output)));
        if (!status) throw std::runtime_error(status.error().Description());
        if (output[0] != 0xdadbaabau) throw std::runtime_error("short-input output canary changed");
        std::cout << "{\"completed_dispatches\":1,\"complete_canary_unchanged\":true,"
                     "\"retained_numerical_path_executed\":false,\"resident_buffer_bytes\":"
                  << device.BufferAllocationBytes() << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
