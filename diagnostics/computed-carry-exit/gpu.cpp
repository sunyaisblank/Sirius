// Focused five-word accumulator witness; no device call without --dispatch.
#include "sirius/backend/device.h"
#include <array>
#include <dlfcn.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<std::uint32_t> Read(const char* path) {
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    const auto size = file.tellg();
    if (!file || size <= 0 || size > 8 * 1024 * 1024 || size % 4 != 0)
        throw std::runtime_error("invalid word file");
    std::vector<std::uint32_t> words(static_cast<std::size_t>(size) / 4);
    file.seekg(0);
    if (!file.read(reinterpret_cast<char*>(words.data()), size))
        throw std::runtime_error("word file read failed");
    return words;
}

int main(int argc, char** argv) {
    try {
        if (argc != 7 || std::string(argv[6]) != "--dispatch")
            throw std::runtime_error("usage: gpu module input expected output mode --dispatch");
        const std::string mode = argv[5];
        if (mode != "integer" && mode != "normal-transport")
            throw std::runtime_error("unknown mode");
        const auto module = Read(argv[1]), input = Read(argv[2]), expected = Read(argv[3]);
        if (module.front() != 0x07230203u || input.front() == 0 || input.front() > 65536 ||
            input.size() != 1u + 6u * input.front() || expected.size() != 5u * input.front() ||
            4u * (input.size() + expected.size()) > 8u * 1024u * 1024u)
            throw std::runtime_error("fixture/module shape or allocation bound mismatch");
        using namespace sirius::backend;
        std::cout << std::unitbuf;
        const auto inventory = EnumerateVulkanDevices();
        if (!inventory) throw std::runtime_error(inventory.error().Description());
        const auto selected = ResolveVulkanDeviceIndex(*inventory);
        if (!selected || *selected >= inventory->size())
            throw std::runtime_error("invalid selected device");
        auto opened = CreateVulkanDevice(*selected);
        if (!opened) throw std::runtime_error(opened.error().Description());
        auto& device = **opened;
        const auto& info = device.Info();
        if (info.kind != DeviceKind::kSoftware || info.driver_id != 13 || info.driver_name != "llvmpipe" ||
            (mode == "normal-transport" && !info.rounds_fp32_to_nearest))
            throw std::runtime_error("actual selected software provider/mode is not admitted");
        void* resident = dlopen("libvulkan_lvp.so", RTLD_NOW | RTLD_NOLOAD);
        Dl_info resident_info{};
        if (!resident || !dladdr(dlsym(resident, "vk_icdGetInstanceProcAddr"), &resident_info) ||
            !resident_info.dli_fname)
            throw std::runtime_error("resident software provider identity missing");
        const std::string resident_path = resident_info.dli_fname;
        dlclose(resident);
        std::cout << "{\"provider\":" << std::quoted(info.driver_name)
                  << ",\"device\":" << std::quoted(info.name)
                  << ",\"driver_id\":" << info.driver_id
                  << ",\"RTE32\":" << info.rounds_fp32_to_nearest
                  << ",\"resident_provider\":" << std::quoted(resident_path)
                  << ",\"mode\":" << std::quoted(mode) << "}\n";
        auto status = device.SetBufferAllocationLimit(8 * 1024 * 1024);
        if (!status) throw std::runtime_error(status.error().Description());
        const auto kernel = device.LoadKernel(module);
        if (!kernel) throw std::runtime_error(kernel.error().Description());
        const auto in = device.CreateBuffer(input.size() * 4, BufferUsage::kStorage);
        const auto out = device.CreateBuffer(expected.size() * 4, BufferUsage::kStorage);
        if (!in || !out) throw std::runtime_error("bounded buffers failed");
        std::vector<std::uint32_t> actual(expected.size());
        for (std::size_t i = 0; i < actual.size(); ++i) actual[i] = ~expected[i];
        status = device.WriteBuffer(*in, std::as_bytes(std::span(input)));
        if (!status) throw std::runtime_error(status.error().Description());
        status = device.WriteBuffer(*out, std::as_bytes(std::span(actual)));
        if (!status) throw std::runtime_error(status.error().Description());
        const std::array buffers{*in, *out};
        std::cout << "{\"dispatch_started\":1,\"cases\":" << input.front()
                  << ",\"complement_initialized_words\":" << actual.size() << "}\n";
        status = device.Dispatch(*kernel, buffers, (input.front() + 63u) / 64u, 1, 1);
        if (!status) throw std::runtime_error(status.error().Description());
        status = device.ReadBuffer(*out, std::as_writable_bytes(std::span(actual)));
        if (!status) throw std::runtime_error(status.error().Description());
        // Preserve complete failed output before deciding numerical acceptance.
        std::ofstream output(argv[4], std::ios::binary);
        output.write(reinterpret_cast<const char*>(actual.data()), static_cast<std::streamsize>(actual.size() * 4));
        if (!output) throw std::runtime_error("output write failed");
        output.close();
        std::size_t mismatches = 0, untouched = 0, malformed = 0;
        for (std::size_t i = 0; i < actual.size(); ++i) {
            untouched += actual[i] == ~expected[i];
            if (actual[i] != expected[i]) {
                if (mismatches < 8)
                    std::cerr << "mismatch row=" << i / 5 << " column=" << i % 5
                              << " expected=" << std::hex << expected[i]
                              << " actual=" << actual[i] << std::dec << '\n';
                ++mismatches;
            }
        }
        for (std::size_t row = 0; row < input.front(); ++row) {
            const auto base = 5 * row;
            malformed += actual[base + 4] > 1;
            for (std::size_t field = 0; field < 4; ++field)
                malformed += actual[base + 4] == 0 ? actual[base + field] != 0
                    : (actual[base + field] & 0x7fffffffu) > 0x7d800000u;
        }
        std::cout << "{\"completed_dispatches\":1,\"records\":" << input.front()
                  << ",\"complete_words\":" << actual.size() << ",\"mismatches\":" << mismatches
                  << ",\"untouched_complement_words\":" << untouched
                  << ",\"malformed_words\":" << malformed
                  << ",\"requested_bytes\":" << 4 * (input.size() + expected.size())
                  << ",\"resident_buffer_bytes\":" << device.BufferAllocationBytes() << "}\n";
        return mismatches == 0 && untouched == 0 && malformed == 0 ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
