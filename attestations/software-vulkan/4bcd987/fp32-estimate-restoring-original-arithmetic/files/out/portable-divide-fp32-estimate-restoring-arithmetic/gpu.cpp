#include "sirius/backend/device.h"
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<std::uint32_t> Read(const char* path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    const auto size = f.tellg();
    if (!f || size <= 0 || size > 8 * 1024 * 1024 || size % 4 != 0)
        throw std::runtime_error("invalid word file");
    std::vector<std::uint32_t> words(static_cast<std::size_t>(size) / 4);
    f.seekg(0);
    if (!f.read(reinterpret_cast<char*>(words.data()), size))
        throw std::runtime_error("word file read failed");
    return words;
}
int main(int argc, char** argv) {
    try {
        if (argc != 5) throw std::runtime_error("usage: gpu module input expected output");
        const auto module = Read(argv[1]), input = Read(argv[2]), expected = Read(argv[3]);
        if (module.front() != 0x07230203u || input.front() == 0 ||
            input.size() != 1u + 4u * input.front() || expected.size() != 4u * input.front())
            throw std::runtime_error("fixture/module shape mismatch");
        using namespace sirius::backend;
        const auto inventory = EnumerateVulkanDevices();
        if (!inventory) throw std::runtime_error(inventory.error().Description());
        const auto selected = ResolveVulkanDeviceIndex(*inventory);
        if (!selected) throw std::runtime_error(selected.error().Description());
        auto opened = CreateVulkanDevice(*selected);
        if (!opened) throw std::runtime_error(opened.error().Description());
        auto& device = **opened;
        const auto& info = device.Info();
        if (info.kind != DeviceKind::kSoftware || info.driver_id != 13 || info.driver_name != "llvmpipe")
            throw std::runtime_error("actual selected provider differs from pinned software route");
        auto status = device.SetBufferAllocationLimit(8 * 1024 * 1024);
        if (!status) throw std::runtime_error(status.error().Description());
        const auto kernel = device.LoadKernel(module);
        if (!kernel) throw std::runtime_error(kernel.error().Description());
        const auto in = device.CreateBuffer(input.size() * 4, BufferUsage::kStorage);
        const auto out = device.CreateBuffer(expected.size() * 4, BufferUsage::kStorage);
        if (!in || !out) throw std::runtime_error("bounded buffers failed");
        std::vector<std::uint32_t> actual(expected.size(), 0xdadbaabau);
        status = device.WriteBuffer(*in, std::as_bytes(std::span(input)));
        if (!status) throw std::runtime_error(status.error().Description());
        status = device.WriteBuffer(*out, std::as_bytes(std::span(actual)));
        if (!status) throw std::runtime_error(status.error().Description());
        const std::array buffers{*in, *out};
        status = device.Dispatch(*kernel, buffers, (input.front() + 63) / 64, 1, 1);
        if (!status) throw std::runtime_error(status.error().Description());
        status = device.ReadBuffer(*out, std::as_writable_bytes(std::span(actual)));
        if (!status) throw std::runtime_error(status.error().Description());
        if (actual != expected) {
            for (std::size_t i = 0; i < actual.size(); ++i)
                if (actual[i] != expected[i]) {
                    std::cerr << "mismatch word " << i << " expected " << std::hex
                              << expected[i] << " actual " << actual[i] << '\n';
                    break;
                }
            throw std::runtime_error("complete independent raw output differs");
        }
        std::ofstream output(argv[4], std::ios::binary);
        output.write(reinterpret_cast<const char*>(actual.data()), actual.size() * 4);
        if (!output) throw std::runtime_error("output write failed");
        std::cout << "{\"records\":" << input.front() << ",\"complete_words\":" << actual.size()
                  << ",\"mismatches\":0,\"device\":" << std::quoted(info.name)
                  << ",\"driver\":" << std::quoted(info.driver_name) << "}\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
