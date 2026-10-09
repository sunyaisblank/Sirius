// Preparation only. No device is enumerated unless --dispatch is supplied.
// Coordinator approval of a noninterfering window is required before running it.
#include "sirius/backend/device.h"

#include <array>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<std::uint32_t> ReadWords(const char* path) {
    std::ifstream stream(path, std::ios::binary | std::ios::ate);
    const auto size = stream.tellg();
    if (!stream || size <= 0 || size > 4 * 1024 * 1024 || size % 4 != 0)
        throw std::runtime_error("invalid/bounded word input");
    std::vector<unsigned char> bytes(static_cast<std::size_t>(size));
    stream.seekg(0);
    if (!stream.read(reinterpret_cast<char*>(bytes.data()), size))
        throw std::runtime_error("input read failed");
    std::vector<std::uint32_t> words(bytes.size() / 4);
    for (std::size_t i = 0; i < words.size(); ++i)
        for (unsigned j = 0; j < 4; ++j)
            words[i] |= std::uint32_t(bytes[4 * i + j]) << (8 * j);
    return words;
}

int main(int argc, char** argv) {
    try {
        if (argc != 4 && argc != 5) {
            std::cerr << "usage: gpu_probe module.spv shader-input.bin shader-expected.bin [--dispatch]\n";
            return 2;
        }
        const auto shader = ReadWords(argv[1]);
        const auto input = ReadWords(argv[2]);
        const auto expected = ReadWords(argv[3]);
        if (shader[0] != 0x07230203u || input[0] == 0 || input[0] > 262144u ||
            input.size() != 1 + 3 * std::size_t(input[0]) || expected.size() != input[0])
            throw std::runtime_error("raw-word shape mismatch");
        std::cout << "scope=raw integer primitive diagnostic; not retained/product qualification\n"
                  << "cases=" << input[0] << " input_bytes=" << input.size() * 4
                  << " output_bytes=" << expected.size() * 4 << '\n';
        if (argc == 4) {
            std::cout << "device_calls=0 dispatches=0 (preflight only)\n";
            return 0;
        }
        if (std::string(argv[4]) != "--dispatch") throw std::runtime_error("unknown argument");
        using namespace sirius::backend;
        const auto inventory = EnumerateVulkanDevices();
        if (!inventory) throw std::runtime_error(inventory.error().Description());
        const auto selected = ResolveVulkanDeviceIndex(*inventory);
        if (!selected) throw std::runtime_error(selected.error().Description());
        auto opened = CreateVulkanDevice(*selected);
        if (!opened) throw std::runtime_error(opened.error().Description());
        auto& device = **opened;
        const auto& info = device.Info();
        std::cout << "device=" << std::quoted(info.name) << " driver=" << std::quoted(info.driver_name)
                  << " driver_info=" << std::quoted(info.driver_info)
                  << " vendor_id=" << info.vendor_id << " device_id=" << info.device_id
                  << " denorm32=" << info.preserves_fp32_denormals
                  << " RTE32=" << info.rounds_fp32_to_nearest << '\n';
        auto status = device.SetBufferAllocationLimit(8 * 1024 * 1024);
        if (!status) throw std::runtime_error(status.error().Description());
        const auto kernel = device.LoadKernel(shader);
        if (!kernel) throw std::runtime_error(kernel.error().Description());
        const auto in = device.CreateBuffer(input.size() * 4, BufferUsage::kStorage);
        if (!in) throw std::runtime_error(in.error().Description());
        const auto out = device.CreateBuffer(expected.size() * 4, BufferUsage::kStorage);
        if (!out) throw std::runtime_error(out.error().Description());
        status = device.WriteBuffer(*in, std::as_bytes(std::span(input)));
        if (!status) throw std::runtime_error(status.error().Description());
        const std::array buffers{*in, *out};
        DispatchTiming timing;
        status = device.Dispatch(*kernel, buffers, (input[0] + 63u) / 64u, 1, 1, &timing);
        if (!status) throw std::runtime_error(status.error().Description());
        std::vector<std::uint32_t> actual(expected.size());
        status = device.ReadBuffer(*out, std::as_writable_bytes(std::span(actual)));
        if (!status) throw std::runtime_error(status.error().Description());
        std::size_t mismatches = 0;
        for (std::size_t i = 0; i < actual.size(); ++i) {
            if (actual[i] == expected[i]) continue;
            if (mismatches < 10)
                std::cout << "mismatch row=" << i << " expected=" << std::hex << expected[i]
                          << " actual=" << actual[i] << std::dec << '\n';
            ++mismatches;
        }
        std::cout << "completed_dispatches=1 observed_words=" << actual.size()
                  << " mismatches=" << mismatches
                  << " actual_resident_bytes=" << device.BufferAllocationBytes() << '\n';
        return mismatches == 0 ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
