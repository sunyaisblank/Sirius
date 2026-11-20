#include "sirius/backend/vulkan/vulkan_device.h"

#include <array>
#include <bit>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <stdexcept>
#include <vector>

using namespace sirius::backend;

template<class T>
T Check(sirius::base::Expected<T> value) {
    if (!value) throw std::runtime_error(value.error().Description());
    if constexpr (!std::is_void_v<T>) return std::move(*value);
}

int main(int argc, char** argv) try {
    if (argc != 3 || std::filesystem::file_size(argv[2]) != 16 * sizeof(float))
        throw std::runtime_error("Usage: replay shader.spv sample.bin");
    std::array<float, 16> sample{};
    if (!std::ifstream(argv[2], std::ios::binary)
            .read(reinterpret_cast<char*>(sample.data()), sizeof(sample)))
        throw std::runtime_error("Read exact sample bytes");
    const auto size = std::filesystem::file_size(argv[1]);
    if (size == 0 || size % sizeof(std::uint32_t))
        throw std::runtime_error("Wrong SPIR-V length");
    std::vector<std::uint32_t> spirv(size / sizeof(std::uint32_t));
    if (!std::ifstream(argv[1], std::ios::binary)
            .read(reinterpret_cast<char*>(spirv.data()), static_cast<std::streamsize>(size)))
        throw std::runtime_error("Read SPIR-V");
    const auto inventory = Check(EnumerateVulkanDevices());
    auto device = Check(CreateVulkanDevice(Check(ResolveVulkanDeviceIndex(inventory))));
    std::cerr << "device=" << device->Info().name << '\n';
    const auto kernel = Check(device->LoadKernel(spirv));
    Check(device->SetBufferAllocationLimit(1024 * 1024));
    std::array<float, 64> results{};
    std::array<float, 4> parameters{8.0f, 0.0f, 16.0f, 64.0f};
    const std::array<BufferHandle, 3> buffers{
        Check(device->CreateBuffer(sizeof(sample), BufferUsage::kStorage)),
        Check(device->CreateBuffer(sizeof(results), BufferUsage::kStorage)),
        Check(device->CreateBuffer(sizeof(parameters), BufferUsage::kStorage))
    };
    Check(device->WriteBuffer(buffers[0], std::as_bytes(std::span{sample})));
    Check(device->WriteBuffer(buffers[1], std::as_bytes(std::span{results})));
    Check(device->WriteBuffer(buffers[2], std::as_bytes(std::span{parameters})));
    DispatchTiming initialization;
    Check(device->Dispatch(kernel, buffers, 1, 1, 1, &initialization));
    Check(device->ReadBuffer(buffers[1], std::as_writable_bytes(std::span{results})));
    for (float value : results) {
        if (value != 0.0f) throw std::runtime_error("Inactive row changed");
    }
    parameters[1] = 1.0f;
    Check(device->WriteBuffer(buffers[2], std::as_bytes(std::span{parameters})));
    DispatchTiming actual;
    Check(device->Dispatch(kernel, buffers, 1, 1, 1, &actual));
    Check(device->ReadBuffer(buffers[1], std::as_writable_bytes(std::span{results})));
    std::cerr << "init_ms=" << initialization.total_ms << " actual_submit_ms=" << actual.submit_wait_ms << '\n';
    if (!(actual.submit_wait_ms > 0.0 && actual.submit_wait_ms <= 1000.0))
        throw std::runtime_error("Actual dispatch outside unchanged governor limit");
    if (results[0] != 1.0f) throw std::runtime_error("Wrong case or shader ABI");
    for (std::size_t i = 0; i != 8; ++i) {
        const std::array<std::size_t, 8> offsets{15, 5, 6, 7, 8, 9, 10, 11};
        if (std::bit_cast<std::uint32_t>(results[20 + i]) !=
            std::bit_cast<std::uint32_t>(sample[offsets[i]]))
            throw std::runtime_error("Initial phase-space readback mismatch");
    }
    std::cout << std::setprecision(9) << '[';
    for (std::size_t i = 0; i != 28; ++i) {
        if (!std::isfinite(results[i])) throw std::runtime_error("Non-finite result");
        std::cout << (i ? ", " : "") << results[i];
    }
    std::cout << "]\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
