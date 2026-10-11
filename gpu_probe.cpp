// Adapted from the isolated portable-arithmetic gpu_probe.cpp. No device call
// occurs without --dispatch; the external runner supervises the synchronous call.
#include "sirius/backend/device.h"
#include "word_io.h"
#include <array>
#include <dlfcn.h>
#include <iomanip>
#include <iostream>
#include <string>
int main(int argc, char** argv) {
    try {
        if (argc != 5 && argc != 6)
            throw std::runtime_error("usage: gpu_probe module input expected actual [--dispatch]");
        const auto shader = ReadWords(argv[1]), input = ReadWords(argv[2]), expected = ReadWords(argv[3]);
        if (input.empty() || input.size() != 1 + 3ULL * input[0] || expected.size() != 5ULL * input[0]) throw std::runtime_error("upward probe shape mismatch");
        if (shader[0] != 0x07230203u) throw std::runtime_error("SPIR-V magic mismatch");
        if (argc == 5) {
            std::cout << "{\"scope\":\"isolated raw-word prototype preflight\",\"cases\":" << input[0]
                      << ",\"device_calls\":0,\"dispatches\":0}\n";
            return 0;
        }
        if (std::string(argv[5]) != "--dispatch") throw std::runtime_error("unknown argument");
        std::cout << std::unitbuf;
        using namespace sirius::backend;
        const auto inventory = EnumerateVulkanDevices();
        if (!inventory) throw std::runtime_error(inventory.error().Description());
        const auto selected = ResolveVulkanDeviceIndex(*inventory);
        if (!selected) throw std::runtime_error(selected.error().Description());
        if (*selected >= inventory->size()) throw std::runtime_error("device selection out of range");
        auto opened = CreateVulkanDevice(*selected);
        if (!opened) throw std::runtime_error(opened.error().Description());
        auto& device = **opened;
        const auto& info = device.Info();
        void* resident_driver = dlopen("libvulkan_lvp.so", RTLD_NOW | RTLD_NOLOAD);
        Dl_info resident_info{};
        if (!resident_driver || !dladdr(dlsym(resident_driver, "vk_icdGetInstanceProcAddr"), &resident_info) ||
            !resident_info.dli_fname)
            throw std::runtime_error("actual resident Lavapipe module identity missing");
        const std::string resident_path = resident_info.dli_fname;
        dlclose(resident_driver);
        std::cout << "{\"scope\":\"isolated raw-word prototype device identity\",\"device\":"
                  << std::quoted(info.name) << ",\"kind\":" << std::quoted(ToString(info.kind))
                  << ",\"driver\":" << std::quoted(info.driver_name)
                  << ",\"driver_info\":" << std::quoted(info.driver_info)
                  << ",\"resident_driver_path\":" << std::quoted(resident_path)
                  << ",\"vendor_id\":" << info.vendor_id << ",\"device_id\":" << info.device_id
                  << ",\"api_version\":" << info.api_version << ",\"driver_id\":" << info.driver_id
                  << ",\"denorm32\":" << info.preserves_fp32_denormals
                  << ",\"RTE32\":" << info.rounds_fp32_to_nearest
                  << ",\"fp64\":" << info.supports_fp64
                  << ",\"RTE64\":" << info.rounds_fp64_to_nearest << "}\n";
        auto status = device.SetBufferAllocationLimit(8 * 1024 * 1024);
        if (!status) throw std::runtime_error(status.error().Description());
        if (!info.rounds_fp32_to_nearest) throw std::runtime_error("upward-multiply probe requires RTE32");
        const auto kernel = device.LoadKernel(shader);
        if (!kernel) throw std::runtime_error(kernel.error().Description());
        const auto in = device.CreateBuffer(input.size() * 4, BufferUsage::kStorage);
        if (!in) throw std::runtime_error(in.error().Description());
        const auto out = device.CreateBuffer(expected.size() * 4, BufferUsage::kStorage);
        if (!out) throw std::runtime_error(out.error().Description());
        status = device.WriteBuffer(*in, std::as_bytes(std::span(input)));
        if (!status) throw std::runtime_error(status.error().Description());
        std::vector<std::uint32_t> actual(expected.size());
        for (std::size_t i = 0; i < expected.size(); ++i) actual[i] = ~expected[i];
        status = device.WriteBuffer(*out, std::as_bytes(std::span(actual)));
        if (!status) throw std::runtime_error(status.error().Description());
        const std::array buffers{*in, *out};
        DispatchTiming timing;
        std::cout << "{\"dispatch_started\":1,\"groups\":" << (input[0] + 63u) / 64u
                  << ",\"cases\":" << input[0] << ",\"complement_initialized_words\":" << actual.size() << "}\n";
        status = device.Dispatch(*kernel, buffers, (input[0] + 63u) / 64u, 1, 1, &timing);
        if (!status) throw std::runtime_error(status.error().Description());
        status = device.ReadBuffer(*out, std::as_writable_bytes(std::span(actual)));
        if (!status) throw std::runtime_error(status.error().Description());
        std::size_t mismatches = 0, untouched = 0;
        for (std::size_t i = 0; i < actual.size(); ++i) {
            untouched += actual[i] == ~expected[i];
            if (actual[i] == expected[i]) continue;
            if (mismatches < 8) std::cerr << "mismatch row=" << i / 5 << " col=" << i % 5
                                         << " expected=" << std::hex << expected[i]
                                         << " actual=" << actual[i] << std::dec << '\n';
            ++mismatches;
        }
        WriteWords(argv[4], actual);
        std::cout << std::setprecision(17)
                  << "{\"completed_dispatches\":1,\"cases\":" << input[0]
                  << ",\"observed_words\":" << actual.size() << ",\"mismatches\":" << mismatches
                  << ",\"untouched_complement_words\":" << untouched
                  << ",\"requested_bytes\":" << (input.size() + expected.size()) * 4
                  << ",\"actual_resident_bytes\":" << device.BufferAllocationBytes()
                  << ",\"host_dispatch_timing_ms\":{\"pipeline\":" << timing.pipeline_setup_ms
                  << ",\"command\":" << timing.command_setup_ms << ",\"submit_wait\":" << timing.submit_wait_ms
                  << ",\"cleanup\":" << timing.cleanup_ms << ",\"total\":" << timing.total_ms << "}}\n";
        return mismatches == 0 && untouched == 0 ? 0 : 1;
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
