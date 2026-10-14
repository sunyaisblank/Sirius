// Setup-only reproduction of the original four-record ordinary-FP32 probe.
// No active record dispatch or scientific/numerical acceptance is performed.
#include "sirius/backend/device.h"
#include "sirius/core/disk/novikov_thorne_disk.h"
#include "support/test_resource.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void Marker(const char* event, const char* phase) {
    const auto ticks = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
    std::fprintf(stderr, "[SetupPhase] {\"event\":\"%s\",\"phase\":\"%s\",\"monotonic_ns\":%lld}\n",
                 event, phase, static_cast<long long>(ticks));
    std::fflush(stderr);
}
struct Scope {
    const char* phase;
    explicit Scope(const char* value) : phase(value) { Marker("begin", phase); }
    ~Scope() { Marker("end", phase); }
};
void Require(bool okay, const char* message) {
    if (!okay) throw std::runtime_error(message);
}
template <typename T>
void RequireResult(const T& value) {
    if (!value.has_value()) throw std::runtime_error(value.error().Description());
}

// Original Sample layout and first original P4 case fixtures. Only zero-active
// params and zero result words are written before the original initialization.
struct Sample {
    float metric_id = 0.0f;
    float p1 = 0.0f, p2 = 0.0f, p3 = 0.0f, p4 = 0.0f;
    float c0 = 0.0f, c1 = 0.0f, c2 = 0.0f;
    float u0 = 0.0f, u1 = 0.0f, u2 = 0.0f, u3 = 0.0f;
    float h = 0.0f, aux0 = 0.0f, aux1 = 0.0f, aux2 = 0.0f;
};
std::vector<Sample> OriginalSamples() {
    Sample s8;
    s8.p1 = 1.0f; s8.p2 = 0.9f; s8.c0 = 8.0f;
    s8.aux0 = 2.32088304f; s8.aux1 = 10000.0f;
    Sample s9;
    s9.p1 = 1.0f; s9.p2 = 0.0f; s9.c0 = 10.0f;
    s9.aux0 = 6.0f; s9.aux1 = 10000.0f;
    const float retrograde_isco = static_cast<float>(sirius::core::AccretionDiskD::ComputeIsco(-0.7));
    Sample s10;
    s10.p1 = 2.0f; s10.p2 = -0.7f; s10.c0 = 20.0f;
    s10.aux0 = 2.0f * retrograde_isco; s10.aux1 = 10000.0f;
    Sample s11;
    s11.p1 = 1.0f; s11.p2 = 0.0f; s11.c0 = 12.0f;
    s11.aux0 = 10.0f; s11.aux1 = 10000.0f;
    return {s8, s9, s10, s11};
}
std::vector<float> Flatten(const std::vector<Sample>& samples) {
    std::vector<float> flat;
    flat.reserve(samples.size() * 16);
    for (const auto& s : samples) {
        flat.insert(flat.end(), {s.metric_id, s.p1, s.p2, s.p3, s.p4, s.c0, s.c1, s.c2,
                                s.u0, s.u1, s.u2, s.u3, s.h, s.aux0, s.aux1, s.aux2});
    }
    return flat;
}
std::vector<std::uint32_t> LoadSpirv(const std::string& path) {
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    if (!file) return {};
    const auto size = static_cast<std::size_t>(file.tellg());
    std::vector<std::uint32_t> words(size / sizeof(std::uint32_t));
    file.seekg(0);
    file.read(reinterpret_cast<char*>(words.data()), static_cast<std::streamsize>(size));
    return words;
}
}  // namespace

int main() {
    const auto wall_begin = std::chrono::steady_clock::now();
    Marker("begin", "driver.whole_setup");
    try {
        // Same OpenProbe public order: enumerate, select, create, load resource,
        // LoadKernel. Failure still returns from its original public operation.
        const auto devices = sirius::backend::EnumerateVulkanDevices();
        RequireResult(devices);
        Require(!devices->empty(), "no Vulkan device");
        const auto selected = sirius::backend::ResolveVulkanDeviceIndex(*devices);
        RequireResult(selected);
        auto created = sirius::backend::CreateVulkanDevice(*selected);
        RequireResult(created);
        std::unique_ptr<sirius::backend::ComputeDevice> device = std::move(*created);
        std::vector<std::uint32_t> spirv;
        {
            Scope phase("driver.read_original_spirv");
            spirv = LoadSpirv(sirius::test::ResourcePath("kernels/parity_probe.spv"));
        }
        Require(spirv.size() == 1390559, "full original module is absent");
        const auto kernel = device->LoadKernel(spirv);
        RequireResult(kernel);

        const auto samples = OriginalSamples();
        constexpr std::size_t capacity = 4, sample_stride = 16, result_stride = 64;
        constexpr std::uint64_t allocation_limit = 1024 * 1024;
        const auto flat = Flatten(samples);
        std::vector<float> completed(samples.size() * result_stride, 0.0f);
        std::vector<float> chunk(capacity * result_stride, 0.0f);
        std::array<float, 4> params = {2.0f, 0.0f, float(sample_stride), float(result_stride)};
        RequireResult(device->SetBufferAllocationLimit(allocation_limit));
        const std::uint64_t resident_before = device->BufferAllocationBytes();
        const auto sbuf = device->CreateBuffer(capacity * sample_stride * sizeof(float), sirius::backend::BufferUsage::kStorage);
        RequireResult(sbuf);
        const auto rbuf = device->CreateBuffer(chunk.size() * sizeof(float), sirius::backend::BufferUsage::kStorage);
        RequireResult(rbuf);
        const auto pbuf = device->CreateBuffer(sizeof(params), sirius::backend::BufferUsage::kStorage);
        RequireResult(pbuf);
        constexpr std::uint64_t requested = capacity * (sample_stride + result_stride) * sizeof(float) + sizeof(params);
        Require(device->BufferAllocationBytes() >= resident_before + requested && device->BufferAllocationBytes() <= allocation_limit,
                "original allocation checks failed");
        const sirius::backend::BufferHandle bindings[] = {*sbuf, *rbuf, *pbuf};
        Require(device->Info().kind == sirius::backend::DeviceKind::kSoftware, "setup-only action requires selected software route");
        RequireResult(device->WriteBuffer(*pbuf, std::as_bytes(std::span<const float>(params))));
        RequireResult(device->WriteBuffer(*rbuf, std::as_bytes(std::span<const float>(chunk))));
        sirius::backend::DispatchTiming initialization;
        {
            Scope phase("driver.zero_record_dispatch");
            RequireResult(device->Dispatch(*kernel, bindings, 1, 1, 1, &initialization));
        }
        Require(std::isfinite(initialization.submit_wait_ms) && initialization.submit_wait_ms > 0,
                "original initialization timing checks failed");
        RequireResult(device->ReadBuffer(*rbuf, std::as_writable_bytes(std::span<float>(chunk))));
        Require(std::all_of(chunk.begin(), chunk.end(), [](float value) { return value == 0.0f; }),
                "original zero-record complement check failed");
        const auto info = device->Info();
        const auto allocated = device->BufferAllocationBytes();
        {
            Scope phase("driver.destroy_device");
            device.reset();
        }
        Marker("end", "driver.whole_setup");
        const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_begin).count();
        std::cout << std::setprecision(17)
                  << "{\"schema\":\"sirius-parity-setup-only-v1\",\"setup_completed\":true,\"active_records\":0,"
                  << "\"fixture_records\":" << samples.size() << ",\"flat_words\":" << flat.size()
                  << ",\"completed_words\":" << completed.size() << ",\"result_words\":" << chunk.size()
                  << ",\"module_words\":" << spirv.size() << ",\"allocated_bytes\":" << allocated
                  << ",\"device_index\":" << *selected << ",\"device_name\":" << std::quoted(info.name)
                  << ",\"driver_name\":" << std::quoted(info.driver_name) << ",\"driver_info\":" << std::quoted(info.driver_info)
                  << ",\"whole_setup_seconds_including_teardown\":" << wall
                  << ",\"pipeline_setup_ms\":" << initialization.pipeline_setup_ms
                  << ",\"submit_wait_ms\":" << initialization.submit_wait_ms
                  << ",\"dispatch_total_ms\":" << initialization.total_ms
                  << ",\"numerical_acceptance\":false,\"qualification_claimed\":false}\n" << std::flush;
        return 0;
    } catch (const std::exception& error) {
        Marker("failure", "driver.whole_setup");
        std::cerr << "setup-only failure: " << error.what() << std::endl;
        return 1;
    }
}
