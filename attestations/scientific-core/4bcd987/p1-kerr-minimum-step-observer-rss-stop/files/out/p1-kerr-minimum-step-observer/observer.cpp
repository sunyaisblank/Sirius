#include "support/test_resource.h"
#include "support/trace_continuation_probe.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace {
using Record = sirius::render::FloatTraceContinuation;
constexpr std::uint32_t kWidth = 64;
constexpr std::uint32_t kHeight = 64;

// Exact original #823 parameter assignments, without a shader/controller change.
std::vector<float> OriginalParameters() {
    std::vector<float> params(sirius::render::kTraceParameterCount, 0.0f);
    params[44] = 0.5f;
    params[45] = 0.5f;
    params[66] = 0.5f;
    params[67] = 0.5f;
    params[51] = 1.0f;
    params[0] = kWidth;   // imageWidth
    params[1] = kHeight;  // imageHeight
    params[2] = 0.0f;     // metricId (Kerr-Schild family)
    params[3] = 1.0f;     // M
    params[4] = 0.9f;     // a
    params[5] = 0.0f;     // Q
    params[6] = 0.0f;     // Lambda
    params[7] = 0.0f;
    params[8] = 0.0f;
    params[9] = -40.0f;  // camera position
    params[10] = 0.0f;
    params[11] = 0.0f;
    params[12] = 1.0f;  // forward
    params[13] = 1.0f;
    params[14] = 0.0f;
    params[15] = 0.0f;  // right
    params[16] = 0.0f;
    params[17] = 1.0f;
    params[18] = 0.0f;            // up
    params[19] = 0.6f;            // fov (radians)
    params[20] = 1.0f;            // aspect
    params[21] = 2000.0f;         // max_steps
    params[22] = 0.06f;           // stepScale
    params[23] = 0.02f;           // min_step
    params[24] = 2.0f;            // max_step
    params[25] = 100.0f;          // escape_radius
    params[26] = 1.12f;           // captureFactor
    params[27] = 0.0f;            // disk disabled
    params[31] = 0.0f;            // tileOriginX
    params[32] = 0.0f;            // tileOriginY
    params[33] = float(kWidth);   // tileWidth
    params[34] = float(kHeight);  // tileHeight
    params[35] = 0.0f;            // starfield disabled (gradient background)
    params[36] = 1.0f;            // starfieldWidth (dummy)
    params[37] = 1.0f;            // starfieldHeight (dummy)
    return params;
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

std::ofstream NewOutput(const std::filesystem::path& path, bool binary = false) {
    if (std::filesystem::exists(path)) throw std::runtime_error("output already exists: " + path.string());
    std::ofstream stream(path, binary ? std::ios::binary : std::ios::out);
    stream.exceptions(std::ios::badbit | std::ios::failbit);
    return stream;
}

void WriteBytes(const std::filesystem::path& path, std::span<const std::byte> bytes) {
    auto stream = NewOutput(path, true);
    stream.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    stream.close();
}

void String(std::ostream& out, std::string_view value) {
    out << '"';
    for (const unsigned char ch : value) {
        if (ch == '"' || ch == '\\') out << '\\' << static_cast<char>(ch);
        else if (ch < 0x20) out << std::format("\\u{:04x}", ch);
        else out << static_cast<char>(ch);
    }
    out << '"';
}

template <typename T, std::size_t N>
void Array(std::ostream& out, const std::array<T, N>& values) {
    out << '[';
    for (std::size_t i = 0; i < N; ++i) {
        if (i) out << ',';
        if constexpr (std::is_floating_point_v<T>) {
            if (std::isfinite(values[i])) out << std::setprecision(std::numeric_limits<T>::max_digits10) << values[i];
            else out << "null";
        } else out << values[i];
    }
    out << ']';
}

void FloatBits(std::ostream& out, const std::array<float, 4>& values) {
    std::array<std::uint32_t, 4> bits{};
    for (std::size_t i = 0; i < 4; ++i) bits[i] = std::bit_cast<std::uint32_t>(values[i]);
    Array(out, bits);
}

bool SaveFailure(const std::filesystem::path& folder,
                 sirius::test::TraceContinuationProbe<float>& probe,
                 const sirius::backend::DeviceInfo& device,
                 std::span<const float> row_parameters,
                 const sirius::base::Error& error, std::uint32_t completed_rows) {
    const auto records = probe.Records();
    WriteBytes(folder / "failed-row-request-parameters.bin", std::as_bytes(row_parameters));
    WriteBytes(folder / "current-records.bin", std::as_bytes(records));
    auto out = NewOutput(folder / "observation.json");
    out << "{\n\"schema\":\"sirius-original-kerr-continuation-observation-v1\",\n"
           "\"original_case\":\"KernelTrace.KerrRenderIsFiniteNonConstantWithBoundedShadow\",\n"
           "\"numerical_acceptance\":false,\"complete_image_published\":false,\n"
           "\"original_failure_preserved\":true,\"failure\":{\"domain\":";
    String(out, sirius::base::ToString(error.domain()));
    out << ",\"operation\":"; String(out, error.operation());
    out << ",\"detail\":"; String(out, error.detail());
    out << "},\n\"device\":{\"name\":"; String(out, device.name);
    out << ",\"driver_name\":"; String(out, device.driver_name);
    out << ",\"driver_info\":"; String(out, device.driver_info);
    out << ",\"kind\":"; String(out, sirius::backend::ToString(device.kind));
    out << ",\"vendor_id\":" << device.vendor_id << ",\"device_id\":" << device.device_id << "},\n"
        << "\"completed_original_rows\":" << completed_rows << ",\n"
        << "\"probe_region_submissions_excluding_prepare_warmup\":" << probe.PhysicalSubmissions() << ",\n"
        << "\"active_record_count\":" << records.size() << ",\n"
           "\"record_stride\":" << sizeof(Record) << ",\n"
           "\"previous_record_available\":false,\"attempted_error_metrics_available\":false,\n"
           "\"last_dt_means_last_accepted_interval\":true,\n"
           "\"parameters_are_request_before_probe_action_and_epoch_assignment\":true,\n"
           "\"extra_observer_reads_or_dispatches\":0,\n"
           "\"failed_region_finalized\":false,\n\"minimum_step_record_indices\":[";
    bool first = true;
    for (std::size_t i = 0; i < records.size(); ++i) {
        if (records[i].control[0] == static_cast<std::uint32_t>(sirius::render::TracePhase::Terminal) &&
            records[i].control[1] == static_cast<std::uint32_t>(sirius::render::TraceTermination::MinimumStep)) {
            if (!first) out << ',';
            out << i;
            first = false;
        }
    }
    out << "],\n\"records\":[\n";
    for (std::size_t i = 0; i < records.size(); ++i) {
        if (i) out << ",\n";
        const auto& record = records[i];
        out << "{\"active_index\":" << i << ",\"identity\":"; Array(out, record.identity);
        out << ",\"control\":"; Array(out, record.control);
        out << ",\"physical\":[";
        for (std::size_t column = 0; column < record.physical.size(); ++column) {
            if (column) out << ',';
            Array(out, record.physical[column]);
        }
        out << "],\"physical_bits\":[";
        for (std::size_t column = 0; column < record.physical.size(); ++column) {
            if (column) out << ',';
            FloatBits(out, record.physical[column]);
        }
        out << "],\"source\":"; Array(out, record.source);
        out << ",\"source_bits\":"; FloatBits(out, record.source);
        out << ",\"volume\":"; Array(out, record.volume);
        out << ",\"volume_bits\":"; FloatBits(out, record.volume);
        out << '}';
    }
    out << "\n]}\n";
    out.close();
    return !first;
}

// Exact original TraceContinuationImage row order and Probe execution; only the
// failure branch observes already-cached Records before returning that error.
sirius::base::Expected<std::vector<float>> ObserveOriginalImage(
    sirius::backend::ComputeDevice& device, sirius::backend::KernelHandle kernel,
    std::span<const float> original, std::uint32_t width, std::uint32_t height,
    const std::filesystem::path& folder, bool& saved_minimum_step_record) {
    sirius::test::TraceContinuationProbe<float> probe(device, kernel);
    auto prepared = probe.Prepare();
    if (!prepared) return std::unexpected(prepared.error());
    if (original.size() != sirius::render::kTraceParameterCount)
        return sirius::base::Fail(sirius::base::ErrorDomain::kKernel, "direct trace image", "wrong parameter count");
    std::vector<float> result(static_cast<std::size_t>(width) * height * 4);
    std::vector<float> packet(original.begin(), original.end());
    for (std::uint32_t y = 0; y < height; ++y) {
        for (std::uint32_t x = 0; x < width; x += 64) {
            const auto count = std::min(64U, width - x);
            packet[31] = original[31] + static_cast<float>(x);
            packet[32] = original[32] + static_cast<float>(y);
            packet[33] = static_cast<float>(count);
            packet[34] = 1.0f;
            auto values = probe.RunRegion(packet);
            if (!values) {
                saved_minimum_step_record = SaveFailure(folder, probe, device.Info(), packet, values.error(), y);
                return std::unexpected(values.error());
            }
            std::copy(values->begin(), values->end(),
                      result.begin() + static_cast<std::ptrdiff_t>(
                                           (static_cast<std::size_t>(y) * width + x) * 4));
        }
    }
    return result;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::runtime_error("one existing owned output directory required");
        const std::filesystem::path folder(argv[1]);
        if (!std::filesystem::is_directory(folder)) throw std::runtime_error("owned output directory absent");
        const auto parameters = OriginalParameters();
        WriteBytes(folder / "original-request-parameters.bin", std::as_bytes(std::span(parameters)));
        const auto devices = sirius::backend::EnumerateVulkanDevices();
        if (!devices) throw std::runtime_error(devices.error().Description());
        if (devices->empty()) throw std::runtime_error("no Vulkan device present");
        const auto selected = sirius::backend::ResolveVulkanDeviceIndex(*devices);
        if (!selected) throw std::runtime_error(selected.error().Description());
        auto device = sirius::backend::CreateVulkanDevice(*selected);
        if (!device) throw std::runtime_error(device.error().Description());
        const auto spirv = LoadSpirv(sirius::test::ResourcePath("kernels/trace.spv"));
        if (spirv.empty()) throw std::runtime_error("trace.spv missing or empty");
        const auto kernel = (*device)->LoadKernel(spirv);
        if (!kernel) throw std::runtime_error(kernel.error().Description());
        bool saved_minimum_step_record = false;
        const auto observed = ObserveOriginalImage(**device, *kernel, parameters, kWidth, kHeight,
                                                   folder, saved_minimum_step_record);
        if (!observed) {
            std::cerr << observed.error().Description() << '\n';
            const bool captured_original_decline = saved_minimum_step_record &&
                observed.error().domain() == sirius::base::ErrorDomain::kKernel &&
                observed.error().operation() == "validate trace continuation" &&
                observed.error().detail() == "trace termination 5 cannot publish radiance";
            // Only a fully saved original MinimumStep refusal satisfies the
            // diagnostic owner. Every other error stops as incomplete.
            return captured_original_decline ? 1 : 2;
        }
        std::cerr << "original image completed without the expected failed record; no qualification claim\n";
        return 0;  // Owner expects the preserved decline, so this stops for review.
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
