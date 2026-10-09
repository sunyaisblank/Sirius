#include "sirius/backend/retained_compute.h"
#include "sirius/core/twofold.h"
#include "retained_kernels.h"
#include "support/retained_transport/reference_cases.h"
#include "word_io.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <dlfcn.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace {
using namespace sirius::backend;
using namespace sirius::backend::retained_program;
std::filesystem::path module_root, evidence_root;

class ModuleSubstitutionDevice final : public ComputeDevice {
  public:
    struct Module {
        std::string kind;
        std::span<const std::uint32_t> original;
        std::vector<std::uint32_t> candidate;
        unsigned loads = 0;
        std::optional<KernelHandle> handle;
    };
    ModuleSubstitutionDevice(std::unique_ptr<ComputeDevice> actual, bool wide,
                             std::filesystem::path directory)
        : actual_(std::move(actual)), directory_(std::move(directory)) {
        const std::string suffix = wide ? "_portable_fp64_products.spv" : "_portable_fp32_products.spv";
        const auto add = [&](const char* kind, const char* stem,
                             std::span<const std::uint32_t> narrow,
                             std::span<const std::uint32_t> wide_words) {
            auto candidate = ReadWords((module_root/(std::string("retained_") + stem + suffix)).c_str());
            if (candidate.size() < 5 || candidate[0] != 0x07230203u)
                throw std::runtime_error("candidate module malformed");
            modules.push_back({kind, wide ? wide_words : narrow, std::move(candidate), 0, std::nullopt});
        };
        add("Camera", "camera", kCameraPortableShader, kCameraPortableFp64Shader);
        add("Transport", "transport", kTransportPortableShader, kTransportPortableFp64Shader);
        add("Endpoint", "endpoint", kEndpointPortableShader, kEndpointPortableFp64Shader);
        add("Dense", "dense", kDensePortableShader, kDensePortableFp64Shader);
        add("Initialize", "initialize", kInitializePortableShader, kInitializePortableFp64Shader);
        add("RayCamera", "ray_camera", kRayCameraPortableShader, kRayCameraPortableFp64Shader);
    }
    const DeviceInfo& Info() const noexcept override { return actual_->Info(); }
    sirius::base::Expected<KernelHandle> LoadKernel(std::span<const std::uint32_t> code) override {
        for (auto& module : modules)
            if (code.size() == module.original.size() && std::equal(code.begin(), code.end(), module.original.begin())) {
                if (module.loads != 0)
                    return sirius::base::Fail(sirius::base::ErrorDomain::kKernel, "candidate substitution", "duplicate source module");
                WriteWords((directory_/("loaded-original-" + module.kind + ".spv")).c_str(), code);
                const auto result = actual_->LoadKernel(module.candidate);
                if (result) { ++module.loads; module.handle = *result; }
                return result;
            }
        return sirius::base::Fail(sirius::base::ErrorDomain::kKernel, "candidate substitution", "unexpected embedded source module");
    }
    sirius::base::Expected<BufferHandle> CreateBuffer(std::uint64_t bytes, BufferUsage usage) override {
        const auto result = actual_->CreateBuffer(bytes, usage);
        if (result) {
            allocations.push_back(bytes);
            std::cout << "{\"allocation\":1,\"handle\":" << result->value << ",\"bytes\":" << bytes << "}\n";
        }
        return result;
    }
    sirius::base::Expected<void> WriteBuffer(BufferHandle buffer, std::span<const std::byte> bytes) override {
        const auto result = actual_->WriteBuffer(buffer, bytes);
        if (result) { ++writes; written_bytes += bytes.size(); }
        return result;
    }
    sirius::base::Expected<void> ReadBuffer(BufferHandle buffer, std::span<std::byte> bytes) override {
        const auto result = actual_->ReadBuffer(buffer, bytes);
        if (result) {
            ++reads; read_bytes += bytes.size();
            std::ofstream file(directory_/"production-readback.bin", std::ios::binary | std::ios::trunc);
            file.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
            file.close();
            if (!file) throw std::runtime_error("passive readback recording failed");
        }
        return result;
    }
    sirius::base::Expected<void> Dispatch(KernelHandle kernel, std::span<const BufferHandle> buffers,
                                          std::uint32_t x, std::uint32_t y, std::uint32_t z,
                                          DispatchTiming* timing) override {
        ++dispatch_attempts;
        if (!modules[1].handle || modules[1].handle->value != kernel.value)
            return sirius::base::Fail(sirius::base::ErrorDomain::kDevice, "candidate transport diagnostic", "unexpected nontransport dispatch");
        std::cout << "{\"dispatch_started\":1,\"groups\":[" << x << ',' << y << ',' << z
                  << "],\"binding_count\":" << buffers.size() << "}\n";
        const auto result = actual_->Dispatch(kernel, buffers, x, y, z, timing);
        if (result) ++completed_dispatches;
        return result;
    }
    sirius::base::Expected<void> SetBufferAllocationLimit(std::uint64_t bytes) override {
        return actual_->SetBufferAllocationLimit(bytes);
    }
    std::uint64_t BufferAllocationBytes() const noexcept override { return actual_->BufferAllocationBytes(); }
    std::vector<Module> modules;
    std::vector<std::uint64_t> allocations;
    std::uint64_t writes = 0, reads = 0, written_bytes = 0, read_bytes = 0;
    std::uint64_t dispatch_attempts = 0, completed_dispatches = 0;
  private:
    std::unique_ptr<ComputeDevice> actual_;
    std::filesystem::path directory_;
};

::testing::AssertionResult StepAgrees(const RetainedStepOutput& output,
                                      const sirius::test::retained_transport::Case& fixture) {
    if (!output.valid || output.stages != 7) return ::testing::AssertionFailure();
    const std::array records{&output.fifth, &output.fourth, &output.increment, &output.error};
    for (std::size_t record = 0; record < records.size(); ++record)
        for (std::size_t group = 0; group < 10; ++group) {
            long double scale = 0, error = 0;
            for (std::size_t axis = 0; axis < 4; ++axis) {
                const auto i = group * 4 + axis, index = record * 40 + i;
                const auto& value = (*records[record])[i];
                if (!value.IsRepresented()) return ::testing::AssertionFailure();
                // Preserve the third device term even on platforms where
                // long double has only binary64 precision. The frozen oracle
                // supplies an independent twofold binary64 conversion.
                const auto reference = fixture.wide_reference[index];
                const sirius::core::Twofold exact(reference.high, reference.low);
                const auto observed = sirius::core::Twofold(value.high) +
                                      sirius::core::Twofold(value.low) +
                                      sirius::core::Twofold(value.tail);
                const double difference = std::abs((observed - exact).Rounded());
                constexpr double epsilon = std::numeric_limits<double>::epsilon();
                const double rounding =
                    128 * epsilon * epsilon * (std::abs(reference.high) + std::abs(value.high));
                if (difference > value.radius + fixture.precision_gap[index] + rounding)
                    return ::testing::AssertionFailure()
                           << "step component " << index << " difference=" << difference
                           << " radius=" << value.radius;
                const long double center =
                    (static_cast<long double>(value.high) + value.low) + value.tail;
                error = std::max(error, std::abs(center - fixture.reference[index]));
                // The small embedded difference is checked on the increment's
                // scale. Its finite arithmetic enclosure remains checked above.
                const auto scale_index = (record == 3 ? 2 : record) * 40 + i;
                scale = std::max(scale, std::abs(fixture.reference[scale_index]));
            }
            if (error > 1e-11L * scale) return ::testing::AssertionFailure();
        }
    return ::testing::AssertionSuccess();
}

class TransportPrototype : public ::testing::TestWithParam<bool> {
  protected:
    void SetUp() override {
        const bool wide = GetParam();
        mode_directory = evidence_root/(wide ? "fp64-products" : "fp32-products");
        std::filesystem::create_directories(mode_directory);
        const auto inventory = EnumerateVulkanDevices();
        ASSERT_TRUE(inventory) << inventory.error().Description();
        ASSERT_FALSE(inventory->empty());
        const auto index = ResolveVulkanDeviceIndex(*inventory);
        ASSERT_TRUE(index) << index.error().Description();
        auto opened = CreateVulkanDevice(*index);
        ASSERT_TRUE(opened) << opened.error().Description();
        const auto& info = (*opened)->Info();
        ASSERT_EQ(info.kind, DeviceKind::kSoftware);
        ASSERT_NE(info.name.find("llvmpipe"), std::string::npos);
        ASSERT_TRUE(info.supports_fp64);
        ASSERT_TRUE(info.rounds_fp64_to_nearest);
        ASSERT_TRUE(RetainedUsesPortableArithmetic(info));
        void* driver = dlopen("libvulkan_lvp.so", RTLD_NOW | RTLD_NOLOAD);
        Dl_info resident{};
        ASSERT_NE(driver, nullptr);
        ASSERT_NE(dladdr(dlsym(driver, "vk_icdGetInstanceProcAddr"), &resident), 0);
        ASSERT_NE(resident.dli_fname, nullptr);
        const std::string driver_path = resident.dli_fname;
        dlclose(driver);
        std::cout << "{\"device_identity\":1,\"product_mode\":" << std::quoted(wide ? "fp64" : "fp32")
                  << ",\"device\":" << std::quoted(info.name) << ",\"driver\":" << std::quoted(info.driver_name)
                  << ",\"driver_info\":" << std::quoted(info.driver_info)
                  << ",\"resident_driver_path\":" << std::quoted(driver_path)
                  << ",\"vendor_id\":" << info.vendor_id << ",\"device_id\":" << info.device_id
                  << ",\"driver_id\":" << info.driver_id << ",\"api_version\":" << info.api_version
                  << ",\"denorm32\":" << info.preserves_fp32_denormals
                  << ",\"RTE32\":" << info.rounds_fp32_to_nearest
                  << ",\"fp64\":" << info.supports_fp64
                  << ",\"RTE64\":" << info.rounds_fp64_to_nearest << "}\n";
        device = std::make_unique<ModuleSubstitutionDevice>(std::move(*opened), wide, mode_directory);
        ASSERT_TRUE(device->SetBufferAllocationLimit(8 * 1024 * 1024));
        auto created = RetainedCompute::Create(*device, 24, wide);
        ASSERT_TRUE(created) << created.error().Description();
        compute = std::move(*created);
        ASSERT_EQ(device->modules.size(), 6U);
        for (const auto& module : device->modules) ASSERT_EQ(module.loads, 1U);
        ASSERT_EQ(device->allocations.size(), 12U);
        RecordProperty("product_mode", wide ? "fp64" : "fp32");
        RecordProperty("scope", "finite scalar-assist transport feasibility; no production or full-stage qualification");
        RecordProperty("fixture_rows", "15");
        RecordProperty("scientific_values_per_row", "160");
        RecordProperty("negative_controls_per_mode", "28");
    }
    void TearDown() override {
        if (!compute || !device) return;
        const auto stats = compute->Statistics();
        const auto& transport = stats[1];
        std::cout << std::setprecision(17)
                  << "{\"mode_terminal\":1,\"product_mode\":" << std::quoted(GetParam() ? "fp64" : "fp32")
                  << ",\"fixture_contract_pass\":" << (!HasFailure())
                  << ",\"completed_dispatches\":" << device->completed_dispatches
                  << ",\"dispatch_attempts\":" << device->dispatch_attempts
                  << ",\"allocation_count\":" << device->allocations.size()
                  << ",\"resident_bytes\":" << device->BufferAllocationBytes()
                  << ",\"required_buffer_bytes\":" << RetainedCompute::RequiredBufferBytes(24)
                  << ",\"writes\":" << device->writes << ",\"written_bytes\":" << device->written_bytes
                  << ",\"reads\":" << device->reads << ",\"read_bytes\":" << device->read_bytes
                  << ",\"host_timing_ms\":{\"pipeline\":" << transport.pipeline_setup_ms
                  << ",\"command\":" << transport.command_setup_ms << ",\"submit_wait\":" << transport.submit_wait_ms
                  << ",\"cleanup\":" << transport.cleanup_ms << ",\"dispatch_total\":" << transport.dispatch_total_ms
                  << ",\"write\":" << transport.write_buffer_ms << ",\"read\":" << transport.read_buffer_ms
                  << "},\"stage_submissions\":[";
        for (std::size_t i = 0; i < stats.size(); ++i) std::cout << (i ? "," : "") << stats[i].submissions;
        std::cout << "]}\n";
        RecordProperty("actual_resident_bytes", std::to_string(device->BufferAllocationBytes()));
        RecordProperty("completed_dispatches", std::to_string(device->completed_dispatches));
        RecordProperty("transport_submit_wait_ms", std::to_string(transport.submit_wait_ms));
        compute.reset();
        device.reset();
    }
    std::filesystem::path mode_directory;
    std::unique_ptr<ModuleSubstitutionDevice> device;
    std::unique_ptr<RetainedCompute> compute;
};

TEST_P(TransportPrototype, JointRkStagesRetainCriticalIncrementsAndEmbeddedError) {

#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    const auto allocation = device->BufferAllocationBytes();
    std::vector<RetainedStepInput> inputs;
    for (const auto& fixture : sirius::test::retained_transport::cases) {
        inputs.push_back(std::bit_cast<RetainedStepInput>(fixture.input));
    }
    const auto outputs = compute->Step(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    ASSERT_EQ(outputs->size(), inputs.size());
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        SCOPED_TRACE(sirius::test::retained_transport::cases[row].name);
        ASSERT_TRUE(StepAgrees((*outputs)[row], sirius::test::retained_transport::cases[row]));
        if (row + 1 < inputs.size()) {
            auto two_term = (*outputs)[row];
            for (auto* record :
                 {&two_term.fifth, &two_term.fourth, &two_term.increment, &two_term.error})
                for (auto& value : *record) value.tail = 0;
            EXPECT_FALSE(StepAgrees(two_term, sirius::test::retained_transport::cases[row]));
            auto narrowed = (*outputs)[row];
            for (auto* record :
                 {&narrowed.fifth, &narrowed.fourth, &narrowed.increment, &narrowed.error})
                for (auto& value : *record) value.low = 0;
            EXPECT_FALSE(StepAgrees(narrowed, sirius::test::retained_transport::cases[row]));
        }
    }
    EXPECT_EQ(device->BufferAllocationBytes(), allocation);
    const auto& flat = outputs->back();
    for (std::size_t i = 0; i < 40; ++i) {
        EXPECT_EQ(flat.error[i].Center(), 0);
        EXPECT_EQ(flat.fifth[i].Center(), flat.fourth[i].Center());
    }
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif

    ASSERT_EQ(device->dispatch_attempts, 1U);
    ASSERT_EQ(device->completed_dispatches, 1U);
    ASSERT_EQ(device->writes, 1U);
    ASSERT_EQ(device->reads, 1U);
    std::vector<std::uint32_t> words{0x50545352u, 1u, static_cast<std::uint32_t>(outputs->size())};
    for (const auto& output : *outputs) {
        words.push_back(output.valid ? 1u : 0u);
        words.push_back(output.stages);
        for (const auto* record : {&output.fifth, &output.fourth, &output.increment, &output.error})
            for (const auto& value : *record) {
                const auto limbs = std::bit_cast<std::array<std::uint32_t, 5>>(value);
                words.insert(words.end(), limbs.begin(), limbs.end());
            }
    }
    WriteWords((mode_directory/"scientific-output.bin").c_str(), words);
}
INSTANTIATE_TEST_SUITE_P(BothProductModes, TransportPrototype, ::testing::Values(false, true),
                         [](const auto& info) { return info.param ? "Fp64Products" : "Fp32Products"; });
}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) return 2;
    module_root = argv[1]; evidence_root = argv[2];
    for (int i = 1; i + 2 < argc; ++i) argv[i] = argv[i + 2];
    argc -= 2;
    argv[argc] = nullptr;
    std::cout << std::unitbuf;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
