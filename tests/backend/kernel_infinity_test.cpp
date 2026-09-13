// Actual bounded device continuation against independent flat/axis/Carter
// oracles. Each dispatch performs at most one adaptive attempt per case.
// Accuracy budgets were declared before the first device comparison; they
// constrain this mathematical helper, not the whole-detector error budget.
#include "sirius/backend/device.h"
#include "sirius/render/dispatch_governor.h"

#include <gtest/gtest.h>

#include "../support/kerr_infinity_oracle.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <span>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace {
using namespace sirius::backend;
std::vector<std::uint32_t> Spirv(const std::string& path) {
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    if (!file) return {};
    const auto bytes = file.tellg();
    if (bytes <= 0 || static_cast<std::size_t>(bytes) % 4 != 0) return {};
    std::vector<std::uint32_t> words(static_cast<std::size_t>(bytes) / 4);
    file.seekg(0);
    file.read(reinterpret_cast<char*>(words.data()), bytes);
    return file ? words : std::vector<std::uint32_t>{};
}

constexpr std::uint64_t kInfinityAllocationLimit = 1024 * 1024;

std::string InfinityNumber(double value) {
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return stream.str();
}

void RecordInfinityDevice(const ComputeDevice& device, std::size_t index,
                          const std::string& prefix) {
    const auto& info = device.Info();
    ::testing::Test::RecordProperty(prefix + "_device_index", std::to_string(index));
    ::testing::Test::RecordProperty(prefix + "_device_name", info.name);
    ::testing::Test::RecordProperty(prefix + "_device_kind", ToString(info.kind));
    ::testing::Test::RecordProperty(prefix + "_vendor_id", std::to_string(info.vendor_id));
    ::testing::Test::RecordProperty(prefix + "_device_id", std::to_string(info.device_id));
    ::testing::Test::RecordProperty(prefix + "_driver_id", std::to_string(info.driver_id));
    ::testing::Test::RecordProperty(prefix + "_driver_name", info.driver_name);
    ::testing::Test::RecordProperty(prefix + "_driver_info", info.driver_info);
    ::testing::Test::RecordProperty(prefix + "_supports_fp64",
                                    info.supports_fp64 ? "true" : "false");
}

void RecordInfinityAllocation(const ComputeDevice& device, const std::string& prefix,
                              std::uint64_t payload_bytes) {
    const auto allocated = device.BufferAllocationBytes();
    ::testing::Test::RecordProperty(prefix + "_payload_bytes", std::to_string(payload_bytes));
    ::testing::Test::RecordProperty(prefix + "_allocated_bytes", std::to_string(allocated));
    ::testing::Test::RecordProperty(prefix + "_allocation_limit_bytes",
                                    std::to_string(kInfinityAllocationLimit));
    // Actual allocations include any adapter-required buffer padding.
    ASSERT_GE(allocated, payload_bytes);
    ASSERT_LE(allocated, kInfinityAllocationLimit);
}

void DispatchInfinityProbe(ComputeDevice& device, KernelHandle kernel,
                           std::span<const BufferHandle> buffers, std::uint32_t count,
                           const std::string& prefix, std::size_t& submissions) {
    DispatchTiming timing;
    const auto dispatched = device.Dispatch(kernel, buffers, count, 1, 1, &timing);
    ASSERT_TRUE(dispatched.has_value()) << dispatched.error().Description();
    const std::string key = prefix + "_dispatch_" + std::to_string(submissions++);
    ::testing::Test::RecordProperty(prefix + "_actual_submissions", std::to_string(submissions));
    ::testing::Test::RecordProperty(key + "_submit_wait_ms", InfinityNumber(timing.submit_wait_ms));
    ::testing::Test::RecordProperty(key + "_pipeline_setup_ms",
                                    InfinityNumber(timing.pipeline_setup_ms));
    ::testing::Test::RecordProperty(key + "_total_ms", InfinityNumber(timing.total_ms));
    ASSERT_TRUE(std::isfinite(timing.submit_wait_ms));
    ASSERT_GT(timing.submit_wait_ms, 0.0);
    ASSERT_LE(timing.submit_wait_ms, sirius::render::kDispatchStopMs);
    ASSERT_TRUE(std::isfinite(timing.pipeline_setup_ms));
    ASSERT_GE(timing.pipeline_setup_ms, 0.0);
    ASSERT_TRUE(std::isfinite(timing.total_ms));
    ASSERT_GE(timing.total_ms, timing.pipeline_setup_ms);
    ASSERT_GE(timing.total_ms, timing.submit_wait_ms);
}

template <class Real>
void RunInfinityProbe(const std::string& artifact) {
    SCOPED_TRACE(artifact);
#ifdef SIRIUS_KERNEL_DIR
    const std::string kernel_dir = SIRIUS_KERNEL_DIR;
#else
    const std::string kernel_dir;
    GTEST_SKIP() << "Vulkan kernels are not configured";
#endif
    // Declared before the first device run. These are smooth-case diagnostic
    // budgets; the final whole-detector accuracy contract is not established.
    constexpr bool wide = std::is_same_v<Real, double>;
    constexpr double direction_budget = wide ? 1e-9 : 2e-5;
    constexpr double jacobian_absolute_budget = wide ? 5e-8 : 3e-5;
    constexpr double jacobian_relative_budget = wide ? 0.0 : 2e-4;
    constexpr double frequency_relative_budget = wide ? 1e-10 : 1e-5;
    constexpr double derivative_absolute_budget = wide ? 5e-8 : 3e-5;
    constexpr Real absolute_tolerance = Real(wide ? 1e-12 : 2e-7);
    constexpr Real relative_tolerance = Real(wide ? 1e-11 : 2e-6);
    constexpr Real input_tolerance =
        Real(wide ? 1e-8 : 32.0 * std::numeric_limits<float>::epsilon());
    const auto cases = sirius::test::kerr_infinity_oracle::Cases();
    ASSERT_EQ(cases.size(), 29u);
    const auto inventory = EnumerateVulkanDevices();
    if (!inventory.has_value() || inventory->empty()) GTEST_SKIP() << "no Vulkan device available";
    const auto index = ResolveVulkanDeviceIndex(*inventory);
    ASSERT_TRUE(index.has_value()) << index.error().Description();
    auto opened = CreateVulkanDevice(*index);
    ASSERT_TRUE(opened.has_value()) << opened.error().Description();
    auto device = std::move(*opened);
    // Artifact-specific keys retain all three rungs in the combined Finish test.
    const std::string prefix = artifact.substr(0, artifact.rfind('.'));
    RecordInfinityDevice(*device, *index, prefix);
    std::size_t submissions = 0;
    const auto& info = device->Info();
    const auto words = Spirv(kernel_dir + "/" + artifact);
    ASSERT_FALSE(words.empty()) << artifact;
    if constexpr (wide) {
        if (!info.supports_fp64) {
            const auto refused = device->LoadKernel(words);
            ASSERT_FALSE(refused.has_value());
            EXPECT_EQ(refused.error().domain(), sirius::base::ErrorDomain::kKernel);
            EXPECT_NE(refused.error().detail().find("shaderFloat64"), std::string::npos);
            ::testing::Test::RecordProperty("fp64_evidence", "unsupported_kernel_declined");
            return;
        }
    }
    auto kernel = device->LoadKernel(words);
    ASSERT_TRUE(kernel.has_value()) << kernel.error().Description();
    std::vector<Real> input(cases.size() * 32), output(cases.size() * 20);
    std::vector<std::byte> opaque(cases.size() * (wide ? 640 : 320));
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const auto& c = cases[i];
        auto* in = input.data() + i * 32;
        in[0] = Real(c.mass);
        in[1] = Real(c.spin);
        in[2] = 0;
        in[3] = Real(c.seed);
        in[4] = absolute_tolerance;
        in[5] = relative_tolerance;
        in[6] = input_tolerance;
        in[7] = Real(c.maximum_attempts);
        for (unsigned j = 0; j < 4; ++j) {
            in[8 + j] = Real(c.x(j));
            in[12 + j] = Real(c.k(j));
            for (unsigned col = 0; col < 2; ++col) {
                in[16 + 4 * col + j] = Real(c.X[col](j));
                in[24 + 4 * col + j] = Real(c.V[col](j));
            }
        }
    }
    ASSERT_TRUE(device->SetBufferAllocationLimit(kInfinityAllocationLimit).has_value());
    std::array<BufferHandle, 3> buffers;
    const std::array<std::size_t, 3> sizes{input.size() * sizeof(Real), opaque.size(),
                                           output.size() * sizeof(Real)};
    for (unsigned j = 0; j < 3; ++j) {
        auto created = device->CreateBuffer(sizes[j], BufferUsage::kStorage);
        ASSERT_TRUE(created.has_value()) << created.error().Description();
        buffers[j] = *created;
    }
    ASSERT_NO_FATAL_FAILURE(
        RecordInfinityAllocation(*device, prefix, sizes[0] + sizes[1] + sizes[2]));
    auto written = device->WriteBuffer(buffers[0], std::as_bytes(std::span(input)));
    ASSERT_TRUE(written.has_value()) << written.error().Description();
    written = device->WriteBuffer(buffers[1], opaque);
    ASSERT_TRUE(written.has_value()) << written.error().Description();
    unsigned max_attempts = 0;
    for (const auto& c : cases) max_attempts = std::max(max_attempts, c.maximum_attempts);
    bool finished = false;
    for (unsigned dispatch = 0; dispatch <= max_attempts + 1; ++dispatch) {
        ASSERT_NO_FATAL_FAILURE(DispatchInfinityProbe(*device, *kernel, buffers,
                                                      static_cast<std::uint32_t>(cases.size()),
                                                      prefix, submissions));
        const auto read = device->ReadBuffer(buffers[2], std::as_writable_bytes(std::span(output)));
        ASSERT_TRUE(read.has_value()) << read.error().Description();
        finished = true;
        for (std::size_t i = 0; i < cases.size(); ++i) {
            const auto* out = output.data() + i * 20;
            for (unsigned j = 0; j < 20; ++j)
                ASSERT_TRUE(std::isfinite(out[j]))
                    << artifact << " " << cases[i].name << " row scalar " << j;
            SCOPED_TRACE(cases[i].name);
            if (dispatch == 0) {
                EXPECT_EQ(out[0], cases[i].expected_status == 2   ? 2
                                  : cases[i].expected_status == 3 ? 3
                                                                  : 0);
                EXPECT_EQ(out[1], 0);
                EXPECT_EQ(out[2], 0);
                EXPECT_EQ(out[3], 0);
                EXPECT_EQ(out[16], 0);
            }
            EXPECT_GE(out[1], out[2]);
            EXPECT_LE(out[1], Real(cases[i].maximum_attempts));
            EXPECT_GE(out[3], 0);
            EXPECT_LE(out[3], 1);
            if (out[0] == 0) finished = false;
        }
        if (finished) break;
        if (dispatch == 0) {
            for (std::size_t i = 0; i < cases.size(); ++i) input[i * 32 + 2] = 1;
            written = device->WriteBuffer(buffers[0], std::as_bytes(std::span(input)));
            ASSERT_TRUE(written.has_value()) << written.error().Description();
        }
    }
    ASSERT_TRUE(finished);
    for (std::size_t i = 0; i < cases.size(); ++i) {
        const auto& c = cases[i];
        const auto* out = output.data() + i * 20;
        SCOPED_TRACE(c.name);
        EXPECT_EQ(out[0], Real(c.expected_status));
        EXPECT_EQ(out[16], c.expected_status == 1 ? 1 : 0);
        if (c.expected_status != 1 || out[0] != 1) continue;
        EXPECT_EQ(out[3], 1);
        EXPECT_LE(out[15], 1);
        EXPECT_GE(out[15], 0);
        for (unsigned j = 0; j < 3; ++j) EXPECT_NEAR(out[4 + j], c.direction[j], direction_budget);
        for (unsigned row = 0; row < 2; ++row)
            for (unsigned col = 0; col < 2; ++col) {
                const double expected = c.jacobian[row][col];
                EXPECT_NEAR(
                    out[8 + 2 * row + col], expected,
                    jacobian_absolute_budget + jacobian_relative_budget * std::abs(expected));
            }
        const double determinant =
            c.jacobian[0][0] * c.jacobian[1][1] - c.jacobian[0][1] * c.jacobian[1][0];
        double determinant_budget = 0;
        for (unsigned row = 0; row < 2; ++row)
            for (unsigned col = 0; col < 2; ++col)
                determinant_budget +=
                    (jacobian_absolute_budget +
                     jacobian_relative_budget * std::abs(c.jacobian[row][col])) *
                    (std::abs(c.jacobian[1 - row][1 - col]) + jacobian_absolute_budget);
        EXPECT_NEAR(out[7], determinant, determinant_budget);
        EXPECT_NEAR(out[12], c.frequency, frequency_relative_budget * std::abs(c.frequency));
        for (unsigned col = 0; col < 2; ++col)
            EXPECT_NEAR(out[13 + col], c.frequency_derivative[col],
                        derivative_absolute_budget * std::max(1.0, std::abs(c.frequency)));
    }
    if constexpr (wide)
        ::testing::Test::RecordProperty("fp64_evidence", "device_infinity_continuations_executed");
    // Every completed or declined continuation must be sticky under attempts.
    const auto terminal = output;
    for (std::size_t i = 0; i < cases.size(); ++i) input[i * 32 + 2] = 1;
    written = device->WriteBuffer(buffers[0], std::as_bytes(std::span(input)));
    ASSERT_TRUE(written.has_value());
    ASSERT_NO_FATAL_FAILURE(DispatchInfinityProbe(
        *device, *kernel, buffers, static_cast<std::uint32_t>(cases.size()), prefix, submissions));
    auto read = device->ReadBuffer(buffers[2], std::as_writable_bytes(std::span(output)));
    ASSERT_TRUE(read.has_value());
    EXPECT_EQ(output, terminal);
    // Explicit initialization replaces previous terminal states with new ones.
    for (std::size_t i = 0; i < cases.size(); ++i) {
        input[i * 32 + 2] = 0;
        input[i * 32 + 3] = 0;
    }
    written = device->WriteBuffer(buffers[0], std::as_bytes(std::span(input)));
    ASSERT_TRUE(written.has_value());
    ASSERT_NO_FATAL_FAILURE(DispatchInfinityProbe(
        *device, *kernel, buffers, static_cast<std::uint32_t>(cases.size()), prefix, submissions));
    read = device->ReadBuffer(buffers[2], std::as_writable_bytes(std::span(output)));
    ASSERT_TRUE(read.has_value());
    for (std::size_t i = 0; i < cases.size(); ++i) {
        EXPECT_EQ(output[i * 20], 2);
        EXPECT_EQ(output[i * 20 + 1], 0);
        EXPECT_EQ(output[i * 20 + 16], 0);
    }
}
TEST(KernelInfinityDevice, Fp32IndependentDirectionsVariationsAndTerminalStates) {
    ASSERT_NO_FATAL_FAILURE(RunInfinityProbe<float>("infinity_probe.spv"));
}
TEST(KernelInfinityDevice, CompensatedIndependentDirectionsVariationsAndTerminalStates) {
    ASSERT_NO_FATAL_FAILURE(RunInfinityProbe<float>("infinity_probe_fp32comp.spv"));
}
TEST(KernelInfinityDevice, Fp64IndependentDirectionsVariationsAndTerminalStates) {
    ASSERT_NO_FATAL_FAILURE(RunInfinityProbe<double>("infinity_probe_fp64.spv"));
}

template <class Real>
void RunFinishProbe(const std::string& artifact) {
    SCOPED_TRACE(artifact);
#ifdef SIRIUS_KERNEL_DIR
    const std::string kernel_dir = SIRIUS_KERNEL_DIR;
#else
    const std::string kernel_dir;
    GTEST_SKIP() << "Vulkan kernels are not configured";
#endif
    const auto inventory = EnumerateVulkanDevices();
    if (!inventory.has_value() || inventory->empty()) GTEST_SKIP() << "no Vulkan device available";
    const auto index = ResolveVulkanDeviceIndex(*inventory);
    ASSERT_TRUE(index.has_value()) << index.error().Description();
    auto opened = CreateVulkanDevice(*index);
    ASSERT_TRUE(opened.has_value()) << opened.error().Description();
    auto device = std::move(*opened);
    // Artifact-specific keys retain all three rungs in the combined Finish test.
    const std::string prefix = artifact.substr(0, artifact.rfind('.'));
    RecordInfinityDevice(*device, *index, prefix);
    std::size_t submissions = 0;
    const auto words = Spirv(kernel_dir + "/" + artifact);
    ASSERT_FALSE(words.empty()) << artifact;
    constexpr bool wide = std::is_same_v<Real, double>;
    if constexpr (wide) {
        if (!device->Info().supports_fp64) {
            const auto refused = device->LoadKernel(words);
            ASSERT_FALSE(refused.has_value());
            EXPECT_EQ(refused.error().domain(), sirius::base::ErrorDomain::kKernel);
            EXPECT_NE(refused.error().detail().find("shaderFloat64"), std::string::npos);
            ::testing::Test::RecordProperty("fp64_evidence", "unsupported_kernel_declined");
            return;
        }
    }
    auto kernel = device->LoadKernel(words);
    ASSERT_TRUE(kernel.has_value()) << kernel.error().Description();
    constexpr std::size_t count = 8;
    const std::array<const char*, count> names{
        "identity_success",  "determinant_overflow", "running",    "invalid_input",
        "radial_unresolved", "work_limit",           "arithmetic", "complete_partial"};
    const std::array<unsigned, count> statuses{1, 1, 0, 2, 3, 4, 5, 1};
    std::vector<Real> input(count * 32, Real(0));
    std::vector<Real> output(count * 20, Real(-123));
    std::vector<std::byte> opaque(count * (wide ? 640 : 320));
    for (std::size_t i = 0; i < count; ++i) {
        auto* in = input.data() + i * 32;
        in[0] = Real(statuses[i]);
        in[1] = i == 7 ? Real(0.5) : Real(1);
        in[2] = Real(2);  // Direct Finish action, without initialization or stepping.
        in[4] = Real(1);  // Frequency and its two angular derivatives.
        in[5] = Real(0.25);
        in[6] = Real(-0.5);
        in[7] = Real(0.5);  // Accepted local error ratio.
        // n=(0,0,1), with independently chosen transverse derivatives. Every
        // input is finite; only the determinant product overflows in case 1.
        const Real amplitude = i == 1 ? Real(wide ? 1e200 : 1e20) : Real(1);
        ASSERT_TRUE(std::isfinite(amplitude));
        in[9] = amplitude;
        in[14] = amplitude;
        in[16] = Real(1);
    }
    ASSERT_TRUE(device->SetBufferAllocationLimit(kInfinityAllocationLimit).has_value());
    std::array<BufferHandle, 3> buffers;
    const std::array<std::size_t, 3> sizes{input.size() * sizeof(Real), opaque.size(),
                                           output.size() * sizeof(Real)};
    for (unsigned j = 0; j < buffers.size(); ++j) {
        auto created = device->CreateBuffer(sizes[j], BufferUsage::kStorage);
        ASSERT_TRUE(created.has_value()) << created.error().Description();
        buffers[j] = *created;
    }
    ASSERT_NO_FATAL_FAILURE(
        RecordInfinityAllocation(*device, prefix, sizes[0] + sizes[1] + sizes[2]));
    auto written = device->WriteBuffer(buffers[0], std::as_bytes(std::span(input)));
    ASSERT_TRUE(written.has_value()) << written.error().Description();
    written = device->WriteBuffer(buffers[1], opaque);
    ASSERT_TRUE(written.has_value()) << written.error().Description();
    written = device->WriteBuffer(buffers[2], std::as_bytes(std::span(output)));
    ASSERT_TRUE(written.has_value()) << written.error().Description();
    ASSERT_NO_FATAL_FAILURE(
        DispatchInfinityProbe(*device, *kernel, buffers, count, prefix, submissions));
    const auto read = device->ReadBuffer(buffers[2], std::as_writable_bytes(std::span(output)));
    ASSERT_TRUE(read.has_value()) << read.error().Description();
    for (std::size_t i = 0; i < count; ++i) {
        SCOPED_TRACE(names[i]);
        const auto* out = output.data() + i * 20;
        for (unsigned j = 0; j < 20; ++j) EXPECT_TRUE(std::isfinite(out[j])) << j;
        EXPECT_EQ(out[0], Real(statuses[i]));
        EXPECT_EQ(out[1], Real(9));
        EXPECT_EQ(out[2], Real(7));
        EXPECT_EQ(out[3], input[i * 32 + 1]);
        EXPECT_EQ(out[16], i == 0 ? Real(1) : Real(0));
        EXPECT_EQ(out[19], Real(1));  // Raw state frequency remains separate from availability.
        if (i == 0) {
            const std::array<Real, 12> expected{0, 0, 1, 1,          1,          0,
                                                0, 1, 1, Real(0.25), Real(-0.5), Real(0.5)};
            for (unsigned j = 0; j < expected.size(); ++j)
                EXPECT_EQ(out[4 + j], expected[j]) << "result scalar " << j;
            EXPECT_EQ(out[17], Real(7));
            EXPECT_EQ(out[18], Real(9));
        } else {
            // Complete status alone does not authorize any result publication:
            // an unfinished fraction or derived overflow must leave all fields zero.
            for (unsigned j = 4; j < 16; ++j) EXPECT_EQ(out[j], Real(0)) << "result scalar " << j;
            EXPECT_EQ(out[17], Real(0));
            EXPECT_EQ(out[18], Real(0));
        }
    }
    if constexpr (wide)
        ::testing::Test::RecordProperty("fp64_evidence", "device_finish_boundary_executed");
}

TEST(KernelInfinityDevice, FinishDeclinesWithoutPublishingPartialResults) {
    ASSERT_NO_FATAL_FAILURE(RunFinishProbe<float>("infinity_probe.spv"));
    ASSERT_NO_FATAL_FAILURE(RunFinishProbe<float>("infinity_probe_fp32comp.spv"));
    ASSERT_NO_FATAL_FAILURE(RunFinishProbe<double>("infinity_probe_fp64.spv"));
}
}  // namespace
