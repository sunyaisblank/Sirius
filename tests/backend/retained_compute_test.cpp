#include "sirius/backend/retained_compute.h"

#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/backend/retained_integrator.h"
#include "sirius/backend/retained_trace_executor.h"
#include "sirius/core/twofold.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstring>
#include <functional>
#include <future>
#include <iomanip>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

#if defined(SIRIUS_HAS_RETAINED_COMPUTE) && defined(SIRIUS_RETAINED_CAMERA_TEST_DIR)
#define SIRIUS_RETAINED_TESTS_AVAILABLE 1
#include "sirius/backend/vulkan/vulkan_device.h"

#include "program_fixture.h"
#include "support/cpu_critical/transport_reference.h"
#include "support/retained_camera/continuous_reference.h"
#include "support/retained_camera/ray_reference.h"
#include "support/retained_transport/dense_reference.h"
#include "support/retained_transport/endpoint_reference.h"
#include "support/retained_transport/reference_cases.h"
#endif

namespace {
using namespace sirius::backend;

// Arithmetic admission must finish before any external shader or buffer work.
// The sentinel gives that boundary an observable positive path without a device.
class AdmissionDevice final : public ComputeDevice {
  public:
    DeviceInfo info;
    unsigned kernel_calls = 0, buffer_calls = 0;
    bool kernel_has_float = false;
    const DeviceInfo& Info() const noexcept override { return info; }
    sirius::base::Expected<KernelHandle> LoadKernel(std::span<const std::uint32_t> code) override {
        ++kernel_calls;
        // Inspect the actual selected module: OpTypeFloat is opcode 22.
        for (std::size_t offset = 5; offset < code.size();) {
            const auto words = code[offset] >> 16;
            if (words == 0 || words > code.size() - offset) break;
            kernel_has_float |= (code[offset] & 0xffffU) == 22U;
            offset += words;
        }
        return sirius::base::Fail(sirius::base::ErrorDomain::kKernel, "admission sentinel",
                                  "arithmetic admitted; external kernel loading reached");
    }
    sirius::base::Expected<BufferHandle> CreateBuffer(std::uint64_t, BufferUsage) override {
        ++buffer_calls;
        return sirius::base::Fail(sirius::base::ErrorDomain::kDevice, "admission sentinel",
                                  "unexpected buffer allocation");
    }
    sirius::base::Expected<void> WriteBuffer(BufferHandle, std::span<const std::byte>) override {
        return Unused();
    }
    sirius::base::Expected<void> ReadBuffer(BufferHandle, std::span<std::byte>) override {
        return Unused();
    }
    sirius::base::Expected<void> Dispatch(KernelHandle, std::span<const BufferHandle>,
                                          std::uint32_t, std::uint32_t, std::uint32_t,
                                          DispatchTiming*) override {
        return Unused();
    }
    sirius::base::Expected<void> SetBufferAllocationLimit(std::uint64_t) override {
        return Unused();
    }
    std::uint64_t BufferAllocationBytes() const noexcept override { return 0; }

  private:
    sirius::base::Expected<void> Unused() {
        return sirius::base::Fail(sirius::base::ErrorDomain::kDevice, "admission sentinel",
                                  "unexpected device work");
    }
};

TEST(RetainedComputeAdmission, ArithmeticRefusalPrecedesKernelLoading) {
#ifdef SIRIUS_HAS_RETAINED_COMPUTE
    const DeviceInfo complete{.supports_fp64 = true,
                              .preserves_fp32_denormals = true,
                              .rounds_fp32_to_nearest = true,
                              .rounds_fp64_to_nearest = true};
    for (const bool wide : {false, true}) {
        // Missing native binary32 controls select actual integer modules. The
        // binary64 rung still refuses absent support before external work.
        for (unsigned missing = 0; missing < 6; ++missing) {
            SCOPED_TRACE(wide);
            SCOPED_TRACE(missing);
            AdmissionDevice device;
            device.info = complete;
            if (missing == 1) device.info.preserves_fp32_denormals = false;
            if (missing == 2) device.info.rounds_fp32_to_nearest = false;
            if (missing == 3) device.info.supports_fp64 = false;
            if (missing == 4) device.info.rounds_fp64_to_nearest = false;
            if (missing == 5) {
                device.info.preserves_fp32_denormals = false;
                device.info.rounds_fp32_to_nearest = false;
            }
            const bool admitted = !wide || (missing != 3 && missing != 4);
            const auto created = RetainedCompute::Create(device, 1, wide);
            ASSERT_FALSE(created);
            EXPECT_EQ(device.kernel_calls, admitted ? 1U : 0U);
            EXPECT_EQ(device.buffer_calls, 0U);
            if (admitted) {
                EXPECT_EQ(device.kernel_has_float, missing != 1 && missing != 2 && missing != 5);
            }
            EXPECT_EQ(created.error().domain(), admitted ? sirius::base::ErrorDomain::kKernel
                                                         : sirius::base::ErrorDomain::kDevice);
            EXPECT_EQ(created.error().operation(),
                      admitted ? "admission sentinel" : "create retained compute stages");
        }
    }
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

class RetainedComputeTest : public ::testing::Test {
  protected:
    void SetUp() override {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
        const auto inventory = EnumerateVulkanDevices();
        ASSERT_TRUE(inventory) << inventory.error().Description();
        if (inventory->empty()) GTEST_SKIP() << "No Vulkan device available";
        const auto index = ResolveVulkanDeviceIndex(*inventory);
        ASSERT_TRUE(index) << index.error().Description();
        auto opened = CreateVulkanDevice(*index);
        ASSERT_TRUE(opened) << opened.error().Description();
        device = std::move(*opened);
        ASSERT_TRUE(device->SetBufferAllocationLimit(8 * 1024 * 1024));
        auto created = RetainedCompute::Create(*device, 24);
        ASSERT_TRUE(created) << created.error().Description();
        compute = std::move(*created);
#else
        GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
    }
    std::unique_ptr<ComputeDevice> device;
    std::unique_ptr<RetainedCompute> compute;
};

#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
// Observe real Vulkan transfers; injected returned errors stop at this seam.
class TransferProbeDevice final : public ComputeDevice {
  public:
    enum class Failure { None, Write, Dispatch, Read };
    struct Allocation {
        BufferHandle handle;
        std::uint64_t bytes;
        BufferUsage usage;
    };
    struct Transfer {
        BufferHandle handle;
        std::size_t bytes;
        std::uint32_t capacity = 0;
    };
    struct Submission {
        std::array<BufferHandle, 2> buffers;
        std::size_t binding_count;
        std::uint32_t x, y, z;
    };
    explicit TransferProbeDevice(ComputeDevice& device) : device_(device) {}
    Failure fail_next = Failure::None;
    std::vector<Allocation> allocations;
    std::vector<Transfer> writes, reads;
    std::vector<Submission> submissions;
    std::uint64_t dispatch_calls = 0, forwarded_dispatches = 0;
    // Optional control-test observations. Physical outputs still come from the
    // actual device; injected timing never sleeps or claims hardware duration.
    std::function<double(BufferHandle, std::uint32_t)> submission_ms;
    std::optional<BufferHandle> capture_output;
    std::vector<std::vector<std::byte>> readbacks;

    const DeviceInfo& Info() const noexcept override { return device_.Info(); }
    sirius::base::Expected<KernelHandle> LoadKernel(std::span<const std::uint32_t> code) override {
        return device_.LoadKernel(code);
    }
    sirius::base::Expected<BufferHandle> CreateBuffer(std::uint64_t bytes,
                                                      BufferUsage usage) override {
        auto buffer = device_.CreateBuffer(bytes, usage);
        if (buffer) allocations.push_back({*buffer, bytes, usage});
        return buffer;
    }
    sirius::base::Expected<void> WriteBuffer(BufferHandle buffer,
                                             std::span<const std::byte> data) override {
        std::uint32_t capacity = 0;
        if (data.size_bytes() >= sizeof(capacity))
            std::memcpy(&capacity, data.data(), sizeof(capacity));
        writes.push_back({buffer, data.size_bytes(), capacity});
        if (fail_next == Failure::Write) return Inject("write");
        return device_.WriteBuffer(buffer, data);
    }
    sirius::base::Expected<void> ReadBuffer(BufferHandle buffer,
                                            std::span<std::byte> data) override {
        reads.push_back({buffer, data.size_bytes()});
        if (fail_next == Failure::Read) return Inject("read");
        auto status = device_.ReadBuffer(buffer, data);
        if (status && capture_output && buffer.value == capture_output->value)
            readbacks.emplace_back(data.begin(), data.end());
        return status;
    }
    sirius::base::Expected<void> Dispatch(KernelHandle kernel,
                                          std::span<const BufferHandle> buffers, std::uint32_t x,
                                          std::uint32_t y, std::uint32_t z,
                                          DispatchTiming* timing) override {
        ++dispatch_calls;
        submissions.push_back({{buffers.empty() ? BufferHandle{} : buffers[0],
                                buffers.size() < 2 ? BufferHandle{} : buffers[1]},
                               buffers.size(),
                               x,
                               y,
                               z});
        if (fail_next == Failure::Dispatch) return Inject("dispatch");
        ++forwarded_dispatches;
        auto status = device_.Dispatch(kernel, buffers, x, y, z, timing);
        if (status && timing && submission_ms) {
            timing->submit_wait_ms = submission_ms(buffers[0], x);
            timing->total_ms = timing->pipeline_setup_ms + timing->command_setup_ms +
                               timing->submit_wait_ms + timing->cleanup_ms;
        }
        return status;
    }
    sirius::base::Expected<void> SetBufferAllocationLimit(std::uint64_t bytes) override {
        return device_.SetBufferAllocationLimit(bytes);
    }
    std::uint64_t BufferAllocationBytes() const noexcept override {
        return device_.BufferAllocationBytes();
    }

  private:
    sirius::base::Expected<void> Inject(const char* operation) {
        fail_next = Failure::None;
        return sirius::base::Fail(sirius::base::ErrorDomain::kDevice,
                                  std::string("retained transfer probe ") + operation,
                                  "injected returned failure");
    }
    ComputeDevice& device_;
};

bool Encloses(const RetainedValue& value, long double expected, long double reference_gap) {
    if (!value.IsRepresented()) return false;
    const long double center = (static_cast<long double>(value.high) + value.low) + value.tail;
    const long double rounding =
        16 * std::numeric_limits<long double>::epsilon() * (std::abs(center) + std::abs(expected));
    return std::abs(center - expected) <= value.radius + reference_gap + rounding;
}

// Compare stored words, not a rounded expansion center or structure padding.
::testing::AssertionResult IntervalBitsAgree(const RetainedIntervalOutput& actual,
                                             const RetainedIntervalOutput& expected) {
    if (actual.admissible != expected.admissible || actual.failure != expected.failure ||
        actual.attempted_stages != expected.attempted_stages ||
        std::bit_cast<std::uint64_t>(actual.error_ratio) !=
            std::bit_cast<std::uint64_t>(expected.error_ratio))
        return ::testing::AssertionFailure() << "interval admission metadata differs";
    const auto values_agree = [](const auto& first, const auto& second,
                                 const std::string& name) -> ::testing::AssertionResult {
        for (std::size_t i = 0; i < first.size(); ++i) {
            const auto a = std::bit_cast<std::array<std::uint32_t, 5>>(first[i]);
            const auto b = std::bit_cast<std::array<std::uint32_t, 5>>(second[i]);
            for (std::size_t word = 0; word < a.size(); ++word)
                if (a[word] != b[word])
                    return ::testing::AssertionFailure()
                           << name << " component " << i << " word " << word << " actual=0x"
                           << std::hex << a[word] << " expected=0x" << b[word];
        }
        return ::testing::AssertionSuccess();
    };
    const std::array first{&actual.full, &actual.lower, &actual.midpoint, &actual.refined};
    const std::array second{&expected.full, &expected.lower, &expected.midpoint, &expected.refined};
    constexpr std::array names{"full", "lower", "midpoint", "refined"};
    for (std::size_t i = 0; i < first.size(); ++i) {
        if (first[i]->valid != second[i]->valid || first[i]->component != second[i]->component)
            return ::testing::AssertionFailure() << names[i] << " projection metadata differs";
        const auto phase =
            values_agree(first[i]->phase, second[i]->phase, std::string(names[i]) + " phase");
        if (!phase) return phase;
        const auto physical = values_agree(first[i]->physical, second[i]->physical,
                                           std::string(names[i]) + " physical");
        if (!physical) return physical;
    }
    const std::array increments{&actual.full_increment, &actual.lower_increment,
                                &actual.midpoint_increment, &actual.refined_increment};
    const std::array expected_increments{&expected.full_increment, &expected.lower_increment,
                                         &expected.midpoint_increment, &expected.refined_increment};
    for (std::size_t i = 0; i < increments.size(); ++i) {
        const auto result = values_agree(*increments[i], *expected_increments[i],
                                         std::string(names[i]) + " increment");
        if (!result) return result;
    }
    return ::testing::AssertionSuccess();
}

template <typename Fixture>
::testing::AssertionResult CameraAgrees(const RetainedCameraOutput& output,
                                        const Fixture& fixture) {
    if (!output.valid) return ::testing::AssertionFailure() << "invalid camera row";
    constexpr std::array<std::size_t, 9> boundaries{0, 4, 8, 24, 40, 56, 72, 88, 104};
    for (std::size_t group = 1; group < boundaries.size(); ++group) {
        long double scale = 0, error = 0;
        for (std::size_t i = boundaries[group - 1]; i < boundaries[group]; ++i) {
            if (!Encloses(output.values[i], fixture.reference[i], fixture.reference_gap[i]))
                return ::testing::AssertionFailure()
                       << "camera component " << i << " center=" << std::setprecision(20)
                       << output.values[i].Center() << " reference=" << fixture.reference[i]
                       << " radius=" << output.values[i].radius;
            const long double center =
                (static_cast<long double>(output.values[i].high) + output.values[i].low) +
                output.values[i].tail;
            scale = std::max(scale, std::abs(fixture.reference[i]));
            error = std::max(error, std::abs(center - fixture.reference[i]));
        }
        if (error > 1e-11L * scale) return ::testing::AssertionFailure();
    }
    return ::testing::AssertionSuccess();
}

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

#endif

TEST_F(RetainedComputeTest, BatchedCameraPreservesPhysicalColumnsAndRejectsInvalidRows) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    const auto initial_allocation = device->BufferAllocationBytes();
    EXPECT_FALSE(RetainedCompute::Create(*device, 0));
    EXPECT_FALSE(RetainedCompute::Create(*device, 65536));
    EXPECT_EQ(device->BufferAllocationBytes(), initial_allocation);
    TransferProbeDevice probe(*device);
    auto observed = RetainedCompute::Create(probe, compute->Capacity());
    ASSERT_TRUE(observed) << observed.error().Description();
    auto& camera_compute = **observed;
    const auto allocation = device->BufferAllocationBytes();
    EXPECT_GT(allocation, initial_allocation);
    EXPECT_LE(allocation, 8ull * 1024 * 1024);
    ASSERT_EQ(probe.allocations.size(), 12U);
    for (const auto& buffer : probe.allocations) EXPECT_EQ(buffer.usage, BufferUsage::kStorage);
    EXPECT_TRUE(probe.writes.empty());
    const auto full_input_bytes = probe.allocations[0].bytes;
    const auto output_row_bytes = probe.allocations[1].bytes / camera_compute.Capacity();
    const auto expect_transfer = [&](std::size_t rows, bool first_upload = false) {
        ASSERT_FALSE(probe.writes.empty());
        EXPECT_EQ(probe.writes.back().handle.value, probe.allocations[0].handle.value);
        EXPECT_EQ(probe.writes.back().capacity, camera_compute.Capacity());
        EXPECT_EQ(probe.writes.back().bytes,
                  first_upload ? full_input_bytes
                               : sizeof(std::uint32_t) + rows * sizeof(RetainedCameraInput));
    };
    const auto expect_read = [&](std::size_t rows) {
        ASSERT_FALSE(probe.reads.empty());
        EXPECT_EQ(probe.reads.back().handle.value, probe.allocations[1].handle.value);
        EXPECT_EQ(probe.reads.back().bytes, rows * output_row_bytes);
    };
    const auto expect_reset_timing = [](const DispatchTiming& timing) {
        EXPECT_EQ(timing.submit_wait_ms, 0);
        EXPECT_EQ(timing.pipeline_setup_ms, 0);
        EXPECT_EQ(timing.command_setup_ms, 0);
        EXPECT_EQ(timing.cleanup_ms, 0);
        EXPECT_EQ(timing.total_ms, 0);
        EXPECT_FALSE(timing.pipeline_created);
    };
    std::vector<RetainedCameraInput> inputs;
    for (const auto& fixture : sirius::test::retained_camera::kCases) {
        RetainedCameraInput input;
        for (std::size_t i = 0; i < 32; ++i)
            input.values[i] = RetainedValue::FromDouble(std::bit_cast<float>(fixture.input[i]));
        inputs.push_back(input);
    }
    const auto original_inputs = inputs;
    DispatchTiming timing{.submit_wait_ms = 1,
                          .pipeline_setup_ms = 2,
                          .command_setup_ms = 3,
                          .cleanup_ms = 4,
                          .total_ms = 10,
                          .pipeline_created = true};
    probe.fail_next = TransferProbeDevice::Failure::Write;
    auto failed = camera_compute.Camera(inputs, &timing);
    ASSERT_FALSE(failed);
    EXPECT_EQ(failed.error().operation(), "retained transfer probe write");
    ASSERT_NO_FATAL_FAILURE(expect_transfer(inputs.size(), true));
    ASSERT_NO_FATAL_FAILURE(expect_reset_timing(timing));
    EXPECT_EQ(probe.dispatch_calls, 0U);
    EXPECT_TRUE(probe.reads.empty());
    auto stats = camera_compute.Statistics()[0];
    EXPECT_EQ(stats.submissions, 0U);
    EXPECT_EQ(stats.write_buffer_bytes, 0U);
    EXPECT_EQ(stats.read_buffer_bytes, 0U);

    // A successful full write seeds the immutable program even if submission
    // then fails. The following retry must use the active prefix.
    probe.fail_next = TransferProbeDevice::Failure::Dispatch;
    failed = camera_compute.Camera(inputs, &timing);
    ASSERT_FALSE(failed);
    EXPECT_EQ(failed.error().operation(), "retained transfer probe dispatch");
    ASSERT_NO_FATAL_FAILURE(expect_transfer(inputs.size(), true));
    ASSERT_NO_FATAL_FAILURE(expect_reset_timing(timing));
    EXPECT_EQ(probe.dispatch_calls, 1U);
    EXPECT_EQ(probe.forwarded_dispatches, 0U);
    EXPECT_TRUE(probe.reads.empty());
    stats = camera_compute.Statistics()[0];
    EXPECT_EQ(stats.submissions, 0U);
    EXPECT_EQ(stats.write_buffer_bytes, full_input_bytes);
    EXPECT_EQ(stats.read_buffer_bytes, 0U);

    auto outputs = camera_compute.Camera(inputs, &timing);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    ASSERT_NO_FATAL_FAILURE(expect_transfer(inputs.size()));
    ASSERT_NO_FATAL_FAILURE(expect_read(inputs.size()));
    stats = camera_compute.Statistics()[0];
    EXPECT_EQ(stats.submissions, 1U);
    EXPECT_EQ(stats.write_buffer_bytes, full_input_bytes + sizeof(std::uint32_t) +
                                            inputs.size() * sizeof(RetainedCameraInput));
    EXPECT_EQ(stats.read_buffer_bytes, inputs.size() * output_row_bytes);
    ASSERT_EQ(outputs->size(), inputs.size());
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        SCOPED_TRACE(row);
        const auto& fixture = sirius::test::retained_camera::kCases[row];
        ASSERT_TRUE(CameraAgrees((*outputs)[row], fixture));
        auto narrowed = (*outputs)[row];
        for (auto& value : narrowed.values) value.low = 0;
        EXPECT_FALSE(CameraAgrees(narrowed, fixture));
    }
    const auto before_read_failure = camera_compute.Statistics()[0];
    const auto forwarded_before = probe.forwarded_dispatches;
    probe.fail_next = TransferProbeDevice::Failure::Read;
    failed = camera_compute.Camera(inputs, &timing);
    ASSERT_FALSE(failed);
    EXPECT_EQ(failed.error().operation(), "retained transfer probe read");
    ASSERT_NO_FATAL_FAILURE(expect_transfer(inputs.size()));
    ASSERT_NO_FATAL_FAILURE(expect_read(inputs.size()));
    EXPECT_EQ(probe.forwarded_dispatches, forwarded_before + 1);
    stats = camera_compute.Statistics()[0];
    EXPECT_EQ(stats.submissions, before_read_failure.submissions + 1);
    EXPECT_EQ(stats.submit_wait_ms, before_read_failure.submit_wait_ms + timing.submit_wait_ms);
    EXPECT_EQ(stats.write_buffer_bytes, before_read_failure.write_buffer_bytes +
                                            sizeof(std::uint32_t) +
                                            inputs.size() * sizeof(RetainedCameraInput));
    EXPECT_EQ(stats.read_buffer_bytes, before_read_failure.read_buffer_bytes);

    inputs[2].values[20].high = std::numeric_limits<float>::quiet_NaN();
    inputs[6].values[29] = RetainedValue::FromDouble(2);
    outputs = camera_compute.Camera(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    ASSERT_NO_FATAL_FAILURE(expect_transfer(inputs.size()));
    ASSERT_NO_FATAL_FAILURE(expect_read(inputs.size()));
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        if (row == 2 || row == 6) {
            EXPECT_FALSE((*outputs)[row].valid);
            for (const auto& value : (*outputs)[row].values) EXPECT_EQ(value.valid, 0U);
        } else {
            EXPECT_TRUE(CameraAgrees((*outputs)[row], sirius::test::retained_camera::kCases[row]));
        }
    }
    inputs.resize(1);
    outputs = camera_compute.Camera(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    ASSERT_NO_FATAL_FAILURE(expect_transfer(inputs.size()));
    ASSERT_NO_FATAL_FAILURE(expect_read(inputs.size()));
    ASSERT_EQ(outputs->size(), 1U);
    EXPECT_TRUE(CameraAgrees(outputs->front(), sirius::test::retained_camera::kCases.front()));
    // Grow after shrinking and replace previously invalid rows, preserving all
    // original independent physical-column expectations.
    inputs = original_inputs;
    outputs = camera_compute.Camera(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    ASSERT_NO_FATAL_FAILURE(expect_transfer(inputs.size()));
    ASSERT_NO_FATAL_FAILURE(expect_read(inputs.size()));
    ASSERT_EQ(outputs->size(), inputs.size());
    for (std::size_t row = 0; row < inputs.size(); ++row)
        EXPECT_TRUE(CameraAgrees((*outputs)[row], sirius::test::retained_camera::kCases[row]));
    inputs.clear();
    for (const auto& fixture : sirius::test::continuous_retained_camera::cases) {
        RetainedCameraInput input;
        for (std::size_t i = 0; i < 32; ++i)
            input.values[i] = RetainedValue::FromDouble(fixture.input[i]);
        inputs.push_back(input);
    }
    ASSERT_EQ(inputs[0].values[25].high, inputs[1].values[25].high);
    ASSERT_NE(inputs[0].values[25].low, inputs[1].values[25].low);
    ASSERT_EQ(inputs[0].values[27].high, inputs[2].values[27].high);
    ASSERT_NE(inputs[0].values[27].low, inputs[2].values[27].low);
    outputs = camera_compute.Camera(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    for (std::size_t row = 0; row < inputs.size(); ++row)
        EXPECT_TRUE(
            CameraAgrees((*outputs)[row], sirius::test::continuous_retained_camera::cases[row]));
    for (auto& input : inputs)
        for (auto& value : input.values) {
            value.low = 0;
            value.tail = 0;
        }
    outputs = camera_compute.Camera(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    EXPECT_FALSE(CameraAgrees((*outputs)[1], sirius::test::continuous_retained_camera::cases[1]));
    EXPECT_FALSE(CameraAgrees((*outputs)[2], sirius::test::continuous_retained_camera::cases[2]));
    EXPECT_EQ(device->BufferAllocationBytes(), allocation);
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST_F(RetainedComputeTest, SmoothRayCameraPreservesPhysicalLensDerivatives) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    const auto allocation = device->BufferAllocationBytes();
    std::vector<RetainedRayCameraInput> inputs;
    for (const auto& fixture : sirius::test::retained_ray_camera::cases)
        inputs.push_back(std::bit_cast<RetainedRayCameraInput>(fixture.input));
    auto outputs = compute->RayCamera(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    ASSERT_EQ(outputs->size(), inputs.size());
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        SCOPED_TRACE(sirius::test::retained_ray_camera::cases[row].name);
        ASSERT_TRUE(CameraAgrees((*outputs)[row], sirius::test::retained_ray_camera::cases[row]));
        auto narrowed = (*outputs)[row];
        for (auto& value : narrowed.values) value.low = 0;
        EXPECT_FALSE(CameraAgrees(narrowed, sirius::test::retained_ray_camera::cases[row]));
    }
    inputs[2].values[22].high = std::numeric_limits<float>::quiet_NaN();
    outputs = compute->RayCamera(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        if (row == 2) {
            EXPECT_FALSE((*outputs)[row].valid);
            for (const auto& value : (*outputs)[row].values) EXPECT_EQ(value.valid, 0U);
        } else {
            EXPECT_TRUE(
                CameraAgrees((*outputs)[row], sirius::test::retained_ray_camera::cases[row]));
        }
    }
    inputs.resize(1);
    outputs = compute->RayCamera(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    EXPECT_TRUE(CameraAgrees(outputs->front(), sirius::test::retained_ray_camera::cases.front()));
    EXPECT_EQ(device->BufferAllocationBytes(), allocation);
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST_F(RetainedComputeTest, JointRkStagesRetainCriticalIncrementsAndEmbeddedError) {
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
}

TEST_F(RetainedComputeTest, ProjectedEndpointsKeepPhysicalColumnsAndRetainedContinuation) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    const auto allocation = device->BufferAllocationBytes();
    std::vector<RetainedEndpointInput> inputs;
    for (const auto& fixture : sirius::test::retained_endpoint::cases)
        inputs.push_back(std::bit_cast<RetainedEndpointInput>(fixture.input));
    const auto outputs = compute->Endpoint(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        const auto& fixture = sirius::test::retained_endpoint::cases[row];
        SCOPED_TRACE(fixture.name);
        const auto& output = (*outputs)[row];
        ASSERT_TRUE(output.valid);
        EXPECT_EQ(output.component, fixture.component);
        for (std::size_t i = 0; i < 80; ++i) {
            SCOPED_TRACE(i);
            const auto& value = i < 40 ? output.phase[i] : output.physical[i - 40];
            const auto oracle = fixture.reference[i];
            const sirius::core::Twofold exact(oracle.high, oracle.low);
            const auto center = sirius::core::Twofold(value.high) +
                                sirius::core::Twofold(value.low) +
                                sirius::core::Twofold(value.tail);
            const double difference = std::abs((center - exact).Rounded());
            EXPECT_LE(difference, value.radius + 1e-29 * (1 + std::abs(oracle.high)));
            EXPECT_LE(difference, 1e-11 * (1 + std::abs(oracle.high)));
        }
    }
    // Reuse the accepted phase without converting through the physical output
    // or a binary64 sum. Projection must not move an already projected state.
    for (std::size_t row = 0; row < inputs.size(); ++row)
        std::copy((*outputs)[row].phase.begin(), (*outputs)[row].phase.end(),
                  inputs[row].values.begin() + 4);
    const auto repeated = compute->Endpoint(inputs);
    ASSERT_TRUE(repeated) << repeated.error().Description();
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        ASSERT_TRUE((*repeated)[row].valid) << row;
        for (std::size_t i = 0; i < 40; ++i)
            EXPECT_NEAR((*repeated)[row].physical[i].Center(), (*outputs)[row].physical[i].Center(),
                        1e-11 * (1 + std::abs((*outputs)[row].physical[i].Center())));
    }
    inputs[0].values[4].valid = 0;
    inputs[1].values[8] = RetainedValue::FromDouble(100);
    const auto invalid = compute->Endpoint(inputs);
    ASSERT_TRUE(invalid) << invalid.error().Description();
    for (std::size_t row = 0; row < 2; ++row) {
        EXPECT_FALSE((*invalid)[row].valid);
        for (const auto& value : (*invalid)[row].physical) EXPECT_EQ(value.valid, 0U);
        for (const auto& value : (*invalid)[row].phase) EXPECT_EQ(value.valid, 0U);
    }
    EXPECT_EQ(device->BufferAllocationBytes(), allocation);
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST_F(RetainedComputeTest, DenseSegmentsPreserveSmallCovariantArrivalDerivatives) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    const auto allocation = device->BufferAllocationBytes();
    const auto& cases = sirius::test::retained_dense::cases;
    for (std::size_t begin = 0; begin < cases.size(); begin += compute->Capacity()) {
        const auto count = std::min(compute->Capacity(), cases.size() - begin);
        std::vector<RetainedDenseInput> inputs;
        for (std::size_t i = 0; i < count; ++i)
            inputs.push_back(std::bit_cast<RetainedDenseInput>(cases[begin + i].input));
        const auto outputs = compute->Dense(inputs);
        ASSERT_TRUE(outputs) << outputs.error().Description();
        for (std::size_t row = 0; row < count; ++row) {
            const auto& fixture = cases[begin + row];
            SCOPED_TRACE(fixture.name);
            ASSERT_TRUE((*outputs)[row].valid);
            for (std::size_t i = 0; i < 40; ++i) {
                SCOPED_TRACE(i);
                const auto& value = (*outputs)[row].physical[i];
                const auto oracle = fixture.reference[i];
                const auto center = sirius::core::Twofold(value.high) +
                                    sirius::core::Twofold(value.low) +
                                    sirius::core::Twofold(value.tail);
                const double difference =
                    std::abs((center - sirius::core::Twofold(oracle.high, oracle.low)).Rounded());
                EXPECT_LE(difference, value.radius + 1e-28 * (1 + std::abs(oracle.high)));
                EXPECT_LE(difference, 1e-10 * (1 + std::abs(oracle.high)));
            }
        }
    }
    // Endpoint ownership evaluates the retained sum directly. Cancellation
    // must promote the surviving sparse term into the leading output limb.
    const std::array<RetainedValue, 3> first{{
        {1, 0x1.000002p-35f, 0x1p-120f, 0, 1},
        {0x1p60f, 1, 0x1p-90f, 0, 1},
        {0x1p120f, -0x1p95f, -0x1p70f, 0, 1},
    }};
    const std::array<RetainedValue, 3> delta{{
        {-1, -0x1.000002p-35f, 0, 0, 1},
        {-0x1p60f, -1, 0x1p-100f, 0, 1},
        {-0x1p120f, 0x1p95f, 0x1p69f, 0, 1},
    }};
    const std::array<float, 3> sums{0x1p-120f, 0x1.004p-90f, -0x1p69f};
    std::vector<RetainedDenseInput> cancellation(4);
    for (auto& row : cancellation) {
        row.values.fill(RetainedValue::FromDouble(0));
        row.values[104] = row.values[105] = row.values[111] = RetainedValue::FromDouble(1);
    }
    for (std::size_t row = 0; row < first.size(); ++row) {
        cancellation[row].values[5] = first[row];
        cancellation[row].values[45] = RetainedValue::FromDouble(sums[row]);
        cancellation[row].values[85] = delta[row];
    }
    // With zero displacement and opposite endpoint slopes, the midpoint is
    // h*v/4. For e=2^-24, (1+e+e^2)*(1-e+e^2)/4 = (1+e^2+e^4)/4.
    auto& product = cancellation.back();
    product.values[9] = {1, 0x1p-24f, 0x1p-48f, 0, 1};
    product.values[49] = {-1, -0x1p-24f, -0x1p-48f, 0, 1};
    product.values[104] = {1, -0x1p-24f, 0x1p-48f, 0, 1};
    product.values[105] = RetainedValue::FromDouble(.5);
    const auto cancelled = compute->Dense(cancellation);
    ASSERT_TRUE(cancelled) << cancelled.error().Description();
    for (std::size_t row = 0; row < first.size(); ++row) {
        SCOPED_TRACE(row);
        ASSERT_TRUE((*cancelled)[row].valid);
        const auto& value = (*cancelled)[row].physical[1];
        EXPECT_EQ(value.high, sums[row]);
        EXPECT_EQ(value.low, 0);
        EXPECT_EQ(value.tail, 0);
        EXPECT_EQ(value.radius, 0);
    }
    ASSERT_TRUE(cancelled->back().valid);
    const auto& value = cancelled->back().physical[1];
    const auto center = sirius::core::Twofold(value.high) + sirius::core::Twofold(value.low) +
                        sirius::core::Twofold(value.tail);
    const auto exact = sirius::core::Twofold(.25) + sirius::core::Twofold(0x1p-50) +
                       sirius::core::Twofold(0x1p-98);
    EXPECT_LE(std::abs((center - exact).Rounded()), value.radius);
    EXPECT_LT(value.radius, 0x1p-60);
    EXPECT_NE(value.low, 0);

    auto input = std::bit_cast<RetainedDenseInput>(cases[1].input);
    for (std::size_t i = 106; i < 110; ++i) input.values[i] = RetainedValue::FromDouble(0);
    auto invalid = compute->Dense(std::span(&input, 1));
    ASSERT_TRUE(invalid) << invalid.error().Description();
    EXPECT_FALSE(invalid->front().valid);
    for (const auto& rejected_value : invalid->front().physical)
        EXPECT_EQ(rejected_value.valid, 0U);
    input = std::bit_cast<RetainedDenseInput>(cases[0].input);
    input.values[105] = RetainedValue::FromDouble(1.01);
    invalid = compute->Dense(std::span(&input, 1));
    ASSERT_TRUE(invalid) << invalid.error().Description();
    EXPECT_FALSE(invalid->front().valid);
    EXPECT_EQ(device->BufferAllocationBytes(), allocation);
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST_F(RetainedComputeTest, PhysicalInitializationRetainsTheHamiltonianResidual) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    const auto& initial = sirius::test::critical_fixture::transport_cases;
    std::vector<RetainedInitializeInput> inputs(initial.size());
    for (std::size_t row = 0; row < initial.size(); ++row) {
        auto& input = inputs[row];
        const auto original =
            std::bit_cast<RetainedStepInput>(sirius::test::retained_transport::cases[row].input);
        std::copy_n(original.values.begin(), 4, input.values.begin());
        for (std::size_t i = 0; i < 40; ++i)
            input.values[4 + i] = RetainedValue::FromDouble(initial[row].initial[i]);
        input.values[44] = RetainedValue::FromDouble(-1);
    }
    const auto outputs = compute->Initialize(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        SCOPED_TRACE(row);
        ASSERT_TRUE((*outputs)[row].valid);
        const auto expected =
            std::bit_cast<RetainedStepInput>(sirius::test::retained_transport::cases[row].input);
        for (std::size_t i = 0; i < 40; ++i) {
            SCOPED_TRACE(i);
            const auto& actual = (*outputs)[row].phase[i];
            const auto& reference = expected.values[i + 4];
            const auto center = [](const RetainedValue& value) {
                return sirius::core::Twofold(value.high) + sirius::core::Twofold(value.low) +
                       sirius::core::Twofold(value.tail);
            };
            const double difference = std::abs((center(actual) - center(reference)).Rounded());
            EXPECT_LE(difference, double(actual.radius) + reference.radius +
                                      1e-29 * (1 + std::abs(reference.high)));
            EXPECT_LE(difference, 1e-11 * (1 + std::abs(reference.high)));
        }
    }
    inputs[0].values[16].valid = 0;
    const auto invalid = compute->Initialize(std::span(inputs.data(), 1));
    ASSERT_TRUE(invalid) << invalid.error().Description();
    EXPECT_FALSE(invalid->front().valid);
    for (const auto& value : invalid->front().phase) EXPECT_EQ(value.valid, 0U);
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST_F(RetainedComputeTest, CoupledIntervalsRequireEmbeddedAndIndependentDenseAgreement) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    using Clock = std::chrono::steady_clock;
    using Stats = std::array<RetainedCompute::StageStats, 6>;
    const auto fixture_started = Clock::now();
    const auto fixture_before = compute->Statistics();
    const auto milliseconds = [](Clock::time_point started, Clock::time_point finished) {
        return std::chrono::duration<double, std::milli>(finished - started).count();
    };
    const auto record_timing = [](const std::string& prefix, const Stats& before,
                                  const Stats& after, double wall_ms,
                                  bool lifetime_maxima = false) {
        constexpr std::array names{"camera", "transport",  "endpoint",
                                   "dense",  "initialize", "ray_camera"};
        std::uint64_t submissions = 0, pipeline_creations = 0, target_overshoots = 0;
        std::uint64_t write_buffer_bytes = 0, read_buffer_bytes = 0;
        double submit_wait_ms = 0, pipeline_setup_ms = 0, maximum_submit_wait_ms = 0;
        double command_setup_ms = 0, cleanup_ms = 0, dispatch_total_ms = 0;
        double write_buffer_ms = 0, read_buffer_ms = 0;
        for (std::size_t stage = 0; stage < after.size(); ++stage) {
            const auto count = after[stage].submissions - before[stage].submissions;
            const auto wait = after[stage].submit_wait_ms - before[stage].submit_wait_ms;
            const auto pipeline = after[stage].pipeline_setup_ms - before[stage].pipeline_setup_ms;
            const auto command = after[stage].command_setup_ms - before[stage].command_setup_ms;
            const auto cleanup = after[stage].cleanup_ms - before[stage].cleanup_ms;
            const auto dispatch = after[stage].dispatch_total_ms - before[stage].dispatch_total_ms;
            const auto write = after[stage].write_buffer_ms - before[stage].write_buffer_ms;
            const auto read = after[stage].read_buffer_ms - before[stage].read_buffer_ms;
            const auto written_bytes =
                after[stage].write_buffer_bytes - before[stage].write_buffer_bytes;
            const auto read_bytes =
                after[stage].read_buffer_bytes - before[stage].read_buffer_bytes;
            const auto creations =
                after[stage].pipeline_creations - before[stage].pipeline_creations;
            const auto overshoots =
                after[stage].target_overshoots - before[stage].target_overshoots;
            submissions += count;
            submit_wait_ms += wait;
            pipeline_setup_ms += pipeline;
            command_setup_ms += command;
            cleanup_ms += cleanup;
            dispatch_total_ms += dispatch;
            write_buffer_ms += write;
            read_buffer_ms += read;
            write_buffer_bytes += written_bytes;
            read_buffer_bytes += read_bytes;
            pipeline_creations += creations;
            target_overshoots += overshoots;
            maximum_submit_wait_ms =
                std::max(maximum_submit_wait_ms, after[stage].maximum_submit_wait_ms);
            const auto key = prefix + "_" + names[stage];
            RecordProperty(key + "_submissions", std::to_string(count));
            RecordProperty(key + "_submit_wait_ms", std::to_string(wait));
            RecordProperty(key + "_pipeline_setup_ms", std::to_string(pipeline));
            RecordProperty(key + "_command_setup_ms", std::to_string(command));
            RecordProperty(key + "_cleanup_ms", std::to_string(cleanup));
            RecordProperty(key + "_dispatch_total_ms", std::to_string(dispatch));
            RecordProperty(key + "_write_buffer_ms", std::to_string(write));
            RecordProperty(key + "_read_buffer_ms", std::to_string(read));
            RecordProperty(key + "_write_buffer_bytes", std::to_string(written_bytes));
            RecordProperty(key + "_read_buffer_bytes", std::to_string(read_bytes));
            RecordProperty(key + "_pipeline_creations", std::to_string(creations));
            RecordProperty(key + "_target_overshoots", std::to_string(overshoots));
            if (lifetime_maxima) {
                RecordProperty(key + "_maximum_submit_wait_ms",
                               std::to_string(after[stage].maximum_submit_wait_ms));
            }
        }
        RecordProperty(prefix + "_submissions", std::to_string(submissions));
        RecordProperty(prefix + "_submit_wait_ms", std::to_string(submit_wait_ms));
        RecordProperty(prefix + "_pipeline_setup_ms", std::to_string(pipeline_setup_ms));
        RecordProperty(prefix + "_command_setup_ms", std::to_string(command_setup_ms));
        RecordProperty(prefix + "_cleanup_ms", std::to_string(cleanup_ms));
        RecordProperty(prefix + "_dispatch_total_ms", std::to_string(dispatch_total_ms));
        RecordProperty(prefix + "_write_buffer_ms", std::to_string(write_buffer_ms));
        RecordProperty(prefix + "_read_buffer_ms", std::to_string(read_buffer_ms));
        RecordProperty(prefix + "_write_buffer_bytes", std::to_string(write_buffer_bytes));
        RecordProperty(prefix + "_read_buffer_bytes", std::to_string(read_buffer_bytes));
        RecordProperty(prefix + "_pipeline_creations", std::to_string(pipeline_creations));
        RecordProperty(prefix + "_target_overshoots", std::to_string(target_overshoots));
        const double measured_interface_ms = write_buffer_ms + dispatch_total_ms + read_buffer_ms;
        RecordProperty(prefix + "_measured_interface_ms", std::to_string(measured_interface_ms));
        RecordProperty(prefix + "_dispatch_total_minus_phases_ms",
                       std::to_string(dispatch_total_ms - pipeline_setup_ms - command_setup_ms -
                                      submit_wait_ms - cleanup_ms));
        RecordProperty(prefix + "_wall_ms", std::to_string(wall_ms));
        RecordProperty(prefix + "_wall_minus_measured_ms",
                       std::to_string(wall_ms - measured_interface_ms));
        // Peaks are lifetime maxima, not additive counters. Only the complete
        // fresh-compute fixture can report them as maxima for its own window.
        if (lifetime_maxima)
            RecordProperty(prefix + "_maximum_submit_wait_ms",
                           std::to_string(maximum_submit_wait_ms));
    };
    RecordProperty("device", device->Info().name);
    RecordProperty("retained_capacity", std::to_string(compute->Capacity()));
    RecordProperty("explicit_buffer_bytes", std::to_string(device->BufferAllocationBytes()));
    auto* vulkan = dynamic_cast<VulkanDevice*>(device.get());
    ASSERT_NE(vulkan, nullptr);
    const auto cache_initial = vulkan->PipelineCacheStatistics();
    RecordProperty("pipeline_cache_enabled", cache_initial.enabled ? "true" : "false");
    RecordProperty("pipeline_cache_creation_result",
                   std::to_string(static_cast<int>(cache_initial.creation_result)));
    RecordProperty("pipeline_cache_imported_bytes", std::to_string(cache_initial.imported_bytes));
    RecordProperty("pipeline_cache_import_discarded",
                   cache_initial.import_discarded ? "true" : "false");
    RecordProperty("pipeline_cache_blob_bound_bytes",
                   std::to_string(detail::kVulkanPipelineCacheBlobLimit));
    RecordProperty("timing_scope",
                   "host wall-clock observations, not GPU execution timestamps; interval fixture "
                   "excludes device/compute SetUp; wall-minus-measured subtracts WriteBuffer, "
                   "ReadBuffer and inclusive Dispatch total once, leaving host packing, decoding, "
                   "coordination, assertions and diagnostic bookkeeping; no frame or interactive "
                   "qualification");
    RecordProperty("timing_failure_scope",
                   "buffer timers include returned errors; buffer byte counters include only "
                   "successful spans; dispatch phases and pipeline creations "
                   "require device Dispatch success with valid submission timing; read failure "
                   "retains successful dispatch observations; failed dispatch partial phases are "
                   "not accumulated");
    RecordProperty(
        "pipeline_cache_scope",
        "initial Endpoint warms projection; first interval first uses transport/dense; "
        "repeated intervals reuse stage kernels; pipeline_creations counts successful dispatches "
        "that created a pipeline object, not driver pipeline-cache hits");
    std::vector<RetainedEndpointInput> initial;
    for (const auto& fixture : sirius::test::retained_transport::cases) {
        const auto step = std::bit_cast<RetainedStepInput>(fixture.input);
        RetainedEndpointInput input;
        std::copy_n(step.values.begin(), 45, input.values.begin());
        initial.push_back(input);
    }
    const auto launches = compute->Endpoint(initial);
    ASSERT_TRUE(launches) << launches.error().Description();
    std::vector<RetainedIntervalInput> inputs(initial.size());
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        ASSERT_TRUE((*launches)[row].valid) << row;
        const auto step =
            std::bit_cast<RetainedStepInput>(sirius::test::retained_transport::cases[row].input);
        auto& input = inputs[row];
        std::copy_n(step.values.begin(), 4, input.metric.begin());
        input.start = (*launches)[row];
        input.chart = step.values[44].Center();
        input.interval = step.values[45].Center();
        input.control.length_scale = 1;
        input.control.frequency_scale = 1;
        input.control.tolerance = 1e-4 / (4 * 30000);
        input.control.integrator.min_step = static_cast<float>(input.interval);
        input.control.integrator.max_step = 1;
        input.control.integrator.abs_tolerance = 1e-9f;
        input.control.integrator.rel_tolerance = 1e-9f;
    }
    const auto original_inputs = inputs;
    const auto first_before = compute->Statistics();
    const auto first_started = Clock::now();
    const auto outputs = AttemptRetainedIntervals(*compute, inputs);
    const auto first_finished = Clock::now();
    record_timing("first_interval", first_before, compute->Statistics(),
                  milliseconds(first_started, first_finished));
    ASSERT_TRUE(outputs) << outputs.error().Description();
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        SCOPED_TRACE(sirius::test::retained_transport::cases[row].name);
        const auto& output = (*outputs)[row];
        EXPECT_TRUE(output.admissible)
            << "error=" << output.error_ratio
            << " failure=" << sirius::core::CoupledStepFailureName(output.failure);
        EXPECT_EQ(output.attempted_stages, 21U);
        EXPECT_LE(output.error_ratio, 1);
    }
    // Exact flat flow remains exact through all three independent candidates.
    const auto& flat = outputs->back();
    ASSERT_TRUE(flat.admissible);
    EXPECT_EQ(flat.error_ratio, 0);
    EXPECT_EQ(flat.full.physical[0].Center(), .75);
    EXPECT_EQ(flat.full.physical[3].Center(), 4.25);
    EXPECT_EQ(flat.midpoint.physical[0].Center(), .875);
    EXPECT_EQ(flat.refined.physical[3].Center(), 4.25);

    // Sparse retained tails must survive the lower-order increment split even
    // when the embedded error is exactly zero in flat space.
    auto sparse = inputs.back();
    sparse.start.phase[13] = sparse.start.physical[13] = {1, 0x1.000002p-35f, 0x1p-120f, 0, 1};
    const auto sparse_before = compute->Statistics();
    const auto sparse_started = Clock::now();
    const auto sparse_output = AttemptRetainedIntervals(*compute, {&sparse, 1});
    const auto sparse_finished = Clock::now();
    record_timing("repeated_sparse_interval", sparse_before, compute->Statistics(),
                  milliseconds(sparse_started, sparse_finished));
    ASSERT_TRUE(sparse_output) << sparse_output.error().Description();
    ASSERT_TRUE(sparse_output->front().admissible);
    const auto& full_increment = sparse_output->front().full_increment[5];
    const auto& lower_increment = sparse_output->front().lower_increment[5];
    EXPECT_EQ(full_increment.tail, 0x1p-122f);
    EXPECT_EQ(lower_increment.high, full_increment.high);
    EXPECT_EQ(lower_increment.low, full_increment.low);
    EXPECT_EQ(lower_increment.tail, full_increment.tail);

    inputs.resize(2);
    inputs[0].control.tolerance = 1e-30;
    inputs[1].start.phase[5].valid = 0;
    const auto rejected_before = compute->Statistics();
    const auto rejected_started = Clock::now();
    const auto rejected = AttemptRetainedIntervals(*compute, inputs);
    const auto rejected_finished = Clock::now();
    record_timing("repeated_rejection_interval", rejected_before, compute->Statistics(),
                  milliseconds(rejected_started, rejected_finished));
    ASSERT_TRUE(rejected) << rejected.error().Description();
    for (const auto& output : *rejected) {
        EXPECT_FALSE(output.admissible);
        EXPECT_NE(output.failure, sirius::core::CoupledStepFailure::None);
        EXPECT_FALSE(output.full.valid);
        EXPECT_FALSE(output.lower.valid);
        EXPECT_FALSE(output.midpoint.valid);
        EXPECT_FALSE(output.refined.valid);
        for (const auto& value : output.full_increment) EXPECT_EQ(value.valid, 0U);
    }
    record_timing("interval_fixture", fixture_before, compute->Statistics(),
                  milliseconds(fixture_started, Clock::now()), true);
    // Keep serialization outside the unchanged interval timing window. A
    // repeated gtest iteration creates a fresh device but can import this blob.
    if (cache_initial.enabled) {
        const auto cache_started = Clock::now();
        const auto cache = vulkan->SnapshotPipelineCache();
        RecordProperty("pipeline_cache_snapshot_wall_ms",
                       std::to_string(milliseconds(cache_started, Clock::now())));
        if (cache) {
            RecordProperty("pipeline_cache_available_bytes",
                           std::to_string(cache->available_bytes));
            RecordProperty("pipeline_cache_exported_bytes", std::to_string(cache->exported_bytes));
            RecordProperty("pipeline_cache_export_discarded",
                           cache->export_discarded ? "true" : "false");
            EXPECT_LE(cache->exported_bytes, detail::kVulkanPipelineCacheBlobLimit);
        } else {
            // Cache export is optional in the product. An export failure leaves
            // its performance benefit unmeasured, not its physics disproved.
            RecordProperty("pipeline_cache_export_error", cache.error().Description());
        }
    }
    RecordProperty("pipeline_cache_memory_scope",
                   "one process-lived serialized blob is bounded; driver-internal pipeline/cache "
                   "allocations are separate from explicit_buffer_bytes; no frame qualification");
    // These bounded equivalence and governor checks are outside the original
    // interval timing window and keep its independent science cases unchanged.
    TransferProbeDevice probe(*device);
    const auto initial_allocation = device->BufferAllocationBytes();
    // An intentionally tiny soft target exercises actual feedback/subdivision;
    // its elapsed timings do not claim a production performance result.
    auto observed = RetainedCompute::Create(probe, 4, false, 1e-12);
    ASSERT_TRUE(observed) << observed.error().Description();
    auto& paired_compute = **observed;
    const auto allocated = device->BufferAllocationBytes();
    EXPECT_GT(allocated, initial_allocation);
    EXPECT_LE(allocated, 8ull * 1024 * 1024);
    ASSERT_EQ(probe.allocations.size(), 12U);
    const auto endpoint_input = probe.allocations[4].handle;
    const auto endpoint_output = probe.allocations[5].handle;
    const auto endpoint_row_bytes = probe.allocations[5].bytes / paired_compute.Capacity();
    const auto check_endpoint_transfers = [&](std::size_t rows, std::size_t calls,
                                              std::size_t writes_begin, std::size_t reads_begin,
                                              std::size_t submissions_begin) {
        std::size_t writes = 0, reads = 0, submissions = 0;
        for (std::size_t i = writes_begin; i < probe.writes.size(); ++i)
            if (probe.writes[i].handle.value == endpoint_input.value) {
                ++writes;
                EXPECT_EQ(probe.writes[i].capacity, paired_compute.Capacity());
                EXPECT_EQ(probe.writes[i].bytes,
                          sizeof(std::uint32_t) + rows * sizeof(RetainedEndpointInput));
            }
        for (std::size_t i = reads_begin; i < probe.reads.size(); ++i)
            if (probe.reads[i].handle.value == endpoint_output.value) {
                ++reads;
                EXPECT_EQ(probe.reads[i].bytes, rows * endpoint_row_bytes);
            }
        for (std::size_t i = submissions_begin; i < probe.submissions.size(); ++i)
            if (probe.submissions[i].buffers[0].value == endpoint_input.value) {
                ++submissions;
                EXPECT_EQ(probe.submissions[i].binding_count, 2U);
                EXPECT_EQ(probe.submissions[i].buffers[1].value, endpoint_output.value);
                // The existing retained shader gives each row one workgroup.
                EXPECT_EQ(probe.submissions[i].x, rows);
                EXPECT_EQ(probe.submissions[i].y, 1U);
                EXPECT_EQ(probe.submissions[i].z, 1U);
            }
        EXPECT_EQ(writes, calls);
        EXPECT_EQ(reads, calls);
        EXPECT_EQ(submissions, calls);
    };
    const auto run = [&](std::span<const RetainedIntervalInput> batch, std::size_t budget,
                         std::uint64_t endpoint_calls, std::size_t endpoint_rows,
                         bool first_upload = false) {
        const auto before = paired_compute.Statistics();
        const auto writes_begin = probe.writes.size(), reads_begin = probe.reads.size(),
                   submissions_begin = probe.submissions.size();
        const auto dispatches = probe.forwarded_dispatches;
        auto result = AttemptRetainedIntervals(paired_compute, batch, budget);
        const auto after = paired_compute.Statistics();
        EXPECT_EQ(after[1].submissions - before[1].submissions, 3U);
        EXPECT_EQ(after[2].submissions - before[2].submissions, endpoint_calls);
        EXPECT_EQ(after[3].submissions - before[3].submissions, 1U);
        EXPECT_EQ(probe.forwarded_dispatches - dispatches, 4 + endpoint_calls);
        EXPECT_EQ(device->BufferAllocationBytes(), allocated);
        EXPECT_EQ(probe.allocations.size(), 12U);
        if (!first_upload)
            check_endpoint_transfers(endpoint_rows, endpoint_calls, writes_begin, reads_begin,
                                     submissions_begin);
        return result;
    };
    // Include a critical curved row and the independently pinned sparse flat
    // row. Four projection rows exactly fit capacity; budget three must fall back.
    const std::array pair_inputs{original_inputs.front(), sparse};
    const auto serialized_before = paired_compute.Statistics();
    const auto serialized = run(pair_inputs, 2, 6, 2, true);
    const auto serialized_endpoint_calls =
        paired_compute.Statistics()[2].submissions - serialized_before[2].submissions;
    ASSERT_TRUE(serialized) << serialized.error().Description();
    for (const auto& output : *serialized) ASSERT_TRUE(output.admissible);
    const auto paired_before = paired_compute.Statistics();
    const auto paired = run(pair_inputs, 4, 3, 4);
    const auto paired_endpoint_calls =
        paired_compute.Statistics()[2].submissions - paired_before[2].submissions;
    ASSERT_TRUE(paired) << paired.error().Description();
    const auto fallback = run(pair_inputs, 3, 6, 2);
    ASSERT_TRUE(fallback) << fallback.error().Description();
    for (std::size_t row = 0; row < pair_inputs.size(); ++row) {
        SCOPED_TRACE(row);
        EXPECT_TRUE(IntervalBitsAgree((*paired)[row], (*serialized)[row]));
        EXPECT_TRUE(IntervalBitsAgree((*fallback)[row], (*serialized)[row]));
    }
    auto mixed = pair_inputs;
    mixed[1].start.phase[5].valid = 0;
    const auto mixed_paired = run(mixed, 4, 3, 4);
    const auto mixed_serialized = run(mixed, 2, 6, 2);
    ASSERT_TRUE(mixed_paired) << mixed_paired.error().Description();
    ASSERT_TRUE(mixed_serialized) << mixed_serialized.error().Description();
    ASSERT_TRUE((*mixed_paired)[0].admissible);
    EXPECT_TRUE(IntervalBitsAgree((*mixed_paired)[0], (*paired)[0]));
    EXPECT_FALSE((*mixed_paired)[1].admissible);
    EXPECT_EQ((*mixed_paired)[1].failure, sirius::core::CoupledStepFailure::InvalidState);
    for (std::size_t row = 0; row < mixed.size(); ++row)
        EXPECT_TRUE(IntervalBitsAgree((*mixed_paired)[row], (*mixed_serialized)[row]));
    auto rejected_pair = pair_inputs;
    rejected_pair[0].control.tolerance = 1e-30;
    const auto rejected_paired = run(rejected_pair, 4, 3, 4);
    const auto rejected_serialized = run(rejected_pair, 2, 6, 2);
    ASSERT_TRUE(rejected_paired) << rejected_paired.error().Description();
    ASSERT_TRUE(rejected_serialized) << rejected_serialized.error().Description();
    ASSERT_FALSE((*rejected_paired)[0].admissible);
    ASSERT_TRUE((*rejected_paired)[1].admissible);
    EXPECT_TRUE(IntervalBitsAgree((*rejected_paired)[1], (*paired)[1]));
    for (std::size_t row = 0; row < rejected_pair.size(); ++row)
        EXPECT_TRUE(IntervalBitsAgree((*rejected_paired)[row], (*rejected_serialized)[row]));
    for (const auto* rejected_row : {&(*mixed_paired)[1], &(*rejected_paired)[0]}) {
        for (const auto* endpoint : {&rejected_row->full, &rejected_row->lower,
                                     &rejected_row->midpoint, &rejected_row->refined}) {
            EXPECT_FALSE(endpoint->valid);
            for (const auto& value : endpoint->phase) EXPECT_EQ(value.valid, 0U);
            for (const auto& value : endpoint->physical) EXPECT_EQ(value.valid, 0U);
        }
        for (const auto* increment :
             {&rejected_row->full_increment, &rejected_row->lower_increment,
              &rejected_row->midpoint_increment, &rejected_row->refined_increment})
            for (const auto& value : *increment) EXPECT_EQ(value.valid, 0U);
    }
    const auto writes_before_invalid = probe.writes.size();
    const auto dispatches_before_invalid = probe.dispatch_calls;
    for (const std::size_t budget : {0U, 1U, 5U}) {
        const auto invalid = AttemptRetainedIntervals(paired_compute, pair_inputs, budget);
        EXPECT_FALSE(invalid);
    }
    EXPECT_EQ(probe.writes.size(), writes_before_invalid);
    EXPECT_EQ(probe.dispatch_calls, dispatches_before_invalid);
    const auto default_before = paired_compute.Statistics();
    const auto capacity_default = AttemptRetainedIntervals(paired_compute, pair_inputs);
    ASSERT_TRUE(capacity_default) << capacity_default.error().Description();
    EXPECT_EQ(paired_compute.Statistics()[2].submissions - default_before[2].submissions, 3U);
    for (std::size_t row = 0; row < pair_inputs.size(); ++row)
        EXPECT_TRUE(IntervalBitsAgree((*capacity_default)[row], (*serialized)[row]));

    // Prove the live coordinator forwards its reduced limit, rather than the
    // static capacity. The first one-row attempt can pair; measured positive
    // submission time exceeds the tiny target and the next attempt's limit is one.
    EXPECT_GT(paired_compute.TakeSubmissionPeakMs(), 0);
    RetainedTraceExecutor bounded_executor(paired_compute);
    sirius::core::KerrSchildFamily flat_metric(sirius::core::KerrSchildParams::Minkowski());
    sirius::core::Lightray ray{};
    ray.position(1) = 5;
    ray.velocity(0) = -1;
    ray.velocity(1) = 1;
    ray.step_size = 1;
    sirius::core::Rk45CoupledState coupled;
    coupled.length_scale = coupled.frequency_scale = 1;
    coupled.tolerance = 1e-9;
    coupled.variations[0].derivative(2) = .001;
    coupled.variations[1].derivative(3) = .001;
    coupled.variations[2].displacement(2) = 1;
    coupled.variations[3].displacement(3) = 1;
    sirius::core::IntegratorConfig control;
    control.min_step = .01f;
    control.max_step = 2;
    sirius::core::Rk45CoupledComparison comparison;
    const auto first_governed = paired_compute.Statistics();
    ASSERT_TRUE(bounded_executor.Step(ray, flat_metric, control, coupled, comparison));
    EXPECT_EQ(paired_compute.Statistics()[2].submissions - first_governed[2].submissions, 4U);
    ASSERT_GT(bounded_executor.Statistics().batch_subdivisions, 0U);
    const auto second_governed = paired_compute.Statistics();
    const auto governed_writes = probe.writes.size(), governed_reads = probe.reads.size(),
               governed_submissions = probe.submissions.size();
    ASSERT_TRUE(bounded_executor.Step(ray, flat_metric, control, coupled, comparison));
    EXPECT_EQ(paired_compute.Statistics()[1].submissions - second_governed[1].submissions, 3U);
    EXPECT_EQ(paired_compute.Statistics()[2].submissions - second_governed[2].submissions, 6U);
    EXPECT_EQ(paired_compute.Statistics()[3].submissions - second_governed[3].submissions, 1U);
    ASSERT_NO_FATAL_FAILURE(
        check_endpoint_transfers(1, 6, governed_writes, governed_reads, governed_submissions));
    EXPECT_EQ(device->BufferAllocationBytes(), allocated);
    std::uint32_t maximum_governed_endpoint_rows = 0;
    for (std::size_t i = governed_submissions; i < probe.submissions.size(); ++i)
        if (probe.submissions[i].buffers[0].value == endpoint_input.value)
            maximum_governed_endpoint_rows =
                std::max(maximum_governed_endpoint_rows, probe.submissions[i].x);
    RecordProperty("endpoint_pairing_paired_endpoint_submissions",
                   std::to_string(paired_endpoint_calls));
    RecordProperty("endpoint_pairing_serialized_endpoint_submissions",
                   std::to_string(serialized_endpoint_calls));
    RecordProperty("endpoint_pairing_governed_endpoint_rows",
                   std::to_string(maximum_governed_endpoint_rows));
    RecordProperty("endpoint_pairing_scope",
                   "stored-word equivalence within one product mode; capacity-four fit/fallback, "
                   "inactive and rejected row positions, and actual tiny-target subdivision; "
                   "outside interval_fixture timing, no frame-speed claim");

#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST_F(RetainedComputeTest, SharedTracerCompletesDeviceIntervalsAndRetainsRollbackState) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    std::atomic<bool> cancelled{false};
    RetainedTraceExecutor executor(*compute, [&] { return cancelled.load(); });
    for (const double mass : {0.0, 1.0}) {
        sirius::core::KerrSchildFamily metric(mass == 0
                                                  ? sirius::core::KerrSchildParams::Minkowski()
                                                  : sirius::core::KerrSchildParams::Kerr(1, .7));
        sirius::core::CameraConfig config;
        config.width = 32;
        config.height = 20;
        config.theta = 1.1;
        config.phi = .2;
        config.beta_x = .1;
        config.beta_y = .8;
        config.beta_z = .01;
        config.focus_distance = 50;
        sirius::core::ThinLensCamera camera(config);
        const auto film = camera.ProjectFilmForObserver(4.125, 7.5, .25f, .6f);
        ASSERT_TRUE(film);
        const auto expected = sirius::core::LaunchCameraRay(metric, mass * .7, film->ray);
        const auto actual = executor.Launch(metric, mass * .7, film->ray);
        ASSERT_TRUE(expected);
        ASSERT_TRUE(actual);
        const auto compare = [](double observed, double reference) {
            EXPECT_NEAR(observed, reference, 1e-10 * (1 + std::abs(reference)));
        };
        for (int axis = 0; axis < 4; ++axis) {
            compare(actual->position(axis), expected->position(axis));
            compare(actual->tangent(axis), expected->tangent(axis));
            compare(actual->observer.time(axis), expected->observer.time(axis));
            for (int basis = 0; basis < 3; ++basis)
                compare(actual->observer.spatial[basis](axis),
                        expected->observer.spatial[basis](axis));
            for (int column = 0; column < 4; ++column) {
                compare(actual->variations[column].displacement(axis),
                        expected->variations[column].displacement(axis));
                compare(actual->variations[column].derivative(axis),
                        expected->variations[column].derivative(axis));
            }
        }
    }
    std::vector<std::future<std::array<TraceResult, 2>>> workers;
    for (int worker = 0; worker < 4; ++worker) {
        workers.push_back(std::async(std::launch::async, [&, worker] {
            sirius::core::KerrSchildFamily metric(
                worker % 2 == 0 ? sirius::core::KerrSchildParams::Minkowski()
                                : sirius::core::KerrSchildParams::Kerr(1, .5));
            TracerConfig config;
            config.enable_disk = false;
            config.enable_ray_bundles = true;
            config.bundle_point_source = true;
            config.escape_radius = 12;
            config.max_steps = 1000;
            config.integrator.min_step = 1e-5f;
            config.integrator.max_step = 2;
            config.integrator.initial_step = 1;
            sirius::core::CameraRay ray;
            ray.origin(1) = 5;
            ray.origin(2) = std::numbers::pi / 2;
            ray.origin(3) = worker * .2;
            ray.direction(1) = 1;
            GeodesicTracer device_trace(&metric, config), reference(&metric, config);
            device_trace.SetStepExecutor(&executor);
            return std::array{device_trace.Trace(ray), reference.Trace(ray)};
        }));
    }
    for (std::size_t worker = 0; worker < workers.size(); ++worker) {
        SCOPED_TRACE(worker);
        const auto results = workers[worker].get();
        const auto& actual = results[0];
        const auto& expected = results[1];
        ASSERT_FALSE(actual.numerical_failure)
            << sirius::core::CoupledStepFailureName(actual.coupled_failure)
            << " termination=" << actual.integrator_termination
            << " attempts=" << actual.steps_taken;
        EXPECT_EQ(actual.outcome, TraceResult::Outcome::Escaped);
        EXPECT_EQ(actual.outcome, expected.outcome);
        EXPECT_GT(actual.central_stages, 0U);
        EXPECT_TRUE(actual.beam.valid);
        EXPECT_NEAR(actual.affine_length, expected.affine_length, 1e-5);
        EXPECT_NEAR(actual.redshift, expected.redshift, 1e-5);
        for (int axis = 0; axis < 4; ++axis) {
            EXPECT_NEAR(actual.final_position(axis), expected.final_position(axis), 1e-5);
            EXPECT_NEAR(actual.final_direction(axis), expected.final_direction(axis), 1e-5);
        }
    }
    ASSERT_FALSE(executor.Error());
    EXPECT_GT(executor.Statistics().interval_batches, 0U);
    EXPECT_GT(executor.Statistics().camera_batches, 0U);
    EXPECT_GT(executor.Statistics().reused_phases, 0U);
    EXPECT_EQ(executor.Statistics().initialized_phases, 4U);

    sirius::core::KerrSchildFamily flat(sirius::core::KerrSchildParams::Minkowski());
    sirius::core::Lightray ray{};
    ray.position(1) = 5;
    ray.velocity(0) = -1;
    ray.velocity(1) = 1;
    ray.step_size = 1;
    sirius::core::Rk45CoupledState coupled;
    coupled.length_scale = 1;
    coupled.frequency_scale = 1;
    coupled.tolerance = 1e-9;
    coupled.variations[0].derivative(2) = .001;
    coupled.variations[1].derivative(3) = .001;
    coupled.variations[2].displacement(2) = 1;
    coupled.variations[3].displacement(3) = 1;
    sirius::core::IntegratorConfig config;
    config.min_step = .01f;
    config.max_step = 2;
    const auto before = ray;
    const auto columns = coupled;
    sirius::core::Rk45CoupledComparison comparison;
    ASSERT_TRUE(executor.Step(ray, flat, config, coupled, comparison));
    const auto midpoint = comparison.midpoint;
    const auto reused = executor.Statistics().reused_phases;
    // Mimic a rejected localized event: restore the original public state and
    // retry half the interval. The future phase must not contaminate this retry.
    executor.RejectLastInterval();
    ray = before;
    coupled = columns;
    ray.step_size = .5f;
    ASSERT_TRUE(executor.Step(ray, flat, config, coupled, comparison));
    EXPECT_GT(executor.Statistics().reused_phases, reused);
    for (int axis = 0; axis < 4; ++axis) {
        EXPECT_EQ(ray.position(axis), midpoint.position(axis));
        EXPECT_EQ(ray.velocity(axis), midpoint.velocity(axis));
    }
    const auto phase_reuses = executor.Statistics().reused_phases;
    ray.proper_time = std::nextafter(ray.proper_time, std::numeric_limits<float>::infinity());
    ASSERT_TRUE(executor.Step(ray, flat, config, coupled, comparison));
    EXPECT_GT(executor.Statistics().reused_phases, phase_reuses);

    // Nested scopes count one live thread. Ending the inner scope must leave
    // that thread registered; standalone calls after the outer end use the
    // bounded fallback rather than falsely claiming all live traces are queued.
    executor.BeginTrace();
    executor.BeginTrace();
    auto trace_ready = executor.Statistics().coalescing_traces_ready;
    EXPECT_TRUE(executor.Step(ray, flat, config, coupled, comparison));
    EXPECT_GT(executor.Statistics().coalescing_traces_ready, trace_ready);
    executor.EndTrace();
    trace_ready = executor.Statistics().coalescing_traces_ready;
    EXPECT_TRUE(executor.Step(ray, flat, config, coupled, comparison));
    EXPECT_GT(executor.Statistics().coalescing_traces_ready, trace_ready);
    executor.EndTrace();
    trace_ready = executor.Statistics().coalescing_traces_ready;
    EXPECT_TRUE(executor.Step(ray, flat, config, coupled, comparison));
    EXPECT_EQ(executor.Statistics().coalescing_traces_ready, trace_ready);

    {
        // One registered worker is deliberately doing host work. A concurrent
        // standalone camera request cannot stand in for that missing worker.
        executor.BeginTrace();
        std::promise<void> registered, release;
        const auto released = release.get_future().share();
        auto held = std::async(std::launch::async, [&] {
            executor.BeginTrace();
            registered.set_value();
            released.wait();
            executor.EndTrace();
        });
        struct ReleaseHeldTrace {
            std::promise<void>& release;
            std::future<void>& held;
            RetainedTraceExecutor& executor;
            ~ReleaseHeldTrace() {
                release.set_value();
                held.wait();
                executor.EndTrace();
            }
        } release_scope{release, held, executor};
        registered.get_future().wait();
        const auto mixed_before = executor.Statistics();
        trace_ready = executor.Statistics().coalescing_traces_ready;
        auto standalone = std::async(std::launch::async, [&] {
            sirius::core::KerrSchildFamily camera_metric(
                sirius::core::KerrSchildParams::Minkowski());
            sirius::core::CameraRay camera;
            camera.origin(1) = 5;
            camera.origin(2) = std::numbers::pi / 2;
            camera.direction(1) = 1;
            return executor.Launch(camera_metric, 0, camera);
        });
        EXPECT_TRUE(executor.Step(ray, flat, config, coupled, comparison));
        EXPECT_TRUE(standalone.get());
        const auto mixed_after = executor.Statistics();
        EXPECT_EQ(mixed_after.coalescing_traces_ready, trace_ready);
        EXPECT_EQ(mixed_after.coalescing_timeouts - mixed_before.coalescing_timeouts,
                  mixed_after.batches - mixed_before.batches);
    }
    const auto accepted = ray;
    const auto submissions = executor.Statistics().interval_batches;
    cancelled = true;
    EXPECT_FALSE(executor.Step(ray, flat, config, coupled, comparison));
    EXPECT_EQ(executor.Statistics().interval_batches, submissions);
    for (int axis = 0; axis < 4; ++axis) {
        EXPECT_EQ(ray.position(axis), accepted.position(axis));
        EXPECT_EQ(ray.velocity(axis), accepted.velocity(axis));
    }
    const auto timing = executor.Statistics();
    std::uint64_t observed_batches = 0, observed_rows = 0;
    for (std::size_t rows = 0; rows < timing.batch_row_counts.size(); ++rows) {
        observed_batches += timing.batch_row_counts[rows];
        observed_rows += rows * timing.batch_row_counts[rows];
    }
    EXPECT_EQ(timing.batch_row_counts.size(), compute->Capacity() + 1);
    EXPECT_EQ(timing.batch_row_counts.front(), 0U);
    EXPECT_EQ(observed_batches, timing.batches);
    EXPECT_EQ(observed_rows, timing.interval_rows + timing.camera_rows);
    EXPECT_EQ(timing.interval_rows, timing.accepted_intervals + timing.rejected_intervals);
    EXPECT_EQ(timing.acceleration_calls, timing.accepted_intervals);
    EXPECT_LE(timing.full_batches, timing.batches);
    EXPECT_LE(timing.coalescing_timeouts, timing.batches);
    EXPECT_LE(timing.coalescing_timeouts, timing.coalescing_underfilled);
    EXPECT_LE(timing.coalescing_underfilled, timing.batches);
    EXPECT_EQ(timing.coalescing_stopped, 0U);
    EXPECT_GT(timing.coalescing_traces_ready, 0U);
    EXPECT_LE(timing.coalescing_timeouts + timing.coalescing_traces_ready, timing.batches);
    RecordProperty("coordinator_batches", std::to_string(timing.batches));
    RecordProperty("coordinator_full_batches", std::to_string(timing.full_batches));
    RecordProperty("coordinator_rows", std::to_string(observed_rows));
    RecordProperty("coalescing_timeouts", std::to_string(timing.coalescing_timeouts));
    RecordProperty("coalescing_underfilled", std::to_string(timing.coalescing_underfilled));
    RecordProperty("coalescing_traces_ready", std::to_string(timing.coalescing_traces_ready));
    RecordProperty("coalescing_wait_ms", std::to_string(timing.coalescing_wait_ms));
    RecordProperty("maximum_coalescing_wait_ms", std::to_string(timing.maximum_coalescing_wait_ms));
    RecordProperty("coordinator_execute_ms", std::to_string(timing.execute_ms));
    RecordProperty("worker_acceleration_calls", std::to_string(timing.acceleration_calls));
    RecordProperty("worker_acceleration_ms", std::to_string(timing.acceleration_ms));

    // Deterministic hard-guard timings at the real ComputeDevice seam. This is
    // a controller witness, not a measured Radeon/software-driver duration.
    TransferProbeDevice guard_probe(*device);
    auto guard_created = RetainedCompute::Create(guard_probe, 4, false, 0);
    ASSERT_TRUE(guard_created) << guard_created.error().Description();
    auto& guard_compute = **guard_created;
    const auto guard_allocation = device->BufferAllocationBytes();
    EXPECT_LE(guard_allocation, 8ull * 1024 * 1024);
    ASSERT_EQ(guard_probe.allocations.size(), 12U);
    const auto transport_input = guard_probe.allocations[2].handle;
    const auto endpoint_input = guard_probe.allocations[4].handle;
    guard_probe.capture_output = guard_probe.allocations[5].handle;
    const auto endpoint_row_bytes = guard_probe.allocations[5].bytes / guard_compute.Capacity();
    enum class GuardTiming {
        Normal,
        Reducible,
        Irreducible,
        RetryIrreducible,
        FailedPair,
        FailedRetry,
        Cancelled,
        CancelledRetry
    };
    GuardTiming guard_timing = GuardTiming::Normal;
    std::size_t paired_endpoint_calls = 0;
    std::atomic<bool> guard_cancelled{false};
    guard_probe.submission_ms = [&](BufferHandle input, std::uint32_t rows) {
        if (input.value == endpoint_input.value && rows == 2) {
            ++paired_endpoint_calls;
            if (guard_timing == GuardTiming::Cancelled) guard_cancelled = true;
            if (guard_timing == GuardTiming::FailedPair)
                guard_probe.fail_next = TransferProbeDevice::Failure::Read;
            return guard_timing == GuardTiming::Normal ? 800.0 : 1200.0;
        }
        if (input.value == endpoint_input.value && rows == 1 && paired_endpoint_calls == 3) {
            if (guard_timing == GuardTiming::FailedRetry)
                guard_probe.fail_next = TransferProbeDevice::Failure::Read;
            if (guard_timing == GuardTiming::CancelledRetry) guard_cancelled = true;
        }
        if (input.value == transport_input.value && rows == 1 &&
            (guard_timing == GuardTiming::Irreducible ||
             (guard_timing == GuardTiming::RetryIrreducible && paired_endpoint_calls == 3)))
            return 1100.0;
        return 600.0;
    };
    const auto make_ray = [] {
        sirius::core::Lightray result{};
        result.position(1) = 5;
        result.velocity(0) = -1;
        result.velocity(1) = 1;
        result.step_size = 1;
        return result;
    };
    const auto make_columns = [] {
        sirius::core::Rk45CoupledState result;
        result.length_scale = result.frequency_scale = 1;
        result.tolerance = 1e-9;
        result.variations[0].derivative(2) = .001;
        result.variations[1].derivative(3) = .001;
        result.variations[2].displacement(2) = 1;
        result.variations[3].displacement(3) = 1;
        return result;
    };
    const auto seed = make_ray();
    const auto seed_columns = make_columns();
    RetainedInitializeInput seed_input;
    for (std::size_t axis = 0; axis < 4; ++axis) {
        seed_input.values[axis] = RetainedValue::FromDouble(0);
        seed_input.values[4 + axis] =
            RetainedValue::FromDouble(seed.position(static_cast<int>(axis)));
        seed_input.values[8 + axis] =
            RetainedValue::FromDouble(seed.velocity(static_cast<int>(axis)));
        for (std::size_t column = 0; column < 4; ++column) {
            seed_input.values[12 + 8 * column + axis] = RetainedValue::FromDouble(
                seed_columns.variations[column].displacement(static_cast<int>(axis)));
            seed_input.values[16 + 8 * column + axis] = RetainedValue::FromDouble(
                seed_columns.variations[column].derivative(static_cast<int>(axis)));
        }
    }
    seed_input.values[44] = RetainedValue::FromDouble(1);
    const auto phase = guard_compute.Initialize({&seed_input, 1});
    ASSERT_TRUE(phase) << phase.error().Description();
    ASSERT_TRUE(phase->front().valid);
    RetainedEndpointInput projection_input;
    std::copy_n(seed_input.values.begin(), 4, projection_input.values.begin());
    std::copy(phase->front().phase.begin(), phase->front().phase.end(),
              projection_input.values.begin() + 4);
    projection_input.values[44] = RetainedValue::FromDouble(1);
    const auto initial = guard_compute.Endpoint({&projection_input, 1});
    ASSERT_TRUE(initial) << initial.error().Description();
    ASSERT_TRUE(initial->front().valid);
    RetainedIntervalInput serial_input;
    std::copy_n(seed_input.values.begin(), 4, serial_input.metric.begin());
    serial_input.start = initial->front();
    serial_input.chart = 1;
    serial_input.interval = 1;
    serial_input.control = {config, seed_columns.length_scale, seed_columns.frequency_scale,
                            seed_columns.tolerance, seed_columns.column_scale};
    guard_probe.readbacks.clear();
    const auto serial = AttemptRetainedIntervals(guard_compute, {&serial_input, 1}, 1);
    ASSERT_TRUE(serial) << serial.error().Description();
    ASSERT_TRUE(serial->front().admissible);
    const auto serial_rows = guard_probe.readbacks;
    ASSERT_EQ(serial_rows.size(), 6U);
    const auto serial_feedback = guard_compute.TakeSubmissionFeedback();
    EXPECT_EQ(serial_feedback.peak_ms, 600);
    EXPECT_EQ(serial_feedback.peak_rows, 1U);
    EXPECT_EQ(serial_feedback.maximum_rows, 1U);
    EXPECT_EQ(serial_feedback.maximum_one_row_ms, 600);
    const auto consumed = guard_compute.TakeSubmissionFeedback();
    EXPECT_EQ(consumed.peak_ms, 0);
    EXPECT_EQ(consumed.maximum_rows, 0U);
    EXPECT_EQ(consumed.maximum_one_row_ms, 0);
    const auto compare_row = [&](const std::vector<std::byte>& actual, std::size_t row,
                                 const std::vector<std::byte>& expected) {
        ASSERT_EQ(expected.size(), endpoint_row_bytes);
        ASSERT_GE(actual.size(), (row + 1) * endpoint_row_bytes);
        const auto first = actual.begin() + row * endpoint_row_bytes;
        EXPECT_TRUE(std::equal(expected.begin(), expected.end(), first));
    };
    for (const auto mode :
         {GuardTiming::Normal, GuardTiming::Reducible, GuardTiming::Irreducible,
          GuardTiming::RetryIrreducible, GuardTiming::FailedPair, GuardTiming::FailedRetry,
          GuardTiming::Cancelled, GuardTiming::CancelledRetry}) {
        SCOPED_TRACE(static_cast<unsigned>(mode));
        guard_timing = mode;
        paired_endpoint_calls = 0;
        guard_cancelled = false;
        guard_probe.readbacks.clear();
        auto guarded_ray = make_ray();
        auto guarded_columns = make_columns();
        sirius::core::Rk45CoupledComparison guarded_comparison;
        RetainedTraceExecutor guarded_executor(
            guard_compute, [&] { return guard_cancelled.load(); }, 1000);
        const bool completed =
            guarded_executor.Step(guarded_ray, flat, config, guarded_columns, guarded_comparison);
        const auto guarded_stats = guarded_executor.Statistics();
        const bool recovered = mode == GuardTiming::Normal || mode == GuardTiming::Reducible;
        EXPECT_EQ(completed, recovered);
        EXPECT_EQ(guarded_stats.accepted_intervals, recovered ? 1U : 0U);
        EXPECT_EQ(guarded_stats.rejected_intervals, recovered ? 0U : 1U);
        EXPECT_EQ(guarded_stats.interval_rows, 1U);
        EXPECT_EQ(guarded_stats.paired_projection_retries,
                  mode == GuardTiming::Reducible || mode == GuardTiming::RetryIrreducible ||
                          mode == GuardTiming::FailedRetry || mode == GuardTiming::CancelledRetry
                      ? 1U
                      : 0U);
        const auto completed_stages = mode == GuardTiming::FailedPair ? 0U
                                      : mode == GuardTiming::Reducible ||
                                              mode == GuardTiming::RetryIrreducible ||
                                              mode == GuardTiming::CancelledRetry
                                          ? 42U
                                          : 21U;
        EXPECT_EQ(guarded_columns.central_stages, completed_stages);
        EXPECT_EQ(guarded_columns.variation_stages, completed_stages);
        EXPECT_EQ(guard_probe.allocations.size(), 12U);
        EXPECT_EQ(device->BufferAllocationBytes(), guard_allocation);
        if (recovered) {
            EXPECT_FALSE(guarded_executor.Error());
            ASSERT_EQ(guard_probe.readbacks.size(), mode == GuardTiming::Normal ? 4U : 10U);
            for (std::size_t part = 0; part < 3; ++part)
                for (std::size_t order = 0; order < 2; ++order)
                    ASSERT_NO_FATAL_FAILURE(compare_row(guard_probe.readbacks[part + 1], order,
                                                        serial_rows[part * 2 + order]));
            EXPECT_EQ(guarded_columns.central_stages, mode == GuardTiming::Normal ? 21U : 42U);
            for (int axis = 0; axis < 4; ++axis) {
                EXPECT_EQ(guarded_ray.position(axis), serial->front().full.physical[axis].Center());
                EXPECT_EQ(guarded_ray.velocity(axis),
                          serial->front().full.physical[4 + axis].Center());
            }
            if (mode == GuardTiming::Reducible) {
                for (std::size_t row = 0; row < serial_rows.size(); ++row)
                    ASSERT_NO_FATAL_FAILURE(
                        compare_row(guard_probe.readbacks[4 + row], 0, serial_rows[row]));
                EXPECT_EQ(guarded_stats.safety_fallbacks, 1U);
                const auto after_recovery = guard_compute.Statistics();
                ASSERT_TRUE(guarded_executor.Step(guarded_ray, flat, config, guarded_columns,
                                                  guarded_comparison));
                EXPECT_EQ(paired_endpoint_calls, 3U);
                EXPECT_EQ(guard_compute.Statistics()[2].submissions - after_recovery[2].submissions,
                          6U);
                EXPECT_EQ(guarded_executor.Statistics().paired_projection_retries, 1U);
                EXPECT_FALSE(guarded_executor.Error());
                RecordProperty("paired_guard_recovery_retries",
                               std::to_string(guarded_stats.paired_projection_retries));
                RecordProperty("paired_guard_recovery_charged_stages",
                               std::to_string(guarded_columns.central_stages - 21));
            }
        } else {
            const bool fatal =
                mode != GuardTiming::Cancelled && mode != GuardTiming::CancelledRetry;
            EXPECT_EQ(guarded_executor.Error().has_value(), fatal);
            EXPECT_EQ(guarded_ray.terminated, 3);
            for (int axis = 0; axis < 4; ++axis) {
                EXPECT_EQ(guarded_ray.position(axis), seed.position(static_cast<int>(axis)));
                EXPECT_EQ(guarded_ray.velocity(axis), seed.velocity(static_cast<int>(axis)));
                for (std::size_t column = 0; column < 4; ++column) {
                    EXPECT_EQ(guarded_columns.variations[column].displacement(axis),
                              seed_columns.variations[column].displacement(static_cast<int>(axis)));
                    EXPECT_EQ(guarded_columns.variations[column].derivative(axis),
                              seed_columns.variations[column].derivative(static_cast<int>(axis)));
                }
            }
            EXPECT_EQ(guarded_ray.proper_time, seed.proper_time);
            EXPECT_EQ(guarded_ray.step_size, seed.step_size);
            EXPECT_EQ(guard_probe.readbacks.size(),
                      mode == GuardTiming::RetryIrreducible || mode == GuardTiming::CancelledRetry
                          ? 10U
                      : mode == GuardTiming::FailedPair ? 1U
                                                        : 4U);
            if (mode == GuardTiming::FailedPair || mode == GuardTiming::FailedRetry) {
                ASSERT_TRUE(guarded_executor.Error());
                EXPECT_EQ(guarded_executor.Error()->operation(), "retained transfer probe read");
            }
            if (fatal) {
                const auto calls = guard_probe.dispatch_calls;
                EXPECT_FALSE(guarded_executor.Step(guarded_ray, flat, config, guarded_columns,
                                                   guarded_comparison));
                EXPECT_EQ(guard_probe.dispatch_calls, calls);
                EXPECT_EQ(guarded_columns.central_stages, completed_stages);
                EXPECT_EQ(guarded_columns.variation_stages, completed_stages);
            }
        }
    }
    RecordProperty("paired_guard_timing_scope",
                   "injected valid seam observations only: single rows 600ms, paired endpoints "
                   "800/1200ms, transport overshoot 1100ms, hard bound 1000ms; no sleeps and no "
                   "hardware-duration claim; real device outputs compare bitwise to serialization");

#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST_F(RetainedComputeTest, RejectedStepRowsCannotExposeOldOrPartialCandidates) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    const auto good =
        std::bit_cast<RetainedStepInput>(sirius::test::retained_transport::cases.front().input);
    std::vector<RetainedStepInput> inputs(9, good);
    auto outputs = compute->Step(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    for (const auto& output : *outputs) ASSERT_TRUE(output.valid);
    inputs[0].values[4].high = std::numeric_limits<float>::quiet_NaN();
    inputs[1].values[5].low = std::numeric_limits<float>::infinity();
    inputs[2].values[6].radius = -1;
    inputs[3].values[7].valid = 0;
    inputs[4].values[44] = RetainedValue::FromDouble(0);
    inputs[5].values[45] = RetainedValue::FromDouble(0);
    inputs[6].values[45] = RetainedValue::FromDouble(-1);
    inputs[7].values[3] = RetainedValue::FromDouble(.01);
    for (std::size_t i = 5; i < 8; ++i) inputs[8].values[i] = RetainedValue::FromDouble(0);
    outputs = compute->Step(inputs);
    ASSERT_TRUE(outputs) << outputs.error().Description();
    for (const auto& output : *outputs) {
        EXPECT_FALSE(output.valid);
        for (const auto* record : {&output.fifth, &output.fourth, &output.increment, &output.error})
            for (const auto& value : *record) EXPECT_EQ(value.valid, 0U);
    }
    EXPECT_EQ(outputs->back().stages, 1U);
    const auto recovered = compute->Step(std::span(&good, 1));
    ASSERT_TRUE(recovered) << recovered.error().Description();
    ASSERT_EQ(recovered->size(), 1U);
    EXPECT_TRUE(StepAgrees(recovered->front(), sirius::test::retained_transport::cases.front()));
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST_F(RetainedComputeTest, Fp64ProductsPreserveIndependentScienceOrDeclineUnsupportedDevices) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    auto wide = RetainedCompute::Create(*device, 24, true);
    if (!device->Info().supports_fp64 || !device->Info().rounds_fp64_to_nearest) {
        ASSERT_FALSE(wide);
        EXPECT_NE(wide.error().detail().find("binary64"), std::string::npos);
        return;
    }
    ASSERT_TRUE(wide) << wide.error().Description();
    std::vector<RetainedRayCameraInput> cameras;
    for (const auto& fixture : sirius::test::retained_ray_camera::cases)
        cameras.push_back(std::bit_cast<RetainedRayCameraInput>(fixture.input));
    const auto launched = (*wide)->RayCamera(cameras);
    ASSERT_TRUE(launched) << launched.error().Description();
    for (std::size_t row = 0; row < cameras.size(); ++row)
        EXPECT_TRUE(CameraAgrees((*launched)[row], sirius::test::retained_ray_camera::cases[row]));
    std::vector<RetainedStepInput> steps;
    std::vector<RetainedEndpointInput> endpoints;
    for (const auto& fixture : sirius::test::retained_transport::cases) {
        steps.push_back(std::bit_cast<RetainedStepInput>(fixture.input));
        RetainedEndpointInput endpoint;
        std::copy_n(steps.back().values.begin(), 45, endpoint.values.begin());
        endpoints.push_back(endpoint);
    }
    const auto stepped = (*wide)->Step(steps);
    ASSERT_TRUE(stepped) << stepped.error().Description();
    for (std::size_t row = 0; row < steps.size(); ++row)
        EXPECT_TRUE(StepAgrees((*stepped)[row], sirius::test::retained_transport::cases[row]));
    const auto projected = (*wide)->Endpoint(endpoints);
    ASSERT_TRUE(projected) << projected.error().Description();
    std::vector<RetainedIntervalInput> intervals(steps.size());
    std::vector<RetainedInitializeInput> physical(steps.size());
    for (std::size_t row = 0; row < steps.size(); ++row) {
        ASSERT_TRUE((*projected)[row].valid);
        auto& input = intervals[row];
        std::copy_n(steps[row].values.begin(), 4, input.metric.begin());
        input.start = (*projected)[row];
        input.chart = steps[row].values[44].Center();
        input.interval = steps[row].values[45].Center();
        input.control.length_scale = 1;
        input.control.frequency_scale = 1;
        input.control.tolerance = 1e-4 / (4 * 30000);
        input.control.integrator.min_step = 1e-15f;
        input.control.integrator.max_step = 1;
        input.control.integrator.abs_tolerance = 1e-9f;
        input.control.integrator.rel_tolerance = 1e-9f;
        std::copy(input.metric.begin(), input.metric.end(), physical[row].values.begin());
        std::copy(input.start.physical.begin(), input.start.physical.end(),
                  physical[row].values.begin() + 4);
        physical[row].values[44] = steps[row].values[44];
    }
    const auto initialized = (*wide)->Initialize(physical);
    ASSERT_TRUE(initialized) << initialized.error().Description();
    for (const auto& row : *initialized) EXPECT_TRUE(row.valid);
    const auto narrow_intervals = AttemptRetainedIntervals(*compute, intervals);
    const auto wide_intervals = AttemptRetainedIntervals(**wide, intervals);
    ASSERT_TRUE(narrow_intervals) << narrow_intervals.error().Description();
    ASSERT_TRUE(wide_intervals) << wide_intervals.error().Description();
    for (std::size_t row = 0; row < intervals.size(); ++row) {
        SCOPED_TRACE(sirius::test::retained_transport::cases[row].name);
        ASSERT_TRUE((*wide_intervals)[row].admissible)
            << sirius::core::CoupledStepFailureName((*wide_intervals)[row].failure)
            << " error=" << (*wide_intervals)[row].error_ratio;
        ASSERT_TRUE((*narrow_intervals)[row].admissible);
        EXPECT_LT(
            RetainedPhysicalError((*wide_intervals)[row].full.physical,
                                  (*narrow_intervals)[row].full.physical, intervals[row].control),
            1e-4);
    }
    const auto pair_inputs = std::span<const RetainedIntervalInput>(intervals).first(2);
    const auto serialized_before = (*wide)->Statistics();
    const auto serialized = AttemptRetainedIntervals(**wide, pair_inputs, 2);
    ASSERT_TRUE(serialized) << serialized.error().Description();
    EXPECT_EQ((*wide)->Statistics()[2].submissions - serialized_before[2].submissions, 6U);
    const auto paired_before = (*wide)->Statistics();
    const auto paired = AttemptRetainedIntervals(**wide, pair_inputs, 4);
    ASSERT_TRUE(paired) << paired.error().Description();
    EXPECT_EQ((*wide)->Statistics()[2].submissions - paired_before[2].submissions, 3U);
    for (std::size_t row = 0; row < pair_inputs.size(); ++row) {
        SCOPED_TRACE(row);
        ASSERT_TRUE((*paired)[row].admissible);
        EXPECT_TRUE(IntervalBitsAgree((*paired)[row], (*serialized)[row]));
    }
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST(RetainedValue, Binary64InputsKeepTheirRepresentationResidual) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    for (const double input : {0.0, -0.0, 1.0 / 3.0, -1.0 / 7.0, 0x1.123456789abcdp80,
                               0x1.123456789abcdp-120, 0x1p-149, 0x1p-150, 0x1p-1022}) {
        const auto value = RetainedValue::FromDouble(input);
        ASSERT_TRUE(value.IsRepresented());
        EXPECT_TRUE(Encloses(value, input, 0));
    }
    EXPECT_FALSE(
        RetainedValue::FromDouble(std::numeric_limits<double>::infinity()).IsRepresented());
    EXPECT_FALSE(RetainedValue::FromDouble(0x1p121).IsRepresented());
    EXPECT_FALSE((RetainedValue{1, 1, 0, 0, 1}.IsRepresented()));

    std::array<RetainedValue, 40> first;
    first.fill(RetainedValue::FromDouble(0));
    first[8] = {1, 0x1.000002p-35f, 0x1p-120f, 0, 1};
    ASSERT_TRUE(first[8].IsRepresented());
    auto second = first;
    second[8].tail = 0;
    RetainedIntervalControl control;
    control.length_scale = control.frequency_scale = 1;
    control.tolerance = 1e-40;
    const double expected = 0x1p-120 / (control.tolerance * (1 + first[8].Center()));
    EXPECT_DOUBLE_EQ(RetainedPhysicalError(first, second, control), expected);
    EXPECT_DOUBLE_EQ(RetainedPhysicalError(second, first, control), expected);
    EXPECT_EQ(RetainedPhysicalError(first, first, control), 0);
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}
}  // namespace
