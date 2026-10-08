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
#include <format>
#include <functional>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

#ifdef SIRIUS_HAS_RETAINED_COMPUTE
#include "retained_kernels.h"
#endif

#if defined(SIRIUS_HAS_RETAINED_COMPUTE) && defined(SIRIUS_TEST_HAS_RETAINED_CAMERA)
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

#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
namespace sirius::backend {
struct RetainedTraceExecutorTestPeer {
    static std::array<std::size_t, 3> WaitForQueued(RetainedTraceExecutor& executor) {
        std::unique_lock lock(executor.mutex_);
        executor.available_.wait(lock, [&] { return !executor.requests_.empty(); });
        return {executor.requests_.size(), executor.queued_registered_,
                executor.active_traces_.size()};
    }
    static void WaitForQueued(RetainedTraceExecutor& executor, std::size_t rows) {
        std::unique_lock lock(executor.mutex_);
        executor.available_.wait(lock, [&] { return executor.requests_.size() >= rows; });
    }
    static std::array<std::size_t, 3> QueueState(RetainedTraceExecutor& executor) {
        std::lock_guard lock(executor.mutex_);
        return {executor.requests_.size(), executor.queued_registered_,
                executor.active_traces_.size()};
    }
};
}  // namespace sirius::backend
#endif

namespace {
using namespace sirius::backend;

// Arithmetic admission must finish before any external shader or buffer work.
// The sentinel gives that boundary an observable positive path without a device.
class AdmissionDevice final : public ComputeDevice {
  public:
    DeviceInfo info;
    unsigned kernel_calls = 0, buffer_calls = 0, query_calls = 0;
    unsigned query_failure_call = 0;
    std::uint64_t query_padding = 0;
    std::optional<std::uint64_t> fixed_requirement;
    std::vector<std::uint64_t> queried_spans;
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
    sirius::base::Expected<std::uint64_t> RequiredBufferAllocationBytes(
        std::uint64_t bytes, BufferUsage usage) override {
        ++query_calls;
        queried_spans.push_back(bytes);
        if (usage != BufferUsage::kStorage || query_calls == query_failure_call)
            return sirius::base::Fail(sirius::base::ErrorDomain::kDevice,
                                      "allocation query sentinel", "injected query failure");
        return fixed_requirement ? *fixed_requirement : bytes + query_padding;
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
            EXPECT_EQ(device.query_calls, 0U);
            if (admitted) {
                EXPECT_EQ(device.kernel_has_float, missing != 1 && missing != 2 && missing != 5);
            }
            EXPECT_EQ(created.error().domain(), admitted ? sirius::base::ErrorDomain::kKernel
                                                         : sirius::base::ErrorDomain::kDevice);
            EXPECT_EQ(created.error().operation(),
                      admitted ? "admission sentinel" : "create retained compute stages");
        }
    }
    // A driver's allocation requirement can exceed every logical shader span.
    // These fixed layout totals are independent of the production planner.
    for (const auto& [capacity, logical] : std::array<std::pair<std::size_t, std::uint64_t>, 3>{
             {{1, 671396}, {24, 3301216}, {64, 7874816}}}) {
        AdmissionDevice padded;
        padded.query_padding = 128;
        const auto required = RetainedCompute::RequiredAllocationBytes(padded, capacity);
        ASSERT_TRUE(required) << required.error().Description();
        EXPECT_EQ(RetainedCompute::RequiredBufferBytes(capacity), logical);
        EXPECT_EQ(*required, logical + 14 * 128);
        EXPECT_EQ(padded.query_calls, 14U);
        EXPECT_EQ(padded.kernel_calls, 0U);
        EXPECT_EQ(padded.buffer_calls, 0U);
        EXPECT_EQ(padded.BufferAllocationBytes(), 0U);
        if (capacity == 24) {
            EXPECT_EQ(padded.queried_spans,
                      (std::vector<std::uint64_t>{122728, 284544, 74468, 451104, 84864, 367584,
                                                  181492, 410784, 61292, 193344, 145452, 308736,
                                                  216520, 398304}));
        }
    }
    AdmissionDevice invalid;
    EXPECT_FALSE(RetainedCompute::RequiredAllocationBytes(invalid, 0));
    EXPECT_FALSE(RetainedCompute::RequiredAllocationBytes(invalid, 65536));
    EXPECT_EQ(invalid.query_calls, 0U);
    AdmissionDevice failing;
    failing.query_failure_call = 3;
    const auto failure = RetainedCompute::RequiredAllocationBytes(failing, 24);
    ASSERT_FALSE(failure);
    EXPECT_EQ(failure.error().operation(), "allocation query sentinel");
    EXPECT_EQ(failing.query_calls, 3U);
    EXPECT_EQ(failing.kernel_calls + failing.buffer_calls, 0U);
    AdmissionDevice overflow;
    overflow.fixed_requirement = std::numeric_limits<std::uint64_t>::max();
    const auto too_large = RetainedCompute::RequiredAllocationBytes(overflow, 24);
    ASSERT_FALSE(too_large);
    EXPECT_NE(too_large.error().detail().find("overflows"), std::string::npos);
    EXPECT_EQ(overflow.query_calls, 2U);
    EXPECT_EQ(overflow.kernel_calls + overflow.buffer_calls, 0U);
    AdmissionDevice undersized;
    undersized.fixed_requirement = 1;
    const auto too_small = RetainedCompute::RequiredAllocationBytes(undersized, 24);
    ASSERT_FALSE(too_small);
    EXPECT_NE(too_small.error().detail().find("smaller"), std::string::npos);
    EXPECT_EQ(undersized.query_calls, 1U);
    EXPECT_EQ(undersized.kernel_calls + undersized.buffer_calls, 0U);
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

#ifdef SIRIUS_HAS_RETAINED_COMPUTE
// Host-only model of the genuine ComputeDevice boundary. It does not claim
// shader execution; the separately selectable device control checks that join.
class PreparationProbeDevice final : public ComputeDevice {
  public:
    struct Write {
        BufferHandle buffer;
        std::size_t bytes;
        std::uint32_t header;
    };
    DeviceInfo info{.kind = DeviceKind::kSoftware};
    std::vector<std::vector<std::byte>> buffers;
    std::vector<Write> writes;
    std::vector<std::uint32_t> kernels;
    std::vector<std::vector<std::uint32_t>> loaded_codes;
    std::vector<std::size_t> failed_writes;
    std::function<void()> after_write, after_dispatch;
    DispatchTiming observation{.submit_wait_ms = 5000,
                               .pipeline_setup_ms = 1,
                               .command_setup_ms = 2,
                               .cleanup_ms = 3,
                               .total_ms = 5006,
                               .pipeline_created = true};
    bool fail_dispatch = false;
    bool independent_pair = false;
    std::uint32_t expected_groups_x = 1;
    unsigned loads = 0, reads = 0;
    const DeviceInfo& Info() const noexcept override { return info; }
    sirius::base::Expected<KernelHandle> LoadKernel(std::span<const std::uint32_t> code) override {
        loaded_codes.emplace_back(code.begin(), code.end());
        return KernelHandle{loads++};
    }
    sirius::base::Expected<BufferHandle> CreateBuffer(std::uint64_t bytes,
                                                      BufferUsage usage) override {
        EXPECT_EQ(usage, BufferUsage::kStorage);
        const BufferHandle result{static_cast<std::uint32_t>(buffers.size())};
        buffers.emplace_back(bytes, std::byte{0xa5});
        return result;
    }
    sirius::base::Expected<std::uint64_t> RequiredBufferAllocationBytes(std::uint64_t bytes,
                                                                        BufferUsage) override {
        return bytes;
    }
    sirius::base::Expected<void> WriteBuffer(BufferHandle buffer,
                                             std::span<const std::byte> data) override {
        if (buffer.value >= buffers.size() || data.size() > buffers[buffer.value].size())
            return Failure("range");
        std::uint32_t header = 0;
        if (data.size() >= sizeof(header)) std::memcpy(&header, data.data(), sizeof(header));
        writes.push_back({buffer, data.size(), header});
        if (std::ranges::find(failed_writes, writes.size()) != failed_writes.end())
            return Failure("write");
        std::copy(data.begin(), data.end(), buffers[buffer.value].begin());
        if (after_write) after_write();
        return {};
    }
    sirius::base::Expected<void> ReadBuffer(BufferHandle buffer,
                                            std::span<std::byte> data) override {
        ++reads;
        if (buffer.value >= buffers.size() || data.size() > buffers[buffer.value].size())
            return Failure("range");
        std::copy_n(buffers[buffer.value].begin(), data.size(), data.begin());
        return {};
    }
    sirius::base::Expected<void> Dispatch(KernelHandle kernel,
                                          std::span<const BufferHandle> bindings, std::uint32_t x,
                                          std::uint32_t y, std::uint32_t z,
                                          DispatchTiming* timing) override {
        EXPECT_EQ(bindings.size(), 2U);
        EXPECT_EQ(x, expected_groups_x);
        EXPECT_EQ(y, 1U);
        EXPECT_EQ(z, 1U);
        kernels.push_back(kernel.value);
        if (timing) *timing = observation;
        if (after_dispatch) after_dispatch();
        if (fail_dispatch) return Failure("dispatch");
        std::uint32_t rows = 0;
        std::memcpy(&rows, buffers[bindings[0].value].data(), sizeof(rows));
        // A physical control returns a declined, zero-stage row; only its
        // actual submission/feedback ownership is tested against this model.
        if (rows != 0) std::ranges::fill(buffers[bindings[1].value], std::byte{0});
        return {};
    }
    sirius::base::Expected<void> SetBufferAllocationLimit(std::uint64_t) override { return {}; }
    bool SupportsIndependentPair() const noexcept override { return independent_pair; }
    sirius::base::Expected<void> DispatchIndependentPair(
        const std::array<sirius::backend::ComputeDispatch, 2>& commands,
        sirius::backend::IndependentPairTiming* timing) override {
        if (!independent_pair) return ComputeDevice::DispatchIndependentPair(commands, timing);
        for (const auto& command : commands) {
            auto status = Dispatch(command.kernel, command.buffers, command.groups_x,
                                   command.groups_y, command.groups_z, nullptr);
            if (!status) return status;
        }
        if (timing) *timing = {.combined = observation, .pipeline_creations = 2};
        return {};
    }
    std::uint64_t BufferAllocationBytes() const noexcept override {
        std::uint64_t bytes = 0;
        for (const auto& buffer : buffers) bytes += buffer.size();
        return bytes;
    }

  private:
    static sirius::base::Expected<void> Failure(const char* operation) {
        return sirius::base::Fail(sirius::base::ErrorDomain::kDevice,
                                  std::string("preparation probe ") + operation,
                                  "injected returned failure");
    }
};

void ExpectPhysicalStatsEqual(
    const std::array<RetainedCompute::StageStats, RetainedCompute::kStageCount>& actual,
    const std::array<RetainedCompute::StageStats, RetainedCompute::kStageCount>& expected) {
    for (std::size_t i = 0; i < actual.size(); ++i) {
        SCOPED_TRACE(i);
        EXPECT_EQ(actual[i].submissions, expected[i].submissions);
        EXPECT_EQ(actual[i].submit_wait_ms, expected[i].submit_wait_ms);
        EXPECT_EQ(actual[i].maximum_submit_wait_ms, expected[i].maximum_submit_wait_ms);
        EXPECT_EQ(actual[i].pipeline_setup_ms, expected[i].pipeline_setup_ms);
        EXPECT_EQ(actual[i].command_setup_ms, expected[i].command_setup_ms);
        EXPECT_EQ(actual[i].cleanup_ms, expected[i].cleanup_ms);
        EXPECT_EQ(actual[i].dispatch_total_ms, expected[i].dispatch_total_ms);
        EXPECT_EQ(actual[i].write_buffer_ms, expected[i].write_buffer_ms);
        EXPECT_EQ(actual[i].read_buffer_ms, expected[i].read_buffer_ms);
        EXPECT_EQ(actual[i].write_buffer_bytes, expected[i].write_buffer_bytes);
        EXPECT_EQ(actual[i].read_buffer_bytes, expected[i].read_buffer_bytes);
        EXPECT_EQ(actual[i].pipeline_creations, expected[i].pipeline_creations);
        EXPECT_EQ(actual[i].target_overshoots, expected[i].target_overshoots);
    }
}
#endif

#ifdef SIRIUS_HAS_RETAINED_COMPUTE
TEST(RetainedComputeAdmission, FmaSelectsOnlyNativeWideProductsAndPreservesAllocation) {
    using namespace sirius::backend::retained_program;
    for (unsigned mask = 0; mask < 16; ++mask) {
        for (const bool wide : {false, true}) {
            SCOPED_TRACE(mask);
            SCOPED_TRACE(wide);
            PreparationProbeDevice control, candidate;
            control.info.supports_fp64 = control.info.rounds_fp64_to_nearest = true;
            control.info.preserves_fp32_denormals = (mask & 1u) != 0;
            control.info.rounds_fp32_to_nearest = (mask & 2u) != 0;
            control.info.preserves_fp32_signed_zero_inf_nan = (mask & 4u) != 0;
            candidate.info = control.info;
            candidate.info.fma_fp32_enabled = (mask & 8u) != 0;
            auto baseline = RetainedCompute::Create(control, 2, wide);
            auto created = RetainedCompute::Create(candidate, 2, wide);
            ASSERT_TRUE(baseline);
            ASSERT_TRUE(created);
            ASSERT_EQ(control.loaded_codes.size(), RetainedCompute::kStageCount);
            ASSERT_EQ(candidate.loaded_codes.size(), RetainedCompute::kStageCount);
            const bool portable = (mask & 3u) != 3u;
            const auto expected = [portable, wide](std::span<const std::uint32_t> native,
                                                   std::span<const std::uint32_t> native_wide,
                                                   std::span<const std::uint32_t> integer,
                                                   std::span<const std::uint32_t> integer_wide) {
                return portable ? (wide ? integer_wide : integer) : (wide ? native_wide : native);
            };
            std::array<std::span<const std::uint32_t>, RetainedCompute::kStageCount> stages{
                expected(kCameraShader, kCameraFp64Shader, kCameraPortableShader,
                         kCameraPortableFp64Shader),
                expected(kTransportShader, kTransportFp64Shader, kTransportPortableShader,
                         kTransportPortableFp64Shader),
                expected(kEndpointShader, kEndpointFp64Shader, kEndpointPortableShader,
                         kEndpointPortableFp64Shader),
                expected(kDenseShader, kDenseFp64Shader, kDensePortableShader,
                         kDensePortableFp64Shader),
                expected(kInitializeShader, kInitializeFp64Shader, kInitializePortableShader,
                         kInitializePortableFp64Shader),
                expected(kRayCameraShader, kRayCameraFp64Shader, kRayCameraPortableShader,
                         kRayCameraPortableFp64Shader),
                expected(kDopriPhaseShader, kDopriPhaseFp64Shader, kDopriPhasePortableShader,
                         kDopriPhasePortableFp64Shader)};
            if ((mask & 3u) == 2u) {
                stages[1] = std::span(kTransportPortableNormalSumShader);
                stages[2] = std::span(kEndpointPortableNormalSumShader);
            }
            const bool eligible = wide && mask == 15u;
            for (std::size_t stage = 0; stage < RetainedCompute::kStageCount; ++stage) {
                // The baseline excludes FMA and independently checks native,
                // pure-integer fallback and RTE32-only Transport/Endpoint selection.
                EXPECT_EQ(control.loaded_codes[stage],
                          (std::vector<std::uint32_t>(stages[stage].begin(), stages[stage].end())));
                if (stage == 1 && eligible && kTransportFmaAvailable) {
                    EXPECT_EQ(candidate.loaded_codes[stage],
                              (std::vector<std::uint32_t>(kTransportFmaShader.begin(),
                                                          kTransportFmaShader.end())));
                    EXPECT_NE(candidate.loaded_codes[stage], control.loaded_codes[stage]);
                } else if (stage == 2 && eligible && kEndpointFmaAvailable) {
                    EXPECT_EQ(candidate.loaded_codes[stage],
                              (std::vector<std::uint32_t>(kEndpointFmaShader.begin(),
                                                          kEndpointFmaShader.end())));
                    EXPECT_NE(candidate.loaded_codes[stage], control.loaded_codes[stage]);
                } else {
                    EXPECT_EQ(candidate.loaded_codes[stage], control.loaded_codes[stage]);
                }
            }
            EXPECT_EQ(candidate.buffers.size(), 2 * RetainedCompute::kStageCount);
            ASSERT_EQ(candidate.buffers.size(), control.buffers.size());
            for (std::size_t i = 0; i < candidate.buffers.size(); ++i)
                EXPECT_EQ(candidate.buffers[i].size(), control.buffers[i].size());
            const auto required = RetainedCompute::RequiredAllocationBytes(candidate, 2);
            ASSERT_TRUE(required);
            EXPECT_EQ(candidate.BufferAllocationBytes(), *required);
        }
    }
    for (const bool has_fp64 : {false, true}) {
        PreparationProbeDevice refused;
        refused.info = {.supports_fp64 = has_fp64,
                        .preserves_fp32_denormals = true,
                        .rounds_fp32_to_nearest = true,
                        .rounds_fp64_to_nearest = !has_fp64,
                        .preserves_fp32_signed_zero_inf_nan = true,
                        .fma_fp32_enabled = true};
        EXPECT_FALSE(RetainedCompute::Create(refused, 2, true));
        EXPECT_EQ(refused.loads, 0U);
        EXPECT_TRUE(refused.buffers.empty());
    }
}
#endif

TEST(RetainedComputeAdmission, SoftwareRendererPreparationPreservesPhysicalAccounting) {
#ifdef SIRIUS_HAS_RETAINED_COMPUTE
    using Stats = RetainedCompute::PreparationStats;
    for (const auto kind :
         {DeviceKind::kIntegratedGpu, DeviceKind::kDiscreteGpu, DeviceKind::kOther}) {
        PreparationProbeDevice probe;
        probe.info.kind = kind;
        probe.info.preserves_fp32_denormals = true;
        probe.info.rounds_fp32_to_nearest = true;
        auto created = RetainedCompute::Create(probe, 2);
        ASSERT_TRUE(created);
        Stats observation;
        observation.stages[0].attempts = 42;
        const auto refused = (*created)->PrepareSoftwareRendererStages(observation);
        ASSERT_FALSE(refused);
        EXPECT_NE(refused.error().detail().find("software"), std::string::npos);
        EXPECT_TRUE(probe.writes.empty());
        EXPECT_TRUE(probe.kernels.empty());
        for (const auto& stage : observation.stages) EXPECT_EQ(stage.attempts, 0U);
    }
    PreparationProbeDevice probe;
    auto created = RetainedCompute::Create(probe, 2);
    ASSERT_TRUE(created);
    auto& compute = **created;
    ASSERT_EQ(probe.loads, RetainedCompute::kStageCount);
    ASSERT_EQ(probe.buffers.size(), 2 * RetainedCompute::kStageCount);
    ASSERT_TRUE(probe.writes.empty());
    ASSERT_TRUE(probe.kernels.empty()) << "Create must remain free of preparation work";
    const auto allocated = probe.BufferAllocationBytes();
    // Preserve an existing physical observation, including an adjacent hard
    // boundary overshoot, rather than merely proving zero stays zero.
    probe.observation.submit_wait_ms = std::nextafter(1000.0, INFINITY);
    const std::array<RetainedStepInput, 1> input{};
    ASSERT_TRUE(compute.Step(input));
    const auto physical = compute.Statistics();
    auto expected_buffers = probe.buffers;
    const auto before_writes = probe.writes.size();
    probe.observation.submit_wait_ms = 5000;  // Controlled initialization observation only.
    Stats observation;
    ASSERT_TRUE(compute.PrepareSoftwareRendererStages(observation));
    ASSERT_NO_FATAL_FAILURE(ExpectPhysicalStatsEqual(compute.Statistics(), physical));
    EXPECT_EQ(probe.BufferAllocationBytes(), allocated);
    EXPECT_EQ(probe.kernels, (std::vector<std::uint32_t>{1, 5, 4, 1, 2, 3, 6}));
    ASSERT_EQ(probe.writes.size(), before_writes + 12);
    for (std::size_t i = 0; i < 6; ++i) {
        const auto& zero = probe.writes[before_writes + 2 * i];
        const auto& restored = probe.writes[before_writes + 2 * i + 1];
        EXPECT_EQ(zero.bytes, 4U);
        EXPECT_EQ(zero.header, 0U);
        EXPECT_EQ(restored.bytes, 4U);
        EXPECT_EQ(restored.header, compute.Capacity());
        EXPECT_EQ(zero.buffer.value, restored.buffer.value);
        const std::uint32_t capacity = static_cast<std::uint32_t>(compute.Capacity());
        std::memcpy(expected_buffers[zero.buffer.value].data(), &capacity, sizeof(capacity));
    }
    EXPECT_EQ(probe.buffers, expected_buffers) << "Only the six capacity headers may change";
    for (std::size_t i = 0; i < observation.stages.size(); ++i) {
        const auto& stage = observation.stages[i];
        const auto count = i == 0 ? 0U : 1U;
        EXPECT_EQ(stage.attempts, count);
        EXPECT_EQ(stage.dispatch_attempts, count);
        EXPECT_EQ(stage.completed_dispatches, count);
        EXPECT_EQ(stage.completed, count);
        EXPECT_EQ(stage.header_restored, i != 0);
        EXPECT_EQ(stage.write_buffer_calls, 2 * count);
        EXPECT_EQ(stage.write_buffer_bytes, 8 * count);
        EXPECT_EQ(stage.timing.submit_wait_ms, i == 0 ? 0 : 5000);
        EXPECT_EQ(stage.timing.pipeline_setup_ms, i == 0 ? 0 : 1);
        EXPECT_EQ(stage.timing.command_setup_ms, i == 0 ? 0 : 2);
        EXPECT_EQ(stage.timing.cleanup_ms, i == 0 ? 0 : 3);
        EXPECT_EQ(stage.timing.total_ms, i == 0 ? 0 : 5006);
        EXPECT_GE(stage.write_buffer_ms, 0);
    }
    const auto feedback = compute.TakeSubmissionFeedback();
    EXPECT_EQ(feedback.peak_ms, std::nextafter(1000.0, INFINITY));
    EXPECT_EQ(feedback.peak_rows, 1U);
    EXPECT_EQ(feedback.maximum_rows, 1U);
    EXPECT_EQ(feedback.maximum_one_row_ms, feedback.peak_ms);
    EXPECT_EQ(feedback.maximum_one_row_stage, RetainedCompute::KernelStage::kTransport);
    // A previously uploaded stage keeps its prefix behavior; an unused Camera
    // still needs its first complete program upload after preparation.
    ASSERT_TRUE(compute.Step(input));
    EXPECT_EQ(probe.writes.back().bytes, 4U + sizeof(RetainedStepInput));
    const std::array<RetainedCameraInput, 1> camera{};
    ASSERT_TRUE(compute.Camera(camera));
    EXPECT_EQ(probe.writes.back().bytes, probe.buffers[0].size());

    // A retained owner starts a fresh observation after all prior synchronous
    // work. Exercise every stage plus the shared pair so old maxima/counts and
    // feedback cannot leak into a later render's governor or report.
    const std::array<sirius::backend::RetainedEndpointInput, 2> endpoint{};
    const std::array<sirius::backend::RetainedDenseInput, 2> dense{};
    const std::array<sirius::backend::RetainedInitializeInput, 1> initialize{};
    const std::array<sirius::backend::RetainedRayCameraInput, 1> ray_camera{};
    probe.expected_groups_x = 2;
    ASSERT_TRUE(compute.Endpoint(endpoint));
    ASSERT_TRUE(compute.Dense(dense));
    probe.expected_groups_x = 1;
    ASSERT_TRUE(compute.Initialize(initialize));
    ASSERT_TRUE(compute.RayCamera(ray_camera));
    const std::array<RetainedDopriPhaseInput, 1> dopri{};
    ASSERT_TRUE(compute.DopriPhase(dopri));
    probe.independent_pair = true;
    probe.expected_groups_x = 2;
    ASSERT_TRUE(compute.EndpointAndDense(endpoint, dense));
    for (const auto& stage : compute.Statistics()) EXPECT_GT(stage.submissions, 0U);
    EXPECT_GT(compute.EndpointDenseStatistics().submissions, 0U);
    const auto unchanged_buffers = probe.buffers;
    compute.ResetStatistics();
    ASSERT_NO_FATAL_FAILURE(ExpectPhysicalStatsEqual(compute.Statistics(), {}));
    std::array<RetainedCompute::StageStats, RetainedCompute::kStageCount> reset_shared{};
    reset_shared[0] = compute.EndpointDenseStatistics();
    ASSERT_NO_FATAL_FAILURE(ExpectPhysicalStatsEqual(reset_shared, {}));
    const auto reset_feedback = compute.TakeSubmissionFeedback();
    EXPECT_EQ(reset_feedback.peak_ms, 0);
    EXPECT_EQ(reset_feedback.peak_rows, 0U);
    EXPECT_EQ(reset_feedback.maximum_rows, 0U);
    EXPECT_EQ(reset_feedback.maximum_one_row_ms, 0);
    EXPECT_FALSE(reset_feedback.maximum_one_row_stage.has_value());
    EXPECT_EQ(probe.buffers, unchanged_buffers);
    EXPECT_EQ(probe.loads, RetainedCompute::kStageCount);
    EXPECT_EQ(probe.buffers.size(), 2 * RetainedCompute::kStageCount);
    probe.observation.submit_wait_ms = 8;
    probe.expected_groups_x = 1;
    ASSERT_TRUE(compute.Step(input));
    EXPECT_EQ(probe.writes.back().bytes, 4U + sizeof(RetainedStepInput));
    EXPECT_EQ(compute.Statistics()[1].submissions, 1U);
    EXPECT_EQ(compute.Statistics()[1].maximum_submit_wait_ms, 8);
    EXPECT_EQ(compute.TakeSubmissionFeedback().peak_ms, 8);

    for (unsigned mode = 0; mode < 8; ++mode) {
        SCOPED_TRACE(mode);
        PreparationProbeDevice failing;
        auto prepared = RetainedCompute::Create(failing, 2);
        ASSERT_TRUE(prepared);
        bool cancelled = mode == 4;
        if (mode == 0) failing.failed_writes = {1};
        if (mode == 1 || mode == 3) failing.fail_dispatch = true;
        if (mode == 2 || mode == 3) failing.failed_writes.push_back(2);
        if (mode == 5) failing.after_write = [&] { cancelled = true; };
        if (mode == 6) failing.after_dispatch = [&] { cancelled = true; };
        if (mode == 7)
            failing.after_write = [&] {
                if (failing.writes.size() == 2) cancelled = true;
            };
        const auto owner = std::this_thread::get_id();
        Stats partial;
        const auto result = (*prepared)->PrepareSoftwareRendererStages(partial, [&] {
            EXPECT_EQ(std::this_thread::get_id(), owner);
            return cancelled;
        });
        ASSERT_FALSE(result);
        EXPECT_EQ(result.error().operation(), mode == 0 || mode == 2 ? "preparation probe write"
                                              : mode == 1 || mode == 3
                                                  ? "preparation probe dispatch"
                                                  : "prepare retained renderer");
        const auto& stage = partial.stages[5];
        EXPECT_EQ(stage.attempts, mode == 4 ? 0U : 1U);
        EXPECT_EQ(stage.dispatch_attempts, mode == 0 || mode == 4 || mode == 5 ? 0U : 1U);
        EXPECT_EQ(stage.completed_dispatches, mode == 2 || mode == 6 || mode == 7 ? 1U : 0U);
        EXPECT_EQ(stage.completed, mode == 6 || mode == 7 ? 1U : 0U);
        EXPECT_EQ(stage.header_restored, mode != 2 && mode != 3 && mode != 4);
        EXPECT_EQ(stage.write_buffer_calls, mode == 4 ? 0U : 2U);
        EXPECT_EQ(stage.write_buffer_bytes, mode == 4                             ? 0U
                                            : mode == 0 || mode == 2 || mode == 3 ? 4U
                                                                                  : 8U);
        if (stage.dispatch_attempts) {
            EXPECT_EQ(stage.timing.total_ms, 5006);
        }
        EXPECT_TRUE(std::isfinite(partial.wall_ms));
        EXPECT_GE(partial.wall_ms, stage.write_buffer_ms);
        for (std::size_t i = 0; i < RetainedCompute::kStageCount; ++i) {
            if (i != 5) {
                EXPECT_EQ(partial.stages[i].attempts, 0U);
            }
        }
        EXPECT_EQ((*prepared)->TakeSubmissionFeedback().peak_ms, 0);
        ASSERT_NO_FATAL_FAILURE(ExpectPhysicalStatsEqual((*prepared)->Statistics(), {}));
        if (stage.header_restored) {
            std::uint32_t header = 0;
            std::memcpy(&header, failing.buffers[10].data(), sizeof(header));
            EXPECT_EQ(header, 2U);
        }
    }
    for (double DispatchTiming::*field :
         {&DispatchTiming::pipeline_setup_ms, &DispatchTiming::command_setup_ms,
          &DispatchTiming::submit_wait_ms, &DispatchTiming::cleanup_ms,
          &DispatchTiming::total_ms}) {
        for (const double invalid : {-1.0, std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::infinity()}) {
            PreparationProbeDevice invalid_timing;
            invalid_timing.observation.*field = invalid;
            auto prepared = RetainedCompute::Create(invalid_timing, 2);
            ASSERT_TRUE(prepared);
            Stats partial;
            const auto result = (*prepared)->PrepareSoftwareRendererStages(partial);
            ASSERT_FALSE(result);
            EXPECT_EQ(result.error().detail(), "invalid submission timing");
            EXPECT_EQ(partial.stages[5].dispatch_attempts, 1U);
            EXPECT_EQ(partial.stages[5].completed_dispatches, 0U);
            EXPECT_EQ(partial.stages[5].completed, 0U);
            EXPECT_TRUE(partial.stages[5].header_restored);
            EXPECT_EQ(partial.stages[5].write_buffer_bytes, 8U);
            for (std::size_t i = 0; i < RetainedCompute::kStageCount; ++i) {
                if (i != 5) {
                    EXPECT_EQ(partial.stages[i].attempts, 0U);
                }
            }
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
        auto required = RetainedCompute::RequiredAllocationBytes(*device, 24);
        ASSERT_TRUE(required) << required.error().Description();
        EXPECT_EQ(device->BufferAllocationBytes(), 0U);
        auto created = RetainedCompute::Create(*device, 24);
        ASSERT_TRUE(created) << created.error().Description();
        EXPECT_EQ(device->BufferAllocationBytes(), *required);
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
    std::function<void()> before_dispatch;
    std::function<void(std::span<const BufferHandle>, const DispatchTiming*)> after_dispatch;
    std::optional<BufferHandle> capture_output;
    std::vector<std::vector<std::byte>> readbacks;
    bool independent_pairs = false;
    std::uint64_t pair_calls = 0;
    std::function<void()> after_pair;
    std::function<void(const IndependentPairTiming*)> after_pair_timing;
    std::optional<double> pair_submit_ms;
    std::function<void(BufferHandle, std::span<std::byte>)> after_read;

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
        if (status && after_read) after_read(buffer, data);
        return status;
    }
    bool SupportsIndependentPair() const noexcept override {
        return independent_pairs && device_.SupportsIndependentPair();
    }
    sirius::base::Expected<void> DispatchIndependentPair(
        const std::array<ComputeDispatch, 2>& commands, IndependentPairTiming* timing) override {
        ++pair_calls;
        if (fail_next == Failure::Dispatch) return Inject("pair dispatch");
        auto status = device_.DispatchIndependentPair(commands, timing);
        if (status && timing && pair_submit_ms) {
            timing->combined.submit_wait_ms = *pair_submit_ms;
            timing->combined.total_ms = timing->combined.pipeline_setup_ms +
                                        timing->combined.command_setup_ms + *pair_submit_ms +
                                        timing->combined.cleanup_ms;
        }
        if (status && after_pair) after_pair();
        if (status && after_pair_timing) after_pair_timing(timing);
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
        if (before_dispatch) before_dispatch();
        if (fail_next == Failure::Dispatch) return Inject("dispatch");
        ++forwarded_dispatches;
        auto status = device_.Dispatch(kernel, buffers, x, y, z, timing);
        if (status && timing && submission_ms) {
            timing->submit_wait_ms = submission_ms(buffers[0], x);
            timing->total_ms = timing->pipeline_setup_ms + timing->command_setup_ms +
                               timing->submit_wait_ms + timing->cleanup_ms;
        }
        if (status && after_dispatch) after_dispatch(buffers, timing);
        return status;
    }
    sirius::base::Expected<void> SetBufferAllocationLimit(std::uint64_t bytes) override {
        return device_.SetBufferAllocationLimit(bytes);
    }
    std::uint64_t BufferAllocationBytes() const noexcept override {
        return device_.BufferAllocationBytes();
    }
    sirius::base::Expected<std::uint64_t> RequiredBufferAllocationBytes(
        std::uint64_t bytes, BufferUsage usage) override {
        return device_.RequiredBufferAllocationBytes(bytes, usage);
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

void CheckStickyErrorDrainsQueuedRequests(ComputeDevice& device) {
    for (const bool queued_camera : {true, false}) {
        SCOPED_TRACE(queued_camera ? "queued camera" : "queued interval");
        TransferProbeDevice probe(device);
        auto created = RetainedCompute::Create(probe, 1, false);
        ASSERT_TRUE(created) << created.error().Description();
        auto& compute = **created;
        RetainedTraceExecutor executor(compute, {}, 1000);
        std::promise<void> first_dispatch, release_dispatch;
        auto entered = first_dispatch.get_future();
        auto release = release_dispatch.get_future();
        probe.before_dispatch = [&] {
            if (probe.dispatch_calls == 1) {
                first_dispatch.set_value();
                release.wait();
            } else if (!queued_camera) {
                // Bound the old implementation's unnecessary interval work:
                // expose the extra call before forwarding any further stages.
                probe.fail_next = TransferProbeDevice::Failure::Dispatch;
            }
        };
        const double overshoot = std::nextafter(1000.0, std::numeric_limits<double>::infinity());
        probe.submission_ms = [overshoot](BufferHandle, std::uint32_t) { return overshoot; };
        sirius::core::CameraRay camera;
        camera.origin(1) = 5;
        camera.origin(2) = std::numbers::pi / 2;
        camera.direction(1) = 1;
        auto first = std::async(std::launch::async, [&] {
            sirius::core::KerrSchildFamily metric(sirius::core::KerrSchildParams::Minkowski());
            executor.BeginTrace();
            const auto result = executor.Launch(metric, 0, camera);
            executor.EndTrace();
            return result.has_value();
        });
        entered.wait();
        sirius::core::Lightray ray{};
        ray.position(1) = 5;
        ray.velocity(0) = -1;
        ray.velocity(1) = 1;
        ray.step_size = 1;
        const auto original = ray;
        sirius::core::Rk45CoupledState columns;
        columns.length_scale = columns.frequency_scale = 1;
        columns.tolerance = 1e-9;
        columns.variations[0].derivative(2) = .001;
        columns.variations[1].derivative(3) = .001;
        columns.variations[2].displacement(2) = 1;
        columns.variations[3].displacement(3) = 1;
        sirius::core::IntegratorConfig config;
        config.min_step = .01f;
        config.max_step = 2;
        sirius::core::Rk45CoupledComparison comparison;
        auto second = std::async(std::launch::async, [&] {
            sirius::core::KerrSchildFamily metric(sirius::core::KerrSchildParams::Minkowski());
            executor.BeginTrace();
            const bool result = queued_camera
                                    ? executor.Launch(metric, 0, camera).has_value()
                                    : executor.Step(ray, metric, config, columns, comparison);
            executor.EndTrace();
            return result;
        });
        // The first dispatch is held, so this observes the second request's
        // actual queue membership before failure; elapsed time is not the proof.
        const auto queued = RetainedTraceExecutorTestPeer::WaitForQueued(executor);
        release_dispatch.set_value();
        EXPECT_FALSE(first.get());
        EXPECT_FALSE(second.get());
        EXPECT_EQ(queued, (std::array<std::size_t, 3>{1, 1, 2}));
        EXPECT_EQ(RetainedTraceExecutorTestPeer::QueueState(executor),
                  (std::array<std::size_t, 3>{0, 0, 0}));
        ASSERT_TRUE(executor.Error());
        EXPECT_EQ(executor.Error()->detail(),
                  "single-row submission exceeded the safety duration: "
                  "stage=ray_camera, active_rows=1, "
                  "submit_wait_ms=1000.0000000000001, limit_ms=1000");
        EXPECT_EQ(probe.dispatch_calls, 1U);
        EXPECT_EQ(probe.forwarded_dispatches, 1U);
        EXPECT_EQ(probe.writes.size(), 1U);
        EXPECT_EQ(probe.reads.size(), 1U);
        ASSERT_EQ(probe.allocations.size(), 2 * RetainedCompute::kStageCount);
        for (std::size_t stage = 0; stage < compute.Statistics().size(); ++stage)
            EXPECT_EQ(compute.Statistics()[stage].submissions, stage == 5 ? 1U : 0U);
        const auto stats = executor.Statistics();
        EXPECT_EQ(stats.batches, 2U);
        EXPECT_EQ(stats.batch_row_counts[1], 2U);
        EXPECT_EQ(stats.camera_rows, queued_camera ? 2U : 1U);
        EXPECT_EQ(stats.interval_rows, queued_camera ? 0U : 1U);
        EXPECT_EQ(stats.accepted_intervals, 0U);
        EXPECT_EQ(stats.rejected_intervals, queued_camera ? 0U : 1U);
        EXPECT_EQ(stats.paired_projection_retries, 0U);
        EXPECT_EQ(stats.safety_fallbacks, 0U);
        EXPECT_EQ(stats.acceleration_calls, 0U);
        if (!queued_camera) {
            EXPECT_EQ(ray.terminated, 3);
            EXPECT_EQ(columns.failure, sirius::core::CoupledStepFailure::InvalidState);
            EXPECT_EQ(columns.central_stages, 0U);
            EXPECT_EQ(columns.variation_stages, 0U);
            for (int axis = 0; axis < 4; ++axis) {
                EXPECT_EQ(ray.position(axis), original.position(axis));
                EXPECT_EQ(ray.velocity(axis), original.velocity(axis));
            }
        }
        ::testing::Test::RecordProperty(queued_camera ? "sticky_error_queued_camera_dispatches"
                                                      : "sticky_error_queued_interval_dispatches",
                                        std::to_string(probe.dispatch_calls));
    }
}

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
            std::bit_cast<std::uint64_t>(expected.error_ratio) ||
        std::bit_cast<std::uint64_t>(actual.embedded_projected_error_ratio) !=
            std::bit_cast<std::uint64_t>(expected.embedded_projected_error_ratio))
        return ::testing::AssertionFailure() << "interval admission metadata differs";
    for (std::size_t i = 0; i < actual.error_checks.size(); ++i) {
        const auto& a = actual.error_checks[i];
        const auto& b = expected.error_checks[i];
        if (std::bit_cast<std::uint64_t>(a.ratio) != std::bit_cast<std::uint64_t>(b.ratio) ||
            a.field != b.field || a.observed != b.observed || a.evaluated != b.evaluated)
            return ::testing::AssertionFailure() << "passive error observation differs: " << i;
    }
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
    const auto required = RetainedCompute::RequiredAllocationBytes(probe, compute->Capacity());
    ASSERT_TRUE(required) << required.error().Description();
    EXPECT_EQ(device->BufferAllocationBytes(), initial_allocation);
    auto observed = RetainedCompute::Create(probe, compute->Capacity());
    ASSERT_TRUE(observed) << observed.error().Description();
    auto& camera_compute = **observed;
    const auto allocation = device->BufferAllocationBytes();
    EXPECT_GT(allocation, initial_allocation);
    EXPECT_EQ(allocation - initial_allocation, *required);
    EXPECT_LE(allocation, 8ull * 1024 * 1024);
    ASSERT_EQ(probe.allocations.size(), 2 * RetainedCompute::kStageCount);
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

TEST_F(RetainedComputeTest, SoftwareRendererPreparationPreservesBuffersAndPhysicalFeedback) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    if (device->Info().kind != DeviceKind::kSoftware) {
        const auto before = compute->Statistics();
        const auto allocated = device->BufferAllocationBytes();
        RetainedCompute::PreparationStats observation;
        const auto result = compute->PrepareSoftwareRendererStages(observation);
        ASSERT_FALSE(result);
        EXPECT_NE(result.error().detail().find("software"), std::string::npos);
        for (const auto& stage : observation.stages) {
            EXPECT_EQ(stage.attempts, 0U);
            EXPECT_EQ(stage.dispatch_attempts, 0U);
            EXPECT_EQ(stage.write_buffer_calls, 0U);
        }
        ASSERT_NO_FATAL_FAILURE(ExpectPhysicalStatsEqual(compute->Statistics(), before));
        EXPECT_EQ(device->BufferAllocationBytes(), allocated);
        EXPECT_EQ(compute->TakeSubmissionFeedback().peak_ms, 0);
        RecordProperty("software_preparation_disposition", "nonsoftware_refused_without_work");
        return;
    }
    TransferProbeDevice probe(*device);
    auto created = RetainedCompute::Create(probe, compute->Capacity());
    ASSERT_TRUE(created) << created.error().Description();
    auto& prepared = **created;
    ASSERT_EQ(probe.allocations.size(), 2 * RetainedCompute::kStageCount);
    EXPECT_TRUE(probe.writes.empty());
    EXPECT_EQ(probe.dispatch_calls, 0U);
    const auto allocated = device->BufferAllocationBytes();
    std::vector<std::vector<std::byte>> expected;
    for (std::size_t i = 0; i < probe.allocations.size(); ++i) {
        expected.emplace_back(probe.allocations[i].bytes, std::byte(i % 2 ? 0x5a : 0xa5));
        ASSERT_TRUE(device->WriteBuffer(probe.allocations[i].handle, expected.back()));
    }
    RetainedCompute::PreparationStats observation;
    const auto result = prepared.PrepareSoftwareRendererStages(observation);
    // Publish all observations even if preparation failed; these are actual
    // initialization host intervals, not physical work or GPU timestamps.
    RecordProperty("software_preparation_scope",
                   "six explicit zero-row initialization submissions; preparation is outside "
                   "physical feedback; subsequent RayCamera is not cold-render qualification");
    RecordProperty("software_preparation_wall_ms", std::to_string(observation.wall_ms));
    for (std::size_t i = 0; i < observation.stages.size(); ++i) {
        const auto& stage = observation.stages[i];
        const auto name = RetainedCompute::StageName(static_cast<RetainedCompute::KernelStage>(i));
        const std::string prefix = std::string("software_preparation_") + name;
        std::string message = "[Retained preparation] " + std::string(name);
        const auto record = [&](const char* key, const auto value) {
            const auto text = std::format("{}", value);
            RecordProperty(prefix + "_" + key, text);
            message += std::format(" {}={}", key, text);
        };
        record("attempts", stage.attempts);
        record("dispatch_attempts", stage.dispatch_attempts);
        record("completed_dispatches", stage.completed_dispatches);
        record("completed", stage.completed);
        record("header_restored", stage.header_restored);
        record("pipeline_setup_ms", stage.timing.pipeline_setup_ms);
        record("command_setup_ms", stage.timing.command_setup_ms);
        record("submit_wait_ms", stage.timing.submit_wait_ms);
        record("cleanup_ms", stage.timing.cleanup_ms);
        record("dispatch_total_ms", stage.timing.total_ms);
        record("pipeline_created", stage.timing.pipeline_created);
        record("write_buffer_calls", stage.write_buffer_calls);
        record("write_buffer_ms", stage.write_buffer_ms);
        record("write_buffer_bytes", stage.write_buffer_bytes);
        std::cerr << message << std::endl;
    }
    ASSERT_TRUE(result) << result.error().Description();
    const std::array<std::size_t, 6> order{5, 4, 1, 2, 3, 6};
    ASSERT_EQ(probe.submissions.size(), order.size());
    ASSERT_EQ(probe.writes.size(), 2 * order.size());
    EXPECT_TRUE(probe.reads.empty());
    for (std::size_t j = 0; j < order.size(); ++j) {
        const auto i = order[j];
        const auto input = probe.allocations[2 * i].handle;
        const auto output = probe.allocations[2 * i + 1].handle;
        EXPECT_EQ(probe.submissions[j].buffers[0].value, input.value);
        EXPECT_EQ(probe.submissions[j].buffers[1].value, output.value);
        EXPECT_EQ(probe.submissions[j].binding_count, 2U);
        EXPECT_EQ(probe.submissions[j].x, 1U);
        EXPECT_EQ(probe.submissions[j].y, 1U);
        EXPECT_EQ(probe.submissions[j].z, 1U);
        EXPECT_EQ(probe.writes[2 * j].bytes, 4U);
        EXPECT_EQ(probe.writes[2 * j].capacity, 0U);
        EXPECT_EQ(probe.writes[2 * j + 1].bytes, 4U);
        EXPECT_EQ(probe.writes[2 * j + 1].capacity, prepared.Capacity());
        const auto capacity = static_cast<std::uint32_t>(prepared.Capacity());
        std::memcpy(expected[2 * i].data(), &capacity, sizeof(capacity));
        const auto& stage = observation.stages[i];
        EXPECT_EQ(stage.attempts, 1U);
        EXPECT_EQ(stage.dispatch_attempts, 1U);
        EXPECT_EQ(stage.completed_dispatches, 1U);
        EXPECT_EQ(stage.completed, 1U);
        EXPECT_TRUE(stage.header_restored);
        EXPECT_EQ(stage.write_buffer_calls, 2U);
        EXPECT_EQ(stage.write_buffer_bytes, 8U);
        EXPECT_GE(stage.timing.submit_wait_ms, 0);
    }
    EXPECT_EQ(observation.stages[0].attempts, 0U);
    for (std::size_t i = 0; i < expected.size(); ++i) {
        SCOPED_TRACE(i);
        std::vector<std::byte> actual(expected[i].size());
        ASSERT_TRUE(device->ReadBuffer(probe.allocations[i].handle, actual));
        EXPECT_EQ(actual, expected[i]) << "zero rows must preserve output and nonheader inputs";
    }
    EXPECT_EQ(device->BufferAllocationBytes(), allocated);
    ASSERT_NO_FATAL_FAILURE(ExpectPhysicalStatsEqual(prepared.Statistics(), {}));
    const auto empty_feedback = prepared.TakeSubmissionFeedback();
    EXPECT_EQ(empty_feedback.peak_ms, 0);
    EXPECT_EQ(empty_feedback.maximum_rows, 0U);
    EXPECT_FALSE(empty_feedback.maximum_one_row_stage);
    const auto& fixture = sirius::test::retained_ray_camera::cases.front();
    const std::array input{std::bit_cast<RetainedRayCameraInput>(fixture.input)};
    DispatchTiming timing;
    const auto physical = prepared.RayCamera(input, &timing);
    RecordProperty("prepared_ray_camera_physical_submit_wait_ms",
                   std::to_string(timing.submit_wait_ms));
    RecordProperty("prepared_ray_camera_physical_pipeline_setup_ms",
                   std::to_string(timing.pipeline_setup_ms));
    RecordProperty("prepared_ray_camera_physical_dispatch_total_ms",
                   std::to_string(timing.total_ms));
    std::cerr << std::format(
                     "[Retained preparation] subsequent physical RayCamera: submit/wait "
                     "{}ms, pipeline {}ms, total {}ms",
                     timing.submit_wait_ms, timing.pipeline_setup_ms, timing.total_ms)
              << std::endl;
    ASSERT_TRUE(physical) << physical.error().Description();
    ASSERT_EQ(physical->size(), 1U);
    ASSERT_TRUE(CameraAgrees(physical->front(), fixture));
    EXPECT_EQ(probe.writes.back().bytes, probe.allocations[10].bytes)
        << "preparation must not claim the immutable program was uploaded";
    const auto feedback = prepared.TakeSubmissionFeedback();
    EXPECT_EQ(feedback.peak_ms, timing.submit_wait_ms);
    EXPECT_EQ(feedback.maximum_rows, 1U);
    EXPECT_EQ(feedback.maximum_one_row_ms, timing.submit_wait_ms);
    EXPECT_EQ(feedback.maximum_one_row_stage, RetainedCompute::KernelStage::kRayCamera);
    EXPECT_EQ(prepared.Statistics()[5].submissions, 1U);
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

    // Observe first use and repetition of exactly the same original row on the
    // fresh fixture compute. These are host intervals, not GPU timestamps or a
    // claim about driver cache state or the operational renderer's hard guard.
    RecordProperty("ray_camera_reuse_device", device->Info().name);
    RecordProperty("ray_camera_reuse_row", sirius::test::retained_ray_camera::cases.front().name);
    RecordProperty("ray_camera_reuse_timing_scope",
                   "same original row twice on fresh fixture compute; host timings exclude "
                   "SetUp; no speed threshold, cache attribution or cold render qualification");
    const auto observe = [&](const std::string& prefix) {
        const auto before =
            compute
                ->Statistics()[static_cast<std::size_t>(RetainedCompute::KernelStage::kRayCamera)];
        DispatchTiming timing;
        std::cerr << "[RayCamera reuse] " << prefix << " started" << std::endl;
        const auto started = std::chrono::steady_clock::now();
        auto result = compute->RayCamera(std::span(inputs.data(), 1), &timing);
        const auto finished = std::chrono::steady_clock::now();
        const auto after =
            compute
                ->Statistics()[static_cast<std::size_t>(RetainedCompute::KernelStage::kRayCamera)];
        auto message = std::format("[RayCamera reuse] {} completed", prefix);
        const auto record = [&](const char* key, const auto& value) {
            const auto text = std::format("{}", value);
            RecordProperty(prefix + "_" + key, text);
            message += std::format(" {}={}", key, text);
        };
        record("success", result.has_value());
        record("wall_ms", std::chrono::duration<double, std::milli>(finished - started).count());
        record("submit_wait_ms", timing.submit_wait_ms);
        record("pipeline_setup_ms", timing.pipeline_setup_ms);
        record("command_setup_ms", timing.command_setup_ms);
        record("cleanup_ms", timing.cleanup_ms);
        record("dispatch_total_ms", timing.total_ms);
        record("pipeline_created", timing.pipeline_created);
        record("stage_submissions", after.submissions - before.submissions);
        record("stage_submit_wait_ms", after.submit_wait_ms - before.submit_wait_ms);
        record("stage_pipeline_setup_ms", after.pipeline_setup_ms - before.pipeline_setup_ms);
        record("stage_command_setup_ms", after.command_setup_ms - before.command_setup_ms);
        record("stage_cleanup_ms", after.cleanup_ms - before.cleanup_ms);
        record("stage_dispatch_total_ms", after.dispatch_total_ms - before.dispatch_total_ms);
        record("stage_write_buffer_ms", after.write_buffer_ms - before.write_buffer_ms);
        record("stage_read_buffer_ms", after.read_buffer_ms - before.read_buffer_ms);
        record("stage_write_buffer_bytes", after.write_buffer_bytes - before.write_buffer_bytes);
        record("stage_read_buffer_bytes", after.read_buffer_bytes - before.read_buffer_bytes);
        record("stage_pipeline_creations", after.pipeline_creations - before.pipeline_creations);
        record("stage_target_overshoots", after.target_overshoots - before.target_overshoots);
        if (!result) message += std::format(" error={}", result.error().Description());
        std::cerr << message << std::endl;
        return result;
    };
    const auto first = observe("ray_camera_first");
    ASSERT_TRUE(first) << first.error().Description();
    ASSERT_EQ(first->size(), 1U);
    ASSERT_TRUE(CameraAgrees(first->front(), sirius::test::retained_ray_camera::cases.front()));
    const auto repeated = observe("ray_camera_repeat");
    ASSERT_TRUE(repeated) << repeated.error().Description();
    ASSERT_EQ(repeated->size(), 1U);
    ASSERT_TRUE(CameraAgrees(repeated->front(), sirius::test::retained_ray_camera::cases.front()));
    EXPECT_EQ(repeated->front().valid, first->front().valid);
    for (std::size_t i = 0; i < first->front().values.size(); ++i) {
        SCOPED_TRACE(i);
        EXPECT_EQ((std::bit_cast<std::array<std::uint32_t, 5>>(repeated->front().values[i])),
                  (std::bit_cast<std::array<std::uint32_t, 5>>(first->front().values[i])));
    }

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
    // Both represented temporal roots near -1 +/- 2^-45 fail the denominator
    // guard: at M=1,r=1,chart=+1 their denominator magnitudes (~2^-45) are
    // below its ~2^-43 cutoff. The accepted root must therefore be spatial.
    {
        SCOPED_TRACE("guarded temporal roots use spatial fallback");
        RetainedEndpointInput fallback_input;
        fallback_input.values.fill(RetainedValue::FromDouble(0));
        fallback_input.values[0] = RetainedValue::FromDouble(1);
        fallback_input.values[5] = RetainedValue::FromDouble(1);
        fallback_input.values[9] = RetainedValue::FromDouble(-.5);
        fallback_input.values[10] = std::bit_cast<RetainedValue>(
            std::array<std::uint32_t, 5>{0x3f000000U, 0x92800000U, 0U, 0U, 1U});
        fallback_input.values[44] = RetainedValue::FromDouble(1);
        const std::array fallback_inputs{fallback_input};
        const auto fallback = compute->Endpoint(fallback_inputs);
        ASSERT_TRUE(fallback) << fallback.error().Description();
        ASSERT_EQ(fallback->size(), 1U);
        const auto& endpoint = fallback->front();
        ASSERT_TRUE(endpoint.valid);
        ASSERT_TRUE(endpoint.component == 1 || endpoint.component == 2);
        // Independent dyadic centers retain delta in low limbs. Component 1
        // approximation error is below 13*delta^2 < 8.5e-54; component 2 is exact.
        constexpr double delta = 0x1p-90;
        using sirius::core::Twofold;
        std::array<Twofold, 40> phase_reference{}, physical_reference{};
        phase_reference[1] = physical_reference[1] = Twofold(1);
        physical_reference[4] = Twofold(-1);
        if (endpoint.component == 1) {
            phase_reference[4] = Twofold(-2 * delta);
            phase_reference[5] = Twofold(-.5, -3 * delta);
            phase_reference[6] = Twofold(.5, -delta);
            physical_reference[5] = physical_reference[6] = Twofold(.5, -delta);
        } else {
            phase_reference[5] = Twofold(-.5);
            phase_reference[6] = Twofold(.5);
            physical_reference[5] = physical_reference[6] = Twofold(.5);
        }
        for (std::size_t i = 0; i < 80; ++i) {
            SCOPED_TRACE(i);
            const auto& value = i < 40 ? endpoint.phase[i] : endpoint.physical[i - 40];
            const auto& oracle = i < 40 ? phase_reference[i] : physical_reference[i - 40];
            const auto center = Twofold(value.high) + Twofold(value.low) + Twofold(value.tail);
            const double difference = std::abs((center - oracle).Rounded());
            EXPECT_LE(difference, value.radius + 1e-29 * (1 + std::abs(oracle.hi)));
            EXPECT_LE(difference, 1e-11 * (1 + std::abs(oracle.hi)));
        }
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

    // Exact endpoints own their values even when the interior secant exceeds
    // the retained divide domain. This synthetic packet tests representation,
    // rather than a unit-speed geodesic with this displacement and duration.
    std::array<RetainedDenseInput, 3> endpoint_masks;
    const std::array<double, 3> fractions{0, 1, .5};
    for (std::size_t row = 0; row < endpoint_masks.size(); ++row) {
        auto& input = endpoint_masks[row];
        input.values.fill(RetainedValue::FromDouble(0));
        input.values[8] = input.values[9] = RetainedValue::FromDouble(1);
        input.values[48] = input.values[49] = RetainedValue::FromDouble(1);
        input.values[104] = RetainedValue::FromDouble(0x1p-100);
        input.values[44] = input.values[84] = input.values[104];
        input.values[45] = input.values[85] = RetainedValue::FromDouble(1);
        input.values[105] = RetainedValue::FromDouble(fractions[row]);
        input.values[111] = RetainedValue::FromDouble(1);
    }
    const auto masked = compute->Dense(endpoint_masks);
    ASSERT_TRUE(masked) << masked.error().Description();
    ASSERT_EQ(masked->size(), endpoint_masks.size());
    for (std::size_t row = 0; row < 2; ++row) {
        SCOPED_TRACE(row);
        ASSERT_TRUE((*masked)[row].valid);
        for (std::size_t i = 0; i < 40; ++i) {
            SCOPED_TRACE(i);
            const auto expected = endpoint_masks[row].values[(row == 0 ? 4 : 44) + i];
            EXPECT_EQ((std::bit_cast<std::array<std::uint32_t, 5>>((*masked)[row].physical[i])),
                      (std::bit_cast<std::array<std::uint32_t, 5>>(expected)));
        }
    }
    EXPECT_FALSE(masked->back().valid);

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

TEST_F(RetainedComputeTest, DeviceTimestampsPreserveOriginalIntervalResults) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    auto* vulkan = dynamic_cast<VulkanDevice*>(device.get());
    ASSERT_NE(vulkan, nullptr);
    const auto properties = vulkan->TimestampProperties();
    RecordProperty("timestamp_valid_bits", std::to_string(properties.valid_bits));
    RecordProperty("timestamp_period_ns", std::format("{:.17g}", properties.period_ns));
    RecordProperty("timestamp_scope",
                   "opt-in device marker span includes original barriers and scheduling; not "
                   "isolated shader time or calibrated host/driver overhead; existing host "
                   "governor unchanged; finite original critical0/critical1/flat intervals, "
                   "not frame, cold, full-science or release qualification");
    const auto enabled = vulkan->SetDispatchTimestampsEnabled(true);
    if (properties.valid_bits == 0) {
        ASSERT_FALSE(enabled);
        RecordProperty("timestamp_observation", "unsupported selected queue explicitly declined");
        return;
    }
    ASSERT_TRUE(enabled) << enabled.error().Description();
    ASSERT_TRUE(vulkan->SetDispatchTimestampsEnabled(false));
    std::vector<RetainedEndpointInput> launch_inputs;
    for (const auto index :
         {std::size_t{0}, std::size_t{1}, sirius::test::retained_transport::cases.size() - 1}) {
        const auto step =
            std::bit_cast<RetainedStepInput>(sirius::test::retained_transport::cases[index].input);
        RetainedEndpointInput endpoint;
        std::copy_n(step.values.begin(), 45, endpoint.values.begin());
        launch_inputs.push_back(endpoint);
    }
    const auto launches = compute->Endpoint(launch_inputs);
    ASSERT_TRUE(launches) << launches.error().Description();
    std::vector<RetainedIntervalInput> inputs(launch_inputs.size());
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        const auto index = row < 2 ? row : sirius::test::retained_transport::cases.size() - 1;
        const auto step =
            std::bit_cast<RetainedStepInput>(sirius::test::retained_transport::cases[index].input);
        auto& input = inputs[row];
        std::copy_n(step.values.begin(), 4, input.metric.begin());
        input.start = (*launches)[row];
        input.chart = step.values[44].Center();
        input.interval = step.values[45].Center();
        input.control.length_scale = input.control.frequency_scale = 1;
        input.control.tolerance = 1e-4 / (4 * 30000);
        input.control.integrator.min_step = static_cast<float>(input.interval);
        input.control.integrator.max_step = 1;
        input.control.integrator.abs_tolerance = input.control.integrator.rel_tolerance = 1e-9f;
    }
    const auto baseline = AttemptRetainedIntervals(*compute, inputs);
    ASSERT_TRUE(baseline) << baseline.error().Description();
    EXPECT_FALSE(vulkan->LastDispatchTimestamp());
    TransferProbeDevice probe(*device);
    probe.independent_pairs = true;
    auto observed = RetainedCompute::Create(probe, 24);
    ASSERT_TRUE(observed) << observed.error().Description();
    ASSERT_EQ(probe.allocations.size(), 2 * RetainedCompute::kStageCount);
    const auto resident = device->BufferAllocationBytes();
    ASSERT_TRUE(vulkan->SetDispatchTimestampsEnabled(true));
    std::size_t samples = 0;
    std::string phase;
    const auto capture = [&](const std::string& stage, const DispatchTiming* host) {
        ASSERT_LT(samples, 64U);
        ASSERT_NE(host, nullptr);
        const auto timestamp = vulkan->LastDispatchTimestamp();
        ASSERT_TRUE(timestamp);
        ASSERT_EQ(timestamp->query_result, VK_SUCCESS);
        ASSERT_NE(timestamp->availability[0], 0U);
        ASSERT_NE(timestamp->availability[1], 0U);
        ASSERT_TRUE(timestamp->device_span_ms);
        EXPECT_TRUE(std::isfinite(*timestamp->device_span_ms));
        EXPECT_GE(*timestamp->device_span_ms, 0);
        EXPECT_NEAR(timestamp->host_submit_ms + timestamp->host_wait_ms, host->submit_wait_ms,
                    1e-9 * std::max(1., host->submit_wait_ms));
        RecordProperty(
            "device_observation_" + std::to_string(samples++),
            std::format("phase={};stage={};begin={};end={};available0={};available1={};"
                        "device_ms={:.17g};host_submit_ms={:.17g};host_wait_ms={:.17g};"
                        "host_completion_ms={:.17g};pipeline_ms={:.17g};total_ms={:.17g}",
                        phase, stage, timestamp->ticks[0], timestamp->ticks[1],
                        timestamp->availability[0], timestamp->availability[1],
                        *timestamp->device_span_ms, timestamp->host_submit_ms,
                        timestamp->host_wait_ms, host->submit_wait_ms, host->pipeline_setup_ms,
                        host->total_ms));
    };
    probe.after_dispatch = [&](std::span<const BufferHandle> bound, const DispatchTiming* host) {
        constexpr std::array names{"camera", "transport",  "endpoint",
                                   "dense",  "initialize", "ray_camera"};
        ASSERT_FALSE(bound.empty());
        for (std::size_t i = 0; i < probe.allocations.size(); i += 2)
            if (bound[0].value == probe.allocations[i].handle.value) {
                capture(names[i / 2], host);
                return;
            }
        FAIL() << "unknown observed stage binding";
    };
    probe.after_pair_timing = [&](const IndependentPairTiming* host) {
        ASSERT_NE(host, nullptr);
        capture("endpoint+dense", &host->combined);
    };
    for (unsigned repeat = 0; repeat < 3; ++repeat) {
        phase = "repeat" + std::to_string(repeat);
        const auto output = AttemptRetainedIntervals(**observed, inputs);
        ASSERT_TRUE(output) << output.error().Description();
        ASSERT_EQ(output->size(), baseline->size());
        for (std::size_t row = 0; row < output->size(); ++row) {
            const auto& actual = (*output)[row];
            const auto& expected = (*baseline)[row];
            ASSERT_TRUE(actual.admissible);
            EXPECT_EQ(actual.admissible, expected.admissible);
            EXPECT_EQ(actual.attempted_stages, 21U);
            EXPECT_EQ(actual.attempted_stages, expected.attempted_stages);
            EXPECT_EQ(actual.failure, expected.failure);
            EXPECT_EQ(std::bit_cast<std::uint64_t>(actual.error_ratio),
                      std::bit_cast<std::uint64_t>(expected.error_ratio));
            EXPECT_EQ(std::bit_cast<std::uint64_t>(actual.embedded_projected_error_ratio),
                      std::bit_cast<std::uint64_t>(expected.embedded_projected_error_ratio));
            for (const auto& endpoints : {std::pair{&actual.full, &expected.full},
                                          {&actual.lower, &expected.lower},
                                          {&actual.midpoint, &expected.midpoint},
                                          {&actual.refined, &expected.refined}}) {
                EXPECT_EQ(endpoints.first->valid, endpoints.second->valid);
                EXPECT_EQ(endpoints.first->component, endpoints.second->component);
                EXPECT_EQ(std::memcmp(endpoints.first->phase.data(), endpoints.second->phase.data(),
                                      sizeof(endpoints.first->phase)),
                          0);
                EXPECT_EQ(
                    std::memcmp(endpoints.first->physical.data(), endpoints.second->physical.data(),
                                sizeof(endpoints.first->physical)),
                    0);
            }
            for (const auto& increments :
                 {std::pair{&actual.full_increment, &expected.full_increment},
                  {&actual.lower_increment, &expected.lower_increment},
                  {&actual.midpoint_increment, &expected.midpoint_increment},
                  {&actual.refined_increment, &expected.refined_increment}})
                EXPECT_EQ(std::memcmp(increments.first->data(), increments.second->data(),
                                      sizeof(*increments.first)),
                          0);
        }
        EXPECT_EQ(device->BufferAllocationBytes(), resident);
    }
    EXPECT_GT(samples, 0U);
    EXPECT_GT(probe.pair_calls, 0U);
    RecordProperty("device_observations", std::to_string(samples));
    ASSERT_TRUE(vulkan->SetDispatchTimestampsEnabled(false));
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST_F(RetainedComputeTest, CoupledIntervalsRequireEmbeddedAndIndependentDenseAgreement) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    using Clock = std::chrono::steady_clock;
    using Stats = std::array<RetainedCompute::StageStats, RetainedCompute::kStageCount>;
    const auto fixture_started = Clock::now();
    const auto fixture_before = compute->Statistics();
    const auto fixture_pair_before = compute->EndpointDenseStatistics();
    const auto milliseconds = [](Clock::time_point started, Clock::time_point finished) {
        return std::chrono::duration<double, std::milli>(finished - started).count();
    };
    auto record_timing = [&, previous_pair = fixture_pair_before](
                             const std::string& prefix, const Stats& before, const Stats& after,
                             double wall_ms, bool lifetime_maxima = false) mutable {
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
        const auto pair = compute->EndpointDenseStatistics();
        const auto pair_before = lifetime_maxima ? fixture_pair_before : previous_pair;
        const auto shared = pair.submissions - pair_before.submissions;
        submit_wait_ms += pair.submit_wait_ms - pair_before.submit_wait_ms;
        pipeline_setup_ms += pair.pipeline_setup_ms - pair_before.pipeline_setup_ms;
        command_setup_ms += pair.command_setup_ms - pair_before.command_setup_ms;
        cleanup_ms += pair.cleanup_ms - pair_before.cleanup_ms;
        dispatch_total_ms += pair.dispatch_total_ms - pair_before.dispatch_total_ms;
        pipeline_creations += pair.pipeline_creations - pair_before.pipeline_creations;
        target_overshoots += pair.target_overshoots - pair_before.target_overshoots;
        maximum_submit_wait_ms = std::max(maximum_submit_wait_ms, pair.maximum_submit_wait_ms);
        previous_pair = pair;
        RecordProperty(prefix + "_shared_endpoint_dense_submissions", std::to_string(shared));
        RecordProperty(prefix + "_shared_endpoint_dense_submit_wait_ms",
                       std::to_string(pair.submit_wait_ms - pair_before.submit_wait_ms));
        RecordProperty(prefix + "_shared_endpoint_dense_pipeline_setup_ms",
                       std::to_string(pair.pipeline_setup_ms - pair_before.pipeline_setup_ms));
        RecordProperty(prefix + "_shared_endpoint_dense_command_setup_ms",
                       std::to_string(pair.command_setup_ms - pair_before.command_setup_ms));
        RecordProperty(prefix + "_shared_endpoint_dense_cleanup_ms",
                       std::to_string(pair.cleanup_ms - pair_before.cleanup_ms));
        RecordProperty(prefix + "_shared_endpoint_dense_dispatch_total_ms",
                       std::to_string(pair.dispatch_total_ms - pair_before.dispatch_total_ms));
        RecordProperty(prefix + "_shared_endpoint_dense_pipeline_creations",
                       std::to_string(pair.pipeline_creations - pair_before.pipeline_creations));
        RecordProperty(prefix + "_shared_endpoint_dense_target_overshoots",
                       std::to_string(pair.target_overshoots - pair_before.target_overshoots));
        RecordProperty(prefix + "_queue_submissions", std::to_string(submissions - shared));
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
    bool positive_embedded_projected_error = false;
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        SCOPED_TRACE(sirius::test::retained_transport::cases[row].name);
        const auto& output = (*outputs)[row];
        EXPECT_TRUE(output.admissible)
            << "error=" << output.error_ratio
            << " failure=" << sirius::core::CoupledStepFailureName(output.failure);
        EXPECT_EQ(output.attempted_stages, 21U);
        EXPECT_LE(output.error_ratio, 1);
        EXPECT_GE(output.embedded_projected_error_ratio, 0);
        EXPECT_LE(output.embedded_projected_error_ratio, output.error_ratio);
        for (const auto& check : output.error_checks) {
            EXPECT_TRUE(check.observed);
            EXPECT_TRUE(check.evaluated);
            EXPECT_TRUE(std::isfinite(check.ratio));
            EXPECT_LE(check.ratio, output.error_ratio);
            EXPECT_GE(check.field, 8U);
            EXPECT_LE(check.field, 40U);
        }
        positive_embedded_projected_error =
            positive_embedded_projected_error || output.embedded_projected_error_ratio > 0;
        RecordProperty("embedded_projected_error_" + std::to_string(row),
                       std::to_string(output.embedded_projected_error_ratio));
        RecordProperty("combined_error_" + std::to_string(row), std::to_string(output.error_ratio));
    }
    EXPECT_TRUE(positive_embedded_projected_error);
    // Exact flat flow remains exact through all three independent candidates.
    const auto& flat = outputs->back();
    ASSERT_TRUE(flat.admissible);
    EXPECT_EQ(flat.error_ratio, 0);
    EXPECT_EQ(flat.embedded_projected_error_ratio, 0);
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
    EXPECT_TRUE((*rejected)[0].error_checks[1].observed);
    EXPECT_TRUE((*rejected)[0].error_checks[1].evaluated);
    double largest_observed = 0;
    for (const auto& check : (*rejected)[0].error_checks)
        if (check.observed) largest_observed = std::max(largest_observed, check.ratio);
    EXPECT_GT(largest_observed, 1);
    EXPECT_EQ(largest_observed, (*rejected)[0].error_ratio);
    for (const auto& check : (*rejected)[1].error_checks) EXPECT_FALSE(check.observed);
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
    ASSERT_EQ(probe.allocations.size(), 2 * RetainedCompute::kStageCount);
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
        EXPECT_EQ(probe.allocations.size(), 2 * RetainedCompute::kStageCount);
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

TEST_F(RetainedComputeTest, SharedEndpointDensePreservesPrivateIntervalsAndSerialRetry) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    TransferProbeDevice probe(*device);
    ASSERT_TRUE(device->SupportsIndependentPair());
    auto created = RetainedCompute::Create(probe, 4, false, 250);
    ASSERT_TRUE(created) << created.error().Description();
    auto& observed = **created;
    ASSERT_EQ(probe.allocations.size(), 2 * RetainedCompute::kStageCount);
    const auto endpoint_output = probe.allocations[5].handle;
    const auto dense_output = probe.allocations[7].handle;
    const auto resident = device->BufferAllocationBytes();
    std::array<RetainedEndpointInput, 2> initial{};
    std::array<RetainedIntervalInput, 2> inputs{};
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        const auto step =
            std::bit_cast<RetainedStepInput>(sirius::test::retained_transport::cases[row].input);
        std::copy_n(step.values.begin(), 45, initial[row].values.begin());
        auto& input = inputs[row];
        std::copy_n(step.values.begin(), 4, input.metric.begin());
        input.chart = step.values[44].Center();
        input.interval = step.values[45].Center();
        input.control.length_scale = input.control.frequency_scale = 1;
        input.control.tolerance = 1e-4 / (4 * 30000);
        input.control.integrator.min_step = static_cast<float>(input.interval);
        input.control.integrator.max_step = 1;
        input.control.integrator.abs_tolerance = input.control.integrator.rel_tolerance = 1e-9f;
    }
    const auto launched = observed.Endpoint(initial);
    ASSERT_TRUE(launched) << launched.error().Description();
    for (std::size_t row = 0; row < inputs.size(); ++row) {
        ASSERT_TRUE((*launched)[row].valid);
        inputs[row].start = (*launched)[row];
    }
    const auto serial = AttemptRetainedIntervals(observed, inputs, 4);
    ASSERT_TRUE(serial) << serial.error().Description();
    for (const auto& row : *serial) ASSERT_TRUE(row.admissible);
    probe.independent_pairs = true;
    ASSERT_TRUE(observed.SupportsIndependentPair());
    const auto paired = AttemptRetainedIntervals(observed, inputs, 4);
    ASSERT_TRUE(paired) << paired.error().Description();
    EXPECT_EQ(probe.pair_calls, 1U);
    for (std::size_t row = 0; row < inputs.size(); ++row)
        EXPECT_TRUE(IntervalBitsAgree((*paired)[row], (*serial)[row]));

    // The third Endpoint is the final projection in both execution schedules.
    // Its rejected row must not acquire a speculative Dense decode error.
    unsigned endpoint_reads = 0;
    bool corrupt_dense = false;
    probe.after_read = [&](BufferHandle buffer, std::span<std::byte> bytes) {
        std::uint32_t flag = 0;
        if (buffer.value == endpoint_output.value && ++endpoint_reads == 3) {
            std::memcpy(bytes.data(), &flag, sizeof(flag));
        } else if (corrupt_dense && buffer.value == dense_output.value) {
            flag = 2;
            std::memcpy(bytes.data(), &flag, sizeof(flag));
        }
    };
    probe.independent_pairs = false;
    const auto rejected_serial = AttemptRetainedIntervals(observed, inputs, 4);
    ASSERT_TRUE(rejected_serial);
    ASSERT_FALSE(rejected_serial->front().admissible);
    ASSERT_TRUE(rejected_serial->back().admissible);
    EXPECT_EQ(rejected_serial->front().failure, sirius::core::CoupledStepFailure::Projection);
    endpoint_reads = 0;
    corrupt_dense = true;
    probe.independent_pairs = true;
    const auto rejected_pair = AttemptRetainedIntervals(observed, inputs, 4);
    ASSERT_TRUE(rejected_pair) << rejected_pair.error().Description();
    for (std::size_t row = 0; row < inputs.size(); ++row)
        EXPECT_TRUE(IntervalBitsAgree((*rejected_pair)[row], (*rejected_serial)[row]));

    // A malformed ACTIVE dense result remains a global error; a zero-input row
    // already inactive before speculation also keeps the original decode veto.
    probe.after_read = [&](BufferHandle buffer, std::span<std::byte> bytes) {
        if (buffer.value == dense_output.value) {
            const std::uint32_t malformed = 2;
            std::memcpy(bytes.data(), &malformed, sizeof(malformed));
        }
    };
    for (const bool initially_inactive : {false, true}) {
        auto malformed_input = inputs;
        if (initially_inactive) malformed_input[0].start.phase[5].valid = 0;
        const auto malformed = AttemptRetainedIntervals(observed, malformed_input, 4);
        ASSERT_FALSE(malformed);
        EXPECT_EQ(malformed.error().operation(), "read retained dense segment");
    }
    probe.after_read = {};
    const auto before_failed_read = observed.EndpointDenseStatistics();
    probe.after_pair = [&] { probe.fail_next = TransferProbeDevice::Failure::Read; };
    const auto failed_read = AttemptRetainedIntervals(observed, inputs, 4);
    ASSERT_FALSE(failed_read);
    EXPECT_EQ(failed_read.error().operation(), "retained transfer probe read");
    EXPECT_EQ(observed.EndpointDenseStatistics().submissions - before_failed_read.submissions, 1U);
    EXPECT_GT(observed.EndpointDenseStatistics().submit_wait_ms - before_failed_read.submit_wait_ms,
              0);
    probe.after_pair = {};
    EXPECT_EQ(device->BufferAllocationBytes(), resident);

    // Synthetic governor timings are explicit controls, not physical duration.
    // The actual device still calculates every private candidate and retry.
    probe.submission_ms = [](BufferHandle, std::uint32_t) { return 100.; };
    probe.pair_submit_ms = 1200;
    (void)observed.TakeSubmissionFeedback();
    RetainedTraceExecutor executor(observed, {}, 1000);
    sirius::core::KerrSchildFamily flat(sirius::core::KerrSchildParams::Minkowski());
    sirius::core::Lightray ray{};
    ray.position(1) = 5;
    ray.velocity(0) = -1;
    ray.velocity(1) = 1;
    ray.step_size = 1;
    sirius::core::Rk45CoupledState columns;
    columns.length_scale = columns.frequency_scale = 1;
    columns.tolerance = 1e-9;
    columns.variations[0].derivative(2) = .001;
    columns.variations[1].derivative(3) = .001;
    columns.variations[2].displacement(2) = 1;
    columns.variations[3].displacement(3) = 1;
    sirius::core::IntegratorConfig control;
    control.min_step = .01f;
    control.max_step = 2;
    sirius::core::Rk45CoupledComparison comparison;
    const auto pairs_before = probe.pair_calls;
    ASSERT_TRUE(executor.Step(ray, flat, control, columns, comparison));
    EXPECT_EQ(probe.pair_calls - pairs_before, 1U);
    EXPECT_EQ(executor.Statistics().paired_projection_retries, 1U);
    EXPECT_EQ(columns.central_stages, 42U);
    ASSERT_TRUE(executor.Step(ray, flat, control, columns, comparison));
    EXPECT_EQ(probe.pair_calls - pairs_before, 1U);  // Sticky budget one is fully serialized.
    EXPECT_EQ(device->BufferAllocationBytes(), resident);
    RecordProperty("device", device->Info().name);
    RecordProperty("scope",
                   "actual pair arithmetic/word equivalence and injected rejection/read/governor "
                   "controls; no frame-speed qualification");
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}

TEST_F(RetainedComputeTest, ProjectionReserveKeepsLogicalCohortsBounded) {
#ifdef SIRIUS_RETAINED_TESTS_AVAILABLE
    const auto inventory = EnumerateVulkanDevices();
    ASSERT_TRUE(inventory);
    const auto index = ResolveVulkanDeviceIndex(*inventory);
    ASSERT_TRUE(index);
    for (const std::size_t rays : {3U, 7U, 8U, 13U}) {
        SCOPED_TRACE(rays);
        // ComputeDevice owns its buffers until destruction. Each finite case
        // gets a separate device so its unchanged 8 MiB bound includes all rows.
        auto opened = CreateVulkanDevice(*index);
        ASSERT_TRUE(opened) << opened.error().Description();
        auto& cohort_device = **opened;
        ASSERT_TRUE(cohort_device.SetBufferAllocationLimit(8 * 1024 * 1024));
        TransferProbeDevice probe(cohort_device);
        const auto resident_before = cohort_device.BufferAllocationBytes();
        const auto required = RetainedCompute::RequiredAllocationBytes(probe, 2 * rays);
        ASSERT_TRUE(required);
        auto created = RetainedCompute::Create(probe, 2 * rays, false, rays == 3 ? 100 : 250);
        ASSERT_TRUE(created) << created.error().Description();
        auto& reserved = **created;
        EXPECT_EQ(cohort_device.BufferAllocationBytes(), resident_before + *required);
        const auto endpoint_input = probe.allocations[4].handle;
        std::promise<void> entered, release;
        auto dispatched = entered.get_future();
        auto released = release.get_future().share();
        probe.before_dispatch = [&] {
            if (probe.dispatch_calls == 1) {
                entered.set_value();
                released.wait();
            }
        };
        // Timing feedback is controlled at the device seam; numerical outputs
        // still come from the real stages. This makes queue limits deterministic.
        bool reduced = false;
        probe.submission_ms = [&](BufferHandle input, std::uint32_t rows) {
            if (rays == 3 && !reduced && input.value == endpoint_input.value && rows == 6) {
                reduced = true;
                return 160.0;
            }
            return 1.0;
        };
        RetainedTraceExecutor executor(reserved, {}, 1000, rays);
        auto held = std::async(std::launch::async, [&] {
            sirius::core::KerrSchildFamily metric(sirius::core::KerrSchildParams::Minkowski());
            sirius::core::CameraRay camera;
            camera.origin(1) = 5;
            camera.origin(2) = std::numbers::pi / 2;
            camera.direction(1) = 1;
            return executor.Launch(metric, 0, camera).has_value();
        });
        dispatched.wait();
        std::vector<std::future<void>> workers;
        // Extra queued rays prove that storage capacity does not also increase
        // logical gathering. Hold a real preceding dispatch while all requests
        // enter the queue; no elapsed-time assumption establishes it.
        const std::size_t row_count = rays == 3 ? 12 : rays + 1;
        for (std::size_t row = 0; row < row_count; ++row)
            workers.push_back(std::async(std::launch::async, [&, row] {
                sirius::core::KerrSchildFamily metric(sirius::core::KerrSchildParams::Minkowski());
                sirius::core::Lightray ray{};
                ray.position(1) = 5 + static_cast<double>(row);
                ray.velocity(0) = -1;
                ray.velocity(1) = 1;
                ray.step_size = 1;
                sirius::core::Rk45CoupledState columns;
                columns.length_scale = columns.frequency_scale = 1;
                columns.tolerance = 1e-9;
                columns.variations[0].derivative(2) = .001;
                columns.variations[1].derivative(3) = .001;
                columns.variations[2].displacement(2) = 1;
                columns.variations[3].displacement(3) = 1;
                sirius::core::IntegratorConfig config;
                config.min_step = .01f;
                config.max_step = 2;
                sirius::core::Rk45CoupledComparison comparison;
                executor.BeginTrace();
                const bool completed = executor.Step(ray, metric, config, columns, comparison);
                executor.EndTrace();
                EXPECT_TRUE(completed);
                EXPECT_EQ(ray.position(0), -1);
                EXPECT_EQ(ray.position(1), 6 + static_cast<double>(row));
                EXPECT_EQ(ray.velocity(0), -1);
                EXPECT_EQ(ray.velocity(1), 1);
                EXPECT_EQ(columns.central_stages, 21U);
                EXPECT_EQ(columns.variation_stages, 21U);
            }));
        RetainedTraceExecutorTestPeer::WaitForQueued(executor, row_count);
        release.set_value();
        EXPECT_TRUE(held.get());
        for (auto& worker : workers) worker.get();
        EXPECT_FALSE(executor.Error());
        const auto stats = executor.Statistics();
        EXPECT_EQ(stats.maximum_batch_rows, rays);
        EXPECT_EQ(stats.full_batches, rays == 3 ? 5U : 1U);
        EXPECT_EQ(stats.batch_row_counts[rays], rays == 3 ? 3U : 1U);
        EXPECT_EQ(stats.batch_row_counts[1], 2U);
        if (rays == 3) {
            EXPECT_EQ(stats.batch_row_counts[2], 1U);
        }
        EXPECT_EQ(stats.camera_rows, 1U);
        EXPECT_EQ(stats.interval_rows, row_count);
        EXPECT_EQ(stats.accepted_intervals, row_count);
        EXPECT_EQ(stats.batch_subdivisions, rays == 3 ? 1U : 0U);
        EXPECT_EQ(stats.safety_fallbacks, 0U);
        std::vector<std::uint32_t> projection_rows;
        for (const auto& submission : probe.submissions)
            if (submission.buffers[0].value == endpoint_input.value)
                projection_rows.push_back(submission.x);
        if (rays == 3) {
            // The first soft overshoot gives budgets 6->1->2->4->6.
            // Logical gathers 3,1,2,3,3 must recover pairing in the last cohort;
            // comparing fullness against budget four would incorrectly stall.
            EXPECT_EQ(projection_rows,
                      (std::vector<std::uint32_t>{3, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
                                                  2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 6, 6, 6}));
        } else {
            EXPECT_EQ(projection_rows,
                      (std::vector<std::uint32_t>{
                          static_cast<std::uint32_t>(rays), static_cast<std::uint32_t>(2 * rays),
                          static_cast<std::uint32_t>(2 * rays),
                          static_cast<std::uint32_t>(2 * rays), 1, 2, 2, 2}));
        }
        EXPECT_EQ(probe.allocations.size(), 2 * RetainedCompute::kStageCount);
        EXPECT_EQ(cohort_device.BufferAllocationBytes(), resident_before + *required);
    }
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

    // Calibrate only this controller regression's budgets against genuine
    // device stages. The original scientific fixtures above remain unchanged.
    // Reinitialization uses the same public state as the actual Step consumer.
    {
        sirius::core::KerrSchildFamily metric(sirius::core::KerrSchildParams::Kerr(1, .5));
        sirius::core::CameraRay camera;
        camera.origin(1) = 5;
        camera.origin(2) = 1.1;
        camera.direction(1) = .7;
        camera.direction(2) = .2;
        camera.direction(3) = .1;
        camera.direction = camera.direction / std::sqrt(.7 * .7 + .2 * .2 + .1 * .1);
        const auto launch = sirius::core::LaunchCameraRay(metric, .5, camera);
        ASSERT_TRUE(launch);
        sirius::core::Lightray candidate{};
        candidate.position = launch->position;
        candidate.velocity = launch->tangent;
        candidate.step_size = .5f;
        sirius::core::Rk45CoupledState columns;
        columns.variations = launch->variations;
        columns.length_scale = columns.frequency_scale = 1;
        columns.tolerance = 1e-3;
        sirius::core::IntegratorConfig control;
        control.min_step = 1e-6f;
        control.max_step = 2;
        control.abs_tolerance = control.rel_tolerance = 1e-3f;
        RetainedInitializeInput initialization;
        initialization.values[0] = RetainedValue::FromDouble(1);
        initialization.values[1] = RetainedValue::FromDouble(.5);
        initialization.values[2] = initialization.values[3] = RetainedValue::FromDouble(0);
        for (int axis = 0; axis < 4; ++axis) {
            initialization.values[4 + axis] = RetainedValue::FromDouble(candidate.position(axis));
            initialization.values[8 + axis] = RetainedValue::FromDouble(candidate.velocity(axis));
            for (std::size_t column = 0; column < 4; ++column) {
                initialization.values[12 + 8 * column + axis] =
                    RetainedValue::FromDouble(columns.variations[column].displacement(axis));
                initialization.values[16 + 8 * column + axis] =
                    RetainedValue::FromDouble(columns.variations[column].derivative(axis));
            }
        }
        initialization.values[44] = RetainedValue::FromDouble(1);
        const auto start = compute->Initialize({&initialization, 1});
        ASSERT_TRUE(start);
        ASSERT_TRUE(start->front().valid);
        RetainedEndpointInput projection;
        std::copy_n(initialization.values.begin(), 4, projection.values.begin());
        std::copy(start->front().phase.begin(), start->front().phase.end(),
                  projection.values.begin() + 4);
        projection.values[44] = RetainedValue::FromDouble(1);
        const auto projected = compute->Endpoint({&projection, 1});
        ASSERT_TRUE(projected);
        ASSERT_TRUE(projected->front().valid);
        RetainedIntervalInput interval;
        std::copy_n(initialization.values.begin(), 4, interval.metric.begin());
        interval.start = projected->front();
        interval.chart = 1;
        interval.interval = candidate.step_size;
        interval.control = {control, 1, 1, columns.tolerance, columns.column_scale};
        const auto calibration = AttemptRetainedIntervals(*compute, {&interval, 1});
        ASSERT_TRUE(calibration);
        const auto& output = calibration->front();
        ASSERT_TRUE(output.admissible);
        ASSERT_GT(output.error_ratio, 0);
        ASSERT_LT(output.embedded_projected_error_ratio, .5 * output.error_ratio);
        const double scale = output.error_ratio / .8;
        control.abs_tolerance *= static_cast<float>(scale);
        control.rel_tolerance *= static_cast<float>(scale);
        columns.tolerance *= scale;
        const float accepted_interval = candidate.step_size;
        sirius::core::Rk45CoupledComparison comparison;
        const auto growth_limits = executor.Statistics().dense_growth_limits;
        ASSERT_TRUE(executor.Step(candidate, metric, control, columns, comparison));
        ASSERT_GT(comparison.error_ratio, .79);
        ASSERT_LT(comparison.error_ratio, .81);
        EXPECT_FLOAT_EQ(candidate.step_size, accepted_interval);
        EXPECT_EQ(executor.Statistics().dense_growth_limits, growth_limits + 1);
        EXPECT_EQ(columns.central_stages, 21U);
        EXPECT_EQ(columns.variation_stages, 21U);
    }

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
    const auto rollbacks = executor.Statistics().tracer_rollbacks;
    // Mimic a rejected localized event: restore the original public state and
    // retry half the interval. The future phase must not contaminate this retry.
    executor.RejectLastInterval();
    EXPECT_EQ(executor.Statistics().tracer_rollbacks, rollbacks + 1);
    executor.RejectLastInterval();
    EXPECT_EQ(executor.Statistics().tracer_rollbacks, rollbacks + 1);
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
        const auto mixed_batches = mixed_after.batches - mixed_before.batches;
        EXPECT_GT(mixed_batches, 0U);
        // Prior feedback can shrink the batch limit until these two requests
        // fill it. The held trace forbids traces-ready, but permits a full batch
        // or the bounded timeout; no stopping occurs inside this window.
        EXPECT_EQ(mixed_after.coalescing_timeouts - mixed_before.coalescing_timeouts +
                      mixed_after.full_batches - mixed_before.full_batches,
                  mixed_batches);
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
    std::uint64_t scaled_rows = 0;
    for (const auto rows : timing.scaled_interval_bins) scaled_rows += rows;
    EXPECT_EQ(scaled_rows, timing.interval_measurements);
    EXPECT_LE(timing.interval_measurements, timing.interval_rows);
    EXPECT_GT(timing.interval_min, 0);
    EXPECT_GT(timing.scaled_interval_min, 0);
    EXPECT_LE(timing.interval_min, timing.interval_max);
    EXPECT_LE(timing.scaled_interval_min, timing.scaled_interval_max);
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
    ASSERT_EQ(guard_probe.allocations.size(), 2 * RetainedCompute::kStageCount);
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
    ASSERT_TRUE(serial_feedback.maximum_one_row_stage);
    EXPECT_EQ(*serial_feedback.maximum_one_row_stage, RetainedCompute::KernelStage::kDense);
    const auto consumed = guard_compute.TakeSubmissionFeedback();
    EXPECT_EQ(consumed.peak_ms, 0);
    EXPECT_EQ(consumed.maximum_rows, 0U);
    EXPECT_EQ(consumed.maximum_one_row_ms, 0);
    EXPECT_FALSE(consumed.maximum_one_row_stage);
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
            guard_compute, [&] { return guard_cancelled.load(); }, 1000, 1);
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
        EXPECT_EQ(guard_probe.allocations.size(), 2 * RetainedCompute::kStageCount);
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
            if (mode == GuardTiming::Irreducible || mode == GuardTiming::RetryIrreducible) {
                ASSERT_TRUE(guarded_executor.Error());
                EXPECT_EQ(guarded_executor.Error()->operation(), "dispatch retained renderer");
                // The paired endpoint peak is 1200ms. The terminal one-row
                // observation is transport at 1100ms, including a retry-only
                // overshoot after a preceding window with different stage identity.
                EXPECT_EQ(guarded_executor.Error()->detail(),
                          "single-row submission exceeded the safety duration: "
                          "stage=transport, active_rows=1, submit_wait_ms=1100, limit_ms=1000");
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

    // A failed camera launch has no accepted interval through which to export
    // timings. Preserve that stage in the error, including an adjacent timing
    // above the unchanged strict 1000ms guard which coarse decimals could hide.
    sirius::core::CameraRay guard_camera;
    guard_camera.origin(1) = 5;
    guard_camera.origin(2) = std::numbers::pi / 2;
    guard_camera.direction(1) = 1;
    const auto ray_camera_input = guard_probe.allocations[10].handle;
    for (const double duration :
         {1000.0, std::nextafter(1000.0, std::numeric_limits<double>::infinity())}) {
        guard_probe.submission_ms = [&, duration](BufferHandle input, std::uint32_t rows) {
            EXPECT_EQ(input.value, ray_camera_input.value);
            EXPECT_EQ(rows, 1U);
            return duration;
        };
        RetainedTraceExecutor guarded_executor(guard_compute, {}, 1000);
        const auto launch = guarded_executor.Launch(flat, 0, guard_camera);
        const bool permitted = duration == 1000;
        EXPECT_EQ(launch.has_value(), permitted);
        EXPECT_EQ(guarded_executor.Statistics().initialized_phases, 0U);
        EXPECT_EQ(guarded_executor.Statistics().interval_rows, 0U);
        if (permitted) {
            EXPECT_FALSE(guarded_executor.Error());
        } else {
            ASSERT_TRUE(guarded_executor.Error());
            EXPECT_EQ(guarded_executor.Error()->detail(),
                      "single-row submission exceeded the safety duration: "
                      "stage=ray_camera, active_rows=1, "
                      "submit_wait_ms=1000.0000000000001, limit_ms=1000");
            const auto calls = guard_probe.dispatch_calls;
            EXPECT_FALSE(guarded_executor.Launch(flat, 0, guard_camera));
            EXPECT_EQ(guard_probe.dispatch_calls, calls);
        }
    }
    ASSERT_NO_FATAL_FAILURE(CheckStickyErrorDrainsQueuedRequests(*device));

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
    RecordProperty("fma_fp32_enabled", static_cast<int>(device->Info().fma_fp32_enabled));
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
    std::size_t limiting_field = 41;
    EXPECT_DOUBLE_EQ(RetainedPhysicalError(first, second, control, &limiting_field), expected);
    EXPECT_EQ(limiting_field, 8U);
    EXPECT_DOUBLE_EQ(RetainedPhysicalError(second, first, control), expected);
    EXPECT_EQ(RetainedPhysicalError(first, first, control), 0);
    EXPECT_EQ(RetainedPhysicalError(first, first, control, &limiting_field), 0);
    EXPECT_EQ(limiting_field, 40U);
    first.fill(RetainedValue::FromDouble(0));
    second = first;
    first[0] = RetainedValue::FromDouble(1);
    control.integrator.abs_tolerance = control.integrator.rel_tolerance = 1;
    // One physical central component contributes 1/2; the norm is its RMS
    // over eight components, rather than that largest single contributor.
    EXPECT_DOUBLE_EQ(RetainedPhysicalError(first, second, control, &limiting_field),
                     .5 / std::sqrt(8.));
    EXPECT_EQ(limiting_field, 40U);
    second[8].valid = 0;
    EXPECT_EQ(RetainedPhysicalError(first, second, control, &limiting_field),
              std::numeric_limits<double>::infinity());
    EXPECT_EQ(limiting_field, 41U);
#else
    GTEST_SKIP() << "Retained compute build tools unavailable";
#endif
}
}  // namespace
