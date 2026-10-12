#include "sirius/backend/retained_compute.h"

#include "sirius/backend/cpu/geodesic_tracer.h"
#include "sirius/backend/retained_integrator.h"
#include "sirius/backend/retained_trace_executor.h"
#include "sirius/base/sha256.h"
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
#include <string_view>
#include <vector>

namespace sirius::backend::retained_program { inline constexpr std::size_t kTransportRowWords = 4699; }
namespace {
using namespace sirius::backend;
class PreparationProbeDevice final : public ComputeDevice {
  public:
    struct Write {
        BufferHandle buffer;
        std::size_t bytes;
        std::uint32_t header;
    };
    DeviceInfo info{.kind = DeviceKind::kSoftware};
    std::vector<std::vector<std::byte>> buffers;
    std::vector<std::uint64_t> epochs;
    std::vector<BufferHandle> read_buffers;
    std::vector<Write> writes;
    std::vector<std::uint32_t> kernels;
    std::vector<std::vector<std::uint32_t>> loaded_codes;
    std::vector<std::size_t> failed_writes;
    std::function<void()> after_write, after_dispatch, after_read;
    DispatchTiming observation{.submit_wait_ms = 5000,
                               .pipeline_setup_ms = 1,
                               .command_setup_ms = 2,
                               .cleanup_ms = 3,
                               .total_ms = 5006,
                               .pipeline_created = true};
    bool fail_dispatch = false;
    bool independent_pair = false;
    bool track_epochs = false;
    std::size_t expected_bindings = 2;
    std::function<void(std::span<const BufferHandle>)> model_output;
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
        epochs.push_back(1);
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
        ++epochs[buffer.value];
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
        read_buffers.push_back(buffer);
        if (after_read) after_read();
        return {};
    }
    sirius::base::Expected<void> Dispatch(KernelHandle kernel,
                                          std::span<const BufferHandle> bindings, std::uint32_t x,
                                          std::uint32_t y, std::uint32_t z,
                                          DispatchTiming* timing) override {
        if (expected_bindings != 0) {
            EXPECT_EQ(bindings.size(), expected_bindings);
        }
        EXPECT_EQ(x, expected_groups_x);
        EXPECT_EQ(y, 1U);
        EXPECT_EQ(z, 1U);
        kernels.push_back(kernel.value);
        for (std::size_t i = 0; i < bindings.size(); ++i)
            if (std::none_of(bindings.begin(), bindings.begin() + i,
                             [&](auto bound) { return bound.value == bindings[i].value; }))
                ++epochs[bindings[i].value];
        if (timing) *timing = observation;
        if (after_dispatch) after_dispatch();
        if (fail_dispatch) return Failure("dispatch");
        std::uint32_t rows = 0;
        std::memcpy(&rows, buffers[bindings[0].value].data(), sizeof(rows));
        // A physical control returns a declined, zero-stage row; only its
        // actual submission/feedback ownership is tested against this model.
        if (rows != 0) std::ranges::fill(buffers[bindings[1].value], std::byte{0});
        if (model_output) model_output(bindings);
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
    std::optional<std::uint64_t> BufferMutationEpoch(BufferHandle buffer) const noexcept override {
        if (!track_epochs || buffer.value >= epochs.size()) return std::nullopt;
        return epochs[buffer.value];
    }

  private:
    static sirius::base::Expected<void> Failure(const char* operation) {
        return sirius::base::Fail(sirius::base::ErrorDomain::kDevice,
                                  std::string("preparation probe ") + operation,
                                  "injected returned failure");
    }
};

void ModelAcceptedStageSevenRows(PreparationProbeDevice& device) {
    device.model_output = [&](std::span<const BufferHandle> bindings) {
        if (bindings[0].value != 2) return;
        auto& bytes = device.buffers[bindings[1].value];
        std::vector<std::uint32_t> words(bytes.size() / 4);
        for (std::size_t row = 0; row < 2; ++row) {
            auto* out = words.data() + row * retained_program::kTransportRowWords;
            out[0] = 1;
            out[1] = 7;
            out[2] = 1;
            for (std::size_t group = 0; group < 4; ++group)
                for (std::size_t i = 0; i < 40; ++i) out[4 + group * 200 + i * 5 + 4] = 1;
            for (std::size_t stage = 0; stage < 7; ++stage)
                for (std::size_t i = 0; i < 40; ++i) out[1004 + stage * 200 + i * 5 + 4] = 1;
        }
        std::memcpy(bytes.data(), words.data(), bytes.size());
    };
}

// Only the host/device seam is modeled here. Shader replay and numerical
// equivalence require separate executions through the real Vulkan provider.
TEST(RetainedComputeAdmission, StageSevenGrantsPreservePayloadAndConsumeTokens) {
    PreparationProbeDevice device;
    device.info.rounds_fp32_to_nearest = true;
    device.expected_bindings = 3;
    device.expected_groups_x = 2;
    auto created = RetainedCompute::Create(device, 2, false, 250, true);
    ASSERT_TRUE(created);
    auto& compute = **created;
    ASSERT_TRUE(compute.StageSevenOwnershipEnabled());
    ASSERT_EQ(device.buffers.size(), 15U);
    EXPECT_EQ(device.loads, 9U);
    const auto extra = RetainedCompute::RequiredStageSevenOwnerAllocationBytes(device, 2);
    ASSERT_TRUE(extra);
    EXPECT_EQ(*extra, device.buffers.back().size());
    EXPECT_EQ(device.BufferAllocationBytes(), RetainedCompute::RequiredBufferBytes(2) + *extra);
    device.model_output = [&](std::span<const BufferHandle> bindings) {
        if (bindings[0].value != 2) return;
        auto& bytes = device.buffers[bindings[1].value];
        std::vector<std::uint32_t> words(bytes.size() / 4);
        for (std::size_t row = 0; row < 2; ++row) {
            auto* out = words.data() + row * retained_program::kTransportRowWords;
            out[0] = 1;
            out[1] = 7;
            out[2] = 1;
            for (std::size_t group = 0; group < 4; ++group)
                for (std::size_t i = 0; i < 40; ++i) out[4 + group * 200 + i * 5 + 4] = 1;
            for (std::size_t stage = 0; stage < 7; ++stage)
                for (std::size_t i = 0; i < 40; ++i) out[1004 + stage * 200 + i * 5 + 4] = 1;
        }
        std::memcpy(bytes.data(), words.data(), bytes.size());
    };
    std::array<RetainedStepInput, 2> packets;
    for (auto& packet : packets) {
        packet.values.fill(RetainedValue::FromDouble(0));
        packet.values[0] = RetainedValue::FromDouble(1);
        packet.values[44] = RetainedValue::FromDouble(1);
        packet.values[45] = RetainedValue::FromDouble(0.01);
    }
    auto step = compute.Step(packets);
    ASSERT_TRUE(step);
    std::array<RetainedEndpointInput, 2> upper;
    for (std::size_t row = 0; row < 2; ++row) {
        std::copy_n(packets[row].values.begin(), 4, upper[row].values.begin());
        std::copy((*step)[row].fifth.begin(), (*step)[row].fifth.end(),
                  upper[row].values.begin() + 4);
        upper[row].values[44] = packets[row].values[44];
    }
    std::array<StageSevenToken, 2> tokens{compute.LastStepOwner(0), {}};
    const auto sidecar_before = device.buffers.back();
    ASSERT_TRUE(compute.Endpoint(upper, nullptr, tokens));
    EXPECT_EQ(device.kernels.back(), 8U);
    std::uint32_t grants = 0;
    std::memcpy(&grants, device.buffers.back().data() + 16 * 4, 4);
    EXPECT_EQ(grants, 1U);
    EXPECT_TRUE(std::equal(sidecar_before.begin() + 17 * 4, sidecar_before.end(),
                           device.buffers.back().begin() + 17 * 4));
    device.expected_bindings = 2;
    ASSERT_TRUE(compute.Endpoint(upper, nullptr, tokens));
    EXPECT_EQ(device.kernels.back(), 2U);

    // A changed momentum leaves all nine reusable geometry keys unchanged.
    // The full high-packet contract must still decline this token.
    device.expected_bindings = 3;
    ASSERT_TRUE(compute.Step(packets));
    tokens = {compute.LastStepOwner(0), {}};
    upper[0].values[8] = RetainedValue::FromDouble(1);
    device.expected_bindings = 2;
    ASSERT_TRUE(compute.Endpoint(upper, nullptr, tokens));
    EXPECT_EQ(device.kernels.back(), 2U);
    upper[0].values[8] = RetainedValue::FromDouble(0);

    // Identical packets cannot turn a token for physical row one into row zero.
    device.expected_bindings = 3;
    ASSERT_TRUE(compute.Step(packets));
    tokens = {compute.LastStepOwner(1), {}};
    device.expected_bindings = 2;
    ASSERT_TRUE(compute.Endpoint(upper, nullptr, tokens));
    EXPECT_EQ(device.kernels.back(), 2U);
}

TEST(RetainedComputeAdmission, StageSevenFailureCannotPublishPartialOrOldIssuance) {
    PreparationProbeDevice device;
    device.info.rounds_fp32_to_nearest = true;
    device.expected_bindings = 3;
    device.expected_groups_x = 2;
    auto created = RetainedCompute::Create(device, 2, false, 250, true);
    ASSERT_TRUE(created);
    auto& compute = **created;
    bool malformed_second = false;
    device.model_output = [&](std::span<const BufferHandle> bindings) {
        if (bindings[0].value != 2) return;
        auto& bytes = device.buffers[bindings[1].value];
        std::vector<std::uint32_t> words(bytes.size() / 4);
        for (std::size_t row = 0; row < 2; ++row) {
            auto* out = words.data() + row * retained_program::kTransportRowWords;
            out[0] = 1;
            out[1] = malformed_second && row == 1 ? 6 : 7;
            out[2] = 1;
            for (std::size_t group = 0; group < 4; ++group)
                for (std::size_t i = 0; i < 40; ++i) out[4 + group * 200 + i * 5 + 4] = 1;
            for (std::size_t stage = 0; stage < 7; ++stage)
                for (std::size_t i = 0; i < 40; ++i) out[1004 + stage * 200 + i * 5 + 4] = 1;
        }
        std::memcpy(bytes.data(), words.data(), bytes.size());
    };
    std::array<RetainedStepInput, 2> packets;
    for (auto& packet : packets) packet.values.fill(RetainedValue::FromDouble(0));
    auto stepped = compute.Step(packets);
    ASSERT_TRUE(stepped);
    const auto old = compute.LastStepOwner(0);
    malformed_second = true;
    EXPECT_FALSE(compute.Step(packets));
    std::array<RetainedEndpointInput, 2> upper;
    for (auto& row : upper) row.values.fill(RetainedValue::FromDouble(0));
    std::array<StageSevenToken, 2> tokens{compute.LastStepOwner(0), {}};
    device.expected_bindings = 2;
    ASSERT_TRUE(compute.Endpoint(upper, nullptr, tokens));
    EXPECT_EQ(device.kernels.back(), 2U);
    tokens = {old, {}};
    ASSERT_TRUE(compute.Endpoint(upper, nullptr, tokens));
    EXPECT_EQ(device.kernels.back(), 2U);

    malformed_second = false;
    device.expected_bindings = 3;
    ASSERT_TRUE(compute.Step(packets));
    tokens = {compute.LastStepOwner(0), {}};
    device.failed_writes.push_back(device.writes.size() + 1);
    DispatchTiming timing = device.observation;
    EXPECT_FALSE(compute.Endpoint(upper, &timing, tokens));
    EXPECT_EQ(timing.total_ms, 0);
    device.expected_bindings = 2;
    ASSERT_TRUE(compute.Endpoint(upper, nullptr, tokens));
    EXPECT_EQ(device.kernels.back(), 2U);
    device.expected_bindings = 3;
    ASSERT_TRUE(compute.Step(packets));
    tokens = {compute.LastStepOwner(0), {}};
    ASSERT_TRUE(compute.Endpoint(upper, nullptr, tokens));
    EXPECT_EQ(device.kernels.back(), 8U);
}

TEST(RetainedComputeAdmission, StageSevenRequestKeepsStockUnsupportedArithmetic) {
    for (const bool rte : {false, true}) {
        PreparationProbeDevice device;
        device.info.rounds_fp32_to_nearest = rte;
        device.info.preserves_fp32_denormals = rte;
        auto created = RetainedCompute::Create(device, 2, false, 250, true);
        ASSERT_TRUE(created);
        EXPECT_FALSE((*created)->StageSevenOwnershipEnabled());
        EXPECT_EQ(device.buffers.size(), 14U);
        EXPECT_EQ(device.loads, 7U);
        auto extra = RetainedCompute::RequiredStageSevenOwnerAllocationBytes(device, 2);
        ASSERT_TRUE(extra);
        EXPECT_EQ(*extra, 0U);
    }
}

TEST(RetainedComputeAdmission, SealedTablesReuseVerificationAndPreserveChangedTails) {
    PreparationProbeDevice device;
    device.info.rounds_fp32_to_nearest = true;
    device.track_epochs = true;
    device.expected_bindings = 0;
    device.expected_groups_x = 2;
    device.independent_pair = true;
    auto created = RetainedCompute::Create(device, 2, false, 250, true);
    ASSERT_TRUE(created);
    auto& compute = **created;
    EXPECT_EQ(device.loads, 11U);
    ModelAcceptedStageSevenRows(device);
    std::array<RetainedStepInput, 2> packets;
    for (auto& packet : packets) packet.values.fill(RetainedValue::FromDouble(0));
    std::array<RetainedEndpointInput, 2> upper;
    for (auto& packet : upper) packet.values.fill(RetainedValue::FromDouble(0));
    ASSERT_TRUE(compute.Step(packets));
    EXPECT_EQ(device.kernels.back(), 9U);
    const auto transport_reads = std::ranges::count(device.read_buffers, 2U, &BufferHandle::value);
    EXPECT_EQ(transport_reads, 1);
    ASSERT_TRUE(compute.Step(packets));
    EXPECT_EQ(std::ranges::count(device.read_buffers, 2U, &BufferHandle::value), transport_reads);
    std::array<StageSevenToken, 2> tokens{compute.LastStepOwner(0), {}};
    ASSERT_TRUE(compute.Endpoint(upper, nullptr, tokens));
    EXPECT_EQ(device.kernels.back(), 10U);
    EXPECT_EQ(std::ranges::count(device.read_buffers, 4U, &BufferHandle::value), 1);
    ASSERT_TRUE(compute.Step(packets));
    tokens = {compute.LastStepOwner(0), {}};
    std::array<RetainedDenseInput, 2> dense{};
    ASSERT_TRUE(compute.EndpointAndDense(upper, dense, tokens));
    ASSERT_GE(device.kernels.size(), 2U);
    EXPECT_EQ(device.kernels[device.kernels.size() - 2], 10U);
    EXPECT_EQ(std::ranges::count(device.read_buffers, 4U, &BufferHandle::value), 1);

    ASSERT_TRUE(compute.Step(packets));
    tokens = {compute.LastStepOwner(0), {}};
    // Even an identical external input write retires the issued value token.
    const auto unchanged = device.buffers[2];
    ASSERT_TRUE(device.WriteBuffer(BufferHandle{2}, unchanged));
    ASSERT_TRUE(compute.Endpoint(upper, nullptr, tokens));
    EXPECT_EQ(device.kernels.back(), 2U);

    // Edit the actual warmed program, using the device API. A later prefix
    // upload must preserve this tail and route that stage through generic code.
    for (const auto input : {BufferHandle{2}, BufferHandle{4}}) {
        auto changed = device.buffers[input.value];
        changed.back() ^= std::byte{1};
        ASSERT_TRUE(device.WriteBuffer(input, changed));
        ASSERT_TRUE(compute.Step(packets));
        if (input.value == 2) {
            EXPECT_EQ(device.kernels.back(), 1U);
            EXPECT_EQ(device.buffers[2].back(), changed.back());
            EXPECT_EQ(std::ranges::count(device.read_buffers, 2U, &BufferHandle::value), 2);
            // Restore through an admitted write; no CPU upload repairs it.
            changed.back() ^= std::byte{1};
            ASSERT_TRUE(device.WriteBuffer(input, changed));
            ASSERT_TRUE(compute.Step(packets));
            EXPECT_EQ(device.kernels.back(), 9U);
        } else {
            tokens = {compute.LastStepOwner(0), {}};
            ASSERT_TRUE(compute.EndpointAndDense(upper, dense, tokens));
            EXPECT_EQ(device.kernels[device.kernels.size() - 2], 2U);
            EXPECT_EQ(device.buffers[4].back(), changed.back());
            EXPECT_EQ(std::ranges::count(device.read_buffers, 4U, &BufferHandle::value), 2);
        }
    }
}

TEST(RetainedComputeAdmission, SealedOwnerMutationRetiresReuseAndLateInputEditsRejectResults) {
    for (const bool mutate_owner : {false, true}) {
        SCOPED_TRACE(mutate_owner);
        PreparationProbeDevice device;
        device.info.rounds_fp32_to_nearest = true;
        device.track_epochs = true;
        device.expected_bindings = 0;
        device.expected_groups_x = 2;
        auto created = RetainedCompute::Create(device, 2, false, 250, true);
        ASSERT_TRUE(created);
        auto& compute = **created;
        ModelAcceptedStageSevenRows(device);
        std::array<RetainedStepInput, 2> packets;
        for (auto& packet : packets) packet.values.fill(RetainedValue::FromDouble(0));
        ASSERT_TRUE(compute.Step(packets));
        const BufferHandle target{mutate_owner ? 14U : 2U};
        bool injected = false;
        device.after_read = [&] {
            if (injected) return;
            injected = true;
            const auto bytes = device.buffers[target.value];
            ASSERT_TRUE(device.WriteBuffer(target, bytes));
        };
        EXPECT_FALSE(compute.Step(packets));
        EXPECT_TRUE(injected);
        device.after_read = {};
        ASSERT_TRUE(compute.Step(packets));
        EXPECT_EQ(device.kernels.back(), mutate_owner ? 1U : 9U);
        EXPECT_EQ(compute.StageSevenOwnershipEnabled(), !mutate_owner);
    }
}

TEST(RetainedComputeAdmission, SealedDispatchFailureRechecksOwnerPrefixBeforeRecovery) {
    for (const unsigned mode : {0U, 1U, 2U}) {
        SCOPED_TRACE(mode);
        const bool corrupt_prefix = mode != 0;
        PreparationProbeDevice device;
        device.info.rounds_fp32_to_nearest = true;
        device.track_epochs = true;
        device.expected_bindings = 0;
        device.expected_groups_x = 2;
        auto created = RetainedCompute::Create(device, 2, false, 250, true);
        ASSERT_TRUE(created);
        auto& compute = **created;
        ModelAcceptedStageSevenRows(device);
        std::array<RetainedStepInput, 2> packets;
        for (auto& packet : packets) packet.values.fill(RetainedValue::FromDouble(0));
        ASSERT_TRUE(compute.Step(packets));
        device.fail_dispatch = true;
        device.after_dispatch = [&] {
            // An admitted failed dispatch may have changed any bound bytes.
            if (corrupt_prefix) device.buffers.back()[17 * 4] ^= std::byte{1};
            if (mode == 2) {
                const auto input = device.buffers[2];
                ASSERT_TRUE(device.WriteBuffer(BufferHandle{2}, input));
                ASSERT_TRUE(device.WriteBuffer(BufferHandle{2}, input));
            }
        };
        EXPECT_FALSE(compute.Step(packets));
        device.fail_dispatch = false;
        device.after_dispatch = {};
        const auto reads = std::ranges::count(device.read_buffers, 14U, &BufferHandle::value);
        ASSERT_TRUE(compute.Step(packets));
        EXPECT_EQ(std::ranges::count(device.read_buffers, 14U, &BufferHandle::value), reads + 1);
        EXPECT_EQ(device.kernels.back(), corrupt_prefix ? 1U : 9U);
        EXPECT_EQ(compute.StageSevenOwnershipEnabled(), !corrupt_prefix);
    }
}


}
