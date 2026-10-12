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
        EXPECT_EQ(bindings.size(), expected_bindings);
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

  private:
    static sirius::base::Expected<void> Failure(const char* operation) {
        return sirius::base::Fail(sirius::base::ErrorDomain::kDevice,
                                  std::string("preparation probe ") + operation,
                                  "injected returned failure");
    }
};

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
            out[0] = 1; out[1] = 7; out[2] = 1;
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
        std::copy((*step)[row].fifth.begin(), (*step)[row].fifth.end(), upper[row].values.begin() + 4);
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
            out[0] = 1; out[1] = malformed_second && row == 1 ? 6 : 7; out[2] = 1;
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


}
