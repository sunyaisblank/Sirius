#include "sirius/backend/retained_compute.h"
#include "sirius/base/sha256.h"
#include "sirius/core/twofold.h"
#include "support/retained_transport/reference_cases.h"
#include "support/retained_transport/schwarzschild_reference_cases.h"
#include "retained_owner_layout.h"
#include "runtime_sensitivity.h"
#include "runtime_module_sha.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <stdexcept>
#include <vector>
namespace sirius::backend::retained_program {
inline constexpr std::size_t kTransportRowWords = 4699;
inline constexpr std::size_t kEndpointRowWords = 3829;
}
namespace {
using namespace sirius::backend;
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
    explicit TransferProbeDevice(ComputeDevice& device, bool decline_rte32 = false)
        : device_(device) {
        if (decline_rte32) {
            info_ = device.Info();
            info_->rounds_fp32_to_nearest = false;
        }
    }
    Failure fail_next = Failure::None;
    std::vector<KernelHandle> loaded_handles, dispatched_handles;
    std::vector<Allocation> allocations;
    std::vector<Transfer> writes, reads;
    std::vector<Submission> submissions;
    std::uint64_t dispatch_calls = 0, forwarded_dispatches = 0;
    // Optional control-test observations. Physical outputs still come from the
    // actual device; injected timing never sleeps or claims hardware duration.
    std::function<void(std::span<const std::uint32_t>)> after_load;
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
    std::function<void(BufferHandle, std::span<std::byte>)> before_write;

    std::optional<std::uint64_t> BufferMutationEpoch(BufferHandle buffer) const noexcept override { return device_.BufferMutationEpoch(buffer); }
    const DeviceInfo& Info() const noexcept override { return info_ ? *info_ : device_.Info(); }
    sirius::base::Expected<KernelHandle> LoadKernel(std::span<const std::uint32_t> code) override {
        auto loaded = device_.LoadKernel(code);
        if (loaded) loaded_handles.push_back(*loaded);
        if (loaded && after_load) after_load(code);
        return loaded;
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
        if (before_write) {
            std::vector<std::byte> changed(data.begin(), data.end());
            before_write(buffer, changed);
            return device_.WriteBuffer(buffer, changed);
        }
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
        dispatched_handles.push_back(kernel);
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
    std::optional<DeviceInfo> info_;
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

::testing::AssertionResult SchwarzschildStepAgrees(
    const RetainedStepOutput& output, const sirius::test::retained_schwarzschild::Case& fixture) {
    if (!output.valid || output.stages != 7) return ::testing::AssertionFailure();
    const bool flat = std::string_view(fixture.name) == "frozen-analytic-flat";
    if (output.rhs_valid == flat) return ::testing::AssertionFailure() << "RHS classification";
    const std::array records{&output.fifth, &output.fourth, &output.increment, &output.error};
    const auto get = [&](std::size_t index) -> const RetainedValue& {
        return index < 160 ? (*records[index / 40])[index % 40]
                           : output.rhs[(index - 160) / 40][index % 40];
    };
    // The unchanged flat oracle has no performed RHS records. Curved cases
    // check all 160 final and 280 independently differentiated RHS fields.
    const std::size_t limit = flat ? 160 : 440;
    for (std::size_t begin = 0; begin < limit; begin += 4) {
        long double error = 0, scale = 0;
        for (std::size_t axis = 0; axis < 4; ++axis) {
            const auto index = begin + axis;
            const auto& value = get(index);
            if (!value.IsRepresented()) return ::testing::AssertionFailure();
            const auto reference = fixture.reference[index];
            const sirius::core::Twofold exact(reference.high, reference.low);
            const auto center = sirius::core::Twofold(value.high) +
                                sirius::core::Twofold(value.low) +
                                sirius::core::Twofold(value.tail);
            const double difference = std::abs((center - exact).Rounded());
            constexpr double epsilon = std::numeric_limits<double>::epsilon();
            const double rounding =
                128 * epsilon * epsilon * (std::abs(reference.high) + std::abs(value.high));
            if (difference > value.radius + fixture.gap[index] + rounding)
                return ::testing::AssertionFailure()
                       << "Schwarzschild component " << index << " difference=" << difference
                       << " radius=" << value.radius;
            error = std::max(
                error, std::abs((static_cast<long double>(value.high) + value.low + value.tail) -
                                (static_cast<long double>(reference.high) + reference.low)));
            const auto scale_index = index >= 120 && index < 160 ? index - 40 : index;
            const auto scale_reference = fixture.reference[scale_index];
            scale = std::max(scale, std::abs(static_cast<long double>(scale_reference.high) +
                                             scale_reference.low));
        }
        if (error > 1e-11L * scale)
            return ::testing::AssertionFailure() << "Schwarzschild group " << begin;
    }
    return ::testing::AssertionSuccess();
}

std::filesystem::path artifacts;
std::vector<std::uint32_t> ReadWords(ComputeDevice& device, BufferHandle handle, std::size_t bytes) {
    std::vector<std::uint32_t> result(bytes / 4);
    auto read = device.ReadBuffer(handle, std::as_writable_bytes(std::span(result)));
    if (!read) throw std::runtime_error(read.error().Description());
    return result;
}
std::vector<std::uint32_t> ReadAuthority(const char* name) {
    const auto path = std::filesystem::path(__FILE__).parent_path() / name;
    const auto bytes = std::filesystem::file_size(path);
    std::vector<std::uint32_t> result(bytes / 4);
    std::ifstream file(path, std::ios::binary);
    file.read(reinterpret_cast<char*>(result.data()), bytes);
    if (!file) throw std::runtime_error("authority read failed");
    return result;
}
void Save(const std::string& name, std::span<const std::byte> bytes) {
    std::ofstream file(artifacts / name, std::ios::binary);
    file.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
    if (!file) throw std::runtime_error("artifact write failed");
}
std::string Hash(std::span<const std::byte> bytes) {
    const auto digest = sirius::base::Sha256Hex(std::span(
        reinterpret_cast<const std::uint8_t*>(bytes.data()), bytes.size()));
    if (!digest) throw std::runtime_error(digest.error());
    return *digest;
}
TEST(SealedOwnerStreamProbe, PositiveReplayAndMutationRetirementPreserveWords) {
    constexpr std::size_t capacity = 24, rows = 5;
    constexpr std::uint64_t limit = 8 * 1024 * 1024;
    const auto inventory = EnumerateVulkanDevices();
    ASSERT_TRUE(inventory) << inventory.error().Description();
    ASSERT_EQ(inventory->size(), 1U);
    const auto opened = CreateVulkanDevice(0);
    ASSERT_TRUE(opened) << opened.error().Description();
    auto& device = **opened;
    auto expected_info = inventory->front();
    expected_info.fma_fp32_enabled = device.Info().fma_fp32_enabled;
    ASSERT_EQ(device.Info(), expected_info);
    ASSERT_EQ(device.Info().kind, DeviceKind::kSoftware);
    ASSERT_TRUE(RetainedUsesPortableArithmetic(device.Info()));
    ASSERT_TRUE(device.Info().rounds_fp32_to_nearest);
    ASSERT_TRUE(device.SetBufferAllocationLimit(limit));
    const auto base = RetainedCompute::RequiredAllocationBytes(device, capacity);
    const auto extra = RetainedCompute::RequiredStageSevenOwnerAllocationBytes(device, capacity);
    ASSERT_TRUE(base); ASSERT_TRUE(extra); ASSERT_GT(*extra, 0U);
    ASSERT_LE(2 * *base + *extra, limit);
    ASSERT_EQ(device.BufferAllocationBytes(), 0U);
    TransferProbeDevice ordinary(device), owned(device);
    std::vector<std::string> ordinary_hashes, owned_hashes;
    ordinary.after_load = [&](auto words) { ordinary_hashes.push_back(Hash(std::as_bytes(words))); };
    owned.after_load = [&](auto words) { owned_hashes.push_back(Hash(std::as_bytes(words))); };
    auto stock = RetainedCompute::Create(ordinary, capacity);
    auto candidate = RetainedCompute::Create(owned, capacity, false, 250, true);
    ASSERT_TRUE(stock) << stock.error().Description();
    ASSERT_TRUE(candidate) << candidate.error().Description();
    ASSERT_TRUE((*candidate)->StageSevenOwnershipEnabled());
    ASSERT_FALSE((*stock)->StageSevenOwnershipEnabled());
    ASSERT_EQ(ordinary.allocations.size(), 14U); ASSERT_EQ(owned.allocations.size(), 15U);
    ASSERT_EQ(ordinary.loaded_handles.size(), 7U); ASSERT_EQ(owned.loaded_handles.size(), 11U);
    ASSERT_EQ(device.BufferAllocationBytes(), 2 * *base + *extra);
    ASSERT_EQ(ordinary_hashes.size(), 7U); ASSERT_EQ(owned_hashes.size(), 11U);
    for (std::size_t i = 0; i < ordinary_hashes.size(); ++i)
        ASSERT_EQ(owned_hashes[i], ordinary_hashes[i]);
    ASSERT_EQ(owned_hashes[7], kOwnerTransportSha);
    ASSERT_EQ(owned_hashes[8], kOwnerEndpointSha);
    ASSERT_EQ(owned_hashes[9], kSealedTransportSha);
    ASSERT_EQ(owned_hashes[10], kSealedEndpointSha);
    RecordProperty("device_name", device.Info().name);
    RecordProperty("actual_allocation_bytes", std::to_string(device.BufferAllocationBytes()));
    RecordProperty("actual_owner_allocation_bytes", std::to_string(*extra));
    RecordProperty("owner_provisioning_ms_outside_call_samples", std::to_string((*candidate)->StageSevenProvisioningMs()));
    for (std::size_t i = 0; i < owned_hashes.size(); ++i)
        RecordProperty("loaded_owner_role_" + std::to_string(i), owned_hashes[i]);
    const auto owner_handle = owned.allocations.back().handle;
    const auto owner_bytes = owned.allocations.back().bytes;
    ASSERT_EQ(owner_bytes, 4 * (retained_owner::kPrefixWords + 1 + capacity * retained_owner::kRowWords));
    std::vector<std::uint32_t> published;
    std::size_t endpoint_bindings = 0;
    owned.after_dispatch = [&](std::span<const BufferHandle> buffers, const DispatchTiming*) {
        if (buffers[0].value == owned.allocations[2].handle.value) {
            EXPECT_EQ(buffers.size(), 3U);
            ASSERT_EQ(buffers[2].value, owner_handle.value);
            EXPECT_EQ(owned.dispatched_handles.back().value, owned.loaded_handles[9].value);
        } else if (buffers[0].value == owned.allocations[4].handle.value) {
            endpoint_bindings = buffers.size();
            if (buffers.size() == 3) {
                EXPECT_EQ(buffers[2].value, owner_handle.value);
                EXPECT_EQ(owned.dispatched_handles.back().value, owned.loaded_handles[10].value);
            } else {
                EXPECT_EQ(buffers.size(), 2U);
                EXPECT_EQ(owned.dispatched_handles.back().value, owned.loaded_handles[2].value);
            }
        }
    };
    std::array<RetainedStepInput, rows> packets{
        std::bit_cast<RetainedStepInput>(sirius::test::retained_transport::cases[0].input),
        std::bit_cast<RetainedStepInput>(sirius::test::retained_transport::cases[1].input),
        std::bit_cast<RetainedStepInput>(sirius::test::retained_schwarzschild::cases[16].input),
        std::bit_cast<RetainedStepInput>(sirius::test::retained_schwarzschild::cases[17].input),
        std::bit_cast<RetainedStepInput>(sirius::test::retained_transport::cases.back().input)};
    ordinary.capture_output = ordinary.allocations[3].handle;
    auto reference = (*stock)->Step(packets);
    ASSERT_TRUE(reference) << reference.error().Description();
    ASSERT_EQ(ordinary.readbacks.size(), 1U);
    const auto stock_step_words = ordinary.readbacks.back();
    std::array<RetainedEndpointInput, rows> upper, lower;
    for (std::size_t row = 0; row < rows; ++row) {
        std::copy_n(packets[row].values.begin(), 4, upper[row].values.begin());
        std::copy((*reference)[row].fifth.begin(), (*reference)[row].fifth.end(), upper[row].values.begin() + 4);
        upper[row].values[44] = packets[row].values[44];
        lower[row] = upper[row];
        std::copy((*reference)[row].fourth.begin(), (*reference)[row].fourth.end(), lower[row].values.begin() + 4);
    }
    const auto fresh = [&]() -> std::array<StageSevenToken, rows> {
        owned.capture_output = owned.allocations[3].handle; owned.readbacks.clear();
        const auto step = (*candidate)->Step(packets);
        if (!step) throw std::runtime_error(step.error().Description());
        EXPECT_EQ(owned.readbacks.size(), 1U);
        EXPECT_EQ(owned.readbacks.back(), stock_step_words);
        EXPECT_TRUE(StepAgrees((*step)[0], sirius::test::retained_transport::cases[0]));
        EXPECT_TRUE(StepAgrees((*step)[1], sirius::test::retained_transport::cases[1]));
        EXPECT_TRUE(SchwarzschildStepAgrees((*step)[2], sirius::test::retained_schwarzschild::cases[16]));
        EXPECT_TRUE(SchwarzschildStepAgrees((*step)[3], sirius::test::retained_schwarzschild::cases[17]));
        EXPECT_TRUE(StepAgrees((*step)[4], sirius::test::retained_transport::cases.back()));
        published = ReadWords(device, owner_handle, owner_bytes);
        EXPECT_EQ(published[0], 1U); EXPECT_EQ(published[1], capacity);
        EXPECT_EQ(published[4], 1U); EXPECT_EQ(published[5], 1U);
        for (std::size_t row = 0; row < rows; ++row) {
            const auto start = retained_owner::kPrefixWords + 1 + row * retained_owner::kRowWords;
            EXPECT_EQ(published[start], row == 4 ? 0U : 1U);
            if (row == 4) continue;
            EXPECT_EQ(published[start + 2], row < 2 ? 0U : 1U);
            EXPECT_EQ(published[start + 3], published[2]);
            EXPECT_EQ(published[start + 4], published[3]);
            const auto count = row < 2 ? retained_owner::kGeneralValues : retained_owner::kSchwarzschildValues;
            for (std::size_t slot = 0; slot < count; ++slot)
                EXPECT_EQ(published[start + retained_owner::kValuesOffset + 5 * slot + 4], 1U);
        }
        std::array<StageSevenToken, rows> tokens;
        for (std::size_t row = 0; row < rows; ++row) tokens[row] = (*candidate)->LastStepOwner(row);
        return tokens;
    };
    auto tokens = fresh();
    const auto t_ref = ReadAuthority("transport-program.bin"), e_ref = ReadAuthority("endpoint-program.bin");
    ASSERT_TRUE(std::equal(t_ref.begin(), t_ref.end(), published.begin() + retained_owner::kTransportReferenceOffset + 1));
    ASSERT_TRUE(std::equal(e_ref.begin(), e_ref.end(), published.begin() + retained_owner::kEndpointReferenceOffset + 1));
    const auto check_map = [&](auto map, std::size_t offset) {
        EXPECT_TRUE(std::equal(map.begin(), map.end(), published.begin() + offset + 1));
    };
    check_map(retained_owner::kGeneralTransportMap, retained_owner::kGeneralTransportMapOffset);
    check_map(retained_owner::kGeneralEndpointMap, retained_owner::kGeneralEndpointMapOffset);
    check_map(retained_owner::kSchwarzschildTransportMap, retained_owner::kSchwarzschildTransportMapOffset);
    check_map(retained_owner::kSchwarzschildEndpointMap, retained_owner::kSchwarzschildEndpointMapOffset);
    Save("published-owner.bin", std::as_bytes(std::span(published)));
    std::size_t calls = 0;
    const auto compare = [&](std::string label, std::span<const RetainedEndpointInput> inputs,
                             std::span<const StageSevenToken> grants, std::size_t expected_bindings) {
        ordinary.capture_output = ordinary.allocations[5].handle; ordinary.readbacks.clear();
        owned.capture_output = owned.allocations[5].handle; owned.readbacks.clear();
        std::cerr << "[OwnerStream] " << label << " begin" << std::endl;
        auto baseline = (*stock)->Endpoint(inputs);
        auto actual = (*candidate)->Endpoint(inputs, nullptr, grants);
        EXPECT_TRUE(baseline); EXPECT_TRUE(actual);
        EXPECT_EQ(endpoint_bindings, expected_bindings);
        EXPECT_EQ(ordinary.readbacks.size(), 1U); EXPECT_EQ(owned.readbacks.size(), 1U);
        if (!baseline || !actual || ordinary.readbacks.empty() || owned.readbacks.empty())
            throw std::runtime_error("endpoint comparison did not complete");
        EXPECT_EQ(ordinary.readbacks.back(), owned.readbacks.back()) << label;
        const auto& raw = owned.readbacks.back();
        EXPECT_EQ(raw.size(), inputs.size() * retained_program::kEndpointRowWords * 4);
        RecordProperty("readback_" + label, Hash(raw));
        Save(label + "-stock.bin", ordinary.readbacks.back()); Save(label + "-owner.bin", raw);
        ++calls;
        std::cerr << "[OwnerStream] " << label << " completed" << std::endl;
        return *actual;
    };
    const auto high = compare("upper-hit", upper, tokens, 3);
    for (std::size_t row = 0; row < rows; ++row) ASSERT_TRUE(high[row].valid);
    const auto after = ReadWords(device, owner_handle, owner_bytes);
    EXPECT_EQ(after[16], 15U);
    EXPECT_TRUE(std::equal(published.begin() + 17, published.end(), after.begin() + 17));
    const auto e_words = ReadWords(device, owned.allocations[5].handle, owned.allocations[5].bytes);
    for (std::size_t row = 0; row < 4; ++row)
        for (std::size_t key = 0; key < 9; ++key) {
            const auto index = key < 8 ? key : 44;
            const auto owner_key = retained_owner::kPrefixWords + 1 + row * retained_owner::kRowWords + retained_owner::kKeyOffset + key * 5;
            const auto endpoint_key = row * retained_program::kEndpointRowWords + 504 + index * 5;
            EXPECT_TRUE(std::equal(published.begin() + owner_key, published.begin() + owner_key + 5,
                                   e_words.begin() + endpoint_key));
        }
    compare("consumed-miss", upper, tokens, 2);
    tokens = fresh(); compare("lower-token-free", lower, {}, 2);
    ASSERT_NE(std::memcmp(upper[0].values.data(), lower[0].values.data(), sizeof(upper[0].values)), 0);
    tokens = fresh();
    std::array<StageSevenToken, rows> lower_grant{tokens[0]};
    compare("lower-context-miss", lower, lower_grant, 2);
    tokens = fresh(); auto changed = upper;
    changed[0].values[8] = RetainedValue::FromDouble(-changed[0].values[8].Center() + 1);
    std::array<StageSevenToken, rows> single{tokens[0]};
    compare("momentum-miss", changed, single, 2);
    tokens = fresh(); auto reordered = upper; std::swap(reordered[0], reordered[1]);
    single = {}; single[0] = tokens[1]; compare("row-miss", reordered, single, 2);
    const auto stale = fresh(); tokens = fresh(); compare("stale-miss", upper, stale, 2);
    tokens = fresh(); single = {}; single[0] = tokens[0];
    owned.fail_next = TransferProbeDevice::Failure::Write;
    EXPECT_FALSE((*candidate)->Endpoint(upper, nullptr, single));
    compare("failed-grant-consumed", upper, single, 2);
    tokens = fresh(); compare("fresh-recovery", upper, tokens, 3);
    // The new host epoch protocol retires reuse after any external owner write.
    // The legacy checked-shader invalid/valid-record sensitivity control remains
    // a separate outstanding acceptance obligation, not a weakened assertion.
    tokens = fresh(); single = {}; single[0] = tokens[0];
    auto invalid = published;
    const auto value_address = retained_owner::kPrefixWords + 1 + retained_owner::kValuesOffset + 5 * kGeneralSensitivitySlot;
    invalid[value_address + 4] = 0;
    ASSERT_TRUE(device.WriteBuffer(owner_handle, std::as_bytes(std::span(invalid))));
    compare("external-owner-write-stock-fallback", upper, single, 2);
    ASSERT_FALSE((*candidate)->StageSevenOwnershipEnabled());
    RecordProperty("equal_endpoint_calls", std::to_string(calls));
    RecordProperty("independent_step_fields_per_fresh_generation", "1360");
    RecordProperty("scope", "Five original reference rows and actual sealed hit/miss/key/record/full-buffer checks, device-failure recovery and external-owner retirement only. Original full nonvacuous sensitivity witness, eight controls, all performance gates and full qualification remain outstanding.");
}
}  // namespace
int main(int argc, char** argv) {
    if (argc < 2) return 2;
    artifacts = argv[1];
    for (int i = 1; i + 1 < argc; ++i) argv[i] = argv[i + 1];
    --argc;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
