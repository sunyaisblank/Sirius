#include "sirius/backend/retained_compute.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstring>
#include <stdexcept>

using namespace sirius::backend;
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
// A wire/failure model only: every submitted Camera row is declined, not solved.
class WireDevice final : public ComputeDevice {
  public:
    struct Write { BufferHandle buffer; std::vector<std::byte> bytes; };
    DeviceInfo info{.kind = DeviceKind::kSoftware};
    std::vector<std::vector<std::byte>> buffers;
    std::vector<Write> writes;
    std::uint32_t loads = 0, dispatches = 0, reads = 0, expected_groups = 1;
    bool fail_first = true;
    const DeviceInfo& Info() const noexcept override { return info; }
    sirius::base::Expected<KernelHandle> LoadKernel(std::span<const std::uint32_t>) override {
        return KernelHandle{loads++};
    }
    sirius::base::Expected<BufferHandle> CreateBuffer(std::uint64_t bytes, BufferUsage usage) override {
        require(usage == BufferUsage::kStorage, "unexpected allocation usage");
        BufferHandle handle{static_cast<std::uint32_t>(buffers.size())};
        buffers.emplace_back(bytes, std::byte{0xa5}); return handle;
    }
    sirius::base::Expected<std::uint64_t> RequiredBufferAllocationBytes(
        std::uint64_t bytes, BufferUsage) override { return bytes; }
    sirius::base::Expected<void> WriteBuffer(BufferHandle buffer, std::span<const std::byte> bytes) override {
        require(buffer.value < buffers.size() && bytes.size() <= buffers[buffer.value].size(), "write range");
        writes.push_back({buffer, {bytes.begin(), bytes.end()}});
        if (fail_first) {
            fail_first = false;
            return sirius::base::Fail(sirius::base::ErrorDomain::kDevice,
                                      "upload boundary write", "injected first returned failure");
        }
        std::copy(bytes.begin(), bytes.end(), buffers[buffer.value].begin()); return {};
    }
    sirius::base::Expected<void> ReadBuffer(BufferHandle buffer, std::span<std::byte> bytes) override {
        require(buffer.value < buffers.size() && bytes.size() <= buffers[buffer.value].size(), "read range");
        ++reads; std::copy_n(buffers[buffer.value].begin(), bytes.size(), bytes.begin()); return {};
    }
    sirius::base::Expected<void> Dispatch(KernelHandle kernel, std::span<const BufferHandle> bindings,
                                         std::uint32_t x, std::uint32_t y, std::uint32_t z,
                                         DispatchTiming* timing) override {
        require(kernel.value == 0 && bindings.size() == 2 && x == expected_groups && y == 1 && z == 1,
                "unexpected Camera dispatch");
        require(bindings[0].value < buffers.size() && bindings[1].value < buffers.size(), "dispatch buffer range");
        ++dispatches;
        std::fill(buffers[bindings[1].value].begin(), buffers[bindings[1].value].end(), std::byte{0});
        if (timing) *timing = {};
        return {};
    }
    sirius::base::Expected<void> SetBufferAllocationLimit(std::uint64_t) override { return {}; }
    std::uint64_t BufferAllocationBytes() const noexcept override {
        std::uint64_t bytes = 0; for (const auto& buffer : buffers) bytes += buffer.size(); return bytes;
    }
};

RetainedCameraInput row(unsigned tag) {
    RetainedCameraInput result;
    for (std::size_t i = 0; i < result.values.size(); ++i)
        result.values[i] = {static_cast<float>(tag * 100 + i + 1), 0, 0, 0, 1};
    return result;
}
int main() {
    try {
        WireDevice device;
        auto created = RetainedCompute::Create(device, 3);
        require(bool(created), "Create failed");
        auto& compute = **created;
        require(device.loads == RetainedCompute::kStageCount && device.writes.empty(), "Create performed upload");
        const auto allocation = device.BufferAllocationBytes();
        const auto allocation_count = device.buffers.size();
        std::array initial{row(1), row(2), row(3)};
        auto failed = compute.Camera(initial);
        require(!failed && failed.error().operation() == "upload boundary write", "first failure not returned");
        require(device.writes.size() == 1 && !device.dispatches && !device.reads, "failed write submitted work");
        const auto handle = device.writes.front().buffer.value;
        const auto full_bytes = device.buffers[handle].size();
        const auto program_start = sizeof(std::uint32_t) + 3 * sizeof(RetainedCameraInput);
        const auto& attempted = device.writes.front().bytes;
        require(attempted.size() == full_bytes && program_start < full_bytes, "first upload not full");
        std::vector<std::byte> program(attempted.begin() + program_start, attempted.end());
        require(std::ranges::any_of(program, [](auto byte) { return byte != std::byte{0}; }), "empty program oracle");
        require(std::ranges::all_of(device.buffers[handle], [](auto byte) { return byte == std::byte{0xa5}; }),
                "failed fake write mutated storage");
        const auto check = [&](std::span<const RetainedCameraInput> rows) {
            const auto& wire = device.buffers[handle];
            std::uint32_t capacity = 0; std::memcpy(&capacity, wire.data(), sizeof(capacity));
            require(capacity == 3, "capacity header changed");
            const auto bytes = std::as_bytes(rows);
            require(std::equal(bytes.begin(), bytes.end(), wire.begin() + sizeof(capacity)), "active bytes differ");
            require(std::equal(program.begin(), program.end(), wire.begin() + program_start), "program tail changed");
        };
        std::array retry{row(4)};
        auto result = compute.Camera(retry);
        require(bool(result) && result->size() == 1 && !result->front().valid, "retry decline protocol");
        require(device.writes.size() == 2 && device.writes.back().bytes.size() == full_bytes, "retry not full");
        check(retry);
        const auto inactive = sizeof(std::uint32_t) + sizeof(RetainedCameraInput);
        require(std::all_of(device.buffers[handle].begin() + inactive,
                            device.buffers[handle].begin() + program_start,
                            [](auto byte) { return byte == std::byte{0}; }), "failed larger batch leaked into first upload");
        std::array grown{row(5), row(6), row(7)};
        device.expected_groups = 3; result = compute.Camera(grown);
        require(bool(result) && result->size() == 3, "grow failed");
        require(device.writes.size() == 3 && device.writes.back().bytes.size() == program_start, "grow not prefix");
        check(grown);
        std::vector<std::byte> previous_inactive(device.buffers[handle].begin() + inactive,
                                                device.buffers[handle].begin() + program_start);
        std::array shrunk{row(8)};
        device.expected_groups = 1; result = compute.Camera(shrunk);
        require(bool(result) && result->size() == 1, "shrink failed");
        require(device.writes.size() == 4 && device.writes.back().bytes.size() == inactive, "shrink not prefix");
        check(shrunk);
        require(std::equal(previous_inactive.begin(), previous_inactive.end(),
                           device.buffers[handle].begin() + inactive), "prefix write changed inactive device rows");
        const auto stats = compute.Statistics()[0];
        require(device.dispatches == 3 && device.reads == 3 && stats.submissions == 3 &&
                stats.command_row_counts[1] == 2 && stats.command_row_counts[3] == 1, "submission accounting changed");
        require(device.BufferAllocationBytes() == allocation && device.buffers.size() == allocation_count,
                "device allocation changed");
        std::printf("{\"passed\":true,\"write_attempts\":4,\"successful_dispatches\":3,"
                    "\"successful_reads\":3,\"first_failed_rows\":3,\"retry_rows\":1,"
                    "\"grow_rows\":3,\"shrink_rows\":1,\"program_bytes\":%zu,"
                    "\"first_success_zero_inactive\":true,\"later_prefix_preserved_inactive\":true}\n", program.size());
        return std::ferror(stdout) ? 1 : 0;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "upload boundary: %s\n", error.what()); return 1;
    }
}
