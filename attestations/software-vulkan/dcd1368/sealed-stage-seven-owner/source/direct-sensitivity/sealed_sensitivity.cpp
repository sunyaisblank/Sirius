#include "sirius/backend/retained_compute.h"
#include "retained_owner_layout.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <bit>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <span>
#include <stdexcept>
#include <vector>
namespace {
using namespace sirius::backend;
std::filesystem::path artifacts;
std::vector<std::uint32_t> Read(const std::filesystem::path& path) {
    const auto bytes = std::filesystem::file_size(path);
    if (bytes % 4) throw std::runtime_error("unaligned fixture");
    std::vector<std::uint32_t> words(bytes / 4);
    std::ifstream file(path, std::ios::binary);
    file.read(reinterpret_cast<char*>(words.data()), bytes);
    if (!file) throw std::runtime_error("fixture read failed");
    return words;
}
TEST(SealedShaderSensitivity, CanonicalFallbackAndLiveRecordsRemainConsequential) {
    constexpr std::size_t capacity = 24, rows = 5, row_words = 3829;
    constexpr std::uint64_t limit = 8 * 1024 * 1024;
    const auto fixtures = std::filesystem::path(__FILE__).parent_path();
    const auto inventory = EnumerateVulkanDevices();
    ASSERT_TRUE(inventory); ASSERT_EQ(inventory->size(), 1U);
    auto opened = CreateVulkanDevice(0);
    ASSERT_TRUE(opened);
    auto& device = **opened;
    ASSERT_EQ(device.Info().kind, DeviceKind::kSoftware);
    ASSERT_TRUE(device.Info().rounds_fp32_to_nearest);
    ASSERT_TRUE(RetainedUsesPortableArithmetic(device.Info()));
    ASSERT_TRUE(device.SetBufferAllocationLimit(limit));
    const auto program = Read(fixtures / "endpoint-program.bin");
    const auto expected = Read(fixtures / "upper-hit-stock.bin");
    auto published = Read(fixtures / "published-owner.bin");
    ASSERT_EQ(expected.size(), rows * row_words);
    ASSERT_EQ(published.size(), retained_owner::kPrefixWords + 1 + capacity * retained_owner::kRowWords);
    ASSERT_TRUE(std::equal(program.begin(), program.end(), published.begin() + retained_owner::kEndpointReferenceOffset + 1));
    ASSERT_TRUE(std::equal(retained_owner::kGeneralEndpointMap.begin(), retained_owner::kGeneralEndpointMap.end(), published.begin() + retained_owner::kGeneralEndpointMapOffset + 1));
    ASSERT_TRUE(std::equal(retained_owner::kSchwarzschildEndpointMap.begin(), retained_owner::kSchwarzschildEndpointMap.end(), published.begin() + retained_owner::kSchwarzschildEndpointMapOffset + 1));
    std::vector<std::uint32_t> input(1 + capacity * 225 + program.size());
    input[0] = capacity;
    for (std::size_t row = 0; row < rows; ++row)
        std::copy_n(expected.begin() + row * row_words + 504, 225, input.begin() + 1 + row * 225);
    std::copy(program.begin(), program.end(), input.begin() + 1 + capacity * 225);
    const auto shader = Read("bin/linux-gcc/sealed-stage-seven-owner/products/retained_endpoint_portable_stage7_sealed.spv");
    auto kernel = device.LoadKernel(shader);
    ASSERT_TRUE(kernel);
    auto in = device.CreateBuffer(input.size() * 4, BufferUsage::kStorage);
    auto out = device.CreateBuffer(capacity * row_words * 4, BufferUsage::kStorage);
    auto sidecar = device.CreateBuffer(published.size() * 4, BufferUsage::kStorage);
    ASSERT_TRUE(in); ASSERT_TRUE(out); ASSERT_TRUE(sidecar);
    ASSERT_TRUE(device.WriteBuffer(*in, std::as_bytes(std::span(input))));
    std::vector<std::uint32_t> actual_input(input.size());
    ASSERT_TRUE(device.ReadBuffer(*in, std::as_writable_bytes(std::span(actual_input))));
    ASSERT_EQ(actual_input, input);
    const std::array bindings{*in, *out, *sidecar};
    const auto execute = [&](std::string name, const std::vector<std::uint32_t>& owner) {
        EXPECT_TRUE(device.WriteBuffer(*sidecar, std::as_bytes(std::span(owner))));
        std::vector<std::uint32_t> verified(owner.size());
        EXPECT_TRUE(device.ReadBuffer(*sidecar, std::as_writable_bytes(std::span(verified))));
        EXPECT_EQ(verified, owner);
        const auto submitted = device.Dispatch(*kernel, bindings, rows, 1, 1);
        EXPECT_TRUE(submitted);
        std::vector<std::uint32_t> result(rows * row_words);
        EXPECT_TRUE(device.ReadBuffer(*out, std::as_writable_bytes(std::span(result))));
        std::ofstream file(artifacts / (name + ".bin"), std::ios::binary);
        file.write(reinterpret_cast<const char*>(result.data()), result.size() * 4);
        EXPECT_TRUE(file);
        return result;
    };
    published[16] = 15;
    EXPECT_EQ(execute("direct-positive", published), expected);
    auto invalid = published;
    const auto first = retained_owner::kPrefixWords + 1 + retained_owner::kValuesOffset;
    invalid[first + 5 * 70 + 4] = 0;
    EXPECT_EQ(execute("direct-invalid-record-fallback", invalid), expected);
    for (const auto row : {std::size_t{0}, std::size_t{2}}) {
        auto altered = published;
        altered[16] = 1U << row;
        const auto address = retained_owner::kPrefixWords + 1 + row * retained_owner::kRowWords + retained_owner::kValuesOffset + 5 * (row == 0 ? 70 : 25);
        const auto previous = std::bit_cast<RetainedValue>(std::array<std::uint32_t, 5>{altered[address], altered[address+1], altered[address+2], altered[address+3], altered[address+4]});
        ASSERT_TRUE(previous.IsRepresented());
        const auto replacement = RetainedValue::FromDouble(previous.Center() == 0 ? 1 : -previous.Center());
        ASSERT_TRUE(replacement.IsRepresented());
        const auto words = std::bit_cast<std::array<std::uint32_t, 5>>(replacement);
        ASSERT_FALSE(std::equal(words.begin(), words.end(), altered.begin() + address));
        std::copy(words.begin(), words.end(), altered.begin() + address);
        const auto result = execute("direct-sensitivity-" + std::to_string(row), altered);
        const auto start = row * row_words;
        ASSERT_EQ(result[start], 1U);
        std::size_t changed = 0;
        EXPECT_TRUE(std::equal(result.begin() + start, result.begin() + start + 104,
                               expected.begin() + start));
        // Input echo, projection packet and public padding are unchanged.
        // This module keeps live register scratch in group-shared storage.
        EXPECT_TRUE(std::equal(result.begin() + start + 504,
                               result.begin() + start + row_words,
                               expected.begin() + start + 504));
        // These are the actual decoded phase/physical fields, not scratch words.
        for (std::size_t i = 104; i < 504; ++i) changed += result[start + i] != expected[start + i];
        EXPECT_GT(changed, 0U);
        RecordProperty("changed_decoded_words_" + std::to_string(row), std::to_string(changed));
        for (std::size_t other = 0; other < rows; ++other) {
            if (other == row) continue;
            EXPECT_TRUE(std::equal(result.begin() + other * row_words, result.begin() + (other + 1) * row_words, expected.begin() + other * row_words));
        }
    }
    EXPECT_EQ(execute("direct-post-sensitivity-recovery", published), expected);
    RecordProperty("scope", "Direct actual sealed shader seam: five original actual source-bound input rows; complete positive/fallback/recovery and valid-one General/Schwarzschild live-record sensitivity. Host API deliberately bypassed; it separately retires external owner edits. No benefit/frame/full qualification claim.");
}
}
int main(int argc, char** argv) {
    if (argc < 2) return 2;
    artifacts = argv[1];
    for (int i = 1; i + 1 < argc; ++i) argv[i] = argv[i + 1];
    --argc;
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
