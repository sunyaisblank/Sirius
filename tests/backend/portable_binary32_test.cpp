#include "support/test_resource.h"

#include "sirius/kernels/portable_binary32.h"

#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <string>

#ifndef SIRIUS_TEST_HAS_PORTABLE_BINARY32_REFERENCE
#error "The build must provide the generated independent portable-binary32 fixture"
#endif

namespace {
using namespace sirius::portable_binary32;

constexpr std::uint32_t kMagic = 0x32334250u;  // Little-endian bytes "PB32".
constexpr std::uint32_t kVersion = 1u;
constexpr std::uint32_t kCorpusRecords = 120744u + 53524u;
constexpr std::uint32_t kProductResidualOperation = 13u;
constexpr std::uint32_t kMaximumMismatches = 8u;

template <std::size_t N>
std::uint32_t Word(const std::array<std::uint8_t, N>& bytes, std::size_t offset) {
    return std::uint32_t(bytes[offset]) | (std::uint32_t(bytes[offset + 1]) << 8) |
           (std::uint32_t(bytes[offset + 2]) << 16) | (std::uint32_t(bytes[offset + 3]) << 24);
}

PB32Product Evaluate(std::uint32_t operation, std::uint32_t a, std::uint32_t b) {
    PB32Product result{0u, 0u, 1u};
    switch (operation) {
        case 0u:
            result.high = PB32Add(a, b);
            break;
        case 1u:
            result.high = PB32Subtract(a, b);
            break;
        case 2u:
            result.high = PB32Multiply(a, b);
            break;
        case 3u:
            result.high = PB32Divide(a, b);
            break;
        case 4u:
            result.high = PB32Sqrt(a);
            break;
        case 5u:
            result.high = PB32Equal(a, b) ? 1u : 0u;
            break;
        case 6u:
            result.high = PB32Less(a, b) ? 1u : 0u;
            break;
        case 7u:
            result.high = PB32LessEqual(a, b) ? 1u : 0u;
            break;
        case 8u:
            result.high = PB32FromUnsigned(a);
            break;
        case 9u:
            result.high = PB32FromSignedBits(a);
            break;
        case 10u:
            result.high = PB32Abs(a);
            break;
        case 11u:
            result.high = PB32Negate(a);
            break;
        case 12u:
            result.high = PB32Classify(a);
            break;
        case kProductResidualOperation:
            result = PB32ProductResidual(a, b);
            break;
        default:
            result.valid = 0u;
            break;  // Rejected by the fixture parser first.
    }
    return result;
}
}  // namespace

TEST(PortableBinary32, MatchesIndependentExactCorpus) {
    const std::string path = sirius::test::ResourcePath(
        "tests/backend/portable_binary32_reference.bin");
    std::ifstream stream(path, std::ios::binary);
    ASSERT_TRUE(stream.is_open()) << "Missing independent fixture: " << path;
    std::array<std::uint8_t, 12> header{};
    stream.read(reinterpret_cast<char*>(header.data()),
                static_cast<std::streamsize>(header.size()));
    ASSERT_EQ(stream.gcount(), static_cast<std::streamsize>(header.size()))
        << "Empty or truncated independent fixture header";
    ASSERT_EQ(Word(header, 0), kMagic) << "Unsupported fixture magic";
    ASSERT_EQ(Word(header, 4), kVersion) << "Unsupported fixture version";
    const auto count = Word(header, 8);
    ASSERT_EQ(count, kCorpusRecords) << "Empty, incomplete or excessive version-1 corpus";

    std::uint32_t mismatches = 0u;
    std::uint32_t checked = 0u;
    for (std::uint32_t row = 0u; row < count; ++row) {
        std::array<std::uint8_t, 24> record{};
        stream.read(reinterpret_cast<char*>(record.data()),
                    static_cast<std::streamsize>(record.size()));
        ASSERT_EQ(stream.gcount(), static_cast<std::streamsize>(record.size()))
            << "Truncated independent record at row " << row;
        const auto operation = Word(record, 0);
        const auto a = Word(record, 4);
        const auto b = Word(record, 8);
        const auto expectedHigh = Word(record, 12);
        const auto expectedLow = Word(record, 16);
        const auto expectedValid = Word(record, 20);
        ASSERT_LE(operation, kProductResidualOperation) << "Unsupported operation at row " << row;
        ASSERT_LE(expectedValid, 1u) << "Unsupported validity field at row " << row;
        if (operation != kProductResidualOperation) {
            ASSERT_EQ(expectedLow, 0u) << "Unsupported scalar record at row " << row;
            ASSERT_EQ(expectedValid, 1u) << "Unsupported scalar validity at row " << row;
        }
        const auto actual = Evaluate(operation, a, b);
        ++checked;
        if (actual.high != expectedHigh || actual.low != expectedLow ||
            actual.valid != expectedValid) {
            ++mismatches;
            ADD_FAILURE() << "Independent exact corpus row " << row << ", operation " << operation
                          << std::hex << std::setfill('0') << ", a=0x" << std::setw(8) << a
                          << ", b=0x" << std::setw(8) << b << ", expected=(0x" << std::setw(8)
                          << expectedHigh << ",0x" << std::setw(8) << expectedLow << ",0x"
                          << std::setw(8) << expectedValid << "), actual=(0x" << std::setw(8)
                          << actual.high << ",0x" << std::setw(8) << actual.low << ",0x"
                          << std::setw(8) << actual.valid << ')';
            if (mismatches == kMaximumMismatches)
                FAIL() << "Stopped after eight independent mismatches";
        }
    }
    char trailing = 0;
    stream.read(&trailing, 1);
    EXPECT_EQ(stream.gcount(), 0) << "Trailing bytes after complete independent corpus";
    EXPECT_TRUE(stream.eof()) << "Independent fixture did not end cleanly";
    EXPECT_EQ(checked, kCorpusRecords);
    EXPECT_EQ(mismatches, 0u);
    RecordProperty("independent_primitive_cases", "120744");
    RecordProperty("independent_product_residual_cases", "53524");
    RecordProperty("independent_exact_records_checked", std::to_string(checked));
    RecordProperty("exact_raw_word_mismatches", std::to_string(mismatches));
}
