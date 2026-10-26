// Focused routing compatibility, not an independent numerical oracle.
// Compile with the original generic Slang typed-batch directory first among
// system include paths, /opt/slang/include, and the original strict FP flags.
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <immintrin.h>
#include <endpoint.program.hpp>

// The namespace and three extern-C entrypoints must all differ from the actual
// scalar worker's generated body. No arithmetic or selector is substituted.
#define sirius_generated_endpoint sirius_endpoint_header_reference
#define sirius_endpointComputeMain_Thread sirius_endpoint_header_reference_Thread
#define sirius_endpointComputeMain_Group sirius_endpoint_header_reference_Group
#define sirius_endpointComputeMain sirius_endpoint_header_reference_Main
#include <endpoint.native.hpp>
#undef sirius_endpointComputeMain
#undef sirius_endpointComputeMain_Group
#undef sirius_endpointComputeMain_Thread
#undef sirius_generated_endpoint

#include "../../tests/support/retained_transport/endpoint_factor_reference.h"
#include "../../tests/support/retained_transport/endpoint_reference.h"

namespace reference = sirius_endpoint_header_reference;
namespace fixture = sirius::test::retained_endpoint_factor;
using InputRow = std::array<std::uint32_t, 225>;
using Value = std::array<std::uint32_t, 5>;
constexpr std::uint32_t kCapacity = 32;
constexpr std::uint32_t kActive = 23;
constexpr std::uint32_t kOutputWords = 3829;
constexpr std::uint32_t kProgramWords = 21143;
constexpr std::uint32_t kProgramBase = 1 + kCapacity * 225;
constexpr std::uint32_t kHeaderWord = 11241;
constexpr std::size_t kCanaryWords = 16;
static_assert(sizeof(reference::RetainedTriple_0) == 20);
static_assert(std::endian::native == std::endian::little);
static_assert(std::size(sirius_cpu_bridge::endpoint_program) == kProgramWords);
static_assert(sirius_cpu_bridge::endpoint_program[1] == 612);
static_assert(sirius_cpu_bridge::endpoint_program[11240] == 1943);
static_assert(sirius_cpu_bridge::endpoint_program[kHeaderWord] == 612);

#define PACKET_TYPES const void*, std::uint32_t, void*, std::uint32_t, \
    std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t, \
    std::uint32_t, std::uint32_t, std::uint32_t
using Packet = bool (*)(PACKET_TYPES) noexcept;
extern "C" bool sirius_native_endpoint_try_packet(PACKET_TYPES) noexcept;
extern "C" bool sirius_native_endpoint_admit_packet(PACKET_TYPES) noexcept;
#undef PACKET_TYPES
extern "C" void sirius_simd_endpoint_packet(std::uint32_t*, std::uint32_t,
    std::uint32_t*, std::uint32_t, std::uint32_t, std::uint32_t);
// Exact original rows_adapter.cpp bridge body: the raw ISPC program calls the
// same scalar libm fmaxf as its original generated C++ reference.
extern "C" float sirius_simd_fmaxf(float a, float b) noexcept { return ::fmaxf(a, b); }

static void require(bool predicate, const char* reason) {
    if (!predicate) throw std::runtime_error(reason);
}

struct Guarded {
    std::size_t count;
    std::vector<std::uint32_t> storage;
    explicit Guarded(std::size_t words) : count(words),
        storage(words + 2 * kCanaryWords, 0xa5a5a5a5U) {
        for (std::size_t i = 0; i < kCanaryWords; ++i) {
            storage[i] = 0x13570000U + static_cast<std::uint32_t>(i);
            storage[kCanaryWords + count + i] = 0x24680000U + static_cast<std::uint32_t>(i);
        }
    }
    std::uint32_t* data() { return storage.data() + kCanaryWords; }
    const std::uint32_t* data() const { return storage.data() + kCanaryWords; }
    std::uint32_t bytes() const { return static_cast<std::uint32_t>(count * 4); }
};

struct Arguments {
    const void* input;
    std::uint32_t input_bytes;
    void* output;
    std::uint32_t output_bytes;
    std::uint32_t first = 0, y = 0, z = 0;
    std::uint32_t grid = kActive, grid_y = 1, grid_z = 1, groups = 1;
};

static Arguments arguments(Guarded& input, Guarded& output) {
    return {input.data(), input.bytes(), output.data(), output.bytes()};
}

static bool invoke(Packet function, const Arguments& a) {
    return function(a.input, a.input_bytes, a.output, a.output_bytes,
        a.first, a.y, a.z, a.grid, a.grid_y, a.grid_z, a.groups);
}

struct FloatingPointScope {
    unsigned saved = _mm_getcsr();
    unsigned required = (saved | 0x1f80U) & ~0xe040U;
    FloatingPointScope() { _mm_setcsr(required); }
    ~FloatingPointScope() { _mm_setcsr(saved); }
    void check() const {
        // Sticky exception flags are not the original controls contract.
        require((_mm_getcsr() & ~0x3fU) == (required & ~0x3fU), "FP controls changed");
    }
};

struct Counts {
    std::uint64_t reference_rows = 0, scalar_packets = 0, raw_packets = 0;
    std::uint64_t pure_calls = 0, other_program_words = 0, boundary_cases = 0;
    std::uint64_t compared_output_words = 0, compared_input_words = 0;
};

static void same_output(const Guarded& actual, const Guarded& expected, Counts& counts) {
    counts.compared_output_words += actual.storage.size();
    require(actual.storage == expected.storage, "Full output/inactive/canary mismatch");
}

static void same_input(const Guarded& input, const std::vector<std::uint32_t>& expected,
                       Counts& counts) {
    counts.compared_input_words += input.storage.size();
    require(input.storage == expected, "Input or input canary changed");
}

static void set_value(InputRow& row, std::size_t field, const Value& value) {
    std::copy(value.begin(), value.end(), row.begin() + field * 5);
}

static std::vector<InputRow> selection() {
    // The original source builder is retained_compute_test.cpp:2737--2778.
    std::vector<InputRow> rows;
    for (std::size_t i = 20; i < 26; ++i) rows.push_back(fixture::cases[i].input);
    auto row = fixture::cases[20].input;
    set_value(row, 1, {0x3f000000U, 0xbf000000U, 0x80000000U, 0, 1});
    rows.push_back(row);
    rows.push_back(fixture::cases[0].input);
    require(std::strcmp(sirius::test::retained_endpoint::cases[14].name, "analytic-flat") == 0,
            "Original flat fixture identity changed");
    rows.push_back(sirius::test::retained_endpoint::cases[14].input);
    for (std::size_t parameter = 1; parameter <= 3; ++parameter) {
        for (unsigned kind = 0; kind < 3; ++kind) {
            row = fixture::cases[20].input;
            if (kind == 0) set_value(row, parameter, {47U << 23, 0, 0, 0, 1});
            if (kind == 1) row[parameter * 5 + 1] = 37U << 23;
            if (kind == 2) row[parameter * 5 + 3] = 37U << 23;
            rows.push_back(row);
        }
    }
    for (std::size_t parameter = 1; parameter <= 3; ++parameter) {
        row = fixture::cases[20].input;
        row[parameter * 5 + 3] = 0x80000000U;
        rows.push_back(row);
    }
    for (const bool uncertain : {false, true}) {
        row = fixture::cases[20].input;
        for (std::size_t field = 5; field < 8; ++field) set_value(row, field, {0, 0, 0, 0, 1});
        if (uncertain) set_value(row, 5, {0x3e800000U, 0, 0, 0x3f000000U, 1});
        rows.push_back(row);
    }
    require(rows.size() == kActive, "Original selection must contain23 rows");
    return rows;
}

static Guarded packet(const std::vector<InputRow>& rows, std::uint32_t header) {
    Guarded input(kProgramBase + kProgramWords);
    input.data()[0] = kCapacity;
    for (std::size_t i = 0; i < rows.size(); ++i)
        std::copy(rows[i].begin(), rows[i].end(), input.data() + 1 + i * 225);
    std::copy_n(sirius_cpu_bridge::endpoint_program, kProgramWords, input.data() + kProgramBase);
    input.data()[kProgramBase + kHeaderWord] = header;
    return input;
}

static Guarded original_reference(Guarded& input, const FloatingPointScope& fp, Counts& counts) {
    Guarded output(kCapacity * kOutputWords);
    reference::GlobalParams_0 params{};
    params.inputs_0.data = input.data();
    params.inputs_0.count = input.count;
    params.outputs_0.data = output.data();
    params.outputs_0.count = output.count;
    for (std::uint32_t row = 0; row < kActive; ++row) {
        ComputeVaryingInput varying{};
        varying.startGroupID.x = row;
        reference::sirius_endpoint_header_reference_Group(&varying, nullptr, &params);
        ++counts.reference_rows;
        fp.check();
    }
    return output;
}

static void routes(Guarded& input, const Guarded& expected,
                   const FloatingPointScope& fp, Counts& counts) {
    const auto original = input.storage;
    Guarded pure(kCapacity * kOutputWords), untouched = pure;
    auto a = arguments(input, pure);
    ++counts.pure_calls;
    require(invoke(sirius_native_endpoint_admit_packet, a), "Expanded pure admission declined");
    same_output(pure, untouched, counts);
    same_input(input, original, counts);

    Guarded scalar(kCapacity * kOutputWords);
    a = arguments(input, scalar);
    for (a.first = 0; a.first < kActive; ++a.first) {
        require(invoke(sirius_native_endpoint_try_packet, a), "Scalar worker declined admitted input");
        ++counts.scalar_packets;
        fp.check();
    }
    same_output(scalar, expected, counts);
    same_input(input, original, counts);

    Guarded mixed(kCapacity * kOutputWords);
    a = arguments(input, mixed);
    for (a.first = 0; a.first + 8 <= kActive; a.first += 8) {
        a.groups = 8;
        ++counts.pure_calls;
        require(invoke(sirius_native_endpoint_admit_packet, a), "Pure admission declined full raw packet");
        sirius_simd_endpoint_packet(input.data(), static_cast<std::uint32_t>(input.count),
            mixed.data(), static_cast<std::uint32_t>(mixed.count), a.first, 8);
        ++counts.raw_packets;
        fp.check();
    }
    a.groups = 1;
    for (; a.first < kActive; ++a.first) {
        require(invoke(sirius_native_endpoint_try_packet, a), "Scalar tail declined admitted input");
        ++counts.scalar_packets;
        fp.check();
    }
    same_output(mixed, expected, counts);
    same_input(input, original, counts);
}

int main() {
    if (!__builtin_cpu_supports("x86-64-v3")) return 2;
    Counts counts;
    const unsigned caller_state = _mm_getcsr();
    try {
        {
            FloatingPointScope fp;
            const auto rows = selection();
            constexpr std::array<std::uint32_t, 6> headers{612, 0, 1, 611, 613, 0xffffffffU};
            Guarded canonical(kCapacity * kOutputWords);
            for (const auto header : headers) {
                auto input = packet(rows, header);
                const auto before = input.storage;
                auto expected = original_reference(input, fp, counts);
                same_input(input, before, counts);
                if (header == 612) {
                    canonical = expected;
                    for (std::size_t row = 0; row < 18; ++row)
                        require(expected.data()[row * kOutputWords] == 1, "Original good selection row declined");
                } else {
                    for (std::size_t row = 0; row < kActive; ++row) {
                        if (row < 7 || row >= 21) {
                            const auto first = expected.data() + row * kOutputWords;
                            require(std::all_of(first, first + 504, [](auto word) { return word == 0; }),
                                    "Specialized-header refusal published public words");
                        } else if (row < 18) {
                            require(std::equal(expected.data() + row * kOutputWords,
                                expected.data() + (row + 1) * kOutputWords,
                                canonical.data() + row * kOutputWords), "General/flat row changed with specialized header");
                        }
                    }
                }
                routes(input, expected, fp, counts);
                // Original per-row permutation relation puts seven valid
                // specialized rows into scalar tails, retaining full raw packets.
                auto reversed = rows;
                std::reverse(reversed.begin(), reversed.end());
                auto reordered_input = packet(reversed, header);
                Guarded reordered_expected(kCapacity * kOutputWords);
                for (std::size_t row = 0; row < kActive; ++row)
                    std::copy_n(expected.data() + (kActive - 1 - row) * kOutputWords,
                                kOutputWords, reordered_expected.data() + row * kOutputWords);
                routes(reordered_input, reordered_expected, fp, counts);
            }

            auto input = packet(rows, 612);
            const auto initial = input.storage;
            Guarded output(kCapacity * kOutputWords), untouched = output;
            auto a = arguments(input, output);
            // Every other program word must retain full immutable authority.
            for (std::uint32_t word = 0; word < kProgramWords; ++word) {
                if (word == kHeaderWord) continue;
                input.data()[kProgramBase + word] ^= 1U;
                ++counts.pure_calls;
                require(!invoke(sirius_native_endpoint_admit_packet, a), "Unknown program passed pure admission");
                same_output(output, untouched, counts);
                require(!invoke(sirius_native_endpoint_try_packet, a), "Unknown program passed scalar admission");
                same_output(output, untouched, counts);
                input.data()[kProgramBase + word] ^= 1U;
                same_input(input, initial, counts);
                fp.check();
                ++counts.other_program_words;
            }
            std::vector<Arguments> invalid;
            auto add = [&](auto modify) { auto bad = a; modify(bad); invalid.push_back(bad); };
            add([](auto& bad) { bad.input = nullptr; });
            add([](auto& bad) { bad.output = nullptr; });
            add([](auto& bad) { bad.input_bytes = 3; });
            add([](auto& bad) { --bad.input_bytes; });
            add([](auto& bad) { bad.input_bytes += 4; });
            add([](auto& bad) { --bad.output_bytes; });
            add([](auto& bad) { bad.output_bytes += 4; });
            add([](auto& bad) { bad.input = static_cast<const std::byte*>(bad.input) + 1; });
            add([](auto& bad) { bad.output = static_cast<std::byte*>(bad.output) + 1; });
            add([](auto& bad) { bad.output = const_cast<void*>(bad.input); bad.output_bytes = 4; });
            add([](auto& bad) { bad.groups = 0; });
            add([](auto& bad) { bad.groups = 9; });
            add([](auto& bad) { bad.first = kActive; });
            add([](auto& bad) { bad.y = 1; });
            add([](auto& bad) { bad.z = 1; });
            add([](auto& bad) { bad.grid = 0; });
            add([](auto& bad) { bad.grid = kCapacity + 1; });
            add([](auto& bad) { bad.grid_y = 2; });
            add([](auto& bad) { bad.grid_z = 2; });
            for (const auto& bad : invalid) {
                ++counts.pure_calls;
                require(!invoke(sirius_native_endpoint_admit_packet, bad), "Original pure view bound weakened");
                require(!invoke(sirius_native_endpoint_try_packet, bad), "Original scalar view bound weakened");
                same_output(output, untouched, counts);
                same_input(input, initial, counts);
                fp.check();
                ++counts.boundary_cases;
            }
            for (const auto capacity : {0U, 65537U}) {
                input.data()[0] = capacity;
                ++counts.pure_calls;
                require(!invoke(sirius_native_endpoint_admit_packet, a), "Invalid capacity passed pure admission");
                require(!invoke(sirius_native_endpoint_try_packet, a), "Invalid capacity passed scalar admission");
                input.data()[0] = kCapacity;
                same_input(input, initial, counts);
                same_output(output, untouched, counts);
                fp.check();
                ++counts.boundary_cases;
            }
        }
        require(_mm_getcsr() == caller_state, "Caller FP state not restored");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "{\"passed\":true,\"selection_rows\":23,\"capacity\":32,\"header_values\":6,\"layouts\":2"
        << ",\"reference_rows\":" << counts.reference_rows
        << ",\"scalar_packets\":" << counts.scalar_packets
        << ",\"raw_packets\":" << counts.raw_packets
        << ",\"pure_calls\":" << counts.pure_calls
        << ",\"other_program_word_refusals\":" << counts.other_program_words
        << ",\"boundary_refusals\":" << counts.boundary_cases
        << ",\"compared_output_words\":" << counts.compared_output_words
        << ",\"compared_input_words\":" << counts.compared_input_words
        << ",\"whole3829_words_and_inactive_capacity_and_canaries\":true"
        << ",\"input_immutable\":true,\"MXCSR_controls_unchanged\":true,\"caller_MXCSR_state_restored\":true"
        << ",\"scope\":\"Finite selector/routing compatibility against unchanged generated scalar shader and checked-in original input builders; direct unchanged raw Endpoint plus scalar tails, not provider integration, independent numerical or performance/scientific/release qualification\"}\n";
}
