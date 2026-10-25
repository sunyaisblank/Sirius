#include "../../tests/support/retained_transport/reference_cases.h"
#include "../../tests/support/retained_transport/schwarzschild_reference_cases.h"
#include "../../bin/linux-gcc/native-cpu-slang-bridge/typed-batch/transport.program.hpp"
#include "../native-cpu-simd-packet-policy/rows_bridge.h"
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstring>
#include <immintrin.h>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" void sirius_simd_transport_packet(std::uint32_t*, std::uint32_t,
    std::uint32_t*, std::uint32_t, std::uint32_t, std::uint32_t);
extern "C" void baseline_sirius_simd_transport_packet(std::uint32_t*, std::uint32_t,
    std::uint32_t*, std::uint32_t, std::uint32_t, std::uint32_t);
extern "C" bool sirius_native_simd_supported() noexcept {
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("x86-64-v3");
}
using Words = std::vector<std::uint32_t>;
constexpr unsigned capacity = 40, active = 38, stride = 230, output_stride = 4699;
constexpr unsigned input_words = 25168, output_words = 187960, sentinel = 0xa5a5a5a5u;
static_assert(sirius::test::retained_transport::cases.size() == 15);
static_assert(sirius::test::retained_schwarzschild::cases.size() == 23);
static_assert(std::size(sirius_cpu_bridge::transport_program) == 15967);
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
struct Scope {
    unsigned saved = _mm_getcsr();
    Scope() { _mm_setcsr((saved | 0x1f80u) & ~0xe040u); }
    ~Scope() { _mm_setcsr(saved); }
};
void dump(const std::string& name, const Words& words) {
    FILE* file = std::fopen(name.c_str(), "wb");
    require(file, "cannot open raw buffer");
    const bool ok = std::fwrite(words.data(), sizeof(std::uint32_t), words.size(), file) == words.size();
    const int closed = std::fclose(file); require(ok && closed == 0, "cannot preserve raw buffer");
}
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const std::string folder = argv[1];
    Words input(input_words), baseline(output_words, sentinel), candidate(output_words, sentinel);
    std::array<const char*, active> names{};
    FILE* calls = nullptr;
    unsigned completed = 0;
    const auto checkpoint = [&] {
        dump(folder + "/input-after.bin", input);
        dump(folder + "/baseline.bin", baseline); dump(folder + "/candidate.bin", candidate);
    };
    try {
        input[0] = capacity;
        unsigned row = 0;
        const auto copy = [&](const auto& cases) {
            for (const auto& item : cases) {
                require(row < active, "too many fixture rows"); names[row] = item.name;
                std::copy(item.input.begin(), item.input.end(), input.begin() + 1 + row * stride); ++row;
            }
        };
        copy(sirius::test::retained_transport::cases); copy(sirius::test::retained_schwarzschild::cases);
        require(row == active, "missing fixture rows");
        std::copy(std::begin(sirius_cpu_bridge::transport_program), std::end(sirius_cpu_bridge::transport_program),
                  input.begin() + 1 + capacity * stride);
        const Words original = input;
        dump(folder + "/input-before.bin", input); checkpoint();
        FILE* cases = std::fopen((folder + "/cases.tsv").c_str(), "wx");
        require(cases, "cannot preserve case order");
        for (unsigned i = 0; i < active; ++i)
            std::fprintf(cases, "%u\t%s\t%s\n", i, i < 15 ? "general" : "schwarzschild", names[i]);
        require(std::fclose(cases) == 0, "cannot close case order");
        if (!sirius_native_simd_supported()) {
            std::puts("{\"passed\":false,\"not_covered\":\"CPU_OS_ISA\",\"raw_calls\":0}"); return 77;
        }
        calls = std::fopen((folder + "/calls.jsonl").c_str(), "wx");
        require(calls, "cannot preserve raw-call journal");
        using Raw = void (*)(std::uint32_t*, std::uint32_t, std::uint32_t*, std::uint32_t, std::uint32_t, std::uint32_t);
        const std::array<Raw, 2> implementations{baseline_sirius_simd_transport_packet, sirius_simd_transport_packet};
        const std::array<std::array<unsigned, 2>, 5> packets{{{0,8},{8,8},{16,8},{24,8},{32,6}}};
        for (unsigned implementation = 0; implementation < 2; ++implementation) {
            auto& output = implementation ? candidate : baseline;
            for (const auto& packet : packets) {
                require(sirius_native_simd_supported(), "CPU_OS_ISA changed");
                const unsigned saved = _mm_getcsr();
                {
                    Scope scope;
                    const unsigned configured = _mm_getcsr();
                    require((configured & 0xe040u) == 0, "IEEE32 scope not installed");
                    const bool admitted = sirius_native_transport_admit_packet(input.data(), input_words * 4,
                        output.data(), output_words * 4, packet[0], 0, 0, active, 1, 1, packet[1]);
                    require(_mm_getcsr() == configured, "pure guard changed MXCSR");
                    std::fprintf(calls, "{\"implementation\":%u,\"first_x\":%u,\"groups\":%u,"
                        "\"admitted\":%s,\"status\":\"before_raw\"}\n",
                        implementation, packet[0], packet[1], admitted ? "true" : "false");
                    std::fflush(calls); require(admitted, "pure F/B admission refused");
                    implementations[implementation](input.data(), input_words, output.data(), output_words, packet[0], packet[1]);
                    require((_mm_getcsr() & ~0x3fu) == (configured & ~0x3fu), "raw call changed MXCSR controls");
                }
                require(_mm_getcsr() == saved, "MXCSR restoration differs");
                ++completed; checkpoint();
                std::fprintf(calls, "{\"implementation\":%u,\"first_x\":%u,\"groups\":%u,\"status\":\"returned\"}\n",
                             implementation, packet[0], packet[1]); std::fflush(calls);
                require(input == original, "input/program mutated");
            }
        }
        const int closed = std::fclose(calls); calls = nullptr;
        require(closed == 0, "cannot close raw-call journal");
        bool inactive_ok = true;
        for (unsigned i = active * output_stride; i < output_words; ++i)
            inactive_ok &= baseline[i] == sentinel && candidate[i] == sentinel;
        unsigned mismatches = 0, first = output_words;
        for (unsigned i = 0; i < output_words; ++i)
            if (baseline[i] != candidate[i]) { if (!mismatches) first = i; ++mismatches; }
        std::printf("{\"passed\":%s,\"case_count\":38,\"active_rows\":38,\"capacity_rows\":40,"
            "\"packets_per_implementation\":5,\"raw_calls\":%u,\"compared_output_words\":187960,"
            "\"input_unchanged\":true,\"inactive_sentinels_preserved\":%s,\"full_word_match\":%s,"
            "\"mismatch_count\":%u,\"first_mismatch\":",
            (!mismatches && inactive_ok && completed == 10) ? "true" : "false", completed,
            inactive_ok ? "true" : "false", mismatches ? "false" : "true", mismatches);
        if (mismatches) {
            const unsigned row_index = first / output_stride;
            std::printf("{\"index\":%u,\"row\":%u,\"word\":%u,\"case\":\"%s\",\"baseline\":%u,\"candidate\":%u}",
                first, row_index, first % output_stride, row_index < active ? names[row_index] : "inactive",
                baseline[first], candidate[first]);
        } else std::printf("null");
        std::puts("}");
        return mismatches || !inactive_ok || completed != 10 || std::ferror(stdout) ? 1 : 0;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "compatibility failed after%u raw calls: %s\n", completed, error.what());
        try { checkpoint(); } catch (...) {}
        if (calls) std::fclose(calls);
        std::printf("{\"passed\":false,\"raw_calls\":%u,\"failure\":\"guard_FP_input_or_IO\"}\n", completed);
        return 1;
    }
}
