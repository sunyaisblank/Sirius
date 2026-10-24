#include <algorithm>
#include <array>
#include <bit>
#include <cerrno>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <immintrin.h>
#include <limits>
#include <stdexcept>
#include <vector>

__attribute__((noinline, noipa)) double reference_spacing(float value) {
    const auto spacing = [](float value) {
        const float magnitude = std::abs(value);
        return double(std::nextafter(magnitude, std::numeric_limits<float>::infinity())) -
               double(magnitude);
    };
    return spacing(value);
}
__attribute__((noinline, noipa)) double candidate_spacing(float value) {
    const auto spacing = [](float value) {
        static_assert(sizeof(float) == sizeof(std::uint32_t) &&
                      std::numeric_limits<float>::is_iec559 &&
                      std::numeric_limits<float>::radix == 2 &&
                      std::numeric_limits<float>::digits == 24 &&
                      std::numeric_limits<float>::min_exponent == -125 &&
                      std::numeric_limits<float>::max_exponent == 128 &&
                      std::bit_cast<std::uint32_t>(1.0f) == 0x3f800000u);
        static_assert(sizeof(double) == sizeof(std::uint64_t) &&
                      std::numeric_limits<double>::is_iec559 &&
                      std::numeric_limits<double>::radix == 2 &&
                      std::numeric_limits<double>::digits == 53 &&
                      std::numeric_limits<double>::min_exponent == -1021 &&
                      std::numeric_limits<double>::max_exponent == 1024 &&
                      std::bit_cast<std::uint64_t>(1.0) == 0x3ff0000000000000ull);
        const float magnitude = std::abs(value);
        const auto bits = std::bit_cast<std::uint32_t>(magnitude);
        const auto exponent = (bits >> 23) & 255u;
        if (exponent != 0 && bits < 0x7f7fffffu)
            return std::bit_cast<double>(std::uint64_t(exponent + 873u) << 52);
        return double(std::nextafter(magnitude, std::numeric_limits<float>::infinity())) -
               double(magnitude);
    };
    return spacing(value);
}

struct Observed {
    std::uint64_t bits;
    int error, flags;
    unsigned csr;
    int rounding;
    bool operator==(const Observed&) const = default;
};
struct Failure {
    std::uint32_t input;
    int rounding, seed;
    unsigned controls;
    Observed reference, candidate;
};
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
__attribute__((noinline, noipa))
Observed observe(double (*function)(float), std::uint32_t input, int rounding,
                 unsigned controls, int seed, const std::fenv_t& baseline) {
    require(std::fesetenv(&baseline) == 0 && std::fesetround(rounding) == 0,
            "cannot install environment");
    require(std::feclearexcept(FE_ALL_EXCEPT) == 0, "cannot clear flags");
    if (seed) require(std::feraiseexcept(FE_INEXACT) == 0, "cannot seed inexact");
    _mm_setcsr((_mm_getcsr() & ~(0x8040u | 0x3fu)) | controls | (seed ? 0x20u : 0u));
    require(std::fetestexcept(FE_ALL_EXCEPT) == (seed ? FE_INEXACT : 0) &&
            std::fegetround() == rounding && (_mm_getcsr() & 0x8040u) == controls,
            "environment precondition differs");
    const float value = std::bit_cast<float>(input);
    int* const error_pointer = &errno;
    *error_pointer = EDOM;  // A nonzero sentinel detects unexpected clearing.
    asm volatile("" ::: "memory");
    const double result = function(value);
    asm volatile("" ::: "memory");
    const int error = *error_pointer;
    const unsigned csr = _mm_getcsr();
    const int flags = std::fetestexcept(FE_ALL_EXCEPT);
    const int mode = std::fegetround();
    return {std::bit_cast<std::uint64_t>(result), error, flags, csr, mode};
}
void print_observed(const Observed& value) {
    std::printf("{\"bits\":\"%016llx\",\"errno\":%d,\"FEflags\":%d,\"MXCSR\":%u,\"round\":%d}",
                static_cast<unsigned long long>(value.bits), value.error, value.flags, value.csr, value.rounding);
}
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    int* const error_pointer = &errno;
    const int saved_errno = *error_pointer;
    std::fenv_t ambient, baseline;
    if (std::fegetenv(&ambient) != 0) return 2;
    const unsigned ambient_csr = _mm_getcsr();
    const int ambient_flags = std::fetestexcept(FE_ALL_EXCEPT);
    const int ambient_round = std::fegetround();
    const auto restore = [&] {
        const bool success = std::fesetenv(&ambient) == 0;
        _mm_setcsr(ambient_csr); *error_pointer = saved_errno;
        return success;
    };
    std::uint64_t comparisons = 0, mismatches = 0;
    std::size_t input_count = 0;
    unsigned mxcsr_mask = 0;
    bool trap_masked = false;
    const char* setup_error = nullptr;
    std::vector<Failure> failures;
    try {
        // Save actual mapped paths; Python later pins them against preparation.
        FILE* maps = std::fopen("/proc/self/maps", "r");
        FILE* copy = std::fopen(argv[1], "wx");
        require(maps && copy, "cannot preserve actual process maps");
        char line[4096];
        while (std::fgets(line, sizeof(line), maps)) std::fputs(line, copy);
        const bool io_ok = !std::ferror(maps) && !std::ferror(copy);
        const int close_maps = std::fclose(maps), close_copy = std::fclose(copy);
        require(io_ok && close_maps == 0 && close_copy == 0, "cannot finish process maps");
        alignas(16) unsigned char state[512]{};
        _fxsave(state); std::memcpy(&mxcsr_mask, state + 28, sizeof(mxcsr_mask));
        if (!mxcsr_mask) mxcsr_mask = 0xffbfu;
        require((mxcsr_mask & 0x8040u) == 0x8040u, "FTZ/DAZ unsupported");
        require(std::feholdexcept(&baseline) == 0, "cannot mask traps");
        trap_masked = true;
        require(std::fegetenv(&baseline) == 0, "cannot preserve masked baseline");
        std::vector<std::uint32_t> inputs;
        const auto add_signed = [&](std::uint32_t bits) {
            inputs.push_back(bits); inputs.push_back(bits | 0x80000000u);
        };
        for (unsigned exponent = 1; exponent <= 254; ++exponent)
            for (unsigned mantissa : {0u, 1u, 0x400000u, 0x7ffffeu, 0x7fffffu})
                add_signed((exponent << 23) | mantissa);
        for (unsigned bits : {0u, 1u, 2u, 0x003fffffu, 0x00400000u, 0x007ffffeu,
                              0x007fffffu, 0x00800000u, 0x00800001u, 0x00800002u,
                              0x7f7ffffeu, 0x7f7fffffu, 0x7f800000u, 0x7fc00000u,
                              0x7fc00001u, 0x7fffffffu, 0x7f800001u, 0x7f800002u,
                              0x7fa00000u, 0x7fbfffffu})
            add_signed(bits);
        std::sort(inputs.begin(), inputs.end());
        inputs.erase(std::unique(inputs.begin(), inputs.end()), inputs.end());
        input_count = inputs.size();
        for (int rounding : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO})
            for (unsigned mode = 0; mode != 4; ++mode)
                for (int seed : {0, 1}) {
                    const unsigned controls = (mode & 1 ? 0x40u : 0u) | (mode & 2 ? 0x8000u : 0u);
                    for (std::uint32_t input : inputs) {
                        const auto a = observe(reference_spacing, input, rounding, controls, seed, baseline);
                        const auto b = observe(candidate_spacing, input, rounding, controls, seed, baseline);
                        ++comparisons;
                        if (!(a == b)) {
                            ++mismatches;
                            if (failures.size() < 32) failures.push_back({input, rounding, seed, controls, a, b});
                        }
                    }
                }
    } catch (const std::exception& error) {
        std::fprintf(stderr, "setup/probe error: %s\n", error.what());
        setup_error = "setup_or_probe_exception";
    } catch (...) {
        setup_error = "unknown_exception";
    }
    const bool restore_ok = restore();
    const bool restored = restore_ok && _mm_getcsr() == ambient_csr &&
                          std::fetestexcept(FE_ALL_EXCEPT) == ambient_flags &&
                          std::fegetround() == ambient_round && *error_pointer == saved_errno;
    const bool passed = !setup_error && comparisons == input_count * 32 && !mismatches && restored;
    std::printf("{\"passed\":%s,\"unique_inputs\":%zu,\"environments\":32,\"comparisons\":%llu,"
                "\"function_calls\":%llu,\"mismatches\":%llu,\"MXCSR_supported_mask\":%u,"
                "\"trap_masked\":%s,\"errno_seed\":%d,\"ambient_restored\":%s,\"setup_error\":%s,"
                "\"first_mismatches\":[", passed ? "true" : "false", input_count,
                static_cast<unsigned long long>(comparisons), static_cast<unsigned long long>(comparisons * 2),
                static_cast<unsigned long long>(mismatches), mxcsr_mask, trap_masked ? "true" : "false",
                EDOM, restored ? "true" : "false",
                setup_error ? "\"setup_or_probe_exception\"" : "null");
    for (std::size_t i = 0; i < failures.size(); ++i) {
        const auto& failure = failures[i];
        std::printf("%s{\"input\":\"%08x\",\"round\":%d,\"FTZ_DAZ\":%u,\"preexisting_inexact\":%d,"
                    "\"reference\":", i ? "," : "", failure.input, failure.rounding, failure.controls, failure.seed);
        print_observed(failure.reference); std::printf(",\"candidate\":");
        print_observed(failure.candidate); std::printf("}");
    }
    std::printf("]}\n");
    const bool output_ok = !std::ferror(stdout);
    const bool final_restore = restore();
    return passed && output_ok && final_restore ? 0 : 1;
}
