#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <endpoint.program.hpp>

namespace sirius_cpu_bridge {

// These checks must finish before entering any generated shader or writing output.
// The caller supplies live Vulkan views and keeps the input immutable for dispatch.
// dispatch_program_admitted is private: the identical live views and immutable
// program have passed complete admission in this synchronous Mesa job.
// It skips only program-byte comparisons; every other guard below remains.
template<bool endpoint_header_domain = false, bool dispatch_program_admitted = false>
inline bool admit_packet(const void* input, std::uint32_t input_bytes,
                         void* output, std::uint32_t output_bytes,
                         std::uint32_t first_x, std::uint32_t y, std::uint32_t z,
                         std::uint32_t grid_x, std::uint32_t grid_y,
                         std::uint32_t grid_z, std::uint32_t groups,
                         std::uint32_t input_row_words,
                         std::uint32_t output_row_words,
                         const std::uint32_t* program, std::size_t program_words) noexcept {
    if (!input || !output || !program || !program_words || program_words > 65536 ||
        input_bytes < 4 || ((input_bytes | output_bytes) & 3u) ||
        y || z || grid_y != 1 || grid_z != 1 || !groups || groups > 8)
        return false;
    const auto in = reinterpret_cast<std::uintptr_t>(input);
    const auto out = reinterpret_cast<std::uintptr_t>(output);
    constexpr auto maximum = std::numeric_limits<std::uintptr_t>::max();
    if (((in | out) & 3u) || in > maximum - input_bytes || out > maximum - output_bytes)
        return false;
    if (in < out + output_bytes && out < in + input_bytes)
        return false;
    std::uint32_t rows_word;
    std::memcpy(&rows_word, input, sizeof(rows_word));
    const std::uint64_t rows = rows_word;
    // Header rows are immutable allocation capacity, not active dispatch size.
    // Program/stride extents keep that capacity; only the submitted prefix runs.
    if (!rows || rows > 65536 || !grid_x || grid_x > rows ||
        std::uint64_t(first_x) + groups > grid_x)
        return false;
    const std::uint64_t program_base = 1 + rows * input_row_words;
    const std::uint64_t required_input_bytes = (program_base + program_words) * 4;
    const std::uint64_t required_output_bytes = rows * output_row_words * 4;
    if (required_input_bytes != input_bytes || required_output_bytes != output_bytes)
        return false;
    const auto* suffix = static_cast<const std::byte*>(input) + program_base * 4;
    if constexpr (endpoint_header_domain) {
        // The unchanged shader compares this word only for the specialized
        // zero-spin/charge/lambda family, refusing that family on a mismatch.
        // It must still execute: internal and public output words are observable.
        constexpr std::size_t specialized_base = 11240;
        constexpr std::size_t register_word = specialized_base + 1;
        static_assert(sizeof(endpoint_program) / sizeof(std::uint32_t) == 21143);
        static_assert(endpoint_program[1] == 612 &&
                      endpoint_program[register_word] == endpoint_program[1]);
        static_assert(endpoint_program[specialized_base] == 1943 &&
                      endpoint_program[specialized_base + 2] == 100);
        if (input_row_words != 225 || output_row_words != 3829 ||
            program != endpoint_program || program_words != 21143)
            return false;
        if constexpr (dispatch_program_admitted) return true;
        return std::memcmp(suffix, program, register_word * sizeof(std::uint32_t)) == 0 &&
               std::memcmp(suffix + (register_word + 1) * sizeof(std::uint32_t),
                           program + register_word + 1,
                           (program_words - register_word - 1) * sizeof(std::uint32_t)) == 0;
    }
    if constexpr (dispatch_program_admitted) return true;
    return std::memcmp(suffix, program, program_words * sizeof(std::uint32_t)) == 0;
}

} // namespace sirius_cpu_bridge
