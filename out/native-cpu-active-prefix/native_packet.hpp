#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

namespace sirius_cpu_bridge {

// These checks must finish before entering any generated shader or writing output.
// The caller supplies live Vulkan views and keeps the input immutable for dispatch.
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
    return std::memcmp(suffix, program, program_words * sizeof(std::uint32_t)) == 0;
}

} // namespace sirius_cpu_bridge
