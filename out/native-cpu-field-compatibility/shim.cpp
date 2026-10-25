#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>

extern "C" bool sirius_native_simd_supported() noexcept;
extern "C" bool sirius_native_transport_field_try_packet(const void*, std::uint32_t,
    void*, std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t,
    std::uint32_t, std::uint32_t, std::uint32_t, std::uint32_t) noexcept;
extern "C" float sirius_simd_fmaxf(float a, float b) noexcept { return ::fmaxf(a, b); }

// The unchanged driver has already admitted this exact packet against active38.
// Its six-argument raw ABI omits grid size; the minimum prefix covering the
// requested packet gives the scalar wrapper its equivalent bounded row extent.
extern "C" void sirius_simd_transport_packet(std::uint32_t* input, std::uint32_t input_words,
    std::uint32_t* output, std::uint32_t output_words, std::uint32_t first, std::uint32_t groups) {
    if (!sirius_native_simd_supported() || groups == 0 || groups > 8 ||
        first > std::numeric_limits<std::uint32_t>::max() - groups ||
        input_words > std::numeric_limits<std::uint32_t>::max() / 4 ||
        output_words > std::numeric_limits<std::uint32_t>::max() / 4 ||
        !sirius_native_transport_field_try_packet(input, input_words * 4, output, output_words * 4,
            first, 0, 0, first + groups, 1, 1, groups))
        std::abort();
}
