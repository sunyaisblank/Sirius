#include "../native-cpu-simd-packet-policy/rows_bridge.h"
#include <cmath>
#include <type_traits>

// A unique external name keeps ISPC from replacing libm fmaxf with its own
// intrinsic. The original Slang/GCC prelude calls this same scalar function.
extern "C" float sirius_simd_fmaxf(float a, float b) noexcept {
    return ::fmaxf(a, b);
}

// ISPC exports take word counts, while the original packet ABI takes bytes.
extern "C" void sirius_simd_transport_packet(uint32_t *input, uint32_t input_words,
    uint32_t *output, uint32_t output_words, uint32_t first_x, uint32_t groups);
extern "C" void sirius_simd_endpoint_packet(uint32_t *input, uint32_t input_words,
    uint32_t *output, uint32_t output_words, uint32_t first_x, uint32_t groups);
extern "C" void sirius_simd_dopri_phase_packet(uint32_t *input, uint32_t input_words,
    uint32_t *output, uint32_t output_words, uint32_t first_x, uint32_t groups);

extern "C" bool sirius_native_v4_supported() noexcept;
extern "C" void sirius_simd_v4_transport_packet(uint32_t *input, uint32_t input_words,
    uint32_t *output, uint32_t output_words, uint32_t first_x, uint32_t groups);
extern "C" void sirius_simd_v4_endpoint_packet(uint32_t *input, uint32_t input_words,
    uint32_t *output, uint32_t output_words, uint32_t first_x, uint32_t groups);
extern "C" void sirius_simd_v4_dopri_phase_packet(uint32_t *input, uint32_t input_words,
    uint32_t *output, uint32_t output_words, uint32_t first_x, uint32_t groups);

static_assert(sizeof(uint32_t) == 4);

#define ROWS_ARGUMENTS \
    const void *input, uint32_t input_bytes, void *output, uint32_t output_bytes, \
    uint32_t first_x, uint32_t y, uint32_t z, uint32_t grid_x, uint32_t grid_y, \
    uint32_t grid_z, uint32_t groups
#define ROWS_VALUES input, input_bytes, output, output_bytes, first_x, y, z, \
    grid_x, grid_y, grid_z, groups

// ISPC's input pointer type is mutable, but the sealed derived shader only
// reads it. Remove const solely at that raw ABI boundary after full admission.
extern "C" bool sirius_native_transport_rows_try_packet(ROWS_ARGUMENTS) noexcept {
    if (groups != 8U || !sirius_native_simd_supported())
        return sirius_native_transport_try_packet(ROWS_VALUES);
    if (!sirius_native_transport_admit_packet(ROWS_VALUES))
        return false;
    if (sirius_native_v4_supported())
        sirius_simd_v4_transport_packet(const_cast<uint32_t*>(static_cast<const uint32_t*>(input)),
        input_bytes / 4U, static_cast<uint32_t*>(output), output_bytes / 4U, first_x, groups);
    else
        sirius_simd_transport_packet(const_cast<uint32_t*>(static_cast<const uint32_t*>(input)),
        input_bytes / 4U, static_cast<uint32_t*>(output), output_bytes / 4U, first_x, groups);
    return true;
}

extern "C" bool sirius_native_endpoint_rows_try_packet(ROWS_ARGUMENTS) noexcept {
    if (groups != 8U || !sirius_native_simd_supported())
        return sirius_native_endpoint_try_packet(ROWS_VALUES);
    if (!sirius_native_endpoint_admit_packet(ROWS_VALUES))
        return false;
    if (sirius_native_v4_supported())
        sirius_simd_v4_endpoint_packet(const_cast<uint32_t*>(static_cast<const uint32_t*>(input)),
        input_bytes / 4U, static_cast<uint32_t*>(output), output_bytes / 4U, first_x, groups);
    else
        sirius_simd_endpoint_packet(const_cast<uint32_t*>(static_cast<const uint32_t*>(input)),
        input_bytes / 4U, static_cast<uint32_t*>(output), output_bytes / 4U, first_x, groups);
    return true;
}

extern "C" bool sirius_native_dopri_phase_rows_try_packet(ROWS_ARGUMENTS) noexcept {
    if (groups != 8U || !sirius_native_simd_supported())
        return sirius_native_dopri_phase_try_packet(ROWS_VALUES);
    if (!sirius_native_dopri_phase_admit_packet(ROWS_VALUES))
        return false;
    if (sirius_native_v4_supported())
        sirius_simd_v4_dopri_phase_packet(const_cast<uint32_t*>(static_cast<const uint32_t*>(input)),
        input_bytes / 4U, static_cast<uint32_t*>(output), output_bytes / 4U, first_x, groups);
    else
        sirius_simd_dopri_phase_packet(const_cast<uint32_t*>(static_cast<const uint32_t*>(input)),
        input_bytes / 4U, static_cast<uint32_t*>(output), output_bytes / 4U, first_x, groups);
    // Arithmetic refusal remains handled execution, represented in output.
    return true;
}

#define CHECK_ROWS_ABI(stage) \
    static_assert(std::is_same_v<decltype(&sirius_native_##stage##_rows_try_packet), \
                                 decltype(&sirius_native_##stage##_try_packet)>); \
    static_assert(std::is_same_v<decltype(&sirius_native_##stage##_admit_packet), \
                                 decltype(&sirius_native_##stage##_try_packet)>);
CHECK_ROWS_ABI(transport)
CHECK_ROWS_ABI(endpoint)
CHECK_ROWS_ABI(dopri_phase)
#undef CHECK_ROWS_ABI
#undef ROWS_VALUES
#undef ROWS_ARGUMENTS

#include "admitted_bridge.h"
#include "native_packet.hpp"
#include <transport.program.hpp>
#include <dopri_phase.program.hpp>

#define ROWS_ARGUMENTS \
    const void *input, uint32_t input_bytes, void *output, uint32_t output_bytes, \
    uint32_t first_x, uint32_t y, uint32_t z, uint32_t grid_x, uint32_t grid_y, \
    uint32_t grid_z, uint32_t groups
#define ROWS_VALUES input, input_bytes, output, output_bytes, first_x, y, z, \
    grid_x, grid_y, grid_z, groups


extern "C" __attribute__((visibility("hidden"))) bool sirius_native_transport_rows_admitted_packet(ROWS_ARGUMENTS) noexcept {
    if (groups != 8U || !sirius_native_simd_supported())
        return sirius_native_transport_admitted_packet(ROWS_VALUES);
    if (!sirius_cpu_bridge::admit_packet<false, true>(ROWS_VALUES, 230U, 4699U,
        sirius_cpu_bridge::transport_program, 15967U))
        return false;
    if (sirius_native_v4_supported())
        sirius_simd_v4_transport_packet(const_cast<uint32_t*>(static_cast<const uint32_t*>(input)),
        input_bytes / 4U, static_cast<uint32_t*>(output), output_bytes / 4U, first_x, groups);
    else
        sirius_simd_transport_packet(const_cast<uint32_t*>(static_cast<const uint32_t*>(input)),
        input_bytes / 4U, static_cast<uint32_t*>(output), output_bytes / 4U, first_x, groups);
    return true;
}
static_assert(std::is_same_v<decltype(&sirius_native_transport_rows_admitted_packet), decltype(&sirius_native_transport_rows_try_packet)>);

extern "C" __attribute__((visibility("hidden"))) bool sirius_native_endpoint_rows_admitted_packet(ROWS_ARGUMENTS) noexcept {
    if (groups != 8U || !sirius_native_simd_supported())
        return sirius_native_endpoint_admitted_packet(ROWS_VALUES);
    if (!sirius_cpu_bridge::admit_packet<true, true>(ROWS_VALUES, 225U, 3829U,
        sirius_cpu_bridge::endpoint_program, 21143U))
        return false;
    if (sirius_native_v4_supported())
        sirius_simd_v4_endpoint_packet(const_cast<uint32_t*>(static_cast<const uint32_t*>(input)),
        input_bytes / 4U, static_cast<uint32_t*>(output), output_bytes / 4U, first_x, groups);
    else
        sirius_simd_endpoint_packet(const_cast<uint32_t*>(static_cast<const uint32_t*>(input)),
        input_bytes / 4U, static_cast<uint32_t*>(output), output_bytes / 4U, first_x, groups);
    return true;
}
static_assert(std::is_same_v<decltype(&sirius_native_endpoint_rows_admitted_packet), decltype(&sirius_native_endpoint_rows_try_packet)>);

extern "C" __attribute__((visibility("hidden"))) bool sirius_native_dopri_phase_rows_admitted_packet(ROWS_ARGUMENTS) noexcept {
    if (groups != 8U || !sirius_native_simd_supported())
        return sirius_native_dopri_phase_admitted_packet(ROWS_VALUES);
    if (!sirius_cpu_bridge::admit_packet<false, true>(ROWS_VALUES, 1810U, 4149U,
        sirius_cpu_bridge::dopri_phase_program, 10689U))
        return false;
    if (sirius_native_v4_supported())
        sirius_simd_v4_dopri_phase_packet(const_cast<uint32_t*>(static_cast<const uint32_t*>(input)),
        input_bytes / 4U, static_cast<uint32_t*>(output), output_bytes / 4U, first_x, groups);
    else
        sirius_simd_dopri_phase_packet(const_cast<uint32_t*>(static_cast<const uint32_t*>(input)),
        input_bytes / 4U, static_cast<uint32_t*>(output), output_bytes / 4U, first_x, groups);
    // Arithmetic refusal remains handled execution, represented in output.
    return true;
}
static_assert(std::is_same_v<decltype(&sirius_native_dopri_phase_rows_admitted_packet), decltype(&sirius_native_dopri_phase_rows_try_packet)>);

#undef ROWS_VALUES
#undef ROWS_ARGUMENTS
