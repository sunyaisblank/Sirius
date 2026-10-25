#include "../native-cpu-simd-packet-policy/rows_bridge.h"
#include <cmath>
#include <type_traits>

extern "C" float sirius_simd_fmaxf(float a, float b) noexcept { return ::fmaxf(a, b); }

#define ROWS_ARGUMENTS \
    const void *input, uint32_t input_bytes, void *output, uint32_t output_bytes, \
    uint32_t first_x, uint32_t y, uint32_t z, uint32_t grid_x, uint32_t grid_y, \
    uint32_t grid_z, uint32_t groups
#define ROWS_VALUES input, input_bytes, output, output_bytes, first_x, y, z, \
    grid_x, grid_y, grid_z, groups

// AVX2 runs within each row's independent field nodes. Keep the existing OS/ISA
// admission ahead of every field entry and the original D scalar fallback.
// Pure F/B admission still precedes all native output writes.
#define DEFINE_FIELD_ROWS(stage) \
    extern "C" bool sirius_native_##stage##_field_try_packet(ROWS_ARGUMENTS) noexcept; \
    extern "C" bool sirius_native_##stage##_rows_try_packet(ROWS_ARGUMENTS) noexcept { \
        if (!sirius_native_simd_supported()) \
            return sirius_native_##stage##_try_packet(ROWS_VALUES); \
        if (!sirius_native_##stage##_admit_packet(ROWS_VALUES)) \
            return false; \
        return sirius_native_##stage##_field_try_packet(ROWS_VALUES); \
    } \
    static_assert(std::is_same_v<decltype(&sirius_native_##stage##_field_try_packet), \
                                 decltype(&sirius_native_##stage##_try_packet)>); \
    static_assert(std::is_same_v<decltype(&sirius_native_##stage##_rows_try_packet), \
                                 decltype(&sirius_native_##stage##_try_packet)>);
DEFINE_FIELD_ROWS(transport)
DEFINE_FIELD_ROWS(endpoint)
DEFINE_FIELD_ROWS(dopri_phase)
#undef DEFINE_FIELD_ROWS
#undef ROWS_VALUES
#undef ROWS_ARGUMENTS
