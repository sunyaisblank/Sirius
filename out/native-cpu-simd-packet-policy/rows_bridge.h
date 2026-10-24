#pragma once

// Reuse the original scalar ABI authority; no packet layout is forked here.
#include "../native-cpu-slang-bridge/native_bridge.h"

#ifdef __cplusplus
extern "C" {
#define SIRIUS_ROWS_NOEXCEPT noexcept
#else
#define SIRIUS_ROWS_NOEXCEPT
#endif

#define SIRIUS_ROWS_PACKET_ARGUMENTS \
    const void *input, uint32_t input_bytes, void *output, uint32_t output_bytes, \
    uint32_t first_x, uint32_t y, uint32_t z, uint32_t grid_x, uint32_t grid_y, \
    uint32_t grid_z, uint32_t groups

// These unchanged F exports perform the pinned B pure admission checks.
bool sirius_native_transport_admit_packet(SIRIUS_ROWS_PACKET_ARGUMENTS) SIRIUS_ROWS_NOEXCEPT;
bool sirius_native_endpoint_admit_packet(SIRIUS_ROWS_PACKET_ARGUMENTS) SIRIUS_ROWS_NOEXCEPT;
bool sirius_native_dopri_phase_admit_packet(SIRIUS_ROWS_PACKET_ARGUMENTS) SIRIUS_ROWS_NOEXCEPT;

// Mesa supplies the CPU+OS AVX2 capability query; direct calls use it too.
bool sirius_native_simd_supported(void) SIRIUS_ROWS_NOEXCEPT;

// Full packets on supported CPUs use eight SIMD rows; other calls use the
// original scalar adapter, including its no-write refusal for invalid counts.
bool sirius_native_transport_rows_try_packet(SIRIUS_ROWS_PACKET_ARGUMENTS) SIRIUS_ROWS_NOEXCEPT;
bool sirius_native_endpoint_rows_try_packet(SIRIUS_ROWS_PACKET_ARGUMENTS) SIRIUS_ROWS_NOEXCEPT;
bool sirius_native_dopri_phase_rows_try_packet(SIRIUS_ROWS_PACKET_ARGUMENTS) SIRIUS_ROWS_NOEXCEPT;

#undef SIRIUS_ROWS_PACKET_ARGUMENTS
#undef SIRIUS_ROWS_NOEXCEPT

#ifdef __cplusplus
}
#endif
