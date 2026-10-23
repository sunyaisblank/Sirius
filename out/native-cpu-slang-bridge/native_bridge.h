#pragma once

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#define SIRIUS_NATIVE_NOEXCEPT noexcept
#else
#define SIRIUS_NATIVE_NOEXCEPT
#endif

// False means no output byte was written. The caller must invoke the ordinary JIT.
// True means every admitted logical group executed all 64 original invocations.
// Requires live view storage, stable input, and the existing IEEE32 callback guard.
#define SIRIUS_NATIVE_PACKET_ARGUMENTS \
    const void *input, uint32_t input_bytes, void *output, uint32_t output_bytes, \
    uint32_t first_x, uint32_t y, uint32_t z, uint32_t grid_x, uint32_t grid_y, \
    uint32_t grid_z, uint32_t groups

bool sirius_native_transport_try_packet(SIRIUS_NATIVE_PACKET_ARGUMENTS) SIRIUS_NATIVE_NOEXCEPT;
bool sirius_native_endpoint_try_packet(SIRIUS_NATIVE_PACKET_ARGUMENTS) SIRIUS_NATIVE_NOEXCEPT;
bool sirius_native_dopri_phase_try_packet(SIRIUS_NATIVE_PACKET_ARGUMENTS) SIRIUS_NATIVE_NOEXCEPT;

#undef SIRIUS_NATIVE_PACKET_ARGUMENTS
#undef SIRIUS_NATIVE_NOEXCEPT

#ifdef __cplusplus
}
#endif
