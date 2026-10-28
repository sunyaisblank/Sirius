#pragma once

// Private entry contract: only lp_native_cpu_packet_mode's positive job-local
// native_preflight_passed branch may call rows_admitted_packet. That preflight
// validated the same stage, live views/extents, immutable capacity/grid and all
// required canonical program bytes before any output write. The input and
// descriptor table remain immutable through the existing synchronous pool wait.
// Output is disjoint from the whole input and descriptor table. Each worker
// still checks all original pointer/extent/header/shape/alias/packet bounds.
// Arithmetic refusals remain handled execution, never late fallback.
// Public try_packet/admit_packet entries remain fully guarded; no token survives
// this job and no pointer/global/TLS cache or caller-supplied bool is introduced.
#include "../native-cpu-simd-packet-policy/rows_bridge.h"

#ifdef __cplusplus
extern "C" {
#define SIRIUS_ADMITTED_NOEXCEPT noexcept
#else
#define SIRIUS_ADMITTED_NOEXCEPT
#endif
#define SIRIUS_ADMITTED_ARGUMENTS \
    const void *input, uint32_t input_bytes, void *output, uint32_t output_bytes, \
    uint32_t first_x, uint32_t y, uint32_t z, uint32_t grid_x, uint32_t grid_y, \
    uint32_t grid_z, uint32_t groups
__attribute__((visibility("hidden")))
bool sirius_native_transport_admitted_packet(SIRIUS_ADMITTED_ARGUMENTS) SIRIUS_ADMITTED_NOEXCEPT;
__attribute__((visibility("hidden")))
bool sirius_native_transport_rows_admitted_packet(SIRIUS_ADMITTED_ARGUMENTS) SIRIUS_ADMITTED_NOEXCEPT;
__attribute__((visibility("hidden")))
bool sirius_native_endpoint_admitted_packet(SIRIUS_ADMITTED_ARGUMENTS) SIRIUS_ADMITTED_NOEXCEPT;
__attribute__((visibility("hidden")))
bool sirius_native_endpoint_rows_admitted_packet(SIRIUS_ADMITTED_ARGUMENTS) SIRIUS_ADMITTED_NOEXCEPT;
__attribute__((visibility("hidden")))
bool sirius_native_dopri_phase_admitted_packet(SIRIUS_ADMITTED_ARGUMENTS) SIRIUS_ADMITTED_NOEXCEPT;
__attribute__((visibility("hidden")))
bool sirius_native_dopri_phase_rows_admitted_packet(SIRIUS_ADMITTED_ARGUMENTS) SIRIUS_ADMITTED_NOEXCEPT;
#undef SIRIUS_ADMITTED_ARGUMENTS
#undef SIRIUS_ADMITTED_NOEXCEPT
#ifdef __cplusplus
}
#endif
