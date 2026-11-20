#include "../native-cpu-dispatch-admission-token/admitted_bridge.h"
#include "../native-cpu-dispatch-admission-token/native_packet.hpp"
#include <dopri_phase.program.hpp>
#include <transport.program.hpp>
#include <type_traits>

using WideWorker = void(uint32_t*, uint32_t, uint32_t*, uint32_t, uint32_t, uint32_t) noexcept;
extern "C" __attribute__((visibility("hidden"))) WideWorker sirius_wide_strict_transport_packet;
extern "C" __attribute__((visibility("hidden"))) WideWorker sirius_wide_strict_endpoint_packet;
extern "C" __attribute__((visibility("hidden"))) WideWorker sirius_wide_strict_dp_packet;

namespace {

// Both halves use the full original views. Validate both before the one
// transactional publication so a short high half cannot commit the low half.
template<bool EndpointHeader, WideWorker& Worker>
bool admit_and_run(const void* input, uint32_t input_bytes, void* output, uint32_t output_bytes,
                   uint32_t first_x, uint32_t y, uint32_t z, uint32_t grid_x, uint32_t grid_y,
                   uint32_t grid_z, uint32_t groups, uint32_t input_row_words,
                   uint32_t output_row_words, const uint32_t* program, std::size_t program_words) noexcept {
    constexpr uint32_t half = 8;
    if (groups != 2 * half ||
        !sirius_cpu_bridge::admit_packet<EndpointHeader, true>(
            input, input_bytes, output, output_bytes, first_x, y, z, grid_x, grid_y, grid_z,
            half, input_row_words, output_row_words, program, program_words) ||
        !sirius_cpu_bridge::admit_packet<EndpointHeader, true>(
            input, input_bytes, output, output_bytes, first_x + half, y, z, grid_x, grid_y, grid_z,
            half, input_row_words, output_row_words, program, program_words))
        return false;

    // The raw ABI uses word counts and only reads input.
    Worker(const_cast<uint32_t*>(static_cast<const uint32_t*>(input)), input_bytes / sizeof(uint32_t),
           static_cast<uint32_t*>(output), output_bytes / sizeof(uint32_t), first_x, groups);
    return true;
}

} // namespace

// Called from the synchronous positive-preflight branch, after the width-16
// CPU/OS gate and inside the existing saved and restored IEEE32 scope.
#define ROWS_ARGUMENTS \
    const void* input, uint32_t input_bytes, void* output, uint32_t output_bytes, \
    uint32_t first_x, uint32_t y, uint32_t z, uint32_t grid_x, uint32_t grid_y, \
    uint32_t grid_z, uint32_t groups
#define ROWS_VALUES input, input_bytes, output, output_bytes, first_x, y, z, grid_x, grid_y, grid_z, groups

extern "C" __attribute__((visibility("hidden")))
bool sirius_native_transport_rows_wide_admitted_packet(ROWS_ARGUMENTS) noexcept {
    return admit_and_run<false, sirius_wide_strict_transport_packet>(
        ROWS_VALUES, 230U, 4699U, sirius_cpu_bridge::transport_program, 15967U);
}

extern "C" __attribute__((visibility("hidden")))
bool sirius_native_endpoint_rows_wide_admitted_packet(ROWS_ARGUMENTS) noexcept {
    return admit_and_run<true, sirius_wide_strict_endpoint_packet>(
        ROWS_VALUES, 225U, 3829U, sirius_cpu_bridge::endpoint_program, 21143U);
}

extern "C" __attribute__((visibility("hidden")))
bool sirius_native_dopri_phase_rows_wide_admitted_packet(ROWS_ARGUMENTS) noexcept {
    return admit_and_run<false, sirius_wide_strict_dp_packet>(
        ROWS_VALUES, 1810U, 4149U, sirius_cpu_bridge::dopri_phase_program, 10689U);
}

static_assert(std::is_same_v<decltype(&sirius_native_transport_rows_wide_admitted_packet),
                             decltype(&sirius_native_transport_rows_admitted_packet)>);
static_assert(std::is_same_v<decltype(&sirius_native_endpoint_rows_wide_admitted_packet),
                             decltype(&sirius_native_endpoint_rows_admitted_packet)>);
static_assert(std::is_same_v<decltype(&sirius_native_dopri_phase_rows_wide_admitted_packet),
                             decltype(&sirius_native_dopri_phase_rows_admitted_packet)>);

#undef ROWS_VALUES
#undef ROWS_ARGUMENTS
