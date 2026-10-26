#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <native_bridge.h>
#include "native_packet.hpp"
#include <transport.program.hpp>
#include <endpoint.program.hpp>
#include <dopri_phase.program.hpp>

static_assert(sizeof(sirius_cpu_bridge::transport_program) / sizeof(std::uint32_t) == 15967);
static_assert(sizeof(sirius_cpu_bridge::endpoint_program) / sizeof(std::uint32_t) == 21143);
static_assert(sizeof(sirius_cpu_bridge::dopri_phase_program) / sizeof(std::uint32_t) == 10689);

#define ADMISSION_ARGUMENTS \
    const void* input, std::uint32_t input_bytes, void* output, std::uint32_t output_bytes, \
    std::uint32_t first_x, std::uint32_t y, std::uint32_t z, std::uint32_t grid_x, \
    std::uint32_t grid_y, std::uint32_t grid_z, std::uint32_t groups
#define ADMISSION_VALUES input, input_bytes, output, output_bytes, first_x, y, z, \
    grid_x, grid_y, grid_z, groups

// Pure predicates: the shared original guard reads input and never writes output.
extern "C" bool sirius_native_transport_admit_packet(ADMISSION_ARGUMENTS) noexcept {
    return sirius_cpu_bridge::admit_packet(ADMISSION_VALUES, 230U, 4699U,
        sirius_cpu_bridge::transport_program, 15967U);
}

extern "C" bool sirius_native_endpoint_admit_packet(ADMISSION_ARGUMENTS) noexcept {
    return sirius_cpu_bridge::admit_packet<true>(ADMISSION_VALUES, 225U, 3829U,
        sirius_cpu_bridge::endpoint_program, 21143U);
}

extern "C" bool sirius_native_dopri_phase_admit_packet(ADMISSION_ARGUMENTS) noexcept {
    return sirius_cpu_bridge::admit_packet(ADMISSION_VALUES, 1810U, 4149U,
        sirius_cpu_bridge::dopri_phase_program, 10689U);
}

static_assert(std::is_same_v<decltype(&sirius_native_transport_admit_packet),
                             decltype(&sirius_native_transport_try_packet)>);
static_assert(std::is_same_v<decltype(&sirius_native_endpoint_admit_packet),
                             decltype(&sirius_native_endpoint_try_packet)>);
static_assert(std::is_same_v<decltype(&sirius_native_dopri_phase_admit_packet),
                             decltype(&sirius_native_dopri_phase_try_packet)>);

#undef ADMISSION_VALUES
#undef ADMISSION_ARGUMENTS
