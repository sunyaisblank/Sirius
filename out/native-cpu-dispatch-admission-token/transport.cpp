#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include "native_bridge.h"
#include "native_packet.hpp"
#include <transport.native.hpp>
#include <transport.program.hpp>
namespace g = sirius_generated_transport;
static_assert(sizeof(void*) == 8 && sizeof(float) == 4);
static_assert(std::endian::native == std::endian::little);
static_assert(std::numeric_limits<float>::is_iec559 && std::numeric_limits<float>::digits == 24);
static_assert(std::is_standard_layout_v<g::RetainedTriple_0> && sizeof(g::RetainedTriple_0) == 20);
static_assert(offsetof(g::RetainedTriple_0, valid_0) == 16);
static_assert(sizeof(StructuredBuffer<uint32_t>) == 16);
static_assert(sizeof(RWStructuredBuffer<uint32_t>) == 16);
static_assert(offsetof(StructuredBuffer<uint32_t>, data) == 0);
static_assert(offsetof(StructuredBuffer<uint32_t>, count) == 8);
static_assert(offsetof(RWStructuredBuffer<uint32_t>, data) == 0);
static_assert(offsetof(RWStructuredBuffer<uint32_t>, count) == 8);
static_assert(sizeof(g::GlobalParams_0) == 32 && offsetof(g::GlobalParams_0, inputs_0) == 0);
static_assert(offsetof(g::GlobalParams_0, outputs_0) == 16);
static_assert(std::is_standard_layout_v<g::KernelContext_0> && sizeof(g::KernelContext_0) == 9696);
static_assert(offsetof(g::KernelContext_0, retainedScratch_0) == 8);
static_assert(offsetof(g::KernelContext_0, retainedStatus_0) == 9188);
static_assert(offsetof(g::KernelContext_0, retainedCoefficients_0) == 9192);
static_assert(std::is_same_v<decltype(&g::sirius_transportComputeMain_Group), void(*)(ComputeVaryingInput*, void*, void*)>);
static_assert(std::is_same_v<decltype(&g::sirius_transportComputeMain_Thread), void(*)(ComputeThreadVaryingInput*, void*, void*)>);
static_assert(sizeof(sirius_cpu_bridge::transport_program) / sizeof(uint32_t) == 15967);
extern "C" bool sirius_native_transport_try_packet(
    const void* input, uint32_t input_bytes, void* output, uint32_t output_bytes,
    uint32_t first_x, uint32_t y, uint32_t z, uint32_t grid_x, uint32_t grid_y,
    uint32_t grid_z, uint32_t groups) noexcept {
    if (!sirius_cpu_bridge::admit_packet(input, input_bytes, output, output_bytes,
            first_x, y, z, grid_x, grid_y, grid_z, groups, 230, 4699,
            sirius_cpu_bridge::transport_program, 15967))
        return false;
    g::GlobalParams_0 params{};
    params.inputs_0.data = const_cast<uint32_t*>(static_cast<const uint32_t*>(input));
    params.inputs_0.count = input_bytes / sizeof(uint32_t);
    params.outputs_0.data = static_cast<uint32_t*>(output);
    params.outputs_0.count = output_bytes / sizeof(uint32_t);
    for (uint32_t i = 0; i < groups; ++i) {
        ComputeVaryingInput varying{};
        varying.startGroupID.x = first_x + i;
        varying.startGroupID.y = y;
        varying.startGroupID.z = z;
        g::sirius_transportComputeMain_Group(&varying, nullptr, &params);
    }
    return true;
}

#include "admitted_bridge.h"

extern "C" __attribute__((visibility("hidden"))) bool sirius_native_transport_admitted_packet(
    const void* input, uint32_t input_bytes, void* output, uint32_t output_bytes,
    uint32_t first_x, uint32_t y, uint32_t z, uint32_t grid_x, uint32_t grid_y,
    uint32_t grid_z, uint32_t groups) noexcept {
    if (!sirius_cpu_bridge::admit_packet<false, true>(input, input_bytes, output, output_bytes,
            first_x, y, z, grid_x, grid_y, grid_z, groups, 230, 4699,
            sirius_cpu_bridge::transport_program, 15967))
        return false;
    g::GlobalParams_0 params{};
    params.inputs_0.data = const_cast<uint32_t*>(static_cast<const uint32_t*>(input));
    params.inputs_0.count = input_bytes / sizeof(uint32_t);
    params.outputs_0.data = static_cast<uint32_t*>(output);
    params.outputs_0.count = output_bytes / sizeof(uint32_t);
    for (uint32_t i = 0; i < groups; ++i) {
        ComputeVaryingInput varying{};
        varying.startGroupID.x = first_x + i;
        varying.startGroupID.y = y;
        varying.startGroupID.z = z;
        g::sirius_transportComputeMain_Group(&varying, nullptr, &params);
    }
    return true;
}

static_assert(std::is_same_v<decltype(&sirius_native_transport_admitted_packet), decltype(&sirius_native_transport_try_packet)>);
