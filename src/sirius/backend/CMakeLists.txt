# Backend module. The CPU tracer is an always-available backend with no device
# loader dependency. Vulkan and its Slang artefacts are optional additions;
# machines without either still build the complete CPU system.

add_library(sirius_backend_cpu STATIC
    cpu/geodesic_tracer.cpp
)

target_link_libraries(sirius_backend_cpu PUBLIC sirius_core)

sirius_configure_target(sirius_backend_cpu)

# Vulkan device seam plus the Slang kernel build (single source authority in
# ../kernels compiled to SPIR-V here). Both degrade with a clear message rather
# than failing configuration.

find_package(Vulkan QUIET)
if(NOT Vulkan_FOUND)
    message(STATUS "Sirius backend: Vulkan development files not found; Vulkan backend skipped")
    return()
endif()

add_library(sirius_backend STATIC
    vulkan/vulkan_device.cpp
)

target_link_libraries(sirius_backend PUBLIC sirius_base Vulkan::Vulkan ${CMAKE_DL_LIBS})

sirius_configure_target(sirius_backend)

# --- Slang kernels -> SPIR-V --------------------------------------------------

find_program(SIRIUS_SLANGC slangc
    HINTS "$ENV{SLANG_ROOT}/bin" "${SLANG_ROOT}/bin" /opt/slang/bin)

if(NOT SIRIUS_SLANGC)
    message(STATUS "Sirius backend: slangc not found; kernels will not be compiled "
                   "(set SLANG_ROOT to enable)")
    return()
endif()

set(SIRIUS_KERNEL_BINARY_DIR "${CMAKE_BINARY_DIR}/kernels" CACHE INTERNAL "")
file(MAKE_DIRECTORY "${SIRIUS_KERNEL_BINARY_DIR}")

set(SIRIUS_KERNEL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../kernels")

# Entry-point kernels (each carries a ComputeMain). The gr_* modules below are
# imported libraries, not standalone entry points, so they are recompilation
# dependencies rather than build targets.
set(SIRIUS_KERNEL_SOURCES
    smoke
    parity_probe
    trace
)

# Shared physics modules imported (via `import`, -I include path) by the
# entry-point kernels; any edit re-triggers the affected kernel compiles.
file(GLOB SIRIUS_KERNEL_LIBS CONFIGURE_DEPENDS "${SIRIUS_KERNEL_DIR}/gr_*.slang")

# spirv-val gates every emitted module when present (verification gate,
# docs/SPECIFICATION.md section 7); when absent the build proceeds and the note
# is printed once.
find_program(SIRIUS_SPIRV_VAL spirv-val)
if(NOT SIRIUS_SPIRV_VAL)
    message(STATUS "Sirius backend: spirv-val not found; emitted kernels will not be validated at build time")
endif()

# Bounded retained programs are embedded with their validated shaders. Their
# bytes are consequently covered by the normal product/receipt hashes and do
# not depend on a test fixture or an extra installed runtime file.
find_program(SIRIUS_SPIRV_AS spirv-as)
find_program(SIRIUS_SPIRV_DIS spirv-dis)
if(SIRIUS_SPIRV_VAL AND SIRIUS_SPIRV_AS AND SIRIUS_SPIRV_DIS)
    set(retained_build "${CMAKE_CURRENT_BINARY_DIR}/retained")
    set(retained_header "${retained_build}/retained_kernels.h")
    add_custom_command(
        OUTPUT "${retained_header}"
        BYPRODUCTS "${retained_build}/retained_camera.spv"
                   "${retained_build}/retained_transport.spv"
                   "${retained_build}/retained_endpoint.spv"
                   "${retained_build}/retained_dense.spv"
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/build-retained-kernels.py"
            --output "${retained_header}" --compiler "${SIRIUS_SLANGC}"
            --assembler "${SIRIUS_SPIRV_AS}" --disassembler "${SIRIUS_SPIRV_DIS}"
            --validator "${SIRIUS_SPIRV_VAL}"
        DEPENDS "${CMAKE_SOURCE_DIR}/scripts/build-retained-kernels.py"
                "${SIRIUS_KERNEL_DIR}/retained_program.py"
                "${SIRIUS_KERNEL_DIR}/retained_pair.slang"
                "${SIRIUS_KERNEL_DIR}/retained_triple.slang"
                "${SIRIUS_KERNEL_DIR}/retained_camera.slang"
                "${SIRIUS_KERNEL_DIR}/retained_transport.slang"
                "${SIRIUS_KERNEL_DIR}/retained_endpoint.slang"
                "${SIRIUS_KERNEL_DIR}/retained_dense.slang"
        COMMENT "Building bounded retained camera and joint transport stages"
        VERBATIM)
    target_sources(sirius_backend PRIVATE retained_compute.cpp "${retained_header}")
    target_include_directories(sirius_backend PRIVATE "${retained_build}")
    target_compile_definitions(sirius_backend PUBLIC SIRIUS_HAS_RETAINED_COMPUTE=1)
endif()

set(SIRIUS_KERNEL_OUTPUTS "")
foreach(kernel ${SIRIUS_KERNEL_SOURCES})
    set(kernel_src "${SIRIUS_KERNEL_DIR}/${kernel}.slang")
    set(kernel_out "${SIRIUS_KERNEL_BINARY_DIR}/${kernel}.spv")
    if(SIRIUS_SPIRV_VAL)
        set(kernel_validate COMMAND "${SIRIUS_SPIRV_VAL}" "${kernel_out}")
    else()
        set(kernel_validate "")
    endif()
    add_custom_command(
        OUTPUT "${kernel_out}"
        COMMAND ${SIRIUS_SLANGC} "${kernel_src}"
            -target spirv -profile spirv_1_5
            -entry ComputeMain -stage compute
            -I "${SIRIUS_KERNEL_DIR}"
            -o "${kernel_out}"
        ${kernel_validate}
        DEPENDS "${kernel_src}" ${SIRIUS_KERNEL_LIBS}
        COMMENT "Compiling Slang kernel ${kernel}.slang -> SPIR-V (+spirv-val)"
        VERBATIM
    )
    list(APPEND SIRIUS_KERNEL_OUTPUTS "${kernel_out}")
endforeach()

# Precision probes exercise the same Cartesian trajectory and adaptive-step
# authorities as trace.slang without producing an image.  They are test-only
# artefacts. Real is selected session-wide; the original parity probe retains
# float I/O while the dedicated probes below preserve Real-width records.
if(BUILD_TESTS)
    # These probes preserve Real-width inputs and outputs at the test boundary.
    foreach(real_probe IN ITEMS infinity_probe metric_consistency_probe camera_frame_probe)
        foreach(probe_rung IN ITEMS fp32 fp32comp fp64)
            set(probe_flags "")
            set(probe_suffix "")
            if(probe_rung STREQUAL "fp32comp")
                set(probe_flags -DSIRIUS_FP32_COMP)
                set(probe_suffix "_fp32comp")
            elseif(probe_rung STREQUAL "fp64")
                set(probe_flags -DSIRIUS_FP64 -warnings-disable 30082)
                set(probe_suffix "_fp64")
            endif()
            set(probe_out "${SIRIUS_KERNEL_BINARY_DIR}/${real_probe}${probe_suffix}.spv")
            if(SIRIUS_SPIRV_VAL)
                set(probe_validate COMMAND "${SIRIUS_SPIRV_VAL}" "${probe_out}")
            else()
                set(probe_validate "")
            endif()
            add_custom_command(
                OUTPUT "${probe_out}"
                COMMAND ${SIRIUS_SLANGC} "${SIRIUS_KERNEL_DIR}/${real_probe}.slang"
                    ${probe_flags}
                    -target spirv -profile spirv_1_5
                    -entry ComputeMain -stage compute
                    -I "${SIRIUS_KERNEL_DIR}"
                    -o "${probe_out}"
                ${probe_validate}
                DEPENDS "${SIRIUS_KERNEL_DIR}/${real_probe}.slang" ${SIRIUS_KERNEL_LIBS}
                COMMENT "Compiling ${real_probe}.slang (${probe_rung}) -> SPIR-V (+spirv-val)"
                VERBATIM
            )
            list(APPEND SIRIUS_KERNEL_OUTPUTS "${probe_out}")
        endforeach()
    endforeach()
    # The existing coupled regression consumes only fp64 and retains float I/O.
    # Isolate its large call graph from the ordinary three-rung parity probe.
    set(coupled_probe_fp64_out "${SIRIUS_KERNEL_BINARY_DIR}/coupled_probe_fp64.spv")
    if(SIRIUS_SPIRV_VAL)
        set(coupled_probe_fp64_validate COMMAND "${SIRIUS_SPIRV_VAL}" "${coupled_probe_fp64_out}")
    else()
        set(coupled_probe_fp64_validate "")
    endif()
    add_custom_command(
        OUTPUT "${coupled_probe_fp64_out}"
        COMMAND ${SIRIUS_SLANGC} "${SIRIUS_KERNEL_DIR}/coupled_probe.slang"
            -DSIRIUS_FP64
            -target spirv -profile spirv_1_5
            -entry ComputeMain -stage compute
            -I "${SIRIUS_KERNEL_DIR}"
            -warnings-disable 30082
            -o "${coupled_probe_fp64_out}"
        ${coupled_probe_fp64_validate}
        DEPENDS "${SIRIUS_KERNEL_DIR}/coupled_probe.slang" ${SIRIUS_KERNEL_LIBS}
        COMMENT "Compiling coupled_probe.slang (fp64) -> SPIR-V (+spirv-val)"
        VERBATIM
    )
    list(APPEND SIRIUS_KERNEL_OUTPUTS "${coupled_probe_fp64_out}")

    set(kernel_probe_fp64_out "${SIRIUS_KERNEL_BINARY_DIR}/parity_probe_fp64.spv")
    if(SIRIUS_SPIRV_VAL)
        set(kernel_probe_fp64_validate COMMAND "${SIRIUS_SPIRV_VAL}" "${kernel_probe_fp64_out}")
    else()
        set(kernel_probe_fp64_validate "")
    endif()
    add_custom_command(
        OUTPUT "${kernel_probe_fp64_out}"
        COMMAND ${SIRIUS_SLANGC} "${SIRIUS_KERNEL_DIR}/parity_probe.slang"
            -DSIRIUS_FP64
            -target spirv -profile spirv_1_5
            -entry ComputeMain -stage compute
            -I "${SIRIUS_KERNEL_DIR}"
            -warnings-disable 30082
            -o "${kernel_probe_fp64_out}"
        ${kernel_probe_fp64_validate}
        DEPENDS "${SIRIUS_KERNEL_DIR}/parity_probe.slang" ${SIRIUS_KERNEL_LIBS}
        COMMENT "Compiling parity_probe.slang (fp64 trajectory probe) -> SPIR-V (+spirv-val)"
        VERBATIM
    )
    list(APPEND SIRIUS_KERNEL_OUTPUTS "${kernel_probe_fp64_out}")

    set(kernel_probe_comp_out "${SIRIUS_KERNEL_BINARY_DIR}/parity_probe_fp32comp.spv")
    if(SIRIUS_SPIRV_VAL)
        set(kernel_probe_comp_validate COMMAND "${SIRIUS_SPIRV_VAL}" "${kernel_probe_comp_out}")
    else()
        set(kernel_probe_comp_validate "")
    endif()
    add_custom_command(
        OUTPUT "${kernel_probe_comp_out}"
        COMMAND ${SIRIUS_SLANGC} "${SIRIUS_KERNEL_DIR}/parity_probe.slang"
            -DSIRIUS_FP32_COMP
            -target spirv -profile spirv_1_5
            -entry ComputeMain -stage compute
            -I "${SIRIUS_KERNEL_DIR}"
            -o "${kernel_probe_comp_out}"
        ${kernel_probe_comp_validate}
        DEPENDS "${SIRIUS_KERNEL_DIR}/parity_probe.slang" ${SIRIUS_KERNEL_LIBS}
        COMMENT "Compiling parity_probe.slang (compensated fp32 trajectory probe) -> SPIR-V (+spirv-val)"
        VERBATIM
    )
    list(APPEND SIRIUS_KERNEL_OUTPUTS "${kernel_probe_comp_out}")
endif()

# The fp64 rung: the same trace source compiled with SIRIUS_FP64, which widens
# the Cartesian trajectory core (gr_types Real) to double. One source, two
# artefacts; the renderer selects trace_fp64.spv when SIRIUS_PRECISION=fp64 and
# the device reports shaderFloat64. Warning 30082 (implicit float-to-double)
# is disabled for this compile only: the float->Real widenings at parameter
# and literal boundaries are the design, not an accident.
set(kernel_fp64_out "${SIRIUS_KERNEL_BINARY_DIR}/trace_fp64.spv")
if(SIRIUS_SPIRV_VAL)
    set(kernel_fp64_validate COMMAND "${SIRIUS_SPIRV_VAL}" "${kernel_fp64_out}")
else()
    set(kernel_fp64_validate "")
endif()
add_custom_command(
    OUTPUT "${kernel_fp64_out}"
    COMMAND ${SIRIUS_SLANGC} "${SIRIUS_KERNEL_DIR}/trace.slang"
        -DSIRIUS_FP64
        -target spirv -profile spirv_1_5
        -entry ComputeMain -stage compute
        -I "${SIRIUS_KERNEL_DIR}"
        -warnings-disable 30082
        -o "${kernel_fp64_out}"
    ${kernel_fp64_validate}
    DEPENDS "${SIRIUS_KERNEL_DIR}/trace.slang" ${SIRIUS_KERNEL_LIBS}
    COMMENT "Compiling Slang kernel trace.slang (fp64 rung) -> SPIR-V (+spirv-val)"
    VERBATIM
)
list(APPEND SIRIUS_KERNEL_OUTPUTS "${kernel_fp64_out}")

# The compensated rung (ladder middle): plain fp32 arithmetic with Kahan
# compensation on the per-step state accumulation (SIRIUS_FP32_COMP). Runs on
# any device; SIRIUS_PRECISION=fp32-comp selects it.
set(kernel_comp_out "${SIRIUS_KERNEL_BINARY_DIR}/trace_fp32comp.spv")
if(SIRIUS_SPIRV_VAL)
    set(kernel_comp_validate COMMAND "${SIRIUS_SPIRV_VAL}" "${kernel_comp_out}")
else()
    set(kernel_comp_validate "")
endif()
add_custom_command(
    OUTPUT "${kernel_comp_out}"
    COMMAND ${SIRIUS_SLANGC} "${SIRIUS_KERNEL_DIR}/trace.slang"
        -DSIRIUS_FP32_COMP
        -target spirv -profile spirv_1_5
        -entry ComputeMain -stage compute
        -I "${SIRIUS_KERNEL_DIR}"
        -o "${kernel_comp_out}"
    ${kernel_comp_validate}
    DEPENDS "${SIRIUS_KERNEL_DIR}/trace.slang" ${SIRIUS_KERNEL_LIBS}
    COMMENT "Compiling Slang kernel trace.slang (compensated fp32 rung) -> SPIR-V (+spirv-val)"
    VERBATIM
)
list(APPEND SIRIUS_KERNEL_OUTPUTS "${kernel_comp_out}")

# --- Portability emissions ----------------------------------------------------
# Hardware agnosticism is a standing gate, not a one-off claim: the same trace
# source that ships as SPIR-V also emits vendor-native kernel source for the
# two ecosystems Vulkan compute does not reach natively — CUDA C++ for the
# nvcc toolchain and Metal Shading Language for Apple silicon. No adapter
# consumes these yet (recorded deferral, docs/ARCHITECTURE.md section 4); the
# build proves single-source portability and the KernelPortability suite pins
# each artefact's entry point, so a kernel edit that breaks a vendor target
# fails here rather than on first contact with that vendor's hardware.
set(SIRIUS_PORTABILITY_DIR "${SIRIUS_KERNEL_BINARY_DIR}/portability")
file(MAKE_DIRECTORY "${SIRIUS_PORTABILITY_DIR}")

set(kernel_cuda_out "${SIRIUS_PORTABILITY_DIR}/trace.cu")
add_custom_command(
    OUTPUT "${kernel_cuda_out}"
    COMMAND ${SIRIUS_SLANGC} "${SIRIUS_KERNEL_DIR}/trace.slang"
        -target cuda
        -entry ComputeMain -stage compute
        -I "${SIRIUS_KERNEL_DIR}"
        -o "${kernel_cuda_out}"
    DEPENDS "${SIRIUS_KERNEL_DIR}/trace.slang" ${SIRIUS_KERNEL_LIBS}
    COMMENT "Emitting Slang kernel trace.slang -> CUDA C++ (portability gate)"
    VERBATIM
)
list(APPEND SIRIUS_KERNEL_OUTPUTS "${kernel_cuda_out}")

set(kernel_metal_out "${SIRIUS_PORTABILITY_DIR}/trace.metal")
add_custom_command(
    OUTPUT "${kernel_metal_out}"
    COMMAND ${SIRIUS_SLANGC} "${SIRIUS_KERNEL_DIR}/trace.slang"
        -target metal
        -entry ComputeMain -stage compute
        -I "${SIRIUS_KERNEL_DIR}"
        -o "${kernel_metal_out}"
    DEPENDS "${SIRIUS_KERNEL_DIR}/trace.slang" ${SIRIUS_KERNEL_LIBS}
    COMMENT "Emitting Slang kernel trace.slang -> Metal Shading Language (portability gate)"
    VERBATIM
)
list(APPEND SIRIUS_KERNEL_OUTPUTS "${kernel_metal_out}")

add_custom_target(sirius_kernels ALL DEPENDS ${SIRIUS_KERNEL_OUTPUTS})
add_dependencies(sirius_backend sirius_kernels)
