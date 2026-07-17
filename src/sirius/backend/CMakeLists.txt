# Backend module: device seam plus the Vulkan compute adapter, and the
# Slang kernel build (single source authority in ../kernels compiled to
# SPIR-V here). Both degrade with a clear message rather than failing the
# configure: machines without Vulkan or slangc build the CPU system only.

find_package(Vulkan QUIET)
if(NOT Vulkan_FOUND)
    message(STATUS "Sirius backend: Vulkan development files not found; Vulkan backend skipped")
    return()
endif()

add_library(sirius_backend STATIC
    vulkan/vulkan_device.cpp
)

target_link_libraries(sirius_backend PUBLIC sirius_base Vulkan::Vulkan)

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

set(SIRIUS_KERNEL_SOURCES
    smoke
)

set(SIRIUS_KERNEL_OUTPUTS "")
foreach(kernel ${SIRIUS_KERNEL_SOURCES})
    set(kernel_src "${CMAKE_CURRENT_SOURCE_DIR}/../kernels/${kernel}.slang")
    set(kernel_out "${SIRIUS_KERNEL_BINARY_DIR}/${kernel}.spv")
    add_custom_command(
        OUTPUT "${kernel_out}"
        COMMAND ${SIRIUS_SLANGC} "${kernel_src}"
            -target spirv -profile spirv_1_5
            -entry ComputeMain -stage compute
            -o "${kernel_out}"
        DEPENDS "${kernel_src}"
        COMMENT "Compiling Slang kernel ${kernel}.slang -> SPIR-V"
        VERBATIM
    )
    list(APPEND SIRIUS_KERNEL_OUTPUTS "${kernel_out}")
endforeach()

add_custom_target(sirius_kernels ALL DEPENDS ${SIRIUS_KERNEL_OUTPUTS})
add_dependencies(sirius_backend sirius_kernels)
