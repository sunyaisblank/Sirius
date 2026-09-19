# CPU transport remains testable without Vulkan development files or Slang.
# Device tests extend the same executable when their backend is available.
add_executable(sirius_backend_tests
    cpu_jacobi_oracle_test.cpp
    cpu_trace_boundary_test.cpp
    coupled_transport_test.cpp
    source_sky_map_trace_test.cpp
    source_sky_infinity_trace_test.cpp
)

target_link_libraries(sirius_backend_tests PRIVATE sirius_backend_cpu GTest::gtest_main)
target_include_directories(sirius_backend_tests PRIVATE "${CMAKE_SOURCE_DIR}/tests")
sirius_configure_target(sirius_backend_tests CONTRACT_MODE 2)

gtest_discover_tests(
    sirius_backend_tests
    PROPERTIES
        # Keep one device lock for this mixed executable on Vulkan builds.
        RESOURCE_LOCK sirius_vulkan_device
        LABELS "Correctness")

if(NOT TARGET sirius_backend)
    message(STATUS "Sirius tests: CPU backend suite enabled; Vulkan tests unavailable")
    return()
endif()

target_sources(sirius_backend_tests PRIVATE
    vulkan_smoke_test.cpp
    kernel_parity_test.cpp
    kernel_metric_consistency_test.cpp
    kernel_camera_frame_test.cpp
    retained_camera_program_test.cpp
    retained_compute_test.cpp
    kernel_infinity_test.cpp
    kernel_trace_test.cpp
    kernel_beam_test.cpp
    kernel_portability_test.cpp
)

# Direct trace/beam probes consume the production per-submission governor.
# This is a test executable dependency; the backend library remains independent
# of render orchestration. Actual renderer/session tests belong to tests/render.
target_link_libraries(sirius_backend_tests PRIVATE
    sirius_backend sirius_render)

if(TARGET sirius_kernels)
    target_compile_definitions(sirius_backend_tests
        PRIVATE SIRIUS_KERNEL_DIR="${SIRIUS_KERNEL_BINARY_DIR}")
    add_dependencies(sirius_backend_tests sirius_kernels)

    find_program(SIRIUS_SPIRV_AS spirv-as)
    find_program(SIRIUS_SPIRV_DIS spirv-dis)
    if(SIRIUS_SPIRV_AS AND SIRIUS_SPIRV_DIS AND SIRIUS_SPIRV_VAL)
        set(retained_source "${CMAKE_SOURCE_DIR}/tests/support/retained_camera")
        set(retained_build "${CMAKE_CURRENT_BINARY_DIR}/retained_camera")
        set(retained_products
            "${retained_build}/program_fixture.h"
            "${retained_build}/program_camera_probe-fp32.spv"
            "${retained_build}/program_camera_probe-fp32comp.spv")
        add_custom_command(
            OUTPUT ${retained_products}
            COMMAND "${Python3_EXECUTABLE}" "${retained_source}/build_fixture.py"
                --output "${retained_build}/program_fixture.h"
            COMMAND "${Python3_EXECUTABLE}" "${retained_source}/check_shaders.py"
                --probe program_camera_probe --output "${retained_build}"
                --compiler "${SIRIUS_SLANGC}"
            DEPENDS
                "${retained_source}/build_fixture.py"
                "${retained_source}/program.py"
                "${CMAKE_SOURCE_DIR}/src/sirius/kernels/retained_program.py"
                "${retained_source}/reference_cases.json"
                "${retained_source}/check_shaders.py"
                "${retained_source}/program_camera_probe.slang"
                "${retained_source}/retained_pair.slang"
            COMMENT "Building bounded retained camera program and independent witnesses"
            VERBATIM)
        add_custom_target(sirius_retained_camera_test_inputs DEPENDS ${retained_products})
        set_property(TARGET sirius_retained_camera_test_inputs PROPERTY
            SIRIUS_TEST_INPUT_ARTIFACTS
            "retained_camera_fixture=${retained_build}/program_fixture.h"
            "retained_camera_fp32_spv=${retained_build}/program_camera_probe-fp32.spv"
            "retained_camera_fp32comp_spv=${retained_build}/program_camera_probe-fp32comp.spv")
        add_dependencies(sirius_backend_tests sirius_retained_camera_test_inputs)
        target_include_directories(sirius_backend_tests PRIVATE "${retained_build}")
        target_compile_definitions(sirius_backend_tests PRIVATE
            SIRIUS_RETAINED_CAMERA_TEST_DIR="${retained_build}")
    endif()
endif()
