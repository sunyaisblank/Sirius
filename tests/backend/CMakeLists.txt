# Backend suite: runs against whatever Vulkan device is present (Lavapipe
# in CI and WSL2), and skips cleanly when none is.

if(NOT TARGET sirius_backend)
    message(STATUS "Sirius tests: backend module absent; backend suite skipped")
    return()
endif()

add_executable(sirius_backend_tests
    vulkan_smoke_test.cpp
    kernel_parity_test.cpp
    kernel_trace_test.cpp
    kernel_beam_test.cpp
    kernel_portability_test.cpp
)

target_link_libraries(sirius_backend_tests PRIVATE sirius_backend GTest::gtest_main)

sirius_configure_target(sirius_backend_tests CONTRACT_MODE 2)

if(TARGET sirius_kernels)
    target_compile_definitions(sirius_backend_tests
        PRIVATE SIRIUS_KERNEL_DIR="${SIRIUS_KERNEL_BINARY_DIR}")
    add_dependencies(sirius_backend_tests sirius_kernels)
endif()

gtest_discover_tests(sirius_backend_tests PROPERTIES LABELS "Correctness")
