# =============================================================================
# Sirius: general-relativistic renderer.
# Governance: docs/SPECIFICATION.md, docs/ARCHITECTURE.md, docs/STYLE.md.
# =============================================================================
cmake_minimum_required(VERSION 3.28)
project(Sirius VERSION 2.0.0 LANGUAGES CXX C)

option(BUILD_TESTS "Build test suites" ON)
option(SIRIUS_MANDATORY_TESTS "Require mandatory tests to pass for every normal build" ON)
option(SIRIUS_REQUIRE_VULKAN_RUNTIME
    "Require compiled Vulkan kernels and a dispatchable Vulkan device in the operational gate"
    OFF)

if(SIRIUS_MANDATORY_TESTS AND NOT BUILD_TESTS)
    message(FATAL_ERROR
        "SIRIUS_MANDATORY_TESTS=ON requires BUILD_TESTS=ON. "
        "An operational build cannot omit its executable invariants.")
endif()

include(cmake/sirius_toolchain.cmake)
include(cmake/sirius_dependencies.cmake)

find_package(Threads REQUIRED)
if(BUILD_TESTS)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
endif()

enable_testing()

if(NOT SIRIUS_MODERN_AVAILABLE)
    message(FATAL_ERROR
        "Sirius requires a C++26-mode compiler: ${SIRIUS_MODERN_UNAVAILABLE_REASON}. "
        "Use a preset from CMakePresets.json to select a supported compiler.")
endif()

add_subdirectory(src/sirius)
if(SIRIUS_REQUIRE_VULKAN_RUNTIME AND
   (NOT TARGET sirius_backend OR NOT TARGET sirius_kernels))
    message(FATAL_ERROR
        "SIRIUS_REQUIRE_VULKAN_RUNTIME=ON requires the Vulkan SDK and slangc. "
        "The required runtime profile may not degrade to a CPU-only build.")
endif()
if(BUILD_TESTS)
    FetchContent_MakeAvailable(googletest)
    add_subdirectory(tests)
    add_custom_target(SiriusSourceGovernance ALL
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/generate-ctest-labels.py"
            --check
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/verify-operating-model.py"
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        COMMENT "Verifying test-label and operating-model governance")
endif()

# === Mandatory test gating (the build gate; label authority is generated) ===
if(BUILD_TESTS AND SIRIUS_MANDATORY_TESTS)
    set(sirius_mandatory_test_targets)
    foreach(test_target
            sirius_base_tests
            sirius_core_tests
            sirius_oracle_tests
            sirius_backend_tests
            sirius_render_tests
            sirius_app_tests)
        if(TARGET ${test_target})
            list(APPEND sirius_mandatory_test_targets ${test_target})
        endif()
    endforeach()

    set(sirius_mandatory_ctest_command "${CMAKE_CTEST_COMMAND}")
    if(NOT SIRIUS_SANITIZERS STREQUAL "none")
        set(sirius_mandatory_ctest_command
            "${CMAKE_COMMAND}" -E env
            "ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:strict_string_checks=1"
            "LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/sanitizers/lsan-vulkan.supp:print_suppressions=1"
            "UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1"
            "${CMAKE_CTEST_COMMAND}")
    endif()

    add_custom_target(RunMandatoryTests ALL
        COMMAND ${sirius_mandatory_ctest_command}
            --test-dir "${CMAKE_BINARY_DIR}"
            --output-on-failure
            --no-tests=error
            -L "Mandatory"
            -C $<CONFIG>
        # gtest_discover_tests writes each suite's CTest include after linking.
        # Depend on every available suite so the gate cannot race discovery and
        # incorrectly report that no Mandatory tests exist.
        DEPENDS sirius ${sirius_mandatory_test_targets}
        COMMENT "=== MANDATORY TEST GATE ==="
    )
    message(STATUS "Mandatory test gating: ENABLED")
endif()

# === Install and packaging ===
if(TARGET sirius)
    install(TARGETS sirius RUNTIME DESTINATION bin)
    install(FILES assets/Starfield.png DESTINATION share/sirius/assets)
    install(FILES tests/operating_model.json DESTINATION share/sirius/model)
    if(SIRIUS_INSTALL_HAS_VIEWER)
        install(DIRECTORY src/sirius/app/viewer/shaders/
            DESTINATION share/sirius/shaders
            FILES_MATCHING PATTERN "*.vert" PATTERN "*.frag")
    endif()
    if(TARGET sirius_kernels)
        install(FILES
            "${SIRIUS_KERNEL_BINARY_DIR}/trace.spv"
            "${SIRIUS_KERNEL_BINARY_DIR}/trace_fp32comp.spv"
            "${SIRIUS_KERNEL_BINARY_DIR}/trace_fp64.spv"
            DESTINATION share/sirius/kernels)
        set(SIRIUS_INSTALL_HAS_VULKAN TRUE)
    else()
        set(SIRIUS_INSTALL_HAS_VULKAN FALSE)
    endif()
    configure_file(
        cmake/verify_sirius_install.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/verify_sirius_install.cmake"
        @ONLY)
    install(SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/verify_sirius_install.cmake")
endif()

set(CPACK_PACKAGE_NAME "sirius")
set(CPACK_PACKAGE_VENDOR "Sirius")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "General-relativistic ray tracing engine")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
include(CPack)
