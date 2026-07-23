# =============================================================================
# Sirius: general-relativistic renderer.
# Governance: docs/SPECIFICATION.md, docs/ARCHITECTURE.md, docs/STYLE.md.
# =============================================================================
cmake_minimum_required(VERSION 3.28)
project(Sirius VERSION 2.0.0 LANGUAGES CXX C)

option(BUILD_TESTS "Build test suites" ON)
option(SIRIUS_MANDATORY_TESTS "Require mandatory tests to pass for build" OFF)
option(SIRIUS_SKIP_MANDATORY "Skip mandatory tests (development override)" OFF)

include(cmake/sirius_toolchain.cmake)
include(cmake/sirius_dependencies.cmake)

find_package(Threads REQUIRED)

enable_testing()

if(NOT SIRIUS_MODERN_AVAILABLE)
    message(FATAL_ERROR
        "Sirius requires a C++26-mode compiler: ${SIRIUS_MODERN_UNAVAILABLE_REASON}. "
        "Use a preset from CMakePresets.json to select a supported compiler.")
endif()

add_subdirectory(src/sirius)
if(BUILD_TESTS)
    FetchContent_MakeAvailable(googletest)
    add_subdirectory(tests)
endif()

# === Mandatory test gating (the build gate; label authority is generated) ===
if(BUILD_TESTS AND SIRIUS_MANDATORY_TESTS AND NOT SIRIUS_SKIP_MANDATORY)
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

    add_custom_target(RunMandatoryTests ALL
        COMMAND ${CMAKE_CTEST_COMMAND}
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
endif()

set(CPACK_PACKAGE_NAME "sirius")
set(CPACK_PACKAGE_VENDOR "Sirius")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "General-relativistic ray tracing engine")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
include(CPack)
