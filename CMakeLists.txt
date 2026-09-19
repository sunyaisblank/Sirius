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
set(SIRIUS_ALIGNMENT_MODE "development" CACHE STRING
    "Alignment enforcement: development is local; qualification is strict evidence production; release requires all attestations")
set_property(CACHE SIRIUS_ALIGNMENT_MODE PROPERTY STRINGS development qualification release)
set(SIRIUS_ATTESTATION_ROOT "" CACHE PATH
    "Root containing revision-bound external attestation bundles")

if(NOT SIRIUS_ALIGNMENT_MODE STREQUAL "development" AND
   NOT SIRIUS_ALIGNMENT_MODE STREQUAL "qualification" AND
   NOT SIRIUS_ALIGNMENT_MODE STREQUAL "release")
    message(FATAL_ERROR
        "SIRIUS_ALIGNMENT_MODE must be development, qualification, or release, not "
        "'${SIRIUS_ALIGNMENT_MODE}'")
endif()

if(SIRIUS_MANDATORY_TESTS AND NOT BUILD_TESTS)
    message(FATAL_ERROR
        "SIRIUS_MANDATORY_TESTS=ON requires BUILD_TESTS=ON. "
        "An operational build cannot omit its executable invariants.")
endif()

include(cmake/sirius_toolchain.cmake)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

include(cmake/sirius_build_policy.cmake)
sirius_verify_configure_policy()

include(cmake/sirius_dependencies.cmake)
find_package(Threads REQUIRED)

include(cmake/sirius_alignment.cmake)

enable_testing()

if(NOT SIRIUS_MODERN_AVAILABLE)
    message(FATAL_ERROR
        "Sirius requires a C++26-mode compiler: ${SIRIUS_MODERN_UNAVAILABLE_REASON}. "
        "Use a preset from CMakePresets.json to select a supported compiler.")
endif()

add_subdirectory(src/sirius)
sirius_verify_product_policy()

include(cmake/sirius_testing.cmake)

include(cmake/sirius_install.cmake)
