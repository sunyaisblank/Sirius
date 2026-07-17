# =============================================================================
# Sirius: general-relativistic renderer.
# Root build: new C++26 tree (src/sirius, tests) plus the verified legacy
# C++17 system (cmake/legacy_build.cmake) during the port. Governance:
# docs/SPECIFICATION.md, docs/ARCHITECTURE.md, docs/STYLE.md.
# =============================================================================
cmake_minimum_required(VERSION 3.28)
project(Sirius VERSION 2.0.0 LANGUAGES CXX C)

option(BUILD_TESTS "Build test suites" ON)
option(SIRIUS_BUILD_LEGACY "Build the pre-rebuild C++17 system during the port" ON)
option(SIRIUS_MANDATORY_TESTS "Require mandatory tests to pass for build" OFF)
option(SIRIUS_SKIP_MANDATORY "Skip mandatory tests (development override)" OFF)

include(cmake/sirius_toolchain.cmake)
include(cmake/sirius_dependencies.cmake)

find_package(Threads REQUIRED)

enable_testing()

if(SIRIUS_BUILD_LEGACY)
    include(cmake/legacy_build.cmake)
endif()

if(SIRIUS_MODERN_AVAILABLE)
    add_subdirectory(src/sirius)
    if(BUILD_TESTS)
        FetchContent_MakeAvailable(googletest)
        add_subdirectory(tests)
    endif()
else()
    message(WARNING
        "New tree skipped: ${SIRIUS_MODERN_UNAVAILABLE_REASON}. "
        "Use a preset from CMakePresets.json to select a supported compiler.")
endif()

# === Install and packaging (ledger item: no install target existed) ===
if(TARGET Sirius)
    install(TARGETS Sirius RUNTIME DESTINATION bin)
    install(TARGETS SiriusCore
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib)
endif()

set(CPACK_PACKAGE_NAME "sirius")
set(CPACK_PACKAGE_VENDOR "Sirius")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "General-relativistic ray tracing engine")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
include(CPack)
