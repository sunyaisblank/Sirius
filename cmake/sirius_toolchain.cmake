# Toolchain gate and C++26 target configuration for the new tree.
#
# CMake 3.28 maps cxx_std_26 for Clang but not for GCC 14, so the standard
# flag is passed explicitly and identically on both (measured 2026-07-17,
# see docs/SPECIFICATION.md section 1.7). MSVC takes /std:c++latest.

set(SIRIUS_MODERN_AVAILABLE OFF)
set(SIRIUS_MODERN_UNAVAILABLE_REASON "")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 14)
        set(SIRIUS_MODERN_AVAILABLE ON)
    else()
        set(SIRIUS_MODERN_UNAVAILABLE_REASON
            "GCC ${CMAKE_CXX_COMPILER_VERSION} lacks C++26 mode; GCC 14 or newer required")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 17)
        set(SIRIUS_MODERN_AVAILABLE ON)
    else()
        set(SIRIUS_MODERN_UNAVAILABLE_REASON
            "Clang ${CMAKE_CXX_COMPILER_VERSION} lacks C++26 mode; Clang 17 or newer required")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 16)
        set(SIRIUS_MODERN_AVAILABLE ON)
    else()
        set(SIRIUS_MODERN_UNAVAILABLE_REASON
            "AppleClang ${CMAKE_CXX_COMPILER_VERSION} lacks C++26 mode; Xcode 16 or newer required")
    endif()
elseif(MSVC)
    if(MSVC_VERSION GREATER_EQUAL 1940)
        set(SIRIUS_MODERN_AVAILABLE ON)
    else()
        set(SIRIUS_MODERN_UNAVAILABLE_REASON
            "MSVC ${MSVC_VERSION} is older than 19.40; Visual Studio 2022 17.10 or newer required")
    endif()
else()
    set(SIRIUS_MODERN_UNAVAILABLE_REASON
        "unrecognised compiler '${CMAKE_CXX_COMPILER_ID}'")
endif()

option(SIRIUS_WERROR "Treat warnings as errors in the new tree" ON)

# Applies the standard mode, strict warnings, and contract semantics to a
# new-tree target. Contract mode: 2 enforce, 1 observe, 0 ignore; the default
# (enforce in Debug, observe otherwise) is in base/contracts.h, and targets
# that need a different mode say so here rather than in source.
function(sirius_configure_target target)
    cmake_parse_arguments(ARG "" "CONTRACT_MODE" "" ${ARGN})

    if(MSVC)
        target_compile_options(${target} PRIVATE /std:c++latest /W4 /permissive- /utf-8)
        if(SIRIUS_WERROR)
            target_compile_options(${target} PRIVATE /WX)
        endif()
    else()
        target_compile_options(${target} PRIVATE -std=c++2c -Wall -Wextra -Wpedantic)
        if(SIRIUS_WERROR)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
    endif()

    if(DEFINED ARG_CONTRACT_MODE)
        target_compile_definitions(${target} PRIVATE SIRIUS_CONTRACT_MODE=${ARG_CONTRACT_MODE})
    endif()
endfunction()
