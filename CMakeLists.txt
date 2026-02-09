# =============================================================================
# Sirius CMake Configuration
# Restructured: Jan 2026
# =============================================================================
cmake_minimum_required(VERSION 3.16)
project(Sirius LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build output directories
set(CMAKE_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/bin/Sirius.Build")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# =============================================================================
# CUDA and OptiX Support
# =============================================================================
include(CheckLanguage)
check_language(CUDA)

if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    find_package(CUDAToolkit 12.0 REQUIRED)
    set(CMAKE_CUDA_ARCHITECTURES "75;86;89" CACHE STRING "CUDA architectures")
    message(STATUS "CUDA Toolkit found: ${CUDAToolkit_VERSION}")
    message(STATUS "CUDA Compiler: ${CMAKE_CUDA_COMPILER}")
    message(STATUS "CUDA Architectures: ${CMAKE_CUDA_ARCHITECTURES}")
    set(SIRIUS_HAS_CUDA ON)
else()
    message(WARNING "CUDA compiler not found - OptiX backend will be disabled")
    set(SIRIUS_HAS_CUDA OFF)
endif()

# Find OptiX SDK
set(OptiX_INSTALL_DIR "" CACHE PATH "Path to OptiX SDK installation")

# Detect WSL2 environment
set(IS_WSL2 OFF)
if(EXISTS "/proc/sys/fs/binfmt_misc/WSLInterop" OR EXISTS "/run/WSL")
    set(IS_WSL2 ON)
    message(STATUS "WSL2 environment detected")
endif()

if(NOT EXISTS "${OptiX_INSTALL_DIR}/include/optix.h")
    if(DEFINED ENV{OptiX_INSTALL_DIR} AND EXISTS "$ENV{OptiX_INSTALL_DIR}/include/optix.h")
        set(OptiX_INSTALL_DIR "$ENV{OptiX_INSTALL_DIR}" CACHE PATH "Path to OptiX SDK installation" FORCE)
    # Local project path (WSL2 or portable builds)
    elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib/optix/include/optix.h")
        set(OptiX_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/optix" CACHE PATH "Path to OptiX SDK installation" FORCE)
    # Common Linux installation paths
    elseif(EXISTS "/opt/nvidia/optix/include/optix.h")
        set(OptiX_INSTALL_DIR "/opt/nvidia/optix" CACHE PATH "Path to OptiX SDK installation" FORCE)
    elseif(EXISTS "$ENV{HOME}/NVIDIA-OptiX-SDK/include/optix.h")
        set(OptiX_INSTALL_DIR "$ENV{HOME}/NVIDIA-OptiX-SDK" CACHE PATH "Path to OptiX SDK installation" FORCE)
    # Windows paths (for native Windows builds)
    elseif(EXISTS "C:/ProgramData/NVIDIA Corporation/OptiX SDK 9.1.0/include/optix.h")
        set(OptiX_INSTALL_DIR "C:/ProgramData/NVIDIA Corporation/OptiX SDK 9.1.0" CACHE PATH "Path to OptiX SDK installation" FORCE)
    elseif(EXISTS "C:/ProgramData/NVIDIA Corporation/OptiX SDK 8.1.0/include/optix.h")
        set(OptiX_INSTALL_DIR "C:/ProgramData/NVIDIA Corporation/OptiX SDK 8.1.0" CACHE PATH "Path to OptiX SDK installation" FORCE)
    elseif(EXISTS "C:/ProgramData/NVIDIA Corporation/OptiX SDK 8.0.0/include/optix.h")
        set(OptiX_INSTALL_DIR "C:/ProgramData/NVIDIA Corporation/OptiX SDK 8.0.0" CACHE PATH "Path to OptiX SDK installation" FORCE)
    endif()
endif()

if(SIRIUS_HAS_CUDA AND EXISTS "${OptiX_INSTALL_DIR}/include/optix.h")
    message(STATUS "OptiX SDK found at: ${OptiX_INSTALL_DIR}")
    set(SIRIUS_HAS_OPTIX ON)
    set(OptiX_INCLUDE_DIR "${OptiX_INSTALL_DIR}/include")
else()
    if(SIRIUS_HAS_CUDA)
        message(WARNING "OptiX SDK not found - Set OptiX_INSTALL_DIR to enable OptiX backend")
    endif()
    set(SIRIUS_HAS_OPTIX OFF)
endif()

# --- Platform Detection (unified cross-platform) ---
if(WIN32)
    set(IS_WINDOWS TRUE)
else()
    set(IS_WINDOWS FALSE)
endif()

# === Find System Dependencies ===
find_package(OpenGL REQUIRED)
find_package(Threads REQUIRED)

# === Add Subdirectories for Dependencies ===
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/glfw)

# =============================================================================
# CLI Dependencies (FTXUI, nlohmann_json)
# =============================================================================
include(FetchContent)

# FTXUI for rich terminal output (progress bars, tables, colors)
FetchContent_Declare(ftxui
    GIT_REPOSITORY https://github.com/ArthurSonzogni/FTXUI.git
    GIT_TAG v5.0.0
)
set(FTXUI_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(FTXUI_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(FTXUI_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(ftxui)

# nlohmann_json for configuration files
FetchContent_Declare(json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
)
set(JSON_BuildTests OFF CACHE INTERNAL "")
FetchContent_MakeAvailable(json)

message(STATUS "CLI dependencies: FTXUI v5.0.0, nlohmann_json v3.11.3")

# =============================================================================
# Core Library
# =============================================================================
add_library(SiriusCore SHARED
    # Geodesic integration
    src/Sirius.Core/Geodesic/PHGD001A.cpp
    src/Sirius.Core/Geodesic/PHGD002A.cpp
    src/Sirius.Core/Metric/PHMT000A.cpp
    
    # Tensor algebra
    src/Sirius.Core/Tensor/MTTN001A.cpp
)

target_include_directories(SiriusCore PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Tensor"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Autodiff"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Metric"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Geodesic"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Coordinate"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Spectral"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Disk"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Symplectic"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Transport"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Camera"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Jet"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Polarisation"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Constants"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/glm"
)

target_compile_features(SiriusCore PUBLIC cxx_std_17)

set_target_properties(SiriusCore PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
    WINDOWS_EXPORT_ALL_SYMBOLS ON
)

foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG_UPPER)
    set_target_properties(SiriusCore PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${OUTPUTCONFIG}"
        LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${OUTPUTCONFIG}"
        ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${OUTPUTCONFIG}"
    )
endforeach()

# === Plugins ===
set(PLUGIN_OUTPUT_PATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/plugins")
file(MAKE_DIRECTORY ${PLUGIN_OUTPUT_PATH})

macro(add_metric_plugin PLUGIN_NAME SOURCE_FILE)
    add_library(${PLUGIN_NAME} SHARED ${SOURCE_FILE})
    target_include_directories(${PLUGIN_NAME} PRIVATE 
        "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Metric"
        "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Tensor"
        "${CMAKE_CURRENT_SOURCE_DIR}/lib/glm"
    )
    target_link_libraries(${PLUGIN_NAME} PRIVATE SiriusCore)
    set_target_properties(${PLUGIN_NAME} PROPERTIES 
        LIBRARY_OUTPUT_DIRECTORY "${PLUGIN_OUTPUT_PATH}" 
        RUNTIME_OUTPUT_DIRECTORY "${PLUGIN_OUTPUT_PATH}"
        ARCHIVE_OUTPUT_DIRECTORY "${PLUGIN_OUTPUT_PATH}"
    )
    foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
        string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG_UPPER)
        set_target_properties(${PLUGIN_NAME} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} "${PLUGIN_OUTPUT_PATH}/${OUTPUTCONFIG}"
            LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} "${PLUGIN_OUTPUT_PATH}/${OUTPUTCONFIG}"
            ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} "${PLUGIN_OUTPUT_PATH}/${OUTPUTCONFIG}"
        )
    endforeach()
endmacro()

# =============================================================================
# OptiX Render Backend
# =============================================================================
if(SIRIUS_HAS_OPTIX)
    message(STATUS "Building OptiX render backend...")

    set(OPTIX_PTX_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Sirius.Render/ptx")
    file(MAKE_DIRECTORY ${OPTIX_PTX_DIR})

    set(OPTIX_DEVICE_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Render/Acceleration/OptiX/RDOP002A.cu")

    # Platform-specific PTX compile definitions (cross-platform)
    if(WIN32)
        set(OPTIX_PTX_PLATFORM_DEFS -DNOMINMAX -D_USE_MATH_DEFINES)
    else()
        set(OPTIX_PTX_PLATFORM_DEFS "")
    endif()

    # ==========================================================================
    # Multi-Architecture PTX Pre-Compilation
    # Generates PTX for: Turing (sm_75), Ampere (sm_80, sm_86), Ada (sm_89)
    # ==========================================================================

    # Option to control which architectures to build
    option(SIRIUS_PTX_ALL_ARCHS "Build PTX for all supported GPU architectures" ON)

    if(SIRIUS_PTX_ALL_ARCHS)
        set(SIRIUS_PTX_ARCHITECTURES 75 80 86 89)
    else()
        # Default to single architecture for faster builds
        set(SIRIUS_PTX_ARCHITECTURES 75)
    endif()

    set(OPTIX_PTX_OUTPUTS "")

    foreach(SM_ARCH ${SIRIUS_PTX_ARCHITECTURES})
        set(PTX_OUTPUT "${OPTIX_PTX_DIR}/RDOP002A_sm${SM_ARCH}.ptx")
        list(APPEND OPTIX_PTX_OUTPUTS ${PTX_OUTPUT})

        add_custom_command(
            OUTPUT ${PTX_OUTPUT}
            COMMAND ${CMAKE_CUDA_COMPILER}
                -ptx
                -I "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Render/Acceleration/OptiX"
                -I "${OptiX_INCLUDE_DIR}"
                -I "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}"
                --use_fast_math
                -lineinfo
                --expt-relaxed-constexpr
                --maxrregcount=96
                -std=c++17
                --machine=64
                -arch compute_${SM_ARCH}
                -DSIRIUS_HAS_OPTIX=1
                ${OPTIX_PTX_PLATFORM_DEFS}
                ${OPTIX_DEVICE_SOURCE}
                -o ${PTX_OUTPUT}
            DEPENDS ${OPTIX_DEVICE_SOURCE}
                "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Render/Acceleration/OptiX/RDOP003A.h"
            COMMENT "Compiling OptiX PTX for SM ${SM_ARCH}"
            VERBATIM
        )
    endforeach()

    # Also generate a generic PTX (compute_75) as fallback
    set(OPTIX_PTX_GENERIC "${OPTIX_PTX_DIR}/RDOP002A.ptx")
    add_custom_command(
        OUTPUT ${OPTIX_PTX_GENERIC}
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${OPTIX_PTX_DIR}/RDOP002A_sm75.ptx"
            ${OPTIX_PTX_GENERIC}
        DEPENDS "${OPTIX_PTX_DIR}/RDOP002A_sm75.ptx"
        COMMENT "Creating generic PTX fallback"
        VERBATIM
    )
    list(APPEND OPTIX_PTX_OUTPUTS ${OPTIX_PTX_GENERIC})

    add_custom_target(SiriusOptiXPTX ALL DEPENDS ${OPTIX_PTX_OUTPUTS})

    message(STATUS "PTX architectures: ${SIRIUS_PTX_ARCHITECTURES}")
    
    add_library(SiriusOptiX STATIC
        src/Sirius.Render/Acceleration/OptiX/RDOP001A.cu
    )
    
    add_dependencies(SiriusOptiX SiriusOptiXPTX)
    
    target_include_directories(SiriusOptiX PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Render/Acceleration/OptiX"
        "${OptiX_INCLUDE_DIR}"
        "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}"
    )
    
    set_target_properties(SiriusOptiX PROPERTIES
        CUDA_SEPARABLE_COMPILATION OFF
        POSITION_INDEPENDENT_CODE ON
    )
    
    target_compile_options(SiriusOptiX PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:
            --use_fast_math
            -lineinfo
            --expt-relaxed-constexpr
            -Xcudafe --diag_suppress=20012
        >
    )
    
    target_link_libraries(SiriusOptiX PUBLIC
        CUDA::cudart
        CUDA::cuda_driver
    )
    
    target_compile_definitions(SiriusOptiX PUBLIC 
        SIRIUS_HAS_OPTIX=1
        OPTIX_PTX_DIR="${OPTIX_PTX_DIR}"
    )
    
    set(SHADER_OUTPUT_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Sirius.Render")
    file(MAKE_DIRECTORY ${SHADER_OUTPUT_DIR})
    
    add_custom_command(TARGET SiriusOptiXPTX POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Render/Shader/RDSD003A.frag"
            "${SHADER_OUTPUT_DIR}/RDSD003A.frag"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Render/Shader/RDSD003A.vert"
            "${SHADER_OUTPUT_DIR}/RDSD003A.vert"
        COMMENT "Copying shader files to build directory"
    )
    
    message(STATUS "OptiX backend enabled with CUDA ${CUDAToolkit_VERSION}")
    message(STATUS "OptiX PTX output: ${OPTIX_PTX_DIR}")
endif()

# =============================================================================
# Sirius Executable
# =============================================================================
add_executable(Sirius
    # Infrastructure - Application (entry point)
    src/Sirius.Infrastructure/Application/CREP001A.cpp

    # Infrastructure - CLI (command-line interface)
    src/Sirius.Infrastructure/Cli/CRCL001A.cpp
    src/Sirius.Infrastructure/Cli/CRCL002A.cpp
    src/Sirius.Infrastructure/Cli/CRCL003A.cpp
    src/Sirius.Infrastructure/Cli/CRCL004A.cpp
    src/Sirius.Infrastructure/Cli/CRCL005A.cpp
    src/Sirius.Infrastructure/Cli/CRCL006A.cpp

    # Infrastructure - Configuration (JSON config loader)
    src/Sirius.Infrastructure/Configuration/CRCF002A.cpp

    # Infrastructure - Platform (cross-platform paths)
    src/Sirius.Infrastructure/Platform/CRPF001A.cpp

    # Infrastructure - Session (FSM-based render pipeline)
    src/Sirius.Infrastructure/Session/SNRS001A.cpp

    # Infrastructure - Interface (ImGui)
    src/Sirius.Infrastructure/Interface/UIMN001A.cpp

    # Infrastructure - Window (GLFW)
    src/Sirius.Infrastructure/Window/CRWN001A.cpp

    # Infrastructure - Plugin (metric management)
    src/Sirius.Infrastructure/Plugin/CRPM001A.cpp

    # Core (camera)
    src/Sirius.Core/Camera/CRGC001A.cpp

    # Render pipeline
    src/Sirius.Render/Pipeline/RDRT001A.cpp
    src/Sirius.Render/Pipeline/RDST001A.cpp

    # Render integration (geodesic tracer)
    src/Sirius.Render/Integration/GTRC001A.cpp

    # Render acceleration backend (GPU)
    src/Sirius.Render/Acceleration/Backend/ACBF001A.cpp

    # External dependencies
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/glad/src/glad.c
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/imgui/imgui.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/imgui/imgui_draw.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/imgui/imgui_tables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/imgui/imgui_widgets.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/imgui/imgui_demo.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/imgui/backends/imgui_impl_glfw.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/imgui/backends/imgui_impl_opengl3.cpp
)

target_include_directories(Sirius PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Tensor"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Autodiff"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Metric"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Geodesic"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Coordinate"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Spectral"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Disk"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Symplectic"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Transport"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Camera"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Jet"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/Polarisation"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/PostProcess"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Render/Pipeline"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Render/Integration"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Render/OptiX"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Render/Acceleration/Backend"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Render/Acceleration/OptiX"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure/Application"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure/Cli"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure/Configuration"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure/Platform"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure/Session"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure/Interface"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure/Window"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure/Plugin"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure/Viewer"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/imgui"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/imgui/backends"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/glm"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/glad/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/stb"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/tinyexr"
)

if(IS_WINDOWS)
    target_compile_options(Sirius PRIVATE /W4 /wd4996 /wd4267 /wd4244)
else()
    target_compile_options(Sirius PRIVATE -Wall -Wextra -O3)
endif()

target_link_libraries(Sirius PRIVATE
    SiriusCore
    glfw
    OpenGL::GL
    Threads::Threads
    ftxui::screen
    ftxui::dom
    ftxui::component
    nlohmann_json::nlohmann_json
)

# Cross-platform compile definitions (replaces global add_definitions)
target_compile_definitions(Sirius PRIVATE
    $<$<PLATFORM_ID:Windows>:NOMINMAX _USE_MATH_DEFINES>
)

if(SIRIUS_HAS_OPTIX)
    target_sources(Sirius PRIVATE
        src/Sirius.Render/Acceleration/OptiX/RDOX001A.cpp
    )
    target_link_libraries(Sirius PRIVATE SiriusOptiX)
    target_include_directories(Sirius PRIVATE
        "${OptiX_INCLUDE_DIR}"
        "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}"
    )
    target_compile_definitions(Sirius PRIVATE SIRIUS_HAS_OPTIX=1)
endif()

set_target_properties(Sirius PROPERTIES
    BUILD_RPATH "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    INSTALL_RPATH "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
)

if(NOT IS_WINDOWS)
    target_link_libraries(Sirius PRIVATE dl)
endif()

# === Resource Copying ===
set(RESOURCE_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Sirius.Render")
file(MAKE_DIRECTORY ${RESOURCE_DIR})

add_custom_command(
    TARGET Sirius POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Render/Texture"
        "${RESOURCE_DIR}/Texture"
    COMMENT "Copying textures to build directory"
)

# === Status Messages ===
message(STATUS "")
message(STATUS "Sirius Build Configuration:")
message(STATUS "=====================================")
if(IS_WINDOWS)
    message(STATUS "Platform: Windows")
    message(STATUS "1. cmake --build bin/Sirius.Build --config Release")
    message(STATUS "2. The executable will be in bin/Sirius.Build/bin/Release/")
else()
    message(STATUS "Platform: Linux")
    message(STATUS "1. cmake -B bin/Sirius.Build && cmake --build bin/Sirius.Build")
    message(STATUS "2. Run: ./bin/Sirius.Build/bin/Sirius")
endif()

message(STATUS "")
message(STATUS "Features:")
message(STATUS "- OpenGL 4.6 Compute Shader GPU acceleration")
if(SIRIUS_HAS_OPTIX)
    message(STATUS "- OptiX ${OptiX_INSTALL_DIR} RTX ray tracing backend [ENABLED]")
else()
    message(STATUS "- OptiX RTX ray tracing backend [DISABLED]")
endif()
message(STATUS "- Geodesic integration with Christoffel symbols")
message(STATUS "- Multiple spacetime metrics (Minkowski, Schwarzschild, Kerr, RN)")
message(STATUS "- Full relativistic redshift calculation")
message(STATUS "- Adaptive step size geodesic integration")
message(STATUS "=========================================")

# === Testing ===
enable_testing()
option(BUILD_TESTS "Build test suite" ON)
if(BUILD_TESTS)
    add_subdirectory(src/Sirius.Test)
endif()

# === Offline Rendering ===
option(BUILD_OFFLINE "Build offline rendering tools" ON)
if(BUILD_OFFLINE)
    add_subdirectory(src/Sirius.Render/Session)
endif()

# =============================================================================
# MANDATORY TEST GATING
# Governance: docs/specification.md
# =============================================================================
option(SIRIUS_MANDATORY_TESTS "Require mandatory tests to pass for build" OFF)
option(SIRIUS_SKIP_MANDATORY "Skip mandatory tests (development override)" OFF)

if(BUILD_TESTS AND SIRIUS_MANDATORY_TESTS AND NOT SIRIUS_SKIP_MANDATORY)
    add_custom_target(RunMandatoryTests ALL
        COMMAND ${CMAKE_CTEST_COMMAND} 
            --test-dir "${CMAKE_BINARY_DIR}"
            --output-on-failure
            --no-tests=error
            -L "Mandatory"
            -C $<CONFIG>
        DEPENDS SiriusTests
        COMMENT "=== MANDATORY TEST GATE: Running mathematical validity tests ==="
    )
    add_dependencies(Sirius RunMandatoryTests)
    
    message(STATUS "")
    message(STATUS "=== MANDATORY TEST GATING: ENABLED ===")
    message(STATUS "  Build will FAIL if tests labeled 'Mandatory' do not pass.")
    message(STATUS "  To skip (development only): -DSIRIUS_SKIP_MANDATORY=ON")
    message(STATUS "")
elseif(SIRIUS_SKIP_MANDATORY)
    message(WARNING "MANDATORY TEST GATING: DISABLED (development override)")
endif()