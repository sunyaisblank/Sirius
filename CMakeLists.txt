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

if(NOT EXISTS "${OptiX_INSTALL_DIR}/include/optix.h")
    if(DEFINED ENV{OptiX_INSTALL_DIR} AND EXISTS "$ENV{OptiX_INSTALL_DIR}/include/optix.h")
        set(OptiX_INSTALL_DIR "$ENV{OptiX_INSTALL_DIR}" CACHE PATH "Path to OptiX SDK installation" FORCE)
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

# --- Platform-Specific Setup ---
if(WIN32)
    set(IS_WINDOWS TRUE)
    add_definitions(-D_WIN32 -DNOMINMAX -D_USE_MATH_DEFINES)
else()
    set(IS_WINDOWS FALSE)
    add_definitions(-D_GNU_SOURCE)
endif()

# === Find System Dependencies ===
find_package(OpenGL REQUIRED)
find_package(Threads REQUIRED)

# === Add Subdirectories for Dependencies ===
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/lib/glfw)

# =============================================================================
# Core Library
# =============================================================================
add_library(SiriusCore SHARED
    # Geodesic integration
    src/Sirius.Core/Geodesic/PHGD001A.cpp
    src/Sirius.Core/Geodesic/PHGD002A.cpp
    src/Sirius.Core/Geodesic/PHMT000A.cpp
    
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
    set(OPTIX_PTX_OUTPUT "${OPTIX_PTX_DIR}/RDOP002A.ptx")
    
    add_custom_command(
        OUTPUT ${OPTIX_PTX_OUTPUT}
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
            -arch compute_75
            -DSIRIUS_HAS_OPTIX=1
            -D_WIN32 -DNOMINMAX -D_USE_MATH_DEFINES
            ${OPTIX_DEVICE_SOURCE}
            -o ${OPTIX_PTX_OUTPUT}
        DEPENDS ${OPTIX_DEVICE_SOURCE}
            "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Render/Acceleration/OptiX/RDOP003A.h"
        COMMENT "Compiling OptiX device kernels to PTX"
        VERBATIM
    )
    
    add_custom_target(SiriusOptiXPTX ALL DEPENDS ${OPTIX_PTX_OUTPUT})
    
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
    # Infrastructure - Application (entry point, orchestration)
    src/Sirius.Infrastructure/Application/CREP001A.cpp
    src/Sirius.Infrastructure/Application/CRAP001A.cpp
    
    # Infrastructure - Session (FSM-based render pipeline)
    src/Sirius.Infrastructure/Session/SNRS001A.cpp
    src/Sirius.Infrastructure/Session/IRRJ001A.cpp
    src/Sirius.Infrastructure/Session/IREX001A.cpp
    
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
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Core/PostProcess"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Render/Pipeline"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Render/OptiX"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure/Application"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure/Session"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure/Interface"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure/Window"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure/Plugin"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/Sirius.Infrastructure/Configuration"
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
)

if(SIRIUS_HAS_OPTIX)
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