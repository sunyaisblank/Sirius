# Sirius.Offline - GPU-Accelerated Offline Rendering Module
# Creates the SiriusRender executable for high-quality offline rendering
# Uses OptiX backend for GPU ray tracing (same as real-time renderer)

cmake_minimum_required(VERSION 3.18)

# Get parent project settings
if(NOT TARGET SiriusCore)
    message(FATAL_ERROR "Sirius.Offline must be built as part of the main Sirius project")
endif()

# =============================================================================
# SiriusRender Executable (GPU-Accelerated)
# =============================================================================

add_executable(SiriusRender
    SREP001A.cpp
    "${CMAKE_SOURCE_DIR}/src/Sirius.Render/Acceleration/Backend/ACBF001A.cpp"
    # EXR output support
    "${CMAKE_SOURCE_DIR}/src/Sirius.Render/Output/OUEW001A.cpp"
    "${CMAKE_SOURCE_DIR}/lib/tinyexr/miniz.c"
)

target_include_directories(SiriusRender PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/lib/stb
    ${CMAKE_SOURCE_DIR}/lib/tinyexr
    # Include directories for new module structure
    "${CMAKE_SOURCE_DIR}/src/Sirius.Render/Acceleration/Backend"
    "${CMAKE_SOURCE_DIR}/src/Sirius.Render/Acceleration/OptiX"
    "${CMAKE_SOURCE_DIR}/src/Sirius.Render/Buffer"
    "${CMAKE_SOURCE_DIR}/src/Sirius.Render/Output"
    "${CMAKE_SOURCE_DIR}/src/Sirius.Render/Session"
    "${CMAKE_SOURCE_DIR}/src/Sirius.Core/PostProcess"
)

# Link core library
target_link_libraries(SiriusRender PRIVATE
    SiriusCore
)

# Link OptiX backend if available (required for GPU rendering)
if(SIRIUS_HAS_OPTIX)
    target_sources(SiriusRender PRIVATE 
        "${CMAKE_SOURCE_DIR}/src/Sirius.Render/Acceleration/OptiX/RDOX001A.cpp"
    )
    target_link_libraries(SiriusRender PRIVATE SiriusOptiX)
    target_include_directories(SiriusRender PRIVATE
        "${OptiX_INCLUDE_DIR}"
        "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}"
    )
    target_compile_definitions(SiriusRender PRIVATE SIRIUS_HAS_OPTIX=1)
    message(STATUS "SiriusRender: OptiX GPU backend enabled")
else()
    message(WARNING "SiriusRender: OptiX not available - GPU rendering disabled")
endif()

# C++20 required
target_compile_features(SiriusRender PRIVATE cxx_std_20)

# MSVC settings
if(MSVC)
    target_compile_options(SiriusRender PRIVATE /W4 /permissive-)
endif()

# Output directory
set_target_properties(SiriusRender PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/bin/Debug
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin/Release
)

# Installation
install(TARGETS SiriusRender
    RUNTIME DESTINATION bin
)
