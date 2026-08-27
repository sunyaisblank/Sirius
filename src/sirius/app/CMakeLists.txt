# sirius_app: the CLI, configuration loader, platform paths, and the OpenGL
# viewer, plus the `sirius` executable (the name build/render scripts invoke).
#
# All app code except the entry point compiles into sirius_app_lib so the
# CLI/config objects can be linked into the test target (July ledger item 8:
# the CLI parsing tests were never linked to the CLI objects). The executable is
# a thin wrapper over that library.
#
option(SIRIUS_BUILD_VIEWER "Build the GLFW/OpenGL interactive window when available" ON)

# Bind the installed capability ledger to the exact reviewed source document.
# CMAKE_CONFIGURE_DEPENDS ensures editing the JSON reruns configuration before
# the executable can be relinked with a stale embedded authority.
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${SIRIUS_OPERATING_MODEL}")
file(READ "${SIRIUS_OPERATING_MODEL}" SIRIUS_OPERATING_MODEL_JSON)
file(SHA256 "${SIRIUS_OPERATING_MODEL}" SIRIUS_OPERATING_MODEL_SHA256)
string(FIND "${SIRIUS_OPERATING_MODEL_JSON}" ")SIRIUS_MODEL\""
    SIRIUS_OPERATING_MODEL_DELIMITER)
if(NOT SIRIUS_OPERATING_MODEL_DELIMITER EQUAL -1)
    message(FATAL_ERROR
        "operating_model.json contains the reserved embedded-string delimiter")
endif()

# MSVC rejects a single narrow string literal above 16,380 characters (C2026).
# Keep every generated literal below that implementation limit and reconstruct
# the authoritative JSON at the runtime validation boundary.
set(SIRIUS_OPERATING_MODEL_CPP_CHUNKS "")
set(SIRIUS_OPERATING_MODEL_CHUNK_COUNT 0)
set(SIRIUS_OPERATING_MODEL_OFFSET 0)
set(SIRIUS_OPERATING_MODEL_CHUNK_SIZE 8000)
string(LENGTH "${SIRIUS_OPERATING_MODEL_JSON}"
    SIRIUS_OPERATING_MODEL_SIZE)
while(SIRIUS_OPERATING_MODEL_OFFSET LESS SIRIUS_OPERATING_MODEL_SIZE)
    math(EXPR SIRIUS_OPERATING_MODEL_REMAINING
        "${SIRIUS_OPERATING_MODEL_SIZE} - ${SIRIUS_OPERATING_MODEL_OFFSET}")
    if(SIRIUS_OPERATING_MODEL_REMAINING GREATER SIRIUS_OPERATING_MODEL_CHUNK_SIZE)
        set(SIRIUS_OPERATING_MODEL_CURRENT_SIZE
            "${SIRIUS_OPERATING_MODEL_CHUNK_SIZE}")
    else()
        set(SIRIUS_OPERATING_MODEL_CURRENT_SIZE
            "${SIRIUS_OPERATING_MODEL_REMAINING}")
    endif()
    string(SUBSTRING "${SIRIUS_OPERATING_MODEL_JSON}"
        "${SIRIUS_OPERATING_MODEL_OFFSET}"
        "${SIRIUS_OPERATING_MODEL_CURRENT_SIZE}"
        SIRIUS_OPERATING_MODEL_CHUNK)
    string(APPEND SIRIUS_OPERATING_MODEL_CPP_CHUNKS
        "    R\"SIRIUS_MODEL(${SIRIUS_OPERATING_MODEL_CHUNK})SIRIUS_MODEL\",\n")
    math(EXPR SIRIUS_OPERATING_MODEL_OFFSET
        "${SIRIUS_OPERATING_MODEL_OFFSET} + ${SIRIUS_OPERATING_MODEL_CURRENT_SIZE}")
    math(EXPR SIRIUS_OPERATING_MODEL_CHUNK_COUNT
        "${SIRIUS_OPERATING_MODEL_CHUNK_COUNT} + 1")
endwhile()

# Embed the alignment receipt beside the operating model. Release mode reaches
# this point only after the root configure gate has admitted every required
# external domain for the exact clean revision.
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${SIRIUS_ALIGNMENT_RECEIPT}")
file(READ "${SIRIUS_ALIGNMENT_RECEIPT}" SIRIUS_ALIGNMENT_RECEIPT_JSON)
string(FIND "${SIRIUS_ALIGNMENT_RECEIPT_JSON}" ")SIRIUS_ALIGN\""
    SIRIUS_ALIGNMENT_RECEIPT_DELIMITER)
if(NOT SIRIUS_ALIGNMENT_RECEIPT_DELIMITER EQUAL -1)
    message(FATAL_ERROR
        "alignment_receipt.json contains the reserved embedded-string delimiter")
endif()
string(JSON SIRIUS_ALIGNMENT_RECEIPT_SATISFIED
    GET "${SIRIUS_ALIGNMENT_RECEIPT_JSON}" ultimate_ideal satisfied)
if(SIRIUS_ALIGNMENT_MODE STREQUAL "release")
    set(SIRIUS_RELEASE_ALIGNMENT_ENFORCED true)
    set(SIRIUS_QUALIFICATION_ALIGNMENT_ENFORCED false)
    if(NOT SIRIUS_ALIGNMENT_RECEIPT_SATISFIED)
        message(FATAL_ERROR
            "Release alignment reached compilation without a satisfied receipt")
    endif()
elseif(SIRIUS_ALIGNMENT_MODE STREQUAL "qualification")
    set(SIRIUS_RELEASE_ALIGNMENT_ENFORCED false)
    set(SIRIUS_QUALIFICATION_ALIGNMENT_ENFORCED true)
else()
    set(SIRIUS_RELEASE_ALIGNMENT_ENFORCED false)
    set(SIRIUS_QUALIFICATION_ALIGNMENT_ENFORCED false)
endif()
if(SIRIUS_ALIGNMENT_RECEIPT_SATISFIED)
    set(SIRIUS_EMBEDDED_ALIGNMENT_SATISFIED true)
else()
    set(SIRIUS_EMBEDDED_ALIGNMENT_SATISFIED false)
endif()

set(SIRIUS_ALIGNMENT_RECEIPT_CPP_CHUNKS "")
set(SIRIUS_ALIGNMENT_RECEIPT_CHUNK_COUNT 0)
set(SIRIUS_ALIGNMENT_RECEIPT_OFFSET 0)
set(SIRIUS_ALIGNMENT_RECEIPT_CHUNK_SIZE 8000)
string(LENGTH "${SIRIUS_ALIGNMENT_RECEIPT_JSON}"
    SIRIUS_ALIGNMENT_RECEIPT_SIZE)
while(SIRIUS_ALIGNMENT_RECEIPT_OFFSET LESS SIRIUS_ALIGNMENT_RECEIPT_SIZE)
    math(EXPR SIRIUS_ALIGNMENT_RECEIPT_REMAINING
        "${SIRIUS_ALIGNMENT_RECEIPT_SIZE} - ${SIRIUS_ALIGNMENT_RECEIPT_OFFSET}")
    if(SIRIUS_ALIGNMENT_RECEIPT_REMAINING GREATER SIRIUS_ALIGNMENT_RECEIPT_CHUNK_SIZE)
        set(SIRIUS_ALIGNMENT_RECEIPT_CURRENT_SIZE
            "${SIRIUS_ALIGNMENT_RECEIPT_CHUNK_SIZE}")
    else()
        set(SIRIUS_ALIGNMENT_RECEIPT_CURRENT_SIZE
            "${SIRIUS_ALIGNMENT_RECEIPT_REMAINING}")
    endif()
    string(SUBSTRING "${SIRIUS_ALIGNMENT_RECEIPT_JSON}"
        "${SIRIUS_ALIGNMENT_RECEIPT_OFFSET}"
        "${SIRIUS_ALIGNMENT_RECEIPT_CURRENT_SIZE}"
        SIRIUS_ALIGNMENT_RECEIPT_CHUNK)
    string(APPEND SIRIUS_ALIGNMENT_RECEIPT_CPP_CHUNKS
        "    R\"SIRIUS_ALIGN(${SIRIUS_ALIGNMENT_RECEIPT_CHUNK})SIRIUS_ALIGN\",\n")
    math(EXPR SIRIUS_ALIGNMENT_RECEIPT_OFFSET
        "${SIRIUS_ALIGNMENT_RECEIPT_OFFSET} + ${SIRIUS_ALIGNMENT_RECEIPT_CURRENT_SIZE}")
    math(EXPR SIRIUS_ALIGNMENT_RECEIPT_CHUNK_COUNT
        "${SIRIUS_ALIGNMENT_RECEIPT_CHUNK_COUNT} + 1")
endwhile()

file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/generated/sirius/base")
configure_file(
    "${CMAKE_SOURCE_DIR}/cmake/operating_model_embedded.h.in"
    "${CMAKE_BINARY_DIR}/generated/sirius/base/operating_model_embedded.h"
    @ONLY)

set(SIRIUS_VIEWER_AVAILABLE FALSE)
if(SIRIUS_BUILD_VIEWER)
    find_package(OpenGL QUIET)
    if(OpenGL_FOUND)
        set(SIRIUS_VIEWER_AVAILABLE TRUE)
        if(NOT TARGET glfw)
            set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
            set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
            set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
            set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
            add_subdirectory("${CMAKE_SOURCE_DIR}/lib/glfw" "${CMAKE_BINARY_DIR}/lib/glfw")
        endif()
    else()
        message(STATUS
            "Sirius viewer: OpenGL development files absent; building the complete "
            "CPU/Vulkan CLI with a fail-closed `view` command")
    endif()
endif()
set(SIRIUS_INSTALL_HAS_VIEWER "${SIRIUS_VIEWER_AVAILABLE}" CACHE INTERNAL
    "Whether the installed Sirius volume contains the interactive viewer" FORCE)

set(SIRIUS_APP_SOURCES
    alignment_authority.cpp
    cli/command_router.cpp
    cli/render_command.cpp
    cli/info_command.cpp
    cli/config_command.cpp
    cli/cli_output.cpp
    config/config_loader.cpp
    config/session_config_adapter.cpp
    platform_paths.cpp
)

if(SIRIUS_VIEWER_AVAILABLE)
    set(SIRIUS_GLAD "${CMAKE_SOURCE_DIR}/lib/glad/src/glad.c")
    list(APPEND SIRIUS_APP_SOURCES
        cli/view_command.cpp
        viewer/interactive_viewer.cpp
        "${SIRIUS_GLAD}")
else()
    list(APPEND SIRIUS_APP_SOURCES cli/view_command_stub.cpp)
endif()

add_library(sirius_app_lib STATIC ${SIRIUS_APP_SOURCES})
add_dependencies(sirius_app_lib SiriusAlignmentGate)

target_include_directories(sirius_app_lib PRIVATE "${CMAKE_BINARY_DIR}/generated")
target_link_libraries(sirius_app_lib PUBLIC
    sirius_render
    ftxui::screen
    ftxui::dom
    ftxui::component
    nlohmann_json::nlohmann_json
    Threads::Threads
)

target_compile_definitions(sirius_app_lib PUBLIC SIRIUS_VERSION="${PROJECT_VERSION}")
if(SIRIUS_VIEWER_AVAILABLE)
    target_link_libraries(sirius_app_lib PUBLIC glfw OpenGL::GL)
    target_include_directories(sirius_app_lib PUBLIC "${CMAKE_SOURCE_DIR}/lib/glad/include")
    target_compile_definitions(sirius_app_lib PUBLIC SIRIUS_HAS_INTERACTIVE_VIEWER=1)
endif()

# Honest backend reporting: when the Vulkan backend target was configured (its
# CMakeLists returns early when Vulkan dev files are absent), link it and let
# `info` enumerate devices through the backend seam.
if(TARGET sirius_backend AND TARGET sirius_kernels)
    target_link_libraries(sirius_app_lib PUBLIC sirius_backend)
    target_compile_definitions(sirius_app_lib PUBLIC SIRIUS_HAS_VULKAN_BACKEND=1)
endif()

sirius_configure_target(sirius_app_lib)

if(SIRIUS_VIEWER_AVAILABLE)
    # Relax warnings on the vendored glad TU only; "-w" is appended after the
    # target options so -Werror has nothing to promote.
    set_source_files_properties("${SIRIUS_GLAD}" PROPERTIES COMPILE_OPTIONS "-w")
endif()

add_executable(sirius main.cpp)
target_link_libraries(sirius PRIVATE sirius_app_lib)
sirius_configure_target(sirius)
set_property(TARGET sirius APPEND PROPERTY LINK_DEPENDS
    "${SIRIUS_OPERATING_MODEL}")

# Stage the exact runtime tree used by installed volumes beside the build-tree
# executable. Resource discovery is therefore identical in both environments.
add_custom_command(TARGET sirius POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
        "$<TARGET_FILE_DIR:sirius>/resources/assets"
    COMMAND ${CMAKE_COMMAND} -E make_directory
        "$<TARGET_FILE_DIR:sirius>/resources/model"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_SOURCE_DIR}/assets/Starfield.png"
        "$<TARGET_FILE_DIR:sirius>/resources/assets/Starfield.png"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${SIRIUS_OPERATING_MODEL}"
        "$<TARGET_FILE_DIR:sirius>/resources/model/operating_model.json"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${SIRIUS_ALIGNMENT_RECEIPT}"
        "$<TARGET_FILE_DIR:sirius>/resources/model/alignment_receipt.json"
    COMMENT "Staging Sirius runtime resources beside the executable"
)

if(SIRIUS_VIEWER_AVAILABLE)
    add_custom_command(TARGET sirius POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
            "$<TARGET_FILE_DIR:sirius>/resources/shaders"
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/viewer/shaders"
            "$<TARGET_FILE_DIR:sirius>/resources/shaders"
        COMMENT "Staging interactive-viewer shaders beside the executable")
endif()

if(TARGET sirius_kernels)
    add_custom_command(TARGET sirius POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
            "$<TARGET_FILE_DIR:sirius>/resources/kernels"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${SIRIUS_KERNEL_BINARY_DIR}/trace.spv"
            "${SIRIUS_KERNEL_BINARY_DIR}/trace_fp32comp.spv"
            "${SIRIUS_KERNEL_BINARY_DIR}/trace_fp64.spv"
            "$<TARGET_FILE_DIR:sirius>/resources/kernels"
        COMMENT "Staging Vulkan runtime kernels beside the executable"
    )
endif()
