# sirius_app: the CLI, configuration loader, platform paths, and the OpenGL
# viewer, plus the `sirius` executable (the name build/render scripts invoke).
#
# All app code except the entry point compiles into sirius_app_lib so the
# CLI/config objects can be linked into the test target (July ledger item 8:
# the CLI parsing tests were never linked to the CLI objects). The executable is
# a thin wrapper over that library.
#
# GLFW is vendored and OpenGL is a system package; the legacy tree also brings
# both in, so each is guarded to coexist when SIRIUS_BUILD_LEGACY is ON. glad.c
# is vendored C: its warnings are relaxed per-file, exactly as the render target
# relaxes miniz.c, so the strict -Werror never promotes them.

find_package(OpenGL REQUIRED)

if(NOT TARGET glfw)
    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
    add_subdirectory("${CMAKE_SOURCE_DIR}/lib/glfw" "${CMAKE_BINARY_DIR}/lib/glfw")
endif()

set(SIRIUS_GLAD "${CMAKE_SOURCE_DIR}/lib/glad/src/glad.c")

add_library(sirius_app_lib STATIC
    cli/command_router.cpp
    cli/render_command.cpp
    cli/info_command.cpp
    cli/config_command.cpp
    cli/cli_output.cpp
    cli/view_command.cpp
    config/config_loader.cpp
    platform_paths.cpp

    # Vendored GL loader (compiled here, warnings relaxed below).
    "${SIRIUS_GLAD}"
)

target_link_libraries(sirius_app_lib PUBLIC
    sirius_render
    glfw
    OpenGL::GL
    ftxui::screen
    ftxui::dom
    ftxui::component
    nlohmann_json::nlohmann_json
    Threads::Threads
)

# glad's headers are consumed by the viewer sources in this library.
target_include_directories(sirius_app_lib PUBLIC "${CMAKE_SOURCE_DIR}/lib/glad/include")

# Honest backend reporting: when the Vulkan backend target was configured (its
# CMakeLists returns early when Vulkan dev files are absent), link it and let
# `info` enumerate devices through the backend seam.
if(TARGET sirius_backend)
    target_link_libraries(sirius_app_lib PUBLIC sirius_backend)
    target_compile_definitions(sirius_app_lib PUBLIC SIRIUS_HAS_VULKAN_BACKEND=1)
endif()

sirius_configure_target(sirius_app_lib)

# Relax warnings on the vendored glad TU only; "-w" is appended after the target
# options so -Werror has nothing to promote.
set_source_files_properties("${SIRIUS_GLAD}" PROPERTIES COMPILE_OPTIONS "-w")

add_executable(sirius main.cpp)
target_link_libraries(sirius PRIVATE sirius_app_lib)
sirius_configure_target(sirius)

# Copy the viewer shaders next to the executable so the runtime resource path
# (PlatformPaths::ResolveResource searches <exe>/shaders) resolves them; this is
# the resource-copy pattern the legacy build used for shader files.
add_custom_command(TARGET sirius POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/viewer/shaders"
        "$<TARGET_FILE_DIR:sirius>/shaders"
    COMMENT "Copying viewer shaders next to the sirius executable"
)
