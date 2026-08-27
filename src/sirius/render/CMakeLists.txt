# sirius_render: render orchestration, output writers, film processing, and
# session lifecycle. Concrete execution belongs to backend targets.
#
# Vendored deps (stb, tinyexr + miniz) remain isolated implementation details.
# Their translation units must not be poisoned by the strict -Werror policy
# (and miniz.c is C, so the target's -std=c++2c flag is deliberately C++-only):
# each such TU gets a per-source "-w" below, relaxing only those files and never
# the whole target.

set(SIRIUS_RENDER_MINIZ "${CMAKE_SOURCE_DIR}/lib/tinyexr/miniz.c")

add_library(sirius_render STATIC
    # Output writers.
    exr_writer.cpp
    png_writer.cpp

    # Vendored stb implementation translation units.
    stb_impl.cpp
    stb_write_impl.cpp

    # Governors (pure arithmetic; no device dependency). Memory bounds tile
    # residency; dispatch bounds per-submission duration under OS watchdogs.
    memory_governor.cpp
    dispatch_governor.cpp

    # Session orchestrator.
    session/render_session.cpp

    # Vendored deflate codec for tinyexr.
    "${SIRIUS_RENDER_MINIZ}"
)

target_link_libraries(sirius_render PUBLIC sirius_core sirius_backend_cpu Threads::Threads)

# Vulkan render path: the trace-kernel dispatch through the ComputeDevice seam.
# The backend target returns early when Vulkan dev files are absent, so this is
# guarded exactly as sirius_app_lib guards its backend link; machines without
# Vulkan build the CPU render library only.
if(TARGET sirius_backend AND TARGET sirius_kernels)
    target_sources(sirius_render PRIVATE vulkan_renderer.cpp)
    target_link_libraries(sirius_render PUBLIC sirius_backend)
    target_compile_definitions(sirius_render PUBLIC SIRIUS_HAS_VULKAN_BACKEND=1)
    add_dependencies(sirius_render sirius_kernels)
endif()

# Vendored headers are included by the writer/decoder TUs only (PRIVATE).
target_include_directories(sirius_render PRIVATE
    "${CMAKE_SOURCE_DIR}/lib/stb"
    "${CMAKE_SOURCE_DIR}/lib/tinyexr"
)

sirius_configure_target(sirius_render)

# Relax warnings on the vendored translation units only. "-w" is appended after
# the target options, so -Werror has nothing to promote.
#   - miniz.c       : vendored C; the target's -std=c++2c is not valid for C.
#   - stb_impl.cpp  : pure vendored STB_IMAGE_IMPLEMENTATION (no first-party code).
#   - stb_write_impl.cpp: pure vendored STB_IMAGE_WRITE_IMPLEMENTATION.
#   - exr_writer.cpp: hosts TINYEXR_IMPLEMENTATION, which trips
#                     -Werror=calloc-transposed-args inside tinyexr's own source.
set_source_files_properties(
    "${SIRIUS_RENDER_MINIZ}"
    stb_impl.cpp
    stb_write_impl.cpp
    exr_writer.cpp
    PROPERTIES COMPILE_OPTIONS "-w"
)
