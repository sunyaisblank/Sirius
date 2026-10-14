# sirius_render: the CPU render library (output writers, film pipeline, session
# orchestrator) plus the CPU geodesic tracer. The tracer's descriptive home is
# src/sirius/backend/cpu (per docs/ARCHITECTURE.md), but the backend target owns
# a separate CMakeLists this workstream must not touch; the tracer source is
# compiled here for now and re-homes into sirius_backend once that seam is ready.
#
# Vendored deps (stb, tinyexr + miniz) are used exactly as the legacy render
# target used them. Their translation units must not be poisoned by the strict
# -Werror the new tree applies (and miniz.c is C, so the target's -std=c++2c is
# not even valid for it): each such TU gets a per-source "-w" below, relaxing
# only those files and never the whole target.

set(SIRIUS_RENDER_MINIZ "${CMAKE_SOURCE_DIR}/lib/tinyexr/miniz.c")

add_library(sirius_render STATIC
    # Tracer (deferred re-home into sirius_backend; see header comment above).
    "${CMAKE_CURRENT_SOURCE_DIR}/../backend/cpu/geodesic_tracer.cpp"

    # Output writers.
    exr_writer.cpp
    png_writer.cpp

    # stb decoder translation unit (encoder impl lives in png_writer.cpp).
    stb_impl.cpp

    # Session orchestrator.
    session/render_session.cpp

    # Vendored deflate codec for tinyexr.
    "${SIRIUS_RENDER_MINIZ}"
)

target_link_libraries(sirius_render PUBLIC sirius_core Threads::Threads)

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
#   - exr_writer.cpp: hosts TINYEXR_IMPLEMENTATION, which trips
#                     -Werror=calloc-transposed-args inside tinyexr's own source.
# png_writer.cpp hosts STB_IMAGE_WRITE_IMPLEMENTATION but compiles clean under
# the strict flags, so it stays checked (its writer logic is first-party).
set_source_files_properties(
    "${SIRIUS_RENDER_MINIZ}"
    stb_impl.cpp
    exr_writer.cpp
    PROPERTIES COMPILE_OPTIONS "-w"
)
