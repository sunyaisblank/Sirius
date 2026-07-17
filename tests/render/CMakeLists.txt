# Render suite: CPU tracer, output writers, film pipeline, spectral emission,
# and the session-level end-to-end probe. Links sirius_render (which brings
# sirius_core and the vendored stb/tinyexr symbols transitively). The vendored
# headers are included by the round-trip decode side of the writer tests and the
# probe, so this target carries the lib/stb and lib/tinyexr include paths.

add_executable(sirius_render_tests
    geodesic_tracer_test.cpp
    exr_writer_test.cpp
    png_writer_test.cpp
    film_pipeline_test.cpp
    spectral_emission_test.cpp
    render_session_probe_test.cpp
    memory_governor_test.cpp
    vulkan_render_test.cpp
)

target_link_libraries(sirius_render_tests PRIVATE sirius_render GTest::gtest_main)

# Decode-side headers (stbi_load, LoadEXR); the implementations are provided by
# sirius_render, so the tests need declarations only.
target_include_directories(sirius_render_tests PRIVATE
    "${CMAKE_SOURCE_DIR}/lib/stb"
    "${CMAKE_SOURCE_DIR}/lib/tinyexr"
)

sirius_configure_target(sirius_render_tests CONTRACT_MODE 2)

# The Vulkan render/parity gate loads the compiled trace kernel; point it at the
# kernel binary dir when slangc produced one (sirius_render carries the same
# SIRIUS_HAS_VULKAN_BACKEND define transitively). The tests skip cleanly when the
# device or the kernel is absent.
if(TARGET sirius_kernels)
    target_compile_definitions(sirius_render_tests
        PRIVATE SIRIUS_KERNEL_DIR="${SIRIUS_KERNEL_BINARY_DIR}")
    add_dependencies(sirius_render_tests sirius_kernels)
endif()

gtest_discover_tests(sirius_render_tests PROPERTIES LABELS "Correctness")
