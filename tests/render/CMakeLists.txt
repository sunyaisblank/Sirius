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
)

target_link_libraries(sirius_render_tests PRIVATE sirius_render GTest::gtest_main)

# Decode-side headers (stbi_load, LoadEXR); the implementations are provided by
# sirius_render, so the tests need declarations only.
target_include_directories(sirius_render_tests PRIVATE
    "${CMAKE_SOURCE_DIR}/lib/stb"
    "${CMAKE_SOURCE_DIR}/lib/tinyexr"
)

sirius_configure_target(sirius_render_tests CONTRACT_MODE 2)

gtest_discover_tests(sirius_render_tests PROPERTIES LABELS "Correctness")
