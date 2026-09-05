# Render suite: every render-capable CLI/viewer case, CPU tracer, output writers,
# film pipeline, spectral emission, and the session-level end-to-end probe.
# Links sirius_app_lib for command/viewer execution and sirius_render for the
# lower render boundaries. The vendored headers are included by the round-trip
# decode side of the writer tests and the probe, so this target carries the
# lib/stb and lib/tinyexr include paths.

add_executable(sirius_render_tests
    app_render_execution_test.cpp
    geodesic_tracer_test.cpp
    exr_writer_test.cpp
    png_writer_test.cpp
    film_pipeline_test.cpp
    spectral_emission_test.cpp
    render_session_probe_test.cpp
    memory_governor_test.cpp
    dispatch_governor_test.cpp
    vulkan_render_test.cpp
    ray_bundle_test.cpp
    starfield_point_test.cpp
    doppler_toggle_test.cpp
    camera_aberration_test.cpp
    pixel_sampling_test.cpp
    shadow_boundary_test.cpp
)

target_link_libraries(sirius_render_tests PRIVATE sirius_render sirius_app_lib GTest::gtest_main)

if(SIRIUS_REQUIRE_VULKAN_RUNTIME)
    target_compile_definitions(sirius_render_tests
        PRIVATE SIRIUS_TEST_REQUIRE_VULKAN_RUNTIME=1)
endif()

# Decode-side headers (stbi_load, LoadEXR); the implementations are provided by
# sirius_render, so the tests need declarations only.
target_include_directories(sirius_render_tests PRIVATE
    "${CMAKE_SOURCE_DIR}/tests"
    "${CMAKE_SOURCE_DIR}/lib/stb"
    "${CMAKE_SOURCE_DIR}/lib/tinyexr"
)

sirius_configure_target(sirius_render_tests CONTRACT_MODE 2)

if(TARGET sirius)
    # Session-level tests consume the exact executable-adjacent runtime volume
    # in every alignment mode. This keeps strict resource locking active rather
    # than substituting a development-only environment-selected tree.
    add_dependencies(sirius_render_tests sirius)
    set_target_properties(sirius_render_tests PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "$<TARGET_FILE_DIR:sirius>")
endif()

# The Vulkan render/parity gate loads the compiled trace kernel; point it at the
# kernel binary dir when slangc produced one (sirius_render carries the same
# SIRIUS_HAS_VULKAN_BACKEND define transitively). The tests skip cleanly when the
# device or the kernel is absent.
if(TARGET sirius_kernels)
    target_compile_definitions(sirius_render_tests
        PRIVATE SIRIUS_KERNEL_DIR="${SIRIUS_KERNEL_BINARY_DIR}")
    add_dependencies(sirius_render_tests sirius_kernels)
endif()

gtest_discover_tests(
    sirius_render_tests
    PROPERTIES
        LABELS "Correctness")
