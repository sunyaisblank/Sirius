# Core suite. The orchestrator owns this file; port workstreams deliver
# .cpp files and report them, and each lands here once its gate passes.

add_executable(sirius_core_tests
    autodiff_dual_test.cpp
    blackbody_planck_test.cpp
    cpu_geodesic_reference_test.cpp
    render_pipeline_integration_test.cpp
    blackbody_test.cpp
    camera_test.cpp
    cie1931_observer_authority_test.cpp
    coordinates_test.cpp
    geodesic_conservation_test.cpp
    geodesic_determinism_test.cpp
    geodesic_equation_test.cpp
    geodesic_live_path_test.cpp
    geodesic_naninf_test.cpp
    geodesic_path_test.cpp
    geodesic_stability_test.cpp
    killing_vector_conservation_test.cpp
    christoffel_test.cpp
    kerr_orbit_authority_test.cpp
    kerr_test.cpp
    metric_derivative_test.cpp
    metric_loader_test.cpp
    metric_registry_test.cpp
    morris_thorne_cartesian_test.cpp
    reissner_nordstrom_test.cpp
    rk45_integrator_test.cpp
    schwarzschild_test.cpp
    walker_penrose_test.cpp
    warp_drive_test.cpp
    stokes_test.cpp
    tensor_inverse_test.cpp
    novikov_thorne_disk_test.cpp
    turbulence_test.cpp
    volumetric_disk_test.cpp
    postprocess_test.cpp
    srgb_transfer_authority_test.cpp
    xyz_srgb_authority_test.cpp
    spectral_radiance_test.cpp
    starfield_test.cpp
    tensor_test.cpp
)

target_link_libraries(sirius_core_tests PRIVATE sirius_core GTest::gtest_main)

sirius_configure_target(sirius_core_tests CONTRACT_MODE 2)

gtest_discover_tests(
    sirius_core_tests
    PROPERTIES
        LABELS "Correctness")
