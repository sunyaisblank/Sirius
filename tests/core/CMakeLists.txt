# Core suite. The orchestrator owns this file; port workstreams deliver
# .cpp files and report them, and each lands here once its gate passes.

add_executable(sirius_core_tests
    autodiff_dual_test.cpp
    camera_test.cpp
    corona_test.cpp
    jet_test.cpp
    novikov_thorne_disk_test.cpp
    postprocess_test.cpp
    spectral_radiance_test.cpp
    starfield_test.cpp
    tensor_test.cpp
)

target_link_libraries(sirius_core_tests PRIVATE sirius_core GTest::gtest_main)

sirius_configure_target(sirius_core_tests CONTRACT_MODE 2)

gtest_discover_tests(sirius_core_tests PROPERTIES LABELS "Correctness")
