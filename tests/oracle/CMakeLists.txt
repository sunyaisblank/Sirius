# Oracle test estate: the pure double-precision Boyer-Lindquist suites whose
# only project dependency is sirius_oracle (plus GoogleTest). Mixed suites
# that also reach into sirius_core (TSDG007A/008A, TSVAL001A) land with the
# core workstream, not here. Contract mode 2 so the constructor preconditions
# enforce under test.

add_executable(sirius_oracle_tests
    analytic_benchmarks_test.cpp
    analytic_solutions_test.cpp
    metric_invariants_test.cpp
    transport_types_test.cpp
    kerr_metric_test.cpp
    symplectic_integrator_test.cpp
    beam_propagation_test.cpp
    einstein_ring_test.cpp
    geodesic_deviation_test.cpp
)

target_link_libraries(sirius_oracle_tests PRIVATE sirius_oracle sirius_core GTest::gtest_main)

sirius_configure_target(sirius_oracle_tests CONTRACT_MODE 2)

gtest_discover_tests(sirius_oracle_tests PROPERTIES LABELS "Correctness")
