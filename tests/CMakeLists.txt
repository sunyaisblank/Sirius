# New-tree test estate. One executable per subject area; names and labels
# per docs/STYLE.md section 10. Labels stay hand-set only until the label
# generator learns the new tree (specification programme 2).

include(GoogleTest)

add_executable(sirius_base_tests
    base/contracts_test.cpp
    base/error_test.cpp
)

target_link_libraries(sirius_base_tests PRIVATE sirius_base GTest::gtest_main)

# Death tests require enforce semantics regardless of build type.
sirius_configure_target(sirius_base_tests CONTRACT_MODE 2)

gtest_discover_tests(sirius_base_tests PROPERTIES LABELS "Correctness")
