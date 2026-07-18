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

# Module suites land with their modules; same one-owner rule as src/sirius.
# Each module directory also receives the generated label file as a ctest
# include processed after its discovery include: CTest scopes test properties
# per directory, so the attachment happens per module, guarded by if(TEST).
foreach(module core oracle backend render app)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${module}/CMakeLists.txt")
        add_subdirectory(${module})
        set_property(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${module}" APPEND
            PROPERTY TEST_INCLUDE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/labels/CTestLabels.cmake")
    endif()
endforeach()

# The base suite lives in this directory; label it here too.
set_property(DIRECTORY APPEND
    PROPERTY TEST_INCLUDE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/labels/CTestLabels.cmake")
