# New-tree test estate. One executable per subject area; names and labels
# per docs/STYLE.md section 10. Labels stay hand-set only until the label
# generator learns the new tree (specification programme 2).

include(GoogleTest)

# A required operational profile may not convert missing capabilities or
# degenerate fixtures into green evidence. GoogleTest discovery normally maps
# GTEST_SKIP() to CTest's successful skip state. Force-include a test-only
# policy header which turns GTEST_SKIP into a fatal assertion in that profile.
if(SIRIUS_REQUIRE_VULKAN_RUNTIME)
    add_compile_options(
        $<$<COMPILE_LANGUAGE:CXX>:-include${CMAKE_CURRENT_SOURCE_DIR}/support/fail_on_skip.h>)
endif()

add_executable(sirius_base_tests
    base/contracts_test.cpp
    base/error_test.cpp
)

target_link_libraries(sirius_base_tests PRIVATE sirius_base GTest::gtest_main)

# Death tests require enforce semantics regardless of build type.
sirius_configure_target(sirius_base_tests CONTRACT_MODE 2)

gtest_discover_tests(
    sirius_base_tests
    PROPERTIES
        LABELS "Correctness")

# Module suites land with their modules; same one-owner rule as src/sirius.
# Each module directory also receives the generated label file as a ctest
# include processed after its discovery include: CTest scopes test properties
# per directory, so the attachment happens per module, guarded by if(TEST).
foreach(module core oracle backend render app)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${module}/CMakeLists.txt")
        add_subdirectory(${module})
        set_property(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${module}" APPEND
            PROPERTY TEST_INCLUDE_FILES
                "${CMAKE_CURRENT_SOURCE_DIR}/labels/CTestLabels.cmake")
    endif()
endforeach()

# The base suite lives in this directory; label it here too.
set_property(DIRECTORY APPEND
    PROPERTY TEST_INCLUDE_FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/labels/CTestLabels.cmake")

# Negative control for source-to-CTest governance: a synthetic TEST_P must be
# rejected because its dynamically discovered names cannot receive the static
# Mandatory policy without explicit generator support.
add_test(
    NAME OperationalEvidence.SourceGovernanceRejectsUngovernedParameterizedTests
    COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_SOURCE_DIR}/scripts/generate-ctest-labels.py" --self-test)
set_tests_properties(
    OperationalEvidence.SourceGovernanceRejectsUngovernedParameterizedTests
    PROPERTIES
        LABELS "Mandatory;Operational"
        TIMEOUT 30)

# External platform/device evidence is admitted only through a strict schema.
# This local negative control does not attest hardware; it proves fabricated,
# software-device, wrong-resolution, and tampered-artifact records are rejected.
add_test(
    NAME OperationalAttestation.FalseExternalEvidenceIsRejected
    COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_SOURCE_DIR}/scripts/verify-attestation.py" --self-test)
set_tests_properties(
    OperationalAttestation.FalseExternalEvidenceIsRejected
    PROPERTIES
        LABELS "Mandatory;Operational"
        TIMEOUT 30)

# Installation is part of the operating model, not a packaging afterthought.
# This test relocates a staged prefix, initialises from a hostile working
# directory, renders through the CPU path, and proves a missing mandatory
# resource declines loudly.
configure_file(
    operational/volume_initialization.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/volume_initialization.cmake"
    @ONLY)
add_test(
    NAME OperationalVolume.InstallRelocateInitialiseAndFailClosed
    COMMAND "${CMAKE_COMMAND}"
        -DSIRIUS_TEST_CONFIG=$<CONFIG>
        -P "${CMAKE_CURRENT_BINARY_DIR}/volume_initialization.cmake")
set_tests_properties(
    OperationalVolume.InstallRelocateInitialiseAndFailClosed
    PROPERTIES
        LABELS "Mandatory;Operational"
        RUN_SERIAL TRUE
        TIMEOUT 600)

# Operator scripts are part of the callable system surface. Their render
# commands must not be masked by a successful logging pipeline, and an exit-zero
# renderer which emits no artefact must not count as success.
if(UNIX)
    find_program(SIRIUS_BASH_EXECUTABLE bash REQUIRED)
    configure_file(
        operational/operator_script_status.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/operator_script_status.cmake"
        @ONLY)
    add_test(
        NAME OperationalScripts.RenderFailuresAndMissingOutputsPropagate
        COMMAND "${CMAKE_COMMAND}"
            -P "${CMAKE_CURRENT_BINARY_DIR}/operator_script_status.cmake")
    set_tests_properties(
        OperationalScripts.RenderFailuresAndMissingOutputsPropagate
        PROPERTIES
            LABELS "Mandatory;Operational"
            TIMEOUT 30)
endif()

# The required Vulkan runtime profile is intentionally opt-in for portable
# source builds and mandatory for the Linux operational CI profile. Unlike
# device-dependent unit tests, it never self-skips: absence of a device,
# kernels, or a successful governed dispatch is a hard failure.
if(SIRIUS_REQUIRE_VULKAN_RUNTIME)
    add_test(
        NAME OperationalEvidence.RequiredProfileRejectsGoogleTestSkip
        COMMAND "${CMAKE_COMMAND}"
            -DSIRIUS_BACKEND_TESTS=$<TARGET_FILE:sirius_backend_tests>
            -P "${CMAKE_CURRENT_SOURCE_DIR}/operational/required_skip_rejection.cmake")
    set_tests_properties(
        OperationalEvidence.RequiredProfileRejectsGoogleTestSkip
        PROPERTIES
            LABELS "Mandatory;Operational"
            TIMEOUT 30)

    add_test(
        NAME OperationalVulkan.RequiredRuntimeDispatch
        COMMAND "${CMAKE_COMMAND}"
            -DSIRIUS_BINARY=$<TARGET_FILE:sirius>
            -DSIRIUS_OUTPUT=${CMAKE_CURRENT_BINARY_DIR}/required-vulkan-runtime.ppm
            -P "${CMAKE_CURRENT_SOURCE_DIR}/operational/vulkan_runtime_profile.cmake")
    set_tests_properties(
        OperationalVulkan.RequiredRuntimeDispatch
        PROPERTIES
            LABELS "Mandatory;Operational"
            RUN_SERIAL TRUE
            TIMEOUT 120)

    add_test(
        NAME OperationalP1.NearExtremalBurnIn
        COMMAND "${CMAKE_COMMAND}"
            -DSIRIUS_RENDER_TESTS=$<TARGET_FILE:sirius_render_tests>
            -P "${CMAKE_CURRENT_SOURCE_DIR}/operational/p1_burn_in.cmake")
    set_tests_properties(
        OperationalP1.NearExtremalBurnIn
        PROPERTIES
            LABELS "Mandatory;Operational;Stability"
            RUN_SERIAL TRUE
            TIMEOUT 180)
endif()
