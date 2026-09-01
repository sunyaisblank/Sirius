# New-tree test estate. One executable per subject area; names and labels
# per docs/STYLE.md section 10. Labels stay hand-set only until the label
# generator learns the new tree (specification programme 2).

include(GoogleTest)

# A required operational profile may not convert missing capabilities or
# degenerate fixtures into green evidence. GoogleTest discovery normally maps
# GTEST_SKIP() to CTest's successful skip state. Force-include a test-only
# policy header which turns GTEST_SKIP into a fatal assertion in that profile.
if(SIRIUS_REQUIRE_VULKAN_RUNTIME)
    if(MSVC)
        add_compile_options(
            "$<$<COMPILE_LANGUAGE:CXX>:/FI${CMAKE_CURRENT_SOURCE_DIR}/support/fail_on_skip.h>")
    else()
        add_compile_options(
            "$<$<COMPILE_LANGUAGE:CXX>:-include${CMAKE_CURRENT_SOURCE_DIR}/support/fail_on_skip.h>")
    endif()
endif()

add_executable(sirius_base_tests
    base/contracts_test.cpp
    base/error_test.cpp
    base/sha256_test.cpp
)

target_link_libraries(sirius_base_tests PRIVATE sirius_base GTest::gtest_main)
target_include_directories(sirius_base_tests PRIVATE "${CMAKE_SOURCE_DIR}/tests")

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

# Source-to-CTest governance rejects unlabelled dynamic registration, and the
# operating-model verifier requires every P/E acceptance criterion as a
# Mandatory-backed obligation.
configure_file(
    operational/evidence_governance.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/evidence_governance.cmake"
    @ONLY)
add_test(
    NAME OperationalEvidence.SourceAndIdealGovernanceRejectUncoveredClaims
    COMMAND "${CMAKE_COMMAND}"
        -DSIRIUS_TEST_CONFIG=$<CONFIG>
        -P "${CMAKE_CURRENT_BINARY_DIR}/evidence_governance.cmake")
set_tests_properties(
    OperationalEvidence.SourceAndIdealGovernanceRejectUncoveredClaims
    PROPERTIES
        LABELS "Mandatory;Operational"
        TIMEOUT 30)

add_test(
    NAME OperationalBuildPolicy.ReleaseCannotWeakenGates
    COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_SOURCE_DIR}/scripts/verify-build-policy.py" --self-test)
set_tests_properties(
    OperationalBuildPolicy.ReleaseCannotWeakenGates
    PROPERTIES
        LABELS "Mandatory;Operational"
        TIMEOUT 10)

add_test(
    NAME OperationalBuildGate.ReleaseInstallRequiresExactPassedEstate
    COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_SOURCE_DIR}/scripts/verify-build-gate.py" --self-test)
set_tests_properties(
    OperationalBuildGate.ReleaseInstallRequiresExactPassedEstate
    PROPERTIES
        LABELS "Mandatory;Operational"
        TIMEOUT 10)

# External platform/device evidence is admitted only through a strict schema.
# This local negative control does not attest hardware; it proves fabricated,
# software-device, wrong-resolution, substituted-model, and tampered-artifact
# records are rejected.
add_test(
    NAME OperationalAttestation.FalseExternalEvidenceIsRejected
    COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_SOURCE_DIR}/scripts/verify-attestation.py" --self-test)
set_tests_properties(
    OperationalAttestation.FalseExternalEvidenceIsRejected
    PROPERTIES
        LABELS "Mandatory;Operational"
        TIMEOUT 30)

# The two native runtime domains have one cross-platform producer with a
# fail-closed host/device policy. Its local self-test is metadata-only: actual
# Windows Vulkan and MoltenVK completion remains external evidence.
add_test(
    NAME OperationalAttestation.PreflightAndNativeRuntimeRejectFalseHostDevice
    COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_SOURCE_DIR}/scripts/validate-native-runtime.py" --self-test)
set_tests_properties(
    OperationalAttestation.PreflightAndNativeRuntimeRejectFalseHostDevice
    PROPERTIES
        LABELS "Mandatory;Operational"
        TIMEOUT 30)

# A valid individual record is still insufficient for release. Admission must
# cover every external domain exactly once, bind the current clean revision,
# and produce a deterministic receipt for compile/runtime comparison.
add_test(
    NAME OperationalAlignment.IncompleteOrAmbiguousAttestationSetIsRejected
    COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_SOURCE_DIR}/scripts/verify-alignment.py" --self-test)
set_tests_properties(
    OperationalAlignment.IncompleteOrAmbiguousAttestationSetIsRejected
    PROPERTIES
        LABELS "Mandatory;Operational"
        TIMEOUT 30)

configure_file(
    operational/alignment_receipt_integrity.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/alignment_receipt_integrity.cmake"
    @ONLY)
add_test(
    NAME OperationalAlignment.InstalledReceiptTamperingBlocksReadiness
    COMMAND "${CMAKE_COMMAND}"
        -DSIRIUS_BINARY=$<TARGET_FILE:sirius>
        -DSIRIUS_RESOURCE_SOURCE=$<TARGET_FILE_DIR:sirius>/resources
        -P "${CMAKE_CURRENT_BINARY_DIR}/alignment_receipt_integrity.cmake")
set_tests_properties(
    OperationalAlignment.InstalledReceiptTamperingBlocksReadiness
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
