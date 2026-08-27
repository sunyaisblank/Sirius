# =============================================================================
# Sirius: general-relativistic renderer.
# Governance: docs/SPECIFICATION.md, docs/ARCHITECTURE.md, docs/STYLE.md.
# =============================================================================
cmake_minimum_required(VERSION 3.28)
project(Sirius VERSION 2.0.0 LANGUAGES CXX C)

option(BUILD_TESTS "Build test suites" ON)
option(SIRIUS_MANDATORY_TESTS "Require mandatory tests to pass for every normal build" ON)
option(SIRIUS_REQUIRE_VULKAN_RUNTIME
    "Require compiled Vulkan kernels and a dispatchable Vulkan device in the operational gate"
    OFF)
set(SIRIUS_ALIGNMENT_MODE "development" CACHE STRING
    "Alignment enforcement: development is local; qualification is strict evidence production; release requires all attestations")
set_property(CACHE SIRIUS_ALIGNMENT_MODE PROPERTY STRINGS development qualification release)
set(SIRIUS_ATTESTATION_ROOT "" CACHE PATH
    "Root containing revision-bound external attestation bundles")

if(NOT SIRIUS_ALIGNMENT_MODE STREQUAL "development" AND
   NOT SIRIUS_ALIGNMENT_MODE STREQUAL "qualification" AND
   NOT SIRIUS_ALIGNMENT_MODE STREQUAL "release")
    message(FATAL_ERROR
        "SIRIUS_ALIGNMENT_MODE must be development, qualification, or release, not "
        "'${SIRIUS_ALIGNMENT_MODE}'")
endif()

if(SIRIUS_MANDATORY_TESTS AND NOT BUILD_TESTS)
    message(FATAL_ERROR
        "SIRIUS_MANDATORY_TESTS=ON requires BUILD_TESTS=ON. "
        "An operational build cannot omit its executable invariants.")
endif()

include(cmake/sirius_toolchain.cmake)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/scripts/verify-build-policy.py")
foreach(_sirius_policy_option
        BUILD_TESTS
        SIRIUS_MANDATORY_TESTS
        SIRIUS_WERROR
        SIRIUS_REQUIRE_VULKAN_RUNTIME)
    if(${_sirius_policy_option})
        set(${_sirius_policy_option}_POLICY_VALUE true)
    else()
        set(${_sirius_policy_option}_POLICY_VALUE false)
    endif()
endforeach()
if(CMAKE_CONFIGURATION_TYPES)
    string(REPLACE ";" "," SIRIUS_POLICY_CONFIGURATIONS
        "${CMAKE_CONFIGURATION_TYPES}")
else()
    set(SIRIUS_POLICY_CONFIGURATIONS "${CMAKE_BUILD_TYPE}")
endif()
execute_process(
    COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_SOURCE_DIR}/scripts/verify-build-policy.py"
        --phase configure
        --alignment-mode "${SIRIUS_ALIGNMENT_MODE}"
        --build-tests "${BUILD_TESTS_POLICY_VALUE}"
        --mandatory-tests "${SIRIUS_MANDATORY_TESTS_POLICY_VALUE}"
        --warnings-as-errors "${SIRIUS_WERROR_POLICY_VALUE}"
        --require-vulkan-runtime "${SIRIUS_REQUIRE_VULKAN_RUNTIME_POLICY_VALUE}"
        --contract-mode "${SIRIUS_CONTRACT_MODE}"
        --configurations "${SIRIUS_POLICY_CONFIGURATIONS}"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    RESULT_VARIABLE SIRIUS_BUILD_POLICY_RESULT
    OUTPUT_VARIABLE SIRIUS_BUILD_POLICY_STDOUT
    ERROR_VARIABLE SIRIUS_BUILD_POLICY_STDERR)
if(NOT SIRIUS_BUILD_POLICY_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Sirius strict build policy failed:\n"
        "${SIRIUS_BUILD_POLICY_STDOUT}${SIRIUS_BUILD_POLICY_STDERR}")
endif()

include(cmake/sirius_dependencies.cmake)
find_package(Threads REQUIRED)

# External evidence becomes an authority only through one revision-bound
# receipt. Development and qualification builds keep absent domains explicit;
# release builds fail configuration unless the complete clean-revision set verifies. The same
# deterministic receipt is checked on every build and embedded into the binary
# for runtime initialisation.
set(SIRIUS_OPERATING_MODEL
    "${CMAKE_SOURCE_DIR}/tests/operating_model.json")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${SIRIUS_OPERATING_MODEL}")
execute_process(
    COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_SOURCE_DIR}/scripts/verify-operating-model.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    RESULT_VARIABLE SIRIUS_OPERATING_MODEL_RESULT
    OUTPUT_VARIABLE SIRIUS_OPERATING_MODEL_STDOUT
    ERROR_VARIABLE SIRIUS_OPERATING_MODEL_STDERR)
if(NOT SIRIUS_OPERATING_MODEL_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Sirius operating-model configuration failed:\n"
        "${SIRIUS_OPERATING_MODEL_STDOUT}${SIRIUS_OPERATING_MODEL_STDERR}")
endif()
find_package(Git QUIET)
set(SIRIUS_SOURCE_REVISION "0000000000000000000000000000000000000000")
set(SIRIUS_SOURCE_TREE_CLEAN false)
set(SIRIUS_HAS_GIT_WORKTREE false)
if(Git_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set(SIRIUS_HAS_GIT_WORKTREE true)
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        RESULT_VARIABLE SIRIUS_GIT_REVISION_RESULT
        OUTPUT_VARIABLE SIRIUS_SOURCE_REVISION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" status --porcelain --untracked-files=normal
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        RESULT_VARIABLE SIRIUS_GIT_STATUS_RESULT
        OUTPUT_VARIABLE SIRIUS_GIT_STATUS
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    if(NOT SIRIUS_GIT_REVISION_RESULT EQUAL 0 OR
       NOT SIRIUS_GIT_STATUS_RESULT EQUAL 0)
        message(FATAL_ERROR "Could not establish the source revision and cleanliness")
    endif()
    if(SIRIUS_GIT_STATUS STREQUAL "")
        set(SIRIUS_SOURCE_TREE_CLEAN true)
    endif()
elseif(NOT SIRIUS_ALIGNMENT_MODE STREQUAL "development")
    message(FATAL_ERROR
        "Qualification/release alignment requires a Git worktree so evidence can bind "
        "the exact clean source revision")
endif()

set(SIRIUS_ALIGNMENT_RECEIPT
    "${CMAKE_BINARY_DIR}/generated/sirius/alignment_receipt.json")
set(SIRIUS_ALIGNMENT_COMMAND
    "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/verify-alignment.py"
    --mode "${SIRIUS_ALIGNMENT_MODE}"
    --source-revision "${SIRIUS_SOURCE_REVISION}"
    --source-tree-clean "${SIRIUS_SOURCE_TREE_CLEAN}"
    --operating-model "${SIRIUS_OPERATING_MODEL}"
    --output "${SIRIUS_ALIGNMENT_RECEIPT}")
if(SIRIUS_HAS_GIT_WORKTREE)
    list(APPEND SIRIUS_ALIGNMENT_COMMAND --source-root "${CMAKE_SOURCE_DIR}")
endif()
if(NOT SIRIUS_ATTESTATION_ROOT STREQUAL "")
    if(EXISTS "${SIRIUS_ATTESTATION_ROOT}")
        file(GLOB_RECURSE SIRIUS_ATTESTATION_INPUTS CONFIGURE_DEPENDS
            "${SIRIUS_ATTESTATION_ROOT}/*.json")
    endif()
    list(APPEND SIRIUS_ALIGNMENT_COMMAND
        --attestation-root "${SIRIUS_ATTESTATION_ROOT}")
endif()
execute_process(
    COMMAND ${SIRIUS_ALIGNMENT_COMMAND}
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    RESULT_VARIABLE SIRIUS_ALIGNMENT_RESULT
    OUTPUT_VARIABLE SIRIUS_ALIGNMENT_STDOUT
    ERROR_VARIABLE SIRIUS_ALIGNMENT_STDERR)
if(NOT SIRIUS_ALIGNMENT_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Sirius alignment configuration failed:\n"
        "${SIRIUS_ALIGNMENT_STDOUT}${SIRIUS_ALIGNMENT_STDERR}")
endif()
string(STRIP "${SIRIUS_ALIGNMENT_STDOUT}" SIRIUS_ALIGNMENT_SUMMARY)
message(STATUS "${SIRIUS_ALIGNMENT_SUMMARY}")

add_custom_target(SiriusAlignmentGate ALL
    COMMAND ${SIRIUS_ALIGNMENT_COMMAND} --check
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    COMMENT "Verifying revision-bound upstream/downstream alignment")

enable_testing()

if(NOT SIRIUS_MODERN_AVAILABLE)
    message(FATAL_ERROR
        "Sirius requires a C++26-mode compiler: ${SIRIUS_MODERN_UNAVAILABLE_REASON}. "
        "Use a preset from CMakePresets.json to select a supported compiler.")
endif()

add_subdirectory(src/sirius)
set(SIRIUS_VULKAN_BACKEND_POLICY_VALUE false)
set(SIRIUS_COMPILED_KERNELS_POLICY_VALUE false)
set(SIRIUS_SPIRV_VALIDATION_POLICY_VALUE false)
set(SIRIUS_NATIVE_VIEWER_POLICY_VALUE false)
if(TARGET sirius_backend)
    set(SIRIUS_VULKAN_BACKEND_POLICY_VALUE true)
endif()
if(TARGET sirius_kernels)
    set(SIRIUS_COMPILED_KERNELS_POLICY_VALUE true)
endif()
if(SIRIUS_SPIRV_VAL)
    set(SIRIUS_SPIRV_VALIDATION_POLICY_VALUE true)
endif()
if(SIRIUS_INSTALL_HAS_VIEWER)
    set(SIRIUS_NATIVE_VIEWER_POLICY_VALUE true)
endif()
execute_process(
    COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_SOURCE_DIR}/scripts/verify-build-policy.py"
        --phase product
        --alignment-mode "${SIRIUS_ALIGNMENT_MODE}"
        --vulkan-backend "${SIRIUS_VULKAN_BACKEND_POLICY_VALUE}"
        --compiled-kernels "${SIRIUS_COMPILED_KERNELS_POLICY_VALUE}"
        --spirv-validation "${SIRIUS_SPIRV_VALIDATION_POLICY_VALUE}"
        --native-viewer "${SIRIUS_NATIVE_VIEWER_POLICY_VALUE}"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    RESULT_VARIABLE SIRIUS_PRODUCT_POLICY_RESULT
    OUTPUT_VARIABLE SIRIUS_PRODUCT_POLICY_STDOUT
    ERROR_VARIABLE SIRIUS_PRODUCT_POLICY_STDERR)
if(NOT SIRIUS_PRODUCT_POLICY_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Sirius strict product policy failed:\n"
        "${SIRIUS_PRODUCT_POLICY_STDOUT}${SIRIUS_PRODUCT_POLICY_STDERR}")
endif()
if(SIRIUS_REQUIRE_VULKAN_RUNTIME AND
   (NOT TARGET sirius_backend OR NOT TARGET sirius_kernels))
    message(FATAL_ERROR
        "SIRIUS_REQUIRE_VULKAN_RUNTIME=ON requires the Vulkan SDK and slangc. "
        "The required runtime profile may not degrade to a CPU-only build.")
endif()
if(BUILD_TESTS)
    FetchContent_MakeAvailable(googletest)
    add_subdirectory(tests)
    add_custom_target(SiriusSourceGovernance ALL
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/generate-ctest-labels.py"
            --check
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/verify-operating-model.py"
        COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_SOURCE_DIR}/scripts/verify-repository-structure.py"
        COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_SOURCE_DIR}/scripts/verify-build-policy.py" --self-test
        COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_SOURCE_DIR}/scripts/verify-build-gate.py" --self-test
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        COMMENT "Verifying test-label, operating-model, and repository-structure governance")
endif()

# === Mandatory test gating (the build gate; label authority is generated) ===
if(BUILD_TESTS AND SIRIUS_MANDATORY_TESTS)
    set(sirius_mandatory_test_targets)
    foreach(test_target
            sirius_base_tests
            sirius_core_tests
            sirius_oracle_tests
            sirius_backend_tests
            sirius_render_tests
            sirius_app_tests)
        if(TARGET ${test_target})
            list(APPEND sirius_mandatory_test_targets ${test_target})
        endif()
    endforeach()

    set(SIRIUS_MANDATORY_GATE_STAMP
        "${CMAKE_BINARY_DIR}/generated/sirius/mandatory_gate.json")
    set(sirius_gate_test_artifacts
        --tested-artifact "sirius=$<TARGET_FILE:sirius>")
    foreach(test_target IN LISTS sirius_mandatory_test_targets)
        list(APPEND sirius_gate_test_artifacts
            --tested-artifact "${test_target}=$<TARGET_FILE:${test_target}>")
    endforeach()
    set(sirius_gate_product_artifacts
        --product-artifact "sirius=$<TARGET_FILE:sirius>"
        --product-artifact "starfield=${CMAKE_SOURCE_DIR}/assets/Starfield.png"
        --product-artifact "operating_model=${SIRIUS_OPERATING_MODEL}"
        --product-artifact "alignment_receipt=${SIRIUS_ALIGNMENT_RECEIPT}")
    if(SIRIUS_INSTALL_HAS_VIEWER)
        list(APPEND sirius_gate_product_artifacts
            --product-artifact
                "viewer_rdsd003a_vertex=${CMAKE_SOURCE_DIR}/src/sirius/app/viewer/shaders/RDSD003A.vert"
            --product-artifact
                "viewer_rdsd003a_fragment=${CMAKE_SOURCE_DIR}/src/sirius/app/viewer/shaders/RDSD003A.frag"
            --product-artifact
                "viewer_rdsd004a_vertex=${CMAKE_SOURCE_DIR}/src/sirius/app/viewer/shaders/RDSD004A.vert"
            --product-artifact
                "viewer_rdsd004a_fragment=${CMAKE_SOURCE_DIR}/src/sirius/app/viewer/shaders/RDSD004A.frag"
            --product-artifact
                "viewer_rdsd005a_vertex=${CMAKE_SOURCE_DIR}/src/sirius/app/viewer/shaders/RDSD005A.vert"
            --product-artifact
                "viewer_rdsd005a_fragment=${CMAKE_SOURCE_DIR}/src/sirius/app/viewer/shaders/RDSD005A.frag")
    endif()
    if(TARGET sirius_kernels)
        list(APPEND sirius_gate_product_artifacts
            --product-artifact "trace_spv=${SIRIUS_KERNEL_BINARY_DIR}/trace.spv"
            --product-artifact
                "trace_fp32comp_spv=${SIRIUS_KERNEL_BINARY_DIR}/trace_fp32comp.spv"
            --product-artifact "trace_fp64_spv=${SIRIUS_KERNEL_BINARY_DIR}/trace_fp64.spv")
    endif()

    set(sirius_mandatory_gate_command
        "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/verify-build-gate.py"
        --action run
        --stamp "${SIRIUS_MANDATORY_GATE_STAMP}"
        --alignment-mode "${SIRIUS_ALIGNMENT_MODE}"
        --source-root "${CMAKE_SOURCE_DIR}"
        --build-dir "${CMAKE_BINARY_DIR}"
        --source-revision "${SIRIUS_SOURCE_REVISION}"
        --source-tree-clean "${SIRIUS_SOURCE_TREE_CLEAN}"
        --ctest "${CMAKE_CTEST_COMMAND}"
        --config "$<CONFIG>"
        ${sirius_gate_test_artifacts}
        ${sirius_gate_product_artifacts})
    if(NOT SIRIUS_SANITIZERS STREQUAL "none")
        set(sirius_mandatory_gate_command
            "${CMAKE_COMMAND}" -E env
            "ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:strict_string_checks=1"
            "LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/sanitizers/lsan-vulkan.supp:print_suppressions=1"
            "UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1"
            ${sirius_mandatory_gate_command})
    endif()

    add_custom_target(RunMandatoryTests ALL
        COMMAND "${CMAKE_COMMAND}" -E rm -f
            "$<TARGET_FILE_DIR:sirius>/resources/model/mandatory_gate.json"
        COMMAND ${sirius_mandatory_gate_command}
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different
            "${SIRIUS_MANDATORY_GATE_STAMP}"
            "$<TARGET_FILE_DIR:sirius>/resources/model/mandatory_gate.json"
        # gtest_discover_tests writes each suite's CTest include after linking.
        # Depend on every available suite so the gate cannot race discovery and
        # incorrectly report that no Mandatory tests exist.
        DEPENDS sirius ${sirius_mandatory_test_targets}
            SiriusAlignmentGate SiriusSourceGovernance
        BYPRODUCTS
            "${SIRIUS_MANDATORY_GATE_STAMP}"
            "${CMAKE_BINARY_DIR}/generated/sirius/mandatory_gate_junit.xml"
            "${CMAKE_BINARY_DIR}/generated/sirius/mandatory_gate_ctest.log"
        COMMENT "=== MANDATORY TEST GATE ==="
    )
    message(STATUS "Mandatory test gating: ENABLED")
endif()

# === Install and packaging ===
if(TARGET sirius)
    install(TARGETS sirius RUNTIME DESTINATION bin)
    install(FILES assets/Starfield.png DESTINATION share/sirius/assets)
    install(FILES tests/operating_model.json DESTINATION share/sirius/model)
    install(FILES "${SIRIUS_ALIGNMENT_RECEIPT}"
        DESTINATION share/sirius/model)
    if(NOT SIRIUS_ALIGNMENT_MODE STREQUAL "development")
        install(FILES "${SIRIUS_MANDATORY_GATE_STAMP}"
            DESTINATION share/sirius/model)
    endif()
    if(SIRIUS_INSTALL_HAS_VIEWER)
        install(DIRECTORY src/sirius/app/viewer/shaders/
            DESTINATION share/sirius/shaders
            FILES_MATCHING PATTERN "*.vert" PATTERN "*.frag")
    endif()
    if(TARGET sirius_kernels)
        install(FILES
            "${SIRIUS_KERNEL_BINARY_DIR}/trace.spv"
            "${SIRIUS_KERNEL_BINARY_DIR}/trace_fp32comp.spv"
            "${SIRIUS_KERNEL_BINARY_DIR}/trace_fp64.spv"
            DESTINATION share/sirius/kernels)
        set(SIRIUS_INSTALL_HAS_VULKAN TRUE)
    else()
        set(SIRIUS_INSTALL_HAS_VULKAN FALSE)
    endif()
    configure_file(
        cmake/verify_sirius_install.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/verify_sirius_install.cmake"
        @ONLY)
    install(SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/verify_sirius_install.cmake")
endif()

set(CPACK_PACKAGE_NAME "sirius")
set(CPACK_PACKAGE_VENDOR "Sirius")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "General-relativistic ray tracing engine")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
if(SIRIUS_ALIGNMENT_MODE STREQUAL "release")
    include(CPack)
else()
    message(STATUS
        "Release packaging: DISABLED (configure SIRIUS_ALIGNMENT_MODE=release "
        "with the complete attestation set)")
endif()
