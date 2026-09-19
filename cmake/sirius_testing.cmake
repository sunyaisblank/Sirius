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
                "viewer_rdsd003a_fragment=${CMAKE_SOURCE_DIR}/src/sirius/app/viewer/shaders/RDSD003A.frag")
    endif()
    set(sirius_gate_test_input_artifacts)
    if(TARGET sirius_retained_camera_test_inputs)
        get_target_property(retained_camera_inputs sirius_retained_camera_test_inputs
            SIRIUS_TEST_INPUT_ARTIFACTS)
        foreach(retained_camera_input IN LISTS retained_camera_inputs)
            list(APPEND sirius_gate_test_input_artifacts
                --test-input-artifact "${retained_camera_input}")
        endforeach()
    endif()
    if(TARGET sirius_kernels)
        list(APPEND sirius_gate_product_artifacts
            --product-artifact "trace_spv=${SIRIUS_KERNEL_BINARY_DIR}/trace.spv"
            --product-artifact
                "trace_fp32comp_spv=${SIRIUS_KERNEL_BINARY_DIR}/trace_fp32comp.spv"
            --product-artifact "trace_fp64_spv=${SIRIUS_KERNEL_BINARY_DIR}/trace_fp64.spv")
        # These generated inputs are consumed by Mandatory tests but are not
        # installed runtime products. Bind their actual test paths separately.
        list(APPEND sirius_gate_test_input_artifacts
            --test-input-artifact "smoke_spv=${SIRIUS_KERNEL_BINARY_DIR}/smoke.spv"
            --test-input-artifact "parity_probe_spv=${SIRIUS_KERNEL_BINARY_DIR}/parity_probe.spv"
            --test-input-artifact
                "parity_probe_fp32comp_spv=${SIRIUS_KERNEL_BINARY_DIR}/parity_probe_fp32comp.spv"
            --test-input-artifact
                "parity_probe_fp64_spv=${SIRIUS_KERNEL_BINARY_DIR}/parity_probe_fp64.spv"
            --test-input-artifact "infinity_probe_spv=${SIRIUS_KERNEL_BINARY_DIR}/infinity_probe.spv"
            --test-input-artifact
                "infinity_probe_fp32comp_spv=${SIRIUS_KERNEL_BINARY_DIR}/infinity_probe_fp32comp.spv"
            --test-input-artifact
                "infinity_probe_fp64_spv=${SIRIUS_KERNEL_BINARY_DIR}/infinity_probe_fp64.spv"
            --test-input-artifact
                "metric_consistency_probe_spv=${SIRIUS_KERNEL_BINARY_DIR}/metric_consistency_probe.spv"
            --test-input-artifact
                "metric_consistency_probe_fp32comp_spv=${SIRIUS_KERNEL_BINARY_DIR}/metric_consistency_probe_fp32comp.spv"
            --test-input-artifact
                "metric_consistency_probe_fp64_spv=${SIRIUS_KERNEL_BINARY_DIR}/metric_consistency_probe_fp64.spv"
            --test-input-artifact
                "camera_frame_probe_spv=${SIRIUS_KERNEL_BINARY_DIR}/camera_frame_probe.spv"
            --test-input-artifact
                "camera_frame_probe_fp32comp_spv=${SIRIUS_KERNEL_BINARY_DIR}/camera_frame_probe_fp32comp.spv"
            --test-input-artifact
                "camera_frame_probe_fp64_spv=${SIRIUS_KERNEL_BINARY_DIR}/camera_frame_probe_fp64.spv"
            --test-input-artifact
                "coupled_probe_fp64_spv=${SIRIUS_KERNEL_BINARY_DIR}/coupled_probe_fp64.spv"
            --test-input-artifact "trace_cuda=${SIRIUS_KERNEL_BINARY_DIR}/portability/trace.cu"
            --test-input-artifact "trace_metal=${SIRIUS_KERNEL_BINARY_DIR}/portability/trace.metal")
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
        ${sirius_gate_product_artifacts}
        ${sirius_gate_test_input_artifacts})
    if(NOT SIRIUS_SANITIZERS STREQUAL "none")
        set(sirius_mandatory_gate_command
            "${CMAKE_COMMAND}" -E env
            "ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:strict_string_checks=1"
            "LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/tests/sanitizers/lsan-vulkan.supp:print_suppressions=1"
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

    # Build-domain evidence is intentionally distinct from the promotion gate.
    # It hashes the complete strict product/test topology and runs only the
    # fixed non-render authority selection. It is never part of ALL, never
    # staged into the runtime volume, and cannot satisfy a runtime/image domain.
    set(SIRIUS_NATIVE_BUILD_GATE_STAMP
        "${CMAKE_BINARY_DIR}/generated/sirius/native_build_gate.json")
    add_custom_target(RunNativeBuildEvidence
        COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_SOURCE_DIR}/scripts/verify-build-gate.py"
            --action run-native-build
            --stamp "${SIRIUS_NATIVE_BUILD_GATE_STAMP}"
            --alignment-mode "${SIRIUS_ALIGNMENT_MODE}"
            --source-root "${CMAKE_SOURCE_DIR}"
            --build-dir "${CMAKE_BINARY_DIR}"
            --source-revision "${SIRIUS_SOURCE_REVISION}"
            --source-tree-clean "${SIRIUS_SOURCE_TREE_CLEAN}"
            --ctest "${CMAKE_CTEST_COMMAND}"
            --config "$<CONFIG>"
            ${sirius_gate_test_artifacts}
            ${sirius_gate_product_artifacts}
            ${sirius_gate_test_input_artifacts}
        DEPENDS sirius ${sirius_mandatory_test_targets}
            SiriusAlignmentGate SiriusSourceGovernance
        BYPRODUCTS
            "${SIRIUS_NATIVE_BUILD_GATE_STAMP}"
            "${CMAKE_BINARY_DIR}/generated/sirius/native_build_gate_junit.xml"
            "${CMAKE_BINARY_DIR}/generated/sirius/native_build_gate_ctest.log"
        COMMENT "=== NATIVE BUILD NON-RENDER EVIDENCE GATE ==="
    )
    message(STATUS "Mandatory test gating: ENABLED")
endif()
