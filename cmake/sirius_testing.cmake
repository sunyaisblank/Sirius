if(BUILD_TESTS)
    FetchContent_MakeAvailable(googletest)
    add_subdirectory(tests)
    set(sirius_generated_test_inputs)
    foreach(input_target sirius_portable_binary32_test_inputs sirius_retained_camera_test_inputs)
        if(TARGET ${input_target})
            get_target_property(generated_inputs ${input_target} SIRIUS_TEST_INPUT_ARTIFACTS)
            foreach(generated_input IN LISTS generated_inputs)
                list(APPEND sirius_generated_test_inputs
                    "${generated_input}")
            endforeach()
        endif()
    endforeach()
    if(TARGET sirius_kernels)
        # Preserve canonical generated artifact identities; stage their exact bytes
        # beside both consumers below.
        list(APPEND sirius_generated_test_inputs
            "smoke_spv=${SIRIUS_KERNEL_BINARY_DIR}/smoke.spv"
            "parity_probe_spv=${SIRIUS_KERNEL_BINARY_DIR}/parity_probe.spv"
                "parity_probe_fp32comp_spv=${SIRIUS_KERNEL_BINARY_DIR}/parity_probe_fp32comp.spv"
                "parity_probe_fp64_spv=${SIRIUS_KERNEL_BINARY_DIR}/parity_probe_fp64.spv"
            "infinity_probe_spv=${SIRIUS_KERNEL_BINARY_DIR}/infinity_probe.spv"
                "infinity_probe_fp32comp_spv=${SIRIUS_KERNEL_BINARY_DIR}/infinity_probe_fp32comp.spv"
                "infinity_probe_fp64_spv=${SIRIUS_KERNEL_BINARY_DIR}/infinity_probe_fp64.spv"
                "metric_consistency_probe_spv=${SIRIUS_KERNEL_BINARY_DIR}/metric_consistency_probe.spv"
                "metric_consistency_probe_fp32comp_spv=${SIRIUS_KERNEL_BINARY_DIR}/metric_consistency_probe_fp32comp.spv"
                "metric_consistency_probe_fp64_spv=${SIRIUS_KERNEL_BINARY_DIR}/metric_consistency_probe_fp64.spv"
                "camera_frame_probe_spv=${SIRIUS_KERNEL_BINARY_DIR}/camera_frame_probe.spv"
                "camera_frame_probe_fp32comp_spv=${SIRIUS_KERNEL_BINARY_DIR}/camera_frame_probe_fp32comp.spv"
                "camera_frame_probe_fp64_spv=${SIRIUS_KERNEL_BINARY_DIR}/camera_frame_probe_fp64.spv"
                "coupled_probe_fp64_spv=${SIRIUS_KERNEL_BINARY_DIR}/coupled_probe_fp64.spv"
            "trace_cuda=${SIRIUS_KERNEL_BINARY_DIR}/portability/trace.cu"
            "trace_metal=${SIRIUS_KERNEL_BINARY_DIR}/portability/trace.metal")
    endif()


    # Refresh input volumes on every consumer build, even when a regenerated
    # shader does not require relinking the executable. The gate checks these
    # consumed copies against the canonical input/product records pre/post CTest.
    set(sirius_test_input_dependencies sirius_portable_binary32_test_inputs)
    foreach(input_target sirius_retained_camera_test_inputs sirius_kernels)
        if(TARGET ${input_target})
            list(APPEND sirius_test_input_dependencies ${input_target})
        endif()
    endforeach()
    foreach(test_target sirius_backend_tests sirius_render_tests)
        set(volume_commands)
        foreach(generated_input IN LISTS sirius_generated_test_inputs)
            string(REGEX REPLACE "^[^=]+=" "" input_path "${generated_input}")
            file(RELATIVE_PATH input_relative "${CMAKE_BINARY_DIR}" "${input_path}")
            get_filename_component(input_directory "${input_relative}" DIRECTORY)
            list(APPEND volume_commands
                COMMAND "${CMAKE_COMMAND}" -E make_directory
                    "$<TARGET_FILE_DIR:${test_target}>/resources/${input_directory}"
                COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${input_path}"
                    "$<TARGET_FILE_DIR:${test_target}>/resources/${input_relative}")
        endforeach()
        if(TARGET sirius_kernels)
            list(APPEND volume_commands
                COMMAND "${CMAKE_COMMAND}" -E make_directory
                    "$<TARGET_FILE_DIR:${test_target}>/resources/kernels"
                COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                    "${SIRIUS_KERNEL_BINARY_DIR}/trace.spv"
                    "${SIRIUS_KERNEL_BINARY_DIR}/trace_fp32comp.spv"
                    "${SIRIUS_KERNEL_BINARY_DIR}/trace_fp64.spv"
                    "$<TARGET_FILE_DIR:${test_target}>/resources/kernels")
        endif()
        add_custom_target(${test_target}_input_volume
            ${volume_commands}
            DEPENDS ${sirius_test_input_dependencies}
            COMMENT "Staging exact qualification inputs beside ${test_target}"
            VERBATIM)
        add_dependencies(${test_target} ${test_target}_input_volume)
    endforeach()
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
    foreach(generated_input IN LISTS sirius_generated_test_inputs)
        list(APPEND sirius_gate_test_input_artifacts
            --test-input-artifact "${generated_input}")
    endforeach()
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
        ${sirius_gate_product_artifacts}
        ${sirius_gate_test_input_artifacts})
    if(SIRIUS_REQUIRE_VULKAN_RUNTIME)
        set(sirius_identity_require_vulkan true)
    else()
        set(sirius_identity_require_vulkan false)
    endif()
    set(sirius_runtime_identity_command
        "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/runtime_identity.py"
        --source-root "${CMAKE_SOURCE_DIR}"
        --source-revision "${SIRIUS_SOURCE_REVISION}"
        --source-tree-clean "${SIRIUS_SOURCE_TREE_CLEAN}"
        --executable "$<TARGET_FILE:sirius>"
        --require-vulkan "${sirius_identity_require_vulkan}"
        --output "${CMAKE_BINARY_DIR}/generated/sirius/mandatory_runtime_identity_gate.json")
    if(NOT SIRIUS_SANITIZERS STREQUAL "none")
        set(sirius_mandatory_gate_command
            "${CMAKE_COMMAND}" -E env
            "ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:strict_string_checks=1"
            "LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/tests/sanitizers/lsan-vulkan.supp:print_suppressions=1"
            "UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1"
            ${sirius_mandatory_gate_command})
        set(sirius_runtime_identity_command
            "${CMAKE_COMMAND}" -E env
            "ASAN_OPTIONS=detect_leaks=1:halt_on_error=1:strict_string_checks=1"
            "LSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/tests/sanitizers/lsan-vulkan.supp:print_suppressions=1"
            "UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1"
            ${sirius_runtime_identity_command})
    endif()

    add_custom_target(RunMandatoryTests ALL
        COMMAND "${CMAKE_COMMAND}" -E rm -f
            "${SIRIUS_MANDATORY_GATE_STAMP}"
            "$<TARGET_FILE_DIR:sirius>/resources/model/mandatory_gate.json"
        COMMAND ${sirius_runtime_identity_command}
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
            "${CMAKE_BINARY_DIR}/generated/sirius/mandatory_runtime_identity_gate.json"
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
