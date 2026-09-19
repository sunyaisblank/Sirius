# Configure-time checks before and after product target creation.

function(sirius_verify_configure_policy)
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
endfunction()

function(sirius_verify_product_policy)
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
endfunction()
