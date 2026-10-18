# Private compile-only hook; pass as CMAKE_PROJECT_Sirius_INCLUDE on Windows.
# Defer until canonical libraries, fixtures, policies and tests exist.
include_guard(GLOBAL)

function(_sirius_native_rounded_error_matched_producer)
    if(NOT WIN32 OR NOT MSVC OR NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
        message(FATAL_ERROR "Private matched producer requires canonical Windows MSVC x64")
    endif()
    if(NOT BUILD_TESTS OR NOT SIRIUS_MANDATORY_TESTS OR NOT SIRIUS_WERROR OR
       NOT SIRIUS_REQUIRE_VULKAN_RUNTIME OR NOT SIRIUS_CONTRACT_MODE STREQUAL "2")
        message(FATAL_ERROR "Private producer may not weaken canonical required-runtime policy")
    endif()
    foreach(required sirius_backend_cpu sirius_backend sirius_render sirius_base
                     GTest::gtest_main sirius_kernels SiriusSourceGovernance
                     sirius_portable_binary32_test_inputs sirius_retained_camera_test_inputs)
        if(NOT TARGET ${required})
            message(FATAL_ERROR "Missing canonical producer dependency: ${required}")
        endif()
    endforeach()
    set(inputs "${CMAKE_SOURCE_DIR}/out/native-rounded-error-producer-inputs")
    set(observer "${inputs}/observer/audited-observer.cpp")
    if(NOT EXISTS "${observer}")
        message(FATAL_ERROR "ROOT must first materialize the admitted unchanged observer")
    endif()
    file(SHA256 "${observer}" observer_sha)
    if(NOT observer_sha STREQUAL "c9485d03ec55b5c1307a34097f15efbc58c4f84749e9f3736ad063a68676b601")
        message(FATAL_ERROR "Observer differs from the exact f851 audited source")
    endif()
    foreach(arm baseline candidate)
        foreach(header retained_kernels.h retained_owner_kernels.h retained_owner_layout.h)
            if(NOT EXISTS "${inputs}/${arm}/${header}")
                message(FATAL_ERROR "Missing ROOT-admitted private header: ${arm}/${header}")
            endif()
        endforeach()
    endforeach()
    foreach(header retained_owner_kernels.h retained_owner_layout.h)
        file(SHA256 "${inputs}/baseline/${header}" baseline_sha)
        file(SHA256 "${inputs}/candidate/${header}" candidate_sha)
        if(NOT baseline_sha STREQUAL candidate_sha)
            message(FATAL_ERROR "Owner header must be byte-identical between both arms")
        endif()
    endforeach()
    set(fixture "${CMAKE_BINARY_DIR}/tests/backend/retained_camera")
    set(canonical_retained "${CMAKE_BINARY_DIR}/src/sirius/backend/retained")

    # Common DP source does not include retained_kernels.h: compile it once.
    add_library(native_rounded_error_dopri OBJECT
        "${CMAKE_SOURCE_DIR}/tests/backend/retained_dopri_test.cpp")
    sirius_configure_target(native_rounded_error_dopri CONTRACT_MODE 2)
    target_link_libraries(native_rounded_error_dopri PRIVATE
        sirius_backend_cpu sirius_backend sirius_render GTest::gtest_main)
    target_include_directories(native_rounded_error_dopri PRIVATE
        "${CMAKE_SOURCE_DIR}/tests" "${canonical_retained}" "${fixture}")
    target_compile_definitions(native_rounded_error_dopri PRIVATE
        SIRIUS_TEST_HAS_KERNELS=1 SIRIUS_TEST_HAS_PORTABLE_BINARY32_REFERENCE=1
        SIRIUS_TEST_HAS_RETAINED_CAMERA=1)
    target_compile_options(native_rounded_error_dopri PRIVATE
        "/FI${CMAKE_SOURCE_DIR}/tests/support/fail_on_skip.h")
    add_dependencies(native_rounded_error_dopri SiriusSourceGovernance
        sirius_portable_binary32_test_inputs sirius_retained_camera_test_inputs)

    foreach(arm baseline candidate)
        # Keep compute's compile usage requirements at canonical backend scope.
        set(compute native_rounded_error_${arm}_compute)
        add_library(${compute} OBJECT
            "${CMAKE_SOURCE_DIR}/src/sirius/backend/retained_compute.cpp")
        sirius_configure_target(${compute} CONTRACT_MODE 2)
        target_link_libraries(${compute} PRIVATE sirius_backend)
        target_include_directories(${compute} BEFORE PRIVATE
            "${inputs}/${arm}" "${canonical_retained}")
        add_dependencies(${compute} SiriusSourceGovernance)

        set(exe native_rounded_error_${arm}_controls)
        add_executable(${exe} "${observer}"
            $<TARGET_OBJECTS:${compute}> $<TARGET_OBJECTS:native_rounded_error_dopri>)
        sirius_configure_target(${exe} CONTRACT_MODE 2)
        target_link_libraries(${exe} PRIVATE
            sirius_backend_cpu sirius_backend sirius_render GTest::gtest_main)
        target_include_directories(${exe} BEFORE PRIVATE
            "${inputs}/${arm}" "${CMAKE_SOURCE_DIR}/tests"
            "${canonical_retained}" "${fixture}")
        target_compile_definitions(${exe} PRIVATE
            SIRIUS_TEST_HAS_KERNELS=1 SIRIUS_TEST_HAS_PORTABLE_BINARY32_REFERENCE=1
            SIRIUS_TEST_HAS_RETAINED_CAMERA=1)
        target_compile_options(${exe} PRIVATE
            "/FI${CMAKE_SOURCE_DIR}/tests/support/fail_on_skip.h")
        set(output "${CMAKE_BINARY_DIR}/private-native-rounded-error/${arm}/$<CONFIG>")
        set_target_properties(${exe} PROPERTIES
            OUTPUT_NAME original_${arm}_controls
            RUNTIME_OUTPUT_DIRECTORY "${output}" PDB_OUTPUT_DIRECTORY "${output}")
        target_link_options(${exe} PRIVATE
            "/MAP:${output}/original_${arm}_controls.map" /VERBOSE:LIB)
        add_dependencies(${exe} SiriusSourceGovernance
            sirius_portable_binary32_test_inputs sirius_retained_camera_test_inputs)
    endforeach()
    # No ALL, gtest_discover_tests, add_test, provider query or runtime command.
    add_custom_target(SiriusNativeRoundedErrorMatchedProducer DEPENDS
        native_rounded_error_baseline_controls native_rounded_error_candidate_controls)
endfunction()

cmake_language(DEFER DIRECTORY "${CMAKE_SOURCE_DIR}"
    CALL _sirius_native_rounded_error_matched_producer)
