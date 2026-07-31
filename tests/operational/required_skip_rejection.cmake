if(NOT DEFINED SIRIUS_BACKEND_TESTS OR NOT EXISTS "${SIRIUS_BACKEND_TESTS}")
    message(FATAL_ERROR "required skip-rejection probe has no backend test binary")
endif()

# Remove every Vulkan ICD from this child process. The selected test would use
# GTEST_SKIP() in a portable profile; the required profile must compile that
# same branch into a failure.
execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env
        "VK_ICD_FILENAMES=/sirius/intentional/missing-vulkan-icd.json"
        "VK_DRIVER_FILES=/sirius/intentional/missing-vulkan-driver.json"
        "VK_ADD_DRIVER_FILES="
        "${SIRIUS_BACKEND_TESTS}"
        "--gtest_filter=VulkanBackend.EnumerationReportsInsteadOfThrowing"
    RESULT_VARIABLE probe_result
    OUTPUT_VARIABLE probe_stdout
    ERROR_VARIABLE probe_stderr)

set(probe_output "${probe_stdout}\n${probe_stderr}")
if(probe_result EQUAL 0)
    message(FATAL_ERROR
        "required profile accepted a GoogleTest capability skip:\n${probe_output}")
endif()
if(probe_output MATCHES "\\[  SKIPPED \\]")
    message(FATAL_ERROR
        "required profile still emitted GoogleTest's successful skip marker:\n${probe_output}")
endif()
if(NOT probe_output MATCHES "\\[  FAILED  \\]")
    message(FATAL_ERROR
        "skip-rejection probe failed for an unrelated reason:\n${probe_output}")
endif()

message(STATUS "required profile converted the missing-device skip into a test failure")
