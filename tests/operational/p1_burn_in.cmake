if(NOT DEFINED SIRIUS_RENDER_TESTS OR NOT EXISTS "${SIRIUS_RENDER_TESTS}")
    message(FATAL_ERROR "P1 burn-in did not receive the render-test executable")
endif()

# Repeat both independent classifiers in one process and recreate every fixture
# between iterations. This is a stability witness in addition to the ordinary
# one-shot correctness tests: device/session teardown, deterministic boundary
# classification, and the near-extremal path are exercised repeatedly.
execute_process(
    COMMAND "${SIRIUS_RENDER_TESTS}"
        "--gtest_filter=ShadowBoundary.KerrNearExtremalMatchesBardeenWithinOnePixelAt1080p:VulkanRenderSession.KerrNearExtremalBardeenBoundaryAt1080p"
        --gtest_repeat=3
        --gtest_recreate_environments_when_repeating
        --gtest_break_on_failure
    RESULT_VARIABLE _burn_in_result
    OUTPUT_VARIABLE _burn_in_stdout
    ERROR_VARIABLE _burn_in_stderr)
if(NOT _burn_in_result EQUAL 0)
    message(FATAL_ERROR
        "P1 near-extremal burn-in failed (${_burn_in_result})\n"
        "${_burn_in_stdout}\n${_burn_in_stderr}")
endif()
