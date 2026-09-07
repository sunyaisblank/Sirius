if(NOT DEFINED SIRIUS_RENDER_TESTS OR NOT EXISTS "${SIRIUS_RENDER_TESTS}")
    message(FATAL_ERROR "required Vulkan profile did not receive the render test executable")
endif()

execute_process(
    COMMAND "${SIRIUS_RENDER_TESTS}"
        --gtest_filter=RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent
        --gtest_color=no
    RESULT_VARIABLE _render_result
    OUTPUT_VARIABLE _render_stdout
    ERROR_VARIABLE _render_stderr)
if(NOT _render_result EQUAL 0)
    message(FATAL_ERROR
        "required direct Vulkan dispatch failed (${_render_result})\n"
        "${_render_stdout}\n${_render_stderr}")
endif()
