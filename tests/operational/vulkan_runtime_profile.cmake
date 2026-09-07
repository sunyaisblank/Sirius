if(NOT DEFINED SIRIUS_RENDER_TESTS OR NOT EXISTS "${SIRIUS_RENDER_TESTS}")
    message(FATAL_ERROR "required Vulkan profile did not receive the render test executable")
endif()

execute_process(
    # Disable Mesa's disk cache to cover cold startup, including software
    # compilation deferred until submission. This also covers cold Dozen on
    # physical runs; real output must succeed under the unchanged ray limits.
    COMMAND "${CMAKE_COMMAND}" -E env MESA_SHADER_CACHE_DISABLE=true "${SIRIUS_RENDER_TESTS}"
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
