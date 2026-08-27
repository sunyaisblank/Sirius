if(NOT DEFINED SIRIUS_BINARY OR NOT EXISTS "${SIRIUS_BINARY}")
    message(FATAL_ERROR "required Vulkan profile did not receive a Sirius executable")
endif()
if(NOT DEFINED SIRIUS_OUTPUT)
    message(FATAL_ERROR "required Vulkan profile did not receive an output path")
endif()

execute_process(
    COMMAND "${SIRIUS_BINARY}" --json info readiness
    RESULT_VARIABLE _readiness_result
    OUTPUT_VARIABLE _readiness_stdout
    ERROR_VARIABLE _readiness_stderr)
string(JSON _evidence_ready ERROR_VARIABLE _evidence_error
    GET "${_readiness_stdout}" evidence_generation_ready)
string(JSON _vulkan_ready ERROR_VARIABLE _vulkan_error
    GET "${_readiness_stdout}" backends vulkan ready)
if((NOT _readiness_result EQUAL 0 AND NOT _readiness_result EQUAL 1) OR
   NOT _evidence_error STREQUAL "NOTFOUND" OR
   NOT _vulkan_error STREQUAL "NOTFOUND" OR
   NOT _evidence_ready OR NOT _vulkan_ready)
    message(FATAL_ERROR
        "required Vulkan runtime is not ready (${_readiness_result})\n"
        "${_readiness_stdout}\n${_readiness_stderr}")
endif()

file(REMOVE "${SIRIUS_OUTPUT}")
execute_process(
    COMMAND "${SIRIUS_BINARY}" render
        --backend vulkan
        --metric Kerr --spin 0.9
        --width 128 --height 128 --samples 3
        --temperature-model SS
        --doppler-beaming off
        --camera-beta 0.08,0.02,0.01
        --lens ThinLens --focal-length 50 --aperture 2.8 --focus-distance 30
        --beams --starfield point
        --no-bloom --output "${SIRIUS_OUTPUT}"
    RESULT_VARIABLE _render_result
    OUTPUT_VARIABLE _render_stdout
    ERROR_VARIABLE _render_stderr)
if(NOT _render_result EQUAL 0 OR NOT EXISTS "${SIRIUS_OUTPUT}")
    message(FATAL_ERROR
        "required Vulkan dispatch failed (${_render_result})\n"
        "${_render_stdout}\n${_render_stderr}")
endif()
file(SIZE "${SIRIUS_OUTPUT}" _render_size)
if(_render_size LESS 1024)
    message(FATAL_ERROR
        "required Vulkan dispatch produced an implausibly small output: ${_render_size} bytes")
endif()
file(REMOVE "${SIRIUS_OUTPUT}")
