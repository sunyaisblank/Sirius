execute_process(COMMAND "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_LIST_DIR}/child.cmake"
    RESULT_VARIABLE baseline_result OUTPUT_VARIABLE baseline_out ERROR_VARIABLE baseline_err)
execute_process(
    COMMAND_ECHO STDOUT
    ECHO_OUTPUT_VARIABLE
    ECHO_ERROR_VARIABLE
    COMMAND "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_LIST_DIR}/child.cmake"
    RESULT_VARIABLE mirrored_result OUTPUT_VARIABLE mirrored_out ERROR_VARIABLE mirrored_err)
if(NOT baseline_result EQUAL 0 OR NOT mirrored_result EQUAL 0 OR
   NOT baseline_out STREQUAL mirrored_out OR NOT baseline_err STREQUAL mirrored_err OR
   NOT mirrored_out STREQUAL "-- stdout: a;b [quoted]\n" OR
   NOT mirrored_err STREQUAL "stderr: c;d [quoted]\n")
    message(FATAL_ERROR "captured stdout/stderr or result changed")
endif()
message(STATUS "CAPTURE_UNCHANGED")
