execute_process(
    COMMAND_ECHO STDOUT
    ECHO_OUTPUT_VARIABLE
    ECHO_ERROR_VARIABLE
    COMMAND "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_LIST_DIR}/early-child.cmake"
    RESULT_VARIABLE child_result OUTPUT_VARIABLE child_out ERROR_VARIABLE child_err)
message(FATAL_ERROR "stop diagnostic child unexpectedly returned: ${child_result}")
