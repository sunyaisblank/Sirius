message(STATUS "EARLY_STDOUT")
message("EARLY_STDERR")
execute_process(COMMAND "${CMAKE_COMMAND}" -E sleep 30)
