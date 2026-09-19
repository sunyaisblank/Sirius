# External evidence becomes an authority only through one revision-bound
# receipt. Development and qualification builds keep absent domains explicit;
# release builds fail configuration unless the complete clean-revision set verifies. The same
# deterministic receipt is checked on every build and embedded into the binary
# for runtime initialisation.
set(SIRIUS_OPERATING_MODEL
    "${CMAKE_SOURCE_DIR}/tests/operating_model.json")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${SIRIUS_OPERATING_MODEL}")
execute_process(
    COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_SOURCE_DIR}/scripts/verify-operating-model.py"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    RESULT_VARIABLE SIRIUS_OPERATING_MODEL_RESULT
    OUTPUT_VARIABLE SIRIUS_OPERATING_MODEL_STDOUT
    ERROR_VARIABLE SIRIUS_OPERATING_MODEL_STDERR)
if(NOT SIRIUS_OPERATING_MODEL_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Sirius operating-model configuration failed:\n"
        "${SIRIUS_OPERATING_MODEL_STDOUT}${SIRIUS_OPERATING_MODEL_STDERR}")
endif()
find_package(Git QUIET)
set(SIRIUS_SOURCE_REVISION "0000000000000000000000000000000000000000")
set(SIRIUS_SOURCE_TREE_CLEAN false)
set(SIRIUS_HAS_GIT_WORKTREE false)
if(Git_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set(SIRIUS_HAS_GIT_WORKTREE true)
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" rev-parse HEAD
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        RESULT_VARIABLE SIRIUS_GIT_REVISION_RESULT
        OUTPUT_VARIABLE SIRIUS_SOURCE_REVISION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" status --porcelain --untracked-files=normal
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        RESULT_VARIABLE SIRIUS_GIT_STATUS_RESULT
        OUTPUT_VARIABLE SIRIUS_GIT_STATUS
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
    if(NOT SIRIUS_GIT_REVISION_RESULT EQUAL 0 OR
       NOT SIRIUS_GIT_STATUS_RESULT EQUAL 0)
        message(FATAL_ERROR "Could not establish the source revision and cleanliness")
    endif()
    if(SIRIUS_GIT_STATUS STREQUAL "")
        set(SIRIUS_SOURCE_TREE_CLEAN true)
    endif()
elseif(NOT SIRIUS_ALIGNMENT_MODE STREQUAL "development")
    message(FATAL_ERROR
        "Qualification/release alignment requires a Git worktree so evidence can bind "
        "the exact clean source revision")
endif()

set(SIRIUS_ALIGNMENT_RECEIPT
    "${CMAKE_BINARY_DIR}/generated/sirius/alignment_receipt.json")
set(SIRIUS_ALIGNMENT_COMMAND
    "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/verify-alignment.py"
    --mode "${SIRIUS_ALIGNMENT_MODE}"
    --source-revision "${SIRIUS_SOURCE_REVISION}"
    --source-tree-clean "${SIRIUS_SOURCE_TREE_CLEAN}"
    --operating-model "${SIRIUS_OPERATING_MODEL}"
    --output "${SIRIUS_ALIGNMENT_RECEIPT}")
if(SIRIUS_HAS_GIT_WORKTREE)
    list(APPEND SIRIUS_ALIGNMENT_COMMAND --source-root "${CMAKE_SOURCE_DIR}")
endif()
if(NOT SIRIUS_ATTESTATION_ROOT STREQUAL "")
    if(EXISTS "${SIRIUS_ATTESTATION_ROOT}")
        file(GLOB_RECURSE SIRIUS_ATTESTATION_INPUTS CONFIGURE_DEPENDS
            "${SIRIUS_ATTESTATION_ROOT}/*.json")
    endif()
    list(APPEND SIRIUS_ALIGNMENT_COMMAND
        --attestation-root "${SIRIUS_ATTESTATION_ROOT}")
endif()
execute_process(
    COMMAND ${SIRIUS_ALIGNMENT_COMMAND}
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    RESULT_VARIABLE SIRIUS_ALIGNMENT_RESULT
    OUTPUT_VARIABLE SIRIUS_ALIGNMENT_STDOUT
    ERROR_VARIABLE SIRIUS_ALIGNMENT_STDERR)
if(NOT SIRIUS_ALIGNMENT_RESULT EQUAL 0)
    message(FATAL_ERROR
        "Sirius alignment configuration failed:\n"
        "${SIRIUS_ALIGNMENT_STDOUT}${SIRIUS_ALIGNMENT_STDERR}")
endif()
string(STRIP "${SIRIUS_ALIGNMENT_STDOUT}" SIRIUS_ALIGNMENT_SUMMARY)
message(STATUS "${SIRIUS_ALIGNMENT_SUMMARY}")

add_custom_target(SiriusAlignmentGate ALL
    COMMAND ${SIRIUS_ALIGNMENT_COMMAND} --check
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    COMMENT "Verifying revision-bound upstream/downstream alignment")
