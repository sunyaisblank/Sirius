# App suite: non-rendering CLI argument parsing, configuration
# loading/validation, JSON schema compatibility, session projection, display
# transfer contracts, input state, and platform paths. Render-session, dispatch, and frame-publication
# cases are owned by sirius_render_tests; shared source governance rejects their
# return.
# Links sirius_app_lib so the CLI/config objects are exercised directly (July
# ledger item 8: the CLI parsing tests were never linked to the CLI objects).

add_executable(sirius_app_tests
    alignment_authority_test.cpp
    config_validation_test.cpp
    config_environment_test.cpp
    config_schema_roundtrip_test.cpp
    render_command_parse_test.cpp
    session_config_adapter_test.cpp
    platform_paths_test.cpp
    config_loading_test.cpp
    command_router_test.cpp
    viewer_display_contract_test.cpp
)

target_link_libraries(sirius_app_tests PRIVATE sirius_app_lib GTest::gtest_main)
target_include_directories(sirius_app_tests PRIVATE "${CMAKE_SOURCE_DIR}/tests")

sirius_configure_target(sirius_app_tests)

if(TARGET sirius)
    # Resource-locked qualification/release tests must consume the exact staged
    # product volume. Place this executable beside sirius instead of relying on
    # a development-only resource override.
    add_dependencies(sirius_app_tests sirius)
    set_target_properties(sirius_app_tests PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "$<TARGET_FILE_DIR:sirius>")
endif()

gtest_discover_tests(
    sirius_app_tests
    PROPERTIES
        LABELS "Correctness")
