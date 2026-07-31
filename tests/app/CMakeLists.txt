# App suite: CLI argument parsing, configuration loading/validation, JSON schema
# compatibility, and platform paths. Links sirius_app_lib so the CLI/config
# objects are exercised directly (July ledger item 8: the CLI parsing tests were
# never linked to the CLI objects; they are linked here).

add_executable(sirius_app_tests
    config_validation_test.cpp
    config_environment_test.cpp
    config_schema_roundtrip_test.cpp
    render_command_parse_test.cpp
    platform_paths_test.cpp
    config_loading_test.cpp
    command_router_test.cpp
)

target_link_libraries(sirius_app_tests PRIVATE sirius_app_lib GTest::gtest_main)
target_include_directories(sirius_app_tests PRIVATE "${CMAKE_SOURCE_DIR}/tests")

sirius_configure_target(sirius_app_tests)

if(TARGET sirius)
    # RenderCommand is exercised through the application library, so ensure its
    # runtime-volume test sees the same staged resources as the real executable.
    add_dependencies(sirius_app_tests sirius)
endif()

gtest_discover_tests(
    sirius_app_tests
    PROPERTIES
        LABELS "Correctness"
        ENVIRONMENT "SIRIUS_RESOURCE_DIR=$<TARGET_FILE_DIR:sirius>/resources")
