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
)

target_link_libraries(sirius_app_tests PRIVATE sirius_app_lib GTest::gtest_main)
target_include_directories(sirius_app_tests PRIVATE "${CMAKE_SOURCE_DIR}/tests")

sirius_configure_target(sirius_app_tests)

gtest_discover_tests(sirius_app_tests PROPERTIES LABELS "Correctness")
