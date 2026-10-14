# First-party layers. sirius_configure_target owns the C++ standard policy.
unset(CMAKE_CXX_STANDARD)
unset(CMAKE_CXX_STANDARD_REQUIRED)

add_library(sirius_base STATIC
    base/contracts.cpp
    base/resource_locator.cpp
)

target_include_directories(sirius_base PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/..")

# This is derived directly from the non-negotiable alignment mode rather than a
# separate cache option. Keeping it on the base target makes resource lookup
# release-locked in every downstream consumer, including render and viewer code.
if(NOT SIRIUS_ALIGNMENT_MODE STREQUAL "development")
    target_compile_definitions(sirius_base PUBLIC SIRIUS_RELEASE_RESOURCE_LOCKED=1)
else()
    target_compile_definitions(sirius_base PUBLIC SIRIUS_RELEASE_RESOURCE_LOCKED=0)
endif()

sirius_configure_target(sirius_base)

# Every layer is required and owns its sources in its own CMakeLists.txt.
foreach(module core oracle backend render app)
    add_subdirectory(${module})
endforeach()
