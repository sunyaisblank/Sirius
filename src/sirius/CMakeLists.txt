# New C++26 tree. The legacy include sets CMAKE_CXX_STANDARD 17 at root
# scope; clear it here so the only standard flag on new targets is the one
# sirius_configure_target adds.
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

# Modules land one at a time during the port; each owns its CMakeLists so
# parallel port workstreams never share a file.
foreach(module core oracle backend render app)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${module}/CMakeLists.txt")
        add_subdirectory(${module})
    endif()
endforeach()
