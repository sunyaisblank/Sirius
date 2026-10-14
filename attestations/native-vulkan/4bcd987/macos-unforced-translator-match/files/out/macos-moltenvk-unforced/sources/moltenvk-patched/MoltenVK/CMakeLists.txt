################################################################################
# MoltenVK::MoltenVK Library
################################################################################

# Find the source files in this directory. We prefer to use GLOB because it is
# easier to maintain in conjunction with the Xcode project file.
file(GLOB_RECURSE SOURCES "MoltenVK/*.c" "MoltenVK/*.cpp" "MoltenVK/*.m" "MoltenVK/*.mm")

add_library(MoltenVK SHARED ${SOURCES})
add_library(MoltenVK::MoltenVK ALIAS MoltenVK)

target_include_directories(MoltenVK
	PUBLIC
		"${CMAKE_CURRENT_SOURCE_DIR}/include"
	PRIVATE
		"${CMAKE_CURRENT_SOURCE_DIR}/MoltenVK/API"
		"${CMAKE_CURRENT_SOURCE_DIR}/MoltenVK/Commands"
		"${CMAKE_CURRENT_SOURCE_DIR}/MoltenVK/GPUObjects"
		"${CMAKE_CURRENT_SOURCE_DIR}/MoltenVK/Layers"
		"${CMAKE_CURRENT_SOURCE_DIR}/MoltenVK/OS"
		"${CMAKE_CURRENT_SOURCE_DIR}/MoltenVK/Utility"
		"${CMAKE_CURRENT_SOURCE_DIR}/MoltenVK/Vulkan"
		"${CMAKE_CURRENT_BINARY_DIR}" # for mvkGitRevDerived.h
)

# Sirius issue-48 diagnostic: this is a source/configuration content discriminator.
# It is NOT an upstream or local Git revision, despite the legacy template name.
if(NOT DEFINED MVK_DIAGNOSTIC_CONTENT_ID)
    message(FATAL_ERROR "The diagnostic provider requires its declared content discriminator")
endif()
string(LENGTH "${MVK_DIAGNOSTIC_CONTENT_ID}" MVK_DIAGNOSTIC_CONTENT_ID_LENGTH)
if(NOT MVK_DIAGNOSTIC_CONTENT_ID_LENGTH EQUAL 64 OR
   NOT MVK_DIAGNOSTIC_CONTENT_ID MATCHES "^[0-9a-f]+$")
    message(FATAL_ERROR "The diagnostic content discriminator must be 64 lowercase hex digits")
endif()
string(SUBSTRING "${MVK_DIAGNOSTIC_CONTENT_ID}" 0 8 MVK_DIAGNOSTIC_CACHE_PREFIX)
if(MVK_DIAGNOSTIC_CACHE_PREFIX STREQUAL "00000000" OR
   MVK_DIAGNOSTIC_CACHE_PREFIX STREQUAL "db445ff2")
    message(FATAL_ERROR "The diagnostic cache prefix is zero or aliases the upstream base")
endif()
set(MVK_GIT_REV "${MVK_DIAGNOSTIC_CONTENT_ID}")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/../Templates/cmake/mvkGitRevDerived.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/mvkGitRevDerived.h"
)
target_compile_definitions(MoltenVK PRIVATE
    "MVK_DIAGNOSTIC_CONTENT_ID=\"${MVK_DIAGNOSTIC_CONTENT_ID}\""
)


target_compile_definitions(MoltenVK PUBLIC
	"MVK_FRAMEWORK_VERSION=${PROJECT_VERSION}"
)

set_target_properties(MoltenVK PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(MoltenVK PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR})

################################################################################
# Dependencies
################################################################################

target_link_libraries(MoltenVK PUBLIC MoltenVK::ShaderConverter)

if (NOT MVK_EXCLUDE_CEREAL)
	include(cereal)
	target_link_libraries(MoltenVK PRIVATE cereal::cereal)
else()
	target_compile_definitions(MoltenVK PUBLIC MVK_USE_CEREAL=0)
endif()

include(Vulkan-Headers)
target_link_libraries(MoltenVK PUBLIC Vulkan::Headers)