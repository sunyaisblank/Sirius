set(SIRIUS_SOFTWARE_VULKAN_BUNDLE "" CACHE PATH
    "Directory containing the qualified Linux software Vulkan driver and retained modules")
if(NOT SIRIUS_SOFTWARE_VULKAN_BUNDLE)
    return()
endif()

if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux" OR NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64)$")
    message(FATAL_ERROR "The packaged software Vulkan bundle requires Linux x86-64")
endif()

file(REAL_PATH "${SIRIUS_SOFTWARE_VULKAN_BUNDLE}" _sirius_software_bundle)
set(SIRIUS_PACKAGED_SOFTWARE_DRIVER "${_sirius_software_bundle}/libvulkan_lvp.so")
set(SIRIUS_PACKAGED_SOFTWARE_DRIVER_SHA256
    "797cba316a7cb3eb438b0ce3810b9351ee54f0406450f7b2d06d8e1e33e85563")
set(SIRIUS_PACKAGED_SOFTWARE_LLVM_SHA256
    "fd419f2d3eddfbdece987a7757da8f66794ab933dde48fa685c2ad1c9478225e")
set(SIRIUS_PACKAGED_SOFTWARE_SUBDIR "providers/linux-x86_64")

if(NOT EXISTS "${SIRIUS_PACKAGED_SOFTWARE_DRIVER}")
    message(FATAL_ERROR "The qualified software Vulkan driver is missing")
endif()
file(REAL_PATH "${SIRIUS_PACKAGED_SOFTWARE_DRIVER}" SIRIUS_PACKAGED_SOFTWARE_DRIVER)
file(SHA256 "${SIRIUS_PACKAGED_SOFTWARE_DRIVER}" _sirius_software_driver_digest)
if(NOT _sirius_software_driver_digest STREQUAL SIRIUS_PACKAGED_SOFTWARE_DRIVER_SHA256)
    message(FATAL_ERROR "The software Vulkan driver differs from the qualified SHA-256")
endif()

function(sirius_embed_packaged_software retained_directory retained_header)
    set(module_file "${retained_directory}/sealed_shaders.h")
    set(ordinary_sources "")
    set(sealed_sources "")
    foreach(stage transport endpoint dopri_phase)
        list(APPEND ordinary_sources "${retained_directory}/retained_${stage}.spv")
        list(APPEND sealed_sources "${_sirius_software_bundle}/${stage}.spv")
    endforeach()
    add_custom_command(
        OUTPUT "${module_file}"
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_SOURCE_DIR}/scripts/embed-packaged-software.py"
            --bundle "${_sirius_software_bundle}"
            --ordinary "${retained_directory}"
            --output "${module_file}"
            --driver-sha256 "${SIRIUS_PACKAGED_SOFTWARE_DRIVER_SHA256}"
            --llvm-sha256 "${SIRIUS_PACKAGED_SOFTWARE_LLVM_SHA256}"
        DEPENDS "${retained_header}" ${ordinary_sources} ${sealed_sources}
            "${CMAKE_SOURCE_DIR}/scripts/embed-packaged-software.py"
        COMMENT "Verifying and embedding the qualified software Vulkan modules"
        VERBATIM)
endfunction()


function(sirius_packaged_software_stage_arguments target result)
    set(directory "$<TARGET_FILE_DIR:${target}>/resources/${SIRIUS_PACKAGED_SOFTWARE_SUBDIR}")
    set(${result}
        COMMAND "${CMAKE_COMMAND}" -E make_directory "${directory}"
        COMMAND "${CMAKE_COMMAND}" -E rm -f "${directory}/libvulkan_lvp.so"
        COMMAND "${CMAKE_COMMAND}" -E create_hardlink
            "${SIRIUS_PACKAGED_SOFTWARE_DRIVER}" "${directory}/libvulkan_lvp.so"
        PARENT_SCOPE)
endfunction()
