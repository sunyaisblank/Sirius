# === Install and packaging ===
if(TARGET sirius)
    install(TARGETS sirius RUNTIME DESTINATION bin)
    install(FILES assets/Starfield.png DESTINATION share/sirius/assets)
    install(FILES tests/operating_model.json DESTINATION share/sirius/model)
    install(FILES "${SIRIUS_ALIGNMENT_RECEIPT}"
        DESTINATION share/sirius/model)
    if(NOT SIRIUS_ALIGNMENT_MODE STREQUAL "development")
        install(FILES "${SIRIUS_MANDATORY_GATE_STAMP}"
            DESTINATION share/sirius/model)
    endif()
    if(SIRIUS_INSTALL_HAS_VIEWER)
        install(FILES
            src/sirius/app/viewer/shaders/RDSD003A.vert
            src/sirius/app/viewer/shaders/RDSD003A.frag
            DESTINATION share/sirius/shaders)
    endif()
    if(TARGET sirius_kernels)
        install(FILES
            "${SIRIUS_KERNEL_BINARY_DIR}/trace.spv"
            "${SIRIUS_KERNEL_BINARY_DIR}/trace_fp32comp.spv"
            "${SIRIUS_KERNEL_BINARY_DIR}/trace_fp64.spv"
            DESTINATION share/sirius/kernels)
        set(SIRIUS_INSTALL_HAS_VULKAN TRUE)
    else()
        set(SIRIUS_INSTALL_HAS_VULKAN FALSE)
    endif()
    configure_file(
        cmake/verify_sirius_install.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/verify_sirius_install.cmake"
        @ONLY)
    install(SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/verify_sirius_install.cmake")
endif()

set(CPACK_PACKAGE_NAME "sirius")
set(CPACK_PACKAGE_VENDOR "Sirius")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "General-relativistic ray tracing engine")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
if(SIRIUS_ALIGNMENT_MODE STREQUAL "release")
    include(CPack)
else()
    message(STATUS
        "Release packaging: DISABLED (configure SIRIUS_ALIGNMENT_MODE=release "
        "with the complete attestation set)")
endif()
