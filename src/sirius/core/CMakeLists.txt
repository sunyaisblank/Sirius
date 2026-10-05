# sirius_core: the CPU-side physics library (tensors, autodiff, coordinates,
# metric families, geodesic integrator, disk, spectral). Everything is
# header-only except the two compiled translation units below. The src/ include
# root and the sirius_base link both arrive PUBLIC from sirius_base.
add_library(sirius_core STATIC
    tensor.cpp
    geodesic_integrator.cpp
)

target_link_libraries(sirius_core PUBLIC sirius_base)

sirius_configure_target(sirius_core)

# The x86-64 GNU variation-stage clones may use hardware for explicit FMA,
# but ordinary multiply/add operations must retain their ordered rounding.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_SYSTEM_NAME STREQUAL "Linux"
   AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$")
    set_property(SOURCE geodesic_integrator.cpp APPEND PROPERTY
        COMPILE_OPTIONS -ffp-contract=off)
endif()
