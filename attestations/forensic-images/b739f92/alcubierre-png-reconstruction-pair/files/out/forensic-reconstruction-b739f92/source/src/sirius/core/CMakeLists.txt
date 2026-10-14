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
