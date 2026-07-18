# Double-precision Boyer-Lindquist validation oracle (sirius::oracle).
# Header-only: the abstract metric interface, the analytic Kerr metric with
# closed-form Christoffel and Riemann tensors, the symplectic reference
# integrator, the transport value types, and the beam integrator. Off the
# render path; it is the fidelity reference the test estate and every
# downstream gate measure against. Self-contained by design: it depends only
# on sirius_base (contracts), never on sirius_core.

add_library(sirius_oracle INTERFACE)

# Include root is src/, so headers resolve as "sirius/oracle/<name>.h".
target_include_directories(sirius_oracle INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/../..")

target_link_libraries(sirius_oracle INTERFACE sirius_base)
