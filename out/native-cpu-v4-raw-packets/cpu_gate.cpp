// A baseline-ISA fixture precondition; no numerical state or packet is changed.
// Production selection belongs to rows_adapter.cpp and its sealed dispatcher.
#include <cstdlib>
extern "C" bool sirius_native_v4_supported() noexcept;
__attribute__((constructor)) static void require_actual_v4() {
    __builtin_cpu_init();
    if (!sirius_native_v4_supported()) std::abort();
}
