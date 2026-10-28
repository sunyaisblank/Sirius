// Diagnostic routing only. Original fixture objects/assertions are not rebuilt.
#include "admitted_bridge.h"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <immintrin.h>

#define ARGS const void* input, std::uint32_t input_bytes, void* output, std::uint32_t output_bytes, \
    std::uint32_t first_x, std::uint32_t y, std::uint32_t z, std::uint32_t grid_x, \
    std::uint32_t grid_y, std::uint32_t grid_z, std::uint32_t groups
#define VALUES input, input_bytes, output, output_bytes, first_x, y, z, grid_x, grid_y, grid_z, groups
using Packet = bool (*)(ARGS) noexcept;

namespace {
struct Counts {
    std::uint64_t public_scalar = 0, public_rows = 0;
    std::uint64_t preflight_pass = 0, preflight_fail = 0;
    std::uint64_t private_scalar = 0, private_rows = 0, private_refusals = 0;
    std::uint64_t small_rows = 0;
} counts;

[[noreturn]] void fail(const char* why) noexcept {
    std::fprintf(stderr, "diagnostic route failed: %s\n", why);
    std::abort();
}

bool admitted_mode() noexcept {
    const char* mode = std::getenv("SIRIUS_TOKEN_FIXTURE_MODE");
    if (!mode) fail("missing fixed route mode");
    if (std::strcmp(mode, "public") == 0) return false;
    if (std::strcmp(mode, "admitted") == 0) return true;
    fail("unknown route mode");
}

bool route(Packet pure, Packet guarded, Packet trusted, bool scalar, ARGS) noexcept {
    if (!admitted_mode()) {
        if (scalar) ++counts.public_scalar; else ++counts.public_rows;
        return guarded(VALUES);
    }
    // Actual unchanged pure preflight, using Mesa's x=y=z=0/groups=1 tuple.
    // Views/extents/grid are exactly the call's. Immutable input is retained
    // through this synchronous return; no capability is cached between calls.
    const unsigned before = _mm_getcsr();
    const bool admitted = pure(input, input_bytes, output, output_bytes,
                              0, 0, 0, grid_x, grid_y, grid_z, 1);
    if (_mm_getcsr() != before) fail("pure preflight changed MXCSR");
    if (!admitted) {
        ++counts.preflight_fail;
        if (scalar) ++counts.public_scalar; else ++counts.public_rows;
        const bool handled = guarded(VALUES);
        if (handled) fail("failed dispatch preflight unexpectedly handled");
        return false;  // No trusted call is reachable on this branch.
    }
    ++counts.preflight_pass;
    if (scalar) ++counts.private_scalar; else ++counts.private_rows;
    if (!scalar && groups < 8) ++counts.small_rows;
    const bool handled = trusted(VALUES);
    if (!handled) ++counts.private_refusals;
    return handled;
}

__attribute__((destructor)) void report_routes() noexcept {
    const char* destination = std::getenv("SIRIUS_TOKEN_FIXTURE_ROUTES");
    if (!destination) fail("missing diagnostic route output");
    FILE* file = std::fopen(destination, "wx");
    if (!file) fail("cannot create diagnostic route output");
    const int n = std::fprintf(file,
        "{\"mode\":\"%s\",\"public_scalar\":%llu,\"public_rows\":%llu,"
        "\"preflight_pass\":%llu,\"preflight_fail\":%llu,"
        "\"private_scalar\":%llu,\"private_rows\":%llu,\"private_refusals\":%llu,"
        "\"admitted_small_row_packets\":%llu,\"trusted_call_after_failed_preflight\":false,"
        "\"scope\":\"Diagnostic pure-preflight routing; actual Mesa job capability branch unobserved\"}\n",
        admitted_mode() ? "admitted" : "public",
        static_cast<unsigned long long>(counts.public_scalar),
        static_cast<unsigned long long>(counts.public_rows),
        static_cast<unsigned long long>(counts.preflight_pass),
        static_cast<unsigned long long>(counts.preflight_fail),
        static_cast<unsigned long long>(counts.private_scalar),
        static_cast<unsigned long long>(counts.private_rows),
        static_cast<unsigned long long>(counts.private_refusals),
        static_cast<unsigned long long>(counts.small_rows));
    const int closed = std::fclose(file);
    if (n < 0 || closed) fail("cannot finish diagnostic route output");
}
}

// The original Transport fixture uses this identical real CPU/OS query.
// Its own definition and internal uses are renamed together, not replaced.
extern "C" bool sirius_native_simd_supported() noexcept {
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("x86-64-v3");
}

extern "C" bool sirius_token_endpoint_scalar(ARGS) noexcept {
    return route(sirius_native_endpoint_admit_packet, sirius_native_endpoint_try_packet,
                 sirius_native_endpoint_admitted_packet, true, VALUES);
}

extern "C" void sirius_token_transport_raw(std::uint32_t* input, std::uint32_t input_words,
    std::uint32_t* output, std::uint32_t output_words, std::uint32_t first_x,
    std::uint32_t groups) noexcept {
    if (input_words > UINT32_MAX / 4 || output_words > UINT32_MAX / 4) fail("byte count overflow");
    // Exact original Transport fixture active prefix, not allocation capacity40.
    if (!route(sirius_native_transport_admit_packet, sirius_native_transport_rows_try_packet,
               sirius_native_transport_rows_admitted_packet, false, input, input_words * 4,
               output, output_words * 4, first_x, 0, 0, 38, 1, 1, groups))
        fail("original admitted Transport packet refused");
}

extern "C" void sirius_token_endpoint_raw(std::uint32_t* input, std::uint32_t input_words,
    std::uint32_t* output, std::uint32_t output_words, std::uint32_t first_x,
    std::uint32_t groups) noexcept {
    if (input_words > UINT32_MAX / 4 || output_words > UINT32_MAX / 4) fail("byte count overflow");
    // Exact original Endpoint fixture active prefix, not allocation capacity32.
    if (!route(sirius_native_endpoint_admit_packet, sirius_native_endpoint_rows_try_packet,
               sirius_native_endpoint_rows_admitted_packet, false, input, input_words * 4,
               output, output_words * 4, first_x, 0, 0, 23, 1, 1, groups))
        fail("original admitted Endpoint packet refused");
}
#undef VALUES
#undef ARGS
