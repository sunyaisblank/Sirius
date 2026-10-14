#pragma once

// Contract macros per docs/STYLE.md section 2. P2900 contracts are the
// intended end state; on the measured toolchain (GCC 14, Clang 21, MSVC
// 19.4x) no release compiler ships them, so the macros provide checked
// semantics today and a mechanical migration path later. The native mapping
// below is contract_assert for all three macros because these are
// body-placed statements, not declaration-placed pre/post; the migration to
// declaration syntax is the scripted rewrite STYLE.md describes.
//
// SIRIUS_CONTRACT_MODE: 2 enforce (violation terminates), 1 observe
// (violation reported once per site, execution continues), 0 ignore
// (expression compiled but never evaluated). The mode is a build-system
// decision (cmake/sirius_toolchain.cmake); sources never set it.

#ifndef SIRIUS_CONTRACT_MODE
#if defined(NDEBUG)
#define SIRIUS_CONTRACT_MODE 1
#else
#define SIRIUS_CONTRACT_MODE 2
#endif
#endif

#if defined(__cpp_contracts) && __cpp_contracts >= 202502L && SIRIUS_CONTRACT_MODE == 2

#define SIRIUS_PRE(expr) contract_assert(expr)
#define SIRIUS_POST(expr) contract_assert(expr)
#define SIRIUS_ASSERT(expr) contract_assert(expr)

#else

#include <atomic>

namespace sirius::base {

// Reports and returns; the observe-mode channel. Reporting is stderr rather
// than a logging seam so that contracts stay usable below every other layer.
void ReportContractViolation(const char* kind, const char* expression, const char* file,
                             long line) noexcept;

// Reports and terminates; the enforce-mode channel.
[[noreturn]] void EnforceContractViolation(const char* kind, const char* expression,
                                           const char* file, long line) noexcept;

}  // namespace sirius::base

#if SIRIUS_CONTRACT_MODE == 2

#define SIRIUS_CONTRACT_CHECK_(kind, expr)                                              \
    do {                                                                                \
        if (!(expr)) {                                                                  \
            ::sirius::base::EnforceContractViolation(kind, #expr, __FILE__, __LINE__);  \
        }                                                                               \
    } while (false)

#elif SIRIUS_CONTRACT_MODE == 1

#define SIRIUS_CONTRACT_CHECK_(kind, expr)                                                \
    do {                                                                                  \
        if (!(expr)) {                                                                    \
            static ::std::atomic<bool> sirius_contract_reported_{false};                  \
            if (!sirius_contract_reported_.exchange(true, ::std::memory_order_relaxed)) { \
                ::sirius::base::ReportContractViolation(kind, #expr, __FILE__, __LINE__); \
            }                                                                             \
        }                                                                                 \
    } while (false)

#else

// The expression still has to compile, so a renamed symbol cannot silently
// strand a contract in ignore-mode builds.
#define SIRIUS_CONTRACT_CHECK_(kind, expr)   \
    do {                                     \
        if (false) {                         \
            static_cast<void>(expr);         \
        }                                    \
    } while (false)

#endif

#define SIRIUS_PRE(expr) SIRIUS_CONTRACT_CHECK_("precondition", expr)
#define SIRIUS_POST(expr) SIRIUS_CONTRACT_CHECK_("postcondition", expr)
#define SIRIUS_ASSERT(expr) SIRIUS_CONTRACT_CHECK_("assertion", expr)

#endif

// Never evaluated in any mode: documentation of a condition too expensive to
// check, kept compiling so it cannot rot.
#define SIRIUS_AXIOM(expr)           \
    do {                             \
        if (false) {                 \
            static_cast<void>(expr); \
        }                            \
    } while (false)
