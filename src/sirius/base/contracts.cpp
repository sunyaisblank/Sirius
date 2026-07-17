#include "sirius/base/contracts.h"

#include <cstdio>
#include <cstdlib>

namespace sirius::base {

namespace {

void PrintViolation(const char* action, const char* kind, const char* expression,
                    const char* file, long line) noexcept {
    std::fprintf(stderr, "sirius: %s violated (%s): %s at %s:%ld\n", kind, action, expression,
                 file, line);
    std::fflush(stderr);
}

}  // namespace

void ReportContractViolation(const char* kind, const char* expression, const char* file,
                             long line) noexcept {
    PrintViolation("observed, continuing", kind, expression, file, line);
}

void EnforceContractViolation(const char* kind, const char* expression, const char* file,
                              long line) noexcept {
    PrintViolation("enforced, terminating", kind, expression, file, line);
    std::abort();
}

}  // namespace sirius::base
