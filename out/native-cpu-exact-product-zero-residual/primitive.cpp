// Two unchanged source prefixes supplied by primitive.py. Only the candidate
// product32 prefix differs. No candidate predicate is used by the fixture.
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <immintrin.h>
#include <limits>
#include <initializer_list>
#include "original-prefix.hpp"
#include "candidate-prefix.hpp"

#define WRAPPER(label) \
extern "C" void label##_product(const std::uint32_t* a, const std::uint32_t* b, \
                              std::uint32_t* high, std::uint32_t* low) noexcept { \
    label::F aa, bb, hi, lo; \
    std::memcpy(&aa, a, sizeof(aa)); std::memcpy(&bb, b, sizeof(bb)); \
    label::product32(aa, bb, hi, lo); \
    std::memcpy(high, &hi, sizeof(hi)); std::memcpy(low, &lo, sizeof(lo)); \
}
WRAPPER(original)
WRAPPER(candidate)
#undef WRAPPER
