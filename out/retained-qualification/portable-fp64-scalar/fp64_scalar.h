#ifndef SIRIUS_ISOLATED_FP64_SCALAR_H
#define SIRIUS_ISOLATED_FP64_SCALAR_H

// Isolated prototype; not a production arithmetic implementation.
#include "sirius/kernels/portable_binary32.h"
#if defined(__cplusplus)
#include <bit>
#include <limits>
namespace sirius::portable_fp64_scalar {
using namespace sirius::portable_binary32;
static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559 &&
              std::numeric_limits<double>::digits == 53);
#define PB64_INLINE inline
#define PB64_PRECISE
#else
#define PB64_INLINE
#define PB64_PRECISE precise
#endif

// Called only for finite nonzero binary32. Decode words directly; no float32
// conversion, arithmetic or denormal reaches the native floating unit.
PB64_INLINE double PB64Decode(uint a) {
    PB32Finite decoded = PB32Unpack(a);
    uint fraction = decoded.significand & 0x007fffffu;
    uint lo = fraction << 29;
    uint hi = (a & 0x80000000u) | (uint(decoded.exponent + 1023) << 20) |
              (fraction >> 3);
#if defined(__cplusplus)
    return std::bit_cast<double>((std::uint64_t(hi) << 32) | lo);
#else
    return asdouble(lo, hi);
#endif
}

// Finite nonzero sums/products of binary32 operands are normal binary64:
// minimum nonzero sum 2^-149, product 2^-298; maximum product below 2^256.
// Recover 53 significand bits as two uint32 words and jam to PB32Round's
// leading-bit-26/three-rounding-bit input. Only integer binary32 rounding.
PB64_INLINE uint PB64Project(double value) {
    uint lo, hi;
#if defined(__cplusplus)
    std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
    lo = uint(bits);
    hi = uint(bits >> 32);
#else
    asuint(value, lo, hi);
#endif
    uint sign = hi >> 31;
    int exponent = int((hi >> 20) & 0x7ffu) - 1023;
    PB32Wide significand;
    significand.lo = lo;
    significand.hi = (hi & 0x000fffffu) | 0x00100000u;
    return PB32Round(sign, exponent, PB32WideShiftJam(significand, 26u));
}

PB64_INLINE uint PB64Add(uint a, uint b) {
    if (PB32Nan(a) || PB32Nan(b)) return 0x7fc00000u;
    if (PB32Infinity(a) || PB32Infinity(b)) {
        if (PB32Infinity(a) && PB32Infinity(b) && ((a ^ b) >> 31) != 0u)
            return 0x7fc00000u;
        return PB32Infinity(a) ? a : b;
    }
    if (PB32Zero(a) && PB32Zero(b)) return (a & b) & 0x80000000u;
    if (PB32Zero(a)) return b;
    if (PB32Zero(b)) return a;
    if (PB32Abs(a) == PB32Abs(b) && ((a ^ b) >> 31) != 0u) return 0u;
    PB64_PRECISE double left = PB64Decode(a);
    PB64_PRECISE double right = PB64Decode(b);
#if defined(__cplusplus)
    double sum = left + right;
#else
    double sum = spirv_asm {
        OpFAdd $$double %v $left $right;
        OpDecorate %v NoContraction;
        OpCopyObject $$double result %v;
    };
#endif
    return PB64Project(sum);
}
PB64_INLINE uint PB64Subtract(uint a, uint b) { return PB64Add(a, PB32Negate(b)); }
PB64_INLINE uint PB64Multiply(uint a, uint b) {
    uint sign = (a ^ b) >> 31;
    if (PB32Nan(a) || PB32Nan(b)) return 0x7fc00000u;
    if (PB32Infinity(a) || PB32Infinity(b))
        return PB32Zero(a) || PB32Zero(b) ? 0x7fc00000u : (sign << 31) | 0x7f800000u;
    if (PB32Zero(a) || PB32Zero(b)) return sign << 31;
    PB64_PRECISE double left = PB64Decode(a);
    PB64_PRECISE double right = PB64Decode(b);
#if defined(__cplusplus)
    double product = left * right;
#else
    double product = spirv_asm {
        OpFMul $$double %v $left $right;
        OpDecorate %v NoContraction;
        OpCopyObject $$double result %v;
    };
#endif
    return PB64Project(product);
}

PB64_INLINE PB32Product PB64Evaluate(uint operation, uint a, uint b) {
    PB32Product result;
    result.high = 0u;
    result.low = 0u;
    result.valid = 1u;
    switch (operation) {
        case 0u: result.high = PB64Add(a, b); break;
        case 1u: result.high = PB64Subtract(a, b); break;
        case 2u: result.high = PB64Multiply(a, b); break;
        case 3u: result.high = PB32Divide(a, b); break;
        case 4u: result.high = PB32Sqrt(a); break;
        case 5u: result.high = PB32Equal(a, b) ? 1u : 0u; break;
        case 6u: result.high = PB32Less(a, b) ? 1u : 0u; break;
        case 7u: result.high = PB32LessEqual(a, b) ? 1u : 0u; break;
        case 8u: result.high = PB32FromUnsigned(a); break;
        case 9u: result.high = PB32FromSignedBits(a); break;
        case 10u: result.high = PB32Abs(a); break;
        case 11u: result.high = PB32Negate(a); break;
        case 12u: result.high = PB32Classify(a); break;
        case 13u: result = PB32ProductResidual(a, b); break;
        default: result.valid = 0u; break;
    }
    return result;
}

#undef PB64_INLINE
#undef PB64_PRECISE
#if defined(__cplusplus)
}  // namespace sirius::portable_fp64_scalar
#endif
#endif
