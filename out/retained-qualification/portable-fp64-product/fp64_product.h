#ifndef SIRIUS_ISOLATED_FP64_PRODUCT_H
#define SIRIUS_ISOLATED_FP64_PRODUCT_H

// SOURCE-ONLY PENDING CANDIDATE. Mathematical/provider review and the root's
// execution window must settle before compilation or dispatch. No production
// integration or retained-stage acceptance is asserted by this source.
#include "fp64_scalar.h"
#if defined(__cplusplus)
namespace sirius::portable_fp64_scalar {
#define PB64_PRODUCT_INLINE inline
#else
#define PB64_PRODUCT_INLINE
#endif

PB64_PRODUCT_INLINE PB32Wide PB64ProductWords(double value) {
    PB32Wide words;
#if defined(__cplusplus)
    std::uint64_t bits = std::bit_cast<std::uint64_t>(value);
    words.lo = uint(bits);
    words.hi = uint(bits >> 32);
#else
    asuint(value, words.lo, words.hi);
#endif
    return words;
}

// Same bounded high/low/valid value contract as canonical PB32ProductResidual.
// Proof obligations pending: correct 53-bit native normal binary64 product,
// exact Sterbenz subtraction, manual projection and all signed-zero joins.
PB64_PRODUCT_INLINE PB32Product PB64ProductResidual(uint a, uint b) {
    PB32Product result;
    result.high = 0u;
    result.low = 0u;
    result.valid = 0u;
    if (PB32Abs(a) > 0x5d800000u || PB32Abs(b) > 0x5d800000u) return result;
    result.valid = 1u;
    if (PB32Zero(a) || PB32Zero(b)) {
        result.high = (a ^ b) & 0x80000000u;
        return result;  // Signed exact high zero; exact residual is +0.
    }
    double left = PB64Decode(a);
    double right = PB64Decode(b);
#if defined(__cplusplus)
    double product = left * right;
#else
    double product = spirv_asm {
        OpFMul $$double %v $left $right;
        OpDecorate %v NoContraction;
        OpCopyObject $$double result %v;
    };
#endif
    result.high = PB64Project(product);
    if (PB32Zero(result.high)) {
        result.low = result.high;  // Nonzero exact product underflows with its sign.
        return result;
    }
    double rounded_high = PB64Decode(result.high);
#if defined(__cplusplus)
    double residual = product - rounded_high;
#else
    double residual = spirv_asm {
        OpFSub $$double %v $product $rounded_high;
        OpDecorate %v NoContraction;
        OpCopyObject $$double result %v;
    };
#endif
    PB32Wide residual_words = PB64ProductWords(residual);
    if ((residual_words.lo | (residual_words.hi & 0x7fffffffu)) == 0u)
        return result;  // Force exact-zero residual +0 independently of native sign.
    result.low = PB64Project(residual);  // Nonzero underflow preserves residual sign.
    return result;
}

PB64_PRODUCT_INLINE PB32Product PB64ProductEvaluate(uint operation, uint a, uint b) {
    // Keep the canonical residual outside the reachable probe call graph.
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
        case 13u: result = PB64ProductResidual(a, b); break;
        default: result.valid = 0u; break;
    }
    return result;
}

#undef PB64_PRODUCT_INLINE
#if defined(__cplusplus)
}  // namespace sirius::portable_fp64_scalar
#endif
#endif
