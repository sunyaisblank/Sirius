// Private candidate. A native estimate has no correctness authority: exact
// two-word inequalities admit it, or the unchanged integer algorithm runs.
struct PB32Estimate {
    uint quotient;
    uint sticky;
    uint valid;
};

SIRIUS_PB32_INLINE PB32Wide PB32WideAddWord(PB32Wide a, uint b) {
    PB32Wide result;
    result.lo = a.lo + b;
    result.hi = a.hi + (result.lo < a.lo ? 1u : 0u);
    return result;
}

SIRIUS_PB32_INLINE uint PB32SeedUnsigned(double seed) {
    // Sanitize before FToU, including inactive SIMD lanes. No NaN, negative,
    // infinite or out-of-range estimate can reach that conversion.
#if defined(__cplusplus)
    uint high = uint(std::bit_cast<std::uint64_t>(seed) >> 32);
#else
    uint low, high;
    asuint(seed, low, high);
#endif
    double safe = high < 0x41b00000u ? seed : double(0);
    return uint(safe);
}

SIRIUS_PB32_INLINE PB32Estimate PB32CheckedDivide(uint numerator, uint denominator, double seed) {
    // D <= R < 2D, 2^23 <= D < 2^24; N=R*2^26 < 2^51.
    // q < 2^28 initially, so all products/corrections fit two words.
    PB32Wide n;
    n.lo = numerator << 26;
    n.hi = numerator >> 6;
    uint q = PB32SeedUnsigned(seed);
    PB32Wide product = PB32WideProduct(q, denominator);
    PB32Wide divisor;
    divisor.lo = denominator;
    divisor.hi = 0u;
    if (PB32WideLess(n, product)) {
        --q;
        product = PB32WideSubtract(product, divisor);
    } else {
        PB32Wide residual = PB32WideSubtract(n, product);
        if (!PB32WideLess(residual, divisor)) {
            ++q;
            product = PB32WideAddWord(product, denominator);
        }
    }
    PB32Wide residual = PB32WideSubtract(n, product);
    PB32Estimate result;
    result.quotient = q;
    result.sticky = (residual.lo | residual.hi) != 0u ? 1u : 0u;
    result.valid = !PB32WideLess(n, product) && PB32WideLess(residual, divisor) ? 1u : 0u;
    // The postcondition is exactly qD <= N < (q+1)D.
    return result;
}

SIRIUS_PB32_INLINE PB32Estimate PB32CheckedSqrt(PB32Wide n, double seed) {
    // 2^52 <= N < 2^54. One exact product and bounded corrections establish
    // q^2 <= N < (q+1)^2 without a native precision/rounding assumption.
    uint q = PB32SeedUnsigned(seed);
    PB32Wide product = PB32WideProduct(q, q);
    if (PB32WideLess(n, product)) {
        PB32Wide decrement;
        decrement.lo = 2u * q - 1u;
        decrement.hi = 0u;
        product = PB32WideSubtract(product, decrement);
        --q;
    } else {
        PB32Wide next = PB32WideAddWord(product, 2u * q + 1u);
        if (!PB32WideLess(n, next)) {
            product = next;
            ++q;
        }
    }
    PB32Wide residual = PB32WideSubtract(n, product);
    PB32Estimate result;
    result.quotient = q;
    result.sticky = (residual.lo | residual.hi) != 0u ? 1u : 0u;
    result.valid = !PB32WideLess(n, product) && residual.hi == 0u && residual.lo < 2u * q + 1u ? 1u : 0u;
    return result;
}
