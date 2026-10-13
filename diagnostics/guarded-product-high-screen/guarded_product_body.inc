// Private compiler screen only. Positive normal hardware products may return
// either enclosing neighbor. Compare on the exact product's integer lattice,
// correct to nearest-even, and retain the original exact residual rounding.
#if defined(__cplusplus)
extern uint PB32HighModel;
extern uint PB32GuardedCalls;
extern uint PB32Corrections;
extern uint PB32DownwardBinades;
#endif
SIRIUS_PB32_INLINE PB32Product PB32ProductResidual(uint a, uint b) {
    uint am = PB32Abs(a), bm = PB32Abs(b);
    if (am < 0x00800000u || bm < 0x00800000u ||
        am > 0x5d800000u || bm > 0x5d800000u)
        return PB32ProductResidualFallback(a, b);
    int sumExponent = int(am >> 23) + int(bm >> 23) - 254;
    if (sumExponent < -126) return PB32ProductResidualFallback(a, b);
    PB32Wide product = PB32WideProduct((am & 0x007fffffu) | 0x00800000u,
                                       (bm & 0x007fffffu) | 0x00800000u);
#if defined(__cplusplus)
    ++PB32GuardedCalls;
    uint high = PB32HighModel;
#else
    float af = asfloat(am), bf = asfloat(bm);
    float estimate = spirv_asm {
        OpFMul $$float %v $af $bf;
        OpDecorate %v NoContraction;
        OpCopyObject $$float result %v;
    };
    uint high = asuint(estimate);
#endif
    // Normal input/product bounds imply shift in {23,24,25}; the enclosing
    // neighbors lie on this integer lattice. No float subtraction is used.
    uint shift = uint(int(high >> 23) - 127 - sumExponent + 23);
    uint significand = (high & 0x007fffffu) | 0x00800000u;
    PB32Wide rounded;
    rounded.lo = significand << shift;
    rounded.hi = significand >> (32u - shift);
    bool roundedUp = PB32WideLess(product, rounded);
    PB32Wide difference = roundedUp ? PB32WideSubtract(rounded, product)
                                   : PB32WideSubtract(product, rounded);
    // At a power of two, the lower neighbor is half the upper spacing.
    bool downwardBinade = roundedUp && (high & 0x007fffffu) == 0u;
#if defined(__cplusplus)
    if (downwardBinade) ++PB32DownwardBinades;
#endif
    uint half = 1u << (shift - 1u - (downwardBinade ? 1u : 0u));
    if (difference.lo > half || (difference.lo == half && (high & 1u) != 0u)) {
        high = roundedUp ? high - 1u : high + 1u;
        difference.lo = (half << 1u) - difference.lo;
        roundedUp = !roundedUp;
#if defined(__cplusplus)
        ++PB32Corrections;
#endif
    }
    uint sign = (a ^ b) >> 31;
    PB32Product result;
    result.high = high | (sign << 31);
    result.low = (difference.lo | difference.hi) == 0u ? 0u :
        PB32RoundWide(sign ^ (roundedUp ? 1u : 0u), sumExponent - 46, difference);
    result.valid = 1u;
    return result;
}
