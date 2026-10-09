#ifndef SIRIUS_KERNELS_PORTABLE_BINARY32_H
#define SIRIUS_KERNELS_PORTABLE_BINARY32_H

// Raw-word IEEE binary32 roundTiesToEven with gradual underflow and signed zero.
// All arithmetic uses uint32, including explicit two-word exact products.
// Arithmetic NaNs are canonical 0x7fc00000; exception flags/payload propagation
// are outside this value API. Bit operations preserve the input payload bits.
// The same source is used by C++ and Slang; no native floating operation is used.
#if defined(__cplusplus)
#include <bit>
#include <cstdint>
namespace sirius::portable_binary32 {
using uint = std::uint32_t;
static_assert(sizeof(uint) == 4 && sizeof(int) == 4);
#define SIRIUS_PB32_INLINE inline
#else
#define SIRIUS_PB32_INLINE
#endif

struct PB32Wide {
    uint lo;
    uint hi;
};
struct PB32Finite {
    uint significand;
    int exponent;
};

// Callers supply a nonzero word. Its leading bit replaces exact one-bit
// normalization loops without changing significands, exponents or rounding.
SIRIUS_PB32_INLINE int PB32LeadingBit(uint a) {
#if defined(__cplusplus)
    return 31 - int(std::countl_zero(a));
#else
    return int(firstbithigh(a));
#endif
}

SIRIUS_PB32_INLINE uint PB32Abs(uint a) { return a & 0x7fffffffu; }
SIRIUS_PB32_INLINE uint PB32Negate(uint a) { return a ^ 0x80000000u; }
SIRIUS_PB32_INLINE bool PB32Nan(uint a) { return PB32Abs(a) > 0x7f800000u; }
SIRIUS_PB32_INLINE bool PB32Infinity(uint a) { return PB32Abs(a) == 0x7f800000u; }
SIRIUS_PB32_INLINE bool PB32Zero(uint a) { return PB32Abs(a) == 0u; }
SIRIUS_PB32_INLINE uint PB32Classify(uint a) {
    uint x = PB32Abs(a);
    return x == 0u ? 0u : x < 0x00800000u ? 1u : x < 0x7f800000u ? 2u : x == 0x7f800000u ? 3u : 4u;
}
SIRIUS_PB32_INLINE bool PB32Equal(uint a, uint b) {
    return !PB32Nan(a) && !PB32Nan(b) && (a == b || (PB32Zero(a) && PB32Zero(b)));
}
SIRIUS_PB32_INLINE bool PB32Less(uint a, uint b) {
    if (PB32Nan(a) || PB32Nan(b) || PB32Equal(a, b)) return false;
    uint sa = a >> 31, sb = b >> 31;
    return sa != sb ? sa != 0u : sa != 0u ? PB32Abs(a) > PB32Abs(b) : PB32Abs(a) < PB32Abs(b);
}
SIRIUS_PB32_INLINE bool PB32LessEqual(uint a, uint b) { return PB32Equal(a, b) || PB32Less(a, b); }

SIRIUS_PB32_INLINE uint PB32ShiftJam(uint a, uint distance) {
    if (distance == 0u) return a;
    if (distance >= 32u) return a != 0u ? 1u : 0u;
    return (a >> distance) | ((a << (32u - distance)) != 0u ? 1u : 0u);
}
SIRIUS_PB32_INLINE PB32Wide PB32WideProduct(uint a, uint b) {
    uint a0 = a & 0xffffu, a1 = a >> 16;
    uint b0 = b & 0xffffu, b1 = b >> 16;
    uint w0 = a0 * b0;
    uint t = a1 * b0 + (w0 >> 16);
    uint w1 = t & 0xffffu, w2 = t >> 16;
    w1 = a0 * b1 + w1;
    PB32Wide r;
    r.hi = a1 * b1 + w2 + (w1 >> 16);
    r.lo = (w1 << 16) | (w0 & 0xffffu);
    return r;
}
SIRIUS_PB32_INLINE uint PB32WideShiftJam(PB32Wide a, uint distance) {
    // Callers ensure that the retained quotient fits in one word.
    if (distance == 0u) return a.lo;
    if (distance < 32u)
        return (a.lo >> distance) | (a.hi << (32u - distance)) |
               ((a.lo << (32u - distance)) != 0u ? 1u : 0u);
    if (distance == 32u) return a.hi | (a.lo != 0u ? 1u : 0u);
    if (distance < 64u)
        return (a.hi >> (distance - 32u)) | (((a.hi << (64u - distance)) | a.lo) != 0u ? 1u : 0u);
    return (a.lo | a.hi) != 0u ? 1u : 0u;
}
SIRIUS_PB32_INLINE PB32Finite PB32Unpack(uint a) {
    PB32Finite r;
    uint e = (a >> 23) & 0xffu;
    r.significand = a & 0x007fffffu;
    r.exponent = int(e) - 127;
    if (e != 0u)
        r.significand |= 0x00800000u;
    else {
        r.exponent = -126;
        if (r.significand != 0u) {
            uint distance = uint(23 - PB32LeadingBit(r.significand));
            r.significand <<= distance;
            r.exponent -= int(distance);
        }
    }
    return r;
}
SIRIUS_PB32_INLINE uint PB32Round(uint sign, int exponent, uint extended) {
    // extended has a leading bit at bit 26 and three rounding bits.
    if (extended == 0u) return sign << 31;
    if (exponent < -126) {
        extended = PB32ShiftJam(extended, uint(-126 - exponent));
        exponent = -126;
    }
    uint residue = extended & 7u;
    uint sig = extended >> 3;
    if (residue > 4u || (residue == 4u && (sig & 1u) != 0u)) ++sig;
    if (sig >= 0x01000000u) {
        sig >>= 1;
        ++exponent;
    }
    if (exponent > 127) return (sign << 31) | 0x7f800000u;
    if (sig == 0u) return sign << 31;
    uint e = exponent == -126 && sig < 0x00800000u ? 0u : uint(exponent + 127);
    return (sign << 31) | (e << 23) | (sig & 0x007fffffu);
}
SIRIUS_PB32_INLINE uint PB32Add(uint a, uint b) {
    if (PB32Nan(a) || PB32Nan(b)) return 0x7fc00000u;
    if (PB32Infinity(a) || PB32Infinity(b)) {
        if (PB32Infinity(a) && PB32Infinity(b) && ((a ^ b) >> 31) != 0u) return 0x7fc00000u;
        return PB32Infinity(a) ? a : b;
    }
    if (PB32Zero(a) && PB32Zero(b)) return (a & b) & 0x80000000u;
    if (PB32Zero(a)) return b;
    if (PB32Zero(b)) return a;
    if (PB32Abs(a) < PB32Abs(b)) {
        uint t = a;
        a = b;
        b = t;
    }
    PB32Finite aa = PB32Unpack(a), bb = PB32Unpack(b);
    uint ax = aa.significand << 3;
    uint bx = PB32ShiftJam(bb.significand << 3, uint(aa.exponent - bb.exponent));
    uint sign = a >> 31, result;
    int exponent = aa.exponent;
    if (((a ^ b) >> 31) == 0u) {
        result = ax + bx;
        if ((result & 0x08000000u) != 0u) {
            result = PB32ShiftJam(result, 1u);
            ++exponent;
        }
    } else {
        result = ax - bx;
        if (result == 0u) return 0u;
        uint distance = uint(26 - PB32LeadingBit(result));
        result <<= distance;
        exponent -= int(distance);
    }
    return PB32Round(sign, exponent, result);
}
SIRIUS_PB32_INLINE uint PB32Subtract(uint a, uint b) { return PB32Add(a, PB32Negate(b)); }
SIRIUS_PB32_INLINE uint PB32Multiply(uint a, uint b) {
    uint sign = (a ^ b) >> 31;
    if (PB32Nan(a) || PB32Nan(b)) return 0x7fc00000u;
    if (PB32Infinity(a) || PB32Infinity(b))
        return PB32Zero(a) || PB32Zero(b) ? 0x7fc00000u : (sign << 31) | 0x7f800000u;
    if (PB32Zero(a) || PB32Zero(b)) return sign << 31;
    PB32Finite aa = PB32Unpack(a), bb = PB32Unpack(b);
    PB32Wide product = PB32WideProduct(aa.significand, bb.significand);
    bool upper = (product.hi & 0x00008000u) != 0u;
    uint extended = PB32WideShiftJam(product, upper ? 21u : 20u);
    return PB32Round(sign, aa.exponent + bb.exponent + (upper ? 1 : 0), extended);
}
SIRIUS_PB32_INLINE uint PB32Divide(uint a, uint b) {
    uint sign = (a ^ b) >> 31;
    if (PB32Nan(a) || PB32Nan(b)) return 0x7fc00000u;
    if ((PB32Infinity(a) && PB32Infinity(b)) || (PB32Zero(a) && PB32Zero(b))) return 0x7fc00000u;
    if (PB32Infinity(a) || PB32Zero(b)) return (sign << 31) | 0x7f800000u;
    if (PB32Zero(a) || PB32Infinity(b)) return sign << 31;
    PB32Finite aa = PB32Unpack(a), bb = PB32Unpack(b);
    int exponent = aa.exponent - bb.exponent;
    uint remainder = aa.significand;
    if (remainder < bb.significand) {
        remainder <<= 1;
        --exponent;
    }
    // Normalization gives d <= remainder < 2*d. Emit the leading bit,
    // then the same 26 fractional bits in exact radix-256/4 chunks.
    uint quotient = 1u;
    remainder -= bb.significand;
    for (uint chunk = 0u; chunk < 4u; ++chunk) {
        const uint bits = chunk == 3u ? 2u : 8u;
        const uint numerator = remainder << bits;
        const uint digit = numerator / bb.significand;
        remainder = numerator - digit * bb.significand;
        quotient = (quotient << bits) | digit;
    }
    if (remainder != 0u) quotient |= 1u;
    return PB32Round(sign, exponent, quotient);
}
SIRIUS_PB32_INLINE uint PB32Sqrt(uint a) {
    if (PB32Nan(a)) return 0x7fc00000u;
    if (PB32Zero(a)) return a;
    if ((a >> 31) != 0u) return 0x7fc00000u;
    if (PB32Infinity(a)) return a;
    PB32Finite aa = PB32Unpack(a);
    uint odd = uint(aa.exponent) & 1u;
    int exponent = (aa.exponent - int(odd)) / 2;
    uint shift = 29u + odd;
    PB32Wide radicand;
    radicand.lo = aa.significand << shift;
    radicand.hi = aa.significand >> (32u - shift);
    uint remainder = 0u, root = 0u;
    for (int i = 26; i >= 0; --i) {
        uint bit = uint(i) * 2u;
        uint pair = bit >= 32u ? (radicand.hi >> (bit - 32u)) & 3u : (radicand.lo >> bit) & 3u;
        remainder = (remainder << 2) | pair;
        uint trial = (root << 2) | 1u;
        root <<= 1;
        if (remainder >= trial) {
            remainder -= trial;
            root |= 1u;
        }
    }
    if (remainder != 0u) root |= 1u;
    return PB32Round(0u, exponent, root);
}
SIRIUS_PB32_INLINE uint PB32FromUnsigned(uint a) {
    if (a == 0u) return 0u;
    int exponent = PB32LeadingBit(a);
    uint extended =
        exponent <= 26 ? a << uint(26 - exponent) : PB32ShiftJam(a, uint(exponent - 26));
    return PB32Round(0u, exponent, extended);
}
SIRIUS_PB32_INLINE uint PB32FromSignedBits(uint a) {
    uint sign = a >> 31;
    uint magnitude = sign != 0u ? 0u - a : a;
    return PB32FromUnsigned(magnitude) | (sign << 31);
}

// The finite product split used by the fp64-products route: high is RTE(a*b),
// low is RTE(exact(a*b)-high). A binary32 product has at most 48 significant
// bits; for these bounds its nonzero value and nonzero exact difference are
// normal binary64 values. Computing that difference as integers preserves the same value split
// without requiring any native float control or a 64-bit shader integer.
struct PB32Product {
    uint high;
    uint low;
    uint valid;
};

SIRIUS_PB32_INLINE bool PB32WideLess(PB32Wide a, PB32Wide b) {
    return a.hi != b.hi ? a.hi < b.hi : a.lo < b.lo;
}
SIRIUS_PB32_INLINE PB32Wide PB32WideSubtract(PB32Wide a, PB32Wide b) {
    // Precondition: a >= b. Unsigned subtraction and its explicit borrow are
    // exact on the two-word magnitude.
    PB32Wide r;
    r.lo = a.lo - b.lo;
    r.hi = a.hi - b.hi - (a.lo < b.lo ? 1u : 0u);
    return r;
}
SIRIUS_PB32_INLINE uint PB32RoundWide(uint sign, int leastBitExponent, PB32Wide magnitude) {
    if ((magnitude.lo | magnitude.hi) == 0u) return sign << 31;
    int leading =
        magnitude.hi != 0u ? 32 + PB32LeadingBit(magnitude.hi) : PB32LeadingBit(magnitude.lo);
    uint extended;
    if (leading > 26)
        extended = PB32WideShiftJam(magnitude, uint(leading - 26));
    else
        extended = magnitude.lo << uint(26 - leading);
    return PB32Round(sign, leastBitExponent + leading, extended);
}
SIRIUS_PB32_INLINE PB32Product PB32ProductResidual(uint a, uint b) {
    PB32Product r;
    r.high = 0u;
    r.low = 0u;
    r.valid = 0u;
    // This is a bounded primitive, not a new policy for out-of-domain inputs.
    if (PB32Abs(a) > 0x5d800000u || PB32Abs(b) > 0x5d800000u) return r;
    if (PB32Zero(a) || PB32Zero(b)) {
        // Exact cancellation in product - high has positive zero under RTE,
        // including (-0)-(-0). The signed high product remains explicit.
        r.high = (a ^ b) & 0x80000000u;
        r.valid = 1u;
        return r;
    }
    PB32Finite aa = PB32Unpack(a), bb = PB32Unpack(b);
    int sumExponent = aa.exponent + bb.exponent;
    // exact |a*b| = product * 2^(sumExponent-46).
    PB32Wide product = PB32WideProduct(aa.significand, bb.significand);
    // Same finite nonzero high projection as PB32Multiply, sharing the exact product.
    bool upper = (product.hi & 0x00008000u) != 0u;
    uint extended = PB32WideShiftJam(product, upper ? 21u : 20u);
    r.high = PB32Round((a ^ b) >> 31, sumExponent + (upper ? 1 : 0), extended);
    PB32Wide rounded;
    rounded.lo = 0u;
    rounded.hi = 0u;
    if (!PB32Zero(r.high)) {
        PB32Finite high = PB32Unpack(r.high);
        // Reconstruct high on the exact product's integer lattice. For a
        // nonzero rounded binary32 product this shift is 23, 24 or 25, also
        // when high is subnormal or a rounding carry crosses a binade.
        int shift = high.exponent - sumExponent + 23;
        if (shift < 23 || shift > 25) {
            r.high = 0u;
            return r;
        }
        rounded.lo = high.significand << uint(shift);
        rounded.hi = high.significand >> uint(32 - shift);
    }
    bool roundedUp = PB32WideLess(product, rounded);
    PB32Wide difference =
        roundedUp ? PB32WideSubtract(rounded, product) : PB32WideSubtract(product, rounded);
    if ((difference.lo | difference.hi) != 0u) {
        uint sign = ((a ^ b) >> 31) ^ (roundedUp ? 1u : 0u);
        r.low = PB32RoundWide(sign, sumExponent - 46, difference);
    }
    // An exact zero residual is +0; nonzero residuals that underflow retain
    // their own sign, which can differ from the high product's sign.
    r.valid = 1u;
    return r;
}

#undef SIRIUS_PB32_INLINE
#if defined(__cplusplus)
}  // namespace sirius::portable_binary32
#endif
#endif  // SIRIUS_KERNELS_PORTABLE_BINARY32_H
