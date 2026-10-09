// Scalar words retain the same storage ABI in both arithmetic implementations.
// Integer arithmetic establishes RTE and gradual underflow independently of
// native float-control properties. Native arithmetic keeps its explicit order.
#ifdef SIRIUS_RETAINED_PORTABLE
#include "portable_binary32.h"
#include "fp64_scalar.h"
typedef uint RPScalar;
RPScalar RPFromBits(uint bits) { return bits; }
uint RPBits(RPScalar value) { return value; }
RPScalar RPFromSigned(int value) { return PB32FromSignedBits(uint(value)); }
RPScalar RPFromUnsigned(uint value) { return PB32FromUnsigned(value); }
bool RPScalarEqual(RPScalar a, RPScalar b) { return PB32Equal(a, b); }
bool RPScalarLess(RPScalar a, RPScalar b) { return PB32Less(a, b); }
bool RPScalarLessEqual(RPScalar a, RPScalar b) { return PB32LessEqual(a, b); }
RPScalar RPScalarNegate(RPScalar value) { return PB32Negate(value); }
RPScalar RPScalarMax(RPScalar a, RPScalar b) { return PB32Less(a, b) ? b : a; }
RPScalar OrderedAdd32(RPScalar a, RPScalar b) { return PB64Add(a, b); }
RPScalar OrderedSubtract32(RPScalar a, RPScalar b) { return PB64Subtract(a, b); }
RPScalar OrderedMultiply32(RPScalar a, RPScalar b) { return PB64Multiply(a, b); }
RPScalar OrderedDivide32(RPScalar a, RPScalar b) { return PB32Divide(a, b); }
RPScalar OrderedSqrt32(RPScalar value) { return PB32Sqrt(value); }
#else
typedef float RPScalar;
RPScalar RPFromBits(uint bits) { return asfloat(bits); }
uint RPBits(RPScalar value) { return asuint(value); }
RPScalar RPFromSigned(int value) { return float(value); }
RPScalar RPFromUnsigned(uint value) { return float(value); }
bool RPScalarEqual(RPScalar a, RPScalar b) { return a == b; }
bool RPScalarLess(RPScalar a, RPScalar b) { return a < b; }
bool RPScalarLessEqual(RPScalar a, RPScalar b) { return a <= b; }
RPScalar RPScalarNegate(RPScalar value) { return asfloat(asuint(value) ^ 0x80000000u); }
RPScalar RPScalarMax(RPScalar a, RPScalar b) { return max(a, b); }
RPScalar OrderedAdd32(RPScalar a, RPScalar b) {
    return spirv_asm {
        OpFAdd $$float %v $a $b;
        OpDecorate %v NoContraction;
        OpCopyObject $$float result %v;
    };
}
RPScalar OrderedSubtract32(RPScalar a, RPScalar b) {
    return spirv_asm {
        OpFSub $$float %v $a $b;
        OpDecorate %v NoContraction;
        OpCopyObject $$float result %v;
    };
}
RPScalar OrderedMultiply32(RPScalar a, RPScalar b) {
    return spirv_asm {
        OpFMul $$float %v $a $b;
        OpDecorate %v NoContraction;
        OpCopyObject $$float result %v;
    };
}
RPScalar OrderedDivide32(RPScalar a, RPScalar b) {
    return spirv_asm {
        OpFDiv $$float %v $a $b;
        OpDecorate %v NoContraction;
        OpCopyObject $$float result %v;
    };
}
RPScalar OrderedSqrt32(RPScalar value) { return sqrt(value); }
#endif

RPScalar RPZero() { return RPFromBits(0u); }
RPScalar RPOne() { return RPFromBits(0x3f800000u); }
bool RPScalarGreater(RPScalar a, RPScalar b) { return RPScalarLess(b, a); }
bool RPScalarGreaterEqual(RPScalar a, RPScalar b) { return RPScalarLessEqual(b, a); }
