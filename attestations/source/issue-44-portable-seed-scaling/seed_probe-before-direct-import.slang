import retained_pair;
[[vk::binding(0, 0)]] StructuredBuffer<uint> words;
[[vk::binding(1, 0)]] RWStructuredBuffer<uint> results;
[numthreads(64, 1, 1)]
void ComputeMain(uint3 id : SV_DispatchThreadID) {
    uint row = id.x;
    if (row >= words[0]) return;
    uint op = words[1u + 4u * row], a = words[2u + 4u * row], b = words[3u + 4u * row];
    uint high = 0u, low = 0u, valid = 1u;
    switch (op) {
        case 0u: high = PB32Add(a, b); break;
        case 1u: high = PB32Subtract(a, b); break;
        case 2u: high = PB32Multiply(a, b); break;
        case 3u: high = PB32Divide(a, b); break;
        case 4u: high = PB32Sqrt(a); break;
        case 5u: high = PB32Equal(a, b) ? 1u : 0u; break;
        case 6u: high = PB32Less(a, b) ? 1u : 0u; break;
        case 7u: high = PB32LessEqual(a, b) ? 1u : 0u; break;
        case 8u: high = PB32FromUnsigned(a); break;
        case 9u: high = PB32FromSignedBits(a); break;
        case 10u: high = PB32Abs(a); break;
        case 11u: high = PB32Negate(a); break;
        case 12u: high = PB32Classify(a); break;
        case 13u: { PB32Product p = PB32ProductResidual(a, b); high=p.high; low=p.low; valid=p.valid; break; }
        case 14u:
#ifdef SIRIUS_RETAINED_SCALE_SEED
            high=PB32ScalePowerOfTwo(a,int(b));
#else
            { int e=int(b); uint factor=e>=-126 ? uint(e+127)<<23 : 1u<<uint(e+149); high=PB32Multiply(a,factor); }
#endif
            break;
        case 15u: high=RPDivideSeed(a,b); break;
        case 16u: high=RPSqrtSeed(a); break;
        default: valid=0u; break;
    }
    results[4u * row] = high;
    results[4u * row + 1u] = low;
    results[4u * row + 2u] = valid;
    results[4u * row + 3u] = 0u;
}
