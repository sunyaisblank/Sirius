#include "portable_binary32.h"
[[vk::binding(0, 0)]] StructuredBuffer<uint> inputs;
[[vk::binding(1, 0)]] RWStructuredBuffer<uint> outputs;
uint Evaluate(uint op, uint a, uint b) {
    if (op == 0u) return PB32Add(a, b);
    if (op == 1u) return PB32Subtract(a, b);
    if (op == 2u) return PB32Multiply(a, b);
    if (op == 3u) return PB32Divide(a, b);
    if (op == 4u) return PB32Sqrt(a);
    if (op == 5u) return PB32Equal(a, b) ? 1u : 0u;
    if (op == 6u) return PB32Less(a, b) ? 1u : 0u;
    if (op == 7u) return PB32LessEqual(a, b) ? 1u : 0u;
    if (op == 8u) return PB32FromUnsigned(a);
    if (op == 9u) return PB32FromSignedBits(a);
    if (op == 10u) return PB32Abs(a);
    if (op == 11u) return PB32Negate(a);
    if (op == 12u) return PB32Classify(a);
    return 0xffffffffu;
}
[shader("compute")][numthreads(64, 1, 1)] void ComputeMain(uint3 id
                                                           : SV_DispatchThreadID) {
    uint ni, si, no, so;
    inputs.GetDimensions(ni, si);
    outputs.GetDimensions(no, so);
    if (ni < 1u || id.y != 0u || id.z != 0u) return;
    uint count = inputs[0];
    if (count > (ni - 1u) / 3u || count > no / 3u || id.x >= count) return;
    uint base = 1u + id.x * 3u;
    uint op = inputs[base];
    PB32Product result;
    if (op == 13u)
        result = PB32ProductResidual(inputs[base + 1u], inputs[base + 2u]);
    else {
        result.high = Evaluate(op, inputs[base + 1u], inputs[base + 2u]);
        result.low = 0u;
        result.valid = 1u;
    }
    uint outBase = id.x * 3u;
    outputs[outBase] = result.high;
    outputs[outBase + 1u] = result.low;
    outputs[outBase + 2u] = result.valid;
}
