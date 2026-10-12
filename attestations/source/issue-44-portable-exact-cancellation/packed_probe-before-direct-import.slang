import retained_triple;
[[vk::binding(0, 0)]] StructuredBuffer<uint> words;
[[vk::binding(1, 0)]] RWStructuredBuffer<uint> results;
[numthreads(64, 1, 1)]
void ComputeMain(uint3 id : SV_DispatchThreadID) {
    uint row = id.x;
    if (row >= words[0]) return;
    uint input = 1u + 8u * row, output = 8u * row;
    uint op = words[input];
    RPScalar a = RPFromBits(words[input + 1u]);
    RPScalar b = RPFromBits(words[input + 2u]);
    RPScalar c = RPFromBits(words[input + 3u]);
    RPScalar error = RPFromBits(words[input + 4u]);
    RetainedTriple value = RTInvalid();
    if (op == 0u) {
        RetainedPair pair = RPPack(a, b, error);
        value.hi = pair.hi;
        value.lo = pair.lo;
        value.tail = RPZero();
        value.error = pair.error;
        value.valid = pair.valid;
    } else if (op == 1u) value = RTPack(a, b, c, error);
    else if (op == 2u) value = RTFinishAccumulation(a, b, c, error);
    results[output] = RPBits(value.hi);
    results[output + 1u] = RPBits(value.lo);
    results[output + 2u] = RPBits(value.tail);
    results[output + 3u] = RPBits(value.error);
    results[output + 4u] = value.valid;
    results[output + 5u] = 0u;
    results[output + 6u] = 0u;
    results[output + 7u] = 0u;
}
