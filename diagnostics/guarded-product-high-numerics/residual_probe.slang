#include "portable_binary32.h"
[[vk::binding(0, 0)]] StructuredBuffer<uint> words;
[[vk::binding(1, 0)]] RWStructuredBuffer<uint> results;
[numthreads(64, 1, 1)]
void ComputeMain(uint3 id : SV_DispatchThreadID) {
    uint row = id.x;
    if (row >= words[0]) return;
    uint offset = 1u + 4u * row;
    PB32Product product;
    product.high = 0u; product.low = 0u; product.valid = 0u;
    if (words[offset] == 13u)
        product = PB32ProductResidual(words[offset + 1u], words[offset + 2u]);
    results[4u * row] = product.high;
    results[4u * row + 1u] = product.low;
    results[4u * row + 2u] = product.valid;
    results[4u * row + 3u] = 0u;
}
