#include "fp64_product.h"
[[vk::binding(0, 0)]] StructuredBuffer<uint> inputWords;
[[vk::binding(1, 0)]] RWStructuredBuffer<uint> outputWords;
[numthreads(64, 1, 1)]
void ComputeMain(uint3 id : SV_DispatchThreadID) {
    uint row = id.x;
    if (row >= inputWords[0]) return;
    uint offset = 1u + 3u * row;
    PB32Product result = PB64ProductEvaluate(inputWords[offset], inputWords[offset + 1u],
                                           inputWords[offset + 2u]);
    outputWords[3u * row] = result.high;
    outputWords[3u * row + 1u] = result.low;
    outputWords[3u * row + 2u] = result.valid;
}
