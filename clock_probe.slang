// Compiler admission only. This module is never loaded or dispatched.
[[vk::binding(0, 0)]] RWStructuredBuffer<uint> observed;
[shader("compute")]
[numthreads(1, 1, 1)]
void ComputeMain() {
    uint2 value = getRealtimeClock();
    observed[0] = value.x;
    observed[1] = value.y;
}
