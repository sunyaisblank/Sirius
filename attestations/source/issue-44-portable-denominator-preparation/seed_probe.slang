import retained_scalar;
import retained_pair;
[[vk::binding(0, 0)]] StructuredBuffer<uint> words;
[[vk::binding(1, 0)]] RWStructuredBuffer<uint> results;
[numthreads(64, 1, 1)]
void ComputeMain(uint3 id : SV_DispatchThreadID) {
    uint row=id.x;
    if(row>=words[0])return;
    uint input=1u+8u*row,output=8u*row;
    RPScalar a=RPFromBits(words[input+1u]),b=RPFromBits(words[input+2u]);
#ifdef SIRIUS_PREPARED_BASELINE
    RPScalar seed=RPDivideSeed(a,b);
#else
    RPDivideSeedDenominator prepared=RPPrepareDivideSeedDenominator(b);
    RPScalar seed=RPDivideSeedPrepared(a,prepared);
#endif
    results[output]=RPBits(seed);
    results[output+1u]=0u;results[output+2u]=0u;results[output+3u]=0u;
    results[output+4u]=1u;
    results[output+5u]=0u;results[output+6u]=0u;results[output+7u]=0u;
}
