import retained_scalar;
import retained_pair;
import retained_triple;
[[vk::binding(0, 0)]] StructuredBuffer<uint> words;
[[vk::binding(1, 0)]] RWStructuredBuffer<uint> results;
[numthreads(64, 1, 1)]
void ComputeMain(uint3 id : SV_DispatchThreadID) {
    uint row=id.x;
    if(row>=words[0])return;
    uint input=1u+16u*row,output=8u*row;
    RetainedTriple a=RTPack(RPFromBits(words[input+1u]),RPFromBits(words[input+2u]),
        RPFromBits(words[input+3u]),RPFromBits(words[input+4u]));
    RetainedTriple b=RTPack(RPFromBits(words[input+5u]),RPFromBits(words[input+6u]),
        RPFromBits(words[input+7u]),RPFromBits(words[input+8u]));
    RetainedTriple value=RTDivide(a,b);
    results[output]=RPBits(value.hi);results[output+1u]=RPBits(value.lo);
    results[output+2u]=RPBits(value.tail);results[output+3u]=RPBits(value.error);
    results[output+4u]=value.valid;
    results[output+5u]=0u;results[output+6u]=0u;results[output+7u]=0u;
}
