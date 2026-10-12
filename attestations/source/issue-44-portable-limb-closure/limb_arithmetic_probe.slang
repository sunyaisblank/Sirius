import retained_scalar;
import retained_pair;
import retained_triple;
[[vk::binding(0, 0)]] StructuredBuffer<uint> words;
[[vk::binding(1, 0)]] RWStructuredBuffer<uint> results;
[numthreads(64, 1, 1)]
void ComputeMain(uint3 id : SV_DispatchThreadID) {
    uint row=id.x;
    if(row>=words[0]) return;
    uint input=1u+16u*row, output=8u*row, op=words[input];
    RetainedTriple value=RTInvalid();
    if(op==3u) {
        RetainedTriple a;
        a.hi=RPFromBits(words[input+1u]); a.lo=RPFromBits(words[input+2u]);
        a.tail=RPFromBits(words[input+3u]); a.error=RPFromBits(words[input+4u]);
        a.valid=words[input+9u];
        value=RTAccumulateTerm(a,RPFromBits(words[input+11u]));
    } else if(op<=2u) {
        RetainedTriple a=RTPack(RPFromBits(words[input+1u]),RPFromBits(words[input+2u]),
            RPFromBits(words[input+3u]),RPFromBits(words[input+4u]));
        RetainedTriple b=RTPack(RPFromBits(words[input+5u]),RPFromBits(words[input+6u]),
            RPFromBits(words[input+7u]),RPFromBits(words[input+8u]));
        if(op==0u) value=RTDivide(a,b);
        else {
            if(words[input+9u]==0u) a.valid=0u;
            if(words[input+10u]==0u) b.valid=0u;
            value=op==1u?RTAdd(a,b):RTMultiply(a,b);
        }
    }
    results[output]=RPBits(value.hi); results[output+1u]=RPBits(value.lo);
    results[output+2u]=RPBits(value.tail); results[output+3u]=RPBits(value.error);
    results[output+4u]=value.valid;
    results[output+5u]=0u; results[output+6u]=0u; results[output+7u]=0u;
}
