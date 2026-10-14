import retained_scalar;
import retained_triple;
[[vk::binding(0,0)]] StructuredBuffer<uint> inputs;
[[vk::binding(1,0)]] RWStructuredBuffer<uint> outputs;
[numthreads(64,1,1)]
void ComputeMain(uint3 id : SV_DispatchThreadID) {
    uint row=id.x;
    if(row>=inputs[0]) return;
    uint input=1u+6u*row, output=5u*row;
    RetainedTriple a;
    a.hi=RPFromBits(inputs[input]);
    a.lo=RPFromBits(inputs[input+1u]);
    a.tail=RPFromBits(inputs[input+2u]);
    a.error=RPFromBits(inputs[input+3u]);
    a.valid=inputs[input+4u];
    RetainedTriple r=RTAccumulateTerm(a,RPFromBits(inputs[input+5u]));
    outputs[output]=RPBits(r.hi);
    outputs[output+1u]=RPBits(r.lo);
    outputs[output+2u]=RPBits(r.tail);
    outputs[output+3u]=RPBits(r.error);
    outputs[output+4u]=r.valid;
}
