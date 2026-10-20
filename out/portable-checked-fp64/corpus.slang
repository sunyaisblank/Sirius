#include "portable_binary32.h"
uint Evaluate(uint op, uint a, uint b) {
    if (op == 0) return PB32Add(a, b);
    if (op == 1) return PB32Subtract(a, b);
    if (op == 2) return PB32Multiply(a, b);
    if (op == 3) return PB32Divide(a, b);
    if (op == 4) return PB32Sqrt(a);
    if (op == 5) return PB32Equal(a, b) ? 1u : 0u;
    if (op == 6) return PB32Less(a, b) ? 1u : 0u;
    if (op == 7) return PB32LessEqual(a, b) ? 1u : 0u;
    if (op == 8) return PB32FromUnsigned(a);
    if (op == 9) return PB32FromSignedBits(a);
    if (op == 10) return PB32Abs(a);
    if (op == 11) return PB32Negate(a);
    if (op == 12) return PB32Classify(a);
    return 0xffffffffu;
}

[[vk::binding(0,0)]] StructuredBuffer<uint> inputs;
[[vk::binding(1,0)]] RWStructuredBuffer<uint> outputs;
[shader("compute")]
[numthreads(64,1,1)]
void main(uint3 id:SV_DispatchThreadID) {
 uint n, stride, m, other; inputs.GetDimensions(n,stride); outputs.GetDimensions(m,other);
 if(n<1u||id.y!=0u||id.z!=0u)return;
 uint count=inputs[0]; if(count>(n-1u)/3u||count>m/3u||id.x>=count)return;
 uint o=1u+3u*id.x; uint op=inputs[o], a=inputs[o+1u], b=inputs[o+2u];
 PB32Product r; if(op==13u)r=PB32ProductResidual(a,b); else {r.high=Evaluate(op,a,b);r.low=0u;r.valid=1u;}
 outputs[3u*id.x]=r.high;outputs[3u*id.x+1u]=r.low;outputs[3u*id.x+2u]=r.valid;
}
