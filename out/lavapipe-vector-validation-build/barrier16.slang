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
// Every lane reaches every original workgroup rendezvous, including tail lanes.
groupshared uint shared_result[48];
groupshared uint joined_mask;
groupshared uint joined_count;
[shader("compute")]
[numthreads(16,1,1)]
void main(uint3 group:SV_GroupID, uint lane:SV_GroupIndex) {
 uint n,stride,m,other;inputs.GetDimensions(n,stride);outputs.GetDimensions(m,other);
 if(n<1u||group.y!=0u||group.z!=0u)return;
 uint count=inputs[0];if(count>(n-1u)/3u||count>m/3u)return;
 uint base=16u*group.x;
 if(lane==0u){joined_mask=0u;joined_count=0u;}
 GroupMemoryBarrierWithGroupSync();
 uint record_lane=(5u*lane+7u)&15u;uint row=base+record_lane;
 PB32Product r;r.high=0u;r.low=0u;r.valid=0u;
 if(row<count){
  uint o=1u+3u*row;uint op=inputs[o],a=inputs[o+1u],b=inputs[o+2u];
  if(op==13u)r=PB32ProductResidual(a,b);else{r.high=Evaluate(op,a,b);r.low=0u;r.valid=1u;}
  uint ignored;
  InterlockedOr(joined_mask,1u<<record_lane,ignored);
  InterlockedAdd(joined_count,1u,ignored);
 }
 shared_result[3u*lane]=r.high;shared_result[3u*lane+1u]=r.low;shared_result[3u*lane+2u]=r.valid;
 GroupMemoryBarrierWithGroupSync();
 uint active=base<count?min(16u,count-base):0u;
 uint expected_mask=(1u<<active)-1u;
 uint source_lane=(13u*(lane+16u-7u))&15u;
 uint high=shared_result[3u*source_lane],low=shared_result[3u*source_lane+1u],valid=shared_result[3u*source_lane+2u];
 if(joined_mask!=expected_mask||joined_count!=active)valid^=0xdeadbeefu;
 // Complete all reads before lanes overwrite scratch in this independently
 // checked second cross-lane rendezvous.
 GroupMemoryBarrierWithGroupSync();
 shared_result[3u*lane]=high;shared_result[3u*lane+1u]=low;shared_result[3u*lane+2u]=valid;
 GroupMemoryBarrierWithGroupSync();
 uint rotated=(lane+9u)&15u;
 uint original=(rotated+7u)&15u;
 if(base+original<count){
  outputs[3u*(base+original)]=shared_result[3u*original];
  outputs[3u*(base+original)+1u]=shared_result[3u*original+1u];
  outputs[3u*(base+original)+2u]=shared_result[3u*original+2u];
 }
}
