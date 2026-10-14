#include "portable_binary32.h"
[[vk::binding(0,0)]] StructuredBuffer<uint> inputs;
[[vk::binding(1,0)]] RWStructuredBuffer<uint> outputs;
[numthreads(64,1,1)]
void ComputeMain(uint3 id : SV_DispatchThreadID) {
    uint row = id.x;
    if (row >= inputs[0]) return;
    uint offset = 1 + 4 * row, op = inputs[offset], a = inputs[offset+1], b = inputs[offset+2];
    PB32Product r; r.high=0; r.low=0; r.valid=1;
    if(op==0)r.high=PB32Add(a,b);
    else if(op==1)r.high=PB32Subtract(a,b);
    else if(op==2)r.high=PB32Multiply(a,b);
    else if(op==3)r.high=PB32Divide(a,b);
    else if(op==4)r.high=PB32Sqrt(a);
    else if(op==5)r.high=PB32Equal(a,b)?1:0;
    else if(op==6)r.high=PB32Less(a,b)?1:0;
    else if(op==7)r.high=PB32LessEqual(a,b)?1:0;
    else if(op==8)r.high=PB32FromUnsigned(a);
    else if(op==9)r.high=PB32FromSignedBits(a);
    else if(op==10)r.high=PB32Abs(a);
    else if(op==11)r.high=PB32Negate(a);
    else if(op==12)r.high=PB32Classify(a);
    else if(op==13)r=PB32ProductResidual(a,b);
    else r.valid=0;
    outputs[4*row]=r.high;outputs[4*row+1]=r.low;outputs[4*row+2]=r.valid;outputs[4*row+3]=0;
}
