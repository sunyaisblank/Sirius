import retained_pair;
[[vk::binding(0,0)]] StructuredBuffer<uint> inputs;
[[vk::binding(1,0)]] RWStructuredBuffer<uint> outputs;
[shader("compute")]
[numthreads(1,1,1)]
void ComputeMain(uint3 id:SV_DispatchThreadID){
    if(any(id!=uint3(0)))return;uint ni,si,no,so;inputs.GetDimensions(ni,si);outputs.GetDimensions(no,so);
    for(uint i=0;i<no;++i)outputs[i]=0u;if(ni!=32||no!=160)return;
    outputs[0]=asuint(271831.0f);outputs[1]=asuint(32.0f);outputs[4]=inputs[30];outputs[5]=asuint(32.0f);outputs[6]=asuint(160.0f);
#ifdef SIRIUS_FP32_COMP
    outputs[2]=asuint(1.0f);
#endif
    for(uint i=0;i<32;++i)outputs[8+i]=inputs[i];
    RetainedPair a=RPPack(asfloat(inputs[0]),asfloat(inputs[1]),asfloat(inputs[2]));
    RetainedPair b=RPPack(asfloat(inputs[3]),asfloat(inputs[4]),asfloat(inputs[5]));
    RetainedPair results[4];results[0]=RPAdd(a,b);results[1]=RPMultiply(a,b);results[2]=RPDivide(a,b);results[3]=RPSqrt(a);
    for(uint i=0;i<4;++i){uint o=40+4*i;outputs[o]=asuint(results[i].hi);outputs[o+1]=asuint(results[i].lo);outputs[o+2]=asuint(results[i].error);outputs[o+3]=asuint(float(results[i].valid));}
    outputs[3]=asuint(1.0f);
}
