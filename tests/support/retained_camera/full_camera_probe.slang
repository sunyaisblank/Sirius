import retained_pair;
import retained_camera_launch;
[[vk::binding(0,0)]] StructuredBuffer<uint> inputs;
[[vk::binding(1,0)]] RWStructuredBuffer<uint> outputs;
[shader("compute")]
[numthreads(1,1,1)]
void ComputeMain(uint3 id:SV_DispatchThreadID){
    if(any(id!=uint3(0)))return;uint ni,si,no,so;inputs.GetDimensions(ni,si);outputs.GetDimensions(no,so);
    for(uint i=0;i<no;++i)outputs[i]=0u;if(ni!=32||no!=384)return;
    outputs[0]=asuint(271832.0f);outputs[1]=asuint(32.0f);outputs[4]=inputs[30];outputs[5]=asuint(32.0f);outputs[6]=asuint(384.0f);
#ifdef SIRIUS_FP32_COMP
    outputs[2]=asuint(1.0f);
#endif
    for(uint i=0;i<32;++i)outputs[8+i]=inputs[i];
    float values[32];
    for(uint i=0;i<32;++i){if((inputs[i]&0x7fffffffu)>=0x7f800000u)return;values[i]=asfloat(inputs[i]);}
    if(values[31]!=0.0f)return;
    for(uint i=20;i<=24;++i)if(values[i]<=0.0f)return;
    if(values[29]!=1.0f&&(values[29]!=0.0f||values[27]!=0.0f||values[28]!=0.0f))return;
    RetainedPair scientific[104];if(!TryRetainedCameraLaunch(values,scientific))return;
    for(uint i=0;i<104;++i){outputs[40+3*i]=asuint(scientific[i].hi);outputs[41+3*i]=asuint(scientific[i].lo);outputs[42+3*i]=asuint(scientific[i].error);}
    outputs[3]=asuint(1.0f);
}
