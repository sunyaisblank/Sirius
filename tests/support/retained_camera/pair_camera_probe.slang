import retained_pair;
import retained_camera_direction;
[[vk::binding(0,0)]] StructuredBuffer<uint> inputs;
[[vk::binding(1,0)]] RWStructuredBuffer<uint> outputs;
void ExportPair(uint offset,RetainedPair p){outputs[offset]=asuint(p.hi);outputs[offset+1]=asuint(p.lo);outputs[offset+2]=asuint(p.error);}
[shader("compute")]
[numthreads(1,1,1)]
void ComputeMain(uint3 id:SV_DispatchThreadID){
    if(any(id!=uint3(0)))return;uint ni,si,no,so;inputs.GetDimensions(ni,si);outputs.GetDimensions(no,so);
    for(uint i=0;i<no;++i)outputs[i]=0u;if(ni!=32||no!=160)return;
    outputs[0]=asuint(271830.0f);outputs[1]=asuint(32.0f);outputs[4]=inputs[30];outputs[5]=asuint(32.0f);outputs[6]=asuint(160.0f);
#ifdef SIRIUS_FP32_COMP
    outputs[2]=asuint(1.0f);
#endif
    for(uint i=0;i<32;++i)outputs[8+i]=inputs[i];
    for(uint i=0;i<32;++i)if((inputs[i]&0x7fffffffu)>=0x7f800000u)return;
    if(asfloat(inputs[31])!=0.0f)return;
    for(uint i=20;i<=24;++i)if(asfloat(inputs[i])<=0.0f)return;
    bool thin=asfloat(inputs[29])==1.0f;
    if(!thin&&(asfloat(inputs[29])!=0.0f||asfloat(inputs[27])!=0.0f||asfloat(inputs[28])!=0.0f))return;
    RetainedPair n[3],dn[4][3],q[3],norm;
    if(!TryRetainedCameraDirection(asfloat(inputs[20]),asfloat(inputs[21]),asfloat(inputs[22]),asfloat(inputs[23]),
        asfloat(inputs[24]),asfloat(inputs[25]),asfloat(inputs[26]),asfloat(inputs[27]),asfloat(inputs[28]),thin,n,dn,q,norm))return;
    for(uint i=0;i<3;++i){ExportPair(40+3*i,n[i]);ExportPair(85+3*i,q[i]);}
    for(uint c=0;c<4;++c)for(uint i=0;i<3;++i)ExportPair(49+9*c+3*i,dn[c][i]);
    ExportPair(94,norm);outputs[3]=asuint(1.0f);
}
