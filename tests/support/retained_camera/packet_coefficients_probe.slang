import gr_camera;
[[vk::binding(0,0)]] StructuredBuffer<float> input;
[[vk::binding(1,0)]] RWStructuredBuffer<float> output;
[shader("compute")]
[numthreads(1,1,1)]
void ComputeMain(uint3 id: SV_DispatchThreadID) {
    if(any(id!=uint3(0)))return;
    for(uint i=0;i<384;++i)output[i]=0;
    for(uint i=0;i<68;++i)output[i]=input[i];
    float row[32];
    for(uint i=0;i<4;++i)row[i]=input[3+i];
    row[4]=0;for(uint i=0;i<3;++i)row[5+i]=input[7+i];
    for(uint i=0;i<3;++i){row[8+i]=input[10+i];row[11+i]=input[16+i];row[14+i]=input[13+i];row[17+i]=input[47+i];}
    row[20]=input[0];row[21]=input[1];row[22]=input[20];
    row[23]=tan(input[19]*0.5f);row[24]=input[55];
    row[25]=input[31]+input[44];row[26]=input[32]+input[45];
    float imageX=(2.0f*row[25]/row[20]-1.0f)*row[22];
    float imageY=1.0f-2.0f*row[26]/row[21];
    ThinLensProjectionSample pupil=ProjectThinLensSample(imageX,imageY,row[23],input[53],input[54],input[55],input[66],input[67]);
    row[27]=input[52]==1?pupil.pupilRight:0;
    row[28]=input[52]==1?pupil.pupilUp:0;
    row[29]=input[52];row[30]=input[46];row[31]=0;
    for(uint i=0;i<32;++i)output[68+i]=row[i];
}
