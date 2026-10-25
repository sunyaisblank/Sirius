typedef uint32 uint32_t;
#define F32_asuint(x) intbits(x)
// One existing private metadata word. Known2 certifies the exact complete
// tuple's normalization and identity/valid facts. Partial writes clear it.
#include "/home/astra/.project/Sirius/out/native-cpu-normalization-tags/normalization_tag.isph"
static inline uint32_t SiriusValueFacts_0(uint32_t hi,uint32_t lo,
    uint32_t tail,uint32_t error,uint32_t valid) {
    uint32_t lower=(lo|tail|error)&2147483647U;
    bool clean=valid!=0U && lower==0U;
    uint32_t facts=2U | (valid!=0U ? 32U : 0U);
    if(clean && (hi&2147483647U)==0U) facts|=4U;
    if(clean && hi==1065353216U) facts|=8U;
    if(clean && hi==3212836864U) facts|=16U;
    return facts|SiriusNormalizationTag_0(hi,lo,tail,error,valid);
}

struct GlobalParams_0 { uniform uint32_t unused; };
struct RetainedTriple_0
{
    float hi_0;
    float lo_0;
    float tail_0;
    float error_0;
    uint32_t valid_0;
};
struct transport_old_KernelContext_0
{
    uniform GlobalParams_0* uniform globalParams_0;
    uint32_t retainedScratch_0[2295];
    uint32_t SiriusNormalized_0[459];
    uint32_t retainedStatus_0;
    uint32_t retainedCoefficients_0[125];
};
static uint32_t transport_old_RPBits_0(float value_0)
{

#line 29
    return (F32_asuint((value_0)));
}
static uint32_t transport_old_ScratchAddress_0(uint32_t index_1, uint32_t component_0, uniform uint32_t rowBase_0, uint32_t prefix_0)
{

#line 8
    return index_1 * 5U + component_0;
}
static uniform uint32_t transport_old_ScratchAddressUniform_0(uniform uint32_t index_1, uniform uint32_t component_0, uniform uint32_t rowBase_0, uniform uint32_t prefix_0)
{

#line 8
    return index_1 * 5U + component_0;
}
static void transport_old_ScratchWrite_0(uint32_t index_2, uint32_t component_1, uniform uint32_t rowBase_1, uint32_t prefix_1, uint32_t value_3, varying transport_old_KernelContext_0 * uniform kernelContext_2)
{

#line 14
    kernelContext_2->retainedScratch_0[transport_old_ScratchAddress_0(index_2, component_1, rowBase_1, prefix_1)] = value_3;
    return;
    kernelContext_2->SiriusNormalized_0[index_2] = 0U;
}
static void transport_old_ScratchWriteUniform_0(uniform uint32_t index_2, uniform uint32_t component_1, uniform uint32_t rowBase_1, uniform uint32_t prefix_1, uint32_t value_3, varying transport_old_KernelContext_0 * uniform kernelContext_2)
{

#line 14
    kernelContext_2->retainedScratch_0[transport_old_ScratchAddressUniform_0(index_2, component_1, rowBase_1, prefix_1)] = value_3;
    return;
    kernelContext_2->SiriusNormalized_0[index_2] = 0U;
}
static void transport_old_WriteValue_0(uint32_t index_10, varying RetainedTriple_0 * uniform r_10, uniform uint32_t rowBase_5, varying transport_old_KernelContext_0 * uniform kernelContext_10)
{

#line 53
    transport_old_ScratchWrite_0(index_10, 0U, rowBase_5, 2404U, transport_old_RPBits_0(r_10->hi_0), kernelContext_10);

#line 53
    transport_old_ScratchWrite_0(index_10, 1U, rowBase_5, 2404U, transport_old_RPBits_0(r_10->lo_0), kernelContext_10);

#line 53
    transport_old_ScratchWrite_0(index_10, 2U, rowBase_5, 2404U, transport_old_RPBits_0(r_10->tail_0), kernelContext_10);

#line 53
    transport_old_ScratchWrite_0(index_10, 3U, rowBase_5, 2404U, transport_old_RPBits_0(r_10->error_0), kernelContext_10);

#line 53
    transport_old_ScratchWrite_0(index_10, 4U, rowBase_5, 2404U, r_10->valid_0, kernelContext_10);

#line 59
    return;
    kernelContext_10->SiriusNormalized_0[index_10] = SiriusNormalizationTag_0(
        transport_old_RPBits_0(r_10->hi_0), transport_old_RPBits_0(r_10->lo_0), transport_old_RPBits_0(r_10->tail_0),
        transport_old_RPBits_0(r_10->error_0), r_10->valid_0);
}
static void transport_old_WriteValueUniform_0(uniform uint32_t index_10, varying RetainedTriple_0 * uniform r_10, uniform uint32_t rowBase_5, varying transport_old_KernelContext_0 * uniform kernelContext_10)
{

#line 53
    transport_old_ScratchWriteUniform_0(index_10, 0U, rowBase_5, 2404U, transport_old_RPBits_0(r_10->hi_0), kernelContext_10);

#line 53
    transport_old_ScratchWriteUniform_0(index_10, 1U, rowBase_5, 2404U, transport_old_RPBits_0(r_10->lo_0), kernelContext_10);

#line 53
    transport_old_ScratchWriteUniform_0(index_10, 2U, rowBase_5, 2404U, transport_old_RPBits_0(r_10->tail_0), kernelContext_10);

#line 53
    transport_old_ScratchWriteUniform_0(index_10, 3U, rowBase_5, 2404U, transport_old_RPBits_0(r_10->error_0), kernelContext_10);

#line 53
    transport_old_ScratchWriteUniform_0(index_10, 4U, rowBase_5, 2404U, r_10->valid_0, kernelContext_10);

#line 59
    return;
    kernelContext_10->SiriusNormalized_0[index_10] = SiriusNormalizationTag_0(
        transport_old_RPBits_0(r_10->hi_0), transport_old_RPBits_0(r_10->lo_0), transport_old_RPBits_0(r_10->tail_0),
        transport_old_RPBits_0(r_10->error_0), r_10->valid_0);
}
export void writer_transport_old(uniform uint32_t* uniform bank,
    uniform uint32_t* uniform tags, uniform uint32_t* uniform value,
    uniform uint32_t operation, uniform uint32_t index, uniform uint32_t component,
    uniform uint32_t mask, uniform uint32_t vary_index) {
    transport_old_KernelContext_0 context;
    for(uniform uint32_t r=0;r<3U;++r) {
        for(uniform uint32_t c=0;c<5U;++c)
            context.retainedScratch_0[r*5U+c]=bank[(r*5U+c)*programCount+programIndex];
        context.SiriusNormalized_0[r]=tags[r*programCount+programIndex];
    }
    if((mask&(1U<<programIndex))!=0U) {
        uint32_t varying_index=vary_index!=0U ? (index+programIndex)%3U : index;
        RetainedTriple_0 v;
        v.hi_0=floatbits(value[programIndex]);v.lo_0=floatbits(value[programCount+programIndex]);
        v.tail_0=floatbits(value[2U*programCount+programIndex]);v.error_0=floatbits(value[3U*programCount+programIndex]);
        v.valid_0=value[4U*programCount+programIndex];
        if(operation==0U) transport_old_ScratchWrite_0(varying_index,component,0U,0U,value[component*programCount+programIndex],&context);
        else if(operation==1U) transport_old_ScratchWriteUniform_0(index,component,0U,0U,value[component*programCount+programIndex],&context);
        else if(operation==2U) transport_old_WriteValue_0(varying_index,&v,0U,&context);
        else transport_old_WriteValueUniform_0(index,&v,0U,&context);
    }
    for(uniform uint32_t r=0;r<3U;++r) {
        for(uniform uint32_t c=0;c<5U;++c)
            bank[(r*5U+c)*programCount+programIndex]=context.retainedScratch_0[r*5U+c];
        tags[r*programCount+programIndex]=context.SiriusNormalized_0[r];
    }
}
struct transport_fixed_KernelContext_0
{
    uniform GlobalParams_0* uniform globalParams_0;
    uint32_t retainedScratch_0[2295];
    uint32_t SiriusNormalized_0[459];
    uint32_t retainedStatus_0;
    uint32_t retainedCoefficients_0[125];
};
static uint32_t transport_fixed_RPBits_0(float value_0)
{

#line 29
    return (F32_asuint((value_0)));
}
static uint32_t transport_fixed_ScratchAddress_0(uint32_t index_1, uint32_t component_0, uniform uint32_t rowBase_0, uint32_t prefix_0)
{

#line 8
    return index_1 * 5U + component_0;
}
static uniform uint32_t transport_fixed_ScratchAddressUniform_0(uniform uint32_t index_1, uniform uint32_t component_0, uniform uint32_t rowBase_0, uniform uint32_t prefix_0)
{

#line 8
    return index_1 * 5U + component_0;
}
static void transport_fixed_ScratchWrite_0(uint32_t index_2, uint32_t component_1, uniform uint32_t rowBase_1, uint32_t prefix_1, uint32_t value_3, varying transport_fixed_KernelContext_0 * uniform kernelContext_2)
{

#line 14
    kernelContext_2->retainedScratch_0[transport_fixed_ScratchAddress_0(index_2, component_1, rowBase_1, prefix_1)] = value_3;
    kernelContext_2->SiriusNormalized_0[index_2] = 0U;
    return;
}
static void transport_fixed_ScratchWriteUniform_0(uniform uint32_t index_2, uniform uint32_t component_1, uniform uint32_t rowBase_1, uniform uint32_t prefix_1, uint32_t value_3, varying transport_fixed_KernelContext_0 * uniform kernelContext_2)
{

#line 14
    kernelContext_2->retainedScratch_0[transport_fixed_ScratchAddressUniform_0(index_2, component_1, rowBase_1, prefix_1)] = value_3;
    kernelContext_2->SiriusNormalized_0[index_2] = 0U;
    return;
}
static void transport_fixed_WriteValue_0(uint32_t index_10, varying RetainedTriple_0 * uniform r_10, uniform uint32_t rowBase_5, varying transport_fixed_KernelContext_0 * uniform kernelContext_10)
{

#line 53
    transport_fixed_ScratchWrite_0(index_10, 0U, rowBase_5, 2404U, transport_fixed_RPBits_0(r_10->hi_0), kernelContext_10);

#line 53
    transport_fixed_ScratchWrite_0(index_10, 1U, rowBase_5, 2404U, transport_fixed_RPBits_0(r_10->lo_0), kernelContext_10);

#line 53
    transport_fixed_ScratchWrite_0(index_10, 2U, rowBase_5, 2404U, transport_fixed_RPBits_0(r_10->tail_0), kernelContext_10);

#line 53
    transport_fixed_ScratchWrite_0(index_10, 3U, rowBase_5, 2404U, transport_fixed_RPBits_0(r_10->error_0), kernelContext_10);

#line 53
    transport_fixed_ScratchWrite_0(index_10, 4U, rowBase_5, 2404U, r_10->valid_0, kernelContext_10);

#line 59
    kernelContext_10->SiriusNormalized_0[index_10] = SiriusValueFacts_0(
        transport_fixed_RPBits_0(r_10->hi_0), transport_fixed_RPBits_0(r_10->lo_0), transport_fixed_RPBits_0(r_10->tail_0),
        transport_fixed_RPBits_0(r_10->error_0), r_10->valid_0);
    return;
}
static void transport_fixed_WriteValueUniform_0(uniform uint32_t index_10, varying RetainedTriple_0 * uniform r_10, uniform uint32_t rowBase_5, varying transport_fixed_KernelContext_0 * uniform kernelContext_10)
{

#line 53
    transport_fixed_ScratchWriteUniform_0(index_10, 0U, rowBase_5, 2404U, transport_fixed_RPBits_0(r_10->hi_0), kernelContext_10);

#line 53
    transport_fixed_ScratchWriteUniform_0(index_10, 1U, rowBase_5, 2404U, transport_fixed_RPBits_0(r_10->lo_0), kernelContext_10);

#line 53
    transport_fixed_ScratchWriteUniform_0(index_10, 2U, rowBase_5, 2404U, transport_fixed_RPBits_0(r_10->tail_0), kernelContext_10);

#line 53
    transport_fixed_ScratchWriteUniform_0(index_10, 3U, rowBase_5, 2404U, transport_fixed_RPBits_0(r_10->error_0), kernelContext_10);

#line 53
    transport_fixed_ScratchWriteUniform_0(index_10, 4U, rowBase_5, 2404U, r_10->valid_0, kernelContext_10);

#line 59
    kernelContext_10->SiriusNormalized_0[index_10] = SiriusValueFacts_0(
        transport_fixed_RPBits_0(r_10->hi_0), transport_fixed_RPBits_0(r_10->lo_0), transport_fixed_RPBits_0(r_10->tail_0),
        transport_fixed_RPBits_0(r_10->error_0), r_10->valid_0);
    return;
}
export void writer_transport_fixed(uniform uint32_t* uniform bank,
    uniform uint32_t* uniform tags, uniform uint32_t* uniform value,
    uniform uint32_t operation, uniform uint32_t index, uniform uint32_t component,
    uniform uint32_t mask, uniform uint32_t vary_index) {
    transport_fixed_KernelContext_0 context;
    for(uniform uint32_t r=0;r<3U;++r) {
        for(uniform uint32_t c=0;c<5U;++c)
            context.retainedScratch_0[r*5U+c]=bank[(r*5U+c)*programCount+programIndex];
        context.SiriusNormalized_0[r]=tags[r*programCount+programIndex];
    }
    if((mask&(1U<<programIndex))!=0U) {
        uint32_t varying_index=vary_index!=0U ? (index+programIndex)%3U : index;
        RetainedTriple_0 v;
        v.hi_0=floatbits(value[programIndex]);v.lo_0=floatbits(value[programCount+programIndex]);
        v.tail_0=floatbits(value[2U*programCount+programIndex]);v.error_0=floatbits(value[3U*programCount+programIndex]);
        v.valid_0=value[4U*programCount+programIndex];
        if(operation==0U) transport_fixed_ScratchWrite_0(varying_index,component,0U,0U,value[component*programCount+programIndex],&context);
        else if(operation==1U) transport_fixed_ScratchWriteUniform_0(index,component,0U,0U,value[component*programCount+programIndex],&context);
        else if(operation==2U) transport_fixed_WriteValue_0(varying_index,&v,0U,&context);
        else transport_fixed_WriteValueUniform_0(index,&v,0U,&context);
    }
    for(uniform uint32_t r=0;r<3U;++r) {
        for(uniform uint32_t c=0;c<5U;++c)
            bank[(r*5U+c)*programCount+programIndex]=context.retainedScratch_0[r*5U+c];
        tags[r*programCount+programIndex]=context.SiriusNormalized_0[r];
    }
}
struct endpoint_old_KernelContext_0
{
    uniform GlobalParams_0* uniform globalParams_0;
    uint32_t nativeProjection_0[7];
    uint32_t retainedStatus_0;
    uint32_t retainedScratch_0[3060];
    uint32_t SiriusNormalized_0[612];
};
static uint32_t endpoint_old_RPBits_0(float value_0)
{

#line 29
    return (F32_asuint((value_0)));
}
static uint32_t endpoint_old_ScratchAddress_0(uint32_t index_0, uint32_t component_1, uniform uint32_t rowBase_4, uint32_t prefix_1)
{

#line 8
    return index_0 * 5U + component_1;
}
static uniform uint32_t endpoint_old_ScratchAddressUniform_0(uniform uint32_t index_0, uniform uint32_t component_1, uniform uint32_t rowBase_4, uniform uint32_t prefix_1)
{

#line 8
    return index_0 * 5U + component_1;
}
static void endpoint_old_ScratchWrite_0(uint32_t index_3, uint32_t component_3, uniform uint32_t rowBase_7, uint32_t prefix_3, uint32_t value_5, varying endpoint_old_KernelContext_0 * uniform kernelContext_9)
{

#line 14
    kernelContext_9->retainedScratch_0[endpoint_old_ScratchAddress_0(index_3, component_3, rowBase_7, prefix_3)] = value_5;
    return;
    kernelContext_9->SiriusNormalized_0[index_3] = 0U;
}
static void endpoint_old_ScratchWriteUniform_0(uniform uint32_t index_3, uniform uint32_t component_3, uniform uint32_t rowBase_7, uniform uint32_t prefix_3, uint32_t value_5, varying endpoint_old_KernelContext_0 * uniform kernelContext_9)
{

#line 14
    kernelContext_9->retainedScratch_0[endpoint_old_ScratchAddressUniform_0(index_3, component_3, rowBase_7, prefix_3)] = value_5;
    return;
    kernelContext_9->SiriusNormalized_0[index_3] = 0U;
}
static void endpoint_old_WriteValue_0(uint32_t index_4, varying RetainedTriple_0 * uniform r_9, uniform uint32_t rowBase_8, varying endpoint_old_KernelContext_0 * uniform kernelContext_10)
{

#line 54
    endpoint_old_ScratchWrite_0(index_4, 0U, rowBase_8, 769U, endpoint_old_RPBits_0(r_9->hi_0), kernelContext_10);

#line 54
    endpoint_old_ScratchWrite_0(index_4, 1U, rowBase_8, 769U, endpoint_old_RPBits_0(r_9->lo_0), kernelContext_10);

#line 54
    endpoint_old_ScratchWrite_0(index_4, 2U, rowBase_8, 769U, endpoint_old_RPBits_0(r_9->tail_0), kernelContext_10);

#line 54
    endpoint_old_ScratchWrite_0(index_4, 3U, rowBase_8, 769U, endpoint_old_RPBits_0(r_9->error_0), kernelContext_10);

#line 54
    endpoint_old_ScratchWrite_0(index_4, 4U, rowBase_8, 769U, r_9->valid_0, kernelContext_10);

#line 60
    return;
    kernelContext_10->SiriusNormalized_0[index_4] = SiriusNormalizationTag_0(
        endpoint_old_RPBits_0(r_9->hi_0), endpoint_old_RPBits_0(r_9->lo_0), endpoint_old_RPBits_0(r_9->tail_0),
        endpoint_old_RPBits_0(r_9->error_0), r_9->valid_0);
}
static void endpoint_old_WriteValueUniform_0(uniform uint32_t index_4, varying RetainedTriple_0 * uniform r_9, uniform uint32_t rowBase_8, varying endpoint_old_KernelContext_0 * uniform kernelContext_10)
{

#line 54
    endpoint_old_ScratchWriteUniform_0(index_4, 0U, rowBase_8, 769U, endpoint_old_RPBits_0(r_9->hi_0), kernelContext_10);

#line 54
    endpoint_old_ScratchWriteUniform_0(index_4, 1U, rowBase_8, 769U, endpoint_old_RPBits_0(r_9->lo_0), kernelContext_10);

#line 54
    endpoint_old_ScratchWriteUniform_0(index_4, 2U, rowBase_8, 769U, endpoint_old_RPBits_0(r_9->tail_0), kernelContext_10);

#line 54
    endpoint_old_ScratchWriteUniform_0(index_4, 3U, rowBase_8, 769U, endpoint_old_RPBits_0(r_9->error_0), kernelContext_10);

#line 54
    endpoint_old_ScratchWriteUniform_0(index_4, 4U, rowBase_8, 769U, r_9->valid_0, kernelContext_10);

#line 60
    return;
    kernelContext_10->SiriusNormalized_0[index_4] = SiriusNormalizationTag_0(
        endpoint_old_RPBits_0(r_9->hi_0), endpoint_old_RPBits_0(r_9->lo_0), endpoint_old_RPBits_0(r_9->tail_0),
        endpoint_old_RPBits_0(r_9->error_0), r_9->valid_0);
}
export void writer_endpoint_old(uniform uint32_t* uniform bank,
    uniform uint32_t* uniform tags, uniform uint32_t* uniform value,
    uniform uint32_t operation, uniform uint32_t index, uniform uint32_t component,
    uniform uint32_t mask, uniform uint32_t vary_index) {
    endpoint_old_KernelContext_0 context;
    for(uniform uint32_t r=0;r<3U;++r) {
        for(uniform uint32_t c=0;c<5U;++c)
            context.retainedScratch_0[r*5U+c]=bank[(r*5U+c)*programCount+programIndex];
        context.SiriusNormalized_0[r]=tags[r*programCount+programIndex];
    }
    if((mask&(1U<<programIndex))!=0U) {
        uint32_t varying_index=vary_index!=0U ? (index+programIndex)%3U : index;
        RetainedTriple_0 v;
        v.hi_0=floatbits(value[programIndex]);v.lo_0=floatbits(value[programCount+programIndex]);
        v.tail_0=floatbits(value[2U*programCount+programIndex]);v.error_0=floatbits(value[3U*programCount+programIndex]);
        v.valid_0=value[4U*programCount+programIndex];
        if(operation==0U) endpoint_old_ScratchWrite_0(varying_index,component,0U,0U,value[component*programCount+programIndex],&context);
        else if(operation==1U) endpoint_old_ScratchWriteUniform_0(index,component,0U,0U,value[component*programCount+programIndex],&context);
        else if(operation==2U) endpoint_old_WriteValue_0(varying_index,&v,0U,&context);
        else endpoint_old_WriteValueUniform_0(index,&v,0U,&context);
    }
    for(uniform uint32_t r=0;r<3U;++r) {
        for(uniform uint32_t c=0;c<5U;++c)
            bank[(r*5U+c)*programCount+programIndex]=context.retainedScratch_0[r*5U+c];
        tags[r*programCount+programIndex]=context.SiriusNormalized_0[r];
    }
}
struct endpoint_fixed_KernelContext_0
{
    uniform GlobalParams_0* uniform globalParams_0;
    uint32_t nativeProjection_0[7];
    uint32_t retainedStatus_0;
    uint32_t retainedScratch_0[3060];
    uint32_t SiriusNormalized_0[612];
};
static uint32_t endpoint_fixed_RPBits_0(float value_0)
{

#line 29
    return (F32_asuint((value_0)));
}
static uint32_t endpoint_fixed_ScratchAddress_0(uint32_t index_0, uint32_t component_1, uniform uint32_t rowBase_4, uint32_t prefix_1)
{

#line 8
    return index_0 * 5U + component_1;
}
static uniform uint32_t endpoint_fixed_ScratchAddressUniform_0(uniform uint32_t index_0, uniform uint32_t component_1, uniform uint32_t rowBase_4, uniform uint32_t prefix_1)
{

#line 8
    return index_0 * 5U + component_1;
}
static void endpoint_fixed_ScratchWrite_0(uint32_t index_3, uint32_t component_3, uniform uint32_t rowBase_7, uint32_t prefix_3, uint32_t value_5, varying endpoint_fixed_KernelContext_0 * uniform kernelContext_9)
{

#line 14
    kernelContext_9->retainedScratch_0[endpoint_fixed_ScratchAddress_0(index_3, component_3, rowBase_7, prefix_3)] = value_5;
    kernelContext_9->SiriusNormalized_0[index_3] = 0U;
    return;
}
static void endpoint_fixed_ScratchWriteUniform_0(uniform uint32_t index_3, uniform uint32_t component_3, uniform uint32_t rowBase_7, uniform uint32_t prefix_3, uint32_t value_5, varying endpoint_fixed_KernelContext_0 * uniform kernelContext_9)
{

#line 14
    kernelContext_9->retainedScratch_0[endpoint_fixed_ScratchAddressUniform_0(index_3, component_3, rowBase_7, prefix_3)] = value_5;
    kernelContext_9->SiriusNormalized_0[index_3] = 0U;
    return;
}
static void endpoint_fixed_WriteValue_0(uint32_t index_4, varying RetainedTriple_0 * uniform r_9, uniform uint32_t rowBase_8, varying endpoint_fixed_KernelContext_0 * uniform kernelContext_10)
{

#line 54
    endpoint_fixed_ScratchWrite_0(index_4, 0U, rowBase_8, 769U, endpoint_fixed_RPBits_0(r_9->hi_0), kernelContext_10);

#line 54
    endpoint_fixed_ScratchWrite_0(index_4, 1U, rowBase_8, 769U, endpoint_fixed_RPBits_0(r_9->lo_0), kernelContext_10);

#line 54
    endpoint_fixed_ScratchWrite_0(index_4, 2U, rowBase_8, 769U, endpoint_fixed_RPBits_0(r_9->tail_0), kernelContext_10);

#line 54
    endpoint_fixed_ScratchWrite_0(index_4, 3U, rowBase_8, 769U, endpoint_fixed_RPBits_0(r_9->error_0), kernelContext_10);

#line 54
    endpoint_fixed_ScratchWrite_0(index_4, 4U, rowBase_8, 769U, r_9->valid_0, kernelContext_10);

#line 60
    kernelContext_10->SiriusNormalized_0[index_4] = SiriusValueFacts_0(
        endpoint_fixed_RPBits_0(r_9->hi_0), endpoint_fixed_RPBits_0(r_9->lo_0), endpoint_fixed_RPBits_0(r_9->tail_0),
        endpoint_fixed_RPBits_0(r_9->error_0), r_9->valid_0);
    return;
}
static void endpoint_fixed_WriteValueUniform_0(uniform uint32_t index_4, varying RetainedTriple_0 * uniform r_9, uniform uint32_t rowBase_8, varying endpoint_fixed_KernelContext_0 * uniform kernelContext_10)
{

#line 54
    endpoint_fixed_ScratchWriteUniform_0(index_4, 0U, rowBase_8, 769U, endpoint_fixed_RPBits_0(r_9->hi_0), kernelContext_10);

#line 54
    endpoint_fixed_ScratchWriteUniform_0(index_4, 1U, rowBase_8, 769U, endpoint_fixed_RPBits_0(r_9->lo_0), kernelContext_10);

#line 54
    endpoint_fixed_ScratchWriteUniform_0(index_4, 2U, rowBase_8, 769U, endpoint_fixed_RPBits_0(r_9->tail_0), kernelContext_10);

#line 54
    endpoint_fixed_ScratchWriteUniform_0(index_4, 3U, rowBase_8, 769U, endpoint_fixed_RPBits_0(r_9->error_0), kernelContext_10);

#line 54
    endpoint_fixed_ScratchWriteUniform_0(index_4, 4U, rowBase_8, 769U, r_9->valid_0, kernelContext_10);

#line 60
    kernelContext_10->SiriusNormalized_0[index_4] = SiriusValueFacts_0(
        endpoint_fixed_RPBits_0(r_9->hi_0), endpoint_fixed_RPBits_0(r_9->lo_0), endpoint_fixed_RPBits_0(r_9->tail_0),
        endpoint_fixed_RPBits_0(r_9->error_0), r_9->valid_0);
    return;
}
export void writer_endpoint_fixed(uniform uint32_t* uniform bank,
    uniform uint32_t* uniform tags, uniform uint32_t* uniform value,
    uniform uint32_t operation, uniform uint32_t index, uniform uint32_t component,
    uniform uint32_t mask, uniform uint32_t vary_index) {
    endpoint_fixed_KernelContext_0 context;
    for(uniform uint32_t r=0;r<3U;++r) {
        for(uniform uint32_t c=0;c<5U;++c)
            context.retainedScratch_0[r*5U+c]=bank[(r*5U+c)*programCount+programIndex];
        context.SiriusNormalized_0[r]=tags[r*programCount+programIndex];
    }
    if((mask&(1U<<programIndex))!=0U) {
        uint32_t varying_index=vary_index!=0U ? (index+programIndex)%3U : index;
        RetainedTriple_0 v;
        v.hi_0=floatbits(value[programIndex]);v.lo_0=floatbits(value[programCount+programIndex]);
        v.tail_0=floatbits(value[2U*programCount+programIndex]);v.error_0=floatbits(value[3U*programCount+programIndex]);
        v.valid_0=value[4U*programCount+programIndex];
        if(operation==0U) endpoint_fixed_ScratchWrite_0(varying_index,component,0U,0U,value[component*programCount+programIndex],&context);
        else if(operation==1U) endpoint_fixed_ScratchWriteUniform_0(index,component,0U,0U,value[component*programCount+programIndex],&context);
        else if(operation==2U) endpoint_fixed_WriteValue_0(varying_index,&v,0U,&context);
        else endpoint_fixed_WriteValueUniform_0(index,&v,0U,&context);
    }
    for(uniform uint32_t r=0;r<3U;++r) {
        for(uniform uint32_t c=0;c<5U;++c)
            bank[(r*5U+c)*programCount+programIndex]=context.retainedScratch_0[r*5U+c];
        tags[r*programCount+programIndex]=context.SiriusNormalized_0[r];
    }
}
struct dopri_phase_old_KernelContext_0
{
    uniform GlobalParams_0* uniform globalParams_0;
    uint32_t retainedStatus_0;
    uint32_t retainedScratch_0[3145];
    uint32_t SiriusNormalized_0[629];
};
static uint32_t dopri_phase_old_RPBits_0(float value_0)
{

#line 29
    return (F32_asuint((value_0)));
}
static uint32_t dopri_phase_old_ScratchAddress_0(uint32_t index_1, uint32_t component_0, uniform uint32_t rowBase_0, uint32_t prefix_0)
{

#line 8
    return index_1 * 5U + component_0;
}
static uniform uint32_t dopri_phase_old_ScratchAddressUniform_0(uniform uint32_t index_1, uniform uint32_t component_0, uniform uint32_t rowBase_0, uniform uint32_t prefix_0)
{

#line 8
    return index_1 * 5U + component_0;
}
static void dopri_phase_old_ScratchWrite_0(uint32_t index_4, uint32_t component_2, uniform uint32_t rowBase_3, uint32_t prefix_2, uint32_t value_5, varying dopri_phase_old_KernelContext_0 * uniform kernelContext_4)
{

#line 14
    kernelContext_4->retainedScratch_0[dopri_phase_old_ScratchAddress_0(index_4, component_2, rowBase_3, prefix_2)] = value_5;
    return;
    kernelContext_4->SiriusNormalized_0[index_4] = 0U;
}
static void dopri_phase_old_ScratchWriteUniform_0(uniform uint32_t index_4, uniform uint32_t component_2, uniform uint32_t rowBase_3, uniform uint32_t prefix_2, uint32_t value_5, varying dopri_phase_old_KernelContext_0 * uniform kernelContext_4)
{

#line 14
    kernelContext_4->retainedScratch_0[dopri_phase_old_ScratchAddressUniform_0(index_4, component_2, rowBase_3, prefix_2)] = value_5;
    return;
    kernelContext_4->SiriusNormalized_0[index_4] = 0U;
}
static void dopri_phase_old_WriteValue_0(uint32_t index_5, varying RetainedTriple_0 * uniform r_7, uniform uint32_t rowBase_4, varying dopri_phase_old_KernelContext_0 * uniform kernelContext_5)
{

#line 39
    dopri_phase_old_ScratchWrite_0(index_5, 0U, rowBase_4, 1004U, dopri_phase_old_RPBits_0(r_7->hi_0), kernelContext_5);

#line 39
    dopri_phase_old_ScratchWrite_0(index_5, 1U, rowBase_4, 1004U, dopri_phase_old_RPBits_0(r_7->lo_0), kernelContext_5);

#line 39
    dopri_phase_old_ScratchWrite_0(index_5, 2U, rowBase_4, 1004U, dopri_phase_old_RPBits_0(r_7->tail_0), kernelContext_5);

#line 39
    dopri_phase_old_ScratchWrite_0(index_5, 3U, rowBase_4, 1004U, dopri_phase_old_RPBits_0(r_7->error_0), kernelContext_5);

#line 39
    dopri_phase_old_ScratchWrite_0(index_5, 4U, rowBase_4, 1004U, r_7->valid_0, kernelContext_5);

#line 45
    return;
    kernelContext_5->SiriusNormalized_0[index_5] = SiriusNormalizationTag_0(
        dopri_phase_old_RPBits_0(r_7->hi_0), dopri_phase_old_RPBits_0(r_7->lo_0), dopri_phase_old_RPBits_0(r_7->tail_0),
        dopri_phase_old_RPBits_0(r_7->error_0), r_7->valid_0);
}
static void dopri_phase_old_WriteValueUniform_0(uniform uint32_t index_5, varying RetainedTriple_0 * uniform r_7, uniform uint32_t rowBase_4, varying dopri_phase_old_KernelContext_0 * uniform kernelContext_5)
{

#line 39
    dopri_phase_old_ScratchWriteUniform_0(index_5, 0U, rowBase_4, 1004U, dopri_phase_old_RPBits_0(r_7->hi_0), kernelContext_5);

#line 39
    dopri_phase_old_ScratchWriteUniform_0(index_5, 1U, rowBase_4, 1004U, dopri_phase_old_RPBits_0(r_7->lo_0), kernelContext_5);

#line 39
    dopri_phase_old_ScratchWriteUniform_0(index_5, 2U, rowBase_4, 1004U, dopri_phase_old_RPBits_0(r_7->tail_0), kernelContext_5);

#line 39
    dopri_phase_old_ScratchWriteUniform_0(index_5, 3U, rowBase_4, 1004U, dopri_phase_old_RPBits_0(r_7->error_0), kernelContext_5);

#line 39
    dopri_phase_old_ScratchWriteUniform_0(index_5, 4U, rowBase_4, 1004U, r_7->valid_0, kernelContext_5);

#line 45
    return;
    kernelContext_5->SiriusNormalized_0[index_5] = SiriusNormalizationTag_0(
        dopri_phase_old_RPBits_0(r_7->hi_0), dopri_phase_old_RPBits_0(r_7->lo_0), dopri_phase_old_RPBits_0(r_7->tail_0),
        dopri_phase_old_RPBits_0(r_7->error_0), r_7->valid_0);
}
export void writer_dopri_phase_old(uniform uint32_t* uniform bank,
    uniform uint32_t* uniform tags, uniform uint32_t* uniform value,
    uniform uint32_t operation, uniform uint32_t index, uniform uint32_t component,
    uniform uint32_t mask, uniform uint32_t vary_index) {
    dopri_phase_old_KernelContext_0 context;
    for(uniform uint32_t r=0;r<3U;++r) {
        for(uniform uint32_t c=0;c<5U;++c)
            context.retainedScratch_0[r*5U+c]=bank[(r*5U+c)*programCount+programIndex];
        context.SiriusNormalized_0[r]=tags[r*programCount+programIndex];
    }
    if((mask&(1U<<programIndex))!=0U) {
        uint32_t varying_index=vary_index!=0U ? (index+programIndex)%3U : index;
        RetainedTriple_0 v;
        v.hi_0=floatbits(value[programIndex]);v.lo_0=floatbits(value[programCount+programIndex]);
        v.tail_0=floatbits(value[2U*programCount+programIndex]);v.error_0=floatbits(value[3U*programCount+programIndex]);
        v.valid_0=value[4U*programCount+programIndex];
        if(operation==0U) dopri_phase_old_ScratchWrite_0(varying_index,component,0U,0U,value[component*programCount+programIndex],&context);
        else if(operation==1U) dopri_phase_old_ScratchWriteUniform_0(index,component,0U,0U,value[component*programCount+programIndex],&context);
        else if(operation==2U) dopri_phase_old_WriteValue_0(varying_index,&v,0U,&context);
        else dopri_phase_old_WriteValueUniform_0(index,&v,0U,&context);
    }
    for(uniform uint32_t r=0;r<3U;++r) {
        for(uniform uint32_t c=0;c<5U;++c)
            bank[(r*5U+c)*programCount+programIndex]=context.retainedScratch_0[r*5U+c];
        tags[r*programCount+programIndex]=context.SiriusNormalized_0[r];
    }
}
struct dopri_phase_fixed_KernelContext_0
{
    uniform GlobalParams_0* uniform globalParams_0;
    uint32_t retainedStatus_0;
    uint32_t retainedScratch_0[3145];
    uint32_t SiriusNormalized_0[629];
};
static uint32_t dopri_phase_fixed_RPBits_0(float value_0)
{

#line 29
    return (F32_asuint((value_0)));
}
static uint32_t dopri_phase_fixed_ScratchAddress_0(uint32_t index_1, uint32_t component_0, uniform uint32_t rowBase_0, uint32_t prefix_0)
{

#line 8
    return index_1 * 5U + component_0;
}
static uniform uint32_t dopri_phase_fixed_ScratchAddressUniform_0(uniform uint32_t index_1, uniform uint32_t component_0, uniform uint32_t rowBase_0, uniform uint32_t prefix_0)
{

#line 8
    return index_1 * 5U + component_0;
}
static void dopri_phase_fixed_ScratchWrite_0(uint32_t index_4, uint32_t component_2, uniform uint32_t rowBase_3, uint32_t prefix_2, uint32_t value_5, varying dopri_phase_fixed_KernelContext_0 * uniform kernelContext_4)
{

#line 14
    kernelContext_4->retainedScratch_0[dopri_phase_fixed_ScratchAddress_0(index_4, component_2, rowBase_3, prefix_2)] = value_5;
    kernelContext_4->SiriusNormalized_0[index_4] = 0U;
    return;
}
static void dopri_phase_fixed_ScratchWriteUniform_0(uniform uint32_t index_4, uniform uint32_t component_2, uniform uint32_t rowBase_3, uniform uint32_t prefix_2, uint32_t value_5, varying dopri_phase_fixed_KernelContext_0 * uniform kernelContext_4)
{

#line 14
    kernelContext_4->retainedScratch_0[dopri_phase_fixed_ScratchAddressUniform_0(index_4, component_2, rowBase_3, prefix_2)] = value_5;
    kernelContext_4->SiriusNormalized_0[index_4] = 0U;
    return;
}
static void dopri_phase_fixed_WriteValue_0(uint32_t index_5, varying RetainedTriple_0 * uniform r_7, uniform uint32_t rowBase_4, varying dopri_phase_fixed_KernelContext_0 * uniform kernelContext_5)
{

#line 39
    dopri_phase_fixed_ScratchWrite_0(index_5, 0U, rowBase_4, 1004U, dopri_phase_fixed_RPBits_0(r_7->hi_0), kernelContext_5);

#line 39
    dopri_phase_fixed_ScratchWrite_0(index_5, 1U, rowBase_4, 1004U, dopri_phase_fixed_RPBits_0(r_7->lo_0), kernelContext_5);

#line 39
    dopri_phase_fixed_ScratchWrite_0(index_5, 2U, rowBase_4, 1004U, dopri_phase_fixed_RPBits_0(r_7->tail_0), kernelContext_5);

#line 39
    dopri_phase_fixed_ScratchWrite_0(index_5, 3U, rowBase_4, 1004U, dopri_phase_fixed_RPBits_0(r_7->error_0), kernelContext_5);

#line 39
    dopri_phase_fixed_ScratchWrite_0(index_5, 4U, rowBase_4, 1004U, r_7->valid_0, kernelContext_5);

#line 45
    kernelContext_5->SiriusNormalized_0[index_5] = SiriusValueFacts_0(
        dopri_phase_fixed_RPBits_0(r_7->hi_0), dopri_phase_fixed_RPBits_0(r_7->lo_0), dopri_phase_fixed_RPBits_0(r_7->tail_0),
        dopri_phase_fixed_RPBits_0(r_7->error_0), r_7->valid_0);
    return;
}
static void dopri_phase_fixed_WriteValueUniform_0(uniform uint32_t index_5, varying RetainedTriple_0 * uniform r_7, uniform uint32_t rowBase_4, varying dopri_phase_fixed_KernelContext_0 * uniform kernelContext_5)
{

#line 39
    dopri_phase_fixed_ScratchWriteUniform_0(index_5, 0U, rowBase_4, 1004U, dopri_phase_fixed_RPBits_0(r_7->hi_0), kernelContext_5);

#line 39
    dopri_phase_fixed_ScratchWriteUniform_0(index_5, 1U, rowBase_4, 1004U, dopri_phase_fixed_RPBits_0(r_7->lo_0), kernelContext_5);

#line 39
    dopri_phase_fixed_ScratchWriteUniform_0(index_5, 2U, rowBase_4, 1004U, dopri_phase_fixed_RPBits_0(r_7->tail_0), kernelContext_5);

#line 39
    dopri_phase_fixed_ScratchWriteUniform_0(index_5, 3U, rowBase_4, 1004U, dopri_phase_fixed_RPBits_0(r_7->error_0), kernelContext_5);

#line 39
    dopri_phase_fixed_ScratchWriteUniform_0(index_5, 4U, rowBase_4, 1004U, r_7->valid_0, kernelContext_5);

#line 45
    kernelContext_5->SiriusNormalized_0[index_5] = SiriusValueFacts_0(
        dopri_phase_fixed_RPBits_0(r_7->hi_0), dopri_phase_fixed_RPBits_0(r_7->lo_0), dopri_phase_fixed_RPBits_0(r_7->tail_0),
        dopri_phase_fixed_RPBits_0(r_7->error_0), r_7->valid_0);
    return;
}
export void writer_dopri_phase_fixed(uniform uint32_t* uniform bank,
    uniform uint32_t* uniform tags, uniform uint32_t* uniform value,
    uniform uint32_t operation, uniform uint32_t index, uniform uint32_t component,
    uniform uint32_t mask, uniform uint32_t vary_index) {
    dopri_phase_fixed_KernelContext_0 context;
    for(uniform uint32_t r=0;r<3U;++r) {
        for(uniform uint32_t c=0;c<5U;++c)
            context.retainedScratch_0[r*5U+c]=bank[(r*5U+c)*programCount+programIndex];
        context.SiriusNormalized_0[r]=tags[r*programCount+programIndex];
    }
    if((mask&(1U<<programIndex))!=0U) {
        uint32_t varying_index=vary_index!=0U ? (index+programIndex)%3U : index;
        RetainedTriple_0 v;
        v.hi_0=floatbits(value[programIndex]);v.lo_0=floatbits(value[programCount+programIndex]);
        v.tail_0=floatbits(value[2U*programCount+programIndex]);v.error_0=floatbits(value[3U*programCount+programIndex]);
        v.valid_0=value[4U*programCount+programIndex];
        if(operation==0U) dopri_phase_fixed_ScratchWrite_0(varying_index,component,0U,0U,value[component*programCount+programIndex],&context);
        else if(operation==1U) dopri_phase_fixed_ScratchWriteUniform_0(index,component,0U,0U,value[component*programCount+programIndex],&context);
        else if(operation==2U) dopri_phase_fixed_WriteValue_0(varying_index,&v,0U,&context);
        else dopri_phase_fixed_WriteValueUniform_0(index,&v,0U,&context);
    }
    for(uniform uint32_t r=0;r<3U;++r) {
        for(uniform uint32_t c=0;c<5U;++c)
            bank[(r*5U+c)*programCount+programIndex]=context.retainedScratch_0[r*5U+c];
        tags[r*programCount+programIndex]=context.SiriusNormalized_0[r];
    }
}
