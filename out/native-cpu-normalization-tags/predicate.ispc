typedef uint32 uint32_t;
#include "normalization_tag.isph"
export void sirius_native_tag_predicate8(uniform uint32_t* uniform scratch,
    uniform uint32_t* uniform tags, uniform uint32_t index, uniform uint32_t mask) {
    if ((mask & (1U << programIndex)) != 0U) {
        uniform uint32_t base=index*5U*programCount;
        tags[index*programCount+programIndex]=SiriusNormalizationTag_0(
            scratch[base+programIndex],scratch[base+programCount+programIndex],
            scratch[base+2U*programCount+programIndex],scratch[base+3U*programCount+programIndex],
            scratch[base+4U*programCount+programIndex]);
    }
}
