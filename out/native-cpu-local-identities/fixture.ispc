typedef uint32 uint32_t;
#include "/home/astra/.project/Sirius/out/native-cpu-local-identities/identity_store.isph"
extern "C" void sirius_native_tagged_add8(uniform uint32_t* uniform bank,uniform uint32_t* uniform tags,uniform uint32_t d,uniform uint32_t a,uniform uint32_t b,uniform uint32_t mask);
export void local_tagged_add8(uniform uint32_t* uniform bank,uniform uint32_t* uniform tags,uniform uint32_t d,uniform uint32_t a,uniform uint32_t b,uniform uint32_t mask) {
    if((mask&(1U<<programIndex))!=0U) {
        if(!SiriusTryIdentity_0(bank,tags,d,a,b,2U))
            sirius_native_tagged_add8(bank,tags,d,a,b,(uniform uint32_t)lanemask());
    }
}
extern "C" void sirius_native_tagged_subtract8(uniform uint32_t* uniform bank,uniform uint32_t* uniform tags,uniform uint32_t d,uniform uint32_t a,uniform uint32_t b,uniform uint32_t mask);
export void local_tagged_subtract8(uniform uint32_t* uniform bank,uniform uint32_t* uniform tags,uniform uint32_t d,uniform uint32_t a,uniform uint32_t b,uniform uint32_t mask) {
    if((mask&(1U<<programIndex))!=0U) {
        if(!SiriusTryIdentity_0(bank,tags,d,a,b,3U))
            sirius_native_tagged_subtract8(bank,tags,d,a,b,(uniform uint32_t)lanemask());
    }
}
extern "C" void sirius_native_tagged_multiply8(uniform uint32_t* uniform bank,uniform uint32_t* uniform tags,uniform uint32_t d,uniform uint32_t a,uniform uint32_t b,uniform uint32_t mask);
export void local_tagged_multiply8(uniform uint32_t* uniform bank,uniform uint32_t* uniform tags,uniform uint32_t d,uniform uint32_t a,uniform uint32_t b,uniform uint32_t mask) {
    if((mask&(1U<<programIndex))!=0U) {
        if(!SiriusTryIdentity_0(bank,tags,d,a,b,4U))
            sirius_native_tagged_multiply8(bank,tags,d,a,b,(uniform uint32_t)lanemask());
    }
}
