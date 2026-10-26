// One positive exact-residual case. No scalar reference is compiled here.
#include <array>
#include <cstdint>
#include <iostream>
#include <immintrin.h>

extern "C" void sirius_native_dense_multiply8(std::uint32_t*,std::uint32_t,
    std::uint32_t,std::uint32_t,std::uint32_t) noexcept;
extern "C" void sirius_native_v3_dense_multiply8(std::uint32_t*,std::uint32_t,
    std::uint32_t,std::uint32_t,std::uint32_t) noexcept;

int main() {
    const bool v3=__builtin_cpu_supports("x86-64-v3");
    const bool v4=__builtin_cpu_supports("x86-64-v4");
    if(!v3 || !v4) {
        std::cout<<"{\"passed\":false,\"unsupported_actual_CPU\":true}\n";
        return 2;
    }
    constexpr unsigned prefix=8,bank_words=120,suffix=8;
    using Storage=std::array<std::uint32_t,prefix+bank_words+suffix>;
    Storage actual;
    actual.fill(0xa5a5a5a5U);
    auto* bank=actual.data()+prefix;
    for(unsigned reg=0;reg<2;++reg) for(unsigned lane=0;lane<8;++lane) {
        bank[reg*40+lane]=0x3f800001U;
        for(unsigned field=1;field<4;++field) bank[reg*40+field*8+lane]=0;
        bank[reg*40+32+lane]=1;
    }
    const Storage input=actual;
    Storage reference=actual;
    const unsigned saved=_mm_getcsr();
    const unsigned scope=(saved|0x1f80U)&~(0xe040U|0x3fU);
    _mm_setcsr(scope);
    sirius_native_dense_multiply8(bank,2,0,1,255);
    const unsigned candidate_mxcsr=_mm_getcsr();
    _mm_setcsr(scope);
    sirius_native_v3_dense_multiply8(reference.data()+prefix,2,0,1,255);
    const unsigned reference_mxcsr=_mm_getcsr();
    bool equal=true,immutable=true,canaries=true,center=true;
    for(unsigned i=0;i<actual.size();++i) equal=equal && actual[i]==reference[i];
    for(unsigned i=0;i<80;++i)
        immutable=immutable && bank[i]==input[prefix+i] && reference[prefix+i]==input[prefix+i];
    for(unsigned i=0;i<prefix;++i)
        canaries=canaries && actual[i]==input[i] && reference[i]==input[i];
    for(unsigned i=prefix+bank_words;i<actual.size();++i)
        canaries=canaries && actual[i]==input[i] && reference[i]==input[i];
    // (1+2^-23)^2 = 1+2^-22+2^-46, already normalized as these two words.
    for(unsigned lane=0;lane<8;++lane)
        center=center && bank[80+lane]==0x3f800002U && bank[88+lane]==0x28800000U &&
               bank[96+lane]==0 && bank[112+lane]==1;
    const bool controls=(candidate_mxcsr&~0x3fU)==scope && (reference_mxcsr&~0x3fU)==scope;
    const bool sticky=candidate_mxcsr==reference_mxcsr;
    _mm_setcsr(saved);
    const bool restored=_mm_getcsr()==saved;
    const bool passed=equal && immutable && canaries && center && controls && sticky && restored;
    std::cout<<std::boolalpha<<"{\"passed\":"<<passed
        <<",\"candidate_calls\":1,\"original_v3_reference_calls\":1,\"active_lane_results\":8"
        <<",\"compared_bank_words\":120,\"compared_canary_words\":16"
        <<",\"whole_bank_equal\":"<<equal<<",\"input_immutable\":"<<immutable
        <<",\"canaries_preserved\":"<<canaries<<",\"independent_center_equal\":"<<center
        <<",\"actual_CPU_v3\":"<<v3<<",\"actual_CPU_v4\":"<<v4
        <<",\"MXCSR_controls_unchanged\":"<<controls<<",\"MXCSR_sticky_equal\":"<<sticky
        <<",\"candidate_MXCSR\":"<<candidate_mxcsr<<",\"reference_MXCSR\":"<<reference_mxcsr
        <<",\"caller_MXCSR_state_restored\":"<<restored<<"}\n";
    return passed ? 0 : 3;
}
