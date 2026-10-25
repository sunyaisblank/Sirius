// One focused ABI/word-compatibility check against the original generated
// scalar arithmetic. This finite comparison is not an independent math oracle.
#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <limits>
#include <iostream>
#include <immintrin.h>
#include <transport.native.hpp>

namespace g = sirius_generated_transport;
using Words = std::array<std::uint32_t,5>;
using Bank = std::array<std::uint32_t,120>;
extern "C" void sirius_native_tagged_add8(std::uint32_t*,std::uint32_t*,std::uint32_t,std::uint32_t,std::uint32_t,std::uint32_t) noexcept;
extern "C" void sirius_native_tagged_multiply8(std::uint32_t*,std::uint32_t*,std::uint32_t,std::uint32_t,std::uint32_t,std::uint32_t) noexcept;
extern "C" void sirius_native_tagged_divide8(std::uint32_t*,std::uint32_t*,std::uint32_t,std::uint32_t,std::uint32_t,std::uint32_t) noexcept;
extern "C" void sirius_native_tagged_sqrt8(std::uint32_t*,std::uint32_t*,std::uint32_t,std::uint32_t,std::uint32_t,std::uint32_t) noexcept;
extern "C" void sirius_native_tagged_subtract8(std::uint32_t*,std::uint32_t*,std::uint32_t,std::uint32_t,std::uint32_t,std::uint32_t) noexcept;
extern "C" __attribute__((weak)) void sirius_native_tag_predicate8(std::uint32_t*,std::uint32_t*,std::uint32_t,std::uint32_t) noexcept;
static_assert(sizeof(g::RetainedTriple_0)==20);
using Tags=std::array<std::uint32_t,24>;
static std::uint32_t expected_tag(const Words& w) {
    const float h=std::bit_cast<float>(w[0]),l=std::bit_cast<float>(w[1]),
                t=std::bit_cast<float>(w[2]),e=std::bit_cast<float>(w[3]);
    std::uint32_t facts=2|(w[4]!=0 ? 32U : 0U);
    if(w[4]!=0 && l==0 && t==0 && e==0) {
        if(h==0) facts|=4;
        if(h==1) facts|=8;
        if(h==-1) facts|=16;
    }
    if(w[4]!=1 || !std::isfinite(h) || !std::isfinite(l) || !std::isfinite(t) ||
       !std::isfinite(e) || std::abs(h)>0x1p120f || std::signbit(e) || e>0x1p120f)
        return facts;
    const auto spacing=[](float x) {
        const float m=std::abs(x);
        return double(std::nextafter(m,std::numeric_limits<float>::infinity()))-double(m);
    };
    return facts|(std::abs(double(l))<=spacing(h) && std::abs(double(t))<=spacing(l) &&
           (h!=0 || (l==0 && t==0)) ? 1U : 0U);
}


int main(int argc,char** argv) {
    if(!__builtin_cpu_supports("x86-64-v3")) return 2;
    const unsigned saved=_mm_getcsr();
    const unsigned controls=(saved|0x1f80U)&~0xe040U;
    _mm_setcsr(controls);
    const std::array<Words,49> cases={{
        {0,0,0,0,1}, {0x80000000U,0,0,0,1},
        {0,0x80000000U,0x80000000U,0x80000000U,1},
        {0x3f800000U,0,0,0,1}, {0xbf800000U,0,0,0,1},
        {0x40000000U,0,0,0,1}, {0xc0000000U,0,0,0,1},
        {0x3f000000U,0x33000000U,0x26800000U,5,1},
        {0xbf000000U,0xb3000000U,0xa6800000U,5,1},
        {1,0,0,0,1}, {0x80000001U,0,0,0,1},
        {0x007fffffU,1,0,0,1}, {0x00800000U,0x80000001U,0,0,1},
        {1568669695U,1,0,0,1}, {1568669696U,0,0,0,1},
        {1568669697U,0,0,0,1}, {2071986175U,0,0,0,1},
        {2071986176U,0,0,0,1}, {2071986177U,0,0,0,1},
        {2105540608U,0,0,0,1}, {2105540609U,0,0,0,1},
        {0x3f800000U,0,0,0x80000000U,1},
        {0x3f800000U,0,0,2071986176U,1},
        {0x3f800000U,0,0,2071986177U,1},
        {0x7f800000U,0,0,0,1}, {0xff800000U,0,0,0,1},
        {0x7fc01234U,0,0,0,1}, {0x7f801234U,0,0,0,1},
        {0x7fc01234U,0xff800000U,0x80000000U,0x7f801234U,0},
        {0x3f800000U,0,0,0,2}, {0,0,0,0,2},
        {0x3f800000U,0x80000000U,0x80000000U,0,1},
        {0xbf800000U,0x80000000U,0,0x80000000U,1},
        {0x5a7fffffU,0x4e800000U,0x43000000U,0,1},
        {0x5a800000U,0x4f000000U,0x43800000U,0,1},
        {0x5a800001U,0x4f000000U,0x43800000U,0,1},
        {0xda800000U,0xcf000000U,0xc3800000U,0,1},
        {0x5a800000U,0xcf000000U,0x43800000U,0,1},
        {0x3f800000U,0x34000000U,0x28800000U,0x7b800000U,1},
        {0x3f800000U,0x34000001U,0,0,1},
        {0x3f800000U,0x34000000U,0x28800001U,0,1},
        {0x7b800000U,0x70000000U,0x64800000U,0,1},
        {0x7b800000U,0,0,0x7b800000U,1},
        {0xfb800000U,0,0,0x7b800000U,1},
        {0x00800000U,1,1,5,1},
        {0x00000001U,1,1,0,1},
        {0,1,0,0,1},
        {0,0,1,0,1},
        {0,0,0,0x7b800000U,1}
    }};
    const std::array<unsigned,6> masks={0,1,15,85,128,255};
    if(argc==2 && std::strcmp(argv[1],"--predicate-only")==0) {
        if(!sirius_native_tag_predicate8) { _mm_setcsr(saved);return 5; }
        std::uint64_t predicate_calls=0, predicate_tag_words=0, predicate_input_words=0;
        for(unsigned pattern=0;pattern<2;++pattern) for(unsigned p=0;p<cases.size();++p)
        for(unsigned index=0;index<3;++index) for(unsigned mask:masks) {
            Bank input;input.fill(0xa5a5a5a5U);
            Tags tags;tags.fill(0xa5a5a5a5U);Tags expected=tags;
            for(unsigned lane=0;lane<8;++lane) {
                const auto& value=cases[(p+(pattern==0 ? lane : 0))%cases.size()];
                for(unsigned f=0;f<5;++f) input[index*40+f*8+lane]=value[f];
                if(mask&(1U<<lane)) expected[index*8+lane]=expected_tag(value);
            }
            const Bank original=input;
            sirius_native_tag_predicate8(input.data(),tags.data(),index,mask);
            ++predicate_calls;
            for(unsigned i=0;i<tags.size();++i) { ++predicate_tag_words;if(tags[i]!=expected[i]) { _mm_setcsr(saved);return 6; } }
            for(unsigned i=0;i<input.size();++i) { ++predicate_input_words;if(input[i]!=original[i]) { _mm_setcsr(saved);return 7; } }
            if((_mm_getcsr()&~0x3fU)!=(controls&~0x3fU)) { _mm_setcsr(saved);return 3; }
        }
        _mm_setcsr(saved);
        std::cout<<"{\"passed\":true,\"predicate_calls\":"<<predicate_calls
                 <<",\"compared_tag_words\":"<<predicate_tag_words
                 <<",\"input_words_checked\":"<<predicate_input_words<<"}\n";
        return 0;
    }
    if(argc!=1) { _mm_setcsr(saved);return 8; }

    std::uint64_t calls=0, compared=0, tag_words=0, lane_results=0;
    const std::array<unsigned,20> normalized_cases={0,1,3,4,5,6,7,8,9,10,11,12,33,34,36,37,38,44,45,48};
    for(unsigned tag_pattern=0;tag_pattern<3;++tag_pattern)
    for(unsigned pattern=0;pattern<3;++pattern)
    for(unsigned p=0;p<cases.size();p+=8) for(unsigned q=0;q<cases.size();++q)
    for(unsigned mask:masks) for(unsigned destination=0;destination<3;++destination)
    for(unsigned b_index=0;b_index<2;++b_index) for(unsigned operation=0;operation<5;++operation) {
        Bank actual; actual.fill(0xa5a5a5a5U);
        for(unsigned lane=0;lane<8;++lane) for(unsigned field=0;field<5;++field) {
            const auto ai=pattern==0 ? (p+lane)%cases.size() : pattern==1 ? (p+q)%cases.size() : normalized_cases[(p+lane)%normalized_cases.size()];
            const auto bi=pattern==0 ? (q+3*lane)%cases.size() : pattern==1 ? q : normalized_cases[(q+3*lane)%normalized_cases.size()];
            actual[field*8+lane]=cases[ai][field];
            actual[40+field*8+lane]=cases[bi][field];
        }
        Tags tags;tags.fill(0xa5a5a5a5U);
        for(unsigned index=0;index<2;++index) for(unsigned lane=0;lane<8;++lane) {
            Words value;
            for(unsigned field=0;field<5;++field) value[field]=actual[index*40+field*8+lane];
            tags[index*8+lane]=tag_pattern==0 ? expected_tag(value) : tag_pattern==1 ? 0 : (expected_tag(value)&1U);
        }
        Tags expected_tags=tags;
        Bank expected=actual;
        for(unsigned lane=0;lane<8;++lane) if(mask&(1U<<lane)) {
            Words aw,bw,rw;
            for(unsigned field=0;field<5;++field) {
                aw[field]=actual[field*8+lane];
                bw[field]=actual[b_index*40+field*8+lane];
            }
            g::RetainedTriple_0 a,b;
            std::memcpy(&a,aw.data(),20); std::memcpy(&b,bw.data(),20);
            auto r=operation==0 ? g::RTAdd_0(&a,&b) : operation==1 ? g::RTMultiply_0(&a,&b) :
                   operation==2 ? g::RTDivide_0(&a,&b) : operation==3 ? g::RTSqrt_0(&a) : g::RTSubtract_0(&a,&b);
            std::memcpy(rw.data(),&r,20);
            expected_tags[destination*8+lane]=expected_tag(rw);
            for(unsigned field=0;field<5;++field) expected[destination*40+field*8+lane]=rw[field];
            ++lane_results;
        }
        if(operation==0) sirius_native_tagged_add8(actual.data(),tags.data(),destination,0,b_index,mask);
        else if(operation==1) sirius_native_tagged_multiply8(actual.data(),tags.data(),destination,0,b_index,mask);
        else if(operation==2) sirius_native_tagged_divide8(actual.data(),tags.data(),destination,0,b_index,mask);
        else if(operation==3) sirius_native_tagged_sqrt8(actual.data(),tags.data(),destination,0,b_index,mask);
        else sirius_native_tagged_subtract8(actual.data(),tags.data(),destination,0,b_index,mask);
        ++calls;
        for(unsigned i=0;i<tags.size();++i) {
            ++tag_words;
            if(tags[i]!=expected_tags[i]) {
                _mm_setcsr(saved);
                std::cerr<<"tag mismatch tag_pattern="<<tag_pattern<<" pattern="<<pattern<<" p="<<p<<" q="<<q
                         <<" mask="<<mask<<" destination="<<destination<<" b="<<b_index<<" op="<<operation
                         <<" word="<<i<<" expected="<<expected_tags[i]<<" actual="<<tags[i]<<'\n';
                return 4;
            }
        }
        if((_mm_getcsr()&~0x3fU)!=(controls&~0x3fU)) { _mm_setcsr(saved); return 3; }
        for(unsigned word=0;word<actual.size();++word) {
            ++compared;
            if(actual[word]!=expected[word]) {
                _mm_setcsr(saved);
                std::cerr<<"mismatch pattern="<<pattern<<" p="<<p<<" q="<<q<<" mask="<<mask<<" destination="<<destination
                         <<" b="<<b_index<<" op="<<operation<<" word="<<word
                         <<" expected="<<expected[word]<<" actual="<<actual[word]<<'\n';
                return 1;
            }
        }
    }
    _mm_setcsr(saved);
    std::cout<<"{\"passed\":true,\"calls\":"<<calls<<",\"compared_words\":"<<compared
             <<",\"tag_words\":"<<tag_words<<",\"active_lane_results\":"<<lane_results<<",\"case_tuples\":"<<cases.size()
             <<",\"input_alias_and_inactive_words_checked\":true,\"FP_controls_unchanged\":true}\n";
}
