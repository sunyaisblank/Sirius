// Actual generated writer bodies, complete word comparisons, and a stale-tag
// arithmetic counterexample against the original generated scalar helper.
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <immintrin.h>
#include <transport.native.hpp>
using Words=std::array<std::uint32_t,5>;
using Bank=std::array<std::uint32_t,120>;
using Tags=std::array<std::uint32_t,24>;
using Value=std::array<std::uint32_t,40>;
using Writer=void(*)(std::uint32_t*,std::uint32_t*,std::uint32_t*,std::uint32_t,std::uint32_t,std::uint32_t,std::uint32_t,std::uint32_t);
#define WRITER(stage,version) extern "C" void writer_##stage##_##version(std::uint32_t*,std::uint32_t*,std::uint32_t*,std::uint32_t,std::uint32_t,std::uint32_t,std::uint32_t,std::uint32_t);
WRITER(transport,old) WRITER(transport,fixed)
WRITER(endpoint,old) WRITER(endpoint,fixed)
WRITER(dopri_phase,old) WRITER(dopri_phase,fixed)
extern "C" void sirius_native_tagged_multiply8(std::uint32_t*,std::uint32_t*,std::uint32_t,std::uint32_t,std::uint32_t,std::uint32_t) noexcept;
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
static Words get(const Bank& b,unsigned index,unsigned lane) {
 Words w;for(unsigned f=0;f<5;++f)w[f]=b[index*40+f*8+lane];return w;
}
int main() {
 if(!__builtin_cpu_supports("x86-64-v3"))return 2;
 const auto saved=_mm_getcsr();const auto controls=(saved|0x1f80U)&~0xe040U;_mm_setcsr(controls);
 const auto fail=[&](const char* what){std::cerr<<what<<'\n';_mm_setcsr(saved);return 1;};
#include "cases.inc"
 const std::array<Writer,3> old={writer_transport_old,writer_endpoint_old,writer_dopri_phase_old};
 const std::array<Writer,3> fixed={writer_transport_fixed,writer_endpoint_fixed,writer_dopri_phase_fixed};
 const std::array<unsigned,6> masks={0,1,15,85,128,255};
 const std::array<unsigned,4> tag_seeds={0,1,0xa5a5a5a5U,63};
 std::uint64_t calls=0,bank_words=0,tag_words=0,counterexamples=0,witness_words=0;
 for(unsigned stage=0;stage<3;++stage) for(unsigned p=0;p<cases.size();++p)
 for(unsigned pattern=0;pattern<2;++pattern) for(unsigned index=0;index<3;++index)
 for(unsigned mask:masks) for(unsigned seed:tag_seeds) for(unsigned op=0;op<4;++op)
 for(unsigned component=0;component<(op<2?5U:1U);++component) {
  Bank bank;Tags tags;Value value;
  for(unsigned r=0;r<3;++r)for(unsigned lane=0;lane<8;++lane) {
   const auto& w=cases[(p+11*r+lane)%cases.size()];
   for(unsigned f=0;f<5;++f)bank[r*40+f*8+lane]=w[f];
   tags[r*8+lane]=seed;
  }
  for(unsigned lane=0;lane<8;++lane)for(unsigned f=0;f<5;++f)value[f*8+lane]=cases[(p+(pattern==0?0:lane))%cases.size()][f];
  Bank expected=bank;Tags expected_tags=tags;const Value input=value;
  for(unsigned lane=0;lane<8;++lane)if(mask&(1U<<lane)) {
   const unsigned dest=(op==0||op==2)&&pattern!=0?(index+lane)%3:index;
   if(op<2){expected[dest*40+component*8+lane]=value[component*8+lane];expected_tags[dest*8+lane]=0;}
   else {Words w;for(unsigned f=0;f<5;++f){w[f]=value[f*8+lane];expected[dest*40+f*8+lane]=w[f];}expected_tags[dest*8+lane]=expected_tag(w);}
  }
  fixed[stage](bank.data(),tags.data(),value.data(),op,index,component,mask,pattern);++calls;
  for(unsigned i=0;i<bank.size();++i){++bank_words;if(bank[i]!=expected[i])return fail("fixed writer payload/mask/index mismatch");}
  for(unsigned i=0;i<tags.size();++i){++tag_words;if(tags[i]!=expected_tags[i])return fail("fixed writer tag/inactive mismatch");}
  if(value!=input)return fail("writer changed source words");
  if((_mm_getcsr()&~0x3fU)!=(controls&~0x3fU))return fail("writer changed FP controls");
 }
 // Start from coherent 1*2, then replace a by (1,+2^61,-2^61,0,valid1).
 // The exact original multiply refuses the out-of-domain intermediate low
 // product. A stale normalized tag instead admits the private fast path.
 const Words bad={0x3f800000U,0x5e000000U,0xde000000U,0,1};
 for(unsigned stage=0;stage<3;++stage)for(unsigned op=0;op<4;++op)for(unsigned contract=0;contract<2;++contract) {
  Bank before;Tags initial;Value value;
  for(unsigned r=0;r<3;++r)for(unsigned lane=0;lane<8;++lane) {
   const Words w={r==0?0x3f800000U:r==1?0x40000000U:0x40400000U,0,0,0,1};
   for(unsigned f=0;f<5;++f)before[r*40+f*8+lane]=w[f];initial[r*8+lane]=contract==0 ? (expected_tag(w)&1U) : expected_tag(w);
  }
  for(unsigned f=0;f<5;++f)for(unsigned lane=0;lane<8;++lane)value[f*8+lane]=bad[f];
  std::array<Bank,2> banks={before,before};std::array<Tags,2> tags={initial,initial};
  for(unsigned version=0;version<2;++version) {
   const auto fn=version==0?old[stage]:fixed[stage];
   if(op<2){fn(banks[version].data(),tags[version].data(),value.data(),op,0,1,255,0);fn(banks[version].data(),tags[version].data(),value.data(),op,0,2,255,0);}
   else fn(banks[version].data(),tags[version].data(),value.data(),op,0,0,255,0);
  }
  if(banks[0]!=banks[1])return fail("old/new writer payload differs");
  for(unsigned lane=0;lane<8;++lane)if(tags[0][lane]!=(contract==0 ? 1U : 43U)||tags[1][lane]!=(op<2 ? 0U : expected_tag(bad)))return fail("stale-tag counterexample not observed");
  Bank expected=banks[1];
  for(unsigned lane=0;lane<8;++lane) {
   auto aw=get(expected,0,lane),bw=get(expected,1,lane);sirius_generated_transport::RetainedTriple_0 a,b;
   std::memcpy(&a,aw.data(),20);std::memcpy(&b,bw.data(),20);auto result=sirius_generated_transport::RTMultiply_0(&a,&b);Words rw;std::memcpy(rw.data(),&result,20);
   if(rw[4]!=0)return fail("original scalar did not refuse witness");
   for(unsigned f=0;f<5;++f)expected[80+f*8+lane]=rw[f];
  }
  for(unsigned version=0;version<2;++version)sirius_native_tagged_multiply8(banks[version].data(),tags[version].data(),2,0,1,255);
  if(banks[0]==expected)return fail("old stale tag did not expose arithmetic mismatch");
  for(unsigned i=0;i<expected.size();++i){++witness_words;if(banks[1][i]!=expected[i])return fail("fixed writer arithmetic disagrees with original scalar");}
  ++counterexamples;
 }
 _mm_setcsr(saved);
 std::cout<<"{\"passed\":true,\"fixed_writer_calls\":"<<calls<<",\"compared_bank_words\":"<<bank_words
          <<",\"compared_tag_words\":"<<tag_words<<",\"old_stale_tag_arithmetic_counterexamples\":"<<counterexamples
          <<",\"witness_bank_words\":"<<witness_words<<",\"case_tuples\":"<<cases.size()<<",\"stages\":3,\"writer_forms_per_stage\":4}\n";
}
