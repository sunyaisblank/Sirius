// Fixed private TwoProduct cases: historical source comparison, independent
// integer-dyadic centers and exact signed-zero/exception-state expectations.
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <immintrin.h>

extern "C" void original_product(const std::uint32_t*, const std::uint32_t*,
                                  std::uint32_t*, std::uint32_t*) noexcept;
extern "C" void candidate_product(const std::uint32_t*, const std::uint32_t*,
                                   std::uint32_t*, std::uint32_t*) noexcept;

namespace {
struct Lane { std::uint32_t a, b, hi, lo; };
struct Case { const char* name; bool independent; std::array<Lane, 8> lane; };
constexpr Case Uniform(const char* name, std::uint32_t a, std::uint32_t b,
                       bool independent=false, std::uint32_t hi=0, std::uint32_t lo=0) {
    Case c{name, independent, {}}; c.lane.fill({a,b,hi,lo}); return c;
}
constexpr auto cases = [] {
    std::array<Case, 24> c{{
        Uniform("all_signed_zero_pairs",0,0,true),
        Uniform("zero_noncoarse_low_positive",0,0x3f800001U,true),
        Uniform("zero_noncoarse_low_negative",0,0x3f800fffU,true,0,0x80000000U),
        Uniform("zero_tie_even",0,0x3f800800U,true),
        Uniform("zero_tie_odd",0,0x3f801800U,true,0,0x80000000U),
        Uniform("zero_split_carry",0,0x3ffff800U,true,0,0x80000000U),
        Uniform("coarse_exact",0x3fc00000U,0x3fe00000U,true,0x40280000U),
        Uniform("exact_noncoarse_power",0x3f800fffU,0x40000000U,true,0x40000fffU,0x80000000U),
        Uniform("exact_tie_even_power",0x3f800800U,0x40000000U,true,0x40000800U),
        Uniform("exact_tie_odd_power",0x3f801800U,0x40000000U,true,0x40001800U,0x80000000U),
        Uniform("exact_carry_power",0x3ffff800U,0x40000000U,true,0x407ff800U,0x80000000U),
        Uniform("nonzero_residual_all_signs",0x3f800001U,0x3f800001U,true,0x3f800002U,0x28800000U),
        Uniform("mixed_nonzero_exact_zero",0x3f800001U,0x3f800001U,true,0x3f800002U,0x28800000U),
        Uniform("mixed_all_zero_residuals",0x3fc00000U,0x3fe00000U,true,0x40280000U),
        Uniform("joint_floor",0x2b800001U,0x2b800001U,true,0x17800002U,0x00800000U),
        Uniform("below_joint_floor_fallback",0x2b000001U,0x2b800001U),
        Uniform("individual_floor",0x0c000fffU,0x5d800000U,true,0x2a000fffU,0x80000000U),
        Uniform("below_individual_floor_DE_fallback",0x0b800001U,0x5d800000U),
        Uniform("maximum_safe_operands",0x5d800000U,0x5d800000U,true,0x7b800000U),
        Uniform("above_operand_bound_fallback",0x5d800001U,0x3f800001U),
        Uniform("one_unsafe_lane_fallback",0x3f800001U,0x3f800001U),
        Uniform("one_below_joint_lane_fallback",0x3f800001U,0x3f800001U),
        Uniform("exact_nonpower",0x3fc00000U,0x3f800ffeU,true,0x3fc017fdU,0x80000000U),
        Uniform("coarse_and_zero_mixture",0x3fc00000U,0x3fe00000U,true,0x40280000U),
    }};
    for (unsigned lane=0; lane<8; ++lane) {
        const auto sa=(lane&1U)<<31, sb=((lane>>1)&1U)<<31;
        c[0].lane[lane]={sa,sb,sa^sb,0};
        for (unsigned index=1; index<=5; ++index) {
            auto& x=c[index].lane[lane];
            x.a=sa; x.b^=sb; x.hi=sa^sb; x.lo^=sb;
            if (lane&4U) { const auto t=x.a; x.a=x.b; x.b=t; }
        }
        for (unsigned index=6; index<=11; ++index) {
            auto& x=c[index].lane[lane]; x.a^=sa; x.b^=sb; x.hi^=sa^sb;
            if (index==6) x.lo=0;
            else if (index==11) x.lo^=sa^sb;
            else x.lo^=sa; // split low belongs to A; power's sign is irrelevant.
            if (lane&4U) { const auto t=x.a; x.a=x.b; x.b=t; }
        }
        c[12].lane[lane]=c[lane%4==0 ? 11 : lane%4==1 ? 7 : lane%4==2 ? 2 : 6].lane[lane];
        c[13].lane[lane]=c[lane%3==0 ? 7 : lane%3==1 ? 2 : 6].lane[lane];
        c[23].lane[lane]=(lane&1U) ? c[0].lane[lane] : c[6].lane[lane];
    }
    c[20].lane[7]={0x0b800001U,0x5d800000U,0,0};
    c[21].lane[7]={0x2b000001U,0x2b800001U,0,0};
    return c;
}();

// Exact dyadic verification of the declared centers, independent of floating
// operations and of every candidate bound, predicate and zero-sign formula.
struct Dyadic { std::int64_t integer; int exponent; };
constexpr Dyadic Decode(std::uint32_t bits) {
    const auto biased=(bits>>23)&255U;
    auto n=(bits&0x7fffffU)|(biased ? 0x800000U : 0U);
    int e=biased ? static_cast<int>(biased)-150 : -149;
    while (n && !(n&1U)) { n>>=1; ++e; }
    return {(bits>>31) ? -static_cast<std::int64_t>(n) : n,e};
}
constexpr bool Scale(Dyadic d,int e,std::int64_t& result) {
    if (!d.integer) { result=0; return true; }
    const auto shift=d.exponent-e;
    if (shift<0 || shift>62) return false;
    const auto mag=static_cast<std::uint64_t>(d.integer<0 ? -d.integer : d.integer);
    if (mag>(static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())>>shift)) return false;
    const auto n=static_cast<std::int64_t>(mag<<shift); result=d.integer<0 ? -n : n; return true;
}
constexpr bool Centers() {
    unsigned count=0;
    for (const auto& c:cases) if(c.independent) {
        ++count;
        for (const auto& x:c.lane) {
            const auto a=Decode(x.a),b=Decode(x.b),h=Decode(x.hi),l=Decode(x.lo);
            const Dyadic p{a.integer*b.integer,a.exponent+b.exponent};
            int e=p.exponent; if(h.integer && h.exponent<e)e=h.exponent; if(l.integer && l.exponent<e)e=l.exponent;
            std::int64_t pp=0,hh=0,ll=0;
            if(!Scale(p,e,pp)||!Scale(h,e,hh)||!Scale(l,e,ll)||pp!=hh+ll)return false;
        }
    }
    return count==19;
}
static_assert(Centers());
}

int main() {
    if (!__builtin_cpu_supports("x86-64-v4")) return 77;
    const unsigned saved=_mm_getcsr();
    const unsigned scope=(saved|0x1f80U)&~(0xe040U|0x3fU);
    bool all=true,denormal_witness=false;
    std::cout<<std::boolalpha<<"{\"cases\":[";
    for(unsigned index=0;index<cases.size();++index) {
        const auto& c=cases[index];
        std::array<std::uint32_t,8> a{},b{};
        for(unsigned i=0;i<8;++i){a[i]=c.lane[i].a;b[i]=c.lane[i].b;}
        const auto aa=a,bb=b;
        std::array<std::uint32_t,12> oh,ol,nh,nl;
        oh.fill(0xa5a5a5a5U);ol=nh=nl=oh;
        _mm_setcsr(scope);original_product(a.data(),b.data(),oh.data()+2,ol.data()+2);const auto old_fp=_mm_getcsr();
        _mm_setcsr(scope);candidate_product(a.data(),b.data(),nh.data()+2,nl.data()+2);const auto new_fp=_mm_getcsr();
        _mm_setcsr(saved);
        const bool restored=_mm_getcsr()==saved;
        const bool equal=oh==nh && ol==nl;
        const bool immutable=a==aa && b==bb;
        bool canary=true,center=true;
        for(unsigned i: {0U,1U,10U,11U})canary &= oh[i]==0xa5a5a5a5U && ol[i]==0xa5a5a5a5U && nh[i]==0xa5a5a5a5U && nl[i]==0xa5a5a5a5U;
        for(unsigned i=0;i<8;++i)if(c.independent)center &= nh[i+2]==c.lane[i].hi && nl[i+2]==c.lane[i].lo;
        const bool controls=(old_fp&~0x3fU)==scope && (new_fp&~0x3fU)==scope;
        const bool flags=(old_fp&0x3fU)==(new_fp&0x3fU);
        if(index==17)denormal_witness=(old_fp&2U)!=0;
        const bool passed=equal&&immutable&&canary&&center&&controls&&flags&&restored;
        all &= passed;
        if(index)std::cout<<',';
        std::cout<<"{\"name\":\""<<c.name<<"\",\"passed\":"<<passed<<",\"independent\":"<<c.independent
                 <<",\"whole_words_equal\":"<<equal<<",\"center_and_zero_signs_equal\":"<<center
                 <<",\"old_MXCSR\":"<<old_fp<<",\"new_MXCSR\":"<<new_fp<<",\"caller_restored\":"<<restored<<'}';
    }
    all &= denormal_witness;
    std::cout<<"],\"passed\":"<<all<<",\"case_count\":24,\"candidate_calls\":24,\"historical_primitive_calls\":24,"
             <<"\"compared_hi_low_words\":384,\"compared_canary_words\":384,\"independent_center_cases\":19,"
             <<"\"denormal_reference_witness\":"<<denormal_witness<<",\"scope\":\"Fixed direct primitive words/zero-signs/MXCSR compatibility only; no speed/provider/frame admission\"}\n";
    return all ? 0 : 1;
}
