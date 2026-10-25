// C++26 AVX2 retained arithmetic with unused-lane operands sanitized.
// Derived from the checked AC helper; a new candidate awaiting qualification.
// The caller owns canonical scratch indices and supplies ISPC lanemask().
// All scratch words are initialized; load both operands before any aliasing write.
// Compile with the existing strict ordered FP flags and FMA disabled, under the
// existing x86-64-v3 guard and saved/masked/restored IEEE32 worker scope.
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <immintrin.h>
#include <limits>

namespace {
using U = std::uint32_t __attribute__((vector_size(32)));
using F = float __attribute__((vector_size(32)));
static_assert(sizeof(U) == 32 && sizeof(F) == 32);
static_assert(sizeof(float) == 4 && std::numeric_limits<float>::is_iec559 &&
              std::numeric_limits<float>::digits == 24);
struct Triple { F hi, lo, tail, error; U valid; };
static_assert(sizeof(Triple) == 160);

__attribute__((always_inline)) inline U u(std::uint32_t x) { U result={x,x,x,x,x,x,x,x}; return result; }
__attribute__((always_inline)) inline U bits(F x) { return std::bit_cast<U>(x); }
__attribute__((always_inline)) inline F from(U x) { return std::bit_cast<F>(x); }
__attribute__((always_inline)) inline F f(std::uint32_t x) { return from(u(x)); }
__attribute__((always_inline)) inline F zero() { return f(0); }
__attribute__((always_inline)) inline U eq(F a, F b) {
    U aa=bits(a), bb=bits(b), am=aa&u(2147483647U), bm=bb&u(2147483647U);
    U same=std::bit_cast<U>(aa==bb) |
           (std::bit_cast<U>(am==u(0)) & std::bit_cast<U>(bm==u(0)));
    return same & std::bit_cast<U>(am<=u(2139095040U)) &
                  std::bit_cast<U>(bm<=u(2139095040U));
}
__attribute__((always_inline)) inline U nz(U x) { return std::bit_cast<U>(x != u(0)); }
__attribute__((always_inline)) inline U bounded(F x, std::uint32_t limit) {
    return std::bit_cast<U>((bits(x) & u(2147483647U)) <= u(limit));
}
__attribute__((always_inline)) inline F choose(U m, F a, F b) { return from((m & bits(a)) | (~m & bits(b))); }
__attribute__((always_inline)) inline U choose(U m, U a, U b) { return (m & a) | (~m & b); }
__attribute__((always_inline)) inline Triple choose(U m, Triple a, Triple b) {
    return {choose(m,a.hi,b.hi), choose(m,a.lo,b.lo),
            choose(m,a.tail,b.tail), choose(m,a.error,b.error), choose(m,a.valid,b.valid)};
}
__attribute__((always_inline)) inline bool any(U m) { return _mm256_movemask_ps(std::bit_cast<__m256>(m)) != 0; }
__attribute__((always_inline)) inline bool all_on(U m, U active) { return !any(active & ~m); }
__attribute__((always_inline)) inline F add(F a, F b) { F result = a + b; return result; }
__attribute__((always_inline)) inline F sub(F a, F b) { F result = a - b; return result; }
__attribute__((always_inline)) inline F mul(F a, F b) { F result = a * b; return result; }
__attribute__((always_inline)) inline F abs32(F a) { return from(bits(a) & u(2147483647U)); }
__attribute__((always_inline)) inline Triple invalid() { return {zero(),zero(),zero(),zero(),u(0)}; }
__attribute__((always_inline)) inline Triple constant_zero() { return {zero(),zero(),zero(),zero(),u(1)}; }
__attribute__((always_inline)) inline Triple negate(Triple a) {
    a.hi=from(bits(a.hi)^u(2147483648U));
    a.lo=from(bits(a.lo)^u(2147483648U));
    a.tail=from(bits(a.tail)^u(2147483648U));
    return a;
}
__attribute__((always_inline)) inline U exact(Triple a, F x) {
    return nz(a.valid) & eq(a.hi,x) & eq(a.lo,zero()) &
           eq(a.tail,zero()) & eq(a.error,zero());
}
__attribute__((always_inline)) inline void sum32(F a, F b, F& hi, F& lo) {
    U swap=std::bit_cast<U>((bits(a)&u(2147483647U)) < (bits(b)&u(2147483647U)));
    F first=choose(swap,b,a), second=choose(swap,a,b);
    F rounded=add(first,second);
    hi=rounded;
    lo=sub(second,sub(rounded,first));
}
__attribute__((always_inline)) inline F up_add(F a, F b) {
    U good=bounded(a,2105540608U)&bounded(b,2105540608U);
    U work=good & ~eq(a,zero()) & ~eq(b,zero());
    F rounded=from(bits(add(choose(work,a,zero()),choose(work,b,zero())))+u(1));
    F nonzero_a=choose(eq(b,zero()),a,rounded);
    return choose(good,choose(eq(a,zero()),b,nonzero_a),f(2130706432U));
}
__attribute__((always_inline)) inline F up_mul(F a, F b) {
    U is_zero=eq(a,zero())|eq(b,zero());
    U good=bounded(a,1585446912U)&bounded(b,1585446912U);
    U work=good & ~is_zero;
    F rounded=from(bits(mul(choose(work,a,zero()),choose(work,b,zero())))+u(1));
    return choose(is_zero,zero(),choose(good,rounded,f(2130706432U)));
}
__attribute__((always_inline)) inline void split32(F a, F& hi, F& lo) {
    F c=mul(f(1166018560U),a);
    F first=sub(c,sub(c,a));
    hi=first; lo=sub(a,first);
}
__attribute__((always_inline)) inline void product32(F a, F b, F& hi, F& lo) {
    F ah,al,bh,bl;
    split32(a,ah,al); split32(b,bh,bl);
    F first=mul(a,b);
    hi=first;
    lo=sub(mul(al,bl),sub(sub(sub(first,mul(ah,bh)),mul(al,bh)),mul(ah,bl)));
}
__attribute__((always_inline)) inline Triple finish(F hi, F lo, F tail, F error) {
    U input_good=bounded(error,2071986176U)&~nz(bits(error)&u(2147483648U));
    F a,b,c,d,e,gap,g;
    sum32(lo,tail,a,b); sum32(hi,a,c,d); sum32(d,b,e,gap);
    Triple r;
    sum32(c,e,r.hi,g); sum32(g,gap,r.lo,r.tail);
    r.error=error; r.valid=u(1);
    return choose(input_good&bounded(r.hi,2071986176U),r,invalid());
}
__attribute__((always_inline)) inline Triple term(Triple a, F x, U active) {
    U is_zero=std::bit_cast<U>((bits(x)&u(2147483647U))==u(0));
    U zero_good=nz(a.valid)&bounded(a.hi,2105540608U)&bounded(a.lo,2105540608U)&
                bounded(a.tail,2105540608U)&bounded(a.error,2105540608U);
    U positive_zero=std::bit_cast<U>(bits(x)==u(0));
    Triple z{
        choose(positive_zero&eq(a.hi,zero()),zero(),a.hi),
        choose(positive_zero&eq(a.lo,zero()),zero(),a.lo),
        choose(positive_zero&eq(a.tail,zero()),zero(),a.tail),
        choose(eq(a.error,zero()),zero(),a.error),u(1)};
    if(all_on(is_zero&zero_good,active)) return z;
    U work=active & ~is_zero & nz(a.valid) & bounded(x,2071986176U);
    a=choose(work,a,invalid()); x=choose(work,x,zero());
    F hi,first,lo,second,tail,last;
    sum32(a.hi,x,hi,first); sum32(a.lo,first,lo,second);
    sum32(a.tail,second,tail,last);
    F error=up_add(a.error,abs32(last));
    U general_good=work&bounded(hi,2105540608U)&
                   bounded(lo,2105540608U)&bounded(tail,2105540608U)&bounded(error,2105540608U);
    U good=(is_zero&zero_good)|(~is_zero&general_good);
    Triple general{hi,lo,tail,error,u(1)};
    return choose(good,choose(is_zero,z,general),invalid());
}
__attribute__((always_inline)) inline Triple product(Triple sum, F a, F b, U active) {
    U good=nz(sum.valid)&bounded(a,1568669696U)&bounded(b,1568669696U);
    U is_zero=eq(a,zero())|eq(b,zero());
    if(all_on(good&is_zero,active)) return sum;
    U work=active & good & ~is_zero;
    a=choose(work,a,zero()); b=choose(work,b,zero());
    F hi,lo; product32(a,b,hi,lo);
    Triple r=term(sum,hi,work); r=term(r,lo,work);
    r.error=choose(nz(r.valid),up_add(r.error,f(5)),r.error);
    return choose(good,choose(is_zero,sum,r),invalid());
}
__attribute__((always_inline)) inline Triple rounded_product(Triple sum, F a, F b, U active) {
    U good=nz(sum.valid)&bounded(a,1568669696U)&bounded(b,1568669696U);
    U is_zero=eq(a,zero())|eq(b,zero());
    if(all_on(good&is_zero,active)) return sum;
    U work=active & good & ~is_zero;
    a=choose(work,a,zero()); b=choose(work,b,zero());
    F p=mul(a,b);
    F error=up_add(up_mul(abs32(p),f(872415232U)),f(5));
    Triple r=term(sum,p,work);
    r.error=choose(nz(r.valid),up_add(r.error,error),r.error);
    return choose(good,choose(is_zero,sum,r),invalid());
}
__attribute__((always_inline)) inline F center_upper(Triple a) { return up_add(up_add(abs32(a.hi),abs32(a.lo)),abs32(a.tail)); }
__attribute__((always_inline)) inline Triple add_triple(Triple a, Triple b, U active) {
    U good=nz(a.valid)&nz(b.valid), az=exact(a,zero()), bz=exact(b,zero());
    U general_active=active&good&~az&~bz;
    Triple general=invalid();
    if(any(general_active)) {
        Triple r=term(a,b.hi,general_active);
        r=term(r,b.lo,general_active); r=term(r,b.tail,general_active);
        general=choose(nz(r.valid),finish(r.hi,r.lo,r.tail,up_add(r.error,b.error)),r);
    }
    return choose(good,choose(az,b,choose(bz,a,general)),invalid());
}
__attribute__((always_inline)) inline Triple multiply_triple(Triple a, Triple b, U active) {
    U good=nz(a.valid)&nz(b.valid);
    U is_zero=exact(a,zero())|exact(b,zero());
    U a_one=exact(a,f(1065353216U)), b_one=exact(b,f(1065353216U));
    U a_minus=exact(a,f(3212836864U)), b_minus=exact(b,f(3212836864U));
    U general_active=active&good&~is_zero&~a_one&~b_one&~a_minus&~b_minus;
    Triple general=invalid();
    if(any(general_active)) {
        Triple r=product(constant_zero(),a.hi,b.hi,general_active);
        r=product(r,a.hi,b.lo,general_active); r=product(r,a.lo,b.hi,general_active);
        r=rounded_product(r,a.lo,b.lo,general_active);
        r=rounded_product(r,a.hi,b.tail,general_active);
        r=rounded_product(r,a.tail,b.hi,general_active);
        F omitted=up_add(up_add(up_mul(abs32(a.lo),abs32(b.tail)),
                                  up_mul(abs32(a.tail),abs32(b.lo))),
                         up_mul(abs32(a.tail),abs32(b.tail)));
        F a_error=a.error, b_center=center_upper(b);
        F first_error=up_mul(a_error,b_center);
        F b_error=b.error, a_center=center_upper(a);
        F error=up_add(up_add(r.error,omitted),
                       up_add(up_add(first_error,up_mul(b_error,a_center)),up_mul(a_error,b_error)));
        general=choose(nz(r.valid),finish(r.hi,r.lo,r.tail,error),r);
    }
    Triple selected=choose(a_minus,negate(b),choose(b_minus,negate(a),general));
    selected=choose(a_one,b,choose(b_one,a,selected));
    selected=choose(is_zero,constant_zero(),selected);
    return choose(good,selected,invalid());
}
__attribute__((always_inline)) inline U load_words(const std::uint32_t* words) {
    U result; std::memcpy(&result,words,sizeof(result)); return result;
}
__attribute__((always_inline)) inline Triple read(const std::uint32_t* scratch,std::uint32_t index) {
    const auto* p=scratch+std::size_t(index)*40;
    return {from(load_words(p)),from(load_words(p+8)),from(load_words(p+16)),
            from(load_words(p+24)),load_words(p+32)};
}
__attribute__((always_inline)) inline void write(std::uint32_t* scratch,std::uint32_t index,Triple value,U active) {
    auto* p=scratch+std::size_t(index)*40;
    U columns[5]={bits(value.hi),bits(value.lo),bits(value.tail),bits(value.error),value.valid};
    for(unsigned c=0;c<5;++c) {
        U merged=choose(active,columns[c],load_words(p+c*8));
        std::memcpy(p+c*8,&merged,sizeof(merged));
    }
}
__attribute__((always_inline)) inline U active_mask(std::uint32_t mask) {
    U lanes={1,2,4,8,16,32,64,128};
    return nz(lanes&u(mask));
}
} // namespace

extern "C" __attribute__((visibility("hidden")))
void sirius_native_dense_add8(std::uint32_t* scratch,std::uint32_t destination,
                             std::uint32_t a_index,std::uint32_t b_index,
                             std::uint32_t mask) noexcept {
    Triple a=read(scratch,a_index), b=read(scratch,b_index);
    U active=active_mask(mask);
    write(scratch,destination,add_triple(a,b,active),active);
}

extern "C" __attribute__((visibility("hidden")))
void sirius_native_dense_multiply8(std::uint32_t* scratch,std::uint32_t destination,
                                  std::uint32_t a_index,std::uint32_t b_index,
                                  std::uint32_t mask) noexcept {
    Triple a=read(scratch,a_index), b=read(scratch,b_index);
    U active=active_mask(mask);
    write(scratch,destination,multiply_triple(a,b,active),active);
}
