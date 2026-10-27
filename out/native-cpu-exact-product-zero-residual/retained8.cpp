// C++26 AVX2 arithmetic closure; original division/sqrt refinement order retained.
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
#include <initializer_list>

namespace {
using U = std::uint32_t __attribute__((vector_size(32)));
using F = float __attribute__((vector_size(32)));
using I = std::int32_t __attribute__((vector_size(32)));
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
    // RNE, DAZ/FTZ off, gradual underflow, masked exceptions, no contraction.
    // Whole-eight domain: two normals with e>=-103, |x|<=2^60, eA+eB>=-80;
    // or signed zero with such a normal; or two signed zeros. Split lows lie
    // on2^(e-23)>=2^-126; nonzero cross-products on2^(eA+eB-46)>=2^-126.
    // Thus the original TwoProduct and the binary64 residual are exact here.
    // For r=low12 bits, split low is+0 if r=0; otherwise its sign is sign(x)
    // XOR(r>2048 || r==2048 && retained bit12 odd), including exponent carry.
    // Exact residual zero has the sole nonzero split-low sign, or+0 otherwise;
    // a zero operand's sign never enters that low sign. Keep original high mul.
    // Old flags are only IX, iff some operand is noncoarse. A nonzero residual
    // makes high mul supply IX; otherwise retain one original4097*x multiply,
    // consumed by the finite guard below. Its emitted dependency must be checked.
    using D=double __attribute__((vector_size(64)));
    static_assert(sizeof(D)==64 && sizeof(double)==8 &&
                  std::numeric_limits<double>::is_iec559 &&
                  std::numeric_limits<double>::digits==53);
    U am=bits(a)&u(2147483647U),bm=bits(b)&u(2147483647U);
    U ae=am>>23,be=bm>>23;
    U az=~nz(am),bz=~nz(bm);
    U an=std::bit_cast<U>(ae>=u(24))&std::bit_cast<U>(am<=u(1568669696U));
    U bn=std::bit_cast<U>(be>=u(24))&std::bit_cast<U>(bm<=u(1568669696U));
    U domain=(az|an)&(bz|bn)&(az|bz|std::bit_cast<U>((ae+be)>=u(174)));
    if(all_on(domain,u(0xffffffffU))) {
        F first=mul(a,b);
        D wide_a=__builtin_convertvector(a,D),wide_b=__builtin_convertvector(b,D);
        D exact_product=wide_a*wide_b;
        D wide_high=__builtin_convertvector(first,D);
        D exact_residual=exact_product-wide_high;
        F residual=__builtin_convertvector(exact_residual,F);
        U nonzero=nz(bits(residual)&u(2147483647U));
        U ar=am&u(4095),br=bm&u(4095),al_nonzero=nz(ar),bl_nonzero=nz(br);
        bool ix_kept=true;
        if(!any(nonzero) && any(al_nonzero|bl_nonzero)) {
            F chosen=choose(al_nonzero,a,choose(bl_nonzero,b,zero()));
            F retained_split=mul(f(1166018560U),chosen);
            //4097*|chosen|<=4097*2^60<2^73; retain the actual value dependency.
            ix_kept=all_on(std::bit_cast<U>((bits(retained_split)&u(2147483647U))<u(0x64000000U)),
                           u(0xffffffffU));
        }
        if(ix_kept) {
            U a_up=std::bit_cast<U>(ar>u(2048))|
                   (std::bit_cast<U>(ar==u(2048))&nz(am&u(4096)));
            U b_up=std::bit_cast<U>(br>u(2048))|
                   (std::bit_cast<U>(br==u(2048))&nz(bm&u(4096)));
            U a_low_sign=(bits(a)^(a_up&u(2147483648U)))&al_nonzero&u(2147483648U);
            U b_low_sign=(bits(b)^(b_up&u(2147483648U)))&bl_nonzero&u(2147483648U);
            U zero_sign=(a_low_sign|b_low_sign)&(al_nonzero^bl_nonzero);
            hi=first;lo=choose(nonzero,residual,from(zero_sign));return;
        }
    }
    // Every other vector executes the literal original17-operation body below.
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
// Fast-domain invariant: valid1, finite canonical-sign radius <=2^120,
// |lo|<=ulp(hi), |tail|<=ulp(lo), and zero hi has zero lower limbs. The
// integer spacing words equal nextafter(abs(x),+inf)-abs(x) for finite x.
// Every other input retains the complete original implementation below.
__attribute__((always_inline)) inline U spacing_word(F a) {
    U exponent=(bits(a)&u(2147483647U))>>23;
    U normal=std::bit_cast<U>(exponent>=u(24));
    U shift=choose(nz(exponent),exponent-u(1),u(0))&u(31);
    return choose(normal,(exponent-u(23))<<23,u(1)<<shift);
}
__attribute__((always_inline)) inline U normalized(Triple a) {
    U h=bits(a.hi)&u(2147483647U),l=bits(a.lo)&u(2147483647U),
      t=bits(a.tail)&u(2147483647U),e=bits(a.error);
    U result=std::bit_cast<U>(a.valid==u(1))&std::bit_cast<U>(h<=u(2071986176U))&
        std::bit_cast<U>(e<=u(2071986176U))&
        std::bit_cast<U>(l<=spacing_word(a.hi))&std::bit_cast<U>(t<=spacing_word(a.lo));
    return result & (nz(h)|~nz(l|t));
}
__attribute__((always_inline)) inline U fast_eq(F a,F b) {return std::bit_cast<U>(a==b);}
__attribute__((always_inline)) inline U fast_exact(Triple a,F x) {
    return fast_eq(a.hi,x)&fast_eq(a.lo,zero())&fast_eq(a.tail,zero())&fast_eq(a.error,zero());
}
// These additions have nonnegative inputs <=2^124 under the private-domain
// bounds below. Keep the original zero precedence and next-word rounding.
__attribute__((always_inline)) inline F fast_up_add(F a,F b) {
    F rounded=from(bits(add(a,b))+u(1));
    return choose(fast_eq(a,zero()),b,choose(fast_eq(b,zero()),a,rounded));
}
__attribute__((always_inline)) inline Triple fast_term(Triple a,F x,U active) {
    U is_zero=fast_eq(x,zero());
    const auto zero_result=[&]() {
        U positive_zero=std::bit_cast<U>(bits(x)==u(0));
        return Triple{choose(positive_zero&fast_eq(a.hi,zero()),zero(),a.hi),
            choose(positive_zero&fast_eq(a.lo,zero()),zero(),a.lo),
            choose(positive_zero&fast_eq(a.tail,zero()),zero(),a.tail),
            choose(fast_eq(a.error,zero()),zero(),a.error),u(1)};
    };
    if(all_on(is_zero,active)) return zero_result();
    F hi,first,lo,second,tail,last;
    sum32(a.hi,x,hi,first);sum32(a.lo,first,lo,second);sum32(a.tail,second,tail,last);
    Triple result{hi,lo,tail,fast_up_add(a.error,abs32(last)),u(1)};
    if(any(active&is_zero)) result=choose(is_zero,zero_result(),result);
    return result;
}
__attribute__((always_inline)) inline Triple fast_product(Triple sum,F a,F b,U active) {
    U is_zero=fast_eq(a,zero())|fast_eq(b,zero());
    if(all_on(is_zero,active)) return sum;
    F hi,lo;product32(a,b,hi,lo);
    Triple r=fast_term(sum,hi,active);r=fast_term(r,lo,active);
    r.error=fast_up_add(r.error,f(5));
    return choose(is_zero,sum,r);
}
__attribute__((always_inline)) inline Triple fast_rounded_product(Triple sum,F a,F b,U active) {
    U is_zero=fast_eq(a,zero())|fast_eq(b,zero());
    if(all_on(is_zero,active)) return sum;
    F p=mul(a,b);
    F error=up_add(up_mul(abs32(p),f(872415232U)),f(5));
    Triple r=fast_term(sum,p,active);r.error=fast_up_add(r.error,error);
    return choose(is_zero,sum,r);
}
// Add: three scalar terms bounded by2^120 keep private limbs/errors below
// 2^124. Multiply additionally bounds |hi|<=2^54: then |lo|<=2^31 and
// |tail|<=2^8, every rounded low-order product <=2^62, all private centers
// <2^112, and accumulated rounding radius <2^44. Thus the omitted private
// validity/bound tests cannot refuse active lanes. Original final error
// propagation/normalization/public2^120 refusal and identity order remain.
__attribute__((always_inline)) inline Triple fast_add_triple(Triple a,Triple b,U active) {
    U az=fast_exact(a,zero()),bz=fast_exact(b,zero());
    U work=active&~az&~bz;Triple general=invalid();
    if(any(work)) {
        Triple r=fast_term(a,b.hi,work);r=fast_term(r,b.lo,work);r=fast_term(r,b.tail,work);
        general=finish(r.hi,r.lo,r.tail,up_add(r.error,b.error));
    }
    return choose(az,b,choose(bz,a,general));
}
__attribute__((always_inline)) inline Triple fast_multiply_triple(Triple a,Triple b,U active) {
    U is_zero=fast_exact(a,zero())|fast_exact(b,zero());
    U a_one=fast_exact(a,f(1065353216U)),b_one=fast_exact(b,f(1065353216U));
    U a_minus=fast_exact(a,f(3212836864U)),b_minus=fast_exact(b,f(3212836864U));
    U work=active&~is_zero&~a_one&~b_one&~a_minus&~b_minus;Triple general=invalid();
    if(any(work)) {
        Triple r=fast_product(constant_zero(),a.hi,b.hi,work);
        r=fast_product(r,a.hi,b.lo,work);r=fast_product(r,a.lo,b.hi,work);
        r=fast_rounded_product(r,a.lo,b.lo,work);
        r=fast_rounded_product(r,a.hi,b.tail,work);
        r=fast_rounded_product(r,a.tail,b.hi,work);
        F omitted=up_add(up_add(up_mul(abs32(a.lo),abs32(b.tail)),up_mul(abs32(a.tail),abs32(b.lo))),
                         up_mul(abs32(a.tail),abs32(b.tail)));
        F a_error=a.error,b_center=center_upper(b),first_error=up_mul(a_error,b_center);
        F b_error=b.error,a_center=center_upper(a);
        F error=up_add(up_add(r.error,omitted),up_add(up_add(first_error,up_mul(b_error,a_center)),up_mul(a_error,b_error)));
        general=finish(r.hi,r.lo,r.tail,error);
    }
    Triple selected=choose(a_minus,negate(b),choose(b_minus,negate(a),general));
    selected=choose(a_one,b,choose(b_one,a,selected));
    return choose(is_zero,constant_zero(),selected);
}
__attribute__((noinline)) Triple slow_add_triple(Triple a, Triple b, U active) {
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
__attribute__((noinline)) Triple slow_multiply_triple(Triple a, Triple b, U active) {
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
__attribute__((always_inline)) inline Triple add_triple(Triple a,Triple b,U active) {
    if(all_on(normalized(a)&normalized(b),active)) return fast_add_triple(a,b,active);
    return slow_add_triple(a,b,active);
}
__attribute__((always_inline)) inline Triple multiply_triple(Triple a,Triple b,U active) {
    if(all_on(normalized(a)&normalized(b)&bounded(a.hi,1518338048U)&bounded(b.hi,1518338048U),active))
        return fast_multiply_triple(a,b,active);
    return slow_multiply_triple(a,b,active);
}
// Integer exponent reconstruction and normal-mantissa seeds follow the original
// retained_pair.slang. Every vector shift count is bounded even in unused lanes.
__attribute__((always_inline)) inline I signed_words(U a) {return std::bit_cast<I>(a);}
__attribute__((always_inline)) inline U unsigned_words(I a) {return std::bit_cast<U>(a);}
__attribute__((always_inline)) inline U less_equal(F a,F b) {return std::bit_cast<U>(a<=b);}
__attribute__((always_inline)) inline U greater(F a,F b) {return std::bit_cast<U>(a>b);}
__attribute__((always_inline)) inline U highest_bit(U a) {
    U result=u(0);
    for(unsigned shift: {16U,8U,4U,2U,1U}) {
        U take=std::bit_cast<U>(a>=u(1U<<shift));
        U amount=take&u(shift);
        a=a>>amount;result=result+amount;
    }
    return result;
}
__attribute__((always_inline)) inline I floor_exponent(F a) {
    U word=bits(a),exponent=word>>23;
    if(all_on(nz(exponent),u(0xffffffffU))) return signed_words(exponent-u(127));
    return signed_words(choose(nz(exponent),exponent-u(127),highest_bit(word)-u(149)));
}
__attribute__((always_inline)) inline F power_two(I e) {
    U low=std::bit_cast<U>(e<signed_words(u(std::uint32_t(-126))));
    U shift=choose(low,unsigned_words(e)+u(149),u(0))&u(31);
    U word=choose(low,u(1)<<shift,(unsigned_words(e)+u(127))<<23);
    word=choose(std::bit_cast<U>(e<signed_words(u(std::uint32_t(-149)))),u(1),word);
    return from(choose(std::bit_cast<U>(e>signed_words(u(120))),u(2130706432U),word));
}
__attribute__((always_inline)) inline F ratio_upper(F n,F d) {
    U good=bounded(n,2071986176U)&bounded(d,2071986176U)&~less_equal(d,zero());
    U nonzero=~eq(n,zero());
    F nn=choose(good&nonzero,n,f(1065353216U));
    F dd=choose(good&nonzero,d,f(1065353216U));
    F ratio=power_two(floor_exponent(nn)-floor_exponent(dd)+signed_words(u(1)));
    return choose(good,choose(nonzero,ratio,zero()),f(2130706432U));
}
__attribute__((always_inline)) inline F mantissa(F a) {
    U word=bits(a),exponent=word>>23;
    if(all_on(nz(exponent),u(0xffffffffU))) return from(u(1065353216U)|(word&u(8388607U)));
    U shift=(u(23)-highest_bit(word))&u(31);
    U fraction=choose(nz(exponent),word,(word<<shift))&u(8388607U);
    return from(u(1065353216U)|fraction);
}
__attribute__((always_inline)) inline F divide_seed(F a,F b,U active) {
    U work=active&~eq(a,zero());
    F aa=choose(work,abs32(a),f(1065353216U));
    F bb=choose(work,abs32(b),f(1065353216U));
    I exponent=floor_exponent(aa)-floor_exponent(bb);
    work=work & ~std::bit_cast<U>(exponent<signed_words(u(std::uint32_t(-151)))) &
                ~std::bit_cast<U>(exponent>signed_words(u(120)));
    F am=mantissa(aa),bm=mantissa(bb);
    F seed=am/bm;
    work=work&bounded(seed,1082130432U);
    U small=std::bit_cast<U>(exponent<signed_words(u(std::uint32_t(-126))));
    // The first ordered scaling is evaluated only for the original small case.
    F scaled=mul(choose(small,seed,zero()),f(8388608U));
    seed=choose(small,scaled,seed);
    exponent=exponent+signed_words(small&u(126));
    seed=mul(seed,power_two(exponent));
    F signed_seed=from(bits(seed)^((bits(a)^bits(b))&u(2147483648U)));
    return choose(work,signed_seed,zero());
}
__attribute__((always_inline)) inline F sqrt_seed(F a,U active) {
    F input=choose(active,a,f(1065353216U));
    I exponent=floor_exponent(input);U odd=unsigned_words(exponent)&u(1);
    F m=mantissa(input);
    m=choose(nz(odd),mul(m,f(1073741824U)),m);
    F seed=std::bit_cast<F>(_mm256_sqrt_ps(std::bit_cast<__m256>(m)));
    U good=bounded(seed,1082130432U);
    F result=mul(seed,power_two((exponent-signed_words(odd))/signed_words(u(2))));
    return choose(good,result,zero());
}
__attribute__((always_inline)) inline Triple constant(F a) {
    Triple r{choose(eq(a,zero()),zero(),a),zero(),zero(),zero(),u(1)};
    return choose(bounded(a,2071986176U),r,invalid());
}
__attribute__((always_inline)) inline Triple center(Triple a) {a.error=zero();return a;}
__attribute__((always_inline)) inline Triple pack(F hi,F lo,F tail,F error) {
    return choose(bounded(hi,2071986176U)&bounded(lo,2071986176U)&bounded(tail,2071986176U),
                  finish(hi,lo,tail,error),invalid());
}
__attribute__((always_inline)) inline F absolute_upper(Triple a) {return up_add(center_upper(a),a.error);}
__attribute__((always_inline)) inline F absolute_lower(Triple a) {
    F tail=up_add(up_add(abs32(a.lo),abs32(a.tail)),a.error);
    U work=nz(a.valid)&~less_equal(abs32(a.hi),tail);
    F difference,correction;
    sum32(abs32(a.hi),from(bits(tail)^u(2147483648U)),difference,correction);
    U word=bits(difference);
    F lower=choose(less_equal(zero(),correction),difference,
                   choose(nz(word),from(word-u(1)),zero()));
    return choose(work,lower,zero());
}
__attribute__((always_inline)) inline Triple subtract_triple(Triple a,Triple b,U active) {
    return add_triple(a,negate(b),active);
}
__attribute__((always_inline)) inline Triple divide_triple(Triple a,Triple b,U active) {
    U eligible=active&nz(a.valid)&nz(b.valid);
    F lower=absolute_lower(b);
    eligible=eligible&~less_equal(lower,zero())&
              ~greater(ratio_upper(absolute_upper(a),lower),f(1577058304U));
    U identity=eligible&exact(b,f(1065353216U));
    U work=eligible&~identity;
    if(!any(work)) return choose(identity,a,invalid());
    Triple q=constant(divide_seed(a.hi,b.hi,work));
    for(unsigned i=0;i<3 && any(work);++i) {
        Triple residual=subtract_triple(a,multiply_triple(b,q,work),work);
        work=work&nz(residual.valid)&
             ~greater(ratio_upper(abs32(residual.hi),abs32(b.hi)),f(1560281088U));
        q=center(add_triple(q,constant(divide_seed(residual.hi,b.hi,work)),work));
        work=work&nz(q.valid);
    }
    if(!any(work)) return choose(identity,a,invalid());
    Triple residual=subtract_triple(a,multiply_triple(b,q,work),work);
    work=work&nz(residual.valid);
    if(!any(work)) return choose(identity,a,invalid());
    Triple result=pack(q.hi,q.lo,q.tail,ratio_upper(absolute_upper(residual),lower));
    return choose(identity,a,choose(work,result,invalid()));
}
__attribute__((always_inline)) inline Triple sqrt_triple(Triple a,U active) {
    U input_invalid=~nz(a.valid);
    U exact_zero=eq(a.hi,zero())&eq(a.lo,zero())&eq(a.tail,zero())&eq(a.error,zero());
    U work=active&~input_invalid&~exact_zero&~less_equal(a.hi,zero())&
           ~less_equal(absolute_lower(a),zero());
    if(!any(work)) return choose(input_invalid,a,choose(exact_zero,constant_zero(),invalid()));
    Triple q=constant(sqrt_seed(a.hi,work));
    for(unsigned i=0;i<3 && any(work);++i) {
        work=work&nz(q.valid)&~less_equal(q.hi,zero())&bounded(q.hi,1560281088U);
        Triple residual=subtract_triple(a,multiply_triple(q,q,work),work);
        F twice=add(q.hi,q.hi);
        work=work&nz(residual.valid)&
             ~greater(ratio_upper(abs32(residual.hi),twice),f(1560281088U));
        q=center(add_triple(q,constant(divide_seed(residual.hi,twice,work)),work));
    }
    F lower=absolute_lower(q);
    work=work&nz(q.valid)&~less_equal(q.hi,zero())&~less_equal(lower,zero());
    if(!any(work)) return choose(input_invalid,a,choose(exact_zero,constant_zero(),invalid()));
    Triple residual=subtract_triple(a,multiply_triple(q,q,work),work);
    work=work&nz(residual.valid);
    if(!any(work)) return choose(input_invalid,a,choose(exact_zero,constant_zero(),invalid()));
    Triple result=pack(q.hi,q.lo,q.tail,ratio_upper(absolute_upper(residual),lower));
    result=choose(work,result,invalid());
    result=choose(exact_zero,constant_zero(),result);
    return choose(input_invalid,a,result);
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
    if(_mm256_movemask_ps(std::bit_cast<__m256>(active))==255) {
        for(unsigned c=0;c<5;++c) std::memcpy(p+c*8,&columns[c],sizeof(U));
    } else {
        for(unsigned c=0;c<5;++c) {
            U merged=choose(active,columns[c],load_words(p+c*8));
            std::memcpy(p+c*8,&merged,sizeof(merged));
        }
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

extern "C" __attribute__((visibility("hidden")))
void sirius_native_dense_divide8(std::uint32_t* scratch,std::uint32_t destination,
                                 std::uint32_t a_index,std::uint32_t b_index,
                                 std::uint32_t mask) noexcept {
    Triple a=read(scratch,a_index),b=read(scratch,b_index);U active=active_mask(mask);
    write(scratch,destination,divide_triple(a,b,active),active);
}
extern "C" __attribute__((visibility("hidden")))
void sirius_native_dense_sqrt8(std::uint32_t* scratch,std::uint32_t destination,
                              std::uint32_t a_index,std::uint32_t b_index,
                              std::uint32_t mask) noexcept {
    (void)b_index;
    Triple a=read(scratch,a_index);U active=active_mask(mask);
    write(scratch,destination,sqrt_triple(a,active),active);
}

// RTSubtract is exactly RTAdd(a, RTNegate(b)); load both complete operands
// before the destination write, including when destination aliases a or b.
extern "C" __attribute__((visibility("hidden")))
void sirius_native_dense_subtract8(std::uint32_t* scratch,std::uint32_t destination,
                                   std::uint32_t a_index,std::uint32_t b_index,
                                   std::uint32_t mask) noexcept {
    Triple a=read(scratch,a_index), b=read(scratch,b_index);
    U active=active_mask(mask);
    write(scratch,destination,add_triple(a,negate(b),active),active);
}

// Cached predicates belong only to the exact private register words. The
// original five-word functions above remain unchanged for all other callers.
namespace {
__attribute__((always_inline)) inline Triple facts_add_general(Triple a,Triple b,U work) {
    Triple r=fast_term(a,b.hi,work);r=fast_term(r,b.lo,work);r=fast_term(r,b.tail,work);
    return finish(r.hi,r.lo,r.tail,up_add(r.error,b.error));
}
__attribute__((always_inline)) inline Triple facts_multiply_general(Triple a,Triple b,U work) {
    Triple r=fast_product(constant_zero(),a.hi,b.hi,work);
    r=fast_product(r,a.hi,b.lo,work);r=fast_product(r,a.lo,b.hi,work);
    r=fast_rounded_product(r,a.lo,b.lo,work);
    r=fast_rounded_product(r,a.hi,b.tail,work);
    r=fast_rounded_product(r,a.tail,b.hi,work);
    F omitted=up_add(up_add(up_mul(abs32(a.lo),abs32(b.tail)),up_mul(abs32(a.tail),abs32(b.lo))),
                     up_mul(abs32(a.tail),abs32(b.tail)));
    F a_error=a.error,b_center=center_upper(b),first_error=up_mul(a_error,b_center);
    F b_error=b.error,a_center=center_upper(a);
    F error=up_add(up_add(r.error,omitted),up_add(up_add(first_error,up_mul(b_error,a_center)),up_mul(a_error,b_error)));
    return finish(r.hi,r.lo,r.tail,error);
}
// Metadata bits: normalized1, known2, exact-zero4, exact-one8,
// exact-minus-one16, nonzero-valid32. Known certifies every fact, including
// the exact normalization bit. Complete writes set known; partial writes clear all bits.
// Legacy/unknown metadata retains only its normalization bit and recomputes
// identity facts from the same complete operand before making any decision.
__attribute__((always_inline)) inline U identity_facts(Triple a) {
    U lower=(bits(a.lo)|bits(a.tail)|bits(a.error))&u(2147483647U);
    U clean=nz(a.valid)&std::bit_cast<U>(lower==u(0));
    U h=bits(a.hi);
    return u(2)|(nz(a.valid)&u(32))|(clean&std::bit_cast<U>((h&u(2147483647U))==u(0))&u(4))|
        (clean&std::bit_cast<U>(h==u(1065353216U))&u(8))|
        (clean&std::bit_cast<U>(h==u(3212836864U))&u(16));
}
__attribute__((always_inline)) inline U operand_facts(Triple a,U metadata) {
    U known=nz(metadata&u(2));
    if(all_on(known,u(0xffffffffU))) return metadata;
    return choose(known,metadata,(metadata&u(1))|identity_facts(a));
}
__attribute__((always_inline)) inline U negated_facts(U metadata) {
    return (metadata&~u(24))|((metadata&u(8))<<1)|((metadata&u(16))>>1);
}
__attribute__((always_inline)) inline U register_tag(const std::uint32_t* tags,std::uint32_t index) {
    return load_words(tags+std::size_t(index)*8);
}
__attribute__((always_inline)) inline void write_tag(std::uint32_t* tags,std::uint32_t index,
                                                    Triple value,U active) {
    auto* p=tags+std::size_t(index)*8;
    U result=(normalized(value)&u(1))|identity_facts(value);
    if(_mm256_movemask_ps(std::bit_cast<__m256>(active))!=255)
        result=choose(active,result,load_words(p));
    std::memcpy(p,&result,sizeof(result));
}
__attribute__((always_inline)) inline Triple cached_add(Triple a,Triple b,U ta,U tb,U active) {
    U good=nz(a.valid)&nz(b.valid);
    U ma=operand_facts(a,ta),mb=operand_facts(b,tb);
    U az=nz(ma&u(4)),bz=nz(mb&u(4));
    U work=active&good&~az&~bz;Triple general=invalid();
    if(any(work)) {
        if(!all_on(nz(ta&u(1))&nz(tb&u(1)),work)) return add_triple(a,b,active);
        Triple aa=choose(work,a,constant_zero()),bb=choose(work,b,constant_zero());
        general=facts_add_general(aa,bb,work);
    }
    return choose(good,choose(az,b,choose(bz,a,general)),invalid());
}
__attribute__((always_inline)) inline Triple cached_multiply(Triple a,Triple b,U ta,U tb,U active) {
    U good=nz(a.valid)&nz(b.valid);
    U ma=operand_facts(a,ta),mb=operand_facts(b,tb);
    U is_zero=nz(ma&u(4))|nz(mb&u(4)),a_one=nz(ma&u(8)),b_one=nz(mb&u(8));
    U a_minus=nz(ma&u(16)),b_minus=nz(mb&u(16));
    U work=active&good&~is_zero&~a_one&~b_one&~a_minus&~b_minus;Triple general=invalid();
    if(any(work)) {
        U domain=nz(ta&u(1))&nz(tb&u(1))&bounded(a.hi,1518338048U)&bounded(b.hi,1518338048U);
        if(!all_on(domain,work)) return multiply_triple(a,b,active);
        Triple aa=choose(work,a,constant_zero()),bb=choose(work,b,constant_zero());
        general=facts_multiply_general(aa,bb,work);
    }
    Triple selected=choose(a_minus,negate(b),choose(b_minus,negate(a),general));
    selected=choose(a_one,b,choose(b_one,a,selected));
    selected=choose(is_zero,constant_zero(),selected);
    return choose(good,selected,invalid());
}
__attribute__((always_inline)) inline void store_known_facts(std::uint32_t* tags,
    std::uint32_t index,U metadata,U active) {
    auto* p=tags+std::size_t(index)*8;
    if(_mm256_movemask_ps(std::bit_cast<__m256>(active))!=255)
        metadata=choose(active,metadata,load_words(p));
    std::memcpy(p,&metadata,sizeof(metadata));
}
// Only complete, coherent metadata enters this path. Every active lane must
// select an original identity/refusal; otherwise no store occurs here.
__attribute__((always_inline)) inline bool store_identity(std::uint32_t* scratch,
    std::uint32_t* tags,std::uint32_t destination,std::uint32_t ai,std::uint32_t bi,
    U ma,U mb,U active,unsigned operation) {
    if(!any(active)) return true;
    if(!all_on(nz(ma&u(2))&nz(mb&u(2)),active)) return false;
    U live=active&nz(ma&u(32))&nz(mb&u(32));
    U az=nz(ma&u(4)),bz=nz(mb&u(4));
    U zero_mask=u(0),copy_a=u(0),copy_b=u(0),negative_a=u(0),negative_b=u(0);
    if(operation==4) {
        zero_mask=live&(az|bz);live=live&~zero_mask;
        copy_b=live&nz(ma&u(8));live=live&~copy_b;
        copy_a=live&nz(mb&u(8));live=live&~copy_a;
        negative_b=live&nz(ma&u(16));live=live&~negative_b;
        negative_a=live&nz(mb&u(16));live=live&~negative_a;
    } else {
        copy_b=live&az;live=live&~copy_b;
        copy_a=live&bz;live=live&~copy_a;
        if(operation==3) {negative_b=copy_b;copy_b=u(0);}
    }
    if(any(live)) return false;
    Triple a=invalid(),b=invalid();
    // Snapshot every selected operand before writing, including destination
    // aliases. Zero/refusal-only packets load no payload operands.
    if(any(copy_a|negative_a)) a=read(scratch,ai);
    if(any(copy_b|negative_b)) b=read(scratch,bi);
    Triple result=choose(zero_mask,constant_zero(),choose(copy_b,b,choose(copy_a,a,
        choose(negative_b,negate(b),choose(negative_a,negate(a),invalid())))));
    U metadata=choose(zero_mask,u(39),choose(copy_b,mb,choose(copy_a,ma,
        choose(negative_b,negated_facts(mb),choose(negative_a,negated_facts(ma),u(2))))));
    write(scratch,destination,result,active);
    store_known_facts(tags,destination,metadata,active);
    return true;
}

}

extern "C" __attribute__((visibility("hidden")))
void sirius_native_tagged_add8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept {
    U active=active_mask(mask),ta=register_tag(tags,a_index),tb=register_tag(tags,b_index);
    if(store_identity(scratch,tags,destination,a_index,b_index,ta,tb,active,2)) return;
    Triple a=read(scratch,a_index), b=read(scratch,b_index);
    Triple result=cached_add(a,b,ta,tb,active);
    write(scratch,destination,result,active);
    write_tag(tags,destination,result,active);
}

extern "C" __attribute__((visibility("hidden")))
void sirius_native_tagged_multiply8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept {
    U active=active_mask(mask),ta=register_tag(tags,a_index),tb=register_tag(tags,b_index);
    if(store_identity(scratch,tags,destination,a_index,b_index,ta,tb,active,4)) return;
    Triple a=read(scratch,a_index), b=read(scratch,b_index);
    Triple result=cached_multiply(a,b,ta,tb,active);
    write(scratch,destination,result,active);
    write_tag(tags,destination,result,active);
}

extern "C" __attribute__((visibility("hidden")))
void sirius_native_tagged_subtract8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept {
    U active=active_mask(mask),ta=register_tag(tags,a_index),tb=register_tag(tags,b_index);
    if(store_identity(scratch,tags,destination,a_index,b_index,ta,tb,active,3)) return;
    Triple a=read(scratch,a_index), b=read(scratch,b_index);
    Triple result=cached_add(a,negate(b),ta,negated_facts(tb),active);
    write(scratch,destination,result,active);
    write_tag(tags,destination,result,active);
}

extern "C" __attribute__((visibility("hidden")))
void sirius_native_tagged_divide8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept {
    Triple a=read(scratch,a_index), b=read(scratch,b_index);
    U active=active_mask(mask);
    Triple result=divide_triple(a,b,active);
    write(scratch,destination,result,active);
    write_tag(tags,destination,result,active);
}

extern "C" __attribute__((visibility("hidden")))
void sirius_native_tagged_sqrt8(std::uint32_t* scratch,std::uint32_t* tags,
    std::uint32_t destination,std::uint32_t a_index,std::uint32_t b_index,
    std::uint32_t mask) noexcept {
    (void)b_index;
    Triple a=read(scratch,a_index);
    U active=active_mask(mask);
    Triple result=sqrt_triple(a,active);
    write(scratch,destination,result,active);
    write_tag(tags,destination,result,active);
}
