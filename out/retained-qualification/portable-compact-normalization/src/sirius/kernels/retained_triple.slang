// Three-term expansion for transport's small covariant components. The camera
// pair remains a separate authority; no conversion discards its low part.
import retained_pair;
import retained_scalar;

struct RetainedTriple { RPScalar hi; RPScalar lo; RPScalar tail; RPScalar error; uint valid; }
RetainedTriple RTInvalid() {
    RetainedTriple r; r.hi=RPZero(); r.lo=RPZero(); r.tail=RPZero(); r.error=RPZero(); r.valid=0; return r;
}
// Callers bound all three private limbs by 2^124 before normalization, so
// every intermediate sum is finite. The returned value retains the public
// 2^120 component/radius limits, including after large-term cancellation.
RetainedTriple RTFinishAccumulation(RPScalar hi,RPScalar lo,RPScalar tail,RPScalar error) {
    if(!RPBounded(error,0x7b800000u)||(RPBits(error)&0x80000000u)!=0u) return RTInvalid();
#ifdef SIRIUS_RETAINED_PORTABLE
    RetainedTriple r=RTInvalid();
    RPScalar x=lo,y=tail,carry=RPZero(),other=RPZero();
    // Same five TwoSum calls and dependencies as the explicit sequence. Reuse
    // one integer body instead of expanding that body at every call site.
    [loop] for(uint i=0;i<5;++i) {
        RPScalar upper,lower;
        RetainedSum32(x,y,upper,lower);
        if(i==0) {x=hi;y=upper;other=lower;}
        else if(i==1) {x=lower;y=other;carry=upper;}
        else if(i==2) {x=carry;y=upper;other=lower;}
        else if(i==3) {r.hi=upper;x=lower;y=other;}
        else {r.lo=upper;r.tail=lower;}
    }
    r.error=error;r.valid=1;
#else
    RPScalar a,b,c,d,e,f,g;
    RetainedSum32(lo,tail,a,b);
    RetainedSum32(hi,a,c,d);
    RetainedSum32(d,b,e,f);
    RetainedTriple r;
    RetainedSum32(c,e,r.hi,g);
    RetainedSum32(g,f,r.lo,r.tail);r.error=error;r.valid=1;
#endif
    if(!RPBounded(r.hi,0x7b800000u)) return RTInvalid();
    return r;
}
RetainedTriple RTPack(RPScalar hi,RPScalar lo,RPScalar tail,RPScalar error) {
    if(!RPBounded(hi,0x7b800000u)||!RPBounded(lo,0x7b800000u)||
       !RPBounded(tail,0x7b800000u)) return RTInvalid();
    return RTFinishAccumulation(hi,lo,tail,error);
}
RetainedTriple RTConstant(RPScalar x) {return RTPack(x,RPZero(),RPZero(),RPZero());}
bool RTExact(RetainedTriple a,RPScalar x) {
    return a.valid!=0&&RPScalarEqual(a.hi,x)&&RPScalarEqual(a.lo,RPZero())&&RPScalarEqual(a.tail,RPZero())&&RPScalarEqual(a.error,RPZero());
}
RetainedTriple RTCenter(RetainedTriple a) {a.error=RPZero();return a;}
RetainedTriple RTNegate(RetainedTriple a) {
    a.hi=RPFromBits(RPBits(a.hi)^0x80000000u);
    a.lo=RPFromBits(RPBits(a.lo)^0x80000000u);
    a.tail=RPFromBits(RPBits(a.tail)^0x80000000u);return a;
}
RPScalar RTCenterUpper(RetainedTriple a) {
    return RPUpAdd(RPUpAdd(RPAbs(a.hi),RPAbs(a.lo)),RPAbs(a.tail));
}
RPScalar RTAbsoluteUpper(RetainedTriple a) {return RPUpAdd(RTCenterUpper(a),a.error);}
RPScalar RTAbsoluteLower(RetainedTriple a) {
    RPScalar tail=RPUpAdd(RPUpAdd(RPAbs(a.lo),RPAbs(a.tail)),a.error);
    if(a.valid==0||RPScalarLessEqual(RPAbs(a.hi),tail)) return RPZero();
    RPScalar difference,correction;RetainedSum32(RPAbs(a.hi),RPScalarNegate(tail),difference,correction);
    uint bits=RPBits(difference);
    return RPScalarGreaterEqual(correction,RPZero())?difference:(bits==0?RPZero():RPFromBits(bits-1));
}
// A private sum need not normalize after every scalar term. Magnitude-ordered
// TwoSum preserves overlapping limbs; each discarded fourth term is enclosed
// immediately. Normalize once before magnitude estimates consume the limbs.
RetainedTriple RTAccumulateTerm(RetainedTriple a,RPScalar x) {
    if(a.valid==0||!RPBounded(x,0x7b800000u)) return RTInvalid();
    RPScalar hi,lo,tail,first,second,last;
    RetainedSum32(a.hi,x,hi,first);
    RetainedSum32(a.lo,first,lo,second);
    RetainedSum32(a.tail,second,tail,last);
    RPScalar error=RPUpAdd(a.error,RPAbs(last));
    // At most nine bounded terms enter a product. This private 2^124 bound
    // leaves room for cancellation while keeping every sum finite. Finishing
    // the accumulation reapplies the unchanged 2^120 public bound.
    if(!RPBounded(hi,0x7d800000u)||!RPBounded(lo,0x7d800000u)||
       !RPBounded(tail,0x7d800000u)||!RPBounded(error,0x7d800000u)) return RTInvalid();
    RetainedTriple r;r.hi=hi;r.lo=lo;r.tail=tail;r.error=error;r.valid=1;return r;
}
RetainedTriple RTAdd(RetainedTriple a,RetainedTriple b) {
    if(a.valid==0||b.valid==0) return RTInvalid();
    if(RTExact(a,RPZero())) return b;
    if(RTExact(b,RPZero())) return a;
#ifdef SIRIUS_RETAINED_PORTABLE
    // Reuse the integer arithmetic body while preserving hi, lo, tail order.
    RetainedTriple r=a;
    [loop] for(uint i=0;i<3;++i)
        r=RTAccumulateTerm(r,i==0?b.hi:(i==1?b.lo:b.tail));
#else
    RetainedTriple r=RTAccumulateTerm(RTAccumulateTerm(RTAccumulateTerm(a,b.hi),b.lo),b.tail);
#endif
    return r.valid==0?r:RTFinishAccumulation(r.hi,r.lo,r.tail,RPUpAdd(r.error,b.error));
}
RetainedTriple RTSubtract(RetainedTriple a,RetainedTriple b) {return RTAdd(a,RTNegate(b));}
RetainedTriple RTAccumulateProduct(RetainedTriple sum,RPScalar a,RPScalar b) {
    if(sum.valid==0||!RPBounded(a,0x5d800000u)||!RPBounded(b,0x5d800000u)) return RTInvalid();
    if(RPScalarEqual(a,RPZero())||RPScalarEqual(b,RPZero())) return sum;
    RPScalar hi,lo;RetainedProduct32(a,b,hi,lo);
#ifdef SIRIUS_RETAINED_PORTABLE
    RetainedTriple r=sum;
    [loop] for(uint i=0;i<2;++i)
        r=RTAccumulateTerm(r,i==0?hi:lo);
#else
    RetainedTriple r=RTAccumulateTerm(RTAccumulateTerm(sum,hi),lo);
#endif
    // Same product transform and underflow enclosure as RPScalarProduct;
    // the private accumulator does not require a normalized scalar pair.
    if(r.valid!=0) r.error=RPUpAdd(r.error,RPFromBits(5u));
    return r;
}
RetainedTriple RTAccumulateRoundedProduct(RetainedTriple sum,RPScalar a,RPScalar b) {
    if(sum.valid==0||!RPBounded(a,0x5d800000u)||!RPBounded(b,0x5d800000u)) return RTInvalid();
    if(RPScalarEqual(a,RPZero())||RPScalarEqual(b,RPZero())) return sum;
    RPScalar product=OrderedMultiply32(a,b);
    // A full binary32 ulp plus the underflow allowance bounds this rounded
    // low-order product. Its omitted residual stays in the arithmetic radius.
    RPScalar error=RPUpAdd(RPUpMultiply(RPAbs(product),RPFromBits(0x34000000u)),RPFromBits(5u));
    RetainedTriple r=RTAccumulateTerm(sum,product);
    if(r.valid!=0) r.error=RPUpAdd(r.error,error);
    return r;
}
RetainedTriple RTMultiply(RetainedTriple a,RetainedTriple b) {
    if(a.valid==0||b.valid==0) return RTInvalid();
    if(RTExact(a,RPZero())||RTExact(b,RPZero())) return RTConstant(RPZero());
    if(RTExact(a,RPOne())) return b;
    if(RTExact(b,RPOne())) return a;
    if(RTExact(a,RPFromSigned(-1))) return RTNegate(b);
    if(RTExact(b,RPFromSigned(-1))) return RTNegate(a);
#ifdef SIRIUS_RETAINED_PORTABLE
    RetainedTriple r=RTConstant(RPZero());
    // Same three full products, then the same three rounded products. Keeping
    // their calls in loop bodies bounds duplicated compiler input after inlining.
    [loop] for(uint i=0;i<3;++i)
        r=RTAccumulateProduct(r,i==2?a.lo:a.hi,i==1?b.lo:b.hi);
    [loop] for(uint i=0;i<3;++i)
        r=RTAccumulateRoundedProduct(r,i==0?a.lo:(i==1?a.hi:a.tail),
                                    i==0?b.lo:(i==1?b.tail:b.hi));
#else
    RetainedTriple r=RTAccumulateProduct(RTConstant(RPZero()),a.hi,b.hi);
    r=RTAccumulateProduct(r,a.hi,b.lo);
    r=RTAccumulateProduct(r,a.lo,b.hi);
    r=RTAccumulateRoundedProduct(r,a.lo,b.lo);
    r=RTAccumulateRoundedProduct(r,a.hi,b.tail);
    r=RTAccumulateRoundedProduct(r,a.tail,b.hi);
#endif
    if(r.valid==0) return r;
    RPScalar omitted=RPUpAdd(RPUpAdd(RPUpMultiply(RPAbs(a.lo),RPAbs(b.tail)),
        RPUpMultiply(RPAbs(a.tail),RPAbs(b.lo))),RPUpMultiply(RPAbs(a.tail),RPAbs(b.tail)));
    RPScalar uncertainty=RPUpAdd(RPUpAdd(RPUpMultiply(a.error,RTCenterUpper(b)),
        RPUpMultiply(b.error,RTCenterUpper(a))),RPUpMultiply(a.error,b.error));
    return RTFinishAccumulation(r.hi,r.lo,r.tail,RPUpAdd(RPUpAdd(r.error,omitted),uncertainty));
}
RetainedTriple RTDivide(RetainedTriple a,RetainedTriple b) {
    if(a.valid==0||b.valid==0) return RTInvalid();
    RPScalar lower=RTAbsoluteLower(b);
    if(RPScalarLessEqual(lower,RPZero())||RPScalarGreater(RPRatioUpper(RTAbsoluteUpper(a),lower),RPFromBits(0x5e000000u))) return RTInvalid();
    if(RTExact(b,RPOne())) return a;
    RetainedTriple q=RTConstant(RPDivideSeed(a.hi,b.hi));
    [loop] for(uint i=0;i<3;++i) {
        RetainedTriple residual=RTSubtract(a,RTMultiply(b,q));
        if(residual.valid==0||RPScalarGreater(RPRatioUpper(RPAbs(residual.hi),RPAbs(b.hi)),RPFromBits(0x5d000000u)))
            return RTInvalid();
        q=RTCenter(RTAdd(q,RTConstant(RPDivideSeed(residual.hi,b.hi))));
        if(q.valid==0) return q;
    }
    RetainedTriple residual=RTSubtract(a,RTMultiply(b,q));
    if(residual.valid==0) return residual;
    return RTPack(q.hi,q.lo,q.tail,RPRatioUpper(RTAbsoluteUpper(residual),lower));
}
RetainedTriple RTSqrt(RetainedTriple a) {
    if(a.valid==0) return a;
    if(RPScalarEqual(a.hi,RPZero())&&RPScalarEqual(a.lo,RPZero())&&RPScalarEqual(a.tail,RPZero())&&RPScalarEqual(a.error,RPZero())) return RTConstant(RPZero());
    if(RPScalarLessEqual(a.hi,RPZero())||RPScalarLessEqual(RTAbsoluteLower(a),RPZero())) return RTInvalid();
    RetainedTriple q=RTConstant(RPSqrtSeed(a.hi));
    [loop] for(uint i=0;i<3;++i) {
        if(q.valid==0||RPScalarLessEqual(q.hi,RPZero())||!RPBounded(q.hi,0x5d000000u)) return RTInvalid();
        RetainedTriple residual=RTSubtract(a,RTMultiply(q,q));
        RPScalar twice=OrderedAdd32(q.hi,q.hi);
        if(residual.valid==0||RPScalarGreater(RPRatioUpper(RPAbs(residual.hi),twice),RPFromBits(0x5d000000u)))
            return RTInvalid();
        q=RTCenter(RTAdd(q,RTConstant(RPDivideSeed(residual.hi,twice))));
    }
    RPScalar lower=RTAbsoluteLower(q);
    if(q.valid==0||RPScalarLessEqual(q.hi,RPZero())||RPScalarLessEqual(lower,RPZero())) return RTInvalid();
    RetainedTriple residual=RTSubtract(a,RTMultiply(q,q));
    if(residual.valid==0) return residual;
    return RTPack(q.hi,q.lo,q.tail,RPRatioUpper(RTAbsoluteUpper(residual),lower));
}
