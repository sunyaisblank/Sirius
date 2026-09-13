// Three-term expansion for transport's small covariant components. The camera
// pair remains a separate authority; no conversion discards its low part.
import retained_pair;

struct RetainedTriple { float hi; float lo; float tail; float error; uint valid; }
RetainedTriple RTInvalid() {
    RetainedTriple r; r.hi=0; r.lo=0; r.tail=0; r.error=0; r.valid=0; return r;
}
// Callers bound all three private limbs by 2^124 before normalization, so
// every intermediate sum is finite. The returned value retains the public
// 2^120 component/radius limits, including after large-term cancellation.
RetainedTriple RTFinishAccumulation(float hi,float lo,float tail,float error) {
    if(!RPBounded(error,0x7b800000u)||(asuint(error)&0x80000000u)!=0u) return RTInvalid();
    float a,b,c,d,e,f,g;
    RetainedSum32(lo,tail,a,b);
    RetainedSum32(hi,a,c,d);
    RetainedSum32(d,b,e,f);
    RetainedTriple r;
    RetainedSum32(c,e,r.hi,g);
    RetainedSum32(g,f,r.lo,r.tail);r.error=error;r.valid=1;
    if(!RPBounded(r.hi,0x7b800000u)) return RTInvalid();
    return r;
}
RetainedTriple RTPack(float hi,float lo,float tail,float error) {
    if(!RPBounded(hi,0x7b800000u)||!RPBounded(lo,0x7b800000u)||
       !RPBounded(tail,0x7b800000u)) return RTInvalid();
    return RTFinishAccumulation(hi,lo,tail,error);
}
RetainedTriple RTConstant(float x) {return RTPack(x,0,0,0);}
bool RTExact(RetainedTriple a,float x) {
    return a.valid!=0&&a.hi==x&&a.lo==0&&a.tail==0&&a.error==0;
}
RetainedTriple RTCenter(RetainedTriple a) {a.error=0;return a;}
RetainedTriple RTNegate(RetainedTriple a) {
    a.hi=asfloat(asuint(a.hi)^0x80000000u);
    a.lo=asfloat(asuint(a.lo)^0x80000000u);
    a.tail=asfloat(asuint(a.tail)^0x80000000u);return a;
}
float RTCenterUpper(RetainedTriple a) {
    return RPUpAdd(RPUpAdd(RPAbs(a.hi),RPAbs(a.lo)),RPAbs(a.tail));
}
float RTAbsoluteUpper(RetainedTriple a) {return RPUpAdd(RTCenterUpper(a),a.error);}
float RTAbsoluteLower(RetainedTriple a) {
    float tail=RPUpAdd(RPUpAdd(RPAbs(a.lo),RPAbs(a.tail)),a.error);
    if(a.valid==0||RPAbs(a.hi)<=tail) return 0;
    float difference,correction;RetainedSum32(RPAbs(a.hi),-tail,difference,correction);
    uint bits=asuint(difference);
    return correction>=0?difference:(bits==0?0:asfloat(bits-1));
}
// A private sum need not normalize after every scalar term. Magnitude-ordered
// TwoSum preserves overlapping limbs; each discarded fourth term is enclosed
// immediately. Normalize once before magnitude estimates consume the limbs.
RetainedTriple RTAccumulateTerm(RetainedTriple a,float x) {
    if(a.valid==0||!RPBounded(x,0x7b800000u)) return RTInvalid();
    float hi,lo,tail,first,second,last;
    RetainedSum32(a.hi,x,hi,first);
    RetainedSum32(a.lo,first,lo,second);
    RetainedSum32(a.tail,second,tail,last);
    float error=RPUpAdd(a.error,RPAbs(last));
    // At most nine bounded terms enter a product. This private 2^124 bound
    // leaves room for cancellation while keeping every sum finite. Finishing
    // the accumulation reapplies the unchanged 2^120 public bound.
    if(!RPBounded(hi,0x7d800000u)||!RPBounded(lo,0x7d800000u)||
       !RPBounded(tail,0x7d800000u)||!RPBounded(error,0x7d800000u)) return RTInvalid();
    RetainedTriple r;r.hi=hi;r.lo=lo;r.tail=tail;r.error=error;r.valid=1;return r;
}
RetainedTriple RTAdd(RetainedTriple a,RetainedTriple b) {
    if(a.valid==0||b.valid==0) return RTInvalid();
    if(RTExact(a,0)) return b;
    if(RTExact(b,0)) return a;
    RetainedTriple r=RTAccumulateTerm(RTAccumulateTerm(RTAccumulateTerm(a,b.hi),b.lo),b.tail);
    return r.valid==0?r:RTFinishAccumulation(r.hi,r.lo,r.tail,RPUpAdd(r.error,b.error));
}
RetainedTriple RTSubtract(RetainedTriple a,RetainedTriple b) {return RTAdd(a,RTNegate(b));}
RetainedTriple RTAccumulateProduct(RetainedTriple sum,float a,float b) {
    if(sum.valid==0||!RPBounded(a,0x5d800000u)||!RPBounded(b,0x5d800000u)) return RTInvalid();
    if(a==0||b==0) return sum;
    float hi,lo;RetainedProduct32(a,b,hi,lo);
    RetainedTriple r=RTAccumulateTerm(RTAccumulateTerm(sum,hi),lo);
    // Same product transform and underflow enclosure as RPScalarProduct;
    // the private accumulator does not require a normalized scalar pair.
    if(r.valid!=0) r.error=RPUpAdd(r.error,asfloat(5u));
    return r;
}
RetainedTriple RTAccumulateRoundedProduct(RetainedTriple sum,float a,float b) {
    if(sum.valid==0||!RPBounded(a,0x5d800000u)||!RPBounded(b,0x5d800000u)) return RTInvalid();
    if(a==0||b==0) return sum;
    float product=OrderedMultiply32(a,b);
    // A full binary32 ulp plus the underflow allowance bounds this rounded
    // low-order product. Its omitted residual stays in the arithmetic radius.
    float error=RPUpAdd(RPUpMultiply(RPAbs(product),1.1920928955078125e-7f),asfloat(5u));
    RetainedTriple r=RTAccumulateTerm(sum,product);
    if(r.valid!=0) r.error=RPUpAdd(r.error,error);
    return r;
}
RetainedTriple RTMultiply(RetainedTriple a,RetainedTriple b) {
    if(a.valid==0||b.valid==0) return RTInvalid();
    if(RTExact(a,0)||RTExact(b,0)) return RTConstant(0);
    if(RTExact(a,1)) return b;
    if(RTExact(b,1)) return a;
    if(RTExact(a,-1)) return RTNegate(b);
    if(RTExact(b,-1)) return RTNegate(a);
    RetainedTriple r=RTAccumulateProduct(RTConstant(0),a.hi,b.hi);
    r=RTAccumulateProduct(r,a.hi,b.lo);
    r=RTAccumulateProduct(r,a.lo,b.hi);
    r=RTAccumulateRoundedProduct(r,a.lo,b.lo);
    r=RTAccumulateRoundedProduct(r,a.hi,b.tail);
    r=RTAccumulateRoundedProduct(r,a.tail,b.hi);
    if(r.valid==0) return r;
    float omitted=RPUpAdd(RPUpAdd(RPUpMultiply(RPAbs(a.lo),RPAbs(b.tail)),
        RPUpMultiply(RPAbs(a.tail),RPAbs(b.lo))),RPUpMultiply(RPAbs(a.tail),RPAbs(b.tail)));
    float uncertainty=RPUpAdd(RPUpAdd(RPUpMultiply(a.error,RTCenterUpper(b)),
        RPUpMultiply(b.error,RTCenterUpper(a))),RPUpMultiply(a.error,b.error));
    return RTFinishAccumulation(r.hi,r.lo,r.tail,RPUpAdd(RPUpAdd(r.error,omitted),uncertainty));
}
RetainedTriple RTDivide(RetainedTriple a,RetainedTriple b) {
    if(a.valid==0||b.valid==0) return RTInvalid();
    float lower=RTAbsoluteLower(b);
    if(lower<=0||RPRatioUpper(RTAbsoluteUpper(a),lower)>asfloat(0x5e000000u)) return RTInvalid();
    if(RTExact(b,1)) return a;
    RetainedTriple q=RTConstant(RPDivideSeed(a.hi,b.hi));
    [loop] for(uint i=0;i<3;++i) {
        RetainedTriple residual=RTSubtract(a,RTMultiply(b,q));
        if(residual.valid==0||RPRatioUpper(RPAbs(residual.hi),RPAbs(b.hi))>asfloat(0x5d000000u))
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
    if(a.hi==0&&a.lo==0&&a.tail==0&&a.error==0) return RTConstant(0);
    if(a.hi<=0||RTAbsoluteLower(a)<=0) return RTInvalid();
    RetainedTriple q=RTConstant(RPSqrtSeed(a.hi));
    [loop] for(uint i=0;i<3;++i) {
        if(q.valid==0||q.hi<=0||!RPBounded(q.hi,0x5d000000u)) return RTInvalid();
        RetainedTriple residual=RTSubtract(a,RTMultiply(q,q));
        float twice=OrderedAdd32(q.hi,q.hi);
        if(residual.valid==0||RPRatioUpper(RPAbs(residual.hi),twice)>asfloat(0x5d000000u))
            return RTInvalid();
        q=RTCenter(RTAdd(q,RTConstant(RPDivideSeed(residual.hi,twice))));
    }
    float lower=RTAbsoluteLower(q);
    if(q.valid==0||q.hi<=0||lower<=0) return RTInvalid();
    RetainedTriple residual=RTSubtract(a,RTMultiply(q,q));
    if(residual.valid==0) return residual;
    return RTPack(q.hi,q.lo,q.tail,RPRatioUpper(RTAbsoluteUpper(residual),lower));
}
