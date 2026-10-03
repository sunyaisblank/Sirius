// Finite bounded transforms; all low parts remain explicit outputs.
import retained_scalar;
// Magnitude ordering establishes FastTwoSum's required |a| >= |b|.
// The longer Knuth chain lost a correction in the combined camera pipeline.
void RetainedSum32(RPScalar a,RPScalar b,out RPScalar hi,out RPScalar lo){
 if((RPBits(a)&0x7fffffffu)<(RPBits(b)&0x7fffffffu)){RPScalar t=a;a=b;b=t;}
 hi=OrderedAdd32(a,b);lo=OrderedSubtract32(b,OrderedSubtract32(hi,a));
}
void Split32(RPScalar a,out RPScalar hi,out RPScalar lo){
 RPScalar c=OrderedMultiply32(RPFromBits(0x45800800u),a);hi=OrderedSubtract32(c,OrderedSubtract32(c,a));lo=OrderedSubtract32(a,hi);
}
void RetainedProduct32(RPScalar a,RPScalar b,out RPScalar hi,out RPScalar lo){
#if defined(SIRIUS_RETAINED_PORTABLE)
 // The exact bounded integer product/residual avoids emulating Dekker's
 // repeated scalar operations. Both binary32 projections satisfy the same
 // five-minimum-subnormal enclosure used by the native product transforms.
 PB32Product product=PB32ProductResidual(RPBits(a),RPBits(b));
 hi=RPFromBits(product.valid!=0u?product.high:0x7f800000u);
 lo=RPFromBits(product.valid!=0u?product.low:0x7f800000u);
#elif defined(SIRIUS_RETAINED_FP64)
 // A binary32 product has at most 48 significant bits and is normal in
 // binary64 even for two binary32 subnormals. Keep its rounded high part
 // and exact residual separately; the caller bounds binary32 underflow.
 double da=double(a),db=double(b);
 double product=spirv_asm{OpFMul $$double %v $da $db;OpDecorate %v NoContraction;OpCopyObject $$double result %v;};
 // Establish the high part through the ordered binary32 operation itself.
 // A native wide-to-narrow conversion did not preserve this residual split
 // in the complete Radeon/Dozen camera program.
 hi=OrderedMultiply32(a,b);double dh=double(hi);
 double residual=spirv_asm{OpFSub $$double %v $product $dh;OpDecorate %v NoContraction;OpCopyObject $$double result %v;};
 lo=RPScalar(residual);
#else
 RPScalar ah,al,bh,bl;Split32(a,ah,al);Split32(b,bh,bl);hi=OrderedMultiply32(a,b);
 RPScalar e1=OrderedSubtract32(hi,OrderedMultiply32(ah,bh));
 RPScalar e2=OrderedSubtract32(e1,OrderedMultiply32(al,bh));
 RPScalar e3=OrderedSubtract32(e2,OrderedMultiply32(ah,bl));
 lo=OrderedSubtract32(OrderedMultiply32(al,bl),e3);
#endif
}

// Every successful value encloses the real input in hi+lo +/- error.
struct RetainedPair { RPScalar hi; RPScalar lo; RPScalar error; uint valid; }
RPScalar RPAbs(RPScalar x) { return RPFromBits(RPBits(x)&0x7fffffffu); }
bool RPBounded(RPScalar x,uint limit) { return (RPBits(x)&0x7fffffffu)<=limit; }
RPScalar RPFailureBound() { return RPFromBits(0x7f000000u); } // finite 2^127 sentinel
RPScalar RPUpAdd(RPScalar a,RPScalar b) {
    if(!RPBounded(a,0x7d800000u)||!RPBounded(b,0x7d800000u))return RPFailureBound();
    if(RPScalarEqual(a,RPZero()))return b;if(RPScalarEqual(b,RPZero()))return a;
    return RPFromBits(RPBits(OrderedAdd32(a,b))+1u);
}
RPScalar RPUpMultiply(RPScalar a,RPScalar b) {
    if(RPScalarEqual(a,RPZero())||RPScalarEqual(b,RPZero()))return RPZero();
    if(!RPBounded(a,0x5e800000u)||!RPBounded(b,0x5e800000u))return RPFailureBound();
    return RPFromBits(RPBits(OrderedMultiply32(a,b))+1u);
}
RetainedPair RPInvalid() { RetainedPair r;r.hi=RPZero();r.lo=RPZero();r.error=RPZero();r.valid=0;return r; }
RetainedPair RPPack(RPScalar hi,RPScalar lo,RPScalar error) {
    if(!RPBounded(hi,0x7b800000u)||!RPBounded(lo,0x7b800000u)||
       !RPBounded(error,0x7b800000u)||(RPBits(error)&0x80000000u)!=0u)return RPInvalid();
    RetainedPair r;RetainedSum32(hi,lo,r.hi,r.lo);if(!RPBounded(r.hi,0x7b800000u))return RPInvalid();r.error=error;r.valid=1;return r;
}
RetainedPair RPConstant(RPScalar x) {return RPPack(x,RPZero(),RPZero());}
RetainedPair RPCenter(RetainedPair a) {a.error=RPZero();return a;}
RetainedPair RPNegate(RetainedPair a) {
    a.hi=RPFromBits(RPBits(a.hi)^0x80000000u);a.lo=RPFromBits(RPBits(a.lo)^0x80000000u);return a;
}
RPScalar RPCenterUpper(RetainedPair a) {return RPUpAdd(RPAbs(a.hi),RPAbs(a.lo));}
RPScalar RPAbsoluteUpper(RetainedPair a) {return RPUpAdd(RPCenterUpper(a),a.error);}
RPScalar RPAbsoluteLower(RetainedPair a) {
    RPScalar tail=RPUpAdd(RPAbs(a.lo),a.error);
    if(a.valid==0u||RPScalarLessEqual(RPAbs(a.hi),tail))return RPZero();
    RPScalar difference,correction;RetainedSum32(RPAbs(a.hi),RPScalarNegate(tail),difference,correction);
    uint bits=RPBits(difference);return RPScalarGreaterEqual(correction,RPZero())?difference:(bits==0u?RPZero():RPFromBits(bits-1u));
}
RetainedPair RPAppend(RetainedPair a,RPScalar x) {
    if(a.valid==0u||!RPBounded(a.hi,0x7b800000u)||!RPBounded(a.lo,0x7b800000u)||
       !RPBounded(x,0x7b800000u))return RPInvalid();
    RPScalar s,e,t,f,h,l;RetainedSum32(a.hi,x,s,e);RetainedSum32(a.lo,e,t,f);RetainedSum32(s,t,h,l);
    return RPPack(h,l,RPUpAdd(a.error,RPAbs(f)));
}
RetainedPair RPAdd(RetainedPair a,RetainedPair b) {
    if(a.valid==0u||b.valid==0u)return RPInvalid();
    RetainedPair r=RPAppend(RPAppend(a,b.hi),b.lo);
    return r.valid==0u?r:RPPack(r.hi,r.lo,RPUpAdd(r.error,b.error));
}
RetainedPair RPSubtract(RetainedPair a,RetainedPair b) {return RPAdd(a,RPNegate(b));}
RetainedPair RPScalarProduct(RPScalar a,RPScalar b) {
    // 2^60 * splitter < 2^73; products and residual intermediates < 2^123.
    if(!RPBounded(a,0x5d800000u)||!RPBounded(b,0x5d800000u))return RPInvalid();
    if(RPScalarEqual(a,RPZero())||RPScalarEqual(b,RPZero()))return RPConstant(RPZero());
    RPScalar h,l;RetainedProduct32(a,b,h,l);return RPPack(h,l,RPFromBits(5u));
}
RetainedPair RPMultiply(RetainedPair a,RetainedPair b) {
    if(a.valid==0u||b.valid==0u)return RPInvalid();
    RetainedPair r=RPAdd(RPAdd(RPScalarProduct(a.hi,b.hi),RPScalarProduct(a.hi,b.lo)),
                        RPAdd(RPScalarProduct(a.lo,b.hi),RPScalarProduct(a.lo,b.lo)));
    RPScalar uncertainty=RPUpAdd(RPUpAdd(RPUpMultiply(a.error,RPCenterUpper(b)),
        RPUpMultiply(b.error,RPCenterUpper(a))),RPUpMultiply(a.error,b.error));
    return r.valid==0u?r:RPPack(r.hi,r.lo,RPUpAdd(r.error,uncertainty));
}
int RPFloorExponent(RPScalar positive) {
    uint bits=RPBits(positive),exponent=bits>>23;
    return exponent!=0u?int(exponent)-127:int(firstbithigh(bits))-149;
}
RPScalar RPPowerOfTwo(int e) {
    if(e>120)return RPFailureBound();if(e<-149)return RPFromBits(1u);
    return e>=-126?RPFromBits(uint(e+127)<<23):RPFromBits(1u<<uint(e+149));
}
// n < 2^(floor(log2(n))+1), d >= 2^floor(log2(d)). No native division.
RPScalar RPRatioUpper(RPScalar numerator,RPScalar denominator) {
    if(!RPBounded(numerator,0x7b800000u)||!RPBounded(denominator,0x7b800000u)||RPScalarLessEqual(denominator,RPZero()))return RPFailureBound();
    if(RPScalarEqual(numerator,RPZero()))return RPZero();
    return RPPowerOfTwo(RPFloorExponent(numerator)-RPFloorExponent(denominator)+1);
}
// Normalize by integer significand reconstruction; native seed operations only
// see normal inputs. Ordered scaling restores the original exponent afterward.
RPScalar RPMantissa(RPScalar positive) {
    uint bits=RPBits(positive),exponent=bits>>23;
    uint fraction=exponent!=0u?(bits&0x7fffffu):((bits<<uint(23-int(firstbithigh(bits))))&0x7fffffu);
    return RPFromBits(0x3f800000u|fraction);
}
RPScalar RPDivideSeed(RPScalar a,RPScalar b) {
    if(RPScalarEqual(a,RPZero()))return RPZero();
    RPScalar aa=RPAbs(a),bb=RPAbs(b);
    int exponent=RPFloorExponent(aa)-RPFloorExponent(bb);
    if(exponent < -151 || exponent > 120)return RPZero();
    RPScalar am=RPMantissa(aa),bm=RPMantissa(bb);
    RPScalar seed=OrderedDivide32(am,bm);
    if(!RPBounded(seed,0x40800000u))return RPZero();
    if(exponent < -126){seed=OrderedMultiply32(seed,RPFromBits(0x00800000u));exponent+=126;}
    seed=OrderedMultiply32(seed,RPPowerOfTwo(exponent));
    return RPFromBits(RPBits(seed)^((RPBits(a)^RPBits(b))&0x80000000u));
}
RPScalar RPSqrtSeed(RPScalar positive) {
    int exponent=RPFloorExponent(positive),odd=exponent&1;
    RPScalar mantissa=RPMantissa(positive);
    if(odd!=0)mantissa=OrderedMultiply32(mantissa,RPFromBits(0x40000000u));
    RPScalar seed=OrderedSqrt32(mantissa);if(!RPBounded(seed,0x40800000u))return RPZero();
    return OrderedMultiply32(seed,RPPowerOfTwo((exponent-odd)/2));
}
RetainedPair RPDivide(RetainedPair a,RetainedPair b) {
    if(a.valid==0u||b.valid==0u)return RPInvalid();
    RPScalar lower=RPAbsoluteLower(b);
    if(RPScalarLessEqual(lower,RPZero())||RPScalarGreater(RPRatioUpper(RPAbsoluteUpper(a),lower),RPFromBits(0x5e000000u)))return RPInvalid();
    RetainedPair q=RPConstant(RPDivideSeed(a.hi,b.hi));if(q.valid==0u)return q;
    RetainedPair residual=RPSubtract(a,RPMultiply(b,q));if(residual.valid==0u)return residual;
    if(RPScalarGreater(RPRatioUpper(RPAbs(residual.hi),RPAbs(b.hi)),RPFromBits(0x5d000000u)))return RPInvalid();
    q=RPCenter(RPAdd(q,RPConstant(RPDivideSeed(residual.hi,b.hi))));if(q.valid==0u)return q;
    residual=RPSubtract(a,RPMultiply(b,q));if(residual.valid==0u)return residual;
    return RPPack(q.hi,q.lo,RPRatioUpper(RPAbsoluteUpper(residual),lower));
}
RetainedPair RPSqrt(RetainedPair a) {
    if(a.valid==0u)return a;
    if(RPScalarEqual(a.hi,RPZero())&&RPScalarEqual(a.lo,RPZero())&&RPScalarEqual(a.error,RPZero()))return RPConstant(RPZero());
    if(RPScalarLessEqual(a.hi,RPZero())||RPScalarLessEqual(RPAbsoluteLower(a),RPZero()))return RPInvalid();
    RetainedPair q=RPConstant(RPSqrtSeed(a.hi));if(q.valid==0u||RPScalarLessEqual(q.hi,RPZero())||!RPBounded(q.hi,0x5d800000u))return RPInvalid();
    RetainedPair residual=RPSubtract(a,RPMultiply(q,q));if(residual.valid==0u)return residual;
    RPScalar twice=OrderedAdd32(q.hi,q.hi);
    if(RPScalarGreater(RPRatioUpper(RPAbs(residual.hi),twice),RPFromBits(0x5d000000u)))return RPInvalid();
    q=RPCenter(RPAdd(q,RPConstant(RPDivideSeed(residual.hi,twice))));
    RPScalar lower=RPAbsoluteLower(q);if(q.valid==0u||RPScalarLessEqual(q.hi,RPZero())||RPScalarLessEqual(lower,RPZero()))return RPInvalid();
    residual=RPSubtract(a,RPMultiply(q,q));if(residual.valid==0u)return residual;
    return RPPack(q.hi,q.lo,RPRatioUpper(RPAbsoluteUpper(residual),lower));
}
