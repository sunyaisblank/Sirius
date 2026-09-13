// Finite bounded transforms; all low parts remain explicit outputs.
float OrderedAdd32(float a,float b){return spirv_asm{OpFAdd $$float %v $a $b;OpDecorate %v NoContraction;OpCopyObject $$float result %v;};}
float OrderedSubtract32(float a,float b){return spirv_asm{OpFSub $$float %v $a $b;OpDecorate %v NoContraction;OpCopyObject $$float result %v;};}
float OrderedMultiply32(float a,float b){return spirv_asm{OpFMul $$float %v $a $b;OpDecorate %v NoContraction;OpCopyObject $$float result %v;};}
// Magnitude ordering establishes FastTwoSum's required |a| >= |b|.
// The longer Knuth chain lost a correction in the combined camera pipeline.
void RetainedSum32(float a,float b,out float hi,out float lo){
 if((asuint(a)&0x7fffffffu)<(asuint(b)&0x7fffffffu)){float t=a;a=b;b=t;}
 hi=OrderedAdd32(a,b);lo=OrderedSubtract32(b,OrderedSubtract32(hi,a));
}
void Split32(float a,out float hi,out float lo){
 float c=OrderedMultiply32(4097.0f,a);hi=OrderedSubtract32(c,OrderedSubtract32(c,a));lo=OrderedSubtract32(a,hi);
}
void RetainedProduct32(float a,float b,out float hi,out float lo){
#ifdef SIRIUS_RETAINED_FP64
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
 lo=float(residual);
#else
 float ah,al,bh,bl;Split32(a,ah,al);Split32(b,bh,bl);hi=OrderedMultiply32(a,b);
 float e1=OrderedSubtract32(hi,OrderedMultiply32(ah,bh));
 float e2=OrderedSubtract32(e1,OrderedMultiply32(al,bh));
 float e3=OrderedSubtract32(e2,OrderedMultiply32(ah,bl));
 lo=OrderedSubtract32(OrderedMultiply32(al,bl),e3);
#endif
}

// Every successful value encloses the real input in hi+lo +/- error.
struct RetainedPair { float hi; float lo; float error; uint valid; }
float RPAbs(float x) { return asfloat(asuint(x)&0x7fffffffu); }
bool RPBounded(float x,uint limit) { return (asuint(x)&0x7fffffffu)<=limit; }
float RPFailureBound() { return asfloat(0x7f000000u); } // finite 2^127 sentinel
float RPUpAdd(float a,float b) {
    if(!RPBounded(a,0x7d800000u)||!RPBounded(b,0x7d800000u))return RPFailureBound();
    if(a==0.0f)return b;if(b==0.0f)return a;
    return asfloat(asuint(OrderedAdd32(a,b))+1u);
}
float RPUpMultiply(float a,float b) {
    if(a==0.0f||b==0.0f)return 0.0f;
    if(!RPBounded(a,0x5e800000u)||!RPBounded(b,0x5e800000u))return RPFailureBound();
    return asfloat(asuint(OrderedMultiply32(a,b))+1u);
}
RetainedPair RPInvalid() { RetainedPair r;r.hi=0;r.lo=0;r.error=0;r.valid=0;return r; }
RetainedPair RPPack(float hi,float lo,float error) {
    if(!RPBounded(hi,0x7b800000u)||!RPBounded(lo,0x7b800000u)||
       !RPBounded(error,0x7b800000u)||(asuint(error)&0x80000000u)!=0u)return RPInvalid();
    RetainedPair r;RetainedSum32(hi,lo,r.hi,r.lo);if(!RPBounded(r.hi,0x7b800000u))return RPInvalid();r.error=error;r.valid=1;return r;
}
RetainedPair RPConstant(float x) {return RPPack(x,0.0f,0.0f);}
RetainedPair RPCenter(RetainedPair a) {a.error=0.0f;return a;}
RetainedPair RPNegate(RetainedPair a) {
    a.hi=asfloat(asuint(a.hi)^0x80000000u);a.lo=asfloat(asuint(a.lo)^0x80000000u);return a;
}
float RPCenterUpper(RetainedPair a) {return RPUpAdd(RPAbs(a.hi),RPAbs(a.lo));}
float RPAbsoluteUpper(RetainedPair a) {return RPUpAdd(RPCenterUpper(a),a.error);}
float RPAbsoluteLower(RetainedPair a) {
    float tail=RPUpAdd(RPAbs(a.lo),a.error);
    if(a.valid==0u||RPAbs(a.hi)<=tail)return 0.0f;
    float difference,correction;RetainedSum32(RPAbs(a.hi),-tail,difference,correction);
    uint bits=asuint(difference);return correction>=0.0f?difference:(bits==0u?0.0f:asfloat(bits-1u));
}
RetainedPair RPAppend(RetainedPair a,float x) {
    if(a.valid==0u||!RPBounded(a.hi,0x7b800000u)||!RPBounded(a.lo,0x7b800000u)||
       !RPBounded(x,0x7b800000u))return RPInvalid();
    float s,e,t,f,h,l;RetainedSum32(a.hi,x,s,e);RetainedSum32(a.lo,e,t,f);RetainedSum32(s,t,h,l);
    return RPPack(h,l,RPUpAdd(a.error,RPAbs(f)));
}
RetainedPair RPAdd(RetainedPair a,RetainedPair b) {
    if(a.valid==0u||b.valid==0u)return RPInvalid();
    RetainedPair r=RPAppend(RPAppend(a,b.hi),b.lo);
    return r.valid==0u?r:RPPack(r.hi,r.lo,RPUpAdd(r.error,b.error));
}
RetainedPair RPSubtract(RetainedPair a,RetainedPair b) {return RPAdd(a,RPNegate(b));}
RetainedPair RPScalarProduct(float a,float b) {
    // 2^60 * splitter < 2^73; products and residual intermediates < 2^123.
    if(!RPBounded(a,0x5d800000u)||!RPBounded(b,0x5d800000u))return RPInvalid();
    if(a==0.0f||b==0.0f)return RPConstant(0.0f);
    float h,l;RetainedProduct32(a,b,h,l);return RPPack(h,l,asfloat(5u));
}
RetainedPair RPMultiply(RetainedPair a,RetainedPair b) {
    if(a.valid==0u||b.valid==0u)return RPInvalid();
    RetainedPair r=RPAdd(RPAdd(RPScalarProduct(a.hi,b.hi),RPScalarProduct(a.hi,b.lo)),
                        RPAdd(RPScalarProduct(a.lo,b.hi),RPScalarProduct(a.lo,b.lo)));
    float uncertainty=RPUpAdd(RPUpAdd(RPUpMultiply(a.error,RPCenterUpper(b)),
        RPUpMultiply(b.error,RPCenterUpper(a))),RPUpMultiply(a.error,b.error));
    return r.valid==0u?r:RPPack(r.hi,r.lo,RPUpAdd(r.error,uncertainty));
}
int RPFloorExponent(float positive) {
    uint bits=asuint(positive),exponent=bits>>23;
    return exponent!=0u?int(exponent)-127:int(firstbithigh(bits))-149;
}
float RPPowerOfTwo(int e) {
    if(e>120)return RPFailureBound();if(e<-149)return asfloat(1u);
    return e>=-126?asfloat(uint(e+127)<<23):asfloat(1u<<uint(e+149));
}
// n < 2^(floor(log2(n))+1), d >= 2^floor(log2(d)). No native division.
float RPRatioUpper(float numerator,float denominator) {
    if(!RPBounded(numerator,0x7b800000u)||!RPBounded(denominator,0x7b800000u)||denominator<=0.0f)return RPFailureBound();
    if(numerator==0.0f)return 0.0f;
    return RPPowerOfTwo(RPFloorExponent(numerator)-RPFloorExponent(denominator)+1);
}
// Normalize by integer significand reconstruction; native seed operations only
// see normal inputs. Ordered scaling restores the original exponent afterward.
float RPMantissa(float positive) {
    uint bits=asuint(positive),exponent=bits>>23;
    uint fraction=exponent!=0u?(bits&0x7fffffu):((bits<<uint(23-int(firstbithigh(bits))))&0x7fffffu);
    return asfloat(0x3f800000u|fraction);
}
float RPDivideSeed(float a,float b) {
    if(a==0.0f)return 0.0f;
    float aa=RPAbs(a),bb=RPAbs(b);
    int exponent=RPFloorExponent(aa)-RPFloorExponent(bb);
    if(exponent < -151 || exponent > 120)return 0.0f;
    float am=RPMantissa(aa),bm=RPMantissa(bb);
    float seed=spirv_asm {OpFDiv $$float %v $am $bm;OpDecorate %v NoContraction;OpCopyObject $$float result %v;};
    if(!RPBounded(seed,0x40800000u))return 0.0f;
    if(exponent < -126){seed=OrderedMultiply32(seed,asfloat(0x00800000u));exponent+=126;}
    seed=OrderedMultiply32(seed,RPPowerOfTwo(exponent));
    return asfloat(asuint(seed)^((asuint(a)^asuint(b))&0x80000000u));
}
float RPSqrtSeed(float positive) {
    int exponent=RPFloorExponent(positive),odd=exponent&1;
    float mantissa=RPMantissa(positive);
    if(odd!=0)mantissa=OrderedMultiply32(mantissa,2.0f);
    float seed=sqrt(mantissa);if(!RPBounded(seed,0x40800000u))return 0.0f;
    return OrderedMultiply32(seed,RPPowerOfTwo((exponent-odd)/2));
}
RetainedPair RPDivide(RetainedPair a,RetainedPair b) {
    if(a.valid==0u||b.valid==0u)return RPInvalid();
    float lower=RPAbsoluteLower(b);
    if(lower<=0.0f||RPRatioUpper(RPAbsoluteUpper(a),lower)>asfloat(0x5e000000u))return RPInvalid();
    RetainedPair q=RPConstant(RPDivideSeed(a.hi,b.hi));if(q.valid==0u)return q;
    RetainedPair residual=RPSubtract(a,RPMultiply(b,q));if(residual.valid==0u)return residual;
    if(RPRatioUpper(RPAbs(residual.hi),RPAbs(b.hi))>asfloat(0x5d000000u))return RPInvalid();
    q=RPCenter(RPAdd(q,RPConstant(RPDivideSeed(residual.hi,b.hi))));if(q.valid==0u)return q;
    residual=RPSubtract(a,RPMultiply(b,q));if(residual.valid==0u)return residual;
    return RPPack(q.hi,q.lo,RPRatioUpper(RPAbsoluteUpper(residual),lower));
}
RetainedPair RPSqrt(RetainedPair a) {
    if(a.valid==0u)return a;
    if(a.hi==0.0f&&a.lo==0.0f&&a.error==0.0f)return RPConstant(0.0f);
    if(a.hi<=0.0f||RPAbsoluteLower(a)<=0.0f)return RPInvalid();
    RetainedPair q=RPConstant(RPSqrtSeed(a.hi));if(q.valid==0u||q.hi<=0.0f||!RPBounded(q.hi,0x5d800000u))return RPInvalid();
    RetainedPair residual=RPSubtract(a,RPMultiply(q,q));if(residual.valid==0u)return residual;
    float twice=OrderedAdd32(q.hi,q.hi);
    if(RPRatioUpper(RPAbs(residual.hi),twice)>asfloat(0x5d000000u))return RPInvalid();
    q=RPCenter(RPAdd(q,RPConstant(RPDivideSeed(residual.hi,twice))));
    float lower=RPAbsoluteLower(q);if(q.valid==0u||q.hi<=0.0f||lower<=0.0f)return RPInvalid();
    residual=RPSubtract(a,RPMultiply(q,q));if(residual.valid==0u)return residual;
    return RPPack(q.hi,q.lo,RPRatioUpper(RPAbsoluteUpper(residual),lower));
}
