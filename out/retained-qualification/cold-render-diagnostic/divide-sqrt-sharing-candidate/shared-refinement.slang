#ifdef SIRIUS_RETAINED_PORTABLE
// Preserve each operation's guards and ordered correction while sharing the
// two residual evaluations in the portable camera interpreter.
RetainedPair RPEvaluateArithmetic(RetainedPair a,RetainedPair b,uint op) {
    bool square_root=op==6;
    RetainedPair q=RPInvalid();
    RPScalar lower=RPZero();
    if(square_root) {
        if(a.valid==0u)return a;
        if(RPScalarEqual(a.hi,RPZero())&&RPScalarEqual(a.lo,RPZero())&&RPScalarEqual(a.error,RPZero()))return RPConstant(RPZero());
        if(RPScalarLessEqual(a.hi,RPZero())||RPScalarLessEqual(RPAbsoluteLower(a),RPZero()))return RPInvalid();
        q=RPConstant(RPSqrtSeed(a.hi));
        if(q.valid==0u||RPScalarLessEqual(q.hi,RPZero())||!RPBounded(q.hi,0x5d800000u))return RPInvalid();
    } else if(op==5) {
        if(a.valid==0u||b.valid==0u)return RPInvalid();
        lower=RPAbsoluteLower(b);
        if(RPScalarLessEqual(lower,RPZero())||RPScalarGreater(RPRatioUpper(RPAbsoluteUpper(a),lower),RPFromBits(0x5e000000u)))return RPInvalid();
        q=RPConstant(RPDivideSeed(a.hi,b.hi));
        if(q.valid==0u)return q;
    }
    RetainedPair residual=RPInvalid();
    RPScalar denominator=RPZero();
    [loop] for(uint phase=0;phase<3;++phase) {
        RetainedPair left=a,right=b;
        if(op==3)right=RPNegate(b);
        if(op>=4) {
            if(phase==1) {
                left=q;
                right=RPConstant(RPDivideSeed(residual.hi,denominator));
            } else {
                RetainedPair product=RPMultiply(op==4?a:(square_root?q:b),op==4?b:q);
                if(op==4)return product;
                right=RPNegate(product);
            }
        }
        RetainedPair sum=RPAdd(left,right);
        if(op==2||op==3)return sum;
        if(phase==1) {
            q=RPCenter(sum);
            if(square_root) {
                lower=RPAbsoluteLower(q);
                if(q.valid==0u||RPScalarLessEqual(q.hi,RPZero())||RPScalarLessEqual(lower,RPZero()))return RPInvalid();
            } else if(q.valid==0u)return q;
        } else {
            residual=sum;
            if(residual.valid==0u)return residual;
            if(phase==0) {
                denominator=b.hi;
                if(square_root)denominator=OrderedAdd32(q.hi,q.hi);
                if(RPScalarGreater(RPRatioUpper(RPAbs(residual.hi),square_root?denominator:RPAbs(b.hi)),RPFromBits(0x5d000000u)))return RPInvalid();
            }
        }
    }
    return RPPackCanonicalCenter(q,RPRatioUpper(RPAbsoluteUpper(residual),lower));
}
#endif
