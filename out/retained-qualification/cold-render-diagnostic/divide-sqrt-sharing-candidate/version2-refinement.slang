#ifdef SIRIUS_RETAINED_PORTABLE
// Preserve each operation's guards and ordered correction while sharing the
// two residual evaluations in the portable camera interpreter.
RetainedPair RPProductOrRefinement(RetainedPair a,RetainedPair b,uint op) {
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
    [loop] for(uint i=0;i<2;++i) {
        RetainedPair product=RPMultiply(op==4?a:(square_root?q:b),op==4?b:q);
        if(op==4)return product;
        residual=RPSubtract(a,product);
        if(residual.valid==0u)return residual;
        if(i==0) {
            RPScalar denominator=b.hi;
            if(square_root)denominator=OrderedAdd32(q.hi,q.hi);
            if(RPScalarGreater(RPRatioUpper(RPAbs(residual.hi),square_root?denominator:RPAbs(b.hi)),RPFromBits(0x5d000000u)))return RPInvalid();
            q=RPCenter(RPAdd(q,RPConstant(RPDivideSeed(residual.hi,denominator))));
            if(square_root) {
                lower=RPAbsoluteLower(q);
                if(q.valid==0u||RPScalarLessEqual(q.hi,RPZero())||RPScalarLessEqual(lower,RPZero()))return RPInvalid();
            } else if(q.valid==0u)return q;
        }
    }
    return RPPackCanonicalCenter(q,RPRatioUpper(RPAbsoluteUpper(residual),lower));
}
#endif
