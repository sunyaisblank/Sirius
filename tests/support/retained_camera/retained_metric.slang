import retained_pair;
import retained_jet;

bool RPExactZero(RetainedPair a){return a.valid!=0u&&a.hi==0.0f&&a.lo==0.0f&&a.error==0.0f;}
bool RPPositive(RetainedPair a){return a.valid!=0u&&a.hi>0.0f&&RPAbsoluteLower(a)>0.0f;}

// The positive Kerr root and its coordinate gradient use the same retained
// scaled quantities. The selected positive scale is a numerical constant.
bool TryRetainedKerrRadius(RetainedPair position[4],float spin,out RetainedJet radius){
    radius=RJConstant(RPInvalid());RetainedPair a=RPConstant(spin);
    float bound=RPAbs(spin);
    [loop] for(int i=1;i<4;++i){if(position[i].valid==0u)return false;bound=max(bound,RPAbsoluteUpper(position[i]));}
    if(bound<=0.0f||!RPBounded(bound,0x7b800000u))return false;
    RetainedPair scale=RPConstant(RPPowerOfTwo(RPFloorExponent(bound)+1));
    RetainedPair x=RPDivide(position[1],scale),y=RPDivide(position[2],scale),z=RPDivide(position[3],scale),s=RPDivide(a,scale);
    RetainedPair a2=RPMultiply(s,s),z2=RPMultiply(z,z);
    RetainedPair R2=RPAdd(RPAdd(RPMultiply(x,x),RPMultiply(y,y)),z2);
    RetainedPair reduced=RPSubtract(R2,a2);
    RetainedPair root=RPSqrt(RPAdd(RPMultiply(reduced,reduced),RPMultiply(RPConstant(4.0f),RPMultiply(a2,z2))));
    if(!RPPositive(root))return false;
    RetainedPair r2;
    if(reduced.hi>=0.0f)r2=RPMultiply(RPConstant(0.5f),RPAdd(reduced,root));
    else r2=RPDivide(RPMultiply(RPConstant(2.0f),RPMultiply(a2,z2)),RPSubtract(root,reduced));
    RetainedPair r=RPSqrt(r2);if(!RPPositive(r))return false;
    radius.value=RPMultiply(scale,r);radius.derivative[0]=RPConstant(0.0f);
    radius.derivative[1]=RPDivide(RPMultiply(x,r),root);
    radius.derivative[2]=RPDivide(RPMultiply(y,r),root);
    radius.derivative[3]=RPDivide(RPMultiply(z,RPAdd(r2,a2)),RPMultiply(r,root));
    return RJValid(radius)&&RPPositive(radius.value);
}

// g and inverse retain coordinate derivatives, before contraction with camera X.
bool TryRetainedMetric(RetainedPair position[4],float M,float a,float Q,float Lambda,
    out RetainedJet g[4][4],out RetainedJet inverse[4][4]){
    [loop] for(int i=0;i<4;++i){
        if(position[i].valid==0u)return false;
        [loop] for(int j=0;j<4;++j){g[i][j]=RJFloat(i==j?(i==0?-1.0f:1.0f):0.0f);inverse[i][j]=g[i][j];}
    }
    if(a!=0.0f&&Lambda!=0.0f)return false;
    if(M==0.0f&&a==0.0f&&Q==0.0f&&Lambda==0.0f)return true;
    // Pure de Sitter: g-eta=O(|x|^2), hence the first metric jet is exactly flat
    // at the origin. This does not assert a second derivative there.
    if(M==0.0f&&a==0.0f&&Q==0.0f&&RPExactZero(position[1])&&RPExactZero(position[2])&&RPExactZero(position[3]))return true;
    RetainedJet r;if(!TryRetainedKerrRadius(position,a,r))return false;
    RetainedJet x=RJCoordinate(position[1],1),y=RJCoordinate(position[2],2),z=RJCoordinate(position[3],3);
    RetainedJet spin=RJFloat(a),a2=RJMultiply(spin,spin),r2=RJMultiply(r,r);
    RetainedJet denominator=RJAdd(r2,a2),ell[4];ell[0]=RJFloat(1.0f);
    ell[1]=RJDivide(RJAdd(RJMultiply(r,x),RJMultiply(spin,y)),denominator);
    ell[2]=RJDivide(RJSubtract(RJMultiply(r,y),RJMultiply(spin,x)),denominator);
    ell[3]=RJDivide(z,r);
    RetainedJet H=RJFloat(0.0f);
    if(M!=0.0f||Q!=0.0f){
        RetainedJet sigma=RJAdd(r2,RJMultiply(a2,RJMultiply(ell[3],ell[3])));
        H=RJDivide(RJSubtract(RJMultiply(RJFloat(2.0f),RJMultiply(RJFloat(M),r)),RJMultiply(RJFloat(Q),RJFloat(Q))),sigma);
    }
    if(a==0.0f&&Lambda!=0.0f)H=RJAdd(H,RJDivide(RJMultiply(RJFloat(Lambda),r2),RJFloat(3.0f)));
    [loop] for(int i=0;i<4;++i)[loop] for(int j=0;j<4;++j){
        RetainedJet correction=RJMultiply(RJMultiply(H,ell[i]),ell[j]);
        g[i][j]=RJAdd(g[i][j],correction);
        if((i==0)!=(j==0))correction=RJNegate(correction);
        inverse[i][j]=RJSubtract(inverse[i][j],correction);
        if(!RJValid(g[i][j])||!RJValid(inverse[i][j]))return false;
    }
    return true;
}
