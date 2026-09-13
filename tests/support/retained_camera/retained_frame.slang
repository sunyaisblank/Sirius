import retained_pair;
import retained_jet;

struct RetainedVector { RetainedJet component[4]; }
RetainedVector RVAdd(RetainedVector a,RetainedVector b){
    [loop] for(int i=0;i<4;++i)a.component[i]=RJAdd(a.component[i],b.component[i]);return a;
}
RetainedVector RVScale(RetainedVector a,RetainedJet b){
    [loop] for(int i=0;i<4;++i)a.component[i]=RJMultiply(a.component[i],b);return a;
}
RetainedJet RVSpatialDot(RetainedJet metric[4][4],RetainedVector a,RetainedVector b){
    RetainedJet sum=RJFloat(0.0f);
    [loop] for(int i=1;i<4;++i)[loop] for(int j=1;j<4;++j)
        sum=RJAdd(sum,RJMultiply(RJMultiply(metric[i][j],a.component[i]),b.component[j]));return sum;
}
bool TryRetainedFrame(RetainedJet metric[4][4],RetainedJet inverse[4][4],
    float axes[3][3],float beta[3],out RetainedVector frame[4]){
    RetainedVector reference[4];RetainedJet lapse=RJDivide(RJFloat(1.0f),RJSqrt(RJNegate(inverse[0][0])));
    [loop] for(int i=0;i<4;++i)reference[0].component[i]=RJMultiply(RJNegate(lapse),inverse[i][0]);
    [loop] for(int axis=1;axis<4;++axis){
        reference[axis].component[0]=RJFloat(0.0f);
        [loop] for(int i=1;i<4;++i)reference[axis].component[i]=RJFloat(axis==2?-axes[axis-1][i-1]:axes[axis-1][i-1]);
        // u_cov=(-lapse,0,0,0); these t=0 axes already have g(axis,u)=0.
        [loop] for(int previous=1;previous<axis;++previous)
            reference[axis]=RVAdd(reference[axis],RVScale(reference[previous],RJNegate(RVSpatialDot(metric,reference[axis],reference[previous]))));
        RetainedJet norm=RJSqrt(RVSpatialDot(metric,reference[axis],reference[axis]));
        reference[axis]=RVScale(reference[axis],RJDivide(RJFloat(1.0f),norm));
    }
    RetainedPair b2=RPConstant(0.0f);
    [loop] for(int i=0;i<3;++i)b2=RPAdd(b2,RPMultiply(RPConstant(beta[i]),RPConstant(beta[i])));
    RetainedPair s=RPSqrt(RPSubtract(RPConstant(1.0f),b2));
    RetainedJet gamma=RJConstant(RPDivide(RPConstant(1.0f),s));
    RetainedJet coefficient=RJConstant(RPDivide(RPConstant(1.0f),RPMultiply(s,RPAdd(RPConstant(1.0f),s))));
    RetainedVector bv;
    [loop] for(int i=0;i<4;++i)bv.component[i]=RJFloat(0.0f);
    [loop] for(int axis=1;axis<4;++axis)bv=RVAdd(bv,RVScale(reference[axis],RJFloat(beta[axis-1])));
    frame[0]=RVScale(RVAdd(reference[0],bv),gamma);
    [loop] for(int axis=1;axis<4;++axis)frame[axis]=RVAdd(RVAdd(reference[axis],
        RVScale(bv,RJMultiply(coefficient,RJFloat(beta[axis-1])))),
        RVScale(reference[0],RJMultiply(gamma,RJFloat(beta[axis-1]))));
    [loop] for(int axis=0;axis<4;++axis)[loop] for(int i=0;i<4;++i)if(!RJValid(frame[axis].component[i]))return false;
    return true;
}
