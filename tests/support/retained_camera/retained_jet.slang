import retained_pair;

struct RetainedJet { RetainedPair value; RetainedPair derivative[4]; }
RetainedJet RJConstant(RetainedPair value) {
    RetainedJet r;r.value=value;
    [loop] for(int c=0;c<4;++c)r.derivative[c]=RPConstant(0.0f);
    return r;
}
RetainedJet RJFloat(float value){return RJConstant(RPConstant(value));}
RetainedJet RJCoordinate(RetainedPair value,int coordinate){
    RetainedJet r=RJConstant(value);r.derivative[coordinate]=RPConstant(1.0f);return r;
}
RetainedJet RJAdd(RetainedJet a,RetainedJet b){
    RetainedJet r;r.value=RPAdd(a.value,b.value);
    [loop] for(int c=0;c<4;++c)r.derivative[c]=RPAdd(a.derivative[c],b.derivative[c]);return r;
}
RetainedJet RJNegate(RetainedJet a){
    a.value=RPNegate(a.value);
    [loop] for(int c=0;c<4;++c)a.derivative[c]=RPNegate(a.derivative[c]);return a;
}
RetainedJet RJSubtract(RetainedJet a,RetainedJet b){return RJAdd(a,RJNegate(b));}
RetainedJet RJMultiply(RetainedJet a,RetainedJet b){
    RetainedJet r;r.value=RPMultiply(a.value,b.value);
    [loop] for(int c=0;c<4;++c)r.derivative[c]=RPAdd(RPMultiply(a.derivative[c],b.value),RPMultiply(a.value,b.derivative[c]));return r;
}
RetainedJet RJDivide(RetainedJet a,RetainedJet b){
    RetainedJet r;r.value=RPDivide(a.value,b.value);
    [loop] for(int c=0;c<4;++c)r.derivative[c]=RPDivide(RPSubtract(a.derivative[c],RPMultiply(r.value,b.derivative[c])),b.value);return r;
}
RetainedJet RJSqrt(RetainedJet a){
    RetainedJet r;r.value=RPSqrt(a.value);RetainedPair denominator=RPMultiply(RPConstant(2.0f),r.value);
    [loop] for(int c=0;c<4;++c)r.derivative[c]=RPDivide(a.derivative[c],denominator);return r;
}
bool RJValid(RetainedJet a){
    if(a.value.valid==0u)return false;
    [loop] for(int c=0;c<4;++c)if(a.derivative[c].valid==0u)return false;return true;
}
