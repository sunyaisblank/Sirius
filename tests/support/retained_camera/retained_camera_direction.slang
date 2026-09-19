import retained_pair;

// Same smooth film/pupil coordinates and units as the complete camera launch.
// Keep the pair and its radius intact when the differentiated frame consumes it.
bool TryRetainedCameraDirection(float width,float height,float aspect,float tanHalfFov,
    float focus,float filmX,float filmY,float pupilRight,float pupilUp,bool thin,
    out RetainedPair direction[3],out RetainedPair derivative[4][3],
    out RetainedPair q[3],out RetainedPair norm) {
    RetainedPair zero=RPConstant(0.0f),one=RPConstant(1.0f),two=RPConstant(2.0f);
    RetainedPair D=RPConstant(thin?focus:1.0f),T=RPConstant(tanHalfFov);
    RetainedPair DT=RPMultiply(D,T);
    RetainedPair fx=RPMultiply(RPSubtract(RPDivide(RPMultiply(two,RPConstant(filmX)),RPConstant(width)),one),RPConstant(aspect));
    RetainedPair fy=RPSubtract(one,RPDivide(RPMultiply(two,RPConstant(filmY)),RPConstant(height)));
    q[0]=D;q[1]=RPSubtract(RPMultiply(DT,fy),RPConstant(thin?pupilUp:0.0f));
    q[2]=RPSubtract(RPMultiply(DT,fx),RPConstant(thin?pupilRight:0.0f));
    RetainedPair square=zero;
    for(int i=0;i<3;++i)square=RPAdd(square,RPMultiply(q[i],q[i]));
    norm=RPSqrt(square);
    for(int i=0;i<3;++i)direction[i]=RPDivide(q[i],norm);
    RetainedPair dq[4][3];
    for(int c=0;c<4;++c)for(int i=0;i<3;++i)dq[c][i]=zero;
    dq[0][2]=RPDivide(RPMultiply(RPMultiply(two,DT),RPConstant(aspect)),RPConstant(width));
    dq[1][1]=RPNegate(RPDivide(RPMultiply(two,DT),RPConstant(height)));
    if(thin){dq[2][2]=RPNegate(one);dq[3][1]=RPNegate(one);}
    for(int c=0;c<4;++c){
        RetainedPair projection=zero;
        for(int i=0;i<3;++i)projection=RPAdd(projection,RPMultiply(direction[i],dq[c][i]));
        for(int i=0;i<3;++i)derivative[c][i]=RPDivide(
            RPSubtract(dq[c][i],RPMultiply(direction[i],projection)),norm);
    }
    if(norm.valid==0u)return false;
    for(int i=0;i<3;++i)if(q[i].valid==0u||direction[i].valid==0u)return false;
    for(int c=0;c<4;++c)for(int i=0;i<3;++i)if(derivative[c][i].valid==0u)return false;
    return true;
}
