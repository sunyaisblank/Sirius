import retained_pair;
import retained_jet;
import retained_metric;
import retained_frame;
import retained_camera_direction;

RetainedJet RJContract(RetainedJet coordinate,RetainedPair X[4][4]){
    RetainedJet result=RJConstant(coordinate.value);
    [loop] for(int c=0;c<4;++c)[loop] for(int mu=0;mu<4;++mu)
        result.derivative[c]=RPAdd(result.derivative[c],RPMultiply(coordinate.derivative[mu],X[c][mu]));
    return result;
}
bool TryRetainedCameraLaunch(float input[32],out RetainedPair scientific[104]){
    [loop] for(int i=0;i<104;++i)scientific[i]=RPInvalid();
    float axes[3][3],beta[3];
    [loop] for(int axis=0;axis<3;++axis){beta[axis]=input[17+axis];
        [loop] for(int i=0;i<3;++i)axes[axis][i]=input[8+3*axis+i];}
    bool thin=input[29]==1.0f;
    RetainedPair position[4],X[4][4];
    [loop] for(int mu=0;mu<4;++mu){position[mu]=RPConstant(input[4+mu]);
        [loop] for(int c=0;c<4;++c)X[c][mu]=RPConstant(0.0f);}
    RetainedJet coordinateMetric[4][4],coordinateInverse[4][4];
    RetainedVector frame[4];
    // One authority builds the fixed central frame, then the displaced frame.
    [loop] for(int phase=0;phase<2;++phase){
        if(!TryRetainedMetric(position,input[0],input[1],input[2],input[3],coordinateMetric,coordinateInverse))return false;
        RetainedJet metric[4][4],inverse[4][4];
        [loop] for(int i=0;i<4;++i)[loop] for(int j=0;j<4;++j){
            metric[i][j]=RJContract(coordinateMetric[i][j],X);
            inverse[i][j]=RJContract(coordinateInverse[i][j],X);
        }
        if(!TryRetainedFrame(metric,inverse,axes,beta,frame))return false;
        if(phase==0&&thin){
            [loop] for(int mu=0;mu<4;++mu){
                X[2][mu]=frame[3].component[mu].value;X[3][mu]=frame[2].component[mu].value;
                position[mu]=RPAdd(RPAdd(position[mu],RPMultiply(X[3][mu],RPConstant(input[28]))),RPMultiply(X[2][mu],RPConstant(input[27])));
            }
        }
    }
    RetainedPair n[3],dn[4][3],q[3],norm;
    if(!TryRetainedCameraDirection(input[20],input[21],input[22],input[23],input[24],input[25],input[26],input[27],input[28],thin,n,dn,q,norm))return false;
    RetainedVector ray=RVScale(frame[0],RJFloat(-1.0f));
    [loop] for(int i=0;i<3;++i){
        RetainedJet direction=RJConstant(n[i]);
        [loop] for(int c=0;c<4;++c)direction.derivative[c]=dn[c][i];
        ray=RVAdd(ray,RVScale(frame[i+1],direction));
    }
    RetainedPair Gamma[4][4][4];
    [loop] for(int mu=0;mu<4;++mu)[loop] for(int alpha=0;alpha<4;++alpha)[loop] for(int b=alpha;b<4;++b){
        RetainedPair sum=RPConstant(0.0f);
        [loop] for(int s=0;s<4;++s){
            RetainedPair bracket=RPSubtract(RPAdd(coordinateMetric[s][b].derivative[alpha],coordinateMetric[s][alpha].derivative[b]),coordinateMetric[alpha][b].derivative[s]);
            sum=RPAdd(sum,RPMultiply(coordinateInverse[mu][s].value,bracket));
        }
        Gamma[mu][alpha][b]=RPMultiply(RPConstant(0.5f),sum);Gamma[mu][b][alpha]=Gamma[mu][alpha][b];
    }
    [loop] for(int mu=0;mu<4;++mu){
        scientific[mu]=position[mu];scientific[4+mu]=ray.component[mu].value;
        [loop] for(int axis=0;axis<4;++axis){
            scientific[8+4*axis+mu]=frame[axis].component[mu].value;
            scientific[88+4*axis+mu]=scientific[8+4*axis+mu];
        }
        [loop] for(int c=0;c<4;++c){
            RetainedPair K=ray.component[mu].derivative[c],V=K;
            [loop] for(int alpha=0;alpha<4;++alpha)[loop] for(int b=0;b<4;++b)
                V=RPAdd(V,RPMultiply(RPMultiply(Gamma[mu][alpha][b],ray.component[alpha].value),X[c][b]));
            scientific[24+4*c+mu]=X[c][mu];scientific[40+4*c+mu]=K;
            scientific[56+4*c+mu]=V;scientific[72+4*c+mu]=frame[0].component[mu].derivative[c];
        }
    }
    [loop] for(int i=0;i<104;++i)if(scientific[i].valid==0u)return false;
    return true;
}
