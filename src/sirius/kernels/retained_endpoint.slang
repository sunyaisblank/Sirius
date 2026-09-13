import retained_triple;
import retained_pair;
[[vk::binding(0, 0)]] StructuredBuffer<uint> inputs;
[[vk::binding(1, 0)]] RWStructuredBuffer<uint> outputs;
RetainedTriple ReadStored(uint base) {
    RetainedTriple r; r.hi=asfloat(outputs[base]); r.lo=asfloat(outputs[base+1]);
    r.tail=asfloat(outputs[base+2]); r.error=asfloat(outputs[base+3]); r.valid=outputs[base+4]; return r;
}
void WriteStored(uint base, RetainedTriple r) {
    outputs[base]=asuint(r.hi); outputs[base+1]=asuint(r.lo);
    outputs[base+2]=asuint(r.tail); outputs[base+3]=asuint(r.error); outputs[base+4]=r.valid;
}
RetainedTriple ReadValue(uint index,uint rowBase) {return ReadStored(rowBase+769+5*index);}
void WriteValue(uint index,RetainedTriple value,uint rowBase) {WriteStored(rowBase+769+5*index,value);}
bool EvaluateProgram(uint row,uint rowBase,uint programBase,uint count,uint registers,uint required) {
    for (uint index = 0; index < count; ++index) {
        uint address = programBase + 103 + 5 * index, op = inputs[address], destination = inputs[address + 1],
             a = inputs[address + 2], b = inputs[address + 3], c = inputs[address + 4];
        if (destination >= registers || op > 11) return false;
        RetainedTriple r = RTInvalid();
        if (op == 0)
            r = RTConstant(asfloat(a));
        else if (op == 1) {
            if (a >= 53) return false;
            r = ReadStored(rowBase + 504 + 5*a);
        } else {
            if (a >= registers) return false;
            RetainedTriple av = ReadValue(a, rowBase);
            if (op == 6)
                r = RTSqrt(av);
            else if (op == 7)
                r = RTNegate(av);
            else if (op == 8) {
                float bound = RTAbsoluteUpper(av);
                if (av.valid != 0 && bound > 0)
                    r = RTConstant(RPPowerOfTwo(RPFloorExponent(bound) + 1));
            } else {
                if (b >= registers) return false;
                RetainedTriple bv = ReadValue(b, rowBase);
                if (op == 2)
                    r = RTAdd(av, bv);
                else if (op == 3)
                    r = RTSubtract(av, bv);
                else if (op == 4)
                    r = RTMultiply(av, bv);
                else if (op == 5)
                    r = RTDivide(av, bv);
                else if (op == 9) {
                    if (av.valid != 0 && bv.valid != 0)
                        r = RTConstant(max(RTAbsoluteUpper(av), RTAbsoluteUpper(bv)));
                } else if (op == 10 || op == 11) {
                    if (c >= registers) return false;
                    RetainedTriple cv = ReadValue(c, rowBase);
                    if (cv.valid != 0) {
                        bool choose =
                            op == 10 ? cv.hi >= 0 : (cv.hi == 0 && cv.lo == 0 && cv.tail == 0 && cv.error == 0);
                        r = choose ? av : bv;
                    }
                }
            }
        }
        WriteValue(destination, r, rowBase);
    }


    for(uint i=0;i<required;++i) {
        uint index=inputs[programBase+3+i];
        if(index>=registers || ReadValue(index,rowBase).valid==0) return false;
    }
    for(uint i=0;i<100;++i)
        WriteStored(rowBase+4+i*5,ReadValue(inputs[programBase+3+i],rowBase));
    return true;
}

RetainedTriple Metric(uint rowBase,uint i,uint j) {return ReadStored(rowBase+4+5*(4*i+j));}
RetainedTriple Tangent(uint rowBase,uint i) {return ReadStored(rowBase+4+5*(16+i));}
bool Project(uint rowBase,out uint component,out RetainedTriple root) {
    component=4;root=RTInvalid();RetainedTriple best=RTInvalid();
    RetainedTriple defect=RTConstant(0);float scale=0;
    for(uint i=0;i<4;++i) for(uint j=0;j<4;++j) {
        RetainedTriple term=RTMultiply(Metric(rowBase,i,j),RTMultiply(Tangent(rowBase,i),Tangent(rowBase,j)));
        defect=RTAdd(defect,term);scale=RPUpAdd(scale,RTAbsoluteUpper(term));
    }
    if(defect.valid==0 || scale<=0 || RPRatioUpper(RTAbsoluteUpper(defect),scale)>1e-6f) return false;
    // Prefer a temporal root; if none is represented, preserve time orientation
    // and choose the smallest normalized spatial correction.
    for(uint c=0;c<4;++c) {
        if(c>0 && component==0) break;
        RetainedTriple a=Metric(rowBase,c,c),b=RTConstant(0),d=RTConstant(0);
        for(uint i=0;i<4;++i) if(i!=c) {
            b=RTAdd(b,RTMultiply(RTAdd(Metric(rowBase,c,i),Metric(rowBase,i,c)),Tangent(rowBase,i)));
            for(uint j=0;j<4;++j) if(j!=c)
                d=RTAdd(d,RTMultiply(Metric(rowBase,i,j),RTMultiply(Tangent(rowBase,i),Tangent(rowBase,j))));
        }
        RetainedTriple first=RTInvalid(),second=RTInvalid();
        if(RTExact(a,0)) first=RTDivide(RTNegate(d),b);
        else {
            RetainedTriple discriminant=RTSubtract(RTMultiply(b,b),RTMultiply(RTConstant(4),RTMultiply(a,d)));
            RetainedTriple square=RTSqrt(discriminant);
            RetainedTriple q=RTMultiply(RTConstant(-0.5f),RTAdd(b,b.hi>=0?square:RTNegate(square)));
            if(RTExact(q,0)) first=RTDivide(RTNegate(b),RTMultiply(RTConstant(2),a));
            else {first=RTDivide(q,a);second=RTDivide(d,q);}
        }
        for(uint r=0;r<2;++r) {
            RetainedTriple candidate=r==0?first:second;
            if(candidate.valid==0) continue;
            RetainedTriple denominator=RTConstant(0);float denominatorScale=0;
            for(uint i=0;i<4;++i) {
                RetainedTriple term=RTMultiply(Metric(rowBase,c,i),i==c?candidate:Tangent(rowBase,i));
                denominator=RTAdd(denominator,term);
                denominatorScale=RPUpAdd(denominatorScale,RTAbsoluteUpper(term));
            }
            // A central root with a singular differentiated constraint cannot
            // own the coupled candidate. Try the represented spatial roots.
            if(denominator.valid==0 || RTAbsoluteLower(denominator)<=
                RPUpMultiply(5.6843418860808015e-14f,denominatorScale)) continue;
            RetainedTriple difference=RTSubtract(candidate,Tangent(rowBase,c));
            RetainedTriple original=RTCenter(Tangent(rowBase,c));
            RetainedTriple absolute=RTCenter(difference.hi<0?RTNegate(difference):difference);
            RetainedTriple scaleValue=RTAdd(RTConstant(1),original.hi<0?RTNegate(original):original);
            RetainedTriple correction=RTDivide(absolute,scaleValue);
            if(correction.valid!=0 && (best.valid==0 || RTSubtract(correction,best).hi<0)) {
                best=correction;component=c;root=candidate;
            }
        }
    }
    return component<4 && root.valid!=0;
}

[shader("compute")][numthreads(64,1,1)]
void ComputeMain(uint3 id:SV_DispatchThreadID) {
    uint ni,si,no,so;inputs.GetDimensions(ni,si);outputs.GetDimensions(no,so);
    if(ni<1||id.y!=0||id.z!=0) return;
    uint rows=inputs[0],row=id.x;
    if(rows==0||rows>65536||row>=rows||ni<1+225*rows+103||no%rows!=0) return;
    uint programBase=1+225*rows,count=inputs[programBase],registers=inputs[programBase+1];
    uint rowWords=no/rows,rowBase=row*rowWords;
    if(rowWords<769) return;
    for(uint i=0;i<rowWords;++i) outputs[rowBase+i]=0;
    if(inputs[programBase+2]!=100||count>65536||registers>8192||
       ni!=programBase+103+5*count||rowWords!=769+5*registers) return;
    for(uint i=0;i<45;++i) {
        uint address=1+225*row+5*i;
        RetainedTriple value=RTPack(asfloat(inputs[address]),asfloat(inputs[address+1]),
                                   asfloat(inputs[address+2]),asfloat(inputs[address+3]));
        if(inputs[address+4]!=1||value.valid==0) return;
        WriteStored(rowBase+504+5*i,value);
    }
    RetainedTriple chart=ReadStored(rowBase+504+5*44);
    if(!RTExact(chart,1)&&!RTExact(chart,-1)) return;
    if(RTAbsoluteUpper(ReadStored(rowBase+504+5))!=0&&
       RTAbsoluteUpper(ReadStored(rowBase+504+15))!=0) return;
    for(uint i=45;i<53;++i) WriteStored(rowBase+504+5*i,RTConstant(i==49?1:0));
    if(!EvaluateProgram(row,rowBase,programBase,count,registers,20)) return;
    uint component;RetainedTriple root;
    if(!Project(rowBase,component,root)) return;
    for(uint i=0;i<4;++i) {
        WriteStored(rowBase+504+5*(45+i),i==component?root:Tangent(rowBase,i));
        WriteStored(rowBase+504+5*(49+i),RTConstant(i==component?1:0));
    }
    if(!EvaluateProgram(row,rowBase,programBase,count,registers,100)) return;
    outputs[rowBase]=1;outputs[rowBase+1]=component;
}
