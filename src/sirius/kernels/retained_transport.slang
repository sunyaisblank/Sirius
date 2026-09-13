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
RetainedTriple ReadOriginal(uint index, uint row) {
    uint base=1+row*230+5*index;
    RetainedTriple r=RTPack(asfloat(inputs[base]),asfloat(inputs[base+1]),asfloat(inputs[base+2]),asfloat(inputs[base+3]));
    if(inputs[base+4]!=1) r.valid=0; return r;
}
RetainedTriple ReadInput(uint index,uint row,uint rowBase) {
    return index>=4 && index<44 ? ReadStored(rowBase+804+5*(index-4)) : ReadOriginal(index,row);
}
RetainedTriple ReadValue(uint index, uint rowBase) {
    uint base = rowBase + 2404 + 5 * index;
    RetainedTriple r;
    r.hi = asfloat(outputs[base]);
    r.lo = asfloat(outputs[base + 1]);
    r.tail = asfloat(outputs[base + 2]);
    r.error = asfloat(outputs[base + 3]);
    r.valid = outputs[base + 4];
    return r;
}
void WriteValue(uint index, RetainedTriple r, uint rowBase) {
    uint base = rowBase + 2404 + 5 * index;
    outputs[base] = asuint(r.hi);
    outputs[base + 1] = asuint(r.lo);
    outputs[base + 2] = asuint(r.tail);
    outputs[base + 3] = asuint(r.error);
    outputs[base + 4] = r.valid;
}


bool EvaluateProgram(uint row,uint rowBase,uint programBase,uint count,uint registers,uint stage) {
    for (uint index = 0; index < count; ++index) {
        uint address = programBase + 43 + 5 * index, op = inputs[address], destination = inputs[address + 1],
             a = inputs[address + 2], b = inputs[address + 3], c = inputs[address + 4];
        if (destination >= registers || op > 11) return false;
        RetainedTriple r = RTInvalid();
        if (op == 0)
            r = RTConstant(asfloat(a));
        else if (op == 1) {
            if (a >= 45) return false;
            r = ReadInput(a, row, rowBase);
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


    for(uint i=0;i<40;++i) {
        uint index=inputs[programBase+3+i];
        if(index>=registers || ReadValue(index,rowBase).valid==0) return false;
    }
    for(uint i=0;i<40;++i)
        WriteStored(rowBase+1004+stage*200+i*5,ReadValue(inputs[programBase+3+i],rowBase));
    return true;
}
RetainedTriple Rational(int numerator,int denominator) {
    return RTDivide(RTConstant(float(numerator)),RTConstant(float(denominator)));
}
static const int cn[7]={0,1,3,4,8,1,1};
static const int cd[7]={1,5,10,5,9,1,1};
static const int an[7][7]={
    {0,0,0,0,0,0,0},{1,0,0,0,0,0,0},{3,9,0,0,0,0,0},
    {44,-56,32,0,0,0,0},{19372,-25360,64448,-212,0,0,0},
    {9017,-355,46732,49,-5103,0,0},{35,0,500,125,-2187,11,0}};
static const int ad[7][7]={
    {1,1,1,1,1,1,1},{5,1,1,1,1,1,1},{40,40,1,1,1,1,1},
    {45,15,9,1,1,1,1},{6561,2187,6561,729,1,1,1},
    {3168,33,5247,176,18656,1,1},{384,1,1113,192,6784,84,1}};
static const int en[7]={71,0,-71,71,-17253,22,-1};
static const int ed[7]={57600,1,16695,1920,339200,525,40};
[shader("compute")][numthreads(64,1,1)]
void ComputeMain(uint3 id : SV_DispatchThreadID) {
    uint ni,si,no,so; inputs.GetDimensions(ni,si); outputs.GetDimensions(no,so);
    if(ni<1 || id.y!=0 || id.z!=0) return;
    uint rows=inputs[0],row=id.x;
    if(rows==0 || rows>65536 || row>=rows || ni<1+230*rows+43 || no%rows!=0) return;
    uint programBase=1+230*rows,count=inputs[programBase],registers=inputs[programBase+1];
    uint rowWords=no/rows,rowBase=row*rowWords;
    if(rowWords<2404) return;
    for(uint i=0;i<rowWords;++i) outputs[rowBase+i]=0;
    if(inputs[programBase+2]!=40 || count>65536 || registers>8192 ||
       ni!=programBase+43+5*count || rowWords!=2404+5*registers) return;
    for(uint i=0;i<46;++i) if(ReadOriginal(i,row).valid==0) return;
    RetainedTriple chart=ReadOriginal(44,row);
    if((chart.hi!=1 && chart.hi!=-1) || chart.lo!=0 || chart.tail!=0 || chart.error!=0) return;
    RetainedTriple spin=ReadOriginal(1,row),lambda=ReadOriginal(3,row);
    if(RTAbsoluteUpper(spin)!=0 && RTAbsoluteUpper(lambda)!=0) return;
    RetainedTriple h=ReadOriginal(45,row);
    if(h.hi<=0 || RTAbsoluteLower(h)<=0) return;
    if(RTExact(ReadOriginal(0,row),0)&&RTExact(ReadOriginal(2,row),0)&&RTExact(lambda,0)) {
        // All seven Hamiltonian RHS stages are the same constant eta*p and
        // eta*P. The defining tableau row sums are exactly one and zero.
        for(uint i=0;i<40;++i) {
            RetainedTriple delta=RTConstant(0);
            if(i%8<4) {
                RetainedTriple slope=ReadOriginal(8+i,row);
                if(i%8==0) slope=RTNegate(slope);
                delta=RTMultiply(h,slope);
            }
            RetainedTriple full=RTAdd(ReadOriginal(4+i,row),delta);
            if(full.valid==0||delta.valid==0) return;
            WriteStored(rowBase+4+i*5,full);
            WriteStored(rowBase+204+i*5,full);
            WriteStored(rowBase+404+i*5,delta);
            WriteStored(rowBase+604+i*5,RTConstant(0));
        }
        outputs[rowBase]=1;outputs[rowBase+1]=7;return;
    }
    for(uint stage=0;stage<7;++stage) {
        for(uint i=0;i<40;++i) {
            RetainedTriple value=ReadOriginal(4+i,row);
            if(stage!=0) {
                RetainedTriple first=ReadStored(rowBase+1004+i*5);
                RetainedTriple weighted=RTMultiply(first,Rational(cn[stage],cd[stage]));
                for(uint j=1;j<stage;++j) if(an[stage][j]!=0) {
                    RetainedTriple difference=RTSubtract(ReadStored(rowBase+1004+j*200+i*5),first);
                    weighted=RTAdd(weighted,RTMultiply(difference,Rational(an[stage][j],ad[stage][j])));
                }
                RetainedTriple increment=RTMultiply(weighted,h);
                value=RTAdd(value,increment);
                if(stage==6) WriteStored(rowBase+4+80*5+i*5,increment);
            }
            if(value.valid==0) return;
            WriteStored(rowBase+804+i*5,value);
        }
        outputs[rowBase+1]=stage+1;
        if(!EvaluateProgram(row,rowBase,programBase,count,registers,stage)) return;
    }
    for(uint i=0;i<40;++i) {
        RetainedTriple first=ReadStored(rowBase+1004+i*5),weighted=RTConstant(0);
        for(uint j=1;j<7;++j) if(en[j]!=0) {
            RetainedTriple difference=RTSubtract(ReadStored(rowBase+1004+j*200+i*5),first);
            weighted=RTAdd(weighted,RTMultiply(difference,Rational(en[j],ed[j])));
        }
        RetainedTriple error=RTMultiply(weighted,h),full=ReadStored(rowBase+804+i*5);
        RetainedTriple lower=RTSubtract(full,error);
        if(error.valid==0 || lower.valid==0) return;
        WriteStored(rowBase+4+i*5,full); WriteStored(rowBase+4+40*5+i*5,lower);
        WriteStored(rowBase+4+120*5+i*5,error);
    }
    outputs[rowBase]=1;
}
