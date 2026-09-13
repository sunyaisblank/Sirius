import retained_scratch;
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
RetainedTriple ReadValue(uint index, uint rowBase) {
    RetainedTriple r;
    r.hi = asfloat(ScratchRead(index, 0, rowBase, 204));
    r.lo = asfloat(ScratchRead(index, 1, rowBase, 204));
    r.tail = asfloat(ScratchRead(index, 2, rowBase, 204));
    r.error = asfloat(ScratchRead(index, 3, rowBase, 204));
    r.valid = ScratchRead(index, 4, rowBase, 204);
    return r;
}
void WriteValue(uint index, RetainedTriple r, uint rowBase) {
    ScratchWrite(index, 0, rowBase, 204, asuint(r.hi));
    ScratchWrite(index, 1, rowBase, 204, asuint(r.lo));
    ScratchWrite(index, 2, rowBase, 204, asuint(r.tail));
    ScratchWrite(index, 3, rowBase, 204, asuint(r.error));
    ScratchWrite(index, 4, rowBase, 204, r.valid);
}
RetainedTriple ReadOriginal(uint index,uint row) {
    uint base=1+225*row+5*index;
    RetainedTriple value=RTPack(asfloat(inputs[base]),asfloat(inputs[base+1]),
                              asfloat(inputs[base+2]),asfloat(inputs[base+3]));
    if(inputs[base+4]!=1) value.valid=0;
    return value;
}
bool EvaluateProgram(uint row,uint rowBase,uint programBase,uint count,uint registers,uint required) {
    for (uint index = 0; index < count; ++index) {
        uint address = programBase + 43 + 5 * index, op = inputs[address], destination = inputs[address + 1],
             a = inputs[address + 2], b = inputs[address + 3], c = inputs[address + 4];
        if (destination >= registers || op > 11) return false;
        RetainedTriple r = RTInvalid();
        if (op == 0)
            r = RTConstant(asfloat(a));
        else if (op == 1) {
            if (a >= 45) return false;
            r = ReadOriginal(a,row);
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
    for(uint i=0;i<40;++i)
        WriteStored(rowBase+4+i*5,ReadValue(inputs[programBase+3+i],rowBase));
    return true;
}

[shader("compute")][numthreads(SIRIUS_RETAINED_LANES,1,1)]
void ComputeMain(uint3 id:SV_DispatchThreadID) {
    uint ni,si,no,so;inputs.GetDimensions(ni,si);outputs.GetDimensions(no,so);
    if(ni<1||id.y!=0||id.z!=0) return;
    uint rows=inputs[0],row=id.x;
    if(rows==0||rows>65536||row>=rows||ni<1+225*rows+43||no%rows!=0) return;
    uint programBase=1+225*rows,count=inputs[programBase],registers=inputs[programBase+1];
    uint rowWords=no/rows,rowBase=row*rowWords;
    if(rowWords<204) return;
    for(uint i=0;i<rowWords;++i) outputs[rowBase+i]=0;
    if(inputs[programBase+2]!=40||count>65536||registers != SIRIUS_RETAINED_REGISTERS||
       ni!=programBase+43+5*count||rowWords!=204+5*registers) return;
    for(uint i=0;i<45;++i) if(ReadOriginal(i,row).valid==0) return;
    RetainedTriple chart=ReadOriginal(44,row);
    if(!RTExact(chart,1)&&!RTExact(chart,-1)) return;
    if(RTAbsoluteUpper(ReadOriginal(1,row))!=0&&RTAbsoluteUpper(ReadOriginal(3,row))!=0) return;
    if(RTExact(ReadOriginal(0,row),0)&&RTExact(ReadOriginal(2,row),0)&&RTExact(ReadOriginal(3,row),0)) {
        // Exact Minkowski lowering, including at the regular Cartesian origin.
        for(uint i=0;i<40;++i) {
            RetainedTriple value=ReadOriginal(4+i,row);
            WriteStored(rowBase+4+i*5,i%8==4?RTNegate(value):value);
        }
        outputs[rowBase]=1;return;
    }
    if(!EvaluateProgram(row,rowBase,programBase,count,registers,40)) return;
    outputs[rowBase]=1;
}
