import retained_scratch;
import retained_scalar;
import retained_triple;
import retained_pair;
[[vk::binding(0, 0)]] StructuredBuffer<uint> inputs;
[[vk::binding(1, 0)]] RWStructuredBuffer<uint> outputs;
RetainedTriple ReadStored(uint base) {
    RetainedTriple r;
    r.hi = RPFromBits(outputs[base]);
    r.lo = RPFromBits(outputs[base + 1]);
    r.tail = RPFromBits(outputs[base + 2]);
    r.error = RPFromBits(outputs[base + 3]);
    r.valid = outputs[base + 4];
    return r;
}
void WriteStored(uint base, RetainedTriple r) {
    outputs[base] = RPBits(r.hi);
    outputs[base + 1] = RPBits(r.lo);
    outputs[base + 2] = RPBits(r.tail);
    outputs[base + 3] = RPBits(r.error);
    outputs[base + 4] = r.valid;
}
RetainedTriple ReadValue(uint index, uint rowBase) {
    RetainedTriple r;
    r.hi = RPFromBits(ScratchRead(index, 0, rowBase, 1004));
    r.lo = RPFromBits(ScratchRead(index, 1, rowBase, 1004));
    r.tail = RPFromBits(ScratchRead(index, 2, rowBase, 1004));
    r.error = RPFromBits(ScratchRead(index, 3, rowBase, 1004));
    r.valid = ScratchRead(index, 4, rowBase, 1004);
    return r;
}
void WriteValue(uint index, RetainedTriple r, uint rowBase) {
    ScratchWrite(index, 0, rowBase, 1004, RPBits(r.hi));
    ScratchWrite(index, 1, rowBase, 1004, RPBits(r.lo));
    ScratchWrite(index, 2, rowBase, 1004, RPBits(r.tail));
    ScratchWrite(index, 3, rowBase, 1004, RPBits(r.error));
    ScratchWrite(index, 4, rowBase, 1004, r.valid);
}
RetainedTriple ReadOriginal(uint index, uint row) {
    uint base = 1 + 1810 * row + 5 * index;
    RetainedTriple value = RTPack(RPFromBits(inputs[base]), RPFromBits(inputs[base + 1]),
                                  RPFromBits(inputs[base + 2]), RPFromBits(inputs[base + 3]));
    if (inputs[base + 4] != 1) value.valid = 0;
    return value;
}
bool EvaluateOne(uint index, uint row, uint rowBase, uint programBase, uint registers) {
    uint address = programBase + 203 + 5 * index, op = inputs[address],
         destination = inputs[address + 1], a = inputs[address + 2], b = inputs[address + 3],
         c = inputs[address + 4];
    if (destination >= registers || op > 11) return false;
    RetainedTriple r = RTInvalid();
    if (op == 0)
        r = RTConstant(RPFromBits(a));
    else if (op == 1) {
        if (a >= 362) return false;
        r = ReadOriginal(a, row);
    } else {
        if (a >= registers) return false;
        RetainedTriple av = ReadValue(a, rowBase);
#ifdef SIRIUS_RETAINED_PORTABLE
        if (op >= 2 && op <= 6) {
            RetainedTriple bv = RTInvalid();
            if (op != 6) {
                if (b >= registers) return false;
                bv = ReadValue(b, rowBase);
            }
            r = RTEvaluateArithmetic(av, bv, op);
        }
#else
        if (op == 6)
            r = RTSqrt(av);
#endif
        else if (op == 7)
            r = RTNegate(av);
        else if (op == 8) {
            RPScalar bound = RTAbsoluteUpper(av);
            if (av.valid != 0 && RPScalarGreater(bound, RPZero()))
                r = RTConstant(RPPowerOfTwo(RPFloorExponent(bound) + 1));
        } else {
            if (b >= registers) return false;
            RetainedTriple bv = ReadValue(b, rowBase);
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_DOPRI_PHASE)
            if (op == 2)
                r = RTAdd(av, bv);
            else if (op == 3)
                r = RTSubtract(av, bv);
            else if (op == 4)
                r = RTMultiply(av, bv);
            else if (op == 5)
                r = RTDivide(av, bv);
            else
#endif
            if (op == 9) {
                if (av.valid != 0 && bv.valid != 0)
                    r = RTConstant(RPScalarMax(RTAbsoluteUpper(av), RTAbsoluteUpper(bv)));
            } else if (op == 10 || op == 11) {
                if (c >= registers) return false;
                RetainedTriple cv = ReadValue(c, rowBase);
                if (cv.valid != 0) {
                    bool choose = op == 10 ? RPScalarGreaterEqual(cv.hi, RPZero())
                                           : (RPScalarEqual(cv.hi, RPZero()) &&
                                              RPScalarEqual(cv.lo, RPZero()) &&
                                              RPScalarEqual(cv.tail, RPZero()) &&
                                              RPScalarEqual(cv.error, RPZero()));
                    r = choose ? av : bv;
                }
            }
        }
    }
    WriteValue(destination, r, rowBase);

    return true;
}
bool EvaluateProgram(uint row, uint rowBase, uint programBase, uint count, uint registers,
                     uint required, uint lane) {
    uint needed = count;
    if (needed == 0 || needed > count) return false;
    if (lane == 0) retainedStatus = 0;
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_DOPRI_PHASE)
    GroupMemoryBarrierWithGroupSync();
#endif
    uint layerBase = programBase + 203 + 5 * count;
    if (inputs[layerBase] != 0 || inputs[layerBase + SIRIUS_RETAINED_LAYERS] != count) return false;
    uint expected = 0;
    for (uint layer = 0; layer < SIRIUS_RETAINED_LAYERS; ++layer) {
        uint start = inputs[layerBase + layer], end = inputs[layerBase + layer + 1];
        if (start >= needed) break;
        if (start != expected || end <= start || end > needed ||
            end - start > SIRIUS_RETAINED_LANES)
            RejectRetainedProgram();
#if defined(SIRIUS_RETAINED_PORTABLE) && !defined(SIRIUS_RETAINED_PARALLEL_DOPRI_PHASE)
        else {
            // Layer destinations are disjoint from every input register and
            // from one another; the generated plan releases inputs afterwards.
            for (uint node = start; node < end; ++node)
                if (!EvaluateOne(node, row, rowBase, programBase, registers)) RejectRetainedProgram();
        }
#else
        else if (start + lane < end &&
                 !EvaluateOne(start + lane, row, rowBase, programBase, registers))
            RejectRetainedProgram();
#endif
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_DOPRI_PHASE)
        GroupMemoryBarrierWithGroupSync();
#endif
        if (retainedStatus != 0) return false;
        expected = end;
    }

    for (uint i = 0; i < required; ++i) {
        uint index = inputs[programBase + 3 + i];
        if (index >= registers || ReadValue(index, rowBase).valid == 0) return false;
    }
    if (lane == 0) {
        for (uint i = 0; i < 200; ++i)
            WriteStored(rowBase + 4 + i * 5, ReadValue(inputs[programBase + 3 + i], rowBase));
    }
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_DOPRI_PHASE)
    AllMemoryBarrierWithGroupSync();
#endif
    return true;
}

[shader("compute")][numthreads(SIRIUS_RETAINED_EXECUTION_LANES, 1, 1)] void ComputeMain(uint3 id
                                                                              : SV_GroupID,
                                                                                uint lane
                                                                              : SV_GroupIndex) {
    uint ni, si, no, so;
    inputs.GetDimensions(ni, si);
    outputs.GetDimensions(no, so);
    if (ni < 1 || id.y != 0 || id.z != 0) return;
    uint rows = inputs[0], row = id.x;
    if (rows == 0 || rows > 65536 || row >= rows || ni < 1 + 1810 * rows + 203 || no % rows != 0)
        return;
    uint programBase = 1 + 1810 * rows, count = inputs[programBase],
         registers = inputs[programBase + 1];
    uint rowWords = no / rows, rowBase = row * rowWords;
    if (rowWords < 1004) return;
    for (uint i = lane; i < rowWords; i += SIRIUS_RETAINED_EXECUTION_LANES) outputs[rowBase + i] = 0;
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_DOPRI_PHASE)
    AllMemoryBarrierWithGroupSync();
#endif
    if (inputs[programBase + 2] != 200 || count > 65536 || registers != SIRIUS_RETAINED_REGISTERS ||
        ni != programBase + 203 + 5 * count + SIRIUS_RETAINED_LAYERS + 1 ||
        rowWords != 1004 + 5 * registers)
        return;
    for (uint i = 0; i < 362; ++i)
        if (ReadOriginal(i, row).valid == 0) return;
    RetainedTriple h = ReadOriginal(360, row), fraction = ReadOriginal(361, row);
    if (RPScalarLessEqual(h.hi, RPZero()) || RPScalarLessEqual(RTAbsoluteLower(h), RPZero()))
        return;
    if (!RTExact(fraction, RPZero()) && !RTExact(fraction, RPOne())) {
        RetainedTriple remaining = RTSubtract(RTConstant(RPOne()), fraction);
        if (RPScalarLessEqual(fraction.hi, RPZero()) ||
            RPScalarLessEqual(RTAbsoluteLower(fraction), RPZero()) ||
            RPScalarLessEqual(remaining.hi, RPZero()) ||
            RPScalarLessEqual(RTAbsoluteLower(remaining), RPZero()))
            return;
    }
    if (!EvaluateProgram(row, rowBase, programBase, count, registers, 200, lane)) return;
    if (lane == 0) outputs[rowBase] = 1;
}
