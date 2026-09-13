import retained_scratch;
import retained_triple;
import retained_pair;
[[vk::binding(0, 0)]] StructuredBuffer<uint> inputs;
[[vk::binding(1, 0)]] RWStructuredBuffer<uint> outputs;
RetainedTriple ReadStored(uint base) {
    RetainedTriple r;
    r.hi = asfloat(outputs[base]);
    r.lo = asfloat(outputs[base + 1]);
    r.tail = asfloat(outputs[base + 2]);
    r.error = asfloat(outputs[base + 3]);
    r.valid = outputs[base + 4];
    return r;
}
void WriteStored(uint base, RetainedTriple r) {
    outputs[base] = asuint(r.hi);
    outputs[base + 1] = asuint(r.lo);
    outputs[base + 2] = asuint(r.tail);
    outputs[base + 3] = asuint(r.error);
    outputs[base + 4] = r.valid;
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
RetainedTriple ReadOriginal(uint index, uint row) {
    uint base = 1 + 560 * row + 5 * index;
    RetainedTriple value = RTPack(asfloat(inputs[base]), asfloat(inputs[base + 1]),
                                  asfloat(inputs[base + 2]), asfloat(inputs[base + 3]));
    if (inputs[base + 4] != 1) value.valid = 0;
    return value;
}
bool EvaluateOne(uint index, uint row, uint rowBase, uint programBase, uint registers) {
    uint address = programBase + 43 + 5 * index, op = inputs[address],
         destination = inputs[address + 1], a = inputs[address + 2], b = inputs[address + 3],
         c = inputs[address + 4];
    if (destination >= registers || op > 11) return false;
    RetainedTriple r = RTInvalid();
    if (op == 0)
        r = RTConstant(asfloat(a));
    else if (op == 1) {
        if (a >= 112) return false;
        r = ReadOriginal(a, row);
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
                    bool choose = op == 10
                                      ? cv.hi >= 0
                                      : (cv.hi == 0 && cv.lo == 0 && cv.tail == 0 && cv.error == 0);
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
    GroupMemoryBarrierWithGroupSync();
    uint layerBase = programBase + 43 + 5 * count;
    if (inputs[layerBase] != 0 || inputs[layerBase + SIRIUS_RETAINED_LAYERS] != count) return false;
    uint expected = 0;
    for (uint layer = 0; layer < SIRIUS_RETAINED_LAYERS; ++layer) {
        uint start = inputs[layerBase + layer], end = inputs[layerBase + layer + 1];
        if (start >= needed) break;
        if (start != expected || end <= start || end > needed ||
            end - start > SIRIUS_RETAINED_LANES)
            RejectRetainedProgram();
        else if (start + lane < end &&
                 !EvaluateOne(start + lane, row, rowBase, programBase, registers))
            RejectRetainedProgram();
        GroupMemoryBarrierWithGroupSync();
        if (retainedStatus != 0) return false;
        expected = end;
    }

    for (uint i = 0; i < required; ++i) {
        uint index = inputs[programBase + 3 + i];
        if (index >= registers || ReadValue(index, rowBase).valid == 0) return false;
    }
    if (lane == 0) {
        for (uint i = 0; i < 40; ++i)
            WriteStored(rowBase + 4 + i * 5, ReadValue(inputs[programBase + 3 + i], rowBase));
    }
    AllMemoryBarrierWithGroupSync();
    return true;
}

bool FlatDense(uint row, uint rowBase, RetainedTriple h, RetainedTriple s, RetainedTriple moving) {
    // The same factored Hermite polynomial with exactly zero connection.
    // Arbitrary endpoint slopes are retained; flat geometry does not imply
    // that a supplied test segment is linear.
    for (uint group = 0; group < 5; ++group)
        for (uint axis = 0; axis < 4; ++axis) {
            uint index = group * 8 + axis;
            RetainedTriple first = ReadOriginal(4 + index, row),
                           slope0 = ReadOriginal(8 + index, row),
                           slope1 = ReadOriginal(48 + index, row),
                           delta = ReadOriginal(84 + group * 4 + axis, row);
            RetainedTriple a = RTSubtract(RTDivide(delta, h), slope0),
                           b = RTSubtract(slope1, slope0);
            RetainedTriple quadratic = RTSubtract(RTMultiply(RTConstant(3), a), b);
            RetainedTriple cubic = RTSubtract(b, RTMultiply(RTConstant(2), a));
            RetainedTriple value = RTAdd(
                first,
                RTMultiply(RTAdd(slope0, RTMultiply(RTAdd(quadratic, RTMultiply(cubic, s)), s)),
                           RTMultiply(h, s)));
            RetainedTriple slope =
                RTAdd(slope0, RTMultiply(RTAdd(RTMultiply(RTConstant(2), quadratic),
                                               RTMultiply(RTMultiply(RTConstant(3), cubic), s)),
                                         s));
            if (RTExact(s, 0)) {
                value = first;
                slope = slope0;
            }
            if (RTExact(s, 1)) {
                value = RTAdd(first, delta);
                slope = slope1;
            }
            if (value.valid == 0 || slope.valid == 0) return false;
            WriteStored(rowBase + 4 + index * 5, value);
            WriteStored(rowBase + 4 + (index + 4) * 5, slope);
        }
    if (RTExact(moving, 1)) {
        RetainedTriple denominator = RTConstant(0);
        for (uint axis = 0; axis < 4; ++axis)
            denominator = RTAdd(denominator, RTMultiply(ReadOriginal(106 + axis, row),
                                                        ReadStored(rowBase + 4 + (4 + axis) * 5)));
        for (uint column = 0; column < 4; ++column) {
            RetainedTriple numerator = RTConstant(0);
            for (uint axis = 0; axis < 4; ++axis)
                numerator = RTAdd(
                    numerator, RTMultiply(ReadOriginal(106 + axis, row),
                                          ReadStored(rowBase + 4 + (8 + 8 * column + axis) * 5)));
            RetainedTriple shift = RTDivide(RTNegate(numerator), denominator);
            for (uint axis = 0; axis < 4; ++axis) {
                uint address = rowBase + 4 + (8 + 8 * column + axis) * 5;
                RetainedTriple value =
                    RTAdd(ReadStored(address),
                          RTMultiply(ReadStored(rowBase + 4 + (4 + axis) * 5), shift));
                if (value.valid == 0) return false;
                WriteStored(address, value);
            }
        }
    }
    return true;
}

[shader("compute")][numthreads(SIRIUS_RETAINED_LANES, 1, 1)] void ComputeMain(uint3 id
                                                                              : SV_GroupID,
                                                                                uint lane
                                                                              : SV_GroupIndex) {
    uint ni, si, no, so;
    inputs.GetDimensions(ni, si);
    outputs.GetDimensions(no, so);
    if (ni < 1 || id.y != 0 || id.z != 0) return;
    uint rows = inputs[0], row = id.x;
    if (rows == 0 || rows > 65536 || row >= rows || ni < 1 + 560 * rows + 43 || no % rows != 0)
        return;
    uint programBase = 1 + 560 * rows, count = inputs[programBase],
         registers = inputs[programBase + 1];
    uint rowWords = no / rows, rowBase = row * rowWords;
    if (rowWords < 204) return;
    for (uint i = lane; i < rowWords; i += SIRIUS_RETAINED_LANES) outputs[rowBase + i] = 0;
    AllMemoryBarrierWithGroupSync();
    if (inputs[programBase + 2] != 40 || count > 65536 || registers != SIRIUS_RETAINED_REGISTERS ||
        ni != programBase + 43 + 5 * count + SIRIUS_RETAINED_LAYERS + 1 ||
        rowWords != 204 + 5 * registers)
        return;
    for (uint i = 0; i < 112; ++i)
        if (ReadOriginal(i, row).valid == 0) return;
    RetainedTriple chart = ReadOriginal(111, row);
    if (!RTExact(chart, 1) && !RTExact(chart, -1)) return;
    if (RTAbsoluteUpper(ReadOriginal(1, row)) != 0 && RTAbsoluteUpper(ReadOriginal(3, row)) != 0)
        return;
    RetainedTriple h = ReadOriginal(104, row), fraction = ReadOriginal(105, row),
                   moving = ReadOriginal(110, row);
    if (h.hi <= 0 || RTAbsoluteLower(h) <= 0) return;
    if (!RTExact(fraction, 0) && !RTExact(fraction, 1)) {
        RetainedTriple remaining = RTSubtract(RTConstant(1), fraction);
        if (fraction.hi <= 0 || RTAbsoluteLower(fraction) <= 0 || remaining.hi <= 0 ||
            RTAbsoluteLower(remaining) <= 0)
            return;
    }
    if (!RTExact(moving, 0) && !RTExact(moving, 1)) return;
    if (RTExact(moving, 0))
        for (uint i = 106; i < 110; ++i)
            if (!RTExact(ReadOriginal(i, row), 0)) return;
    if (RTExact(ReadOriginal(0, row), 0) && RTExact(ReadOriginal(2, row), 0) &&
        RTExact(ReadOriginal(3, row), 0)) {
        if (lane != 0) return;
        if (FlatDense(row, rowBase, h, fraction, moving))
            if (lane == 0) outputs[rowBase] = 1;
        return;
    }
    if (!EvaluateProgram(row, rowBase, programBase, count, registers, 40, lane)) return;
    if (lane == 0) outputs[rowBase] = 1;
}
