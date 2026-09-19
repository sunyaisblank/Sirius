import retained_scratch;
// Batched retained camera launch. All rows share an immutable arithmetic program.
import retained_pair;
[[vk::binding(0, 0)]] StructuredBuffer<uint> inputs;
[[vk::binding(1, 0)]] RWStructuredBuffer<uint> outputs;
uint ReadInput(uint index, uint row, uint rows) {
    return inputs[index < 32 ? 1 + 160 * row + 5 * index : 1 + 160 * rows + index - 32];
}
RetainedPair ReadOriginal(uint index, uint row) {
    uint base = 1 + 160 * row + 5 * index;
    if (inputs[base + 4] != 1) return RPInvalid();
    RetainedPair r =
        RPPack(asfloat(inputs[base]), asfloat(inputs[base + 1]), asfloat(inputs[base + 3]));
    // Retain the third input term through the camera pair's rounding radius.
    return RPAdd(r, RPConstant(asfloat(inputs[base + 2])));
}
bool Exact(RetainedPair value, float x) {
    return value.valid != 0 && value.hi == x && value.lo == 0 && value.error == 0;
}
RetainedPair ReadValue(uint index, uint rowBase) {
    RetainedPair r;
    r.hi = asfloat(ScratchRead(index, 0, rowBase, 512));
    r.lo = asfloat(ScratchRead(index, 1, rowBase, 512));
    r.error = asfloat(ScratchRead(index, 2, rowBase, 512));
    r.valid = ScratchRead(index, 3, rowBase, 512);
    return r;
}
void WriteValue(uint index, RetainedPair r, uint rowBase) {
    ScratchWrite(index, 0, rowBase, 512, asuint(r.hi));
    ScratchWrite(index, 1, rowBase, 512, asuint(r.lo));
    ScratchWrite(index, 2, rowBase, 512, asuint(r.error));
    ScratchWrite(index, 3, rowBase, 512, r.valid);
}
bool EvaluateOne(uint index, uint row, uint rowBase, uint rows, uint registers) {
    uint address = 138 + 5 * index, op = ReadInput(address, row, rows),
         destination = ReadInput(address + 1, row, rows), a = ReadInput(address + 2, row, rows),
         b = ReadInput(address + 3, row, rows), c = ReadInput(address + 4, row, rows);
    if (destination >= registers || op > 11) return false;
    RetainedPair r = RPInvalid();
    if (op == 0)
        r = RPConstant(asfloat(a));
    else if (op == 1) {
        if (a >= 32) return false;
        r = ReadOriginal(a, row);
    } else {
        if (a >= registers) return false;
        RetainedPair av = ReadValue(a, rowBase);
        if (op == 6)
            r = RPSqrt(av);
        else if (op == 7)
            r = RPNegate(av);
        else if (op == 8) {
            float bound = RPAbsoluteUpper(av);
            if (av.valid != 0 && bound > 0)
                r = RPConstant(RPPowerOfTwo(RPFloorExponent(bound) + 1));
        } else {
            if (b >= registers) return false;
            RetainedPair bv = ReadValue(b, rowBase);
            if (op == 2)
                r = RPAdd(av, bv);
            else if (op == 3)
                r = RPSubtract(av, bv);
            else if (op == 4)
                r = RPMultiply(av, bv);
            else if (op == 5)
                r = RPDivide(av, bv);
            else if (op == 9) {
                if (av.valid != 0 && bv.valid != 0)
                    r = RPConstant(max(RPAbsoluteUpper(av), RPAbsoluteUpper(bv)));
            } else if (op == 10 || op == 11) {
                if (c >= registers) return false;
                RetainedPair cv = ReadValue(c, rowBase);
                if (cv.valid != 0) {
                    bool choose =
                        op == 10 ? cv.hi >= 0 : (cv.hi == 0 && cv.lo == 0 && cv.error == 0);
                    r = choose ? av : bv;
                }
            }
        }
    }
    WriteValue(destination, r, rowBase);

    return true;
}
[shader("compute")][numthreads(SIRIUS_RETAINED_LANES, 1, 1)] void ComputeMain(uint3 id
                                                                              : SV_GroupID,
                                                                                uint lane
                                                                              : SV_GroupIndex) {
    if (id.y != 0 || id.z != 0) return;
    uint rawNi, si, rawNo, so;
    inputs.GetDimensions(rawNi, si);
    outputs.GetDimensions(rawNo, so);
    if (rawNi < 1) return;
    uint rows = inputs[0];
    if (rows == 0 || rows > 65536 || id.x >= rows || rawNi < 1 + 160 * rows || rawNo % rows != 0)
        return;
    uint row = id.x, no = rawNo / rows, ni = rawNi - 1 - 160 * rows + 32, rowBase = row * no;
    if (no < 512) return;
    for (uint i = lane; i < no; i += SIRIUS_RETAINED_LANES) outputs[rowBase + i] = 0;
    AllMemoryBarrierWithGroupSync();
    if (ni < 138) return;
    uint count = ReadInput(32, row, rows), registers = ReadInput(33, row, rows);
    if (count > 65536 || registers != SIRIUS_RETAINED_REGISTERS ||
        ni != 138 + 5 * count + SIRIUS_RETAINED_LAYERS + 1 || no != 512 + 4 * registers)
        return;
    if (lane == 0) {
        outputs[rowBase + 0] = asuint(271832.0f);
        outputs[rowBase + 1] = asuint(32.0f);
#ifdef SIRIUS_FP32_COMP
        outputs[rowBase + 2] = asuint(1.0f);
#endif
        outputs[rowBase + 4] = ReadInput(30, row, rows);
        outputs[rowBase + 5] = asuint(float(ni));
        outputs[rowBase + 6] = asuint(float(no));
        for (uint i = 0; i < 160; ++i) outputs[rowBase + 8 + i] = inputs[1 + 160 * row + i];
    }
    for (uint i = 0; i < 32; ++i)
        if (ReadOriginal(i, row).valid == 0) return;
    if (!Exact(ReadOriginal(31, row), 0)) return;
    for (uint i = 20; i <= 24; ++i) {
        RetainedPair value = ReadOriginal(i, row);
        if (value.hi <= 0 || RPAbsoluteLower(value) <= 0) return;
    }
    RetainedPair thin = ReadOriginal(29, row);
    if (!Exact(thin, 1) &&
        (!Exact(thin, 0) || !Exact(ReadOriginal(27, row), 0) || !Exact(ReadOriginal(28, row), 0)))
        return;
    if (RPAbsoluteUpper(ReadOriginal(1, row)) != 0 && RPAbsoluteUpper(ReadOriginal(3, row)) != 0)
        return;
    if (lane == 0) retainedStatus = 0;
    AllMemoryBarrierWithGroupSync();
    uint layerBase = 138 + 5 * count;
    if (ReadInput(layerBase, row, rows) != 0 ||
        ReadInput(layerBase + SIRIUS_RETAINED_LAYERS, row, rows) != count)
        return;
    uint expected = 0;
    for (uint layer = 0; layer < SIRIUS_RETAINED_LAYERS; ++layer) {
        uint start = ReadInput(layerBase + layer, row, rows),
             end = ReadInput(layerBase + layer + 1, row, rows);
        if (start != expected || end <= start || end > count || end - start > SIRIUS_RETAINED_LANES)
            RejectRetainedProgram();
        else if (start + lane < end && !EvaluateOne(start + lane, row, rowBase, rows, registers))
            RejectRetainedProgram();
        GroupMemoryBarrierWithGroupSync();
        if (retainedStatus != 0) return;
        expected = end;
    }

    for (uint i = 0; i < 104; ++i) {
        uint index = ReadInput(34 + i, row, rows);
        if (index >= registers || ReadValue(index, rowBase).valid == 0) return;
    }
    if (lane != 0) return;
    for (uint i = 0; i < 104; ++i) {
        uint index = ReadInput(34 + i, row, rows);
        RetainedPair r = ReadValue(index, rowBase);
        outputs[rowBase + 168 + 3 * i] = asuint(r.hi);
        outputs[rowBase + 169 + 3 * i] = asuint(r.lo);
        outputs[rowBase + 170 + 3 * i] = asuint(r.error);
    }
    outputs[rowBase + 3] = asuint(1.0f);
}
