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
RetainedTriple ReadOriginal(uint index, uint row) {
    uint base = 1 + row * 230 + 5 * index;
    RetainedTriple r = RTPack(RPFromBits(inputs[base]), RPFromBits(inputs[base + 1]),
                              RPFromBits(inputs[base + 2]), RPFromBits(inputs[base + 3]));
    if (inputs[base + 4] != 1) r.valid = 0;
    return r;
}
RetainedTriple ReadInput(uint index, uint row, uint rowBase) {
    return index >= 4 && index < 44 ? ReadStored(rowBase + 804 + 5 * (index - 4))
                                    : ReadOriginal(index, row);
}
RetainedTriple ReadValue(uint index, uint rowBase) {
    RetainedTriple r;
    r.hi = RPFromBits(ScratchRead(index, 0, rowBase, 2404));
    r.lo = RPFromBits(ScratchRead(index, 1, rowBase, 2404));
    r.tail = RPFromBits(ScratchRead(index, 2, rowBase, 2404));
    r.error = RPFromBits(ScratchRead(index, 3, rowBase, 2404));
    r.valid = ScratchRead(index, 4, rowBase, 2404);
    return r;
}
void WriteValue(uint index, RetainedTriple r, uint rowBase) {
    ScratchWrite(index, 0, rowBase, 2404, RPBits(r.hi));
    ScratchWrite(index, 1, rowBase, 2404, RPBits(r.lo));
    ScratchWrite(index, 2, rowBase, 2404, RPBits(r.tail));
    ScratchWrite(index, 3, rowBase, 2404, RPBits(r.error));
    ScratchWrite(index, 4, rowBase, 2404, r.valid);
}

bool EvaluateOne(uint index, uint row, uint rowBase, uint programBase, uint registers) {
    uint address = programBase + 43 + 5 * index, op = inputs[address],
         destination = inputs[address + 1], a = inputs[address + 2], b = inputs[address + 3],
         c = inputs[address + 4];
    if (destination >= registers || op > 11) return false;
    RetainedTriple r = RTInvalid();
    if (op == 0)
        r = RTConstant(RPFromBits(a));
    else if (op == 1) {
        if (a >= 45) return false;
        r = ReadInput(a, row, rowBase);
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
#ifndef SIRIUS_RETAINED_PORTABLE
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
                     uint stage, uint lane) {
    uint needed = count;
    if (needed == 0 || needed > count) return false;
    if (lane == 0) retainedStatus = 0;
#ifndef SIRIUS_RETAINED_PORTABLE
    GroupMemoryBarrierWithGroupSync();
#endif
    uint layerBase = programBase + 43 + 5 * count;
    if (inputs[layerBase] != 0 || inputs[layerBase + SIRIUS_RETAINED_LAYERS] != count) return false;
    uint expected = 0;
    for (uint layer = 0; layer < SIRIUS_RETAINED_LAYERS; ++layer) {
        uint start = inputs[layerBase + layer], end = inputs[layerBase + layer + 1];
        if (start >= needed) break;
        if (start != expected || end <= start || end > needed ||
            end - start > SIRIUS_RETAINED_LANES)
            RejectRetainedProgram();
#ifdef SIRIUS_RETAINED_PORTABLE
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
#ifndef SIRIUS_RETAINED_PORTABLE
        GroupMemoryBarrierWithGroupSync();
#endif
        if (retainedStatus != 0) return false;
        expected = end;
    }

    for (uint i = 0; i < 40; ++i) {
        uint index = inputs[programBase + 3 + i];
        if (index >= registers || ReadValue(index, rowBase).valid == 0) return false;
    }
    if (lane == 0) {
        for (uint i = 0; i < 40; ++i)
            WriteStored(rowBase + 1004 + stage * 200 + i * 5,
                        ReadValue(inputs[programBase + 3 + i], rowBase));
    }
#ifndef SIRIUS_RETAINED_PORTABLE
    AllMemoryBarrierWithGroupSync();
#endif
    return true;
}
RetainedTriple Rational(int numerator, int denominator) {
    return RTDivide(RTConstant(RPFromSigned(numerator)), RTConstant(RPFromSigned(denominator)));
}
static const int cn[7] = {0, 1, 3, 4, 8, 1, 1};
static const int cd[7] = {1, 5, 10, 5, 9, 1, 1};
static const int an[7][7] = {{0, 0, 0, 0, 0, 0, 0},
                             {1, 0, 0, 0, 0, 0, 0},
                             {3, 9, 0, 0, 0, 0, 0},
                             {44, -56, 32, 0, 0, 0, 0},
                             {19372, -25360, 64448, -212, 0, 0, 0},
                             {9017, -355, 46732, 49, -5103, 0, 0},
                             {35, 0, 500, 125, -2187, 11, 0}};
static const int ad[7][7] = {{1, 1, 1, 1, 1, 1, 1},
                             {5, 1, 1, 1, 1, 1, 1},
                             {40, 40, 1, 1, 1, 1, 1},
                             {45, 15, 9, 1, 1, 1, 1},
                             {6561, 2187, 6561, 729, 1, 1, 1},
                             {3168, 33, 5247, 176, 18656, 1, 1},
                             {384, 1, 1113, 192, 6784, 84, 1}};
static const int en[7] = {71, 0, -71, 71, -17253, 22, -1};
static const int ed[7] = {57600, 1, 16695, 1920, 339200, 525, 40};
#ifndef SIRIUS_RETAINED_PORTABLE
// Six row sums, fourteen nonzero corrections and five embedded-error weights.
// Preserve the complete retained value, including its integer validity word.
groupshared uint retainedCoefficients[SIRIUS_RETAINED_COEFFICIENTS * 5];
void CacheCoefficient(uint lane) {
    if (lane >= SIRIUS_RETAINED_COEFFICIENTS) return;
    RetainedTriple value = RTInvalid();
    if (lane < 6)
        value = Rational(cn[lane + 1], cd[lane + 1]);
    else if (lane >= 20)
        value = Rational(en[lane - 18], ed[lane - 18]);
    else {
        uint index = 6;
        for (uint stage = 2; stage < 7; ++stage)
            for (uint previous = 1; previous < stage; ++previous)
                if (an[stage][previous] != 0) {
                    if (lane == index) value = Rational(an[stage][previous], ad[stage][previous]);
                    ++index;
                }
    }
    uint base = lane * 5;
    retainedCoefficients[base] = RPBits(value.hi);
    retainedCoefficients[base + 1] = RPBits(value.lo);
    retainedCoefficients[base + 2] = RPBits(value.tail);
    retainedCoefficients[base + 3] = RPBits(value.error);
    retainedCoefficients[base + 4] = value.valid;
}
RetainedTriple ReadCoefficient(uint index) {
    uint base = index * 5;
    RetainedTriple value;
    value.hi = RPFromBits(retainedCoefficients[base]);
    value.lo = RPFromBits(retainedCoefficients[base + 1]);
    value.tail = RPFromBits(retainedCoefficients[base + 2]);
    value.error = RPFromBits(retainedCoefficients[base + 3]);
    value.valid = retainedCoefficients[base + 4];
    return value;
}
#endif
[shader("compute")][numthreads(SIRIUS_RETAINED_EXECUTION_LANES, 1, 1)] void ComputeMain(uint3 id
                                                                              : SV_GroupID,
                                                                                uint lane
                                                                              : SV_GroupIndex) {
    uint ni, si, no, so;
    inputs.GetDimensions(ni, si);
    outputs.GetDimensions(no, so);
    if (ni < 1 || id.y != 0 || id.z != 0) return;
    uint rows = inputs[0], row = id.x;
    if (rows == 0 || rows > 65536 || row >= rows || ni < 1 + 230 * rows + 43 || no % rows != 0)
        return;
    uint programBase = 1 + 230 * rows, count = inputs[programBase],
         registers = inputs[programBase + 1];
    uint rowWords = no / rows, rowBase = row * rowWords;
    if (rowWords < 2404) return;
    for (uint i = lane; i < rowWords; i += SIRIUS_RETAINED_EXECUTION_LANES) outputs[rowBase + i] = 0;
#ifndef SIRIUS_RETAINED_PORTABLE
    AllMemoryBarrierWithGroupSync();
#endif
    if (inputs[programBase + 2] != 40 || count > 65536 || registers != SIRIUS_RETAINED_REGISTERS ||
        ni != programBase + 43 + 5 * count + SIRIUS_RETAINED_LAYERS + 1 ||
        rowWords != 2404 + 5 * registers)
        return;
    for (uint i = 0; i < 46; ++i)
        if (ReadOriginal(i, row).valid == 0) return;
    RetainedTriple chart = ReadOriginal(44, row);
    if ((!RPScalarEqual(chart.hi, RPOne()) && !RPScalarEqual(chart.hi, RPScalarNegate(RPOne()))) ||
        !RPScalarEqual(chart.lo, RPZero()) || !RPScalarEqual(chart.tail, RPZero()) ||
        !RPScalarEqual(chart.error, RPZero()))
        return;
    RetainedTriple spin = ReadOriginal(1, row), lambda = ReadOriginal(3, row);
    if (!RPScalarEqual(RTAbsoluteUpper(spin), RPZero()) &&
        !RPScalarEqual(RTAbsoluteUpper(lambda), RPZero()))
        return;
    RetainedTriple h = ReadOriginal(45, row);
    if (RPScalarLessEqual(h.hi, RPZero()) || RPScalarLessEqual(RTAbsoluteLower(h), RPZero()))
        return;
    if (RTExact(ReadOriginal(0, row), RPZero()) && RTExact(ReadOriginal(2, row), RPZero()) &&
        RTExact(lambda, RPZero())) {
        if (lane != 0) return;
        // All seven Hamiltonian RHS stages are the same constant eta*p and
        // eta*P. The defining tableau row sums are exactly one and zero.
        for (uint i = 0; i < 40; ++i) {
            RetainedTriple delta = RTConstant(RPZero());
            if (i % 8 < 4) {
                RetainedTriple slope = ReadOriginal(8 + i, row);
                if (i % 8 == 0) slope = RTNegate(slope);
                delta = RTMultiply(h, slope);
            }
            RetainedTriple full = RTAdd(ReadOriginal(4 + i, row), delta);
            if (full.valid == 0 || delta.valid == 0) return;
            if (lane == 0) WriteStored(rowBase + 4 + i * 5, full);
            if (lane == 0) WriteStored(rowBase + 204 + i * 5, full);
            if (lane == 0) WriteStored(rowBase + 404 + i * 5, delta);
            if (lane == 0) WriteStored(rowBase + 604 + i * 5, RTConstant(RPZero()));
        }
        if (lane == 0) outputs[rowBase] = 1;
        if (lane == 0) outputs[rowBase + 1] = 7;
        return;
    }
    for (uint stage = 0; stage < 7; ++stage) {
        if (lane == 0) retainedStatus = 0;
#ifndef SIRIUS_RETAINED_PORTABLE
        // Stage zero has succeeded. Publish once at the existing stage barrier;
        // expression register reuse never writes this separate immutable cache.
        if (stage == 1) CacheCoefficient(lane);
        GroupMemoryBarrierWithGroupSync();
#endif
        for (uint i = lane; i < 40; i += SIRIUS_RETAINED_EXECUTION_LANES) {
            RetainedTriple value = ReadOriginal(4 + i, row);
            if (stage != 0) {
                RetainedTriple first = ReadStored(rowBase + 1004 + i * 5);
#ifdef SIRIUS_RETAINED_PORTABLE
                RetainedTriple weighted = RTMultiply(first, Rational(cn[stage], cd[stage]));
#else
                RetainedTriple weighted = RTMultiply(first, ReadCoefficient(stage - 1));
#endif
                for (uint j = 1; j < stage; ++j)
                    if (an[stage][j] != 0) {
                        RetainedTriple difference =
                            RTSubtract(ReadStored(rowBase + 1004 + j * 200 + i * 5), first);
#ifdef SIRIUS_RETAINED_PORTABLE
                        weighted = RTAdd(
                            weighted, RTMultiply(difference, Rational(an[stage][j], ad[stage][j])));
#else
                        uint index = 6 + (stage - 2) * (stage - 1) / 2 + j - 1;
                        if (stage == 6) --index; // The preceding correction is zero.
                        weighted = RTAdd(weighted, RTMultiply(difference, ReadCoefficient(index)));
#endif
                    }
                RetainedTriple increment = RTMultiply(weighted, h);
                value = RTAdd(value, increment);
                if (stage == 6) WriteStored(rowBase + 4 + 80 * 5 + i * 5, increment);
            }
            if (value.valid == 0) RejectRetainedProgram();
            WriteStored(rowBase + 804 + i * 5, value);
        }
#ifndef SIRIUS_RETAINED_PORTABLE
        AllMemoryBarrierWithGroupSync();
#endif
        if (retainedStatus != 0) return;
        if (lane == 0) outputs[rowBase + 1] = stage + 1;
        if (!EvaluateProgram(row, rowBase, programBase, count, registers, stage, lane)) return;
    }
    if (lane == 0) retainedStatus = 0;
#ifndef SIRIUS_RETAINED_PORTABLE
    GroupMemoryBarrierWithGroupSync();
#endif
    for (uint i = lane; i < 40; i += SIRIUS_RETAINED_EXECUTION_LANES) {
        RetainedTriple first = ReadStored(rowBase + 1004 + i * 5), weighted = RTConstant(RPZero());
        for (uint j = 1; j < 7; ++j)
            if (en[j] != 0) {
                RetainedTriple difference =
                    RTSubtract(ReadStored(rowBase + 1004 + j * 200 + i * 5), first);
#ifdef SIRIUS_RETAINED_PORTABLE
                weighted = RTAdd(weighted, RTMultiply(difference, Rational(en[j], ed[j])));
#else
                weighted = RTAdd(weighted, RTMultiply(difference, ReadCoefficient(18 + j)));
#endif
            }
        RetainedTriple error = RTMultiply(weighted, h), full = ReadStored(rowBase + 804 + i * 5);
        RetainedTriple lower = RTSubtract(full, error);
        if (error.valid == 0 || lower.valid == 0) RejectRetainedProgram();
        WriteStored(rowBase + 4 + i * 5, full);
        WriteStored(rowBase + 4 + 40 * 5 + i * 5, lower);
        WriteStored(rowBase + 4 + 120 * 5 + i * 5, error);
    }
#ifndef SIRIUS_RETAINED_PORTABLE
    AllMemoryBarrierWithGroupSync();
#endif
    if (retainedStatus != 0) return;
    if (lane == 0) outputs[rowBase] = 1;
}
