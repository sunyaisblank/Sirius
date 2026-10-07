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
    r.hi = RPFromBits(ScratchRead(index, 0, rowBase, 204));
    r.lo = RPFromBits(ScratchRead(index, 1, rowBase, 204));
    r.tail = RPFromBits(ScratchRead(index, 2, rowBase, 204));
    r.error = RPFromBits(ScratchRead(index, 3, rowBase, 204));
    r.valid = ScratchRead(index, 4, rowBase, 204);
    return r;
}
void WriteValue(uint index, RetainedTriple r, uint rowBase) {
    ScratchWrite(index, 0, rowBase, 204, RPBits(r.hi));
    ScratchWrite(index, 1, rowBase, 204, RPBits(r.lo));
    ScratchWrite(index, 2, rowBase, 204, RPBits(r.tail));
    ScratchWrite(index, 3, rowBase, 204, RPBits(r.error));
    ScratchWrite(index, 4, rowBase, 204, r.valid);
}
RetainedTriple ReadOriginal(uint index, uint row) {
    uint base = 1 + 560 * row + 5 * index;
    RetainedTriple value = RTPack(RPFromBits(inputs[base]), RPFromBits(inputs[base + 1]),
                                  RPFromBits(inputs[base + 2]), RPFromBits(inputs[base + 3]));
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
        r = RTConstant(RPFromBits(a));
    else if (op == 1) {
        if (a >= 112) return false;
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
                     uint required, uint lane) {
    uint needed = count;
    if (needed == 0 || needed > count) return false;
    if (lane == 0) retainedStatus = 0;
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_DENSE)
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
#if defined(SIRIUS_RETAINED_PORTABLE) && !defined(SIRIUS_RETAINED_PARALLEL_DENSE)
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
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_DENSE)
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
        for (uint i = 0; i < 40; ++i)
            WriteStored(rowBase + 4 + i * 5, ReadValue(inputs[programBase + 3 + i], rowBase));
    }
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_DENSE)
    AllMemoryBarrierWithGroupSync();
#endif
    return true;
}

bool FlatDense(uint row, uint rowBase, RetainedTriple h, RetainedTriple s, RetainedTriple moving) {
#ifdef SIRIUS_RETAINED_PORTABLE
    // The same primitive dependencies as the native Hermite expression. The
    // last phase owns endpoint overrides even when an interior phase refuses.
    // One call site shares the existing portable arithmetic authority.
    RetainedTriple first = RTInvalid(), slope0 = RTInvalid(), slope1 = RTInvalid(),
                   delta = RTInvalid(), a = RTInvalid(), b = RTInvalid(),
                   quadratic = RTInvalid(), cubic = RTInvalid(), temporary = RTInvalid(),
                   hs = RTInvalid(), left = RTInvalid(), right = RTInvalid(),
                   value = RTInvalid(), slope = RTInvalid(), denominator = RTInvalid(),
                   numerator = RTInvalid(), shift = RTInvalid();
    const uint hermitePhases = 21, hermiteComponents = 20;
    const uint arrivalPhases = 76;
    [loop] for (uint phase = 0; phase < hermiteComponents * hermitePhases + arrivalPhases;
                ++phase) {
        bool hermite = phase < hermiteComponents * hermitePhases;
        uint step = 0, index = 0, column = 0, axis = 0;
        uint op = 0;
        RetainedTriple lhs = RTInvalid(), rhs = RTInvalid();
        if (hermite) {
            uint component = phase / hermitePhases;
            step = phase % hermitePhases;
            index = (component / 4) * 8 + component % 4;
            if (step == 0) {
                first = ReadOriginal(4 + index, row);
                slope0 = ReadOriginal(8 + index, row);
                slope1 = ReadOriginal(48 + index, row);
                delta = ReadOriginal(84 + component, row);
                lhs = delta; rhs = h; op = 5;
            } else if (step == 1) {
                lhs = temporary; rhs = slope0; op = 3;
            } else if (step == 2) {
                lhs = slope1; rhs = slope0; op = 3;
            } else if (step == 3) {
                lhs = RTConstant(RPFromUnsigned(3u)); rhs = a; op = 4;
            } else if (step == 4) {
                lhs = temporary; rhs = b; op = 3;
            } else if (step == 5) {
                lhs = RTConstant(RPFromUnsigned(2u)); rhs = a; op = 4;
            } else if (step == 6) {
                lhs = b; rhs = temporary; op = 3;
            } else if (step == 7) {
                lhs = cubic; rhs = s; op = 4;
            } else if (step == 8) {
                lhs = quadratic; rhs = temporary; op = 2;
            } else if (step == 9) {
                lhs = temporary; rhs = s; op = 4;
            } else if (step == 10) {
                lhs = slope0; rhs = temporary; op = 2;
            } else if (step == 11) {
                lhs = h; rhs = s; op = 4;
            } else if (step == 12) {
                lhs = temporary; rhs = hs; op = 4;
            } else if (step == 13) {
                lhs = first; rhs = value; op = 2;
            } else if (step == 14) {
                lhs = RTConstant(RPFromUnsigned(2u)); rhs = quadratic; op = 4;
            } else if (step == 15) {
                lhs = RTConstant(RPFromUnsigned(3u)); rhs = cubic; op = 4;
            } else if (step == 16) {
                lhs = right; rhs = s; op = 4;
            } else if (step == 17) {
                lhs = left; rhs = right; op = 2;
            } else if (step == 18) {
                lhs = temporary; rhs = s; op = 4;
            } else if (step == 19) {
                lhs = slope0; rhs = temporary; op = 2;
            } else {
                if (RTExact(s, RPZero())) {
                    value = first;
                    slope = slope0;
                }
                if (RTExact(s, RPOne())) {
                    lhs = first; rhs = delta; op = 2;
                    slope = slope1;
                }
            }
        } else {
            uint arrival = phase - hermiteComponents * hermitePhases;
            if (arrival == 0) {
                if (!RTExact(moving, RPOne())) return true;
                denominator = RTConstant(RPZero());
            }
            if (arrival < 8) {
                axis = arrival / 2;
                step = arrival % 2;
                if (step == 0) {
                    lhs = ReadOriginal(106 + axis, row);
                    rhs = ReadStored(rowBase + 4 + (4 + axis) * 5); op = 4;
                } else {
                    lhs = denominator; rhs = temporary; op = 2;
                }
            } else {
                column = (arrival - 8) / 17;
                step = (arrival - 8) % 17;
                if (step < 8) {
                    axis = step / 2;
                    if (step == 0) numerator = RTConstant(RPZero());
                    if ((step & 1u) == 0) {
                        lhs = ReadOriginal(106 + axis, row);
                        rhs = ReadStored(rowBase + 4 + (8 + 8 * column + axis) * 5); op = 4;
                    } else {
                        lhs = numerator; rhs = temporary; op = 2;
                    }
                } else if (step == 8) {
                    lhs = RTNegate(numerator); rhs = denominator; op = 5;
                } else {
                    axis = (step - 9) / 2;
                    index = 8 + 8 * column + axis;
                    if (((step - 9) & 1u) == 0) {
                        lhs = ReadStored(rowBase + 4 + (4 + axis) * 5);
                        rhs = shift; op = 4;
                    } else {
                        lhs = ReadStored(rowBase + 4 + index * 5);
                        rhs = temporary; op = 2;
                    }
                }
            }
        }
        if (op != 0) {
            RetainedTriple result = RTEvaluateArithmetic(lhs, rhs, op);
            if (hermite) {
                if (step == 1) a = result;
                else if (step == 2) b = result;
                else if (step == 4) quadratic = result;
                else if (step == 6) cubic = result;
                else if (step == 11) hs = result;
                else if (step == 12 || step == 13 || step == 20) value = result;
                else if (step == 14) left = result;
                else if (step == 15 || step == 16) right = result;
                else if (step == 19) slope = result;
                else temporary = result;
            } else {
                uint arrival = phase - hermiteComponents * hermitePhases;
                if (arrival < 8) {
                    if (step == 0) temporary = result;
                    else denominator = result;
                } else if (step < 8) {
                    if ((step & 1u) == 0) temporary = result;
                    else numerator = result;
                } else if (step == 8) {
                    shift = result;
                } else if (((step - 9) & 1u) == 0) {
                    temporary = result;
                } else {
                    value = result;
                    if (value.valid == 0) return false;
                    WriteStored(rowBase + 4 + index * 5, value);
                }
            }
        }
        if (hermite && step == 20) {
            if (value.valid == 0 || slope.valid == 0) return false;
            WriteStored(rowBase + 4 + index * 5, value);
            WriteStored(rowBase + 4 + (index + 4) * 5, slope);
        }
    }
    return true;
#else
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
            RetainedTriple quadratic = RTSubtract(RTMultiply(RTConstant(RPFromUnsigned(3u)), a), b);
            RetainedTriple cubic = RTSubtract(b, RTMultiply(RTConstant(RPFromUnsigned(2u)), a));
            RetainedTriple value = RTAdd(
                first,
                RTMultiply(RTAdd(slope0, RTMultiply(RTAdd(quadratic, RTMultiply(cubic, s)), s)),
                           RTMultiply(h, s)));
            RetainedTriple slope = RTAdd(
                slope0,
                RTMultiply(RTAdd(RTMultiply(RTConstant(RPFromUnsigned(2u)), quadratic),
                                 RTMultiply(RTMultiply(RTConstant(RPFromUnsigned(3u)), cubic), s)),
                           s));
            if (RTExact(s, RPZero())) {
                value = first;
                slope = slope0;
            }
            if (RTExact(s, RPOne())) {
                value = RTAdd(first, delta);
                slope = slope1;
            }
            if (value.valid == 0 || slope.valid == 0) return false;
            WriteStored(rowBase + 4 + index * 5, value);
            WriteStored(rowBase + 4 + (index + 4) * 5, slope);
        }
    if (RTExact(moving, RPOne())) {
        RetainedTriple denominator = RTConstant(RPZero());
        for (uint axis = 0; axis < 4; ++axis)
            denominator = RTAdd(denominator, RTMultiply(ReadOriginal(106 + axis, row),
                                                        ReadStored(rowBase + 4 + (4 + axis) * 5)));
        for (uint column = 0; column < 4; ++column) {
            RetainedTriple numerator = RTConstant(RPZero());
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
#endif
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
    if (rows == 0 || rows > 65536 || row >= rows || ni < 1 + 560 * rows + 43 || no % rows != 0)
        return;
    uint programBase = 1 + 560 * rows, count = inputs[programBase],
         registers = inputs[programBase + 1];
    uint rowWords = no / rows, rowBase = row * rowWords;
    if (rowWords < 204) return;
    for (uint i = lane; i < rowWords; i += SIRIUS_RETAINED_EXECUTION_LANES) outputs[rowBase + i] = 0;
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_DENSE)
    AllMemoryBarrierWithGroupSync();
#endif
    if (inputs[programBase + 2] != 40 || count > 65536 || registers != SIRIUS_RETAINED_REGISTERS ||
        ni != programBase + 43 + 5 * count + SIRIUS_RETAINED_LAYERS + 1 ||
        rowWords != 204 + 5 * registers)
        return;
    for (uint i = 0; i < 112; ++i)
        if (ReadOriginal(i, row).valid == 0) return;
    RetainedTriple chart = ReadOriginal(111, row);
    if (!RTExact(chart, RPOne()) && !RTExact(chart, RPScalarNegate(RPOne()))) return;
    if (!RPScalarEqual(RTAbsoluteUpper(ReadOriginal(1, row)), RPZero()) &&
        !RPScalarEqual(RTAbsoluteUpper(ReadOriginal(3, row)), RPZero()))
        return;
    RetainedTriple h = ReadOriginal(104, row), fraction = ReadOriginal(105, row),
                   moving = ReadOriginal(110, row);
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
    if (!RTExact(moving, RPZero()) && !RTExact(moving, RPOne())) return;
    if (RTExact(moving, RPZero()))
        for (uint i = 106; i < 110; ++i)
            if (!RTExact(ReadOriginal(i, row), RPZero())) return;
    if (RTExact(ReadOriginal(0, row), RPZero()) && RTExact(ReadOriginal(2, row), RPZero()) &&
        RTExact(ReadOriginal(3, row), RPZero())) {
        if (lane != 0) return;
        if (FlatDense(row, rowBase, h, fraction, moving))
            if (lane == 0) outputs[rowBase] = 1;
        return;
    }
    if (!EvaluateProgram(row, rowBase, programBase, count, registers, 40, lane)) return;
    if (lane == 0) outputs[rowBase] = 1;
}
