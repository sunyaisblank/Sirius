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
    r.hi = RPFromBits(ScratchRead(index, 0, rowBase, 769));
    r.lo = RPFromBits(ScratchRead(index, 1, rowBase, 769));
    r.tail = RPFromBits(ScratchRead(index, 2, rowBase, 769));
    r.error = RPFromBits(ScratchRead(index, 3, rowBase, 769));
    r.valid = ScratchRead(index, 4, rowBase, 769);
    return r;
}
void WriteValue(uint index, RetainedTriple r, uint rowBase) {
    ScratchWrite(index, 0, rowBase, 769, RPBits(r.hi));
    ScratchWrite(index, 1, rowBase, 769, RPBits(r.lo));
    ScratchWrite(index, 2, rowBase, 769, RPBits(r.tail));
    ScratchWrite(index, 3, rowBase, 769, RPBits(r.error));
    ScratchWrite(index, 4, rowBase, 769, r.valid);
}
bool EvaluateOne(uint index, uint row, uint rowBase, uint programBase, uint registers) {
    uint address = programBase + 103 + 5 * index, op = inputs[address],
         destination = inputs[address + 1], a = inputs[address + 2], b = inputs[address + 3],
         c = inputs[address + 4];
    if (destination >= registers || op > 11) return false;
    RetainedTriple r = RTInvalid();
    if (op == 0)
        r = RTConstant(RPFromBits(a));
    else if (op == 1) {
        if (a >= 53) return false;
        r = ReadStored(rowBase + 504 + 5 * a);
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
    uint needed = required == 20 ? SIRIUS_RETAINED_PREFIX : count;
    // The generated prefix owns every metric/tangent dependency until the
    // continuation consumes it. Projected inputs are read only in the suffix.
    uint first = required == 20 ? 0 : SIRIUS_RETAINED_PREFIX;
    if (needed == 0 || needed > count || first > needed) return false;
    if (lane == 0) retainedStatus = 0;
#ifndef SIRIUS_RETAINED_PORTABLE
    GroupMemoryBarrierWithGroupSync();
#endif
    uint layerBase = programBase + 103 + 5 * count;
    if (inputs[layerBase] != 0 || inputs[layerBase + SIRIUS_RETAINED_LAYERS] != count) return false;
    uint expected = 0;
    for (uint layer = 0; layer < SIRIUS_RETAINED_LAYERS; ++layer) {
        uint start = inputs[layerBase + layer], end = inputs[layerBase + layer + 1];
        if (start >= needed) break;
        if (start != expected || end <= start || end > needed ||
            end - start > SIRIUS_RETAINED_LANES || (start < first && end > first))
            RejectRetainedProgram();
        else if (start < first) {
            // Only validated prefix layers skip their now-unneeded barrier.
            // A malformed layer still synchronizes and consumes rejection.
            expected = end;
            continue;
        }
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

    if (expected != needed) return false;
    for (uint i = 0; i < required; ++i) {
        uint index = inputs[programBase + 3 + i];
        if (index >= registers || ReadValue(index, rowBase).valid == 0) return false;
    }
    if (lane == 0) {
        for (uint i = 0; i < required; ++i)
            WriteStored(rowBase + 4 + i * 5, ReadValue(inputs[programBase + 3 + i], rowBase));
    }
#ifndef SIRIUS_RETAINED_PORTABLE
    AllMemoryBarrierWithGroupSync();
#endif
    return true;
}

RetainedTriple Metric(uint rowBase, uint i, uint j) {
    return ReadStored(rowBase + 4 + 5 * (4 * i + j));
}
RetainedTriple Tangent(uint rowBase, uint i) { return ReadStored(rowBase + 4 + 5 * (16 + i)); }
#ifdef SIRIUS_RETAINED_PORTABLE
bool Project(uint rowBase, out uint component, out RetainedTriple root) {
    component = 4;
    root = RTInvalid();
    RetainedTriple best = RTInvalid();
    RetainedTriple defect = RTConstant(RPZero());
    RPScalar scale = RPZero();
    for (uint i = 0; i < 4; ++i)
        for (uint j = 0; j < 4; ++j) {
            RetainedTriple term = RTInvalid();
            // Preserve ti*tj followed by metric*term, including invalid terms.
            [loop] for (uint product = 0; product < 2; ++product)
                term = RTMultiply(product == 0 ? Tangent(rowBase, i) : Metric(rowBase, i, j),
                                  product == 0 ? Tangent(rowBase, j) : term);
            defect = RTAdd(defect, term);
            scale = RPUpAdd(scale, RTAbsoluteUpper(term));
        }
    if (defect.valid == 0 || RPScalarLessEqual(scale, RPZero()) ||
        RPScalarGreater(RPRatioUpper(RTAbsoluteUpper(defect), scale), RPFromBits(0x358637bdu)))
        return false;
    // Prefer a temporal root; if none is represented, preserve time orientation
    // and choose the smallest normalized spatial correction.
    for (uint c = 0; c < 4; ++c) {
        if (c > 0 && component == 0) break;
        RetainedTriple a = Metric(rowBase, c, c), b = RTConstant(RPZero()),
                       d = RTConstant(RPZero());
        for (uint i = 0; i < 4; ++i)
            if (i != c) {
                b = RTAdd(b, RTMultiply(RTAdd(Metric(rowBase, c, i), Metric(rowBase, i, c)),
                                        Tangent(rowBase, i)));
                for (uint j = 0; j < 4; ++j)
                    if (j != c)
                    {
                        RetainedTriple term = RTInvalid();
                        [loop] for (uint product = 0; product < 2; ++product)
                            term = RTMultiply(product == 0 ? Tangent(rowBase, i) : Metric(rowBase, i, j),
                                              product == 0 ? Tangent(rowBase, j) : term);
                        d = RTAdd(d, term);
                    }
            }
        RetainedTriple first = RTInvalid(), second = RTInvalid();
        if (RTExact(a, RPZero()))
            first = RTDivide(RTNegate(d), b);
        else {
            RetainedTriple square_b = RTInvalid(), product_term = RTInvalid();
            [loop] for (uint product = 0; product < 3; ++product) {
                product_term = RTMultiply(product == 0 ? b : (product == 1 ? a : RTConstant(RPFromUnsigned(4u))),
                                          product == 0 ? b : (product == 1 ? d : product_term));
                if (product == 0) square_b = product_term;
            }
            RetainedTriple discriminant = RTSubtract(square_b, product_term);
            RetainedTriple square = RTSqrt(discriminant);
            RetainedTriple q = RTMultiply(
                RTConstant(RPFromBits(0xbf000000u)),
                RTAdd(b, RPScalarGreaterEqual(b.hi, RPZero()) ? square : RTNegate(square)));
            if (RTExact(q, RPZero()))
                first = RTDivide(RTNegate(b), RTMultiply(RTConstant(RPFromUnsigned(2u)), a));
            else {
                // Both divisions execute in order even if the first is invalid.
                [loop] for (uint division = 0; division < 2; ++division) {
                    RetainedTriple quotient = RTDivide(division == 0 ? q : d, division == 0 ? a : q);
                    if (division == 0) first = quotient;
                    else second = quotient;
                }
            }
        }
        for (uint r = 0; r < 2; ++r) {
            RetainedTriple candidate = r == 0 ? first : second;
            if (candidate.valid == 0) continue;
            RetainedTriple denominator = RTConstant(RPZero());
            RPScalar denominatorScale = RPZero();
            for (uint i = 0; i < 4; ++i) {
                RetainedTriple term =
                    RTMultiply(Metric(rowBase, c, i), i == c ? candidate : Tangent(rowBase, i));
                denominator = RTAdd(denominator, term);
                denominatorScale = RPUpAdd(denominatorScale, RTAbsoluteUpper(term));
            }
            // A central root with a singular differentiated constraint cannot
            // own the coupled candidate. Try the represented spatial roots.
            if (denominator.valid == 0 ||
                RPScalarLessEqual(RTAbsoluteLower(denominator),
                                  RPUpMultiply(RPFromBits(0x29800000u), denominatorScale)))
                continue;
            RetainedTriple difference = RTSubtract(candidate, Tangent(rowBase, c));
            RetainedTriple original = RTCenter(Tangent(rowBase, c));
            RetainedTriple absolute =
                RTCenter(RPScalarLess(difference.hi, RPZero()) ? RTNegate(difference) : difference);
            RetainedTriple scaleValue =
                RTAdd(RTConstant(RPOne()),
                      RPScalarLess(original.hi, RPZero()) ? RTNegate(original) : original);
            RetainedTriple correction = RTDivide(absolute, scaleValue);
            if (correction.valid != 0 &&
                (best.valid == 0 || RPScalarLess(RTSubtract(correction, best).hi, RPZero()))) {
                best = correction;
                component = c;
                root = candidate;
            }
        }
    }
    return component < 4 && root.valid != 0;
}
#else
bool Project(uint rowBase, out uint component, out RetainedTriple root) {
    component = 4;
    root = RTInvalid();
    RetainedTriple best = RTInvalid();
    RetainedTriple defect = RTConstant(RPZero());
    RPScalar scale = RPZero();
    for (uint i = 0; i < 4; ++i)
        for (uint j = 0; j < 4; ++j) {
            RetainedTriple term = RTMultiply(Metric(rowBase, i, j),
                                             RTMultiply(Tangent(rowBase, i), Tangent(rowBase, j)));
            defect = RTAdd(defect, term);
            scale = RPUpAdd(scale, RTAbsoluteUpper(term));
        }
    if (defect.valid == 0 || RPScalarLessEqual(scale, RPZero()) ||
        RPScalarGreater(RPRatioUpper(RTAbsoluteUpper(defect), scale), RPFromBits(0x358637bdu)))
        return false;
    // Prefer a temporal root; if none is represented, preserve time orientation
    // and choose the smallest normalized spatial correction.
    for (uint c = 0; c < 4; ++c) {
        if (c > 0 && component == 0) break;
        RetainedTriple a = Metric(rowBase, c, c), b = RTConstant(RPZero()),
                       d = RTConstant(RPZero());
        for (uint i = 0; i < 4; ++i)
            if (i != c) {
                b = RTAdd(b, RTMultiply(RTAdd(Metric(rowBase, c, i), Metric(rowBase, i, c)),
                                        Tangent(rowBase, i)));
                for (uint j = 0; j < 4; ++j)
                    if (j != c)
                        d = RTAdd(d,
                                  RTMultiply(Metric(rowBase, i, j),
                                             RTMultiply(Tangent(rowBase, i), Tangent(rowBase, j))));
            }
        RetainedTriple first = RTInvalid(), second = RTInvalid();
        if (RTExact(a, RPZero()))
            first = RTDivide(RTNegate(d), b);
        else {
            RetainedTriple discriminant = RTSubtract(
                RTMultiply(b, b), RTMultiply(RTConstant(RPFromUnsigned(4u)), RTMultiply(a, d)));
            RetainedTriple square = RTSqrt(discriminant);
            RetainedTriple q = RTMultiply(
                RTConstant(RPFromBits(0xbf000000u)),
                RTAdd(b, RPScalarGreaterEqual(b.hi, RPZero()) ? square : RTNegate(square)));
            if (RTExact(q, RPZero()))
                first = RTDivide(RTNegate(b), RTMultiply(RTConstant(RPFromUnsigned(2u)), a));
            else {
                first = RTDivide(q, a);
                second = RTDivide(d, q);
            }
        }
        for (uint r = 0; r < 2; ++r) {
            RetainedTriple candidate = r == 0 ? first : second;
            if (candidate.valid == 0) continue;
            RetainedTriple denominator = RTConstant(RPZero());
            RPScalar denominatorScale = RPZero();
            for (uint i = 0; i < 4; ++i) {
                RetainedTriple term =
                    RTMultiply(Metric(rowBase, c, i), i == c ? candidate : Tangent(rowBase, i));
                denominator = RTAdd(denominator, term);
                denominatorScale = RPUpAdd(denominatorScale, RTAbsoluteUpper(term));
            }
            // A central root with a singular differentiated constraint cannot
            // own the coupled candidate. Try the represented spatial roots.
            if (denominator.valid == 0 ||
                RPScalarLessEqual(RTAbsoluteLower(denominator),
                                  RPUpMultiply(RPFromBits(0x29800000u), denominatorScale)))
                continue;
            RetainedTriple difference = RTSubtract(candidate, Tangent(rowBase, c));
            RetainedTriple original = RTCenter(Tangent(rowBase, c));
            RetainedTriple absolute =
                RTCenter(RPScalarLess(difference.hi, RPZero()) ? RTNegate(difference) : difference);
            RetainedTriple scaleValue =
                RTAdd(RTConstant(RPOne()),
                      RPScalarLess(original.hi, RPZero()) ? RTNegate(original) : original);
            RetainedTriple correction = RTDivide(absolute, scaleValue);
            if (correction.valid != 0 &&
                (best.valid == 0 || RPScalarLess(RTSubtract(correction, best).hi, RPZero()))) {
                best = correction;
                component = c;
                root = candidate;
            }
        }
    }
    return component < 4 && root.valid != 0;
}
#endif

#ifdef SIRIUS_RETAINED_PORTABLE
[shader("compute")][numthreads(SIRIUS_RETAINED_EXECUTION_LANES, 1, 1)] void ComputeMain(uint3 id
                                                                              : SV_GroupID,
                                                                                uint lane
                                                                              : SV_GroupIndex) {
    uint ni, si, no, so;
    inputs.GetDimensions(ni, si);
    outputs.GetDimensions(no, so);
    if (ni < 1 || id.y != 0 || id.z != 0) return;
    uint rows = inputs[0], row = id.x;
    if (rows == 0 || rows > 65536 || row >= rows || ni < 1 + 225 * rows + 103 || no % rows != 0)
        return;
    uint programBase = 1 + 225 * rows, count = inputs[programBase],
         registers = inputs[programBase + 1];
    uint rowWords = no / rows, rowBase = row * rowWords;
    if (rowWords < 769) return;
    for (uint i = lane; i < rowWords; i += SIRIUS_RETAINED_EXECUTION_LANES) outputs[rowBase + i] = 0;
#ifndef SIRIUS_RETAINED_PORTABLE
    AllMemoryBarrierWithGroupSync();
#endif
    if (inputs[programBase + 2] != 100 || count > 65536 || registers != SIRIUS_RETAINED_REGISTERS ||
        ni != programBase + 103 + 5 * count + SIRIUS_RETAINED_LAYERS + 1 ||
        rowWords != 769 + 5 * registers)
        return;
    for (uint i = 0; i < 45; ++i) {
        uint address = 1 + 225 * row + 5 * i;
        RetainedTriple value =
            RTPack(RPFromBits(inputs[address]), RPFromBits(inputs[address + 1]),
                   RPFromBits(inputs[address + 2]), RPFromBits(inputs[address + 3]));
        if (inputs[address + 4] != 1 || value.valid == 0) return;
        if (lane == 0) WriteStored(rowBase + 504 + 5 * i, value);
    }
#ifndef SIRIUS_RETAINED_PORTABLE
    AllMemoryBarrierWithGroupSync();
#endif
    RetainedTriple chart = ReadStored(rowBase + 504 + 5 * 44);
    if (!RTExact(chart, RPOne()) && !RTExact(chart, RPScalarNegate(RPOne()))) return;
    if (!RPScalarEqual(RTAbsoluteUpper(ReadStored(rowBase + 504 + 5)), RPZero()) &&
        !RPScalarEqual(RTAbsoluteUpper(ReadStored(rowBase + 504 + 15)), RPZero()))
        return;
    for (uint i = 45; i < 53; ++i)
        if (lane == 0) WriteStored(rowBase + 504 + 5 * i, RTConstant(i == 49 ? RPOne() : RPZero()));
#ifndef SIRIUS_RETAINED_PORTABLE
    AllMemoryBarrierWithGroupSync();
#endif
    bool flat = RTExact(ReadStored(rowBase + 504), RPZero()) &&
        RTExact(ReadStored(rowBase + 514), RPZero()) &&
        RTExact(ReadStored(rowBase + 519), RPZero());
    uint component = 4;
    RetainedTriple root = RTInvalid();
    // Preserve flat initialization or the curved metric prefix, projection,
    // and finally the curved continuation, sharing each large call body once.
    [loop] for (uint phase_index = 0; phase_index < 3; ++phase_index) {
        if (phase_index == 1) {
            if (!Project(rowBase, component, root)) return;
            if (flat) {
                if (component != 0) return;
                for (uint i = 0; i < 40; ++i) {
                    RetainedTriple phase = ReadStored(rowBase + 504 + 5 * (4 + i));
                    RetainedTriple physical = i % 8 == 4 ? RTNegate(phase) : phase;
                    if (i == 4) {
                        physical = root;
                        phase = RTNegate(root);
                    }
                    if (i >= 8 && i % 8 == 4) {
                        RetainedTriple constraint = RTConstant(RPZero());
                        for (uint axis = 1; axis < 4; ++axis)
                            constraint = RTAdd(constraint,
                                               RTMultiply(Tangent(rowBase, axis),
                                                          ReadStored(rowBase + 504 + 5 * (4 + i + axis))));
                        physical = RTDivide(constraint, root);
                        phase = RTNegate(physical);
                    }
                    if (phase.valid == 0 || physical.valid == 0) return;
                    if (lane == 0) WriteStored(rowBase + 104 + i * 5, phase);
                    if (lane == 0) WriteStored(rowBase + 304 + i * 5, physical);
                }
                if (lane == 0) outputs[rowBase] = 1;
                if (lane == 0) outputs[rowBase + 1] = 0;
                return;
            }
            for (uint i = 0; i < 4; ++i) {
                if (lane == 0)
                    WriteStored(rowBase + 504 + 5 * (45 + i), i == component ? root : Tangent(rowBase, i));
                if (lane == 0)
                    WriteStored(rowBase + 504 + 5 * (49 + i),
                                RTConstant(i == component ? RPOne() : RPZero()));
            }
        } else if (flat) {
            if (lane != 0) return;
            for (uint i = 0; i < 4; ++i) {
                for (uint j = 0; j < 4; ++j)
                    if (lane == 0)
                        WriteStored(rowBase + 4 + 5 * (i * 4 + j),
                                    RTConstant(i == j ? (i == 0 ? RPScalarNegate(RPOne()) : RPOne())
                                                      : RPZero()));
                RetainedTriple p = ReadStored(rowBase + 504 + 5 * (8 + i));
                if (lane == 0) WriteStored(rowBase + 4 + 5 * (16 + i), i == 0 ? RTNegate(p) : p);
            }
        } else if (!EvaluateProgram(row, rowBase, programBase, count, registers,
                                    phase_index == 0 ? 20 : 100, lane)) return;
    }
    if (lane == 0) outputs[rowBase] = 1;
    if (lane == 0) outputs[rowBase + 1] = component;
}
#else
[shader("compute")][numthreads(SIRIUS_RETAINED_EXECUTION_LANES, 1, 1)] void ComputeMain(uint3 id
                                                                              : SV_GroupID,
                                                                                uint lane
                                                                              : SV_GroupIndex) {
    uint ni, si, no, so;
    inputs.GetDimensions(ni, si);
    outputs.GetDimensions(no, so);
    if (ni < 1 || id.y != 0 || id.z != 0) return;
    uint rows = inputs[0], row = id.x;
    if (rows == 0 || rows > 65536 || row >= rows || ni < 1 + 225 * rows + 103 || no % rows != 0)
        return;
    uint programBase = 1 + 225 * rows, count = inputs[programBase],
         registers = inputs[programBase + 1];
    uint rowWords = no / rows, rowBase = row * rowWords;
    if (rowWords < 769) return;
    for (uint i = lane; i < rowWords; i += SIRIUS_RETAINED_EXECUTION_LANES) outputs[rowBase + i] = 0;
#ifndef SIRIUS_RETAINED_PORTABLE
    AllMemoryBarrierWithGroupSync();
#endif
    if (inputs[programBase + 2] != 100 || count > 65536 || registers != SIRIUS_RETAINED_REGISTERS ||
        ni != programBase + 103 + 5 * count + SIRIUS_RETAINED_LAYERS + 1 ||
        rowWords != 769 + 5 * registers)
        return;
    for (uint i = 0; i < 45; ++i) {
        uint address = 1 + 225 * row + 5 * i;
        RetainedTriple value =
            RTPack(RPFromBits(inputs[address]), RPFromBits(inputs[address + 1]),
                   RPFromBits(inputs[address + 2]), RPFromBits(inputs[address + 3]));
        if (inputs[address + 4] != 1 || value.valid == 0) return;
        if (lane == 0) WriteStored(rowBase + 504 + 5 * i, value);
    }
#ifndef SIRIUS_RETAINED_PORTABLE
    AllMemoryBarrierWithGroupSync();
#endif
    RetainedTriple chart = ReadStored(rowBase + 504 + 5 * 44);
    if (!RTExact(chart, RPOne()) && !RTExact(chart, RPScalarNegate(RPOne()))) return;
    if (!RPScalarEqual(RTAbsoluteUpper(ReadStored(rowBase + 504 + 5)), RPZero()) &&
        !RPScalarEqual(RTAbsoluteUpper(ReadStored(rowBase + 504 + 15)), RPZero()))
        return;
    for (uint i = 45; i < 53; ++i)
        if (lane == 0) WriteStored(rowBase + 504 + 5 * i, RTConstant(i == 49 ? RPOne() : RPZero()));
#ifndef SIRIUS_RETAINED_PORTABLE
    AllMemoryBarrierWithGroupSync();
#endif
    if (RTExact(ReadStored(rowBase + 504), RPZero()) &&
        RTExact(ReadStored(rowBase + 514), RPZero()) &&
        RTExact(ReadStored(rowBase + 519), RPZero())) {
        if (lane != 0) return;
        for (uint i = 0; i < 4; ++i) {
            for (uint j = 0; j < 4; ++j)
                if (lane == 0)
                    WriteStored(rowBase + 4 + 5 * (i * 4 + j),
                                RTConstant(i == j ? (i == 0 ? RPScalarNegate(RPOne()) : RPOne())
                                                  : RPZero()));
            RetainedTriple p = ReadStored(rowBase + 504 + 5 * (8 + i));
            if (lane == 0) WriteStored(rowBase + 4 + 5 * (16 + i), i == 0 ? RTNegate(p) : p);
        }
        uint selected;
        RetainedTriple temporal;
        if (!Project(rowBase, selected, temporal) || selected != 0) return;
        for (uint i = 0; i < 40; ++i) {
            RetainedTriple phase = ReadStored(rowBase + 504 + 5 * (4 + i));
            RetainedTriple physical = i % 8 == 4 ? RTNegate(phase) : phase;
            if (i == 4) {
                physical = temporal;
                phase = RTNegate(temporal);
            }
            if (i >= 8 && i % 8 == 4) {
                RetainedTriple constraint = RTConstant(RPZero());
                for (uint axis = 1; axis < 4; ++axis)
                    constraint = RTAdd(constraint,
                                       RTMultiply(Tangent(rowBase, axis),
                                                  ReadStored(rowBase + 504 + 5 * (4 + i + axis))));
                physical = RTDivide(constraint, temporal);
                phase = RTNegate(physical);
            }
            if (phase.valid == 0 || physical.valid == 0) return;
            if (lane == 0) WriteStored(rowBase + 104 + i * 5, phase);
            if (lane == 0) WriteStored(rowBase + 304 + i * 5, physical);
        }
        if (lane == 0) outputs[rowBase] = 1;
        if (lane == 0) outputs[rowBase + 1] = 0;
        return;
    }
    if (!EvaluateProgram(row, rowBase, programBase, count, registers, 20, lane)) return;
    uint component;
    RetainedTriple root;
    if (!Project(rowBase, component, root)) return;
    for (uint i = 0; i < 4; ++i) {
        if (lane == 0)
            WriteStored(rowBase + 504 + 5 * (45 + i), i == component ? root : Tangent(rowBase, i));
        if (lane == 0)
            WriteStored(rowBase + 504 + 5 * (49 + i),
                        RTConstant(i == component ? RPOne() : RPZero()));
    }
#ifndef SIRIUS_RETAINED_PORTABLE
    AllMemoryBarrierWithGroupSync();
#endif
    if (!EvaluateProgram(row, rowBase, programBase, count, registers, 100, lane)) return;
    if (lane == 0) outputs[rowBase] = 1;
    if (lane == 0) outputs[rowBase + 1] = component;
}
#endif
