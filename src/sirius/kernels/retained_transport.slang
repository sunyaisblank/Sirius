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

bool EvaluateOne(uint index, uint row, uint rowBase, uint programBase, uint registers, uint op) {
    uint address = programBase + 43 + 5 * index,
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
groupshared uint retainedHomogeneousLayer[SIRIUS_RETAINED_LAYERS];
void ClassifyLayers(uint programBase, uint count, uint layers, uint lane) {
    uint layerBase = programBase + 43 + 5 * count;
    for (uint layer = lane; layer < layers; layer += SIRIUS_RETAINED_EXECUTION_LANES) {
        uint start = inputs[layerBase + layer], end = inputs[layerBase + layer + 1];
        uint homogeneous = 0;
        if (start < end && end <= count && end - start <= SIRIUS_RETAINED_LANES) {
            uint opcode = inputs[programBase + 43 + 5 * start];
            homogeneous = 1;
            for (uint node = start + 1; node < end; ++node)
                if (inputs[programBase + 43 + 5 * node] != opcode) homogeneous = 0;
        }
        retainedHomogeneousLayer[layer] = homogeneous;
    }
}
bool EvaluateNode(uint index, uint row, uint rowBase, uint programBase, uint registers,
                  uint layerOpcode, uint homogeneous) {
    if (homogeneous != 0)
        return EvaluateOne(index, row, rowBase, programBase, registers, layerOpcode);
    return EvaluateOne(index, row, rowBase, programBase, registers,
                       inputs[programBase + 43 + 5 * index]);
}
bool EvaluateProgram(uint row, uint rowBase, uint programBase, uint count, uint registers,
                     uint layers, uint stage, uint lane) {
    uint needed = count;
    if (needed == 0 || needed > count || layers == 0 || layers > SIRIUS_RETAINED_LAYERS) return false;
    if (lane == 0) retainedStatus = 0;
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
    GroupMemoryBarrierWithGroupSync();
#endif
    uint layerBase = programBase + 43 + 5 * count;
    if (inputs[layerBase] != 0 || inputs[layerBase + layers] != count) return false;
    uint expected = 0;
    for (uint layer = 0; layer < layers; ++layer) {
        uint start = inputs[layerBase + layer], end = inputs[layerBase + layer + 1];
        if (start >= needed) break;
        if (start != expected || end <= start || end > needed ||
            end - start > SIRIUS_RETAINED_LANES)
            RejectRetainedProgram();
#if defined(SIRIUS_RETAINED_PORTABLE) && !defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
        else {
            uint layerOpcode = inputs[programBase + 43 + 5 * start];
            // Layer destinations are disjoint from every input register and
            // from one another; the generated plan releases inputs afterwards.
            for (uint node = start; node < end; ++node)
                if (!EvaluateNode(node, row, rowBase, programBase, registers, layerOpcode, retainedHomogeneousLayer[layer])) RejectRetainedProgram();
        }
#elif defined(SIRIUS_RETAINED_PORTABLE) && defined(SIRIUS_RETAINED_NORMAL_SUM32) && defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
        else {
            uint layerOpcode = inputs[programBase + 43 + 5 * start];
            // Cover the original complete layer before its publication barrier.
            for (uint node = start + lane; node < end;
                 node += SIRIUS_RETAINED_EXECUTION_LANES)
                if (!EvaluateNode(node, row, rowBase, programBase, registers, layerOpcode, retainedHomogeneousLayer[layer]))
                    RejectRetainedProgram();
        }
#else
        else if (start + lane < end &&
                 !EvaluateNode(start + lane, row, rowBase, programBase, registers,
                               inputs[programBase + 43 + 5 * start], retainedHomogeneousLayer[layer]))
            RejectRetainedProgram();
#endif
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
        GroupMemoryBarrierWithGroupSync();
#endif
        if (retainedStatus != 0) return false;
        expected = end;
    }

#if defined(SIRIUS_RETAINED_PORTABLE) && !defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
    for (uint i = 0; i < 40; ++i) {
        uint index = inputs[programBase + 3 + i];
        if (index >= registers || ReadValue(index, rowBase).valid == 0) return false;
    }
    if (lane == 0) {
        for (uint i = 0; i < 40; ++i)
            WriteStored(rowBase + 1004 + stage * 200 + i * 5,
                        ReadValue(inputs[programBase + 3 + i], rowBase));
    }
#else
    // Every root must pass before any stage value is published.
    for (uint i = lane; i < 40; i += SIRIUS_RETAINED_EXECUTION_LANES) {
        uint index = inputs[programBase + 3 + i];
        if (index >= registers)
            RejectRetainedProgram();
        else if (ScratchRead(index, 4, rowBase, 2404) == 0)
            RejectRetainedProgram();
    }
    GroupMemoryBarrierWithGroupSync();
    if (retainedStatus != 0) return false;
    for (uint i = lane; i < 40; i += SIRIUS_RETAINED_EXECUTION_LANES)
        WriteStored(rowBase + 1004 + stage * 200 + i * 5,
                    ReadValue(inputs[programBase + 3 + i], rowBase));
#endif
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
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
// Six row sums, fourteen nonzero corrections and five embedded-error weights.
// Preserve the complete retained value, including its integer validity word.
groupshared uint retainedCoefficients[SIRIUS_RETAINED_COEFFICIENTS * 5];
void CacheCoefficient(uint lane) {
    if (lane >= SIRIUS_RETAINED_COEFFICIENTS) return;
#ifdef SIRIUS_RETAINED_PORTABLE
    // Complete integer-arithmetic quotient words; retain the radius and status.
    // The kernel builder checks every exact rational enclosure before compile.
    uint base = lane * 5;
    if (lane == 0u) {
        retainedCoefficients[base + 0u] = 0x3e4ccccdu;
        retainedCoefficients[base + 1u] = 0xb14ccccdu;
        retainedCoefficients[base + 2u] = 0x244ccccdu;
        retainedCoefficients[base + 3u] = 0x19800000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 1u) {
        retainedCoefficients[base + 0u] = 0x3e99999au;
        retainedCoefficients[base + 1u] = 0xb24ccccdu;
        retainedCoefficients[base + 2u] = 0x254ccccdu;
        retainedCoefficients[base + 3u] = 0x1a800000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 2u) {
        retainedCoefficients[base + 0u] = 0x3f4ccccdu;
        retainedCoefficients[base + 1u] = 0xb24ccccdu;
        retainedCoefficients[base + 2u] = 0x254ccccdu;
        retainedCoefficients[base + 3u] = 0x1a800000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 3u) {
        retainedCoefficients[base + 0u] = 0x3f638e39u;
        retainedCoefficients[base + 1u] = 0xb1e38e39u;
        retainedCoefficients[base + 2u] = 0x24638e39u;
        retainedCoefficients[base + 3u] = 0x19800000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 4u) {
        retainedCoefficients[base + 0u] = 0x3f800000u;
        retainedCoefficients[base + 1u] = 0x00000000u;
        retainedCoefficients[base + 2u] = 0x00000000u;
        retainedCoefficients[base + 3u] = 0x00000000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 5u) {
        retainedCoefficients[base + 0u] = 0x3f800000u;
        retainedCoefficients[base + 1u] = 0x00000000u;
        retainedCoefficients[base + 2u] = 0x00000000u;
        retainedCoefficients[base + 3u] = 0x00000000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 6u) {
        retainedCoefficients[base + 0u] = 0x3e666666u;
        retainedCoefficients[base + 1u] = 0x31cccccdu;
        retainedCoefficients[base + 2u] = 0xa4cccccdu;
        retainedCoefficients[base + 3u] = 0x1a000000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 7u) {
        retainedCoefficients[base + 0u] = 0xc06eeeefu;
        retainedCoefficients[base + 1u] = 0x32888889u;
        retainedCoefficients[base + 2u] = 0xa66eeeefu;
        retainedCoefficients[base + 3u] = 0x1b800000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 8u) {
        retainedCoefficients[base + 0u] = 0x40638e39u;
        retainedCoefficients[base + 1u] = 0xb2e38e39u;
        retainedCoefficients[base + 2u] = 0x25638e39u;
        retainedCoefficients[base + 3u] = 0x1a800000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 9u) {
        retainedCoefficients[base + 0u] = 0xc139885fu;
        retainedCoefficients[base + 1u] = 0x34d6adf2u;
        retainedCoefficients[base + 2u] = 0xa80fbe73u;
        retainedCoefficients[base + 3u] = 0x1d000000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 10u) {
        retainedCoefficients[base + 0u] = 0x411d2a92u;
        retainedCoefficients[base + 1u] = 0xb49c3ad7u;
        retainedCoefficients[base + 2u] = 0x2589584du;
        retainedCoefficients[base + 3u] = 0x1b000000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 11u) {
        retainedCoefficients[base + 0u] = 0xbe94e4f6u;
        retainedCoefficients[base + 1u] = 0x31bc3997u;
        retainedCoefficients[base + 2u] = 0x25347fd3u;
        retainedCoefficients[base + 3u] = 0x1a800000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 12u) {
        retainedCoefficients[base + 0u] = 0xc12c1f08u;
        retainedCoefficients[base + 1u] = 0x34783e10u;
        retainedCoefficients[base + 2u] = 0xa7f83e10u;
        retainedCoefficients[base + 3u] = 0x1d000000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 13u) {
        retainedCoefficients[base + 0u] = 0x410e80b5u;
        retainedCoefficients[base + 1u] = 0x33dc87a3u;
        retainedCoefficients[base + 2u] = 0xa7127789u;
        retainedCoefficients[base + 3u] = 0x1c000000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 14u) {
        retainedCoefficients[base + 0u] = 0x3e8e8ba3u;
        retainedCoefficients[base + 1u] = 0xb13a2e8cu;
        retainedCoefficients[base + 2u] = 0x24ba2e8cu;
        retainedCoefficients[base + 3u] = 0x1a000000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 15u) {
        retainedCoefficients[base + 0u] = 0xbe8c0c4cu;
        retainedCoefficients[base + 1u] = 0x3272633du;
        retainedCoefficients[base + 2u] = 0x24138a52u;
        retainedCoefficients[base + 3u] = 0x19000000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 16u) {
        retainedCoefficients[base + 0u] = 0x3ee6024du;
        retainedCoefficients[base + 1u] = 0xb1b53f41u;
        retainedCoefficients[base + 2u] = 0x253c25b9u;
        retainedCoefficients[base + 3u] = 0x1a000000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 17u) {
        retainedCoefficients[base + 0u] = 0x3f26aaabu;
        retainedCoefficients[base + 1u] = 0xb2aaaaabu;
        retainedCoefficients[base + 2u] = 0x262aaaabu;
        retainedCoefficients[base + 3u] = 0x1b800000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 18u) {
        retainedCoefficients[base + 0u] = 0xbea50e7eu;
        retainedCoefficients[base + 1u] = 0x3254873fu;
        retainedCoefficients[base + 2u] = 0xa554873fu;
        retainedCoefficients[base + 3u] = 0x1a800000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 19u) {
        retainedCoefficients[base + 0u] = 0x3e061862u;
        retainedCoefficients[base + 1u] = 0xb1f3cf3du;
        retainedCoefficients[base + 2u] = 0x23c30c31u;
        retainedCoefficients[base + 3u] = 0x19000000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 20u) {
        retainedCoefficients[base + 0u] = 0xbb8b5ad3u;
        retainedCoefficients[base + 1u] = 0x2f2f5404u;
        retainedCoefficients[base + 2u] = 0xa2845ccdu;
        retainedCoefficients[base + 3u] = 0x17800000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 21u) {
        retainedCoefficients[base + 0u] = 0x3d177777u;
        retainedCoefficients[base + 1u] = 0x30eeeeefu;
        retainedCoefficients[base + 2u] = 0xa3088889u;
        retainedCoefficients[base + 3u] = 0x18800000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 22u) {
        retainedCoefficients[base + 0u] = 0xbd50568fu;
        retainedCoefficients[base + 1u] = 0x30d35e7au;
        retainedCoefficients[base + 2u] = 0x247cb729u;
        retainedCoefficients[base + 3u] = 0x19800000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 23u) {
        retainedCoefficients[base + 0u] = 0x3d2ba454u;
        retainedCoefficients[base + 1u] = 0xb004a1e4u;
        retainedCoefficients[base + 2u] = 0x239055d2u;
        retainedCoefficients[base + 3u] = 0x18800000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
    if (lane == 24u) {
        retainedCoefficients[base + 0u] = 0xbccccccdu;
        retainedCoefficients[base + 1u] = 0x2fcccccdu;
        retainedCoefficients[base + 2u] = 0xa2cccccdu;
        retainedCoefficients[base + 3u] = 0x18000000u;
        retainedCoefficients[base + 4u] = 0x00000001u;
        return;
    }
#else
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
#endif
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
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
    AllMemoryBarrierWithGroupSync();
#endif
    uint layers = SIRIUS_RETAINED_LAYERS;
    if (inputs[programBase + 2] != 40 || count > 65536 || registers != SIRIUS_RETAINED_REGISTERS ||
        43 + 5 * count + layers + 1 != SIRIUS_RETAINED_GENERAL_PROGRAM_WORDS ||
        ni != programBase + SIRIUS_RETAINED_GENERAL_PROGRAM_WORDS +
                            SIRIUS_RETAINED_SCHWARZSCHILD_PROGRAM_WORDS ||
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
    // The parameters are immutable for the complete seven-stage trial. Every
    // lane in this row chooses the same table, so all existing joins remain
    // uniform. Nonzero centres or parameter radii retain the general family.
    if (RTExact(spin, RPZero()) && RTExact(ReadOriginal(2, row), RPZero()) &&
        RTExact(lambda, RPZero())) {
        programBase += SIRIUS_RETAINED_GENERAL_PROGRAM_WORDS;
        count = inputs[programBase];
        layers = SIRIUS_RETAINED_SCHWARZSCHILD_LAYERS;
        if (inputs[programBase + 1] != registers || inputs[programBase + 2] != 40 ||
            count > 65536 ||
            43 + 5 * count + layers + 1 != SIRIUS_RETAINED_SCHWARZSCHILD_PROGRAM_WORDS)
            return;
    }
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
    // The immutable plan leaves never-assigned slots invalid across all stages.
    // The existing stage-zero barrier publishes this disjoint validity clear.
    for (uint index = lane; index < SIRIUS_RETAINED_REGISTERS;
         index += SIRIUS_RETAINED_EXECUTION_LANES)
        ScratchWrite(index, 4, rowBase, 2404, 0);
#endif
    ClassifyLayers(programBase, count, layers, lane);
    for (uint stage = 0; stage < 7; ++stage) {
        if (lane == 0) retainedStatus = 0;
#if defined(SIRIUS_RETAINED_PORTABLE) && !defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
        // One invocation owns this row; cache all complete coefficients after
        // stage zero succeeds, before their first use. No inter-invocation join.
        if (stage == 1)
            for (uint index = 0; index < SIRIUS_RETAINED_COEFFICIENTS; ++index)
                CacheCoefficient(index);
#else
        // Stage zero has succeeded. Publish once at the existing stage barrier;
        // expression register reuse never writes this separate immutable cache.
#if defined(SIRIUS_RETAINED_PORTABLE) && defined(SIRIUS_RETAINED_NORMAL_SUM32) && defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
        if (stage == 1)
            for (uint index = lane; index < SIRIUS_RETAINED_COEFFICIENTS;
                 index += SIRIUS_RETAINED_EXECUTION_LANES)
                CacheCoefficient(index);
#else
        if (stage == 1) CacheCoefficient(lane);
#endif
        GroupMemoryBarrierWithGroupSync();
#endif
        for (uint i = lane; i < 40; i += SIRIUS_RETAINED_EXECUTION_LANES) {
            RetainedTriple value = ReadOriginal(4 + i, row);
            if (stage != 0) {
                RetainedTriple first = ReadStored(rowBase + 1004 + i * 5);
                RetainedTriple weighted = RTMultiply(first, ReadCoefficient(stage - 1));
                for (uint j = 1; j < stage; ++j)
                    if (an[stage][j] != 0) {
                        RetainedTriple difference =
                            RTSubtract(ReadStored(rowBase + 1004 + j * 200 + i * 5), first);
                        uint index = 6 + (stage - 2) * (stage - 1) / 2 + j - 1;
                        if (stage == 6) --index; // The preceding correction is zero.
                        weighted = RTAdd(weighted, RTMultiply(difference, ReadCoefficient(index)));
                    }
                RetainedTriple increment = RTMultiply(weighted, h);
                value = RTAdd(value, increment);
                if (stage == 6) WriteStored(rowBase + 4 + 80 * 5 + i * 5, increment);
            }
            if (value.valid == 0) RejectRetainedProgram();
            WriteStored(rowBase + 804 + i * 5, value);
        }
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
        AllMemoryBarrierWithGroupSync();
#endif
        if (retainedStatus != 0) return;
        if (lane == 0) outputs[rowBase + 1] = stage + 1;
        if (!EvaluateProgram(row, rowBase, programBase, count, registers, layers, stage, lane)) return;
    }
    if (lane == 0) retainedStatus = 0;
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
    GroupMemoryBarrierWithGroupSync();
#endif
    for (uint i = lane; i < 40; i += SIRIUS_RETAINED_EXECUTION_LANES) {
        RetainedTriple first = ReadStored(rowBase + 1004 + i * 5), weighted = RTConstant(RPZero());
        for (uint j = 1; j < 7; ++j)
            if (en[j] != 0) {
                RetainedTriple difference =
                    RTSubtract(ReadStored(rowBase + 1004 + j * 200 + i * 5), first);
                weighted = RTAdd(weighted, RTMultiply(difference, ReadCoefficient(18 + j)));
            }
        RetainedTriple error = RTMultiply(weighted, h), full = ReadStored(rowBase + 804 + i * 5);
        RetainedTriple lower = RTSubtract(full, error);
        if (error.valid == 0 || lower.valid == 0) RejectRetainedProgram();
        WriteStored(rowBase + 4 + i * 5, full);
        WriteStored(rowBase + 4 + 40 * 5 + i * 5, lower);
        WriteStored(rowBase + 4 + 120 * 5 + i * 5, error);
    }
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
    AllMemoryBarrierWithGroupSync();
#endif
    if (retainedStatus != 0) return;
    if (lane == 0) outputs[rowBase] = 1;
    // Unlike the analytic flat path, this row owns seven complete RHS records.
    if (lane == 0) outputs[rowBase + 2] = 1;
}
