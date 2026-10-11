#if !defined(SIRIUS_RETAINED_PORTABLE) || !defined(SIRIUS_RETAINED_NORMAL_SUM32) || !defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT) || SIRIUS_RETAINED_EXECUTION_LANES != 16
#error This compiler-only diagnostic requires the existing portable normal Transport layout.
#endif
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
bool EvaluateProgramGeneric(uint row, uint rowBase, uint programBase, uint count, uint registers,
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
            // Layer destinations are disjoint from every input register and
            // from one another; the generated plan releases inputs afterwards.
            for (uint node = start; node < end; ++node)
                if (!EvaluateOne(node, row, rowBase, programBase, registers)) RejectRetainedProgram();
        }
#elif defined(SIRIUS_RETAINED_PORTABLE) && defined(SIRIUS_RETAINED_NORMAL_SUM32) && defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
        else {
            // Cover the original complete layer before its publication barrier.
            for (uint node = start + lane; node < end;
                 node += SIRIUS_RETAINED_EXECUTION_LANES)
                if (!EvaluateOne(node, row, rowBase, programBase, registers))
                    RejectRetainedProgram();
        }
#else
        else if (start + lane < end &&
                 !EvaluateOne(start + lane, row, rowBase, programBase, registers))
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
static uint LiteralGeneralTable[13096] = {
2587u, 459u, 40u, 14u, 21u, 2u, 60u, 4u, 413u, 253u, 158u, 145u, 412u, 300u, 343u, 41u,
231u, 182u, 142u, 372u, 280u, 232u, 324u, 41u, 385u, 209u, 258u, 196u, 178u, 92u, 289u, 41u,
282u, 243u, 377u, 387u, 143u, 216u, 394u, 41u, 210u, 165u, 6u, 1u, 0u, 0u, 0u, 0u,
1u, 1u, 1u, 0u, 0u, 1u, 2u, 2u, 0u, 0u, 1u, 3u, 3u, 0u, 0u, 1u,
4u, 5u, 0u, 0u, 1u, 5u, 6u, 0u, 0u, 1u, 6u, 7u, 0u, 0u, 1u, 7u,
8u, 0u, 0u, 1u, 8u, 9u, 0u, 0u, 1u, 9u, 10u, 0u, 0u, 1u, 10u, 11u,
0u, 0u, 1u, 11u, 13u, 0u, 0u, 1u, 12u, 14u, 0u, 0u, 1u, 13u, 15u, 0u,
0u, 1u, 14u, 16u, 0u, 0u, 1u, 15u, 17u, 0u, 0u, 1u, 16u, 18u, 0u, 0u,
1u, 17u, 19u, 0u, 0u, 1u, 18u, 21u, 0u, 0u, 1u, 19u, 22u, 0u, 0u, 1u,
20u, 23u, 0u, 0u, 1u, 21u, 24u, 0u, 0u, 1u, 22u, 25u, 0u, 0u, 1u, 23u,
26u, 0u, 0u, 1u, 24u, 27u, 0u, 0u, 1u, 25u, 29u, 0u, 0u, 1u, 26u, 30u,
0u, 0u, 1u, 27u, 31u, 0u, 0u, 1u, 28u, 32u, 0u, 0u, 1u, 29u, 33u, 0u,
0u, 1u, 30u, 34u, 0u, 0u, 1u, 31u, 35u, 0u, 0u, 1u, 32u, 37u, 0u, 0u,
1u, 33u, 38u, 0u, 0u, 1u, 34u, 39u, 0u, 0u, 1u, 35u, 40u, 0u, 0u, 1u,
36u, 41u, 0u, 0u, 1u, 37u, 42u, 0u, 0u, 1u, 38u, 43u, 0u, 0u, 1u, 39u,
44u, 0u, 0u, 0u, 40u, 1065353216u, 0u, 0u, 0u, 41u, 0u, 0u, 0u, 0u, 42u, 1082130432u,
0u, 0u, 0u, 43u, 1073741824u, 0u, 0u, 0u, 44u, 1056964608u, 0u, 0u, 0u, 45u, 1077936128u, 0u,
0u, 0u, 46u, 3212836864u, 0u, 0u, 4u, 47u, 4u, 40u, 0u, 4u, 48u, 5u, 39u, 0u,
4u, 49u, 6u, 40u, 0u, 4u, 50u, 1u, 1u, 0u, 4u, 51u, 0u, 43u, 0u, 4u,
52u, 2u, 2u, 0u, 4u, 53u, 40u, 39u, 0u, 4u, 54u, 7u, 46u, 0u, 4u, 55u,
8u, 40u, 0u, 4u, 56u, 9u, 40u, 0u, 4u, 57u, 10u, 40u, 0u, 4u, 58u, 40u,
40u, 0u, 4u, 59u, 14u, 46u, 0u, 4u, 60u, 11u, 40u, 0u, 4u, 61u, 12u, 39u,
0u, 4u, 62u, 13u, 40u, 0u, 4u, 63u, 15u, 40u, 0u, 4u, 64u, 16u, 40u, 0u,
4u, 65u, 17u, 40u, 0u, 4u, 66u, 21u, 46u, 0u, 4u, 67u, 18u, 40u, 0u, 4u,
68u, 19u, 39u, 0u, 4u, 69u, 20u, 40u, 0u, 4u, 70u, 22u, 40u, 0u, 4u, 71u,
23u, 40u, 0u, 4u, 72u, 24u, 40u, 0u, 4u, 73u, 28u, 46u, 0u, 4u, 74u, 25u,
40u, 0u, 4u, 75u, 26u, 39u, 0u, 4u, 76u, 27u, 40u, 0u, 4u, 77u, 29u, 40u,
0u, 4u, 78u, 30u, 40u, 0u, 4u, 79u, 31u, 40u, 0u, 4u, 80u, 35u, 46u, 0u,
4u, 81u, 32u, 40u, 0u, 4u, 82u, 33u, 39u, 0u, 4u, 83u, 34u, 40u, 0u, 4u,
84u, 36u, 40u, 0u, 4u, 85u, 37u, 40u, 0u, 4u, 86u, 38u, 40u, 0u, 4u, 34u,
1u, 48u, 0u, 4u, 33u, 1u, 47u, 0u, 4u, 32u, 53u, 46u, 0u, 4u, 27u, 1u,
58u, 0u, 4u, 26u, 1u, 53u, 0u, 4u, 25u, 1u, 61u, 0u, 4u, 20u, 1u, 60u,
0u, 4u, 19u, 1u, 68u, 0u, 4u, 18u, 1u, 67u, 0u, 4u, 13u, 1u, 75u, 0u,
4u, 12u, 1u, 74u, 0u, 4u, 11u, 1u, 82u, 0u, 4u, 6u, 1u, 81u, 0u, 7u,
46u, 33u, 0u, 0u, 7u, 5u, 52u, 0u, 0u, 7u, 4u, 41u, 0u, 0u, 7u, 87u,
27u, 0u, 0u, 7u, 88u, 20u, 0u, 0u, 7u, 89u, 18u, 0u, 0u, 7u, 90u, 12u,
0u, 0u, 7u, 91u, 6u, 0u, 0u, 4u, 6u, 32u, 7u, 0u, 4u, 12u, 32u, 14u,
0u, 4u, 18u, 32u, 21u, 0u, 4u, 20u, 32u, 28u, 0u, 4u, 27u, 32u, 35u, 0u,
9u, 35u, 47u, 48u, 0u, 9u, 28u, 49u, 1u, 0u, 9u, 21u, 0u, 1u, 0u, 9u,
0u, 35u, 28u, 0u, 9u, 14u, 35u, 49u, 0u, 2u, 28u, 41u, 6u, 0u, 4u, 6u,
3u, 14u, 0u, 8u, 14u, 0u, 0u, 0u, 5u, 0u, 47u, 14u, 0u, 5u, 35u, 48u,
14u, 0u, 5u, 7u, 49u, 14u, 0u, 5u, 33u, 1u, 14u, 0u, 5u, 52u, 58u, 14u,
0u, 5u, 92u, 53u, 14u, 0u, 5u, 93u, 60u, 14u, 0u, 5u, 94u, 61u, 14u, 0u,
5u, 95u, 62u, 14u, 0u, 5u, 96u, 67u, 14u, 0u, 5u, 97u, 68u, 14u, 0u, 5u,
98u, 69u, 14u, 0u, 5u, 99u, 74u, 14u, 0u, 5u, 100u, 75u, 14u, 0u, 5u, 101u,
76u, 14u, 0u, 5u, 102u, 81u, 14u, 0u, 5u, 103u, 82u, 14u, 0u, 5u, 104u, 83u,
14u, 0u, 4u, 1u, 33u, 33u, 0u, 4u, 105u, 7u, 7u, 0u, 4u, 106u, 0u, 0u,
0u, 4u, 107u, 35u, 35u, 0u, 4u, 108u, 52u, 0u, 0u, 4u, 109u, 0u, 52u, 0u,
4u, 110u, 92u, 35u, 0u, 4u, 111u, 35u, 92u, 0u, 4u, 112u, 52u, 7u, 0u, 4u,
113u, 7u, 52u, 0u, 4u, 114u, 93u, 0u, 0u, 4u, 115u, 0u, 93u, 0u, 4u, 116u,
94u, 35u, 0u, 4u, 117u, 35u, 94u, 0u, 4u, 118u, 95u, 7u, 0u, 4u, 119u, 7u,
95u, 0u, 4u, 120u, 52u, 93u, 0u, 4u, 121u, 93u, 52u, 0u, 4u, 122u, 92u, 94u,
0u, 4u, 123u, 94u, 92u, 0u, 4u, 124u, 52u, 95u, 0u, 4u, 125u, 95u, 52u, 0u,
4u, 126u, 96u, 0u, 0u, 4u, 127u, 0u, 96u, 0u, 4u, 128u, 97u, 35u, 0u, 4u,
129u, 35u, 97u, 0u, 4u, 130u, 98u, 7u, 0u, 4u, 131u, 7u, 98u, 0u, 4u, 132u,
52u, 96u, 0u, 4u, 133u, 96u, 52u, 0u, 4u, 134u, 92u, 97u, 0u, 4u, 135u, 97u,
92u, 0u, 4u, 136u, 52u, 98u, 0u, 4u, 137u, 98u, 52u, 0u, 4u, 138u, 99u, 0u,
0u, 4u, 139u, 0u, 99u, 0u, 4u, 140u, 100u, 35u, 0u, 4u, 141u, 35u, 100u, 0u,
4u, 142u, 101u, 7u, 0u, 4u, 143u, 7u, 101u, 0u, 4u, 144u, 52u, 99u, 0u, 4u,
145u, 99u, 52u, 0u, 4u, 146u, 92u, 100u, 0u, 4u, 147u, 100u, 92u, 0u, 4u, 148u,
52u, 101u, 0u, 4u, 149u, 101u, 52u, 0u, 4u, 150u, 102u, 0u, 0u, 4u, 151u, 0u,
102u, 0u, 4u, 152u, 103u, 35u, 0u, 4u, 153u, 35u, 103u, 0u, 4u, 154u, 104u, 7u,
0u, 4u, 155u, 7u, 104u, 0u, 4u, 156u, 52u, 102u, 0u, 4u, 157u, 102u, 52u, 0u,
4u, 158u, 92u, 103u, 0u, 4u, 159u, 103u, 92u, 0u, 4u, 160u, 52u, 104u, 0u, 4u,
161u, 104u, 52u, 0u, 2u, 104u, 106u, 107u, 0u, 2u, 103u, 108u, 109u, 0u, 2u, 102u,
110u, 111u, 0u, 2u, 101u, 112u, 113u, 0u, 2u, 100u, 114u, 115u, 0u, 2u, 99u, 116u,
117u, 0u, 2u, 98u, 118u, 119u, 0u, 2u, 97u, 120u, 121u, 0u, 2u, 96u, 122u, 123u,
0u, 2u, 95u, 124u, 125u, 0u, 2u, 94u, 126u, 127u, 0u, 2u, 93u, 128u, 129u, 0u,
2u, 92u, 130u, 131u, 0u, 2u, 52u, 132u, 133u, 0u, 2u, 33u, 134u, 135u, 0u, 2u,
7u, 136u, 137u, 0u, 2u, 35u, 138u, 139u, 0u, 2u, 0u, 140u, 141u, 0u, 2u, 162u,
142u, 143u, 0u, 2u, 163u, 144u, 145u, 0u, 2u, 164u, 146u, 147u, 0u, 2u, 165u, 148u,
149u, 0u, 2u, 166u, 150u, 151u, 0u, 2u, 167u, 152u, 153u, 0u, 2u, 168u, 154u, 155u,
0u, 2u, 169u, 156u, 157u, 0u, 2u, 170u, 158u, 159u, 0u, 2u, 171u, 160u, 161u, 0u,
7u, 161u, 1u, 0u, 0u, 7u, 160u, 103u, 0u, 0u, 7u, 159u, 102u, 0u, 0u, 7u,
158u, 101u, 0u, 0u, 7u, 157u, 97u, 0u, 0u, 7u, 156u, 96u, 0u, 0u, 7u, 155u,
95u, 0u, 0u, 7u, 154u, 52u, 0u, 0u, 7u, 153u, 33u, 0u, 0u, 7u, 152u, 7u,
0u, 0u, 7u, 151u, 163u, 0u, 0u, 7u, 150u, 164u, 0u, 0u, 7u, 149u, 165u, 0u,
0u, 7u, 148u, 169u, 0u, 0u, 7u, 147u, 170u, 0u, 0u, 7u, 146u, 171u, 0u, 0u,
2u, 145u, 104u, 105u, 0u, 2u, 144u, 100u, 99u, 0u, 2u, 143u, 94u, 93u, 0u, 2u,
142u, 35u, 0u, 0u, 2u, 141u, 166u, 167u, 0u, 2u, 167u, 145u, 161u, 0u, 2u, 166u,
144u, 98u, 0u, 2u, 0u, 143u, 92u, 0u, 2u, 35u, 142u, 162u, 0u, 2u, 93u, 141u,
168u, 0u, 4u, 141u, 167u, 167u, 0u, 4u, 142u, 1u, 42u, 0u, 4u, 143u, 1u, 43u,
0u, 4u, 144u, 103u, 167u, 0u, 4u, 145u, 167u, 103u, 0u, 4u, 161u, 102u, 167u, 0u,
4u, 94u, 167u, 102u, 0u, 4u, 99u, 101u, 167u, 0u, 4u, 100u, 167u, 101u, 0u, 4u,
104u, 166u, 167u, 0u, 4u, 140u, 167u, 166u, 0u, 4u, 139u, 97u, 167u, 0u, 4u, 138u,
103u, 166u, 0u, 4u, 137u, 166u, 103u, 0u, 4u, 136u, 167u, 97u, 0u, 4u, 135u, 96u,
167u, 0u, 4u, 134u, 102u, 166u, 0u, 4u, 133u, 166u, 102u, 0u, 4u, 132u, 167u, 96u,
0u, 4u, 131u, 95u, 167u, 0u, 4u, 130u, 101u, 166u, 0u, 4u, 129u, 166u, 101u, 0u,
4u, 128u, 167u, 95u, 0u, 4u, 127u, 0u, 167u, 0u, 4u, 126u, 167u, 0u, 0u, 4u,
125u, 52u, 167u, 0u, 4u, 124u, 103u, 0u, 0u, 4u, 123u, 0u, 103u, 0u, 4u, 122u,
167u, 52u, 0u, 4u, 121u, 33u, 167u, 0u, 4u, 120u, 102u, 0u, 0u, 4u, 119u, 0u,
102u, 0u, 4u, 118u, 167u, 33u, 0u, 4u, 117u, 7u, 167u, 0u, 4u, 116u, 101u, 0u,
0u, 4u, 115u, 0u, 101u, 0u, 4u, 114u, 167u, 7u, 0u, 4u, 113u, 35u, 167u, 0u,
4u, 112u, 167u, 35u, 0u, 4u, 111u, 163u, 167u, 0u, 4u, 110u, 103u, 35u, 0u, 4u,
109u, 35u, 103u, 0u, 4u, 108u, 167u, 163u, 0u, 4u, 107u, 164u, 167u, 0u, 4u, 106u,
102u, 35u, 0u, 4u, 172u, 35u, 102u, 0u, 4u, 173u, 167u, 164u, 0u, 4u, 174u, 165u,
167u, 0u, 4u, 175u, 101u, 35u, 0u, 4u, 176u, 35u, 101u, 0u, 4u, 177u, 167u, 165u,
0u, 4u, 178u, 93u, 167u, 0u, 4u, 179u, 167u, 93u, 0u, 4u, 180u, 169u, 167u, 0u,
4u, 181u, 103u, 93u, 0u, 4u, 182u, 93u, 103u, 0u, 4u, 183u, 167u, 169u, 0u, 4u,
184u, 170u, 167u, 0u, 4u, 185u, 102u, 93u, 0u, 4u, 186u, 93u, 102u, 0u, 4u, 187u,
167u, 170u, 0u, 4u, 188u, 171u, 167u, 0u, 4u, 189u, 101u, 93u, 0u, 4u, 190u, 93u,
101u, 0u, 2u, 1u, 144u, 145u, 0u, 2u, 42u, 161u, 94u, 0u, 2u, 191u, 99u, 100u,
0u, 2u, 192u, 104u, 140u, 0u, 2u, 193u, 139u, 138u, 0u, 2u, 194u, 137u, 136u, 0u,
2u, 195u, 135u, 134u, 0u, 2u, 196u, 133u, 132u, 0u, 2u, 197u, 131u, 130u, 0u, 2u,
198u, 129u, 128u, 0u, 2u, 199u, 127u, 126u, 0u, 2u, 200u, 125u, 124u, 0u, 2u, 201u,
123u, 122u, 0u, 2u, 202u, 121u, 120u, 0u, 2u, 203u, 119u, 118u, 0u, 2u, 204u, 117u,
116u, 0u, 2u, 205u, 115u, 114u, 0u, 2u, 206u, 113u, 112u, 0u, 2u, 207u, 111u, 110u,
0u, 2u, 208u, 109u, 108u, 0u, 2u, 209u, 107u, 106u, 0u, 2u, 210u, 172u, 173u, 0u,
2u, 211u, 174u, 175u, 0u, 2u, 212u, 176u, 177u, 0u, 2u, 213u, 178u, 179u, 0u, 2u,
214u, 180u, 181u, 0u, 2u, 215u, 182u, 183u, 0u, 2u, 216u, 184u, 185u, 0u, 2u, 217u,
186u, 187u, 0u, 2u, 218u, 188u, 189u, 0u, 4u, 189u, 142u, 105u, 0u, 4u, 188u, 143u,
105u, 0u, 4u, 187u, 142u, 101u, 0u, 4u, 186u, 143u, 101u, 0u, 4u, 185u, 142u, 98u,
0u, 4u, 184u, 143u, 98u, 0u, 4u, 183u, 142u, 95u, 0u, 4u, 182u, 143u, 95u, 0u,
4u, 181u, 142u, 92u, 0u, 4u, 180u, 143u, 92u, 0u, 4u, 179u, 142u, 7u, 0u, 4u,
178u, 143u, 7u, 0u, 4u, 177u, 142u, 162u, 0u, 4u, 176u, 143u, 162u, 0u, 4u, 175u,
142u, 165u, 0u, 4u, 174u, 143u, 165u, 0u, 4u, 173u, 142u, 168u, 0u, 4u, 172u, 143u,
168u, 0u, 4u, 106u, 167u, 171u, 0u, 4u, 107u, 142u, 171u, 0u, 4u, 108u, 143u, 171u,
0u, 2u, 143u, 141u, 189u, 0u, 2u, 142u, 191u, 187u, 0u, 2u, 168u, 192u, 185u, 0u,
2u, 162u, 193u, 194u, 0u, 2u, 92u, 195u, 196u, 0u, 2u, 98u, 197u, 198u, 0u, 2u,
105u, 199u, 181u, 0u, 2u, 109u, 200u, 201u, 0u, 2u, 110u, 202u, 203u, 0u, 2u, 111u,
204u, 205u, 0u, 2u, 112u, 206u, 177u, 0u, 2u, 113u, 207u, 208u, 0u, 2u, 114u, 209u,
210u, 0u, 2u, 115u, 211u, 212u, 0u, 2u, 116u, 213u, 173u, 0u, 2u, 117u, 214u, 215u,
0u, 2u, 118u, 216u, 217u, 0u, 2u, 119u, 190u, 106u, 0u, 7u, 106u, 167u, 0u, 0u,
7u, 173u, 166u, 0u, 0u, 7u, 177u, 0u, 0u, 0u, 7u, 181u, 35u, 0u, 0u, 7u,
185u, 93u, 0u, 0u, 2u, 187u, 98u, 183u, 0u, 2u, 189u, 111u, 179u, 0u, 2u, 217u,
115u, 175u, 0u, 2u, 216u, 218u, 119u, 0u, 2u, 119u, 216u, 107u, 0u, 6u, 216u, 143u,
0u, 0u, 2u, 143u, 167u, 216u, 0u, 2u, 107u, 216u, 106u, 0u, 4u, 106u, 43u, 216u,
0u, 4u, 115u, 143u, 44u, 0u, 5u, 143u, 188u, 107u, 0u, 5u, 216u, 1u, 106u, 0u,
5u, 111u, 42u, 106u, 0u, 5u, 98u, 142u, 106u, 0u, 5u, 175u, 168u, 106u, 0u, 5u,
179u, 105u, 106u, 0u, 5u, 183u, 112u, 106u, 0u, 5u, 218u, 116u, 106u, 0u, 2u, 116u,
103u, 216u, 0u, 2u, 112u, 216u, 160u, 0u, 2u, 105u, 102u, 111u, 0u, 2u, 168u, 111u,
159u, 0u, 2u, 142u, 101u, 98u, 0u, 2u, 188u, 98u, 158u, 0u, 2u, 42u, 166u, 175u,
0u, 2u, 1u, 175u, 173u, 0u, 2u, 215u, 0u, 179u, 0u, 2u, 214u, 179u, 177u, 0u,
2u, 213u, 35u, 183u, 0u, 2u, 212u, 183u, 181u, 0u, 2u, 211u, 93u, 218u, 0u, 2u,
210u, 218u, 185u, 0u, 4u, 185u, 116u, 44u, 0u, 4u, 181u, 143u, 112u, 0u, 4u, 177u,
105u, 44u, 0u, 4u, 173u, 143u, 168u, 0u, 4u, 93u, 142u, 44u, 0u, 4u, 35u, 143u,
188u, 0u, 4u, 0u, 42u, 44u, 0u, 4u, 166u, 143u, 1u, 0u, 4u, 158u, 43u, 175u,
0u, 4u, 159u, 215u, 44u, 0u, 4u, 160u, 143u, 214u, 0u, 4u, 101u, 43u, 179u, 0u,
4u, 102u, 213u, 44u, 0u, 4u, 103u, 143u, 212u, 0u, 4u, 209u, 43u, 183u, 0u, 4u,
208u, 211u, 44u, 0u, 4u, 207u, 143u, 210u, 0u, 4u, 206u, 43u, 218u, 0u, 4u, 211u,
216u, 158u, 0u, 4u, 213u, 111u, 158u, 0u, 4u, 215u, 98u, 158u, 0u, 4u, 42u, 216u,
101u, 0u, 4u, 142u, 111u, 101u, 0u, 4u, 105u, 98u, 101u, 0u, 4u, 116u, 216u, 209u,
0u, 4u, 218u, 111u, 209u, 0u, 4u, 183u, 98u, 209u, 0u, 4u, 179u, 216u, 206u, 0u,
4u, 175u, 111u, 206u, 0u, 4u, 205u, 98u, 206u, 0u, 7u, 206u, 181u, 0u, 0u, 7u,
209u, 173u, 0u, 0u, 7u, 101u, 35u, 0u, 0u, 7u, 158u, 166u, 0u, 0u, 7u, 98u,
211u, 0u, 0u, 7u, 111u, 213u, 0u, 0u, 7u, 216u, 215u, 0u, 0u, 7u, 204u, 160u,
0u, 0u, 7u, 203u, 42u, 0u, 0u, 7u, 202u, 142u, 0u, 0u, 7u, 201u, 105u, 0u,
0u, 7u, 200u, 103u, 0u, 0u, 7u, 199u, 116u, 0u, 0u, 7u, 198u, 218u, 0u, 0u,
7u, 197u, 183u, 0u, 0u, 7u, 196u, 207u, 0u, 0u, 7u, 195u, 179u, 0u, 0u, 7u,
194u, 175u, 0u, 0u, 7u, 193u, 205u, 0u, 0u, 2u, 205u, 186u, 101u, 0u, 2u, 175u,
184u, 158u, 0u, 2u, 179u, 162u, 98u, 0u, 2u, 183u, 92u, 111u, 0u, 2u, 218u, 187u,
216u, 0u, 2u, 116u, 180u, 204u, 0u, 2u, 105u, 109u, 203u, 0u, 2u, 142u, 110u, 202u,
0u, 2u, 42u, 189u, 201u, 0u, 2u, 215u, 176u, 200u, 0u, 2u, 213u, 113u, 199u, 0u,
2u, 211u, 114u, 198u, 0u, 2u, 207u, 217u, 197u, 0u, 2u, 103u, 172u, 196u, 0u, 2u,
160u, 117u, 195u, 0u, 2u, 166u, 118u, 194u, 0u, 2u, 35u, 119u, 193u, 0u, 5u, 193u,
206u, 107u, 0u, 5u, 194u, 209u, 107u, 0u, 5u, 195u, 205u, 107u, 0u, 5u, 196u, 175u,
107u, 0u, 5u, 197u, 179u, 106u, 0u, 5u, 198u, 183u, 106u, 0u, 5u, 199u, 218u, 106u,
0u, 5u, 200u, 116u, 107u, 0u, 5u, 201u, 105u, 106u, 0u, 5u, 202u, 142u, 106u, 0u,
5u, 203u, 42u, 106u, 0u, 5u, 204u, 215u, 107u, 0u, 5u, 216u, 213u, 106u, 0u, 5u,
111u, 211u, 106u, 0u, 5u, 98u, 207u, 106u, 0u, 5u, 158u, 103u, 107u, 0u, 5u, 101u,
160u, 106u, 0u, 5u, 119u, 166u, 106u, 0u, 5u, 217u, 35u, 106u, 0u, 2u, 35u, 97u,
197u, 0u, 2u, 166u, 197u, 157u, 0u, 2u, 160u, 96u, 198u, 0u, 2u, 103u, 198u, 156u,
0u, 2u, 207u, 95u, 199u, 0u, 2u, 211u, 199u, 155u, 0u, 2u, 213u, 52u, 201u, 0u,
2u, 215u, 201u, 154u, 0u, 2u, 42u, 33u, 202u, 0u, 2u, 142u, 202u, 153u, 0u, 2u,
105u, 7u, 203u, 0u, 2u, 116u, 203u, 152u, 0u, 2u, 218u, 163u, 216u, 0u, 2u, 183u,
216u, 151u, 0u, 2u, 179u, 164u, 111u, 0u, 2u, 175u, 111u, 150u, 0u, 2u, 205u, 165u,
98u, 0u, 2u, 209u, 98u, 149u, 0u, 2u, 206u, 169u, 101u, 0u, 2u, 106u, 101u, 148u,
0u, 2u, 189u, 170u, 119u, 0u, 2u, 187u, 119u, 147u, 0u, 2u, 118u, 171u, 217u, 0u,
2u, 117u, 217u, 146u, 0u, 4u, 217u, 35u, 44u, 0u, 4u, 119u, 196u, 112u, 0u, 4u,
101u, 143u, 166u, 0u, 4u, 98u, 193u, 1u, 0u, 4u, 111u, 160u, 44u, 0u, 4u, 216u,
196u, 168u, 0u, 4u, 203u, 143u, 103u, 0u, 4u, 202u, 194u, 1u, 0u, 4u, 201u, 207u,
44u, 0u, 4u, 199u, 196u, 188u, 0u, 4u, 198u, 143u, 211u, 0u, 4u, 197u, 195u, 1u,
0u, 4u, 146u, 213u, 44u, 0u, 4u, 147u, 200u, 112u, 0u, 4u, 148u, 143u, 215u, 0u,
4u, 149u, 193u, 214u, 0u, 4u, 150u, 42u, 44u, 0u, 4u, 151u, 200u, 168u, 0u, 4u,
152u, 143u, 142u, 0u, 4u, 153u, 194u, 214u, 0u, 4u, 154u, 105u, 44u, 0u, 4u, 155u,
200u, 188u, 0u, 4u, 156u, 143u, 116u, 0u, 4u, 157u, 195u, 214u, 0u, 4u, 171u, 218u,
44u, 0u, 4u, 170u, 204u, 112u, 0u, 4u, 169u, 143u, 183u, 0u, 4u, 165u, 193u, 212u,
0u, 4u, 164u, 179u, 44u, 0u, 4u, 163u, 204u, 168u, 0u, 4u, 7u, 143u, 175u, 0u,
4u, 33u, 194u, 212u, 0u, 4u, 52u, 205u, 44u, 0u, 4u, 95u, 204u, 188u, 0u, 4u,
96u, 143u, 209u, 0u, 4u, 97u, 195u, 212u, 0u, 4u, 114u, 206u, 44u, 0u, 4u, 113u,
158u, 112u, 0u, 4u, 110u, 143u, 106u, 0u, 4u, 109u, 193u, 210u, 0u, 4u, 92u, 189u,
44u, 0u, 4u, 162u, 158u, 168u, 0u, 4u, 172u, 143u, 187u, 0u, 4u, 176u, 194u, 210u,
0u, 4u, 180u, 118u, 44u, 0u, 4u, 184u, 158u, 188u, 0u, 4u, 186u, 143u, 117u, 0u,
4u, 173u, 195u, 210u, 0u, 2u, 117u, 119u, 101u, 0u, 2u, 118u, 216u, 203u, 0u, 2u,
187u, 199u, 198u, 0u, 2u, 189u, 147u, 148u, 0u, 2u, 106u, 151u, 152u, 0u, 2u, 206u,
155u, 156u, 0u, 2u, 209u, 170u, 169u, 0u, 2u, 205u, 163u, 7u, 0u, 2u, 175u, 95u,
96u, 0u, 2u, 179u, 113u, 110u, 0u, 2u, 183u, 162u, 172u, 0u, 2u, 218u, 184u, 186u,
0u, 7u, 186u, 117u, 0u, 0u, 7u, 184u, 98u, 0u, 0u, 7u, 172u, 118u, 0u, 0u,
7u, 162u, 202u, 0u, 0u, 7u, 110u, 187u, 0u, 0u, 7u, 113u, 197u, 0u, 0u, 7u,
96u, 189u, 0u, 0u, 7u, 95u, 149u, 0u, 0u, 7u, 7u, 106u, 0u, 0u, 7u, 163u,
153u, 0u, 0u, 7u, 169u, 206u, 0u, 0u, 7u, 170u, 157u, 0u, 0u, 7u, 156u, 209u,
0u, 0u, 7u, 155u, 165u, 0u, 0u, 7u, 152u, 205u, 0u, 0u, 7u, 151u, 33u, 0u,
0u, 7u, 148u, 175u, 0u, 0u, 7u, 147u, 97u, 0u, 0u, 7u, 198u, 179u, 0u, 0u,
7u, 199u, 109u, 0u, 0u, 7u, 203u, 183u, 0u, 0u, 7u, 216u, 176u, 0u, 0u, 7u,
101u, 218u, 0u, 0u, 7u, 119u, 173u, 0u, 0u, 2u, 218u, 186u, 184u, 0u, 2u, 183u,
172u, 162u, 0u, 2u, 179u, 182u, 110u, 0u, 2u, 175u, 96u, 95u, 0u, 2u, 205u, 7u,
163u, 0u, 2u, 209u, 178u, 169u, 0u, 2u, 206u, 156u, 155u, 0u, 2u, 106u, 152u, 151u,
0u, 2u, 189u, 174u, 148u, 0u, 2u, 187u, 198u, 199u, 0u, 2u, 118u, 203u, 216u, 0u,
2u, 117u, 108u, 101u, 0u, 5u, 101u, 218u, 107u, 0u, 5u, 216u, 183u, 107u, 0u, 5u,
203u, 175u, 107u, 0u, 5u, 199u, 205u, 107u, 0u, 5u, 198u, 206u, 107u, 0u, 5u, 148u,
106u, 107u, 0u, 5u, 151u, 187u, 107u, 0u, 5u, 152u, 118u, 107u, 0u, 10u, 118u, 115u,
143u, 167u, 10u, 187u, 185u, 193u, 167u, 10u, 106u, 177u, 194u, 167u, 10u, 206u, 93u, 195u,
167u, 10u, 205u, 0u, 196u, 167u, 10u, 175u, 217u, 101u, 167u, 10u, 183u, 111u, 216u, 167u,
10u, 218u, 159u, 200u, 167u, 10u, 155u, 146u, 203u, 167u, 10u, 156u, 150u, 199u, 167u, 10u,
169u, 102u, 204u, 167u, 10u, 163u, 171u, 198u, 167u, 10u, 7u, 164u, 148u, 167u, 10u, 95u,
208u, 158u, 167u, 10u, 96u, 114u, 151u, 167u, 10u, 110u, 92u, 152u, 167u, 2u, 152u, 179u,
113u, 0u, 2u, 151u, 209u, 170u, 0u, 2u, 148u, 189u, 147u, 0u, 2u, 198u, 117u, 119u,
0u, 5u, 117u, 152u, 107u, 0u, 5u, 189u, 151u, 107u, 0u, 5u, 209u, 148u, 107u, 0u,
5u, 179u, 198u, 107u, 0u, 10u, 198u, 201u, 117u, 167u, 10u, 148u, 154u, 189u, 167u, 10u,
151u, 52u, 209u, 167u, 10u, 152u, 180u, 179u, 167u, 6u, 179u, 118u, 0u, 0u, 4u, 118u,
43u, 179u, 0u, 4u, 209u, 14u, 179u, 0u, 4u, 179u, 209u, 209u, 0u, 4u, 189u, 209u,
47u, 0u, 4u, 117u, 209u, 48u, 0u, 4u, 180u, 51u, 209u, 0u, 4u, 52u, 3u, 209u,
0u, 4u, 154u, 209u, 58u, 0u, 4u, 201u, 209u, 53u, 0u, 4u, 167u, 209u, 60u, 0u,
4u, 107u, 209u, 61u, 0u, 4u, 119u, 209u, 67u, 0u, 4u, 147u, 209u, 68u, 0u, 4u,
170u, 209u, 74u, 0u, 4u, 113u, 209u, 75u, 0u, 4u, 199u, 209u, 81u, 0u, 4u, 203u,
209u, 82u, 0u, 5u, 216u, 49u, 209u, 0u, 5u, 101u, 187u, 118u, 0u, 5u, 92u, 106u,
118u, 0u, 5u, 114u, 206u, 118u, 0u, 5u, 164u, 205u, 118u, 0u, 5u, 171u, 218u, 118u,
0u, 5u, 150u, 169u, 118u, 0u, 5u, 146u, 95u, 118u, 0u, 4u, 95u, 50u, 216u, 0u,
4u, 169u, 52u, 209u, 0u, 4u, 218u, 216u, 40u, 0u, 4u, 205u, 14u, 101u, 0u, 4u,
206u, 14u, 92u, 0u, 4u, 106u, 14u, 114u, 0u, 4u, 187u, 14u, 164u, 0u, 4u, 49u,
43u, 164u, 0u, 4u, 111u, 14u, 171u, 0u, 4u, 217u, 43u, 171u, 0u, 4u, 158u, 14u,
150u, 0u, 4u, 204u, 43u, 150u, 0u, 4u, 200u, 14u, 146u, 0u, 4u, 196u, 43u, 146u,
0u, 4u, 146u, 95u, 216u, 0u, 4u, 150u, 218u, 40u, 0u, 4u, 171u, 51u, 205u, 0u,
4u, 164u, 205u, 209u, 0u, 4u, 195u, 209u, 205u, 0u, 4u, 194u, 216u, 205u, 0u, 4u,
193u, 3u, 205u, 0u, 4u, 208u, 52u, 205u, 0u, 4u, 102u, 205u, 47u, 0u, 4u, 159u,
205u, 48u, 0u, 4u, 0u, 51u, 206u, 0u, 4u, 93u, 206u, 209u, 0u, 4u, 177u, 209u,
206u, 0u, 4u, 185u, 216u, 206u, 0u, 4u, 143u, 3u, 206u, 0u, 4u, 115u, 52u, 206u,
0u, 4u, 162u, 206u, 47u, 0u, 4u, 172u, 206u, 48u, 0u, 4u, 184u, 51u, 106u, 0u,
4u, 186u, 106u, 209u, 0u, 4u, 108u, 209u, 106u, 0u, 4u, 174u, 216u, 106u, 0u, 4u,
178u, 3u, 106u, 0u, 4u, 182u, 52u, 106u, 0u, 4u, 173u, 106u, 47u, 0u, 4u, 176u,
106u, 48u, 0u, 4u, 109u, 51u, 187u, 0u, 4u, 97u, 187u, 209u, 0u, 4u, 33u, 209u,
187u, 0u, 4u, 165u, 216u, 187u, 0u, 4u, 157u, 3u, 187u, 0u, 4u, 153u, 52u, 187u,
0u, 4u, 149u, 187u, 47u, 0u, 4u, 197u, 187u, 48u, 0u, 4u, 202u, 101u, 49u, 0u,
4u, 98u, 205u, 187u, 0u, 4u, 116u, 187u, 205u, 0u, 4u, 105u, 205u, 60u, 0u, 4u,
142u, 187u, 58u, 0u, 4u, 42u, 205u, 61u, 0u, 4u, 215u, 92u, 49u, 0u, 4u, 213u,
206u, 187u, 0u, 4u, 211u, 187u, 206u, 0u, 4u, 207u, 206u, 60u, 0u, 4u, 103u, 206u,
61u, 0u, 4u, 160u, 187u, 53u, 0u, 4u, 166u, 114u, 49u, 0u, 4u, 35u, 106u, 187u,
0u, 4u, 210u, 187u, 106u, 0u, 4u, 212u, 106u, 60u, 0u, 4u, 214u, 106u, 61u, 0u,
4u, 1u, 51u, 111u, 0u, 4u, 188u, 111u, 209u, 0u, 4u, 168u, 209u, 111u, 0u, 4u,
112u, 216u, 111u, 0u, 4u, 44u, 3u, 111u, 0u, 4u, 181u, 52u, 111u, 0u, 4u, 192u,
111u, 47u, 0u, 4u, 191u, 111u, 48u, 0u, 4u, 190u, 101u, 217u, 0u, 4u, 141u, 205u,
111u, 0u, 4u, 120u, 111u, 205u, 0u, 4u, 121u, 205u, 67u, 0u, 4u, 122u, 111u, 58u,
0u, 4u, 49u, 150u, 10u, 0u, 4u, 218u, 193u, 209u, 0u, 4u, 61u, 143u, 209u, 0u,
4u, 60u, 178u, 209u, 0u, 4u, 123u, 157u, 209u, 0u, 4u, 124u, 150u, 17u, 0u, 4u,
125u, 193u, 187u, 0u, 4u, 126u, 157u, 205u, 0u, 4u, 127u, 143u, 187u, 0u, 4u, 128u,
157u, 206u, 0u, 4u, 129u, 178u, 187u, 0u, 4u, 130u, 157u, 106u, 0u, 4u, 131u, 44u,
209u, 0u, 4u, 132u, 150u, 24u, 0u, 4u, 133u, 193u, 111u, 0u, 4u, 134u, 44u, 205u,
0u, 4u, 135u, 205u, 68u, 0u, 4u, 136u, 92u, 217u, 0u, 4u, 137u, 206u, 111u, 0u,
4u, 138u, 111u, 206u, 0u, 4u, 139u, 143u, 111u, 0u, 4u, 140u, 44u, 206u, 0u, 4u,
104u, 206u, 67u, 0u, 4u, 100u, 206u, 68u, 0u, 4u, 99u, 111u, 53u, 0u, 4u, 94u,
114u, 217u, 0u, 4u, 161u, 106u, 111u, 0u, 4u, 145u, 111u, 106u, 0u, 4u, 144u, 178u,
111u, 0u, 4u, 219u, 44u, 106u, 0u, 4u, 220u, 106u, 67u, 0u, 4u, 221u, 106u, 68u,
0u, 4u, 222u, 51u, 158u, 0u, 4u, 223u, 158u, 209u, 0u, 4u, 224u, 209u, 158u, 0u,
4u, 225u, 216u, 158u, 0u, 4u, 226u, 3u, 158u, 0u, 4u, 227u, 52u, 158u, 0u, 4u,
228u, 158u, 47u, 0u, 4u, 229u, 158u, 48u, 0u, 4u, 230u, 150u, 31u, 0u, 4u, 231u,
101u, 204u, 0u, 4u, 232u, 205u, 158u, 0u, 4u, 233u, 158u, 205u, 0u, 4u, 234u, 193u,
158u, 0u, 4u, 235u, 205u, 74u, 0u, 4u, 236u, 158u, 58u, 0u, 4u, 237u, 205u, 75u,
0u, 4u, 238u, 92u, 204u, 0u, 4u, 239u, 206u, 158u, 0u, 4u, 240u, 158u, 206u, 0u,
4u, 241u, 143u, 158u, 0u, 4u, 242u, 206u, 74u, 0u, 4u, 243u, 206u, 75u, 0u, 4u,
244u, 158u, 53u, 0u, 4u, 245u, 114u, 204u, 0u, 4u, 246u, 106u, 158u, 0u, 4u, 247u,
158u, 106u, 0u, 4u, 248u, 178u, 158u, 0u, 4u, 249u, 106u, 74u, 0u, 4u, 250u, 106u,
75u, 0u, 4u, 251u, 51u, 200u, 0u, 4u, 252u, 200u, 209u, 0u, 4u, 253u, 209u, 200u,
0u, 4u, 44u, 226u, 209u, 0u, 4u, 157u, 226u, 205u, 0u, 4u, 204u, 226u, 206u, 0u,
4u, 217u, 226u, 106u, 0u, 4u, 75u, 216u, 200u, 0u, 4u, 74u, 3u, 200u, 0u, 4u,
68u, 52u, 200u, 0u, 4u, 67u, 200u, 47u, 0u, 4u, 254u, 200u, 48u, 0u, 4u, 255u,
150u, 38u, 0u, 4u, 256u, 101u, 196u, 0u, 4u, 257u, 205u, 200u, 0u, 4u, 258u, 200u,
205u, 0u, 4u, 259u, 193u, 200u, 0u, 4u, 260u, 205u, 81u, 0u, 4u, 261u, 200u, 58u,
0u, 4u, 262u, 205u, 82u, 0u, 4u, 263u, 92u, 196u, 0u, 4u, 264u, 206u, 200u, 0u,
4u, 265u, 200u, 206u, 0u, 4u, 266u, 143u, 200u, 0u, 4u, 267u, 206u, 81u, 0u, 4u,
268u, 206u, 82u, 0u, 4u, 269u, 200u, 53u, 0u, 4u, 270u, 114u, 196u, 0u, 4u, 271u,
106u, 200u, 0u, 4u, 272u, 200u, 106u, 0u, 4u, 273u, 178u, 200u, 0u, 4u, 274u, 106u,
81u, 0u, 4u, 275u, 106u, 82u, 0u, 2u, 226u, 179u, 50u, 0u, 2u, 178u, 189u, 34u,
0u, 2u, 143u, 117u, 46u, 0u, 2u, 193u, 180u, 5u, 0u, 2u, 196u, 179u, 146u, 0u,
2u, 114u, 164u, 195u, 0u, 2u, 92u, 218u, 208u, 0u, 2u, 101u, 102u, 154u, 0u, 2u,
82u, 159u, 87u, 0u, 2u, 81u, 93u, 177u, 0u, 2u, 53u, 61u, 115u, 0u, 2u, 276u,
162u, 26u, 0u, 2u, 277u, 172u, 201u, 0u, 2u, 278u, 186u, 108u, 0u, 2u, 279u, 60u,
182u, 0u, 2u, 280u, 97u, 33u, 0u, 2u, 281u, 123u, 153u, 0u, 2u, 282u, 149u, 167u,
0u, 2u, 283u, 197u, 107u, 0u, 2u, 284u, 188u, 168u, 0u, 2u, 285u, 131u, 181u, 0u,
2u, 286u, 192u, 119u, 0u, 2u, 287u, 191u, 147u, 0u, 2u, 288u, 223u, 224u, 0u, 2u,
289u, 44u, 227u, 0u, 2u, 290u, 228u, 170u, 0u, 2u, 291u, 229u, 113u, 0u, 2u, 292u,
252u, 253u, 0u, 2u, 293u, 67u, 199u, 0u, 2u, 294u, 254u, 203u, 0u, 7u, 254u, 194u,
0u, 0u, 7u, 67u, 185u, 0u, 0u, 7u, 44u, 174u, 0u, 0u, 7u, 253u, 165u, 0u,
0u, 7u, 252u, 202u, 0u, 0u, 7u, 229u, 215u, 0u, 0u, 7u, 228u, 166u, 0u, 0u,
7u, 227u, 112u, 0u, 0u, 7u, 224u, 190u, 0u, 0u, 7u, 223u, 136u, 0u, 0u, 7u,
131u, 94u, 0u, 0u, 7u, 123u, 225u, 0u, 0u, 7u, 60u, 231u, 0u, 0u, 7u, 61u,
238u, 0u, 0u, 7u, 218u, 245u, 0u, 0u, 7u, 191u, 75u, 0u, 0u, 7u, 192u, 256u,
0u, 0u, 7u, 181u, 263u, 0u, 0u, 7u, 168u, 270u, 0u, 0u, 2u, 270u, 58u, 44u,
0u, 2u, 263u, 62u, 253u, 0u, 2u, 256u, 282u, 25u, 0u, 2u, 75u, 283u, 88u, 0u,
2u, 245u, 175u, 252u, 0u, 2u, 238u, 183u, 229u, 0u, 2u, 231u, 198u, 228u, 0u, 2u,
225u, 69u, 227u, 0u, 2u, 94u, 286u, 19u, 0u, 2u, 136u, 287u, 89u, 0u, 2u, 190u,
155u, 224u, 0u, 2u, 112u, 156u, 223u, 0u, 2u, 166u, 148u, 131u, 0u, 2u, 215u, 76u,
123u, 0u, 2u, 202u, 290u, 13u, 0u, 2u, 165u, 291u, 90u, 0u, 2u, 174u, 163u, 60u,
0u, 2u, 185u, 7u, 61u, 0u, 2u, 194u, 151u, 218u, 0u, 2u, 188u, 83u, 191u, 0u,
2u, 197u, 293u, 11u, 0u, 2u, 149u, 294u, 91u, 0u, 2u, 153u, 96u, 192u, 0u, 2u,
33u, 110u, 181u, 0u, 2u, 97u, 152u, 168u, 0u, 5u, 168u, 178u, 226u, 0u, 5u, 181u,
143u, 226u, 0u, 5u, 192u, 193u, 196u, 0u, 5u, 191u, 169u, 45u, 0u, 5u, 218u, 254u,
209u, 0u, 5u, 61u, 92u, 45u, 0u, 5u, 60u, 67u, 209u, 0u, 5u, 123u, 53u, 45u,
0u, 5u, 131u, 270u, 209u, 0u, 5u, 223u, 279u, 45u, 0u, 5u, 224u, 263u, 209u, 0u,
5u, 227u, 281u, 45u, 0u, 5u, 228u, 245u, 118u, 0u, 5u, 229u, 238u, 118u, 0u, 5u,
252u, 231u, 118u, 0u, 5u, 253u, 225u, 209u, 0u, 5u, 44u, 285u, 45u, 0u, 5u, 294u,
190u, 118u, 0u, 5u, 293u, 112u, 118u, 0u, 5u, 291u, 166u, 118u, 0u, 5u, 290u, 215u,
209u, 0u, 5u, 287u, 289u, 45u, 0u, 5u, 286u, 174u, 118u, 0u, 5u, 283u, 185u, 118u,
0u, 5u, 282u, 194u, 118u, 0u, 5u, 152u, 188u, 209u, 0u, 5u, 151u, 153u, 118u, 0u,
5u, 148u, 33u, 118u, 0u, 5u, 198u, 97u, 118u, 0u, 4u, 97u, 168u, 40u, 0u, 4u,
33u, 181u, 39u, 0u, 4u, 153u, 50u, 218u, 0u, 4u, 188u, 95u, 218u, 0u, 4u, 194u,
168u, 114u, 0u, 4u, 185u, 181u, 114u, 0u, 4u, 174u, 218u, 40u, 0u, 4u, 215u, 50u,
60u, 0u, 4u, 166u, 95u, 60u, 0u, 4u, 112u, 168u, 81u, 0u, 4u, 190u, 181u, 81u,
0u, 4u, 225u, 60u, 40u, 0u, 4u, 231u, 50u, 131u, 0u, 4u, 238u, 95u, 131u, 0u,
4u, 245u, 168u, 278u, 0u, 4u, 263u, 181u, 278u, 0u, 4u, 270u, 131u, 40u, 0u, 4u,
67u, 50u, 224u, 0u, 4u, 254u, 95u, 224u, 0u, 4u, 289u, 168u, 280u, 0u, 4u, 285u,
181u, 280u, 0u, 4u, 281u, 224u, 40u, 0u, 4u, 279u, 14u, 228u, 0u, 4u, 53u, 224u,
205u, 0u, 4u, 92u, 218u, 187u, 0u, 4u, 193u, 14u, 229u, 0u, 4u, 143u, 224u, 206u,
0u, 4u, 178u, 60u, 187u, 0u, 4u, 169u, 14u, 252u, 0u, 4u, 118u, 224u, 106u, 0u,
4u, 110u, 131u, 187u, 0u, 4u, 96u, 50u, 253u, 0u, 4u, 7u, 95u, 253u, 0u, 4u,
163u, 168u, 284u, 0u, 4u, 156u, 181u, 284u, 0u, 4u, 155u, 253u, 40u, 0u, 4u, 183u,
14u, 294u, 0u, 4u, 175u, 253u, 205u, 0u, 4u, 91u, 218u, 111u, 0u, 4u, 90u, 14u,
293u, 0u, 4u, 89u, 253u, 206u, 0u, 4u, 88u, 60u, 111u, 0u, 4u, 11u, 14u, 291u,
0u, 4u, 13u, 253u, 106u, 0u, 4u, 19u, 131u, 111u, 0u, 4u, 25u, 50u, 290u, 0u,
4u, 83u, 95u, 290u, 0u, 4u, 76u, 168u, 288u, 0u, 4u, 69u, 181u, 288u, 0u, 4u,
62u, 290u, 40u, 0u, 4u, 58u, 14u, 286u, 0u, 4u, 182u, 290u, 205u, 0u, 4u, 108u,
218u, 158u, 0u, 4u, 186u, 14u, 283u, 0u, 4u, 172u, 290u, 206u, 0u, 4u, 162u, 60u,
158u, 0u, 4u, 115u, 14u, 282u, 0u, 4u, 177u, 290u, 106u, 0u, 4u, 93u, 131u, 158u,
0u, 4u, 159u, 50u, 152u, 0u, 4u, 102u, 95u, 152u, 0u, 4u, 208u, 74u, 209u, 0u,
4u, 195u, 168u, 292u, 0u, 4u, 164u, 181u, 292u, 0u, 4u, 282u, 97u, 40u, 0u, 4u,
283u, 33u, 40u, 0u, 4u, 286u, 153u, 216u, 0u, 4u, 291u, 174u, 40u, 0u, 4u, 293u,
215u, 216u, 0u, 4u, 294u, 225u, 40u, 0u, 4u, 252u, 231u, 216u, 0u, 4u, 229u, 270u,
40u, 0u, 4u, 228u, 67u, 216u, 0u, 4u, 158u, 281u, 40u, 0u, 4u, 111u, 51u, 279u,
0u, 4u, 187u, 279u, 209u, 0u, 4u, 146u, 209u, 279u, 0u, 4u, 203u, 216u, 279u, 0u,
4u, 199u, 153u, 224u, 0u, 4u, 113u, 67u, 218u, 0u, 4u, 170u, 3u, 279u, 0u, 4u,
147u, 52u, 279u, 0u, 4u, 119u, 279u, 47u, 0u, 4u, 107u, 279u, 48u, 0u, 4u, 167u,
51u, 193u, 0u, 4u, 201u, 193u, 209u, 0u, 4u, 154u, 209u, 193u, 0u, 4u, 180u, 216u,
193u, 0u, 4u, 117u, 215u, 224u, 0u, 4u, 189u, 67u, 60u, 0u, 4u, 179u, 3u, 193u,
0u, 4u, 87u, 52u, 193u, 0u, 4u, 5u, 193u, 47u, 0u, 4u, 46u, 193u, 48u, 0u,
4u, 26u, 51u, 169u, 0u, 4u, 34u, 169u, 209u, 0u, 4u, 295u, 209u, 169u, 0u, 4u,
296u, 216u, 169u, 0u, 4u, 297u, 231u, 224u, 0u, 4u, 298u, 67u, 131u, 0u, 4u, 299u,
3u, 169u, 0u, 4u, 300u, 52u, 169u, 0u, 4u, 301u, 169u, 47u, 0u, 4u, 302u, 169u,
48u, 0u, 4u, 303u, 96u, 216u, 0u, 4u, 304u, 155u, 40u, 0u, 4u, 305u, 51u, 183u,
0u, 4u, 306u, 183u, 209u, 0u, 4u, 307u, 209u, 183u, 0u, 4u, 308u, 216u, 183u, 0u,
4u, 309u, 153u, 253u, 0u, 4u, 310u, 96u, 218u, 0u, 4u, 311u, 3u, 183u, 0u, 4u,
312u, 52u, 183u, 0u, 4u, 313u, 183u, 47u, 0u, 4u, 314u, 183u, 48u, 0u, 4u, 315u,
51u, 90u, 0u, 4u, 316u, 90u, 209u, 0u, 4u, 317u, 209u, 90u, 0u, 4u, 318u, 216u,
90u, 0u, 4u, 319u, 215u, 253u, 0u, 4u, 320u, 96u, 60u, 0u, 4u, 321u, 3u, 90u,
0u, 4u, 322u, 52u, 90u, 0u, 4u, 323u, 90u, 47u, 0u, 4u, 324u, 90u, 48u, 0u,
4u, 325u, 51u, 11u, 0u, 4u, 326u, 11u, 209u, 0u, 4u, 90u, 282u, 8u, 0u, 4u,
183u, 283u, 9u, 0u, 4u, 155u, 291u, 10u, 0u, 4u, 169u, 294u, 10u, 0u, 4u, 193u,
229u, 10u, 0u, 4u, 279u, 282u, 15u, 0u, 4u, 281u, 283u, 16u, 0u, 4u, 67u, 158u,
10u, 0u, 4u, 270u, 170u, 209u, 0u, 4u, 225u, 291u, 17u, 0u, 4u, 174u, 179u, 209u,
0u, 4u, 33u, 294u, 17u, 0u, 4u, 97u, 299u, 209u, 0u, 4u, 224u, 229u, 17u, 0u,
4u, 327u, 282u, 22u, 0u, 4u, 328u, 283u, 23u, 0u, 4u, 329u, 304u, 10u, 0u, 4u,
330u, 311u, 209u, 0u, 4u, 331u, 291u, 24u, 0u, 4u, 332u, 321u, 209u, 0u, 4u, 333u,
294u, 24u, 0u, 4u, 334u, 209u, 11u, 0u, 4u, 335u, 216u, 11u, 0u, 4u, 336u, 231u,
253u, 0u, 4u, 337u, 96u, 131u, 0u, 4u, 338u, 3u, 11u, 0u, 4u, 339u, 52u, 11u,
0u, 4u, 340u, 11u, 47u, 0u, 4u, 341u, 11u, 48u, 0u, 4u, 342u, 229u, 24u, 0u,
4u, 343u, 25u, 216u, 0u, 4u, 344u, 282u, 29u, 0u, 4u, 345u, 283u, 30u, 0u, 4u,
346u, 62u, 40u, 0u, 4u, 347u, 51u, 58u, 0u, 4u, 348u, 58u, 209u, 0u, 4u, 349u,
209u, 58u, 0u, 4u, 350u, 216u, 58u, 0u, 4u, 351u, 153u, 290u, 0u, 4u, 352u, 25u,
218u, 0u, 4u, 353u, 3u, 58u, 0u, 4u, 354u, 52u, 58u, 0u, 4u, 355u, 58u, 47u,
0u, 4u, 356u, 58u, 48u, 0u, 4u, 357u, 291u, 31u, 0u, 4u, 358u, 51u, 186u, 0u,
4u, 359u, 186u, 209u, 0u, 4u, 360u, 209u, 186u, 0u, 4u, 361u, 216u, 186u, 0u, 4u,
362u, 215u, 290u, 0u, 4u, 363u, 25u, 60u, 0u, 4u, 364u, 3u, 186u, 0u, 4u, 365u,
52u, 186u, 0u, 4u, 366u, 186u, 47u, 0u, 4u, 367u, 186u, 48u, 0u, 4u, 368u, 294u,
31u, 0u, 4u, 369u, 51u, 115u, 0u, 4u, 370u, 115u, 209u, 0u, 4u, 371u, 209u, 115u,
0u, 4u, 372u, 216u, 115u, 0u, 4u, 373u, 231u, 290u, 0u, 4u, 374u, 25u, 131u, 0u,
4u, 375u, 3u, 115u, 0u, 4u, 376u, 52u, 115u, 0u, 2u, 321u, 192u, 191u, 0u, 2u,
311u, 28u, 90u, 0u, 2u, 299u, 286u, 188u, 0u, 2u, 179u, 293u, 166u, 0u, 2u, 170u,
252u, 238u, 0u, 2u, 186u, 228u, 254u, 0u, 2u, 58u, 67u, 124u, 0u, 2u, 62u, 187u,
98u, 0u, 2u, 25u, 116u, 146u, 0u, 2u, 11u, 53u, 203u, 0u, 2u, 96u, 270u, 125u,
0u, 2u, 290u, 126u, 147u, 0u, 2u, 253u, 119u, 105u, 0u, 2u, 24u, 107u, 42u, 0u,
2u, 17u, 201u, 213u, 0u, 2u, 377u, 211u, 154u, 0u, 2u, 378u, 143u, 180u, 0u, 2u,
379u, 174u, 127u, 0u, 2u, 380u, 128u, 87u, 0u, 2u, 381u, 5u, 207u, 0u, 2u, 382u,
46u, 103u, 0u, 2u, 383u, 34u, 35u, 0u, 2u, 384u, 210u, 295u, 0u, 2u, 385u, 118u,
296u, 0u, 2u, 386u, 97u, 129u, 0u, 2u, 387u, 130u, 300u, 0u, 2u, 388u, 301u, 212u,
0u, 2u, 389u, 302u, 214u, 0u, 2u, 390u, 303u, 7u, 0u, 2u, 391u, 329u, 132u, 0u,
2u, 392u, 306u, 141u, 0u, 2u, 393u, 120u, 307u, 0u, 2u, 394u, 175u, 308u, 0u, 2u,
395u, 330u, 133u, 0u, 2u, 396u, 134u, 312u, 0u, 2u, 397u, 313u, 121u, 0u, 2u, 398u,
314u, 135u, 0u, 2u, 399u, 316u, 137u, 0u, 2u, 400u, 138u, 317u, 0u, 2u, 401u, 89u,
318u, 0u, 2u, 402u, 332u, 139u, 0u, 2u, 403u, 140u, 322u, 0u, 2u, 404u, 323u, 104u,
0u, 2u, 405u, 324u, 100u, 0u, 2u, 406u, 326u, 161u, 0u, 2u, 407u, 145u, 334u, 0u,
2u, 408u, 13u, 335u, 0u, 2u, 409u, 219u, 339u, 0u, 2u, 410u, 340u, 220u, 0u, 2u,
411u, 341u, 221u, 0u, 2u, 412u, 343u, 83u, 0u, 2u, 413u, 348u, 232u, 0u, 2u, 414u,
233u, 349u, 0u, 2u, 415u, 182u, 350u, 0u, 2u, 416u, 157u, 354u, 0u, 2u, 417u, 355u,
235u, 0u, 2u, 418u, 356u, 237u, 0u, 2u, 419u, 359u, 239u, 0u, 2u, 420u, 240u, 360u,
0u, 2u, 421u, 172u, 361u, 0u, 2u, 422u, 204u, 365u, 0u, 2u, 423u, 366u, 242u, 0u,
2u, 424u, 367u, 243u, 0u, 2u, 425u, 370u, 246u, 0u, 4u, 370u, 321u, 32u, 0u, 4u,
367u, 321u, 282u, 0u, 4u, 366u, 321u, 283u, 0u, 4u, 365u, 321u, 150u, 0u, 4u, 361u,
321u, 158u, 0u, 4u, 360u, 321u, 304u, 0u, 4u, 359u, 338u, 209u, 0u, 4u, 356u, 346u,
10u, 0u, 4u, 355u, 321u, 346u, 0u, 4u, 354u, 353u, 209u, 0u, 4u, 350u, 364u, 209u,
0u, 4u, 349u, 375u, 209u, 0u, 4u, 348u, 115u, 47u, 0u, 4u, 343u, 115u, 48u, 0u,
4u, 341u, 229u, 31u, 0u, 4u, 340u, 159u, 216u, 0u, 4u, 339u, 282u, 36u, 0u, 4u,
335u, 283u, 37u, 0u, 4u, 334u, 152u, 40u, 0u, 4u, 332u, 14u, 151u, 0u, 4u, 330u,
152u, 205u, 0u, 4u, 329u, 218u, 200u, 0u, 4u, 97u, 153u, 152u, 0u, 4u, 174u, 159u,
218u, 0u, 4u, 270u, 74u, 205u, 0u, 4u, 67u, 291u, 38u, 0u, 4u, 90u, 14u, 148u,
0u, 4u, 326u, 152u, 206u, 0u, 4u, 324u, 60u, 200u, 0u, 4u, 323u, 215u, 152u, 0u,
4u, 322u, 159u, 60u, 0u, 4u, 318u, 74u, 206u, 0u, 4u, 317u, 294u, 38u, 0u, 4u,
316u, 14u, 198u, 0u, 4u, 314u, 152u, 106u, 0u, 4u, 313u, 131u, 200u, 0u, 4u, 312u,
231u, 152u, 0u, 4u, 308u, 159u, 131u, 0u, 4u, 307u, 74u, 106u, 0u, 4u, 306u, 229u,
38u, 0u, 2u, 375u, 311u, 183u, 0u, 2u, 364u, 114u, 299u, 0u, 2u, 353u, 81u, 179u,
0u, 2u, 346u, 278u, 170u, 0u, 2u, 338u, 280u, 186u, 0u, 2u, 304u, 62u, 25u, 0u,
2u, 158u, 96u, 290u, 0u, 2u, 229u, 253u, 142u, 0u, 2u, 294u, 17u, 377u, 0u, 2u,
291u, 379u, 380u, 0u, 2u, 159u, 382u, 160u, 0u, 2u, 115u, 383u, 384u, 0u, 2u, 231u,
386u, 387u, 0u, 2u, 215u, 284u, 390u, 0u, 2u, 153u, 392u, 393u, 0u, 2u, 198u, 395u,
396u, 0u, 2u, 148u, 397u, 122u, 0u, 2u, 151u, 399u, 400u, 0u, 2u, 152u, 402u, 403u,
0u, 2u, 131u, 405u, 99u, 0u, 2u, 60u, 406u, 407u, 0u, 2u, 218u, 359u, 144u, 0u,
2u, 74u, 288u, 412u, 0u, 2u, 200u, 356u, 230u, 0u, 2u, 106u, 413u, 414u, 0u, 2u,
206u, 354u, 234u, 0u, 2u, 205u, 417u, 236u, 0u, 2u, 14u, 419u, 420u, 0u, 2u, 38u,
350u, 241u, 0u, 2u, 31u, 424u, 244u, 0u, 2u, 303u, 247u, 371u, 0u, 2u, 302u, 177u,
372u, 0u, 2u, 301u, 349u, 248u, 0u, 2u, 300u, 217u, 376u, 0u, 2u, 296u, 348u, 249u,
0u, 2u, 295u, 343u, 250u, 0u, 2u, 34u, 340u, 102u, 0u, 2u, 46u, 208u, 68u, 0u,
4u, 340u, 192u, 364u, 0u, 4u, 343u, 192u, 353u, 0u, 4u, 348u, 192u, 346u, 0u, 4u,
349u, 192u, 338u, 0u, 4u, 350u, 168u, 304u, 0u, 4u, 354u, 181u, 304u, 0u, 4u, 356u,
168u, 294u, 0u, 4u, 359u, 181u, 294u, 0u, 4u, 424u, 168u, 115u, 0u, 4u, 420u, 181u,
115u, 0u, 4u, 419u, 192u, 215u, 0u, 4u, 417u, 168u, 153u, 0u, 4u, 414u, 181u, 153u,
0u, 4u, 413u, 168u, 151u, 0u, 4u, 412u, 181u, 151u, 0u, 4u, 407u, 168u, 60u, 0u,
4u, 406u, 181u, 60u, 0u, 4u, 405u, 192u, 74u, 0u, 4u, 403u, 168u, 106u, 0u, 4u,
402u, 181u, 106u, 0u, 4u, 400u, 168u, 14u, 0u, 4u, 399u, 181u, 14u, 0u, 4u, 397u,
334u, 40u, 0u, 4u, 396u, 51u, 332u, 0u, 4u, 395u, 332u, 209u, 0u, 4u, 393u, 209u,
332u, 0u, 4u, 392u, 216u, 332u, 0u, 4u, 390u, 3u, 332u, 0u, 4u, 387u, 52u, 332u,
0u, 4u, 386u, 332u, 47u, 0u, 4u, 384u, 332u, 48u, 0u, 4u, 383u, 51u, 90u, 0u,
4u, 382u, 90u, 209u, 0u, 4u, 380u, 209u, 90u, 0u, 4u, 379u, 216u, 90u, 0u, 4u,
377u, 3u, 90u, 0u, 4u, 17u, 52u, 90u, 0u, 4u, 253u, 90u, 47u, 0u, 4u, 290u,
90u, 48u, 0u, 4u, 96u, 51u, 316u, 0u, 4u, 25u, 316u, 209u, 0u, 4u, 62u, 209u,
316u, 0u, 4u, 186u, 216u, 316u, 0u, 4u, 170u, 3u, 316u, 0u, 4u, 179u, 52u, 316u,
0u, 4u, 299u, 316u, 47u, 0u, 4u, 311u, 316u, 48u, 0u, 7u, 316u, 340u, 0u, 0u,
7u, 90u, 194u, 0u, 0u, 7u, 332u, 185u, 0u, 0u, 7u, 334u, 343u, 0u, 0u, 7u,
52u, 112u, 0u, 0u, 7u, 51u, 190u, 0u, 0u, 7u, 48u, 348u, 0u, 0u, 7u, 47u,
245u, 0u, 0u, 7u, 3u, 263u, 0u, 0u, 7u, 376u, 349u, 0u, 0u, 7u, 372u, 289u,
0u, 0u, 7u, 371u, 285u, 0u, 0u, 7u, 183u, 11u, 0u, 0u, 7u, 208u, 92u, 0u,
0u, 7u, 102u, 378u, 0u, 0u, 7u, 177u, 178u, 0u, 0u, 7u, 68u, 385u, 0u, 0u,
7u, 217u, 110u, 0u, 0u, 7u, 250u, 419u, 0u, 0u, 7u, 249u, 163u, 0u, 0u, 7u,
248u, 156u, 0u, 0u, 7u, 247u, 394u, 0u, 0u, 7u, 244u, 91u, 0u, 0u, 7u, 241u,
401u, 0u, 0u, 7u, 236u, 88u, 0u, 0u, 7u, 234u, 408u, 0u, 0u, 7u, 230u, 19u,
0u, 0u, 7u, 144u, 405u, 0u, 0u, 7u, 99u, 76u, 0u, 0u, 7u, 122u, 69u, 0u,
0u, 7u, 160u, 415u, 0u, 0u, 7u, 142u, 108u, 0u, 0u, 7u, 5u, 421u, 0u, 0u,
7u, 87u, 162u, 0u, 0u, 7u, 180u, 302u, 0u, 0u, 7u, 154u, 93u, 0u, 0u, 7u,
201u, 195u, 0u, 0u, 7u, 107u, 164u, 0u, 0u, 7u, 119u, 329u, 0u, 0u, 7u, 147u,
324u, 0u, 0u, 7u, 203u, 313u, 0u, 0u, 2u, 405u, 375u, 49u, 0u, 2u, 419u, 171u,
316u, 0u, 2u, 349u, 101u, 90u, 0u, 2u, 348u, 82u, 332u, 0u, 2u, 343u, 0u, 334u,
0u, 2u, 340u, 276u, 52u, 0u, 2u, 302u, 277u, 51u, 0u, 2u, 313u, 184u, 48u, 0u,
2u, 324u, 173u, 47u, 0u, 2u, 329u, 176u, 3u, 0u, 2u, 421u, 109u, 376u, 0u, 2u,
415u, 256u, 372u, 0u, 2u, 408u, 75u, 371u, 0u, 2u, 401u, 183u, 208u, 0u, 2u, 394u,
102u, 177u, 0u, 2u, 385u, 68u, 217u, 0u, 2u, 378u, 1u, 250u, 0u, 2u, 11u, 94u,
249u, 0u, 2u, 164u, 136u, 248u, 0u, 2u, 195u, 247u, 244u, 0u, 2u, 93u, 241u, 236u,
0u, 2u, 162u, 234u, 230u, 0u, 2u, 108u, 218u, 409u, 0u, 2u, 69u, 222u, 144u, 0u,
2u, 76u, 202u, 99u, 0u, 2u, 19u, 165u, 122u, 0u, 2u, 88u, 160u, 142u, 0u, 2u,
91u, 206u, 416u, 0u, 2u, 156u, 5u, 87u, 0u, 2u, 163u, 38u, 422u, 0u, 2u, 110u,
425u, 303u, 0u, 2u, 178u, 180u, 154u, 0u, 2u, 92u, 301u, 300u, 0u, 2u, 285u, 292u,
34u, 0u, 2u, 289u, 197u, 201u, 0u, 2u, 263u, 149u, 107u, 0u, 2u, 245u, 395u, 257u,
0u, 2u, 190u, 258u, 393u, 0u, 2u, 112u, 330u, 392u, 0u, 2u, 185u, 270u, 387u, 0u,
2u, 194u, 386u, 260u, 0u, 2u, 146u, 384u, 262u, 0u, 2u, 187u, 382u, 264u, 0u, 2u,
228u, 265u, 380u, 0u, 2u, 252u, 326u, 379u, 0u, 2u, 293u, 318u, 17u, 0u, 2u, 286u,
253u, 267u, 0u, 2u, 172u, 290u, 268u, 0u, 2u, 182u, 25u, 271u, 0u, 2u, 83u, 272u,
62u, 0u, 2u, 13u, 314u, 186u, 0u, 2u, 89u, 307u, 179u, 0u, 2u, 175u, 299u, 274u,
0u, 2u, 7u, 311u, 275u, 0u, 5u, 107u, 419u, 196u, 0u, 5u, 201u, 349u, 226u, 0u,
5u, 154u, 348u, 226u, 0u, 5u, 180u, 343u, 196u, 0u, 5u, 87u, 340u, 226u, 0u, 5u,
5u, 302u, 226u, 0u, 5u, 142u, 313u, 196u, 0u, 5u, 160u, 324u, 226u, 0u, 5u, 122u,
329u, 226u, 0u, 5u, 99u, 421u, 196u, 0u, 5u, 144u, 415u, 226u, 0u, 5u, 230u, 408u,
226u, 0u, 5u, 234u, 401u, 209u, 0u, 5u, 236u, 158u, 45u, 0u, 5u, 241u, 394u, 209u,
0u, 5u, 244u, 291u, 45u, 0u, 5u, 247u, 385u, 209u, 0u, 5u, 248u, 231u, 45u, 0u,
5u, 249u, 378u, 196u, 0u, 5u, 250u, 11u, 226u, 0u, 5u, 217u, 164u, 226u, 0u, 5u,
68u, 195u, 209u, 0u, 5u, 177u, 198u, 45u, 0u, 5u, 102u, 93u, 209u, 0u, 5u, 208u,
152u, 45u, 0u, 5u, 183u, 162u, 209u, 0u, 5u, 371u, 108u, 45u, 0u, 5u, 372u, 69u,
196u, 0u, 5u, 376u, 76u, 226u, 0u, 5u, 3u, 19u, 226u, 0u, 5u, 47u, 88u, 209u,
0u, 5u, 48u, 91u, 45u, 0u, 5u, 51u, 156u, 209u, 0u, 5u, 52u, 163u, 45u, 0u,
5u, 334u, 178u, 209u, 0u, 5u, 332u, 92u, 45u, 0u, 5u, 90u, 46u, 45u, 0u, 5u,
316u, 289u, 226u, 0u, 5u, 311u, 263u, 226u, 0u, 4u, 263u, 370u, 405u, 0u, 4u, 289u,
367u, 405u, 0u, 4u, 92u, 366u, 405u, 0u, 4u, 178u, 365u, 405u, 0u, 4u, 163u, 321u,
405u, 0u, 4u, 156u, 201u, 40u, 0u, 4u, 91u, 154u, 39u, 0u, 4u, 88u, 87u, 40u,
0u, 4u, 19u, 5u, 39u, 0u, 4u, 76u, 160u, 40u, 0u, 4u, 69u, 122u, 39u, 0u,
4u, 108u, 144u, 40u, 0u, 4u, 162u, 230u, 39u, 0u, 4u, 93u, 99u, 364u, 0u, 4u,
195u, 50u, 234u, 0u, 4u, 164u, 95u, 234u, 0u, 4u, 11u, 107u, 338u, 0u, 4u, 378u,
144u, 114u, 0u, 4u, 385u, 201u, 280u, 0u, 4u, 394u, 230u, 114u, 0u, 4u, 401u, 154u,
280u, 0u, 4u, 408u, 234u, 40u, 0u, 4u, 415u, 99u, 353u, 0u, 4u, 421u, 50u, 241u,
0u, 4u, 329u, 95u, 241u, 0u, 4u, 324u, 180u, 338u, 0u, 4u, 313u, 144u, 81u, 0u,
4u, 302u, 87u, 280u, 0u, 4u, 340u, 230u, 81u, 0u, 4u, 343u, 5u, 280u, 0u, 4u,
348u, 241u, 40u, 0u, 4u, 349u, 99u, 346u, 0u, 4u, 419u, 50u, 247u, 0u, 4u, 46u,
95u, 247u, 0u, 4u, 152u, 142u, 338u, 0u, 4u, 198u, 144u, 278u, 0u, 4u, 231u, 160u,
280u, 0u, 4u, 291u, 230u, 278u, 0u, 4u, 158u, 122u, 280u, 0u, 4u, 299u, 247u, 40u,
0u, 4u, 179u, 250u, 40u, 0u, 4u, 186u, 217u, 39u, 0u, 4u, 62u, 249u, 364u, 0u,
4u, 25u, 50u, 68u, 0u, 4u, 290u, 95u, 68u, 0u, 4u, 253u, 107u, 215u, 0u, 4u,
17u, 250u, 114u, 0u, 4u, 379u, 201u, 284u, 0u, 4u, 380u, 217u, 114u, 0u, 4u, 382u,
154u, 284u, 0u, 4u, 384u, 68u, 40u, 0u, 4u, 386u, 249u, 353u, 0u, 4u, 387u, 50u,
102u, 0u, 4u, 392u, 95u, 102u, 0u, 4u, 393u, 180u, 215u, 0u, 4u, 395u, 250u, 81u,
0u, 4u, 34u, 87u, 284u, 0u, 4u, 300u, 217u, 81u, 0u, 4u, 301u, 5u, 284u, 0u,
4u, 303u, 102u, 40u, 0u, 4u, 38u, 249u, 346u, 0u, 4u, 206u, 50u, 183u, 0u, 4u,
218u, 95u, 183u, 0u, 4u, 375u, 142u, 215u, 0u, 4u, 102u, 156u, 40u, 0u, 4u, 68u,
91u, 40u, 0u, 4u, 247u, 88u, 40u, 0u, 4u, 241u, 19u, 40u, 0u, 4u, 234u, 76u,
40u, 0u, 4u, 230u, 69u, 40u, 0u, 4u, 144u, 108u, 40u, 0u, 4u, 215u, 162u, 40u,
0u, 4u, 338u, 195u, 216u, 0u, 4u, 280u, 408u, 40u, 0u, 4u, 307u, 421u, 216u, 0u,
4u, 314u, 348u, 40u, 0u, 4u, 318u, 419u, 216u, 0u, 4u, 326u, 299u, 40u, 0u, 4u,
270u, 179u, 40u, 0u, 4u, 330u, 186u, 40u, 0u, 4u, 425u, 25u, 216u, 0u, 4u, 422u,
384u, 40u, 0u, 4u, 416u, 387u, 216u, 0u, 4u, 409u, 303u, 40u, 0u, 4u, 149u, 206u,
216u, 0u, 4u, 197u, 250u, 278u, 0u, 4u, 165u, 160u, 284u, 0u, 4u, 202u, 217u, 278u,
0u, 4u, 136u, 122u, 284u, 0u, 4u, 94u, 183u, 40u, 0u, 4u, 75u, 376u, 40u, 0u,
4u, 256u, 3u, 39u, 0u, 4u, 277u, 372u, 364u, 0u, 4u, 276u, 50u, 47u, 0u, 4u,
82u, 95u, 47u, 0u, 4u, 101u, 107u, 74u, 0u, 4u, 275u, 376u, 114u, 0u, 4u, 274u,
201u, 288u, 0u, 4u, 272u, 3u, 114u, 0u, 4u, 271u, 154u, 288u, 0u, 4u, 268u, 47u,
40u, 0u, 4u, 267u, 372u, 353u, 0u, 4u, 265u, 50u, 51u, 0u, 4u, 264u, 95u, 51u,
0u, 4u, 262u, 180u, 74u, 0u, 4u, 260u, 376u, 81u, 0u, 4u, 258u, 87u, 288u, 0u,
4u, 257u, 3u, 81u, 0u, 4u, 222u, 5u, 288u, 0u, 4u, 49u, 51u, 40u, 0u, 4u,
1u, 372u, 346u, 0u, 4u, 109u, 50u, 334u, 0u, 4u, 176u, 95u, 334u, 0u, 4u, 173u,
142u, 74u, 0u, 4u, 184u, 376u, 278u, 0u, 4u, 0u, 168u, 110u, 0u, 4u, 171u, 160u,
288u, 0u, 4u, 118u, 3u, 278u, 0u, 4u, 143u, 181u, 110u, 0u, 4u, 53u, 122u, 288u,
0u, 4u, 254u, 334u, 40u, 0u, 4u, 238u, 192u, 285u, 0u, 4u, 166u, 316u, 40u, 0u,
4u, 188u, 311u, 39u, 0u, 4u, 191u, 397u, 10u, 0u, 4u, 204u, 321u, 397u, 0u, 4u,
157u, 107u, 285u, 0u, 4u, 246u, 390u, 209u, 0u, 4u, 206u, 102u, 8u, 0u, 4u, 303u,
68u, 9u, 0u, 4u, 387u, 247u, 8u, 0u, 4u, 384u, 241u, 9u, 0u, 4u, 25u, 234u,
8u, 0u, 4u, 186u, 230u, 9u, 0u, 4u, 179u, 144u, 8u, 0u, 4u, 299u, 215u, 9u,
0u, 4u, 419u, 321u, 144u, 0u, 4u, 348u, 321u, 215u, 0u, 4u, 421u, 102u, 15u, 0u,
4u, 408u, 68u, 16u, 0u, 4u, 195u, 280u, 10u, 0u, 4u, 162u, 247u, 15u, 0u, 4u,
108u, 241u, 16u, 0u, 4u, 69u, 314u, 10u, 0u, 4u, 76u, 234u, 15u, 0u, 4u, 19u,
230u, 16u, 0u, 4u, 88u, 326u, 10u, 0u, 4u, 91u, 270u, 8u, 0u, 4u, 156u, 330u,
9u, 0u, 4u, 334u, 321u, 270u, 0u, 4u, 51u, 321u, 330u, 0u, 4u, 47u, 102u, 22u,
0u, 4u, 3u, 68u, 23u, 0u, 4u, 376u, 422u, 10u, 0u, 4u, 183u, 247u, 22u, 0u,
4u, 217u, 241u, 23u, 0u, 4u, 250u, 409u, 10u, 0u, 4u, 390u, 234u, 22u, 0u, 4u,
397u, 230u, 23u, 0u, 4u, 74u, 94u, 40u, 0u, 4u, 288u, 75u, 40u, 0u, 4u, 284u,
256u, 40u, 0u, 4u, 243u, 276u, 216u, 0u, 4u, 242u, 102u, 29u, 0u, 4u, 240u, 68u,
30u, 0u, 4u, 239u, 268u, 40u, 0u, 4u, 237u, 265u, 216u, 0u, 4u, 235u, 247u, 29u,
0u, 4u, 233u, 241u, 30u, 0u, 4u, 232u, 49u, 40u, 0u, 4u, 221u, 109u, 216u, 0u,
4u, 220u, 234u, 29u, 0u, 4u, 219u, 230u, 30u, 0u, 4u, 145u, 254u, 40u, 0u, 4u,
161u, 166u, 40u, 0u, 4u, 100u, 188u, 40u, 0u, 4u, 104u, 316u, 114u, 0u, 4u, 140u,
201u, 292u, 0u, 4u, 139u, 102u, 36u, 0u, 4u, 138u, 311u, 114u, 0u, 4u, 137u, 154u,
292u, 0u, 4u, 135u, 68u, 37u, 0u, 4u, 134u, 180u, 285u, 0u, 4u, 133u, 377u, 209u,
0u, 4u, 132u, 316u, 81u, 0u, 4u, 130u, 87u, 292u, 0u, 4u, 129u, 247u, 36u, 0u,
4u, 128u, 311u, 81u, 0u, 4u, 127u, 5u, 292u, 0u, 4u, 126u, 241u, 37u, 0u, 4u,
125u, 142u, 285u, 0u, 4u, 124u, 170u, 209u, 0u, 2u, 188u, 107u, 61u, 0u, 2u, 166u,
206u, 303u, 0u, 2u, 254u, 180u, 123u, 0u, 2u, 109u, 387u, 384u, 0u, 2u, 49u, 142u,
223u, 0u, 2u, 265u, 25u, 186u, 0u, 2u, 268u, 99u, 227u, 0u, 2u, 276u, 179u, 279u,
0u, 2u, 256u, 299u, 281u, 0u, 2u, 75u, 338u, 199u, 0u, 2u, 94u, 113u, 164u, 0u,
2u, 409u, 378u, 350u, 0u, 2u, 422u, 394u, 354u, 0u, 2u, 330u, 195u, 225u, 0u, 2u,
270u, 307u, 117u, 0u, 2u, 326u, 189u, 329u, 0u, 2u, 314u, 313u, 356u, 0u, 2u, 280u,
340u, 359u, 0u, 2u, 215u, 69u, 33u, 0u, 2u, 144u, 318u, 297u, 0u, 2u, 241u, 298u,
46u, 0u, 2u, 247u, 198u, 424u, 0u, 2u, 68u, 291u, 420u, 0u, 2u, 102u, 88u, 224u,
0u, 2u, 5u, 249u, 44u, 0u, 2u, 87u, 91u, 327u, 0u, 2u, 154u, 156u, 328u, 0u,
2u, 201u, 425u, 309u, 0u, 2u, 285u, 310u, 290u, 0u, 2u, 170u, 17u, 417u, 0u, 2u,
377u, 380u, 414u, 0u, 2u, 81u, 376u, 331u, 0u, 2u, 114u, 416u, 319u, 0u, 2u, 30u,
320u, 392u, 0u, 2u, 29u, 395u, 413u, 0u, 2u, 23u, 300u, 412u, 0u, 2u, 22u, 250u,
333u, 0u, 2u, 16u, 149u, 336u, 0u, 2u, 15u, 337u, 218u, 0u, 2u, 121u, 197u, 407u,
0u, 2u, 120u, 202u, 406u, 0u, 2u, 141u, 372u, 287u, 0u, 2u, 214u, 243u, 351u, 0u,
2u, 212u, 352u, 82u, 0u, 2u, 210u, 275u, 403u, 0u, 2u, 35u, 272u, 402u, 0u, 2u,
103u, 237u, 362u, 0u, 2u, 207u, 363u, 264u, 0u, 2u, 211u, 260u, 400u, 0u, 2u, 213u,
257u, 399u, 0u, 2u, 42u, 221u, 373u, 0u, 2u, 105u, 374u, 176u, 0u, 2u, 116u, 184u,
0u, 0u, 2u, 98u, 118u, 143u, 0u, 2u, 28u, 191u, 255u, 0u, 2u, 426u, 245u, 190u,
0u, 2u, 427u, 246u, 259u, 0u, 2u, 428u, 194u, 261u, 0u, 2u, 429u, 187u, 228u, 0u,
2u, 430u, 133u, 266u, 0u, 2u, 431u, 172u, 269u, 0u, 2u, 432u, 182u, 83u, 0u, 2u,
433u, 124u, 273u, 0u, 7u, 124u, 11u, 0u, 0u, 7u, 133u, 409u, 0u, 0u, 7u, 221u,
385u, 0u, 0u, 7u, 237u, 422u, 0u, 0u, 7u, 243u, 401u, 0u, 0u, 7u, 250u, 324u,
0u, 0u, 7u, 376u, 314u, 0u, 0u, 7u, 156u, 302u, 0u, 0u, 7u, 91u, 280u, 0u,
0u, 7u, 88u, 343u, 0u, 0u, 7u, 69u, 152u, 0u, 0u, 7u, 195u, 247u, 0u, 0u,
7u, 299u, 231u, 0u, 0u, 7u, 179u, 68u, 0u, 0u, 7u, 186u, 158u, 0u, 0u, 7u,
25u, 253u, 0u, 0u, 7u, 384u, 170u, 0u, 0u, 7u, 387u, 379u, 0u, 0u, 7u, 303u,
377u, 0u, 0u, 7u, 206u, 382u, 0u, 0u, 7u, 246u, 393u, 0u, 0u, 7u, 191u, 29u,
0u, 0u, 7u, 143u, 34u, 0u, 0u, 7u, 118u, 23u, 0u, 0u, 7u, 0u, 301u, 0u,
0u, 7u, 184u, 375u, 0u, 0u, 7u, 176u, 121u, 0u, 0u, 7u, 257u, 165u, 0u, 0u,
7u, 260u, 120u, 0u, 0u, 7u, 264u, 136u, 0u, 0u, 7u, 272u, 101u, 0u, 0u, 7u,
275u, 210u, 0u, 0u, 7u, 82u, 274u, 0u, 0u, 7u, 202u, 35u, 0u, 0u, 7u, 197u,
271u, 0u, 0u, 7u, 149u, 262u, 0u, 0u, 7u, 416u, 211u, 0u, 0u, 7u, 425u, 258u,
0u, 0u, 7u, 318u, 213u, 0u, 0u, 7u, 307u, 222u, 0u, 0u, 7u, 338u, 173u, 0u,
0u, 7u, 218u, 116u, 0u, 0u, 7u, 300u, 171u, 0u, 0u, 7u, 395u, 98u, 0u, 0u,
7u, 392u, 53u, 0u, 0u, 7u, 380u, 238u, 0u, 0u, 7u, 17u, 112u, 0u, 0u, 7u,
290u, 157u, 0u, 0u, 7u, 291u, 140u, 0u, 0u, 7u, 198u, 137u, 0u, 0u, 7u, 46u,
252u, 0u, 0u, 7u, 340u, 134u, 0u, 0u, 7u, 313u, 130u, 0u, 0u, 7u, 329u, 127u,
0u, 0u, 7u, 394u, 13u, 0u, 0u, 7u, 378u, 125u, 0u, 0u, 4u, 98u, 188u, 405u,
0u, 4u, 116u, 254u, 405u, 0u, 4u, 213u, 49u, 405u, 0u, 4u, 211u, 268u, 32u, 0u,
4u, 35u, 268u, 282u, 0u, 4u, 210u, 268u, 283u, 0u, 4u, 120u, 268u, 150u, 0u, 4u,
121u, 268u, 405u, 0u, 4u, 23u, 5u, 32u, 0u, 4u, 29u, 5u, 282u, 0u, 4u, 377u,
5u, 283u, 0u, 4u, 170u, 5u, 150u, 0u, 4u, 68u, 5u, 405u, 0u, 4u, 247u, 74u,
10u, 0u, 4u, 280u, 141u, 32u, 0u, 4u, 314u, 288u, 8u, 0u, 4u, 422u, 284u, 9u,
0u, 4u, 409u, 141u, 282u, 0u, 4u, 125u, 321u, 288u, 0u, 4u, 127u, 141u, 283u, 0u,
4u, 130u, 321u, 284u, 0u, 4u, 134u, 141u, 150u, 0u, 4u, 137u, 239u, 10u, 0u, 4u,
140u, 141u, 405u, 0u, 4u, 157u, 232u, 10u, 0u, 4u, 238u, 145u, 10u, 0u, 4u, 53u,
161u, 8u, 0u, 4u, 171u, 100u, 9u, 0u, 4u, 173u, 321u, 161u, 0u, 4u, 222u, 321u,
100u, 0u, 4u, 258u, 168u, 426u, 0u, 4u, 262u, 181u, 426u, 0u, 4u, 271u, 168u, 429u,
0u, 4u, 274u, 181u, 429u, 0u, 4u, 101u, 316u, 278u, 0u, 4u, 136u, 168u, 432u, 0u,
4u, 165u, 160u, 292u, 0u, 4u, 375u, 234u, 36u, 0u, 4u, 301u, 311u, 278u, 0u, 4u,
34u, 181u, 432u, 0u, 4u, 393u, 122u, 292u, 0u, 4u, 382u, 230u, 37u, 0u, 2u, 100u,
166u, 155u, 0u, 2u, 161u, 109u, 169u, 0u, 2u, 145u, 265u, 193u, 0u, 2u, 232u, 12u,
276u, 0u, 2u, 239u, 35u, 419u, 0u, 2u, 284u, 210u, 348u, 0u, 2u, 288u, 120u, 361u,
0u, 2u, 74u, 75u, 94u, 0u, 2u, 230u, 229u, 133u, 0u, 2u, 234u, 24u, 237u, 0u,
2u, 311u, 270u, 326u, 0u, 2u, 316u, 381u, 376u, 0u, 2u, 122u, 159u, 91u, 0u, 2u,
160u, 144u, 241u, 0u, 2u, 181u, 388u, 195u, 0u, 2u, 168u, 389u, 179u, 0u, 2u, 292u,
18u, 87u, 0u, 2u, 278u, 29u, 334u, 0u, 2u, 37u, 377u, 51u, 0u, 2u, 36u, 170u,
360u, 0u, 2u, 379u, 201u, 285u, 0u, 2u, 253u, 148u, 384u, 0u, 2u, 158u, 398u, 303u,
0u, 2u, 231u, 114u, 30u, 0u, 2u, 152u, 404u, 191u, 0u, 2u, 343u, 131u, 118u, 0u,
2u, 302u, 16u, 15u, 0u, 2u, 324u, 410u, 176u, 0u, 2u, 401u, 411u, 260u, 0u, 2u,
385u, 247u, 342u, 0u, 2u, 11u, 314u, 344u, 0u, 2u, 13u, 422u, 345u, 0u, 2u, 252u,
409u, 125u, 0u, 2u, 112u, 127u, 130u, 0u, 2u, 164u, 134u, 355u, 0u, 2u, 372u, 214u,
212u, 0u, 2u, 249u, 205u, 275u, 0u, 2u, 99u, 418u, 202u, 0u, 2u, 142u, 137u, 357u,
0u, 2u, 180u, 103u, 207u, 0u, 2u, 107u, 423u, 416u, 0u, 2u, 83u, 31u, 318u, 0u,
2u, 182u, 157u, 368u, 0u, 2u, 172u, 42u, 105u, 0u, 2u, 228u, 296u, 218u, 0u, 2u,
187u, 295u, 395u, 0u, 2u, 194u, 238u, 341u, 0u, 2u, 190u, 251u, 380u, 0u, 2u, 245u,
53u, 339u, 0u, 2u, 399u, 171u, 335u, 0u, 2u, 400u, 17u, 119u, 0u, 2u, 402u, 427u,
185u, 0u, 2u, 403u, 104u, 258u, 0u, 2u, 406u, 138u, 262u, 0u, 2u, 407u, 46u, 147u,
0u, 2u, 412u, 430u, 293u, 0u, 2u, 413u, 132u, 271u, 0u, 2u, 414u, 128u, 274u, 0u,
2u, 417u, 394u, 203u, 0u, 2u, 420u, 433u, 89u, 0u, 2u, 424u, 101u, 136u, 0u, 2u,
359u, 301u, 34u, 0u, 2u, 34u, 232u, 256u, 0u, 2u, 301u, 304u, 74u, 0u, 2u, 136u,
230u, 221u, 0u, 2u, 101u, 234u, 243u, 0u, 2u, 274u, 294u, 311u, 0u, 2u, 271u, 316u,
156u, 0u, 2u, 262u, 122u, 88u, 0u, 2u, 258u, 115u, 160u, 0u, 2u, 171u, 181u, 299u,
0u, 2u, 53u, 168u, 186u, 0u, 2u, 238u, 292u, 154u, 0u, 2u, 157u, 153u, 379u, 0u,
2u, 137u, 253u, 387u, 0u, 2u, 134u, 158u, 206u, 0u, 2u, 130u, 151u, 231u, 0u, 2u,
127u, 152u, 143u, 0u, 2u, 125u, 343u, 0u, 0u, 2u, 409u, 60u, 302u, 0u, 2u, 422u,
324u, 257u, 0u, 2u, 314u, 401u, 264u, 0u, 2u, 247u, 20u, 11u, 0u, 2u, 170u, 106u,
372u, 0u, 2u, 377u, 249u, 82u, 0u, 2u, 29u, 99u, 197u, 0u, 2u, 120u, 14u, 180u,
0u, 2u, 210u, 107u, 425u, 0u, 2u, 35u, 83u, 307u, 0u, 2u, 394u, 110u, 172u, 0u,
2u, 46u, 228u, 300u, 0u, 2u, 17u, 187u, 392u, 0u, 2u, 380u, 27u, 245u, 0u, 4u,
245u, 321u, 100u, 0u, 4u, 187u, 163u, 100u, 0u, 4u, 228u, 321u, 161u, 0u, 4u, 172u,
163u, 161u, 0u, 4u, 83u, 321u, 145u, 0u, 4u, 107u, 163u, 145u, 0u, 4u, 180u, 211u,
405u, 0u, 4u, 99u, 239u, 405u, 0u, 4u, 249u, 284u, 405u, 0u, 4u, 372u, 288u, 405u,
0u, 4u, 11u, 192u, 301u, 0u, 4u, 401u, 268u, 100u, 0u, 4u, 324u, 192u, 274u, 0u,
4u, 302u, 268u, 161u, 0u, 4u, 343u, 192u, 258u, 0u, 4u, 152u, 268u, 145u, 0u, 4u,
231u, 23u, 405u, 0u, 4u, 158u, 278u, 405u, 0u, 4u, 253u, 37u, 405u, 0u, 4u, 379u,
36u, 405u, 0u, 4u, 292u, 192u, 157u, 0u, 4u, 168u, 5u, 100u, 0u, 4u, 181u, 192u,
130u, 0u, 4u, 160u, 5u, 161u, 0u, 4u, 122u, 192u, 409u, 0u, 4u, 316u, 5u, 145u,
0u, 4u, 311u, 280u, 405u, 0u, 4u, 234u, 252u, 405u, 0u, 4u, 230u, 112u, 405u, 0u,
4u, 74u, 164u, 405u, 0u, 4u, 232u, 192u, 170u, 0u, 4u, 392u, 141u, 100u, 0u, 4u,
300u, 192u, 120u, 0u, 4u, 307u, 141u, 161u, 0u, 4u, 425u, 192u, 394u, 0u, 4u, 197u,
141u, 145u, 0u, 5u, 394u, 136u, 226u, 0u, 5u, 120u, 101u, 226u, 0u, 5u, 170u, 271u,
226u, 0u, 5u, 409u, 262u, 226u, 0u, 5u, 130u, 171u, 226u, 0u, 5u, 157u, 53u, 226u,
0u, 5u, 258u, 137u, 226u, 0u, 5u, 274u, 134u, 226u, 0u, 5u, 301u, 127u, 226u, 0u,
5u, 164u, 125u, 226u, 0u, 5u, 112u, 422u, 226u, 0u, 5u, 252u, 314u, 226u, 0u, 5u,
36u, 377u, 226u, 0u, 5u, 37u, 29u, 226u, 0u, 5u, 278u, 210u, 226u, 0u, 5u, 288u,
35u, 226u, 0u, 5u, 284u, 46u, 226u, 0u, 5u, 239u, 17u, 226u, 0u, 5u, 280u, 190u,
196u, 0u, 5u, 23u, 400u, 209u, 0u, 5u, 211u, 402u, 45u, 0u, 5u, 141u, 407u, 209u,
0u, 5u, 5u, 412u, 45u, 0u, 5u, 268u, 417u, 209u, 0u, 5u, 82u, 420u, 45u, 0u,
4u, 17u, 394u, 40u, 0u, 4u, 46u, 120u, 39u, 0u, 4u, 35u, 170u, 40u, 0u, 4u,
210u, 409u, 39u, 0u, 4u, 29u, 130u, 40u, 0u, 4u, 377u, 157u, 39u, 0u, 4u, 314u,
258u, 40u, 0u, 4u, 422u, 274u, 39u, 0u, 4u, 125u, 301u, 40u, 0u, 4u, 127u, 164u,
39u, 0u, 4u, 134u, 112u, 40u, 0u, 4u, 137u, 252u, 39u, 0u, 4u, 53u, 36u, 40u,
0u, 4u, 171u, 37u, 39u, 0u, 4u, 262u, 278u, 40u, 0u, 4u, 271u, 288u, 39u, 0u,
4u, 101u, 284u, 40u, 0u, 4u, 136u, 239u, 39u, 0u, 4u, 420u, 280u, 364u, 0u, 4u,
417u, 50u, 23u, 0u, 4u, 412u, 95u, 23u, 0u, 4u, 407u, 23u, 40u, 0u, 4u, 402u,
280u, 353u, 0u, 4u, 400u, 50u, 141u, 0u, 4u, 190u, 95u, 141u, 0u, 4u, 209u, 141u,
40u, 0u, 4u, 45u, 280u, 346u, 0u, 4u, 264u, 50u, 268u, 0u, 4u, 257u, 95u, 268u,
0u, 4u, 0u, 268u, 40u, 0u, 4u, 268u, 17u, 40u, 0u, 4u, 141u, 46u, 40u, 0u,
4u, 23u, 35u, 40u, 0u, 4u, 239u, 210u, 40u, 0u, 4u, 284u, 29u, 40u, 0u, 4u,
288u, 377u, 40u, 0u, 4u, 278u, 314u, 40u, 0u, 4u, 37u, 422u, 40u, 0u, 4u, 36u,
125u, 40u, 0u, 4u, 252u, 127u, 40u, 0u, 4u, 112u, 134u, 40u, 0u, 4u, 164u, 137u,
40u, 0u, 4u, 301u, 53u, 40u, 0u, 4u, 274u, 171u, 40u, 0u, 4u, 258u, 262u, 40u,
0u, 4u, 157u, 271u, 40u, 0u, 4u, 130u, 101u, 40u, 0u, 4u, 409u, 136u, 40u, 0u,
4u, 170u, 417u, 216u, 0u, 4u, 120u, 407u, 40u, 0u, 4u, 394u, 400u, 216u, 0u, 4u,
346u, 209u, 40u, 0u, 4u, 353u, 264u, 216u, 0u, 4u, 364u, 0u, 40u, 0u, 2u, 0u,
98u, 245u, 0u, 2u, 264u, 116u, 228u, 0u, 2u, 209u, 213u, 83u, 0u, 2u, 400u, 34u,
58u, 0u, 2u, 407u, 93u, 11u, 0u, 2u, 417u, 415u, 324u, 0u, 2u, 136u, 349u, 343u,
0u, 2u, 101u, 238u, 391u, 0u, 2u, 271u, 62u, 292u, 0u, 2u, 262u, 386u, 181u, 0u,
2u, 171u, 38u, 122u, 0u, 2u, 53u, 247u, 13u, 0u, 2u, 137u, 277u, 232u, 0u, 2u,
134u, 267u, 300u, 0u, 2u, 127u, 1u, 425u, 0u, 2u, 125u, 280u, 90u, 0u, 2u, 422u,
380u, 399u, 0u, 2u, 314u, 170u, 97u, 0u, 2u, 377u, 174u, 412u, 0u, 2u, 29u, 394u,
323u, 0u, 2u, 210u, 322u, 190u, 0u, 2u, 35u, 353u, 312u, 0u, 2u, 46u, 308u, 257u,
0u, 4u, 353u, 0u, 405u, 0u, 4u, 394u, 264u, 405u, 0u, 4u, 170u, 209u, 405u, 0u,
4u, 257u, 370u, 400u, 0u, 4u, 190u, 367u, 400u, 0u, 4u, 412u, 366u, 400u, 0u, 4u,
280u, 365u, 400u, 0u, 4u, 425u, 188u, 400u, 0u, 4u, 300u, 268u, 8u, 0u, 4u, 232u,
141u, 9u, 0u, 4u, 122u, 0u, 400u, 0u, 4u, 181u, 321u, 400u, 0u, 4u, 292u, 254u,
400u, 0u, 4u, 343u, 23u, 8u, 0u, 4u, 324u, 239u, 9u, 0u, 4u, 11u, 264u, 400u,
0u, 4u, 83u, 49u, 400u, 0u, 4u, 228u, 284u, 8u, 0u, 4u, 245u, 288u, 9u, 0u,
4u, 380u, 209u, 400u, 0u, 4u, 247u, 370u, 101u, 0u, 4u, 238u, 367u, 101u, 0u, 4u,
34u, 366u, 101u, 0u, 4u, 399u, 365u, 101u, 0u, 4u, 13u, 188u, 101u, 0u, 4u, 213u,
278u, 8u, 0u, 4u, 116u, 37u, 9u, 0u, 4u, 98u, 0u, 101u, 0u, 4u, 1u, 321u,
101u, 0u, 4u, 267u, 254u, 101u, 0u, 4u, 277u, 36u, 8u, 0u, 4u, 38u, 252u, 9u,
0u, 4u, 386u, 264u, 101u, 0u, 4u, 62u, 49u, 101u, 0u, 4u, 349u, 112u, 8u, 0u,
4u, 415u, 164u, 9u, 0u, 4u, 93u, 209u, 101u, 0u, 4u, 90u, 301u, 8u, 0u, 4u,
308u, 274u, 9u, 0u, 4u, 312u, 258u, 8u, 0u, 4u, 322u, 157u, 9u, 0u, 4u, 323u,
130u, 8u, 0u, 4u, 174u, 409u, 9u, 0u, 4u, 97u, 125u, 32u, 0u, 4u, 391u, 125u,
282u, 0u, 4u, 58u, 125u, 283u, 0u, 4u, 17u, 125u, 150u, 0u, 4u, 216u, 125u, 100u,
0u, 4u, 95u, 120u, 10u, 0u, 4u, 50u, 125u, 405u, 0u, 4u, 143u, 125u, 161u, 0u,
4u, 206u, 346u, 10u, 0u, 4u, 387u, 125u, 145u, 0u, 4u, 186u, 364u, 10u, 0u, 2u,
125u, 353u, 187u, 0u, 2u, 101u, 394u, 172u, 0u, 2u, 400u, 170u, 107u, 0u, 2u, 364u,
180u, 257u, 0u, 2u, 346u, 99u, 190u, 0u, 2u, 120u, 249u, 412u, 0u, 2u, 409u, 372u,
280u, 0u, 2u, 130u, 300u, 421u, 0u, 2u, 157u, 232u, 408u, 0u, 2u, 258u, 121u, 181u,
0u, 2u, 274u, 343u, 162u, 0u, 2u, 301u, 324u, 108u, 0u, 2u, 164u, 228u, 76u, 0u,
2u, 112u, 245u, 19u, 0u, 2u, 252u, 231u, 247u, 0u, 2u, 36u, 158u, 238u, 0u, 2u,
37u, 253u, 34u, 0u, 2u, 278u, 379u, 399u, 0u, 2u, 288u, 213u, 47u, 0u, 2u, 284u,
116u, 3u, 0u, 2u, 239u, 68u, 1u, 0u, 2u, 23u, 277u, 183u, 0u, 2u, 141u, 38u,
217u, 0u, 2u, 268u, 349u, 390u, 0u, 2u, 283u, 415u, 397u, 0u, 2u, 282u, 53u, 200u,
0u, 2u, 150u, 90u, 242u, 0u, 2u, 32u, 308u, 240u, 0u, 2u, 10u, 312u, 235u, 0u,
2u, 299u, 322u, 233u, 0u, 2u, 88u, 323u, 220u, 0u, 2u, 156u, 174u, 219u, 0u, 2u,
243u, 422u, 28u, 0u, 2u, 221u, 391u, 173u, 0u, 2u, 154u, 58u, 222u, 0u, 2u, 256u,
17u, 204u, 0u, 2u, 110u, 314u, 377u, 0u, 2u, 14u, 95u, 67u, 0u, 2u, 106u, 29u,
210u, 0u, 2u, 60u, 206u, 317u, 0u, 2u, 151u, 35u, 46u, 0u, 2u, 153u, 186u, 306u,
0u, 4u, 186u, 258u, 100u, 0u, 4u, 206u, 258u, 161u, 0u, 4u, 95u, 258u, 145u, 0u,
4u, 17u, 239u, 100u, 0u, 4u, 58u, 239u, 161u, 0u, 4u, 391u, 239u, 145u, 0u, 4u,
174u, 370u, 282u, 0u, 4u, 323u, 367u, 282u, 0u, 4u, 322u, 366u, 282u, 0u, 4u, 312u,
365u, 282u, 0u, 4u, 308u, 188u, 282u, 0u, 4u, 90u, 0u, 282u, 0u, 4u, 415u, 321u,
282u, 0u, 4u, 349u, 254u, 282u, 0u, 4u, 38u, 264u, 282u, 0u, 4u, 277u, 49u, 282u,
0u, 4u, 1u, 209u, 282u, 0u, 4u, 116u, 97u, 405u, 0u, 4u, 213u, 370u, 243u, 0u,
4u, 399u, 221u, 405u, 0u, 4u, 34u, 367u, 243u, 0u, 4u, 238u, 154u, 405u, 0u, 4u,
247u, 366u, 243u, 0u, 4u, 245u, 256u, 405u, 0u, 4u, 228u, 365u, 243u, 0u, 4u, 324u,
188u, 243u, 0u, 4u, 343u, 0u, 243u, 0u, 4u, 181u, 321u, 243u, 0u, 4u, 232u, 254u,
243u, 0u, 4u, 300u, 264u, 243u, 0u, 4u, 280u, 49u, 243u, 0u, 4u, 412u, 209u, 243u,
0u, 7u, 256u, 364u, 0u, 0u, 7u, 154u, 346u, 0u, 0u, 7u, 221u, 120u, 0u, 0u,
7u, 243u, 409u, 0u, 0u, 7u, 282u, 407u, 0u, 0u, 7u, 239u, 417u, 0u, 0u, 7u,
258u, 136u, 0u, 0u, 7u, 97u, 252u, 0u, 0u, 7u, 209u, 36u, 0u, 0u, 7u, 264u,
37u, 0u, 0u, 7u, 0u, 278u, 0u, 0u, 7u, 49u, 271u, 0u, 0u, 7u, 254u, 262u,
0u, 0u, 7u, 188u, 171u, 0u, 0u, 7u, 365u, 137u, 0u, 0u, 7u, 366u, 134u, 0u,
0u, 7u, 367u, 127u, 0u, 0u, 7u, 370u, 403u, 0u, 0u, 7u, 190u, 406u, 0u, 0u,
7u, 257u, 413u, 0u, 0u, 7u, 170u, 414u, 0u, 0u, 7u, 394u, 424u, 0u, 0u, 7u,
353u, 165u, 0u, 0u, 7u, 46u, 359u, 0u, 0u, 7u, 35u, 393u, 0u, 0u, 2u, 278u,
59u, 256u, 0u, 2u, 37u, 63u, 154u, 0u, 2u, 36u, 64u, 221u, 0u, 2u, 252u, 65u,
243u, 0u, 2u, 409u, 111u, 282u, 0u, 2u, 120u, 130u, 157u, 0u, 2u, 346u, 167u, 239u,
0u, 2u, 364u, 274u, 301u, 0u, 2u, 127u, 26u, 258u, 0u, 2u, 134u, 164u, 112u, 0u,
2u, 137u, 66u, 97u, 0u, 2u, 171u, 70u, 209u, 0u, 2u, 262u, 71u, 264u, 0u, 2u,
271u, 72u, 0u, 0u, 2u, 136u, 305u, 49u, 0u, 2u, 417u, 288u, 284u, 0u, 2u, 407u,
315u, 254u, 0u, 2u, 359u, 23u, 141u, 0u, 2u, 424u, 325u, 188u, 0u, 2u, 414u, 268u,
283u, 0u, 2u, 413u, 311u, 174u, 0u, 2u, 406u, 234u, 323u, 0u, 2u, 403u, 230u, 322u,
0u, 2u, 393u, 74u, 312u, 0u, 2u, 165u, 347u, 365u, 0u, 2u, 210u, 150u, 32u, 0u,
2u, 29u, 140u, 415u, 0u, 2u, 377u, 358u, 366u, 0u, 2u, 314u, 10u, 299u, 0u, 2u,
422u, 369u, 367u, 0u, 2u, 53u, 88u, 156u, 0u, 2u, 379u, 116u, 213u, 0u, 2u, 253u,
399u, 34u, 0u, 2u, 158u, 238u, 247u, 0u, 2u, 231u, 245u, 228u, 0u, 2u, 372u, 426u,
110u, 0u, 2u, 249u, 428u, 370u, 0u, 2u, 99u, 146u, 190u, 0u, 2u, 180u, 50u, 181u,
0u, 2u, 107u, 429u, 106u, 0u, 2u, 172u, 286u, 257u, 0u, 2u, 187u, 431u, 170u, 0u,
2u, 222u, 432u, 151u, 0u, 2u, 173u, 175u, 394u, 0u, 2u, 68u, 7u, 46u, 0u, 2u,
46u, 409u, 124u, 0u, 2u, 394u, 120u, 330u, 0u, 2u, 170u, 346u, 250u, 0u, 2u, 257u,
364u, 215u, 0u, 2u, 190u, 127u, 69u, 0u, 2u, 370u, 134u, 102u, 0u, 2u, 367u, 136u,
25u, 0u, 2u, 366u, 417u, 81u, 0u, 2u, 365u, 407u, 246u, 0u, 2u, 188u, 359u, 22u,
0u, 2u, 254u, 424u, 184u, 0u, 2u, 49u, 414u, 385u, 0u, 2u, 0u, 165u, 272u, 0u,
2u, 264u, 210u, 142u, 0u, 2u, 209u, 377u, 149u, 0u, 2u, 97u, 314u, 182u, 0u, 2u,
258u, 422u, 338u, 0u, 2u, 239u, 53u, 194u, 0u, 2u, 282u, 249u, 291u, 0u, 2u, 243u,
99u, 198u, 0u, 2u, 221u, 172u, 313u, 0u, 2u, 154u, 187u, 329u, 0u, 2u, 256u, 173u,
353u, 0u, 2u, 181u, 68u, 35u, 0u, 4u, 68u, 321u, 394u, 0u, 4u, 173u, 163u, 394u,
0u, 4u, 187u, 321u, 257u, 0u, 4u, 172u, 163u, 257u, 0u, 4u, 99u, 321u, 370u, 0u,
4u, 249u, 163u, 370u, 0u, 4u, 53u, 321u, 366u, 0u, 4u, 422u, 163u, 366u, 0u, 4u,
314u, 321u, 188u, 0u, 4u, 377u, 163u, 188u, 0u, 4u, 210u, 321u, 49u, 0u, 4u, 165u,
163u, 49u, 0u, 4u, 414u, 321u, 264u, 0u, 4u, 424u, 29u, 100u, 0u, 4u, 359u, 163u,
264u, 0u, 4u, 407u, 321u, 97u, 0u, 4u, 417u, 29u, 161u, 0u, 4u, 136u, 163u, 97u,
0u, 4u, 134u, 321u, 239u, 0u, 4u, 127u, 29u, 145u, 0u, 4u, 364u, 163u, 239u, 0u,
4u, 346u, 192u, 372u, 0u, 4u, 120u, 180u, 100u, 0u, 4u, 409u, 192u, 107u, 0u, 4u,
35u, 180u, 161u, 0u, 4u, 353u, 192u, 222u, 0u, 4u, 194u, 180u, 145u, 0u, 2u, 239u,
401u, 68u, 0u, 2u, 97u, 186u, 173u, 0u, 2u, 264u, 302u, 187u, 0u, 2u, 49u, 206u,
172u, 0u, 2u, 188u, 152u, 99u, 0u, 2u, 366u, 95u, 249u, 0u, 2u, 370u, 168u, 53u,
0u, 2u, 257u, 17u, 422u, 0u, 2u, 394u, 160u, 314u, 0u, 2u, 222u, 58u, 377u, 0u,
2u, 107u, 316u, 210u, 0u, 2u, 180u, 391u, 165u, 0u, 2u, 372u, 392u, 414u, 0u, 2u,
29u, 424u, 359u, 0u, 2u, 145u, 307u, 407u, 0u, 2u, 161u, 417u, 136u, 0u, 2u, 100u,
197u, 134u, 0u, 2u, 192u, 127u, 364u, 0u, 2u, 182u, 420u, 346u, 0u, 2u, 142u, 402u,
409u, 0u, 2u, 385u, 45u, 353u, 0u, 5u, 353u, 46u, 196u, 0u, 5u, 409u, 170u, 196u,
0u, 5u, 346u, 190u, 196u, 0u, 5u, 364u, 367u, 196u, 0u, 5u, 127u, 365u, 196u, 0u,
5u, 134u, 254u, 196u, 0u, 5u, 136u, 0u, 196u, 0u, 5u, 417u, 209u, 196u, 0u, 5u,
407u, 258u, 196u, 0u, 5u, 359u, 282u, 226u, 0u, 5u, 424u, 243u, 226u, 0u, 5u, 414u,
221u, 226u, 0u, 5u, 165u, 154u, 226u, 0u, 5u, 210u, 256u, 226u, 0u, 5u, 377u, 181u,
226u, 0u, 7u, 181u, 413u, 0u, 0u, 7u, 256u, 406u, 0u, 0u, 7u, 154u, 403u, 0u,
0u, 7u, 221u, 393u, 0u, 0u, 7u, 243u, 379u, 0u, 0u, 7u, 282u, 253u, 0u, 0u,
7u, 258u, 158u, 0u, 0u, 7u, 209u, 231u, 0u, 0u, 7u, 0u, 182u, 0u, 0u, 7u,
254u, 142u, 0u, 0u, 7u, 365u, 385u, 0u, 0u, 2u, 385u, 353u, 236u, 0u, 2u, 142u,
409u, 244u, 0u, 2u, 182u, 346u, 248u, 0u, 2u, 231u, 364u, 177u, 0u, 2u, 158u, 127u,
208u, 0u, 2u, 253u, 134u, 371u, 0u, 2u, 379u, 73u, 181u, 0u, 2u, 393u, 77u, 256u,
0u, 2u, 403u, 78u, 154u, 0u, 2u, 406u, 79u, 221u, 0u, 2u, 413u, 136u, 48u, 0u,
2u, 367u, 417u, 52u, 0u, 2u, 190u, 407u, 332u, 0u, 2u, 170u, 80u, 243u, 0u, 2u,
46u, 84u, 282u, 0u, 2u, 226u, 85u, 258u, 0u, 2u, 314u, 86u, 209u, 0u, 2u, 422u,
396u, 0u, 0u, 2u, 53u, 383u, 254u, 0u, 2u, 249u, 96u, 365u, 0u, 4u, 365u, 385u,
405u, 0u, 4u, 254u, 142u, 405u, 0u, 4u, 0u, 182u, 405u, 0u, 4u, 209u, 231u, 405u,
0u, 4u, 258u, 158u, 405u, 0u, 4u, 282u, 253u, 405u, 0u, 4u, 243u, 413u, 405u, 0u,
4u, 221u, 367u, 405u, 0u, 4u, 154u, 190u, 405u, 0u, 4u, 256u, 359u, 40u, 0u, 4u,
181u, 424u, 39u, 0u, 4u, 407u, 414u, 40u, 0u, 4u, 417u, 165u, 39u, 0u, 4u, 136u,
210u, 40u, 0u, 4u, 134u, 377u, 39u, 0u, 2u, 190u, 365u, 425u, 0u, 2u, 367u, 254u,
292u, 0u, 2u, 413u, 0u, 83u, 0u, 2u, 253u, 209u, 13u, 0u, 2u, 158u, 258u, 267u,
0u, 2u, 231u, 282u, 62u, 0u, 2u, 182u, 243u, 308u, 0u, 2u, 142u, 221u, 349u, 0u,
2u, 385u, 154u, 277u, 0u, 2u, 377u, 422u, 290u, 0u, 2u, 210u, 53u, 340u, 0u, 2u,
165u, 249u, 378u, 0u, 2u, 154u, 190u, 239u, 0u, 2u, 221u, 367u, 264u, 0u, 2u, 243u,
413u, 188u, 0u, 2u, 282u, 253u, 370u, 0u, 2u, 258u, 158u, 394u, 0u, 2u, 209u, 231u,
107u, 0u, 2u, 0u, 182u, 372u, 0u, 2u, 254u, 142u, 145u, 0u, 2u, 365u, 385u, 100u,
0u, 4u, 385u, 154u, 405u, 0u, 4u, 142u, 221u, 405u, 0u, 4u, 182u, 243u, 405u, 0u,
4u, 231u, 282u, 405u, 0u, 4u, 158u, 258u, 405u, 0u, 4u, 253u, 209u, 405u, 0u, 4u,
413u, 0u, 405u, 0u, 4u, 367u, 254u, 405u, 0u, 4u, 190u, 365u, 405u, 0u, 4u, 100u,
256u, 40u, 0u, 4u, 145u, 181u, 40u, 0u, 4u, 372u, 407u, 40u, 0u, 4u, 107u, 417u,
40u, 0u, 4u, 394u, 136u, 40u, 0u, 4u, 370u, 134u, 40u, 0u, 2u, 365u, 385u, 122u,
0u, 2u, 254u, 142u, 11u, 0u, 2u, 0u, 182u, 380u, 0u, 2u, 209u, 231u, 98u, 0u,
2u, 258u, 158u, 386u, 0u, 2u, 282u, 253u, 93u, 0u, 2u, 243u, 413u, 90u, 0u, 2u,
221u, 367u, 38u, 0u, 2u, 154u, 190u, 1u, 0u, 2u, 190u, 365u, 97u, 0u, 2u, 367u,
254u, 49u, 0u, 2u, 413u, 0u, 366u, 0u, 2u, 253u, 209u, 257u, 0u, 2u, 158u, 258u,
222u, 0u, 2u, 231u, 282u, 180u, 0u, 2u, 182u, 243u, 29u, 0u, 2u, 142u, 221u, 161u,
0u, 2u, 385u, 154u, 192u, 0u, 4u, 154u, 100u, 8u, 0u, 4u, 221u, 145u, 9u, 0u,
4u, 243u, 372u, 8u, 0u, 4u, 282u, 107u, 9u, 0u, 4u, 258u, 394u, 8u, 0u, 4u,
209u, 370u, 9u, 0u, 2u, 370u, 154u, 139u, 0u, 2u, 394u, 221u, 135u, 0u, 2u, 107u,
243u, 129u, 0u, 2u, 372u, 282u, 126u, 0u, 2u, 145u, 258u, 375u, 0u, 2u, 100u, 209u,
382u, 0u, 3u, 209u, 54u, 263u, 0u, 3u, 258u, 55u, 289u, 0u, 3u, 282u, 56u, 92u,
0u, 3u, 243u, 57u, 178u, 0u, 2u, 178u, 370u, 394u, 0u, 2u, 92u, 107u, 372u, 0u,
2u, 289u, 145u, 100u, 0u, 2u, 100u, 178u, 14u, 0u, 2u, 145u, 92u, 60u, 0u, 2u,
372u, 289u, 153u, 0u, 4u, 289u, 321u, 100u, 0u, 4u, 92u, 163u, 100u, 0u, 4u, 178u,
321u, 145u, 0u, 4u, 153u, 163u, 145u, 0u, 4u, 60u, 321u, 372u, 0u, 4u, 14u, 163u,
372u, 0u, 2u, 372u, 216u, 289u, 0u, 2u, 145u, 120u, 92u, 0u, 2u, 100u, 143u, 178u,
0u, 2u, 163u, 35u, 153u, 0u, 2u, 321u, 387u, 60u, 0u, 2u, 107u, 194u, 14u, 0u,
5u, 14u, 377u, 196u, 0u, 5u, 60u, 210u, 196u, 0u, 5u, 153u, 165u, 196u, 0u, 2u,
165u, 14u, 211u, 0u, 2u, 210u, 60u, 5u, 0u, 2u, 377u, 153u, 82u, 0u, 4u, 153u,
165u, 405u, 0u, 4u, 60u, 210u, 405u, 0u, 4u, 14u, 377u, 405u, 0u, 2u, 377u, 153u,
324u, 0u, 2u, 210u, 60u, 232u, 0u, 2u, 165u, 14u, 280u, 0u, 2u, 14u, 377u, 372u,
0u, 2u, 60u, 210u, 100u, 0u, 2u, 153u, 165u, 321u, 0u, 4u, 165u, 14u, 405u, 0u,
4u, 210u, 60u, 405u, 0u, 4u, 377u, 153u, 405u, 0u, 2u, 153u, 165u, 343u, 0u, 2u,
60u, 210u, 300u, 0u, 2u, 14u, 377u, 412u, 0u, 2u, 377u, 153u, 145u, 0u, 2u, 210u,
60u, 163u, 0u, 2u, 165u, 14u, 107u, 0u, 9u, 14u, 2u, 6u, 0u, 9u, 6u, 21u,
14u, 0u, 11u, 14u, 54u, 209u, 6u, 11u, 21u, 55u, 258u, 6u, 11u, 2u, 56u, 282u,
6u, 11u, 60u, 57u, 243u, 6u, 11u, 153u, 41u, 125u, 6u, 11u, 107u, 41u, 101u, 6u,
11u, 163u, 41u, 400u, 6u, 11u, 145u, 59u, 278u, 6u, 11u, 412u, 63u, 37u, 6u, 11u,
300u, 64u, 36u, 6u, 11u, 343u, 65u, 252u, 6u, 11u, 405u, 41u, 190u, 6u, 11u, 321u,
41u, 367u, 6u, 11u, 100u, 41u, 413u, 6u, 11u, 372u, 66u, 137u, 6u, 11u, 280u, 70u,
171u, 6u, 11u, 232u, 71u, 262u, 6u, 11u, 324u, 72u, 271u, 6u, 11u, 82u, 41u, 253u,
6u, 11u, 5u, 41u, 158u, 6u, 11u, 211u, 41u, 231u, 6u, 11u, 196u, 73u, 379u, 6u,
11u, 178u, 77u, 393u, 6u, 11u, 92u, 78u, 403u, 6u, 11u, 289u, 79u, 406u, 6u, 11u,
194u, 41u, 182u, 6u, 11u, 35u, 41u, 142u, 6u, 11u, 120u, 41u, 385u, 6u, 11u, 387u,
80u, 170u, 6u, 11u, 143u, 84u, 46u, 6u, 11u, 216u, 85u, 226u, 6u, 11u, 394u, 86u,
314u, 6u, 11u, 370u, 41u, 377u, 6u, 11u, 263u, 41u, 210u, 6u, 11u, 221u, 41u, 165u,
6u, 7u, 6u, 153u, 0u, 0u, 7u, 165u, 107u, 0u, 0u, 7u, 210u, 163u, 0u, 0u,
7u, 377u, 405u, 0u, 0u, 7u, 243u, 321u, 0u, 0u, 7u, 282u, 100u, 0u, 0u, 7u,
258u, 82u, 0u, 0u, 7u, 209u, 5u, 0u, 0u, 7u, 385u, 211u, 0u, 0u, 7u, 142u,
194u, 0u, 0u, 7u, 182u, 35u, 0u, 0u, 7u, 231u, 120u, 0u, 0u, 7u, 158u, 370u,
0u, 0u, 7u, 253u, 263u, 0u, 0u, 7u, 413u, 221u, 0u, 0u, 5u, 221u, 6u, 43u,
0u, 5u, 263u, 165u, 43u, 0u, 5u, 370u, 210u, 43u, 0u, 5u, 120u, 377u, 43u, 0u,
5u, 35u, 243u, 43u, 0u, 5u, 194u, 282u, 43u, 0u, 5u, 211u, 258u, 43u, 0u, 5u,
5u, 209u, 43u, 0u, 5u, 82u, 385u, 43u, 0u, 5u, 100u, 142u, 43u, 0u, 5u, 321u,
182u, 43u, 0u, 5u, 405u, 231u, 43u, 0u, 5u, 163u, 158u, 43u, 0u, 5u, 107u, 253u,
43u, 0u, 5u, 153u, 413u, 43u, 0u, 7u, 413u, 221u, 0u, 0u, 7u, 253u, 263u, 0u,
0u, 7u, 158u, 370u, 0u, 0u, 7u, 231u, 120u, 0u, 0u, 7u, 182u, 35u, 0u, 0u,
7u, 142u, 194u, 0u, 0u, 7u, 385u, 211u, 0u, 0u, 7u, 209u, 5u, 0u, 0u, 7u,
258u, 82u, 0u, 0u, 7u, 282u, 100u, 0u, 0u, 7u, 243u, 321u, 0u, 0u, 7u, 377u,
405u, 0u, 0u, 7u, 210u, 163u, 0u, 0u, 7u, 165u, 107u, 0u, 0u, 7u, 6u, 153u,
0u, 0u, 0u, 40u, 47u, 87u, 100u, 108u, 113u, 116u, 118u, 119u, 120u, 121u, 139u, 197u,
225u, 241u, 246u, 251u, 315u, 345u, 366u, 384u, 389u, 393u, 394u, 395u, 397u, 399u, 407u, 421u,
439u, 451u, 470u, 487u, 506u, 530u, 578u, 590u, 614u, 626u, 634u, 650u, 654u, 658u, 662u, 663u,
665u, 680u, 688u, 702u, 766u, 830u, 860u, 890u, 909u, 934u, 963u, 1027u, 1091u, 1155u, 1219u, 1259u,
1297u, 1344u, 1385u, 1439u, 1478u, 1542u, 1606u, 1670u, 1733u, 1789u, 1831u, 1893u, 1924u, 1960u, 1985u, 2015u,
2039u, 2062u, 2116u, 2158u, 2190u, 2215u, 2260u, 2284u, 2311u, 2332u, 2347u, 2358u, 2378u, 2393u, 2405u, 2414u,
2429u, 2438u, 2447u, 2453u, 2459u, 2463u, 2466u, 2469u, 2475u, 2481u, 2484u, 2487u, 2490u, 2493u, 2496u, 2499u,
2502u, 2505u, 2506u, 2507u, 2542u, 2557u, 2572u, 2587u,
};
static uint LiteralSchwarzschildTable[2871] = {
557u, 459u, 40u, 135u, 0u, 10u, 86u, 44u, 23u, 123u, 117u, 103u, 104u, 52u, 54u, 44u,
47u, 127u, 122u, 85u, 90u, 91u, 92u, 44u, 58u, 51u, 17u, 29u, 37u, 105u, 39u, 44u,
128u, 132u, 35u, 70u, 40u, 38u, 32u, 44u, 34u, 33u, 28u, 1u, 0u, 0u, 0u, 0u,
1u, 1u, 5u, 0u, 0u, 1u, 2u, 6u, 0u, 0u, 1u, 3u, 7u, 0u, 0u, 1u,
4u, 8u, 0u, 0u, 1u, 5u, 9u, 0u, 0u, 1u, 6u, 10u, 0u, 0u, 1u, 7u,
11u, 0u, 0u, 1u, 8u, 13u, 0u, 0u, 1u, 9u, 14u, 0u, 0u, 1u, 10u, 15u,
0u, 0u, 1u, 11u, 16u, 0u, 0u, 1u, 12u, 17u, 0u, 0u, 1u, 13u, 18u, 0u,
0u, 1u, 14u, 19u, 0u, 0u, 1u, 15u, 21u, 0u, 0u, 1u, 16u, 22u, 0u, 0u,
1u, 17u, 23u, 0u, 0u, 1u, 18u, 24u, 0u, 0u, 1u, 19u, 25u, 0u, 0u, 1u,
20u, 26u, 0u, 0u, 1u, 21u, 27u, 0u, 0u, 1u, 22u, 29u, 0u, 0u, 1u, 23u,
30u, 0u, 0u, 1u, 24u, 31u, 0u, 0u, 1u, 25u, 32u, 0u, 0u, 1u, 26u, 33u,
0u, 0u, 1u, 27u, 34u, 0u, 0u, 1u, 28u, 35u, 0u, 0u, 1u, 29u, 37u, 0u,
0u, 1u, 30u, 38u, 0u, 0u, 1u, 31u, 39u, 0u, 0u, 1u, 32u, 40u, 0u, 0u,
1u, 33u, 41u, 0u, 0u, 1u, 34u, 42u, 0u, 0u, 1u, 35u, 43u, 0u, 0u, 1u,
36u, 44u, 0u, 0u, 4u, 37u, 2u, 36u, 0u, 4u, 38u, 36u, 4u, 0u, 4u, 39u,
9u, 36u, 0u, 4u, 40u, 36u, 11u, 0u, 4u, 41u, 16u, 36u, 0u, 4u, 42u, 36u,
18u, 0u, 4u, 43u, 23u, 36u, 0u, 4u, 44u, 36u, 25u, 0u, 4u, 45u, 30u, 36u,
0u, 4u, 46u, 36u, 32u, 0u, 7u, 30u, 38u, 0u, 0u, 7u, 23u, 4u, 0u, 0u,
7u, 16u, 11u, 0u, 0u, 7u, 9u, 40u, 0u, 0u, 7u, 2u, 18u, 0u, 0u, 7u,
47u, 42u, 0u, 0u, 7u, 48u, 25u, 0u, 0u, 7u, 49u, 44u, 0u, 0u, 7u, 50u,
32u, 0u, 0u, 7u, 51u, 46u, 0u, 0u, 0u, 46u, 1065353216u, 0u, 0u, 0u, 44u, 0u,
0u, 0u, 0u, 42u, 1073741824u, 0u, 0u, 4u, 40u, 1u, 46u, 0u, 4u, 38u, 3u, 46u,
0u, 4u, 32u, 0u, 42u, 0u, 4u, 25u, 8u, 46u, 0u, 4u, 18u, 10u, 46u, 0u,
4u, 11u, 15u, 46u, 0u, 4u, 4u, 17u, 46u, 0u, 4u, 52u, 22u, 46u, 0u, 4u,
53u, 24u, 46u, 0u, 4u, 54u, 29u, 46u, 0u, 4u, 55u, 31u, 46u, 0u, 9u, 31u,
40u, 37u, 0u, 9u, 29u, 38u, 44u, 0u, 9u, 24u, 31u, 29u, 0u, 8u, 29u, 24u,
0u, 0u, 5u, 24u, 40u, 29u, 0u, 5u, 31u, 37u, 29u, 0u, 5u, 22u, 38u, 29u,
0u, 5u, 17u, 25u, 29u, 0u, 5u, 15u, 39u, 29u, 0u, 5u, 10u, 18u, 29u, 0u,
5u, 8u, 11u, 29u, 0u, 5u, 3u, 41u, 29u, 0u, 5u, 1u, 4u, 29u, 0u, 5u,
0u, 52u, 29u, 0u, 5u, 56u, 43u, 29u, 0u, 5u, 57u, 53u, 29u, 0u, 5u, 58u,
54u, 29u, 0u, 5u, 59u, 45u, 29u, 0u, 5u, 60u, 55u, 29u, 0u, 4u, 61u, 24u,
24u, 0u, 4u, 62u, 31u, 31u, 0u, 4u, 63u, 22u, 22u, 0u, 4u, 64u, 17u, 24u,
0u, 4u, 65u, 24u, 17u, 0u, 4u, 66u, 15u, 31u, 0u, 4u, 67u, 31u, 15u, 0u,
4u, 68u, 10u, 22u, 0u, 4u, 69u, 22u, 10u, 0u, 4u, 70u, 8u, 24u, 0u, 4u,
71u, 24u, 8u, 0u, 4u, 72u, 3u, 31u, 0u, 4u, 73u, 31u, 3u, 0u, 4u, 74u,
1u, 22u, 0u, 4u, 75u, 22u, 1u, 0u, 4u, 76u, 0u, 24u, 0u, 4u, 77u, 24u,
0u, 0u, 4u, 78u, 56u, 31u, 0u, 4u, 79u, 31u, 56u, 0u, 4u, 80u, 57u, 22u,
0u, 4u, 81u, 22u, 57u, 0u, 4u, 82u, 58u, 24u, 0u, 4u, 83u, 24u, 58u, 0u,
4u, 84u, 59u, 31u, 0u, 4u, 85u, 31u, 59u, 0u, 4u, 86u, 60u, 22u, 0u, 4u,
87u, 22u, 60u, 0u, 2u, 60u, 61u, 62u, 0u, 2u, 59u, 64u, 65u, 0u, 2u, 58u,
66u, 67u, 0u, 2u, 57u, 68u, 69u, 0u, 2u, 56u, 70u, 71u, 0u, 2u, 0u, 72u,
73u, 0u, 2u, 1u, 74u, 75u, 0u, 2u, 3u, 76u, 77u, 0u, 2u, 8u, 78u, 79u,
0u, 2u, 10u, 80u, 81u, 0u, 2u, 15u, 82u, 83u, 0u, 2u, 17u, 84u, 85u, 0u,
2u, 22u, 86u, 87u, 0u, 2u, 87u, 60u, 63u, 0u, 2u, 86u, 59u, 58u, 0u, 2u,
85u, 56u, 0u, 0u, 2u, 84u, 3u, 8u, 0u, 2u, 83u, 15u, 17u, 0u, 2u, 17u,
86u, 57u, 0u, 2u, 15u, 85u, 1u, 0u, 2u, 8u, 84u, 10u, 0u, 2u, 3u, 83u,
22u, 0u, 6u, 83u, 87u, 0u, 0u, 4u, 87u, 42u, 83u, 0u, 4u, 84u, 29u, 83u,
0u, 5u, 83u, 40u, 84u, 0u, 5u, 85u, 37u, 84u, 0u, 5u, 86u, 38u, 84u, 0u,
5u, 22u, 32u, 84u, 0u, 5u, 10u, 17u, 87u, 0u, 5u, 1u, 15u, 87u, 0u, 5u,
57u, 8u, 87u, 0u, 5u, 0u, 3u, 87u, 0u, 4u, 87u, 83u, 46u, 0u, 4u, 3u,
85u, 36u, 0u, 4u, 8u, 86u, 46u, 0u, 4u, 15u, 36u, 22u, 0u, 4u, 17u, 29u,
10u, 0u, 4u, 32u, 29u, 1u, 0u, 4u, 38u, 29u, 57u, 0u, 4u, 40u, 29u, 0u,
0u, 4u, 0u, 87u, 5u, 0u, 4u, 57u, 3u, 6u, 0u, 4u, 1u, 8u, 7u, 0u,
4u, 10u, 22u, 87u, 0u, 4u, 29u, 22u, 3u, 0u, 4u, 37u, 22u, 8u, 0u, 4u,
56u, 22u, 17u, 0u, 4u, 58u, 83u, 17u, 0u, 4u, 59u, 87u, 12u, 0u, 4u, 60u,
85u, 17u, 0u, 4u, 63u, 3u, 13u, 0u, 4u, 82u, 86u, 17u, 0u, 4u, 81u, 8u,
14u, 0u, 4u, 80u, 22u, 32u, 0u, 4u, 79u, 83u, 32u, 0u, 4u, 78u, 87u, 19u,
0u, 4u, 77u, 85u, 32u, 0u, 4u, 76u, 3u, 20u, 0u, 4u, 75u, 86u, 32u, 0u,
4u, 74u, 8u, 21u, 0u, 4u, 73u, 22u, 38u, 0u, 4u, 72u, 83u, 38u, 0u, 4u,
71u, 87u, 26u, 0u, 4u, 70u, 85u, 38u, 0u, 4u, 69u, 3u, 27u, 0u, 4u, 68u,
86u, 38u, 0u, 4u, 67u, 8u, 28u, 0u, 4u, 66u, 22u, 40u, 0u, 4u, 65u, 83u,
40u, 0u, 4u, 64u, 87u, 33u, 0u, 4u, 62u, 85u, 40u, 0u, 4u, 61u, 3u, 34u,
0u, 4u, 31u, 86u, 40u, 0u, 4u, 24u, 8u, 35u, 0u, 7u, 86u, 56u, 0u, 0u,
7u, 85u, 58u, 0u, 0u, 7u, 83u, 60u, 0u, 0u, 7u, 88u, 82u, 0u, 0u, 7u,
89u, 80u, 0u, 0u, 7u, 90u, 79u, 0u, 0u, 7u, 91u, 77u, 0u, 0u, 7u, 92u,
75u, 0u, 0u, 7u, 93u, 73u, 0u, 0u, 7u, 94u, 72u, 0u, 0u, 7u, 95u, 70u,
0u, 0u, 7u, 96u, 68u, 0u, 0u, 7u, 97u, 66u, 0u, 0u, 7u, 98u, 65u, 0u,
0u, 7u, 99u, 62u, 0u, 0u, 7u, 100u, 31u, 0u, 0u, 2u, 31u, 44u, 0u, 0u,
2u, 62u, 25u, 85u, 0u, 2u, 65u, 39u, 83u, 0u, 2u, 66u, 18u, 88u, 0u, 2u,
68u, 11u, 90u, 0u, 2u, 70u, 41u, 91u, 0u, 2u, 72u, 4u, 92u, 0u, 2u, 73u,
52u, 94u, 0u, 2u, 75u, 43u, 95u, 0u, 2u, 77u, 53u, 96u, 0u, 2u, 79u, 54u,
98u, 0u, 2u, 80u, 45u, 99u, 0u, 2u, 82u, 55u, 100u, 0u, 5u, 100u, 22u, 84u,
0u, 5u, 99u, 86u, 84u, 0u, 5u, 98u, 62u, 84u, 0u, 5u, 96u, 65u, 84u, 0u,
5u, 95u, 66u, 84u, 0u, 5u, 94u, 89u, 84u, 0u, 5u, 92u, 68u, 84u, 0u, 5u,
91u, 70u, 84u, 0u, 5u, 90u, 72u, 84u, 0u, 5u, 88u, 93u, 84u, 0u, 5u, 83u,
73u, 84u, 0u, 5u, 85u, 75u, 84u, 0u, 5u, 0u, 77u, 84u, 0u, 5u, 55u, 97u,
84u, 0u, 5u, 54u, 79u, 84u, 0u, 5u, 53u, 80u, 84u, 0u, 5u, 52u, 82u, 84u,
0u, 4u, 82u, 36u, 99u, 0u, 4u, 80u, 98u, 46u, 0u, 4u, 79u, 96u, 36u, 0u,
4u, 77u, 95u, 46u, 0u, 4u, 75u, 99u, 87u, 0u, 4u, 73u, 99u, 3u, 0u, 4u,
72u, 99u, 8u, 0u, 4u, 70u, 100u, 17u, 0u, 4u, 68u, 36u, 94u, 0u, 4u, 66u,
92u, 46u, 0u, 4u, 65u, 91u, 36u, 0u, 4u, 62u, 90u, 46u, 0u, 4u, 97u, 94u,
87u, 0u, 4u, 93u, 94u, 3u, 0u, 4u, 89u, 94u, 8u, 0u, 4u, 86u, 100u, 32u,
0u, 4u, 4u, 36u, 88u, 0u, 4u, 11u, 83u, 46u, 0u, 4u, 18u, 85u, 36u, 0u,
4u, 25u, 0u, 46u, 0u, 4u, 45u, 88u, 87u, 0u, 4u, 43u, 88u, 3u, 0u, 4u,
41u, 88u, 8u, 0u, 4u, 39u, 100u, 38u, 0u, 4u, 60u, 36u, 55u, 0u, 4u, 58u,
54u, 46u, 0u, 4u, 56u, 53u, 36u, 0u, 4u, 101u, 52u, 46u, 0u, 4u, 102u, 55u,
87u, 0u, 4u, 103u, 55u, 3u, 0u, 4u, 104u, 55u, 8u, 0u, 4u, 105u, 100u, 40u,
0u, 4u, 52u, 80u, 5u, 0u, 4u, 53u, 79u, 6u, 0u, 4u, 54u, 77u, 7u, 0u,
4u, 0u, 22u, 80u, 0u, 4u, 85u, 22u, 79u, 0u, 4u, 83u, 22u, 77u, 0u, 4u,
90u, 66u, 5u, 0u, 4u, 91u, 65u, 6u, 0u, 4u, 92u, 62u, 7u, 0u, 4u, 95u,
22u, 66u, 0u, 4u, 96u, 22u, 65u, 0u, 4u, 98u, 22u, 62u, 0u, 4u, 40u, 11u,
5u, 0u, 4u, 38u, 18u, 6u, 0u, 4u, 32u, 25u, 7u, 0u, 4u, 17u, 22u, 11u,
0u, 4u, 46u, 22u, 18u, 0u, 4u, 36u, 22u, 25u, 0u, 4u, 106u, 58u, 5u, 0u,
4u, 107u, 56u, 6u, 0u, 4u, 108u, 101u, 7u, 0u, 4u, 109u, 22u, 58u, 0u, 4u,
110u, 22u, 56u, 0u, 4u, 111u, 22u, 101u, 0u, 2u, 22u, 31u, 57u, 0u, 2u, 112u,
52u, 59u, 0u, 2u, 113u, 53u, 63u, 0u, 2u, 114u, 54u, 81u, 0u, 2u, 115u, 75u,
0u, 0u, 2u, 116u, 73u, 85u, 0u, 2u, 117u, 72u, 83u, 0u, 2u, 118u, 90u, 78u,
0u, 2u, 119u, 91u, 76u, 0u, 2u, 120u, 92u, 74u, 0u, 2u, 121u, 97u, 95u, 0u,
2u, 122u, 93u, 96u, 0u, 2u, 123u, 89u, 98u, 0u, 2u, 124u, 40u, 71u, 0u, 2u,
125u, 38u, 69u, 0u, 2u, 126u, 32u, 67u, 0u, 2u, 127u, 45u, 17u, 0u, 2u, 128u,
43u, 46u, 0u, 2u, 129u, 41u, 36u, 0u, 2u, 130u, 106u, 64u, 0u, 2u, 131u, 107u,
61u, 0u, 2u, 132u, 108u, 24u, 0u, 2u, 133u, 102u, 109u, 0u, 2u, 134u, 103u, 110u,
0u, 2u, 135u, 104u, 111u, 0u, 2u, 111u, 22u, 1u, 0u, 2u, 110u, 112u, 113u, 0u,
2u, 109u, 118u, 119u, 0u, 2u, 108u, 124u, 125u, 0u, 2u, 107u, 130u, 131u, 0u, 2u,
131u, 30u, 111u, 0u, 2u, 130u, 110u, 114u, 0u, 2u, 125u, 109u, 120u, 0u, 2u, 124u,
108u, 126u, 0u, 2u, 119u, 107u, 132u, 0u, 4u, 107u, 15u, 131u, 0u, 4u, 108u, 10u,
131u, 0u, 4u, 109u, 29u, 131u, 0u, 4u, 110u, 37u, 131u, 0u, 4u, 132u, 100u, 131u,
0u, 4u, 126u, 82u, 131u, 0u, 4u, 120u, 115u, 131u, 0u, 4u, 114u, 116u, 131u, 0u,
4u, 30u, 117u, 131u, 0u, 4u, 118u, 68u, 131u, 0u, 4u, 113u, 121u, 131u, 0u, 4u,
112u, 122u, 131u, 0u, 4u, 22u, 123u, 131u, 0u, 4u, 1u, 4u, 131u, 0u, 4u, 106u,
127u, 131u, 0u, 4u, 36u, 128u, 131u, 0u, 4u, 46u, 129u, 131u, 0u, 4u, 17u, 60u,
131u, 0u, 4u, 32u, 133u, 131u, 0u, 4u, 38u, 134u, 131u, 0u, 4u, 40u, 135u, 131u,
0u, 2u, 135u, 23u, 107u, 0u, 2u, 134u, 9u, 130u, 0u, 2u, 133u, 47u, 125u, 0u,
2u, 129u, 49u, 124u, 0u, 2u, 128u, 51u, 119u, 0u, 4u, 107u, 15u, 134u, 0u, 4u,
51u, 10u, 134u, 0u, 4u, 49u, 29u, 134u, 0u, 4u, 47u, 37u, 134u, 0u, 4u, 9u,
100u, 134u, 0u, 4u, 23u, 15u, 133u, 0u, 4u, 127u, 10u, 133u, 0u, 4u, 123u, 29u,
133u, 0u, 4u, 122u, 37u, 133u, 0u, 4u, 121u, 100u, 133u, 0u, 4u, 117u, 15u, 129u,
0u, 4u, 116u, 10u, 129u, 0u, 4u, 115u, 29u, 129u, 0u, 4u, 60u, 37u, 129u, 0u,
4u, 4u, 100u, 129u, 0u, 4u, 68u, 15u, 128u, 0u, 4u, 82u, 10u, 128u, 0u, 4u,
98u, 29u, 128u, 0u, 4u, 96u, 37u, 128u, 0u, 4u, 95u, 100u, 128u, 0u, 2u, 100u,
126u, 107u, 0u, 2u, 37u, 120u, 51u, 0u, 2u, 29u, 114u, 49u, 0u, 2u, 10u, 30u,
47u, 0u, 2u, 15u, 118u, 23u, 0u, 2u, 92u, 113u, 127u, 0u, 2u, 91u, 112u, 123u,
0u, 2u, 90u, 22u, 122u, 0u, 2u, 83u, 1u, 117u, 0u, 2u, 85u, 106u, 116u, 0u,
2u, 0u, 36u, 115u, 0u, 2u, 54u, 46u, 60u, 0u, 2u, 53u, 17u, 68u, 0u, 2u,
52u, 32u, 82u, 0u, 2u, 104u, 38u, 98u, 0u, 2u, 103u, 40u, 96u, 0u, 7u, 96u,
37u, 0u, 0u, 7u, 98u, 29u, 0u, 0u, 7u, 82u, 10u, 0u, 0u, 7u, 68u, 70u,
0u, 0u, 7u, 60u, 92u, 0u, 0u, 7u, 115u, 91u, 0u, 0u, 7u, 116u, 90u, 0u,
0u, 7u, 117u, 86u, 0u, 0u, 7u, 122u, 85u, 0u, 0u, 7u, 123u, 0u, 0u, 0u,
7u, 127u, 54u, 0u, 0u, 7u, 23u, 39u, 0u, 0u, 7u, 47u, 52u, 0u, 0u, 7u,
49u, 104u, 0u, 0u, 7u, 51u, 103u, 0u, 0u, 7u, 107u, 105u, 0u, 0u, 2u, 103u,
16u, 100u, 0u, 2u, 104u, 12u, 96u, 0u, 2u, 52u, 13u, 98u, 0u, 2u, 54u, 14u,
82u, 0u, 2u, 0u, 99u, 68u, 0u, 2u, 85u, 2u, 15u, 0u, 2u, 90u, 19u, 60u,
0u, 2u, 91u, 20u, 115u, 0u, 2u, 92u, 21u, 116u, 0u, 2u, 10u, 94u, 117u, 0u,
2u, 29u, 48u, 83u, 0u, 2u, 37u, 26u, 122u, 0u, 2u, 105u, 27u, 123u, 0u, 2u,
39u, 28u, 127u, 0u, 2u, 86u, 88u, 23u, 0u, 2u, 70u, 50u, 53u, 0u, 2u, 40u,
33u, 47u, 0u, 2u, 38u, 34u, 49u, 0u, 2u, 32u, 35u, 51u, 0u, 2u, 17u, 55u,
107u, 0u, 5u, 107u, 131u, 42u, 0u, 5u, 51u, 0u, 84u, 0u, 5u, 49u, 134u, 42u,
0u, 5u, 47u, 10u, 84u, 0u, 5u, 23u, 133u, 42u, 0u, 5u, 127u, 86u, 84u, 0u,
5u, 123u, 129u, 42u, 0u, 5u, 122u, 17u, 84u, 0u, 5u, 117u, 128u, 42u, 0u, 2u,
17u, 111u, 107u, 0u, 2u, 86u, 130u, 49u, 0u, 2u, 10u, 125u, 23u, 0u, 2u, 0u,
124u, 123u, 0u, 2u, 128u, 119u, 117u, 0u, 4u, 117u, 17u, 87u, 0u, 4u, 123u, 17u,
3u, 0u, 4u, 23u, 17u, 8u, 0u, 4u, 49u, 51u, 131u, 0u, 4u, 107u, 86u, 87u,
0u, 4u, 119u, 17u, 80u, 0u, 4u, 124u, 86u, 3u, 0u, 4u, 125u, 17u, 79u, 0u,
4u, 130u, 86u, 8u, 0u, 4u, 111u, 17u, 77u, 0u, 4u, 129u, 47u, 131u, 0u, 4u,
133u, 10u, 87u, 0u, 4u, 134u, 17u, 66u, 0u, 4u, 84u, 10u, 3u, 0u, 4u, 42u,
17u, 65u, 0u, 4u, 116u, 10u, 8u, 0u, 4u, 115u, 17u, 62u, 0u, 4u, 60u, 127u,
131u, 0u, 4u, 68u, 0u, 87u, 0u, 4u, 82u, 17u, 11u, 0u, 4u, 98u, 0u, 3u,
0u, 4u, 96u, 17u, 18u, 0u, 4u, 53u, 0u, 8u, 0u, 4u, 83u, 17u, 25u, 0u,
4u, 15u, 122u, 131u, 0u, 4u, 100u, 128u, 87u, 0u, 4u, 55u, 17u, 58u, 0u, 4u,
88u, 128u, 3u, 0u, 4u, 94u, 17u, 56u, 0u, 4u, 99u, 128u, 8u, 0u, 4u, 50u,
17u, 101u, 0u, 2u, 128u, 49u, 9u, 0u, 2u, 0u, 107u, 119u, 0u, 2u, 10u, 124u,
125u, 0u, 2u, 86u, 130u, 111u, 0u, 2u, 17u, 129u, 121u, 0u, 2u, 122u, 133u, 134u,
0u, 2u, 127u, 84u, 42u, 0u, 2u, 47u, 116u, 115u, 0u, 2u, 51u, 60u, 4u, 0u,
2u, 131u, 68u, 82u, 0u, 2u, 101u, 98u, 96u, 0u, 2u, 56u, 53u, 83u, 0u, 2u,
58u, 15u, 95u, 0u, 2u, 25u, 100u, 55u, 0u, 2u, 18u, 88u, 94u, 0u, 2u, 11u,
99u, 50u, 0u, 7u, 50u, 0u, 0u, 0u, 7u, 99u, 10u, 0u, 0u, 7u, 94u, 86u,
0u, 0u, 7u, 88u, 122u, 0u, 0u, 7u, 55u, 127u, 0u, 0u, 7u, 100u, 47u, 0u,
0u, 7u, 15u, 131u, 0u, 0u, 7u, 83u, 101u, 0u, 0u, 7u, 53u, 56u, 0u, 0u,
7u, 96u, 25u, 0u, 0u, 7u, 98u, 18u, 0u, 0u, 7u, 82u, 11u, 0u, 0u, 2u,
11u, 12u, 50u, 0u, 2u, 18u, 13u, 99u, 0u, 2u, 25u, 14u, 94u, 0u, 2u, 56u,
19u, 88u, 0u, 2u, 101u, 20u, 55u, 0u, 2u, 131u, 21u, 100u, 0u, 2u, 47u, 26u,
15u, 0u, 2u, 127u, 27u, 83u, 0u, 2u, 122u, 28u, 53u, 0u, 2u, 86u, 33u, 96u,
0u, 2u, 10u, 34u, 98u, 0u, 2u, 0u, 35u, 82u, 0u, 4u, 82u, 132u, 11u, 0u,
4u, 98u, 132u, 18u, 0u, 4u, 96u, 132u, 25u, 0u, 4u, 53u, 132u, 56u, 0u, 4u,
83u, 132u, 101u, 0u, 4u, 15u, 132u, 131u, 0u, 4u, 100u, 132u, 47u, 0u, 4u, 55u,
132u, 127u, 0u, 4u, 88u, 132u, 122u, 0u, 4u, 94u, 132u, 86u, 0u, 4u, 99u, 132u,
10u, 0u, 4u, 50u, 132u, 0u, 0u, 3u, 0u, 5u, 108u, 0u, 3u, 10u, 6u, 109u,
0u, 3u, 86u, 7u, 110u, 0u, 3u, 122u, 5u, 117u, 0u, 3u, 127u, 6u, 123u, 0u,
3u, 47u, 7u, 23u, 0u, 4u, 23u, 132u, 122u, 0u, 4u, 123u, 132u, 127u, 0u, 4u,
117u, 132u, 47u, 0u, 4u, 110u, 128u, 122u, 0u, 4u, 109u, 128u, 127u, 0u, 4u, 108u,
128u, 47u, 0u, 4u, 7u, 17u, 122u, 0u, 4u, 6u, 17u, 127u, 0u, 4u, 5u, 17u,
47u, 0u, 4u, 131u, 51u, 122u, 0u, 4u, 101u, 51u, 127u, 0u, 4u, 56u, 51u, 47u,
0u, 4u, 25u, 58u, 122u, 0u, 4u, 18u, 58u, 127u, 0u, 4u, 11u, 58u, 47u, 0u,
2u, 47u, 110u, 82u, 0u, 2u, 127u, 109u, 98u, 0u, 2u, 122u, 108u, 96u, 0u, 2u,
58u, 7u, 53u, 0u, 2u, 51u, 6u, 83u, 0u, 2u, 17u, 5u, 15u, 0u, 2u, 128u,
131u, 100u, 0u, 2u, 132u, 101u, 55u, 0u, 2u, 35u, 56u, 88u, 0u, 2u, 34u, 25u,
94u, 0u, 2u, 33u, 18u, 99u, 0u, 2u, 28u, 11u, 50u, 0u, 0u, 37u, 47u, 57u,
60u, 71u, 73u, 74u, 75u, 90u, 117u, 130u, 135u, 139u, 140u, 142u, 150u, 158u, 192u, 208u,
221u, 238u, 270u, 294u, 319u, 324u, 329u, 350u, 355u, 375u, 391u, 407u, 427u, 436u, 441u, 472u,
488u, 500u, 512u, 524u, 530u, 545u, 557u,
};
groupshared uint literalTransportMatch;
bool RecognizeLiteralTransport(uint programBase, bool schwarzschild, uint lane) {
if (lane == 0) literalTransportMatch = 1;
GroupMemoryBarrierWithGroupSync();
if (schwarzschild) {
for (uint i = lane; i < 2871u; i += 16u)
    if (inputs[programBase + i] != LiteralSchwarzschildTable[i])
        InterlockedAnd(literalTransportMatch, 0u);
} else {
for (uint i = lane; i < 13096u; i += 16u)
    if (inputs[programBase + i] != LiteralGeneralTable[i])
        InterlockedAnd(literalTransportMatch, 0u);
}
GroupMemoryBarrierWithGroupSync();
return literalTransportMatch != 0;
}
uint CanonicalTransportWord(uint offset, bool schwarzschild) {
if (schwarzschild) return LiteralSchwarzschildTable[offset];
return LiteralGeneralTable[offset];
}
bool EvaluateOneCanonical(uint index, uint row, uint rowBase, bool schwarzschild, uint registers) {
    uint address = 43 + 5 * index, op = CanonicalTransportWord(address, schwarzschild),
         destination = CanonicalTransportWord(address + 1, schwarzschild), a = CanonicalTransportWord(address + 2, schwarzschild), b = CanonicalTransportWord(address + 3, schwarzschild),
         c = CanonicalTransportWord(address + 4, schwarzschild);
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

bool EvaluateProgramCanonical(uint row, uint rowBase, bool schwarzschild, uint count, uint registers,
                     uint layers, uint stage, uint lane) {
    uint needed = count;
    if (needed == 0 || needed > count || layers == 0 || layers > SIRIUS_RETAINED_LAYERS) return false;
    if (lane == 0) retainedStatus = 0;
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
    GroupMemoryBarrierWithGroupSync();
#endif
    uint layerBase = 43 + 5 * count;
    if (CanonicalTransportWord(layerBase, schwarzschild) != 0 || CanonicalTransportWord(layerBase + layers, schwarzschild) != count) return false;
    uint expected = 0;
    for (uint layer = 0; layer < layers; ++layer) {
        uint start = CanonicalTransportWord(layerBase + layer, schwarzschild), end = CanonicalTransportWord(layerBase + layer + 1, schwarzschild);
        if (start >= needed) break;
        if (start != expected || end <= start || end > needed ||
            end - start > SIRIUS_RETAINED_LANES)
            RejectRetainedProgram();
#if defined(SIRIUS_RETAINED_PORTABLE) && !defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
        else {
            // Layer destinations are disjoint from every input register and
            // from one another; the generated plan releases inputs afterwards.
            for (uint node = start; node < end; ++node)
                if (!EvaluateOneCanonical(node, row, rowBase, schwarzschild, registers)) RejectRetainedProgram();
        }
#elif defined(SIRIUS_RETAINED_PORTABLE) && defined(SIRIUS_RETAINED_NORMAL_SUM32) && defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
        else {
            // Cover the original complete layer before its publication barrier.
            for (uint node = start + lane; node < end;
                 node += SIRIUS_RETAINED_EXECUTION_LANES)
                if (!EvaluateOneCanonical(node, row, rowBase, schwarzschild, registers))
                    RejectRetainedProgram();
        }
#else
        else if (start + lane < end &&
                 !EvaluateOneCanonical(start + lane, row, rowBase, schwarzschild, registers))
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
        uint index = CanonicalTransportWord(3 + i, schwarzschild);
        if (index >= registers || ReadValue(index, rowBase).valid == 0) return false;
    }
    if (lane == 0) {
        for (uint i = 0; i < 40; ++i)
            WriteStored(rowBase + 1004 + stage * 200 + i * 5,
                        ReadValue(CanonicalTransportWord(3 + i, schwarzschild), rowBase));
    }
#else
    // Every root must pass before any stage value is published.
    for (uint i = lane; i < 40; i += SIRIUS_RETAINED_EXECUTION_LANES) {
        uint index = CanonicalTransportWord(3 + i, schwarzschild);
        if (index >= registers)
            RejectRetainedProgram();
        else if (ScratchRead(index, 4, rowBase, 2404) == 0)
            RejectRetainedProgram();
    }
    GroupMemoryBarrierWithGroupSync();
    if (retainedStatus != 0) return false;
    for (uint i = lane; i < 40; i += SIRIUS_RETAINED_EXECUTION_LANES)
        WriteStored(rowBase + 1004 + stage * 200 + i * 5,
                    ReadValue(CanonicalTransportWord(3 + i, schwarzschild), rowBase));
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
bool literalSchwarzschild = programBase != 1u + 230u * rows;
    bool literalProgram = RecognizeLiteralTransport(programBase, literalSchwarzschild, lane);
#if !defined(SIRIUS_RETAINED_PORTABLE) || defined(SIRIUS_RETAINED_PARALLEL_TRANSPORT)
    // The immutable plan leaves never-assigned slots invalid across all stages.
    // The existing stage-zero barrier publishes this disjoint validity clear.
    for (uint index = lane; index < SIRIUS_RETAINED_REGISTERS;
         index += SIRIUS_RETAINED_EXECUTION_LANES)
        ScratchWrite(index, 4, rowBase, 2404, 0);
#endif
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
        if (literalProgram) {
            if (!EvaluateProgramCanonical(row, rowBase, literalSchwarzschild, count, registers, layers, stage, lane)) return;
        } else {
            if (!EvaluateProgramGeneric(row, rowBase, programBase, count, registers, layers, stage, lane)) return;
        }
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
