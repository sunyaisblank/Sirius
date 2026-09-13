// Bounded retained arithmetic stream; development probe, not production launch.
import retained_pair;
[[vk::binding(0, 0)]] StructuredBuffer<uint> inputs;
[[vk::binding(1, 0)]] RWStructuredBuffer<uint> outputs;
RetainedPair ReadValue(uint index) {
    uint base = 384 + 4 * index;
    RetainedPair r;
    r.hi = asfloat(outputs[base]);
    r.lo = asfloat(outputs[base + 1]);
    r.error = asfloat(outputs[base + 2]);
    r.valid = outputs[base + 3];
    return r;
}
void WriteValue(uint index, RetainedPair r) {
    uint base = 384 + 4 * index;
    outputs[base] = asuint(r.hi);
    outputs[base + 1] = asuint(r.lo);
    outputs[base + 2] = asuint(r.error);
    outputs[base + 3] = r.valid;
}
[shader("compute")][numthreads(1, 1, 1)] void ComputeMain(uint3 id
                                                          : SV_DispatchThreadID) {
    if (any(id != uint3(0))) return;
    uint ni, si, no, so;
    inputs.GetDimensions(ni, si);
    outputs.GetDimensions(no, so);
    if (no < 384) return;
    for (uint i = 0; i < no; ++i) outputs[i] = 0;
    if (ni < 138) return;
    uint count = inputs[32], registers = inputs[33];
    if (count > 65536 || registers > 8192 || ni != 138 + 5 * count || no != 384 + 4 * registers)
        return;
    outputs[0] = asuint(271832.0f);
    outputs[1] = asuint(32.0f);
#ifdef SIRIUS_FP32_COMP
    outputs[2] = asuint(1.0f);
#endif
    outputs[4] = inputs[30];
    outputs[5] = asuint(float(ni));
    outputs[6] = asuint(float(no));
    for (uint i = 0; i < 32; ++i) {
        outputs[8 + i] = inputs[i];
        if ((inputs[i] & 0x7fffffffu) >= 0x7f800000u) return;
    }
    if (asfloat(inputs[31]) != 0) return;
    for (uint i = 20; i <= 24; ++i)
        if (asfloat(inputs[i]) <= 0) return;
    float thin = asfloat(inputs[29]);
    if (thin != 1 && (thin != 0 || asfloat(inputs[27]) != 0 || asfloat(inputs[28]) != 0)) return;
    if (asfloat(inputs[1]) != 0 && asfloat(inputs[3]) != 0) return;
    for (uint index = 0; index < count; ++index) {
        uint address = 138 + 5 * index, op = inputs[address], destination = inputs[address + 1],
             a = inputs[address + 2], b = inputs[address + 3], c = inputs[address + 4];
        if (destination >= registers || op > 11) return;
        RetainedPair r = RPInvalid();
        if (op == 0)
            r = RPConstant(asfloat(a));
        else if (op == 1) {
            if (a >= 32) return;
            r = RPConstant(asfloat(inputs[a]));
        } else {
            if (a >= registers) return;
            RetainedPair av = ReadValue(a);
            if (op == 6)
                r = RPSqrt(av);
            else if (op == 7)
                r = RPNegate(av);
            else if (op == 8) {
                float bound = RPAbsoluteUpper(av);
                if (av.valid != 0 && bound > 0)
                    r = RPConstant(RPPowerOfTwo(RPFloorExponent(bound) + 1));
            } else {
                if (b >= registers) return;
                RetainedPair bv = ReadValue(b);
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
                    if (c >= registers) return;
                    RetainedPair cv = ReadValue(c);
                    if (cv.valid != 0) {
                        bool choose =
                            op == 10 ? cv.hi >= 0 : (cv.hi == 0 && cv.lo == 0 && cv.error == 0);
                        r = choose ? av : bv;
                    }
                }
            }
        }
        WriteValue(destination, r);
    }
    for (uint i = 0; i < 104; ++i) {
        uint index = inputs[34 + i];
        if (index >= registers || ReadValue(index).valid == 0) return;
    }
    for (uint i = 0; i < 104; ++i) {
        uint index = inputs[34 + i];
        RetainedPair r = ReadValue(index);
        outputs[40 + 3 * i] = asuint(r.hi);
        outputs[41 + 3 * i] = asuint(r.lo);
        outputs[42 + 3 * i] = asuint(r.error);
    }
    outputs[3] = asuint(1.0f);
}
