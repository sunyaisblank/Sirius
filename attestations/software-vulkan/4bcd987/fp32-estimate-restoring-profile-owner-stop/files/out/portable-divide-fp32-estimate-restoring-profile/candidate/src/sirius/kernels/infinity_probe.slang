// Test-only bounded Kerr continuation entry. Each dispatched lane owns one
// eight-row input, opaque continuation, and five-row result. The host must
// dispatch exactly the allocated ray count and preserve state between attempts.
import gr_types;
import gr_kerr_infinity;

[[vk::binding(0, 0)]] StructuredBuffer<Real4> inputs;
[[vk::binding(1, 0)]] RWStructuredBuffer<KerrInfinityStateGPU> continuations;
[[vk::binding(2, 0)]] RWStructuredBuffer<Real4> outputs;

Vec4Cart ReadVector(uint offset) {
    Real4 v = inputs[offset];
    return MakeVec4Cart(Real(v.x), Real(v.y), Real(v.z), Real(v.w));
}

[shader("compute")]
[numthreads(1, 1, 1)]
void ComputeMain(uint3 id : SV_DispatchThreadID) {
    KerrInfinityStateGPU state;
    uint inputOffset = id.x * 8;
    bool directFinish = inputs[inputOffset].z == Real(2);
    if (directFinish) {
        // Direct publication-boundary controls construct finite completed or
        // unavailable states without requiring a tail to overflow first.
        state = (KerrInfinityStateGPU)0;
        state.status = uint(inputs[inputOffset].x);
        state.fraction = inputs[inputOffset].y;
        state.frequency = inputs[inputOffset + 1].xyz;
        state.maximumLocalErrorRatio = inputs[inputOffset + 1].w;
        state.acceptedSteps = 7;
        state.attemptedSteps = 9;
        for (uint i = 0; i < 6; ++i) state.values[i] = inputs[inputOffset + 2 + i].xyz;
    } else if (inputs[inputOffset].z == 0) {
        KerrInfinityConfigGPU config;
        config.absoluteTolerance = inputs[inputOffset + 1].x;
        config.relativeTolerance = inputs[inputOffset + 1].y;
        config.inputConstraintTolerance = inputs[inputOffset + 1].z;
        config.maximumAttempts = uint(inputs[inputOffset + 1].w);
        Vec4Cart X[2] = {ReadVector(inputOffset + 4), ReadVector(inputOffset + 5)};
        Vec4Cart V[2] = {ReadVector(inputOffset + 6), ReadVector(inputOffset + 7)};
        InitKerrInfinity(float(inputs[inputOffset].x), float(inputs[inputOffset].y),
                         ReadVector(inputOffset + 2), ReadVector(inputOffset + 3),
                         X, V, inputs[inputOffset].w, config, state);
    } else {
        state = continuations[id.x];
        AttemptKerrInfinityStep(state);
    }
    continuations[id.x] = state;
    KerrInfinityResultGPU result;
    bool complete = FinishKerrInfinity(state, result);
    uint offset = id.x * 5;
    outputs[offset] = Real4(Real(state.status), Real(state.attemptedSteps), Real(state.acceptedSteps), state.fraction);
    outputs[offset + 1] = Real4(result.direction, result.determinant);
    outputs[offset + 2] = Real4(result.jacobian[0], result.jacobian[1]);
    outputs[offset + 3] = Real4(result.frequency, result.frequencyDerivative, result.maximumLocalErrorRatio);
    outputs[offset + 4] = directFinish
        ? Real4(complete ? Real(1) : Real(0), Real(result.acceptedSteps), Real(result.attemptedSteps), state.frequency.x)
        : Real4(complete ? Real(1) : Real(0), state.step, Real(state.config.maximumAttempts), state.frequency.x);
}
