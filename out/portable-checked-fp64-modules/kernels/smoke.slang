// Backend smoke kernel: f(r) = 1 - 2M/r per element. Exercises the Slang
// interface and generic machinery the physics kernels rely on, and pins
// the host/kernel ABI (storage buffers bound 0..N-1 in declaration order,
// entry point renamed to "main" in SPIR-V).

interface IMetricFactor {
    float Factor(float r);
}

struct Schwarzschild : IMetricFactor {
    float mass;
    float Factor(float r) { return 1.0f - 2.0f * mass / r; }
}

float ApplyFactor<T : IMetricFactor>(T metric, float r) {
    return metric.Factor(r);
}

StructuredBuffer<float> radii;       // binding 0
RWStructuredBuffer<float> factors;   // binding 1
StructuredBuffer<float> params;      // binding 2: [0] mass, [1] element count

[shader("compute")]
[numthreads(64, 1, 1)]
void ComputeMain(uint3 tid: SV_DispatchThreadID) {
    uint count = uint(params[1]);
    if (tid.x >= count) {
        return;
    }
    Schwarzschild metric;
    metric.mass = params[0];
    factors[tid.x] = ApplyFactor(metric, radii[tid.x]);
}
