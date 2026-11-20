import gr_types;
import gr_trace_event;

StructuredBuffer<Real> inputs;
RWStructuredBuffer<Real> outputs;
StructuredBuffer<uint> configuration;

[shader("compute")]
[numthreads(8,1,1)]
void ComputeMain(uint3 thread:SV_DispatchThreadID) {
    if (thread.x >= configuration[0]) return;
    uint offset = 10 * thread.x;
    uint output = 20 * thread.x;
    Real radius = inputs[offset], spin = inputs[offset + 1];
    Real aOverR = spin / radius;
    Real equatorial = radius * sqrt(Real(1) + aOverR * aOverR);
    Vec4Cart start = MakeVec4Cart(Real(0), inputs[offset + 2] * equatorial, Real(0), Real(0));
    Vec4Cart finish = MakeVec4Cart(Real(1), inputs[offset + 3] * equatorial, Real(0), Real(0));
    Vec4Cart velocity = AddVec4Cart(finish, ScaleVec4Cart(start, Real(-1)));
    Real fraction;
    Vec4Cart position, tangent;
    bool found = FindKerrCaptureEvent(start, velocity, finish, velocity, Real(1), radius,
                                     spin, fraction, position, tangent);
    outputs[output] = found ? Real(1) : Real(0);
    outputs[output + 1] = fraction;
    outputs[output + 2] = start.x;
    outputs[output + 3] = finish.x;
    outputs[output + 4] = equatorial;
    outputs[output + 5] = position.x;
}
