import gr_types;
import gr_metrics;
import gr_trace_event;

StructuredBuffer<Real> inputs;
RWStructuredBuffer<Real> outputs;
StructuredBuffer<uint> configuration;

void Store(uint base, Vec4Cart vector) {
    outputs[base] = vector.t;
    outputs[base + 1] = vector.x;
    outputs[base + 2] = vector.y;
    outputs[base + 3] = vector.z;
}

[shader("compute")]
[numthreads(8,1,1)]
void ComputeMain(uint3 thread:SV_DispatchThreadID) {
    if (thread.x >= configuration[0]) return;
    uint at = 10 * thread.x;
    uint outAt = 20 * thread.x;
    Real mass = inputs[at], spin = inputs[at + 1];
    Vec4Cart x = MakeVec4Cart(inputs[at+2],inputs[at+3],inputs[at+4],inputs[at+5]);
    Vec4Cart u = MakeVec4Cart(inputs[at+6],inputs[at+7],inputs[at+8],inputs[at+9]);
    KerrPastHorizonMap map, reverse;
    if (!BuildKerrPastHorizonMap(x, mass, spin, map)) return;
    Vec4Cart mapped = MapKerrPastHorizonVector(map,u);
    outputs[outAt] = Real(1);
    Store(outAt+1,map.position);
    Store(outAt+5,mapped);
    if (!BuildKerrPastHorizonMap(map.position,mass,spin,reverse)) return;
    outputs[outAt+9] = Real(1);
    Store(outAt+10,reverse.position);
    Store(outAt+14,MapKerrPastHorizonVector(reverse,mapped));

    // Noncentral curved Hermite segment with coincident exterior endpoints
    // and two contacts, on either side of a central interior interval.
    Real radius = mass * Real(2);
    Real spinRatio = spin / radius;
    Real equatorial = radius * sqrt(Real(1) + spinRatio*spinRatio);
    Vec4Cart start = MakeVec4Cart(Real(6),Real(1.25f)*equatorial,Real(0),Real(0.5f)*radius);
    Vec4Cart v0 = MakeVec4Cart(Real(2),-Real(3)*equatorial,Real(0),Real(0));
    Vec4Cart v1 = MakeVec4Cart(Real(2),Real(3)*equatorial,Real(0),Real(0));
    Vec4Cart hit, tangent;
    Real fraction;
    if (FindKerrCaptureEvent(start,v0,start,v1,Real(1),radius,spin,fraction,hit,tangent)) {
        outputs[outAt+18] = Real(1);
        outputs[outAt+19] = fraction;
    }
}
