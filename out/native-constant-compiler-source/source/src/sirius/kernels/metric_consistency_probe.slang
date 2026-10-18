// Pointwise production metric/connection observation, with Real-width I/O.
// One invocation owns one event; no integration or radiance publication.
import gr_types;
import gr_metrics;

[[vk::binding(0, 0)]] StructuredBuffer<Real> inputs;
[[vk::binding(1, 0)]] RWStructuredBuffer<Real> outputs;

[shader("compute")]
[numthreads(1, 1, 1)]
void ComputeMain(uint3 id : SV_DispatchThreadID) {
    if (id.x != 0 || id.y != 0 || id.z != 0) return;
    outputs[0] = Real(0);
    float M = float(inputs[0]), a = float(inputs[1]);
    float Q = float(inputs[2]), Lambda = float(inputs[3]);
    Vec4Cart position = MakeVec4Cart(inputs[4], inputs[5], inputs[6], inputs[7]);
    if (!IsKerrSchildFamilyEventRepresented(position, M, a, Q, Lambda)) return;
    Real metric[4][4], inverse[4][4], connection[4][4][4];
    GetKerrSchildFamilyMetric(position, M, a, Q, Lambda, metric, inverse);
    GetKerrSchildFamilyChristoffel(position, M, a, Q, Lambda, connection);
    Real4 columns[4];
    for (int i = 0; i < 4; ++i) {
        columns[i] = Real4(Real(0));
        columns[i][i] = Real(1);
    }
    CoupledJet4 jetMetric[4][4], jetConnection[4][4][4];
    if (!TryKerrSchildConnectionJet(position, columns, M, a, Q, Lambda,
                                   jetMetric, jetConnection)) return;
    Real metricDifference = Real(0), connectionDifference = Real(0);
    for (int mu = 0; mu < 4; ++mu) for (int nu = 0; nu < 4; ++nu) {
        outputs[1 + 4*mu + nu] = metric[mu][nu];
        outputs[17 + 4*mu + nu] = inverse[mu][nu];
        metricDifference = max(metricDifference, abs(metric[mu][nu] - jetMetric[mu][nu].value));
        for (int column = 0; column < 4; ++column)
            outputs[97 + 16*column + 4*mu + nu] = jetMetric[mu][nu].derivative[column];
        for (int rho = 0; rho < 4; ++rho) {
            outputs[33 + 16*mu + 4*nu + rho] = connection[mu][nu][rho];
            connectionDifference = max(connectionDifference,
                abs(connection[mu][nu][rho] - jetConnection[mu][nu][rho].value));
            for (int column = 0; column < 4; ++column)
                outputs[161 + 64*column + 16*mu + 4*nu + rho] =
                    jetConnection[mu][nu][rho].derivative[column];
        }
    }
    outputs[417] = metricDifference;
    outputs[418] = connectionDifference;
    outputs[419] = Real(0);
    outputs[0] = Real(1);
}
