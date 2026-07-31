// Parity probe: evaluates the ported device physics at host-supplied sample
// points and writes the results back for comparison against reference values
// computed from the legacy RDOP002A.cu functions (tests/backend/
// kernel_parity_test.cpp). One kernel, several opcodes, so the host drives
// every quantity through the same buffer ABI.
//
// Buffer ABI (storage buffers, bound 0..2 in declaration order):
//   0 samples : flat float array, sampleStride floats per sample.
//   1 results : flat float array, resultStride floats per sample.
//   2 params  : [0] opcode, [1] sampleCount, [2] sampleStride, [3] resultStride.
//
// Sample record (sampleStride = 16):
//   [0] metricId   [1..4] family params (Kerr-Schild M,a,Q,Lambda; Ellis m,n;
//                          warp vs,sigma,R)
//   [5..7] position coords (Cartesian x,y,z or spherical r,theta,phi)
//   [8..11] 4-velocity u (symplectic step)
//   [12] step h   [13] aux0   [14] aux1   [15] aux2
import gr_types;
import gr_metrics;
import gr_integrator;
import gr_disk;
import gr_deviation;

StructuredBuffer<float> samples;     // binding 0
RWStructuredBuffer<float> results;   // binding 1
StructuredBuffer<float> params;      // binding 2

static const uint OP_METRIC = 0;
static const uint OP_CHRISTOFFEL = 1;
static const uint OP_SYMPLECTIC = 2;
static const uint OP_DISK_TEMP = 3;
static const uint OP_BLACKBODY = 4;
static const uint OP_SYMPLECTIC_S4 = 5;
static const uint OP_DEVIATION = 6;

void MetricById(uint metricId, float p1, float p2, float p3, float p4,
                float c0, float c1, float c2,
                out float g[4][4], out float g_inv[4][4]) {
    if (metricId == kMetricEllis) {
        GetEllisDrainholeMetric(MakeVec4(0.0f, c0, c1, c2), p1, p2, g, g_inv);
    } else if (metricId == kMetricWarpDrive) {
        // family params: p1=vs, p2=sigma, p3=R; bubble at origin.
        GetWarpDriveMetric(MakeVec4Cart(0.0f, c0, c1, c2), p1, p2, p3, 0.0f, 0.0f, 0.0f, g, g_inv);
    } else if (metricId == kMetricMorrisThorne) {
        // family params: p1=b0, p2=Phi0, p3=shapeType.
        GetMorrisThorneMetric(MakeVec4(0.0f, c0, c1, c2), p1, p2, int(p3), g, g_inv);
    } else {
        GetKerrSchildFamilyMetric(MakeVec4Cart(0.0f, c0, c1, c2), p1, p2, p3, p4, g, g_inv);
    }
}

void ChristoffelById(uint metricId, float p1, float p2, float p3, float p4,
                     float c0, float c1, float c2,
                     out float Gamma[4][4][4]) {
    if (metricId == kMetricWarpDrive) {
        GetWarpDriveChristoffel(MakeVec4Cart(0.0f, c0, c1, c2), p1, p2, p3, 0.0f, 0.0f, 0.0f, Gamma);
    } else if (metricId == kMetricMorrisThorne) {
        GetMorrisThorneChristoffel(MakeVec4(0.0f, c0, c1, c2), p1, p2, int(p3), Gamma);
    } else {
        GetKerrSchildFamilyChristoffel(MakeVec4Cart(0.0f, c0, c1, c2), p1, p2, p3, p4, Gamma);
    }
}

[shader("compute")]
[numthreads(64, 1, 1)]
void ComputeMain(uint3 tid: SV_DispatchThreadID) {
    uint opcode = uint(params[0]);
    uint count = uint(params[1]);
    uint sStride = uint(params[2]);
    uint rStride = uint(params[3]);
    uint i = tid.x;
    if (i >= count) return;

    uint sb = i * sStride;
    uint rb = i * rStride;

    uint metricId = uint(samples[sb + 0]);
    float p1 = samples[sb + 1];
    float p2 = samples[sb + 2];
    float p3 = samples[sb + 3];
    float p4 = samples[sb + 4];
    float c0 = samples[sb + 5];
    float c1 = samples[sb + 6];
    float c2 = samples[sb + 7];

    if (opcode == OP_METRIC) {
        float g[4][4];
        float g_inv[4][4];
        MetricById(metricId, p1, p2, p3, p4, c0, c1, c2, g, g_inv);
        [unroll] for (int m = 0; m < 4; m++)
            [unroll] for (int n = 0; n < 4; n++) {
                results[rb + uint(m * 4 + n)] = g[m][n];
                results[rb + uint(16 + m * 4 + n)] = g_inv[m][n];
            }
        // Inverse-metric identity error: max |g^{mk} g_{kn} - delta^m_n|.
        float maxErr = 0.0f;
        for (int m = 0; m < 4; m++)
            for (int n = 0; n < 4; n++) {
                float prod = 0.0f;
                for (int k = 0; k < 4; k++) prod += g[m][k] * g_inv[k][n];
                float target = (m == n) ? 1.0f : 0.0f;
                maxErr = max(maxErr, abs(prod - target));
            }
        results[rb + 32] = maxErr;
    } else if (opcode == OP_CHRISTOFFEL) {
        float Gamma[4][4][4];
        ChristoffelById(metricId, p1, p2, p3, p4, c0, c1, c2, Gamma);
        [unroll] for (int m = 0; m < 4; m++)
            [unroll] for (int n = 0; n < 4; n++)
                [unroll] for (int k = 0; k < 4; k++)
                    results[rb + uint(m * 16 + n * 4 + k)] = Gamma[m][n][k];
    } else if (opcode == OP_SYMPLECTIC) {
        // Kerr Boyer-Lindquist symplectic step: p1=M, p2=a.
        GeodesicState st;
        st.x = MakeVec4(0.0f, c0, c1, c2);
        st.u = MakeVec4(samples[sb + 8], samples[sb + 9], samples[sb + 10], samples[sb + 11]);
        float h = samples[sb + 12];
        float tolerance = samples[sb + 13];
        float minStep = samples[sb + 14];
        float maxStep = samples[sb + 15];
        bool accepted;
        GeodesicState res = IntegrateGeodesicSymplectic(st, h, p1, p2, tolerance, minStep, maxStep, accepted);
        results[rb + 0] = res.x.t;
        results[rb + 1] = res.x.r;
        results[rb + 2] = res.x.theta;
        results[rb + 3] = res.x.phi;
        results[rb + 4] = res.u.t;
        results[rb + 5] = res.u.r;
        results[rb + 6] = res.u.theta;
        results[rb + 7] = res.u.phi;
        results[rb + 8] = h;
        results[rb + 9] = accepted ? 1.0f : 0.0f;
    } else if (opcode == OP_SYMPLECTIC_S4) {
        // One Yoshida 4th-order step on the Hamiltonian state built from the
        // input geodesic: p1=M, p2=a; h=sample[12].
        GeodesicState st;
        st.x = MakeVec4(0.0f, c0, c1, c2);
        st.u = MakeVec4(samples[sb + 8], samples[sb + 9], samples[sb + 10], samples[sb + 11]);
        float h = samples[sb + 12];
        HamiltonianState hs = GeodesicToHamiltonian(st, p1, p2);
        IntegrateSymplecticS4(hs, h, p1, p2);
        [unroll] for (int k = 0; k < 4; k++) results[rb + uint(k)] = hs.q[k];
        [unroll] for (int k = 0; k < 4; k++) results[rb + uint(4 + k)] = hs.p[k];
    } else if (opcode == OP_DISK_TEMP) {
        // p1=M, p2=a_star; c0=r; aux0=r_isco, aux1=innerTemperature.
        float r = c0;
        float r_isco = samples[sb + 13];
        float innerT = samples[sb + 14];
        float fluxShape = FullPageThorneFluxShape(r / p1, r_isco / p1, p2);
        float T = FullPageThorneTemperature(innerT, r, r_isco, p1, p2);
        results[rb + 0] = T;
        results[rb + 1] = fluxShape;
        results[rb + 2] = ComputeKerrISCO(p1, p2);
    } else if (opcode == OP_BLACKBODY) {
        // c0 = temperature (K).
        float3 rgb = BlackbodyColor(c0);
        results[rb + 0] = rgb.x;
        results[rb + 1] = rgb.y;
        results[rb + 2] = rgb.z;
    } else if (opcode == OP_DEVIATION) {
        // Ray-bundle deviation at a Cartesian point: full-Riemann and tidal
        // acceleration on a unit deviation vector. p1=M, p2=a, p3=Q, p4=Lambda;
        // c=(x,y,z); u=(k^t,k^x,k^y,k^z); aux0=step.
        Vec4Cart pos = MakeVec4Cart(0.0f, c0, c1, c2);
        Vec4Cart k = MakeVec4Cart(samples[sb + 8], samples[sb + 9], samples[sb + 10], samples[sb + 11]);
        float4 xi = float4(0.0f, 1.0f, 0.0f, 0.0f);
        float step = samples[sb + 13];
        float rKerr = ComputeKerrRadiusCart(c0, c1, c2, p2);
        float4 full = GeodesicDeviationAccelFullRiemann(pos, k, xi, p1, p2, p3, p4);
        float4 tidal = TidalDeviationAccel(xi, p1, rKerr, step);
        results[rb + 0] = full.x; results[rb + 1] = full.y;
        results[rb + 2] = full.z; results[rb + 3] = full.w;
        results[rb + 4] = tidal.x; results[rb + 5] = tidal.y;
        results[rb + 6] = tidal.z; results[rb + 7] = tidal.w;
        results[rb + 8] = length(float3(full.y, full.z, full.w));
    }
}
