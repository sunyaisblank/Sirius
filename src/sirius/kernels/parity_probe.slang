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
static const uint OP_DISK_REDSHIFT = 7;
static const uint OP_LIVE_CART_CONSERVATION = 8;
static const uint OP_BEAM_ELLIPSE = 9;

// Future-directed null root used by the render kernel, repeated here only for
// the diagnostic launch boundary. The integration below calls the production
// Cartesian RK4 authority directly.
Real ProbeNullTimeComponent(Real g[4][4], Real ux, Real uy, Real uz) {
    Real A = g[0][0];
    Real B = 2.0f * (g[0][1] * ux + g[0][2] * uy + g[0][3] * uz);
    Real C = g[1][1] * ux * ux + g[2][2] * uy * uy + g[3][3] * uz * uz
           + 2.0f * (g[1][2] * ux * uy + g[1][3] * ux * uz + g[2][3] * uy * uz);
    Real discriminant = B * B - 4.0f * A * C;
    if (discriminant >= 0.0f) {
        Real sqrtD = sqrt(discriminant);
        Real root1 = (-B + sqrtD) / (2.0f * A);
        Real root2 = (-B - sqrtD) / (2.0f * A);
        return (root1 > 0.0f) ? root1 : ((root2 > 0.0f) ? root2 : Real(1.0f));
    }
    return sqrt(max(-C / A, Real(0.0001f)));
}

// Stationary-axisymmetric Kerr invariants for a Cartesian Kerr-Schild state.
// E and L_z are evaluated directly with the two Killing vectors. Carter Q
// requires p_theta, obtained from the exact inverse oblate spatial Jacobian;
// this is the Slang counterpart of coordinates::
// TransformVectorKerrSchildCartToBl used by the CPU live-path monitor.
void KerrCartConservedQuantities(GeodesicStateCart state, float M, float a,
                                 out Real energy, out Real angularMomentum,
                                 out Real carter, out Real nullResidual) {
    Real g[4][4];
    Real gInv[4][4];
    GetKerrSchildFamilyMetric(state.x, M, a, 0.0f, 0.0f, g, gInv);
    Real u[4] = { state.u.t, state.u.x, state.u.y, state.u.z };

    energy = 0.0f;
    [unroll] for (int nu = 0; nu < 4; nu++) energy -= g[0][nu] * u[nu];

    Real xi[4] = { 0.0f, -state.x.y, state.x.x, 0.0f };
    angularMomentum = 0.0f;
    nullResidual = 0.0f;
    [unroll] for (int mu = 0; mu < 4; mu++) {
        [unroll] for (int nu = 0; nu < 4; nu++) {
            angularMomentum += xi[mu] * g[mu][nu] * u[nu];
            nullResidual += g[mu][nu] * u[mu] * u[nu];
        }
    }

    Real x = state.x.x;
    Real y = state.x.y;
    Real z = state.x.z;
    Real spin = Real(a);
    Real a2 = spin * spin;
    Real radius2 = x * x + y * y + z * z;
    Real reducedRadius2 = radius2 - a2;
    Real discriminant = max(reducedRadius2 * reducedRadius2 + 4.0f * a2 * z * z,
                            Real(1.0e-20f));
    Real discriminantRoot = sqrt(discriminant);
    Real r2 = max((reducedRadius2 + discriminantRoot) / 2.0f, Real(1.0e-12f));
    Real r = sqrt(r2);

    Real discriminantDx = 4.0f * x * reducedRadius2;
    Real discriminantDy = 4.0f * y * reducedRadius2;
    Real discriminantDz = 4.0f * z * (radius2 + a2);
    Real drDx = (x + discriminantDx / (4.0f * discriminantRoot)) / (2.0f * r);
    Real drDy = (y + discriminantDy / (4.0f * discriminantRoot)) / (2.0f * r);
    Real drDz = (z + discriminantDz / (4.0f * discriminantRoot)) / (2.0f * r);

    Real cosTheta = clamp(z / r, Real(-1.0f), Real(1.0f));
    Real sin2Theta = max(Real(1.0f) - cosTheta * cosTheta, Real(1.0e-20f));
    Real sinTheta = sqrt(sin2Theta);
    Real dThetaDx = z * drDx / (r2 * sinTheta);
    Real dThetaDy = z * drDy / (r2 * sinTheta);
    Real dThetaDz = -(r - z * drDz) / (r2 * sinTheta);
    Real uTheta = dThetaDx * state.u.x + dThetaDy * state.u.y + dThetaDz * state.u.z;
    Real sigma = r2 + a2 * cosTheta * cosTheta;
    Real pTheta = sigma * uTheta;
    carter = pTheta * pTheta + cosTheta * cosTheta *
             (-a2 * energy * energy + angularMomentum * angularMomentum / sin2Theta);
}

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
    } else if (opcode == OP_DISK_REDSHIFT) {
        // c0 = emitted bolometric source, c1 = redshift g.
        results[rb + 0] = ObservedDiskBolometricSource(c0, c1);
    } else if (opcode == OP_LIVE_CART_CONSERVATION) {
        // Actual Vulkan render integrator burn-in. p1=M, p2=a; c=(x,y,z);
        // u1..u3 are the spatial launch tangent; h=stepScale, aux0=maxSteps,
        // aux1=minStep, aux2=maxStep. Termination matches the production Kerr
        // trace (exact outer horizon and the governed 200M escape surface).
        GeodesicStateCart state;
        state.x = MakeVec4Cart(0.0f, c0, c1, c2);
        state.u = MakeVec4Cart(0.0f, samples[sb + 9], samples[sb + 10], samples[sb + 11]);
        Real launchMetric[4][4];
        Real launchMetricInv[4][4];
        GetKerrSchildFamilyMetric(state.x, p1, p2, 0.0f, 0.0f,
                                  launchMetric, launchMetricInv);
        state.u.t = ProbeNullTimeComponent(launchMetric, state.u.x, state.u.y, state.u.z);

        Real energy0, angularMomentum0, carter0, null0;
        KerrCartConservedQuantities(state, p1, p2, energy0, angularMomentum0, carter0, null0);
        Real worstEnergy = 0.0f;
        Real worstAngularMomentum = 0.0f;
        Real worstCarter = 0.0f;
        Real worstNull = abs(null0);
        int completedSteps = 0;
        int maxSteps = int(samples[sb + 13]);
        Real stepScale = Real(samples[sb + 12]);
        Real minStep = Real(samples[sb + 14]);
        Real maxStep = Real(samples[sb + 15]);
        Real horizon = Real(p1 + sqrt(max(p1 * p1 - p2 * p2, 0.0f)));
        const Real escapeRadius = Real(200.0f);
        int termination = 0;

        for (int step = 0; step < maxSteps; step++) {
            Real rKerr = ComputeKerrRadiusCart(state.x.x, state.x.y, state.x.z, Real(p2));
            Real rCartesian = sqrt(state.x.x * state.x.x + state.x.y * state.x.y +
                                   state.x.z * state.x.z);
            if (rKerr <= horizon) {
                termination = 1;
                break;
            }
            if (rCartesian >= escapeRadius) {
                termination = 2;
                break;
            }
            Real dt = clamp(stepScale * rKerr, minStep, maxStep);
            state = IntegrateGeodesicRK4Cart(state, dt, p1, p2, 0.0f, 0.0f);
            completedSteps++;

            Real energy, angularMomentum, carter, nullResidual;
            KerrCartConservedQuantities(state, p1, p2, energy, angularMomentum, carter, nullResidual);
            worstEnergy = max(worstEnergy,
                              abs(energy - energy0) / max(abs(energy0), Real(1.0e-12f)));
            worstAngularMomentum =
                max(worstAngularMomentum,
                    abs(angularMomentum - angularMomentum0) /
                    max(abs(angularMomentum0), Real(1.0e-12f)));
            worstCarter = max(worstCarter,
                              abs(carter - carter0) / max(abs(carter0), Real(1.0e-12f)));
            worstNull = max(worstNull, abs(nullResidual));
        }

        results[rb + 0] = float(worstEnergy);
        results[rb + 1] = float(worstAngularMomentum);
        results[rb + 2] = float(worstCarter);
        results[rb + 3] = float(worstNull);
        results[rb + 4] = float(completedSteps);
        results[rb + 5] = float(sqrt(state.x.x * state.x.x + state.x.y * state.x.y +
                                     state.x.z * state.x.z));
        results[rb + 6] = float(energy0);
        results[rb + 7] = float(angularMomentum0);
        results[rb + 8] = float(carter0);
        results[rb + 9] = float(termination);
    } else if (opcode == OP_BEAM_ELLIPSE) {
        // c0,c1,c2,u0 are the row-major transverse Jacobi map.
        float semiMajor, semiMinor, orientation, determinant;
        BeamEllipseFromTransverseMap(c0, c1, c2, samples[sb + 8], semiMajor, semiMinor,
                                     orientation, determinant);
        results[rb + 0] = semiMajor;
        results[rb + 1] = semiMinor;
        results[rb + 2] = orientation;
        results[rb + 3] = determinant;
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
