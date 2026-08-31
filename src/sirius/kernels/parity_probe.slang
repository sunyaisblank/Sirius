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
static const uint OP_DISK_TEMP = 2;
static const uint OP_BLACKBODY = 3;
static const uint OP_DEVIATION = 4;
static const uint OP_DISK_REDSHIFT = 5;
static const uint OP_LIVE_CART_CONSERVATION = 6;
static const uint OP_BEAM_ELLIPSE = 7;
static const uint OP_ADAPTIVE_EVENT_DOMAIN = 8;
static const uint OP_VOLUME_OPACITY = 9;
static const uint OP_NULL_PROJECTION = 10;
static const uint OP_CELESTIAL_TANGENT_BASIS = 11;
static const uint OP_JACOBI_RADIAL_CONGRUENCE = 12;

// Future-directed null root used by the render kernel, repeated here only for
// the diagnostic launch boundary. The integration below calls the production
// Cartesian RK4 authority directly.
Real ProbeNullTimeComponent(Real g[4][4], Real ux, Real uy, Real uz) {
    Real A = g[0][0];
    Real B = 2.0f * (g[0][1] * ux + g[0][2] * uy + g[0][3] * uz);
    Real C = g[1][1] * ux * ux + g[2][2] * uy * uy + g[3][3] * uz * uz
           + 2.0f * (g[1][2] * ux * uy + g[1][3] * ux * uz + g[2][3] * uy * uz);
    Real discriminant = B * B - 4.0f * A * C;
    if (discriminant >= 0.0f && abs(A) > Real(1.0e-20f)) {
        Real sqrtD = sqrt(discriminant);
        Real root1 = (-B + sqrtD) / (2.0f * A);
        Real root2 = (-B - sqrtD) / (2.0f * A);
        return (root1 > 0.0f) ? root1 : ((root2 > 0.0f) ? root2 : Real(0.0f));
    }
    if (abs(B) > Real(1.0e-20f)) {
        Real linearRoot = -C / B;
        return linearRoot > Real(0.0f) ? linearRoot : Real(0.0f);
    }
    return Real(0.0f);
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
    Real r = ComputeKerrRadiusCart(x, y, z, spin);
    Real r2 = r * r;
    Real discriminantRoot = 2.0f * r2 - reducedRadius2;

    Real drDx = x * r / discriminantRoot;
    Real drDy = y * r / discriminantRoot;
    Real drDz = z * (r2 + a2) / (r * discriminantRoot);

    Real cosTheta = clamp(z / r, Real(-1.0f), Real(1.0f));
    Real sin2Theta = Real(1.0f) - cosTheta * cosTheta;
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
                out Real g[4][4], out Real g_inv[4][4]) {
    if (metricId == kMetricEllis) {
        float sphericalG[4][4];
        float sphericalInverse[4][4];
        GetEllisDrainholeMetric(MakeVec4(0.0f, c0, c1, c2), p1, p2,
                                sphericalG, sphericalInverse);
        [unroll] for (int m = 0; m < 4; ++m)
            [unroll] for (int n = 0; n < 4; ++n) {
                g[m][n] = Real(sphericalG[m][n]);
                g_inv[m][n] = Real(sphericalInverse[m][n]);
            }
    } else if (metricId == kMetricWarpDrive) {
        // family params: p1=vs, p2=sigma, p3=R; bubble at origin.
        GetWarpDriveMetric(MakeVec4Cart(0.0f, c0, c1, c2), p1, p2, p3, 0.0f, 0.0f, 0.0f, g, g_inv);
    } else if (metricId == kMetricMorrisThorne) {
        // family params: p1=b0, p2=Phi0, p3=shapeType.
        float sphericalG[4][4];
        float sphericalInverse[4][4];
        GetMorrisThorneMetric(MakeVec4(0.0f, c0, c1, c2), p1, p2, int(p3),
                              sphericalG, sphericalInverse);
        [unroll] for (int m = 0; m < 4; ++m)
            [unroll] for (int n = 0; n < 4; ++n) {
                g[m][n] = Real(sphericalG[m][n]);
                g_inv[m][n] = Real(sphericalInverse[m][n]);
            }
    } else {
        GetKerrSchildFamilyMetric(MakeVec4Cart(0.0f, c0, c1, c2), p1, p2, p3, p4, g, g_inv);
    }
}

void ChristoffelById(uint metricId, float p1, float p2, float p3, float p4,
                     float c0, float c1, float c2,
                     out Real Gamma[4][4][4]) {
    if (metricId == kMetricWarpDrive) {
        GetWarpDriveChristoffel(MakeVec4Cart(0.0f, c0, c1, c2), p1, p2, p3, 0.0f, 0.0f, 0.0f, Gamma);
    } else if (metricId == kMetricMorrisThorne) {
        float sphericalGamma[4][4][4];
        GetMorrisThorneChristoffel(MakeVec4(0.0f, c0, c1, c2), p1, p2, int(p3),
                                   sphericalGamma);
        [unroll] for (int m = 0; m < 4; ++m)
            [unroll] for (int n = 0; n < 4; ++n)
                [unroll] for (int k = 0; k < 4; ++k)
                    Gamma[m][n][k] = Real(sphericalGamma[m][n][k]);
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
        Real g[4][4];
        Real g_inv[4][4];
        MetricById(metricId, p1, p2, p3, p4, c0, c1, c2, g, g_inv);
        [unroll] for (int m = 0; m < 4; m++)
            [unroll] for (int n = 0; n < 4; n++) {
                results[rb + uint(m * 4 + n)] = float(g[m][n]);
                results[rb + uint(16 + m * 4 + n)] = float(g_inv[m][n]);
            }
        // Inverse-metric identity error: max |g^{mk} g_{kn} - delta^m_n|.
        Real maxErr = Real(0.0f);
        for (int m = 0; m < 4; m++)
            for (int n = 0; n < 4; n++) {
                Real prod = Real(0.0f);
                for (int k = 0; k < 4; k++) prod += g[m][k] * g_inv[k][n];
                Real target = (m == n) ? Real(1.0f) : Real(0.0f);
                maxErr = max(maxErr, abs(prod - target));
            }
        results[rb + 32] = float(maxErr);
    } else if (opcode == OP_CHRISTOFFEL) {
        Real Gamma[4][4][4];
        ChristoffelById(metricId, p1, p2, p3, p4, c0, c1, c2, Gamma);
        [unroll] for (int m = 0; m < 4; m++)
            [unroll] for (int n = 0; n < 4; n++)
                [unroll] for (int k = 0; k < 4; k++)
                    results[rb + uint(m * 16 + n * 4 + k)] = float(Gamma[m][n][k]);
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
        results[rb + 3] = ShakuraSunyaevTemperature(innerT, r, r_isco);
    } else if (opcode == OP_BLACKBODY) {
        // c0 = temperature (K).
        float3 rgb = BlackbodyColor(c0);
        results[rb + 0] = rgb.x;
        results[rb + 1] = rgb.y;
        results[rb + 2] = rgb.z;
    } else if (opcode == OP_DISK_REDSHIFT) {
        // c0 = emitted bolometric source, c1 = redshift g.
        results[rb + 0] = ObservedDiskBolometricSource(c0, c1);
    } else if (opcode == OP_VOLUME_OPACITY) {
        // p1=reference radius, p2=H, p3=tau_ref; c0=radius, c1=z.
        results[rb + 0] =
            TruncatedGaussianOpacityDensity(p3, c0, p1, p2, c1);
    } else if (opcode == OP_NULL_PROJECTION) {
        // Direct witness for the production accepted-state constraint method.
        // The sample supplies a metric family/event and u=(t,x,y,z).
        Real metric[4][4];
        Real metricInverse[4][4];
        MetricById(metricId, p1, p2, p3, p4, c0, c1, c2,
                   metric, metricInverse);
        Vec4Cart tangent = MakeVec4Cart(Real(samples[sb + 8]), Real(samples[sb + 9]),
                                        Real(samples[sb + 10]), Real(samples[sb + 11]));
        Vec4Cart projected;
        bool represented =
            ProjectNullTangentPreservingBranch(metric, tangent, projected);
        results[rb + 0] = represented ? 1.0f : 0.0f;
        results[rb + 1] = float(projected.t);
        results[rb + 2] = float(projected.x);
        results[rb + 3] = float(projected.y);
        results[rb + 4] = float(projected.z);
        results[rb + 5] = float(RelativeNullConstraint(metric, tangent));
        results[rb + 6] = represented
                              ? float(RelativeNullConstraint(metric, projected))
                              : 1.0e30f;
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
        if (state.u.t <= Real(0.0f)) {
            [unroll] for (int component = 0; component < 9; ++component) {
                results[rb + uint(component)] = 0.0f;
            }
            results[rb + 9] = -1.0f;
            return;
        }

        Real energy0, angularMomentum0, carter0, null0;
        KerrCartConservedQuantities(state, p1, p2, energy0, angularMomentum0, carter0, null0);
        Real worstEnergy = 0.0f;
        Real worstAngularMomentum = 0.0f;
        Real worstCarter = 0.0f;
        Real worstNull = RelativeNullConstraint(launchMetric, state.u);
        int completedSteps = 0;
        int maxSteps = int(samples[sb + 13]);
        Real stepScale = Real(samples[sb + 12]);
        Real minStep = Real(samples[sb + 14]);
        Real maxStep = Real(samples[sb + 15]);
        Real adaptiveStep = maxStep;
        Real horizon = Real(p1 + sqrt(max(p1 * p1 - p2 * p2, 0.0f)));
        const Real escapeRadius = Real(200.0f);
        int termination = 0;
        Real lastRejectedError = Real(0.0f);
        Real lastRejectedNull = Real(0.0f);
#ifdef SIRIUS_FP32_COMP
        Vec4Cart compX = MakeVec4Cart(0.0f, 0.0f, 0.0f, 0.0f);
        Vec4Cart compU = MakeVec4Cart(0.0f, 0.0f, 0.0f, 0.0f);
#endif

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
            Real geometricStep = clamp(stepScale * rKerr, minStep, maxStep);
            Real dt = min(geometricStep, adaptiveStep);
            AdaptiveGeodesicStep attempted = AttemptAdaptiveGeodesicStep(
                state, dt, false, false, p1, p2, 0.0f, 0.0f,
                0.0f, 0.0f, 0.0f, 0.0f, 0,
                AdaptiveIntegrationTolerance(), AdaptiveNullTolerance());
            if (!attempted.accepted) {
                lastRejectedError = attempted.errorRatio;
                lastRejectedNull = attempted.nullResidual;
                if (dt <= minStep * Real(1.0001f)) {
                    termination = 3;
                    break;
                }
                adaptiveStep = max(minStep, dt * attempted.stepFactor);
                continue;
            }
            adaptiveStep = clamp(dt * attempted.stepFactor, minStep, maxStep);
            GeodesicStateCart accepted = attempted.candidate;
#ifdef SIRIUS_FP32_COMP
            Vec4Cart acceptedDX = AddVec4Cart(
                accepted.x, ScaleVec4Cart(state.x, Real(-1.0f)));
            Vec4Cart acceptedDU = AddVec4Cart(
                accepted.u, ScaleVec4Cart(state.u, Real(-1.0f)));
            accepted.x = KahanAddVec4Cart(state.x, compX, acceptedDX);
            accepted.u = KahanAddVec4Cart(state.u, compU, acceptedDU);
#endif
            state = accepted;
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
            worstNull = max(worstNull, attempted.nullResidual);
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
        results[rb + 10] = float(lastRejectedError);
        results[rb + 11] = float(lastRejectedNull);
    } else if (opcode == OP_ADAPTIVE_EVENT_DOMAIN) {
        // Direct event-domain witness for the production adaptive RK4
        // authority. c=(x,y,z), u=(t,x,y,z), h=attempted step.
        GeodesicStateCart state;
        state.x = MakeVec4Cart(Real(0.0f), Real(c0), Real(c1), Real(c2));
        state.u = MakeVec4Cart(Real(samples[sb + 8]), Real(samples[sb + 9]),
                               Real(samples[sb + 10]), Real(samples[sb + 11]));
        Vec4Cart secondStage = AddVec4Cart(state.x, ScaleVec4Cart(state.u, Real(0.5f * samples[sb + 12])));
        AdaptiveGeodesicStep attempted = AttemptAdaptiveGeodesicStep(
            state, Real(samples[sb + 12]), false, false, p1, p2, p3, p4,
            0.0f, 0.0f, 0.0f, 0.0f, 0, Real(1.0f), Real(1.0f));
        results[rb + 0] = IsKerrSchildFamilyEventRepresented(state.x, p1, p2, p3, p4) ? 1.0f : 0.0f;
        results[rb + 1] = IsKerrSchildFamilyEventRepresented(secondStage, p1, p2, p3, p4) ? 1.0f : 0.0f;
        results[rb + 2] = attempted.accepted ? 1.0f : 0.0f;
        results[rb + 3] = float(attempted.candidate.x.x);
        results[rb + 4] = float(attempted.candidate.x.y);
        results[rb + 5] = float(attempted.candidate.x.z);
        results[rb + 6] = float(attempted.errorRatio);
        results[rb + 7] = float(attempted.stepFactor);
    } else if (opcode == OP_BEAM_ELLIPSE) {
        // c0,c1,c2,u0 are the row-major transverse Jacobi map.
        float semiMajor, semiMinor, orientation, determinant;
        BeamEllipseFromTransverseMap(c0, c1, c2, samples[sb + 8], semiMajor, semiMinor,
                                     orientation, determinant);
        results[rb + 0] = semiMajor;
        results[rb + 1] = semiMinor;
        results[rb + 2] = orientation;
        results[rb + 3] = determinant;
    } else if (opcode == OP_CELESTIAL_TANGENT_BASIS) {
        // c0,c1,c2 are a finite nonzero celestial direction. This is the same
        // helper called by the terminal Sachs screen and point-star filter.
        float3 first;
        float3 second;
        MakeCelestialTangentBasis(float3(c0, c1, c2), first, second);
        results[rb + 0] = first.x;
        results[rb + 1] = first.y;
        results[rb + 2] = first.z;
        results[rb + 3] = second.x;
        results[rb + 4] = second.y;
        results[rb + 5] = second.z;
    } else if (opcode == OP_JACOBI_RADIAL_CONGRUENCE) {
        // Integrate a two-axis point-source bundle along the exact affinely
        // parameterised ingoing Schwarzschild principal null ray
        // (T,x)=(lambda,r0-lambda), k=(1,-1,0,0). The analytic Sachs axes are
        // both seed*lambda. p1=M; c0=start radius, c1=end radius, c2=max step;
        // h=seed. This is scalar/tiny-compute evidence and creates no image.
        Real startRadius = Real(c0);
        Real endRadius = Real(c1);
        Real maxSegment = Real(c2);
        Real seed = Real(samples[sb + 12]);
        Real4 xiA = Real4(0.0f);
        Real4 xiB = Real4(0.0f);
        Real4 VA = Real4(Real(0.0), Real(0.0), seed, Real(0.0));
        Real4 VB = Real4(Real(0.0), Real(0.0), Real(0.0), seed);
        Real radius = startRadius;
        Real affine = Real(0.0);
        int steps = 0;
        if (p1 > 0.0f && endRadius > Real(2.0f * p1) &&
            startRadius > endRadius && maxSegment > Real(0.0)) {
            for (int segment = 0; segment < 128 && radius > endRadius; ++segment) {
                Real dLambda = min(maxSegment, radius - endRadius);
                Vec4Cart startPosition = MakeVec4Cart(affine, radius, Real(0.0), Real(0.0));
                Vec4Cart endPosition = MakeVec4Cart(
                    affine + dLambda, radius - dLambda, Real(0.0), Real(0.0));
                Vec4Cart tangent = MakeVec4Cart(Real(1.0), Real(-1.0), Real(0.0), Real(0.0));
                AdvanceJacobiPairAcceptedSegment(
                    startPosition, tangent, endPosition, tangent, dLambda,
                    p1, 0.0f, 0.0f, 0.0f, xiA, VA, xiB, VB);
                affine += dLambda;
                radius -= dLambda;
                steps += 1;
            }
        }
        results[rb + 0] = float(xiA.z);
        results[rb + 1] = float(xiB.w);
        results[rb + 2] = float(xiA.x);
        results[rb + 3] = float(xiA.y);
        results[rb + 4] = float(xiB.x);
        results[rb + 5] = float(xiB.y);
        results[rb + 6] = float(VA.z);
        results[rb + 7] = float(VB.w);
        results[rb + 8] = float(affine);
        results[rb + 9] = float(steps);
    } else if (opcode == OP_DEVIATION) {
        // Full-Riemann ray-bundle deviation at a Cartesian point.
        // p1=M, p2=a, p3=Q, p4=Lambda; c=(x,y,z);
        // u=(k^t,k^x,k^y,k^z); h,aux0..2=(xi^t,xi^x,xi^y,xi^z).
        Vec4Cart pos = MakeVec4Cart(0.0f, c0, c1, c2);
        Vec4Cart k = MakeVec4Cart(samples[sb + 8], samples[sb + 9], samples[sb + 10], samples[sb + 11]);
        Real4 xi = Real4(Real(samples[sb + 12]), Real(samples[sb + 13]),
                          Real(samples[sb + 14]), Real(samples[sb + 15]));
        Real4 full = GeodesicDeviationAccelFullRiemann(pos, k, xi, p1, p2, p3, p4);
        results[rb + 0] = float(full.x); results[rb + 1] = float(full.y);
        results[rb + 2] = float(full.z); results[rb + 3] = float(full.w);
        results[rb + 4] = 0.0f; results[rb + 5] = 0.0f;
        results[rb + 6] = 0.0f; results[rb + 7] = 0.0f;
        results[rb + 8] = length(float3(float(full.y), float(full.z), float(full.w)));
    }
}
