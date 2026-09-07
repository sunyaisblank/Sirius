// Cartesian RK4 geodesic increments for the represented live metric families.
// trace.slang owns adaptive step doubling, null-residual acceptance, retry, and
// precision-rung accumulation around these local tableaux. The independent
// canonical symplectic method belongs exclusively to the double-precision C++
// oracle; no disconnected device candidate is retained for legacy parity.
module gr_integrator;

import gr_types;
import gr_metrics;

// du^mu/dlambda = -Gamma^mu_nu_rho u^nu u^rho, Kerr-Schild Cartesian
// Christoffels. (M, a, Q, Lambda) carry the represented family parameters.
public bool TryGeodesicAccelerationCart(GeodesicStateCart state, float M, float a,
                                        float Q, float Lambda, out Vec4Cart acceleration) {
    acceleration = MakeVec4Cart(Real(0.0f), Real(0.0f), Real(0.0f), Real(0.0f));
    if (!IsKerrSchildFamilyEventRepresented(state.x, M, a, Q, Lambda)) return false;
    Real Gamma[4][4][4];
    GetKerrSchildFamilyChristoffel(state.x, M, a, Q, Lambda, Gamma);
    Real u[4] = { state.u.t, state.u.x, state.u.y, state.u.z };
    Real acc[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
    [unroll] for (int mu = 0; mu < 4; mu++)
        [unroll] for (int nu = 0; nu < 4; nu++)
            [unroll] for (int rho = 0; rho < 4; rho++)
                acc[mu] -= Gamma[mu][nu][rho] * u[nu] * u[rho];
    acceleration = MakeVec4Cart(acc[0], acc[1], acc[2], acc[3]);
    return true;
}

// Increment form leaves acceptance and precision-rung commit to trace.slang.
public bool IntegrateGeodesicRK4CartDelta(GeodesicStateCart state, Real dt,
                                          float M, float a, float Q, float Lambda,
                                          out Vec4Cart dx, out Vec4Cart du) {
    dx = MakeVec4Cart(Real(0.0f), Real(0.0f), Real(0.0f), Real(0.0f));
    du = MakeVec4Cart(Real(0.0f), Real(0.0f), Real(0.0f), Real(0.0f));
    Vec4Cart acceleration;
    if (!TryGeodesicAccelerationCart(state, M, a, Q, Lambda, acceleration)) return false;
    Vec4Cart dx1 = ScaleVec4Cart(state.u, dt);
    Vec4Cart du1 = ScaleVec4Cart(acceleration, dt);

    GeodesicStateCart s2;
    s2.x = AddVec4Cart(state.x, ScaleVec4Cart(dx1, 0.5f));
    s2.u = AddVec4Cart(state.u, ScaleVec4Cart(du1, 0.5f));
    if (!TryGeodesicAccelerationCart(s2, M, a, Q, Lambda, acceleration)) return false;
    Vec4Cart dx2 = ScaleVec4Cart(s2.u, dt);
    Vec4Cart du2 = ScaleVec4Cart(acceleration, dt);

    GeodesicStateCart s3;
    s3.x = AddVec4Cart(state.x, ScaleVec4Cart(dx2, 0.5f));
    s3.u = AddVec4Cart(state.u, ScaleVec4Cart(du2, 0.5f));
    if (!TryGeodesicAccelerationCart(s3, M, a, Q, Lambda, acceleration)) return false;
    Vec4Cart dx3 = ScaleVec4Cart(s3.u, dt);
    Vec4Cart du3 = ScaleVec4Cart(acceleration, dt);

    GeodesicStateCart s4;
    s4.x = AddVec4Cart(state.x, dx3);
    s4.u = AddVec4Cart(state.u, du3);
    if (!TryGeodesicAccelerationCart(s4, M, a, Q, Lambda, acceleration)) return false;
    Vec4Cart dx4 = ScaleVec4Cart(s4.u, dt);
    Vec4Cart du4 = ScaleVec4Cart(acceleration, dt);

    Real sixth = 1.0f / 6.0f;
    dx = ScaleVec4Cart(AddVec4Cart(AddVec4Cart(dx1, ScaleVec4Cart(dx2, 2.0f)),
                                   AddVec4Cart(ScaleVec4Cart(dx3, 2.0f), dx4)), sixth);
    du = ScaleVec4Cart(AddVec4Cart(AddVec4Cart(du1, ScaleVec4Cart(du2, 2.0f)),
                                   AddVec4Cart(ScaleVec4Cart(du3, 2.0f), du4)), sixth);
    Vec4Cart candidate = AddVec4Cart(state.x, dx);
    return IsKerrSchildFamilyEventRepresented(candidate, M, a, Q, Lambda);
}

public bool IntegrateGeodesicRK4Cart(GeodesicStateCart state, Real dt,
                                     float M, float a, float Q, float Lambda,
                                     out GeodesicStateCart result) {
    result = state;
    Vec4Cart dx, du;
    if (!IntegrateGeodesicRK4CartDelta(state, dt, M, a, Q, Lambda, dx, du)) return false;
    result.x = AddVec4Cart(state.x, dx);
    result.u = AddVec4Cart(state.u, du);
    return true;
}

// The live Morris-Thorne dispatch is the exact isotropic Ellis chart shared
// with the CPU authority.  Positive-rho second-sheet stages are represented so
// the accepted segment can be clipped at the regular one-sheet throat; rho=0
// and non-Ellis shape selectors fail before metric evaluation.
public bool TryGeodesicAccelerationMorrisThorneCart(
    GeodesicStateCart state, float b0, float Phi0, int shapeType,
    out Vec4Cart acceleration) {
    acceleration = MakeVec4Cart(Real(0.0f), Real(0.0f), Real(0.0f), Real(0.0f));
    if (!IsMorrisThorneCartesianEventRepresented(state.x, b0, Phi0, shapeType)) return false;
    Real Gamma[4][4][4];
    GetMorrisThorneCartesianChristoffel(state.x, b0, Phi0, shapeType, Gamma);
    Real u[4] = { state.u.t, state.u.x, state.u.y, state.u.z };
    Real acc[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
    [unroll] for (int mu = 0; mu < 4; mu++)
        [unroll] for (int nu = 0; nu < 4; nu++)
            [unroll] for (int rho = 0; rho < 4; rho++)
                acc[mu] -= Gamma[mu][nu][rho] * u[nu] * u[rho];
    acceleration = MakeVec4Cart(acc[0], acc[1], acc[2], acc[3]);
    return true;
}

public bool IntegrateGeodesicRK4MorrisThorneCartDelta(
    GeodesicStateCart state, Real dt, float b0, float Phi0, int shapeType,
    out Vec4Cart dx, out Vec4Cart du) {
    dx = MakeVec4Cart(Real(0.0f), Real(0.0f), Real(0.0f), Real(0.0f));
    du = MakeVec4Cart(Real(0.0f), Real(0.0f), Real(0.0f), Real(0.0f));
    Vec4Cart acceleration;
    if (!TryGeodesicAccelerationMorrisThorneCart(state, b0, Phi0, shapeType,
                                                 acceleration)) return false;
    Vec4Cart dx1 = ScaleVec4Cart(state.u, dt);
    Vec4Cart du1 = ScaleVec4Cart(acceleration, dt);

    GeodesicStateCart s2;
    s2.x = AddVec4Cart(state.x, ScaleVec4Cart(dx1, 0.5f));
    s2.u = AddVec4Cart(state.u, ScaleVec4Cart(du1, 0.5f));
    if (!TryGeodesicAccelerationMorrisThorneCart(s2, b0, Phi0, shapeType,
                                                 acceleration)) return false;
    Vec4Cart dx2 = ScaleVec4Cart(s2.u, dt);
    Vec4Cart du2 = ScaleVec4Cart(acceleration, dt);

    GeodesicStateCart s3;
    s3.x = AddVec4Cart(state.x, ScaleVec4Cart(dx2, 0.5f));
    s3.u = AddVec4Cart(state.u, ScaleVec4Cart(du2, 0.5f));
    if (!TryGeodesicAccelerationMorrisThorneCart(s3, b0, Phi0, shapeType,
                                                 acceleration)) return false;
    Vec4Cart dx3 = ScaleVec4Cart(s3.u, dt);
    Vec4Cart du3 = ScaleVec4Cart(acceleration, dt);

    GeodesicStateCart s4;
    s4.x = AddVec4Cart(state.x, dx3);
    s4.u = AddVec4Cart(state.u, du3);
    if (!TryGeodesicAccelerationMorrisThorneCart(s4, b0, Phi0, shapeType,
                                                 acceleration)) return false;
    Vec4Cart dx4 = ScaleVec4Cart(s4.u, dt);
    Vec4Cart du4 = ScaleVec4Cart(acceleration, dt);

    Real sixth = 1.0f / 6.0f;
    dx = ScaleVec4Cart(AddVec4Cart(AddVec4Cart(dx1, ScaleVec4Cart(dx2, 2.0f)),
                                   AddVec4Cart(ScaleVec4Cart(dx3, 2.0f), dx4)), sixth);
    du = ScaleVec4Cart(AddVec4Cart(AddVec4Cart(du1, ScaleVec4Cart(du2, 2.0f)),
                                   AddVec4Cart(ScaleVec4Cart(du3, 2.0f), du4)), sixth);
    Vec4Cart candidate = AddVec4Cart(state.x, dx);
    return IsMorrisThorneCartesianEventRepresented(candidate, b0, Phi0, shapeType);
}

public bool IntegrateGeodesicRK4MorrisThorneCart(
    GeodesicStateCart state, Real dt, float b0, float Phi0, int shapeType,
    out GeodesicStateCart result) {
    result = state;
    Vec4Cart dx, du;
    if (!IntegrateGeodesicRK4MorrisThorneCartDelta(
            state, dt, b0, Phi0, shapeType, dx, du)) return false;
    result.x = AddVec4Cart(state.x, dx);
    result.u = AddVec4Cart(state.u, du);
    return true;
}

// Alcubierre uses its analytic Cartesian Christoffels and is lensing-only. As
// with every other live family, every RK stage must remain inside the complete
// represented parameter/event domain; an invalid profile cannot degrade to
// zero connection and masquerade as flat space.
public bool TryGeodesicAccelerationWarpCart(
    GeodesicStateCart state, float vs, float sigma, float R,
    float xs, float ys, float zs, out Vec4Cart acceleration) {
    acceleration = MakeVec4Cart(Real(0.0f), Real(0.0f), Real(0.0f), Real(0.0f));
    if (!IsWarpDriveEventRepresented(state.x, vs, sigma, R, xs, ys, zs)) return false;
    Real Gamma[4][4][4];
    GetWarpDriveChristoffel(state.x, vs, sigma, R, xs, ys, zs, Gamma);
    Real u[4] = { state.u.t, state.u.x, state.u.y, state.u.z };
    Real acc[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
    [unroll] for (int mu = 0; mu < 4; mu++)
        [unroll] for (int nu = 0; nu < 4; nu++)
            [unroll] for (int rho = 0; rho < 4; rho++)
                acc[mu] -= Gamma[mu][nu][rho] * u[nu] * u[rho];
    acceleration = MakeVec4Cart(acc[0], acc[1], acc[2], acc[3]);
    return isfinite(acceleration.t) && isfinite(acceleration.x) &&
           isfinite(acceleration.y) && isfinite(acceleration.z);
}

public bool IntegrateGeodesicRK4WarpCartDelta(
    GeodesicStateCart state, Real dt, float vs, float sigma, float R,
    float xs, float ys, float zs, out Vec4Cart dx, out Vec4Cart du) {
    dx = MakeVec4Cart(Real(0.0f), Real(0.0f), Real(0.0f), Real(0.0f));
    du = MakeVec4Cart(Real(0.0f), Real(0.0f), Real(0.0f), Real(0.0f));
    Vec4Cart acceleration;
    if (!TryGeodesicAccelerationWarpCart(state, vs, sigma, R, xs, ys, zs,
                                         acceleration)) return false;
    Vec4Cart dx1 = ScaleVec4Cart(state.u, dt);
    Vec4Cart du1 = ScaleVec4Cart(acceleration, dt);

    GeodesicStateCart s2;
    s2.x = AddVec4Cart(state.x, ScaleVec4Cart(dx1, 0.5f));
    s2.u = AddVec4Cart(state.u, ScaleVec4Cart(du1, 0.5f));
    if (!TryGeodesicAccelerationWarpCart(s2, vs, sigma, R, xs, ys, zs,
                                         acceleration)) return false;
    Vec4Cart dx2 = ScaleVec4Cart(s2.u, dt);
    Vec4Cart du2 = ScaleVec4Cart(acceleration, dt);

    GeodesicStateCart s3;
    s3.x = AddVec4Cart(state.x, ScaleVec4Cart(dx2, 0.5f));
    s3.u = AddVec4Cart(state.u, ScaleVec4Cart(du2, 0.5f));
    if (!TryGeodesicAccelerationWarpCart(s3, vs, sigma, R, xs, ys, zs,
                                         acceleration)) return false;
    Vec4Cart dx3 = ScaleVec4Cart(s3.u, dt);
    Vec4Cart du3 = ScaleVec4Cart(acceleration, dt);

    GeodesicStateCart s4;
    s4.x = AddVec4Cart(state.x, dx3);
    s4.u = AddVec4Cart(state.u, du3);
    if (!TryGeodesicAccelerationWarpCart(s4, vs, sigma, R, xs, ys, zs,
                                         acceleration)) return false;
    Vec4Cart dx4 = ScaleVec4Cart(s4.u, dt);
    Vec4Cart du4 = ScaleVec4Cart(acceleration, dt);

    Real sixth = 1.0f / 6.0f;
    dx = ScaleVec4Cart(AddVec4Cart(AddVec4Cart(dx1, ScaleVec4Cart(dx2, 2.0f)),
                                   AddVec4Cart(ScaleVec4Cart(dx3, 2.0f), dx4)), sixth);
    du = ScaleVec4Cart(AddVec4Cart(AddVec4Cart(du1, ScaleVec4Cart(du2, 2.0f)),
                                   AddVec4Cart(ScaleVec4Cart(du3, 2.0f), du4)), sixth);
    Vec4Cart candidate = AddVec4Cart(state.x, dx);
    Vec4Cart candidateTangent = AddVec4Cart(state.u, du);
    return IsWarpDriveEventRepresented(candidate, vs, sigma, R, xs, ys, zs) &&
           isfinite(candidateTangent.t) && isfinite(candidateTangent.x) &&
           isfinite(candidateTangent.y) && isfinite(candidateTangent.z);
}

public bool IntegrateGeodesicRK4WarpCart(
    GeodesicStateCart state, Real dt, float vs, float sigma, float R,
    float xs, float ys, float zs, out GeodesicStateCart result) {
    result = state;
    Vec4Cart dx, du;
    if (!IntegrateGeodesicRK4WarpCartDelta(
            state, dt, vs, sigma, R, xs, ys, zs, dx, du)) return false;
    result.x = AddVec4Cart(state.x, dx);
    result.u = AddVec4Cart(state.u, du);
    return true;
}

public struct AdaptiveGeodesicStep {
    public GeodesicStateCart candidate;
    public Real errorRatio;
    // Defect produced by the unprojected RK4 candidate.  This is the local
    // constraint error judged by the adaptive controller.
    public Real nullResidual;
    // Residual of candidate.u after the accepted-state null projection.
    public Real projectedNullResidual;
    public Real stepFactor;
    public bool accepted;
}

bool AdvanceGeodesicRK4(GeodesicStateCart state, Real step,
                        bool isWarp, bool isWormhole,
                        float M, float a, float Q, float Lambda,
                        float warpVs, float warpSigma, float warpR,
                        float throatB0, int wormShape,
                        out GeodesicStateCart result) {
    if (isWarp) {
        return IntegrateGeodesicRK4WarpCart(state, step, warpVs, warpSigma, warpR,
                                            0.0f, 0.0f, 0.0f, result);
    }
    if (isWormhole) {
        return IntegrateGeodesicRK4MorrisThorneCart(
            state, step, throatB0, 0.0f, wormShape, result);
    }
    return IntegrateGeodesicRK4Cart(state, step, M, a, Q, Lambda, result);
}

Real RelativeStepDifference(GeodesicStateCart start, GeodesicStateCart first,
                            GeodesicStateCart second) {
    Real startValues[8] = { start.x.t, start.x.x, start.x.y, start.x.z,
                            start.u.t, start.u.x, start.u.y, start.u.z };
    Real firstValues[8] = { first.x.t, first.x.x, first.x.y, first.x.z,
                            first.u.t, first.u.x, first.u.y, first.u.z };
    Real secondValues[8] = { second.x.t, second.x.x, second.x.y, second.x.z,
                             second.u.t, second.u.x, second.u.y, second.u.z };
    Real error = Real(0.0f);
    [unroll] for (int component = 0; component < 8; ++component) {
        Real scale = max(Real(1.0f), max(abs(startValues[component]),
                                        abs(secondValues[component])));
        error = max(error, abs(secondValues[component] - firstValues[component]) / scale);
    }
    // One full RK4 step versus two half steps differs by 15 times the leading
    // local truncation term.
    return error / Real(15.0f);
}

public Real RelativeNullConstraint(Real metric[4][4], Vec4Cart tangent) {
    Real components[4] = { tangent.t, tangent.x, tangent.y, tangent.z };
    Real contraction = Real(0.0f);
    Real scale = Real(0.0f);
    [unroll] for (int mu = 0; mu < 4; ++mu) {
        [unroll] for (int nu = 0; nu < 4; ++nu) {
            Real term = metric[mu][nu] * components[mu] * components[nu];
            contraction += term;
            scale += abs(term);
        }
    }
    return abs(contraction) / max(scale, Real(1.0e-20f));
}

Real DeviceRealEpsilon() {
#ifdef SIRIUS_FP64
    return Real(2.2204460492503131e-16);
#else
    return Real(1.1920928955078125e-7f);
#endif
}

// Project one contravariant tangent onto the metric null cone using the nearest
// temporal root when represented. A coordinate-time root need not exist inside
// a stationary limit; in that case the smallest normalized spatial correction
// leaves k^0 unchanged and preserves the incoming cone orientation.
public bool ProjectNullTangentPreservingBranch(Real metric[4][4], Vec4Cart tangent,
                                               out Vec4Cart projected) {
    projected = tangent;
    if (!isfinite(tangent.t) || !isfinite(tangent.x) ||
        !isfinite(tangent.y) || !isfinite(tangent.z)) {
        return false;
    }

    Real components[4] = { tangent.t, tangent.x, tangent.y, tangent.z };
    Real metricScale = Real(0.0f);
    [unroll] for (int mu = 0; mu < 4; ++mu) {
        [unroll] for (int nu = 0; nu < 4; ++nu) {
            if (!isfinite(metric[mu][nu])) return false;
            metricScale = max(metricScale, abs(metric[mu][nu]));
        }
    }
    if (!(metricScale > Real(0.0f))) return false;

    Real roundoff = Real(64.0f) * DeviceRealEpsilon();
    Real bestComponents[4] = { components[0], components[1], components[2], components[3] };
    Real bestCorrection = Real(1.0e30f);
    bool bestRepresented = false;
    bool temporalRepresented = false;
    [unroll] for (int component = 0; component < 4; ++component) {
        if (component > 0 && temporalRepresented) continue;
        Real A = metric[component][component];
        Real B = Real(0.0f);
        Real C = Real(0.0f);
        [unroll] for (int i = 0; i < 4; ++i) {
            if (i == component) continue;
            B += (metric[component][i] + metric[i][component]) * components[i];
            [unroll] for (int j = 0; j < 4; ++j) {
                if (j == component) continue;
                C += metric[i][j] * components[i] * components[j];
            }
        }
        if (!isfinite(A) || !isfinite(B) || !isfinite(C)) continue;

        Real roots[2] = { Real(0.0f), Real(0.0f) };
        int rootCount = 0;
        if (A == Real(0.0f)) {
            if (B == Real(0.0f)) continue;
            roots[rootCount++] = -C / B;
        } else {
            Real fourAC = Real(4.0f) * A * C;
            Real discriminant = B * B - fourAC;
            Real discriminantScale = B * B + abs(fourAC);
            if (!isfinite(discriminant) ||
                discriminant < -roundoff * discriminantScale) {
                continue;
            }
            discriminant = max(discriminant, Real(0.0f));
            Real sqrtDiscriminant = sqrt(discriminant);
            Real q = Real(-0.5f) *
                     (B + (B >= Real(0.0f) ? sqrtDiscriminant : -sqrtDiscriminant));
            if (q == Real(0.0f)) {
                roots[rootCount++] = -B / (Real(2.0f) * A);
                roots[rootCount++] = roots[0];
            } else {
                roots[rootCount++] = q / A;
                roots[rootCount++] = C / q;
            }
        }

        [unroll] for (int rootIndex = 0; rootIndex < 2; ++rootIndex) {
            if (rootIndex >= rootCount) continue;
            Real root = roots[rootIndex];
            if (!isfinite(root)) continue;
            Real correction =
                abs(root - components[component]) /
                (Real(1.0f) + abs(components[component]));
            if (correction >= bestCorrection) continue;
            [unroll] for (int i = 0; i < 4; ++i) bestComponents[i] = components[i];
            bestComponents[component] = root;
            bestCorrection = correction;
            bestRepresented = true;
            if (component == 0) temporalRepresented = true;
        }
    }
    if (!bestRepresented) return false;
    projected = MakeVec4Cart(bestComponents[0], bestComponents[1],
                             bestComponents[2], bestComponents[3]);
    return true;
}

// One device-precision policy consumed by both the production trace and its P1
// live parity probe. Keeping these here prevents a probe from silently judging
// a different adaptive/null controller envelope than the renderer executes.
public Real AdaptiveIntegrationTolerance() {
#ifdef SIRIUS_FP64
    return Real(1.0e-9);
#else
    return Real(2.0e-5f);
#endif
}

public Real AdaptiveNullTolerance() {
#ifdef SIRIUS_FP64
    // Keep a small guard band above the 1e-9 integration target.  Requiring
    // the independently accumulated null contraction to equal that target at
    // the last bit can reject an otherwise converged minimum step (observed
    // residual 1.00002e-9) and strand the trajectory before termination.
    return Real(5.0e-9);
#else
    return Real(5.0e-4f);
#endif
}

// Shared adaptive candidate authority for trace.slang and the P1 parity probe.
// The accepted state is the two-half-step solution; the caller owns retry and
// precision-rung commit so a rejected attempt cannot mutate downstream state.
public AdaptiveGeodesicStep AttemptAdaptiveGeodesicStep(
    GeodesicStateCart state, Real step, bool isWarp, bool isWormhole,
    float M, float a, float Q, float Lambda,
    float warpVs, float warpSigma, float warpR,
    float throatB0, int wormShape, Real integrationTolerance, Real nullTolerance) {
    AdaptiveGeodesicStep result;
    result.candidate = state;
    GeodesicStateCart fullStep = state;
    bool fullRepresented = AdvanceGeodesicRK4(
        state, step, isWarp, isWormhole, M, a, Q, Lambda,
        warpVs, warpSigma, warpR, throatB0, wormShape, fullStep);
    GeodesicStateCart firstHalf = state;
    bool firstHalfRepresented = AdvanceGeodesicRK4(
        state, step * Real(0.5f), isWarp, isWormhole, M, a, Q, Lambda,
        warpVs, warpSigma, warpR, throatB0, wormShape, firstHalf);
    bool candidateRepresented = firstHalfRepresented && AdvanceGeodesicRK4(
        firstHalf, step * Real(0.5f), isWarp, isWormhole, M, a, Q, Lambda,
        warpVs, warpSigma, warpR, throatB0, wormShape, result.candidate);
    if (!fullRepresented || !firstHalfRepresented || !candidateRepresented) {
        result.candidate = state;
        result.errorRatio = Real(1.0e30f);
        result.nullResidual = Real(1.0e30f);
        result.projectedNullResidual = Real(1.0e30f);
        result.stepFactor = Real(0.5f);
        result.accepted = false;
        return result;
    }
    result.errorRatio = RelativeStepDifference(state, fullStep, result.candidate) /
                        integrationTolerance;

    Real candidateMetric[4][4];
    Real candidateInverse[4][4];
    if (isWarp) {
        GetWarpDriveMetric(result.candidate.x, warpVs, warpSigma, warpR,
                           0.0f, 0.0f, 0.0f, candidateMetric, candidateInverse);
    } else if (isWormhole) {
        GetMorrisThorneCartesianMetric(result.candidate.x, throatB0, 0.0f, wormShape,
                                       candidateMetric, candidateInverse);
    } else {
        GetKerrSchildFamilyMetric(result.candidate.x, M, a, Q, Lambda,
                                  candidateMetric, candidateInverse);
    }
    result.nullResidual = RelativeNullConstraint(candidateMetric, result.candidate.u);
    Vec4Cart projectedTangent;
    bool projectionRepresented = ProjectNullTangentPreservingBranch(
        candidateMetric, result.candidate.u, projectedTangent);
    if (projectionRepresented) {
        result.candidate.u = projectedTangent;
        result.projectedNullResidual =
            RelativeNullConstraint(candidateMetric, result.candidate.u);
    } else {
        result.projectedNullResidual = Real(1.0e30f);
    }
    float controllerError = max(float(result.errorRatio), 1.0e-12f);
    result.stepFactor =
        Real(clamp(0.9f * pow(1.0f / controllerError, 0.2f), 0.2f, 2.0f));
    result.accepted = projectionRepresented && result.errorRatio <= Real(1.0f) &&
                      result.nullResidual <= nullTolerance &&
                      result.projectedNullResidual <= nullTolerance;
    if (!result.accepted &&
        (!projectionRepresented || result.nullResidual > nullTolerance ||
         result.projectedNullResidual > nullTolerance)) {
        result.stepFactor = min(result.stepFactor, Real(0.5f));
    }
    return result;
}
