// Two-ended topology authority for the exact isotropic Ellis geometry.  The
// metric itself remains horizonless; this module owns only the selected trace
// boundary, the inversion-related second infinity, and its celestial basis.
module gr_ellis_topology;

import gr_types;
import gr_integrator;
import gr_trace_event;

public static const int kEllisTraceUnresolved = 0;
public static const int kEllisTraceObserverInfinity = 1;
public static const int kEllisTraceOppositeInfinity = 2;

public struct EllisTopologyTraceResult {
    public GeodesicStateCart state;
    public int termination;
    public int attempts;
    public Real affineLength;
};

// Inversion rho' = a^2/rho, a=b0/2, exchanges the two asymptotic ends and
// preserves areal radius.  Matching both ends at the same outer isotropic
// radius therefore places the second boundary at a^2/escapeRadius.
public Real EllisOppositeEscapeRadius(Real b0, Real escapeRadius) {
    return (b0 * b0) / (Real(4.0f) * escapeRadius);
}

// At the second end, convert chart-orthonormal direction components to the
// inversion-related asymptotically Cartesian frame.  The conformal scale
// cancels, leaving the radial Householder reflection I-2nn^T.
public bool MapEllisSecondSheetSkyDirection(
    Vec4Cart position, float3 chartDirection, out float3 mappedDirection) {
    mappedDirection = float3(0.0f, 0.0f, 1.0f);
    Real radius = sqrt(position.x * position.x + position.y * position.y +
                       position.z * position.z);
    if (!(isfinite(radius) && radius > Real(0.0f)) ||
        !all(isfinite(chartDirection))) return false;
    Real3 normal = Real3(position.x, position.y, position.z) / radius;
    Real3 direction = Real3(Real(chartDirection.x), Real(chartDirection.y),
                            Real(chartDirection.z));
    Real norm = sqrt(dot(direction, direction));
    if (!(isfinite(norm) && norm > Real(0.0f))) return false;
    Real3 reflected = direction - Real(2.0f) * dot(normal, direction) * normal;
    mappedDirection = float3(float(reflected.x), float(reflected.y), float(reflected.z));
    return all(isfinite(mappedDirection));
}

// Integrate the production two-sheet central ray without shading.  This exact
// function is consumed by both trace.slang and the compute-only parity probe,
// so the tested trajectory is the one that reaches the renderer.  The Ellis
// profile has no disk, volume, polarisation, or bundle model in Sirius; those
// combinations decline at the typed boundary before this path is entered.
public EllisTopologyTraceResult TraceEllisTwoSheet(
    GeodesicStateCart launch, int maxAttempts, Real stepScale, Real minStep,
    Real maxStep, Real escapeRadius, Real b0, int wormShape) {
    EllisTopologyTraceResult result;
    result.state = launch;
    result.termination = kEllisTraceUnresolved;
    result.attempts = 0;
    result.affineLength = Real(0.0f);

    Real throatRadius = Real(0.5f) * b0;
    Real oppositeRadius = EllisOppositeEscapeRadius(b0, escapeRadius);
    if (!(maxAttempts > 0 && isfinite(stepScale) && stepScale > Real(0.0f) &&
          isfinite(minStep) && minStep > Real(0.0f) && isfinite(maxStep) &&
          maxStep >= minStep && isfinite(escapeRadius) &&
          escapeRadius > throatRadius && isfinite(b0) && b0 > Real(0.0f) &&
          isfinite(oppositeRadius) && oppositeRadius > Real(0.0f) &&
          oppositeRadius < throatRadius && wormShape == 0)) {
        return result;
    }

    Real adaptiveStep = maxStep;
#ifdef SIRIUS_FP32_COMP
    Vec4Cart compensationX = MakeVec4Cart(0.0f, 0.0f, 0.0f, 0.0f);
    Vec4Cart compensationU = MakeVec4Cart(0.0f, 0.0f, 0.0f, 0.0f);
#endif

    for (int attempt = 0; attempt < maxAttempts; ++attempt) {
        result.attempts = attempt + 1;
        Real radius = sqrt(result.state.x.x * result.state.x.x +
                           result.state.x.y * result.state.x.y +
                           result.state.x.z * result.state.x.z);
        Real radialNumerator = result.state.x.x * result.state.u.x +
                               result.state.x.y * result.state.u.y +
                               result.state.x.z * result.state.u.z;
        if (radius >= escapeRadius && radialNumerator > Real(0.0f)) {
            result.termination = kEllisTraceObserverInfinity;
            return result;
        }
        if (radius <= oppositeRadius && radialNumerator < Real(0.0f)) {
            result.termination = kEllisTraceOppositeInfinity;
            return result;
        }

        Real geometricStep = clamp(stepScale * max(radius, b0), minStep, maxStep);
        Real dt = min(geometricStep, adaptiveStep);
        AdaptiveGeodesicStep attempted = AttemptAdaptiveGeodesicStep(
            result.state, dt, false, true, 0.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 0.0f, float(b0), wormShape,
            AdaptiveIntegrationTolerance(), AdaptiveNullTolerance());
        if (!attempted.accepted) {
            if (dt <= minStep * Real(1.0001f)) return result;
            adaptiveStep = max(minStep, dt * attempted.stepFactor);
            continue;
        }
        adaptiveStep = clamp(dt * attempted.stepFactor, minStep, maxStep);

        GeodesicStateCart candidate = attempted.candidate;
#ifdef SIRIUS_FP32_COMP
        Vec4Cart deltaX = AddVec4Cart(
            candidate.x, ScaleVec4Cart(result.state.x, Real(-1.0f)));
        Vec4Cart deltaU = AddVec4Cart(
            candidate.u, ScaleVec4Cart(result.state.u, Real(-1.0f)));
        candidate.x = KahanAddVec4Cart(result.state.x, compensationX, deltaX);
        candidate.u = KahanAddVec4Cart(result.state.u, compensationU, deltaU);
#endif

        Real observerFraction;
        Vec4Cart observerPosition;
        Vec4Cart observerTangent;
        bool observerEscape = FindSphericalBoundaryEvent(
            result.state.x, result.state.u, candidate.x, candidate.u, dt,
            escapeRadius, kSphericalBoundaryIncreasingRadius,
            observerFraction, observerPosition, observerTangent);
        Real oppositeFraction;
        Vec4Cart oppositePosition;
        Vec4Cart oppositeTangent;
        bool oppositeEscape = FindSphericalBoundaryEvent(
            result.state.x, result.state.u, candidate.x, candidate.u, dt,
            oppositeRadius, kSphericalBoundaryDecreasingRadius,
            oppositeFraction, oppositePosition, oppositeTangent);

        if (observerEscape && (!oppositeEscape || observerFraction <= oppositeFraction)) {
            result.state.x = observerPosition;
            result.state.u = observerTangent;
            result.affineLength += dt * observerFraction;
            result.termination = kEllisTraceObserverInfinity;
            return result;
        }
        if (oppositeEscape) {
            result.state.x = oppositePosition;
            result.state.u = oppositeTangent;
            result.affineLength += dt * oppositeFraction;
            result.termination = kEllisTraceOppositeInfinity;
            return result;
        }

        result.state = candidate;
        result.affineLength += dt;
    }
    return result;
}
