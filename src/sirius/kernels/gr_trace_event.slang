// Accepted-segment geometry shared by the production trace and compute-only
// parity probe.  In particular, a finite capture event is a property of the
// cubic-Hermite segment, not merely of its endpoint.
module gr_trace_event;

import gr_types;

bool IsFiniteAcceptedSegmentValue(Vec4Cart value) {
    return isfinite(value.t) && isfinite(value.x) &&
           isfinite(value.y) && isfinite(value.z);
}

public Vec4Cart AcceptedSegmentPosition(Vec4Cart start, Vec4Cart startVelocity,
                                        Vec4Cart finish, Vec4Cart finishVelocity,
                                        Real affineLength, Real fraction) {
    Real s2 = fraction * fraction;
    Real s3 = s2 * fraction;
    Real h00 = Real(2.0f) * s3 - Real(3.0f) * s2 + Real(1.0f);
    Real h10 = s3 - Real(2.0f) * s2 + fraction;
    Real h01 = Real(-2.0f) * s3 + Real(3.0f) * s2;
    Real h11 = s3 - s2;
    return AddVec4Cart(
        AddVec4Cart(ScaleVec4Cart(start, h00),
                    ScaleVec4Cart(startVelocity, affineLength * h10)),
        AddVec4Cart(ScaleVec4Cart(finish, h01),
                    ScaleVec4Cart(finishVelocity, affineLength * h11)));
}

public Vec4Cart AcceptedSegmentTangent(Vec4Cart start, Vec4Cart startVelocity,
                                       Vec4Cart finish, Vec4Cart finishVelocity,
                                       Real affineLength, Real fraction) {
    Real s2 = fraction * fraction;
    Real dh00 = Real(6.0f) * s2 - Real(6.0f) * fraction;
    Real dh10 = Real(3.0f) * s2 - Real(4.0f) * fraction + Real(1.0f);
    Real dh01 = Real(-6.0f) * s2 + Real(6.0f) * fraction;
    Real dh11 = Real(3.0f) * s2 - Real(2.0f) * fraction;
    return AddVec4Cart(
        AddVec4Cart(ScaleVec4Cart(start, dh00 / affineLength),
                    ScaleVec4Cart(startVelocity, dh10)),
        AddVec4Cart(ScaleVec4Cart(finish, dh01 / affineLength),
                    ScaleVec4Cart(finishVelocity, dh11)));
}

Real PolynomialEpsilon() {
#ifdef SIRIUS_FP64
    return Real(2.2204460492503131e-16);
#else
    return Real(1.1920928955078125e-7f);
#endif
}

Real EvaluatePolynomial(Real coefficients[7], int degree, Real argument) {
    Real value = coefficients[degree];
    for (int power = degree - 1; power >= 0; --power) {
        value = value * argument + coefficients[power];
    }
    return value;
}

Real PolynomialValueTolerance(Real coefficients[7], int degree) {
    Real scale = Real(0.0f);
    for (int power = 0; power <= degree; ++power) {
        scale = max(scale, abs(coefficients[power]));
    }
    return Real(256.0f) * PolynomialEpsilon() * max(scale, Real(1.0e-30f));
}

void AppendUnitRoot(inout Real roots[6], inout int rootCount, Real root,
                    Real fractionTolerance) {
    root = clamp(root, Real(0.0f), Real(1.0f));
    if (rootCount > 0 && abs(root - roots[rootCount - 1]) <= fractionTolerance) return;
    if (rootCount < 6) roots[rootCount++] = root;
}

int FindLinearRoots(Real coefficients[7], out Real roots[6]) {
    [unroll] for (int index = 0; index < 6; ++index) roots[index] = Real(0.0f);
    Real tolerance = PolynomialValueTolerance(coefficients, 1);
    Real fractionTolerance = Real(512.0f) * PolynomialEpsilon();
    if (abs(coefficients[1]) <= tolerance) return 0;
    Real root = -coefficients[0] / coefficients[1];
    int count = 0;
    if (isfinite(root) && root >= -fractionTolerance &&
        root <= Real(1.0f) + fractionTolerance) {
        AppendUnitRoot(roots, count, root, fractionTolerance);
    }
    return count;
}

int FindRootsFromCriticalPoints(Real coefficients[7], int degree,
                                Real criticalRoots[6], int criticalCount,
                                out Real roots[6]) {
    [unroll] for (int index = 0; index < 6; ++index) roots[index] = Real(0.0f);
    Real valueTolerance = PolynomialValueTolerance(coefficients, degree);
    Real fractionTolerance = Real(512.0f) * PolynomialEpsilon();
    Real boundaries[8];
    [unroll] for (int index = 0; index < 8; ++index) boundaries[index] = Real(0.0f);
    int boundaryCount = 1;
    boundaries[0] = Real(0.0f);
    for (int index = 0; index < criticalCount; ++index) {
        Real point = criticalRoots[index];
        if (point > fractionTolerance && point < Real(1.0f) - fractionTolerance) {
            boundaries[boundaryCount++] = point;
        }
    }
    boundaries[boundaryCount++] = Real(1.0f);

    int rootCount = 0;
    for (int interval = 0; interval + 1 < boundaryCount; ++interval) {
        Real lower = boundaries[interval];
        Real upper = boundaries[interval + 1];
        Real lowerValue = EvaluatePolynomial(coefficients, degree, lower);
        Real upperValue = EvaluatePolynomial(coefficients, degree, upper);
        if (abs(lowerValue) <= valueTolerance) {
            AppendUnitRoot(roots, rootCount, lower, fractionTolerance);
        }
        bool signChange = (lowerValue < -valueTolerance && upperValue > valueTolerance) ||
                          (lowerValue > valueTolerance && upperValue < -valueTolerance);
        if (signChange) {
            Real rootLower = lower;
            Real rootUpper = upper;
            Real rootLowerValue = lowerValue;
            for (int iteration = 0; iteration < 48; ++iteration) {
                Real midpoint = (rootLower + rootUpper) * Real(0.5f);
                Real midpointValue = EvaluatePolynomial(coefficients, degree, midpoint);
                if ((rootLowerValue < Real(0.0f)) == (midpointValue < Real(0.0f))) {
                    rootLower = midpoint;
                    rootLowerValue = midpointValue;
                } else {
                    rootUpper = midpoint;
                }
            }
            AppendUnitRoot(roots, rootCount, (rootLower + rootUpper) * Real(0.5f),
                           fractionTolerance);
        }
    }
    if (abs(EvaluatePolynomial(coefficients, degree, Real(1.0f))) <= valueTolerance) {
        AppendUnitRoot(roots, rootCount, Real(1.0f), fractionTolerance);
    }
    return rootCount;
}

int FindQuadraticRoots(Real coefficients[7], out Real roots[6]) {
    if (abs(coefficients[2]) <= PolynomialValueTolerance(coefficients, 2)) {
        return FindLinearRoots(coefficients, roots);
    }
    Real derivative[7];
    [unroll] for (int index = 0; index < 7; ++index) derivative[index] = Real(0.0f);
    derivative[0] = coefficients[1];
    derivative[1] = Real(2.0f) * coefficients[2];
    Real critical[6];
    int criticalCount = FindLinearRoots(derivative, critical);
    return FindRootsFromCriticalPoints(coefficients, 2, critical, criticalCount, roots);
}

int FindCubicRoots(Real coefficients[7], out Real roots[6]) {
    if (abs(coefficients[3]) <= PolynomialValueTolerance(coefficients, 3)) {
        return FindQuadraticRoots(coefficients, roots);
    }
    Real derivative[7];
    [unroll] for (int index = 0; index < 7; ++index) derivative[index] = Real(0.0f);
    for (int power = 1; power <= 3; ++power)
        derivative[power - 1] = Real(power) * coefficients[power];
    Real critical[6];
    int criticalCount = FindQuadraticRoots(derivative, critical);
    return FindRootsFromCriticalPoints(coefficients, 3, critical, criticalCount, roots);
}

int FindQuarticRoots(Real coefficients[7], out Real roots[6]) {
    if (abs(coefficients[4]) <= PolynomialValueTolerance(coefficients, 4)) {
        return FindCubicRoots(coefficients, roots);
    }
    Real derivative[7];
    [unroll] for (int index = 0; index < 7; ++index) derivative[index] = Real(0.0f);
    for (int power = 1; power <= 4; ++power)
        derivative[power - 1] = Real(power) * coefficients[power];
    Real critical[6];
    int criticalCount = FindCubicRoots(derivative, critical);
    return FindRootsFromCriticalPoints(coefficients, 4, critical, criticalCount, roots);
}

int FindQuinticRoots(Real coefficients[7], out Real roots[6]) {
    if (abs(coefficients[5]) <= PolynomialValueTolerance(coefficients, 5)) {
        return FindQuarticRoots(coefficients, roots);
    }
    Real derivative[7];
    [unroll] for (int index = 0; index < 7; ++index) derivative[index] = Real(0.0f);
    for (int power = 1; power <= 5; ++power)
        derivative[power - 1] = Real(power) * coefficients[power];
    Real critical[6];
    int criticalCount = FindQuarticRoots(derivative, critical);
    return FindRootsFromCriticalPoints(coefficients, 5, critical, criticalCount, roots);
}

int FindSexticRoots(Real coefficients[7], out Real roots[6]) {
    if (abs(coefficients[6]) <= PolynomialValueTolerance(coefficients, 6)) {
        return FindQuinticRoots(coefficients, roots);
    }
    Real derivative[7];
    [unroll] for (int index = 0; index < 7; ++index) derivative[index] = Real(0.0f);
    for (int power = 1; power <= 6; ++power)
        derivative[power - 1] = Real(power) * coefficients[power];
    Real critical[6];
    int criticalCount = FindQuinticRoots(derivative, critical);
    return FindRootsFromCriticalPoints(coefficients, 6, critical, criticalCount, roots);
}

void SphericalSegmentPolynomial(Vec4Cart start, Vec4Cart startVelocity,
                                Vec4Cart finish, Vec4Cart finishVelocity,
                                Real affineLength, Real captureRadius,
                                out Real coefficients[7]) {
    [unroll] for (int index = 0; index < 7; ++index) coefficients[index] = Real(0.0f);
    coefficients[0] = -captureRadius * captureRadius;
    Real startComponents[3] = { start.x, start.y, start.z };
    Real finishComponents[3] = { finish.x, finish.y, finish.z };
    Real startTangents[3] = { startVelocity.x, startVelocity.y, startVelocity.z };
    Real finishTangents[3] = { finishVelocity.x, finishVelocity.y, finishVelocity.z };
    for (int axis = 0; axis < 3; ++axis) {
        Real startDelta = affineLength * startTangents[axis];
        Real finishDelta = affineLength * finishTangents[axis];
        Real component[4];
        component[0] = startComponents[axis];
        component[1] = startDelta;
        component[2] = Real(-3.0f) * startComponents[axis]
                     + Real(3.0f) * finishComponents[axis]
                     - Real(2.0f) * startDelta - finishDelta;
        component[3] = Real(2.0f) * startComponents[axis]
                     - Real(2.0f) * finishComponents[axis]
                     + startDelta + finishDelta;
        for (int left = 0; left <= 3; ++left)
            for (int right = 0; right <= 3; ++right)
                coefficients[left + right] += component[left] * component[right];
    }
}

public static const int kSphericalBoundaryAnyContact = 0;
public static const int kSphericalBoundaryIncreasingRadius = 1;
public static const int kSphericalBoundaryDecreasingRadius = -1;

// Returns the first selected contact with |x|=boundaryRadius in
// observer-to-scene order. Squared radius minus R^2 is a sextic; derivative
// roots partition it into monotone intervals, so equal-sign endpoints cannot
// conceal an interior contact. Directional modes exclude tangencies.
public bool FindSphericalBoundaryEvent(
    Vec4Cart start, Vec4Cart startVelocity, Vec4Cart finish, Vec4Cart finishVelocity,
    Real affineLength, Real boundaryRadius, int radialSense, out Real fraction,
    out Vec4Cart eventPosition, out Vec4Cart eventTangent) {
    fraction = Real(0.0f);
    eventPosition = start;
    eventTangent = startVelocity;
    if (!(isfinite(affineLength) && affineLength > Real(0.0f) &&
          isfinite(boundaryRadius) && boundaryRadius > Real(0.0f) &&
          (radialSense == kSphericalBoundaryAnyContact ||
           radialSense == kSphericalBoundaryIncreasingRadius ||
           radialSense == kSphericalBoundaryDecreasingRadius) &&
          IsFiniteAcceptedSegmentValue(start) &&
          IsFiniteAcceptedSegmentValue(startVelocity) &&
          IsFiniteAcceptedSegmentValue(finish) &&
          IsFiniteAcceptedSegmentValue(finishVelocity))) return false;

    Real coefficients[7];
    SphericalSegmentPolynomial(start, startVelocity, finish, finishVelocity,
                               affineLength, boundaryRadius, coefficients);
    Real roots[6];
    int rootCount = FindSexticRoots(coefficients, roots);
    if (rootCount == 0) return false;

    for (int root = 0; root < rootCount; ++root) {
        Real candidateFraction = roots[root];
        Vec4Cart candidatePosition = AcceptedSegmentPosition(
            start, startVelocity, finish, finishVelocity, affineLength,
            candidateFraction);
        Vec4Cart candidateTangent = AcceptedSegmentTangent(
            start, startVelocity, finish, finishVelocity, affineLength,
            candidateFraction);
        Real sampledRadius = sqrt(candidatePosition.x * candidatePosition.x +
                                  candidatePosition.y * candidatePosition.y +
                                  candidatePosition.z * candidatePosition.z);
        if (!(isfinite(sampledRadius) && sampledRadius > Real(0.0f))) continue;
        Real radialNumerator = candidatePosition.x * candidateTangent.x +
                               candidatePosition.y * candidateTangent.y +
                               candidatePosition.z * candidateTangent.z;
        if (!isfinite(radialNumerator)) continue;
        if (radialSense == kSphericalBoundaryIncreasingRadius &&
            !(radialNumerator > Real(0.0f))) continue;
        if (radialSense == kSphericalBoundaryDecreasingRadius &&
            !(radialNumerator < Real(0.0f))) continue;
        fraction = candidateFraction;
        eventPosition = candidatePosition;
        eventTangent = candidateTangent;
        Real surfaceScale = boundaryRadius / sampledRadius;
        eventPosition.x *= surfaceScale;
        eventPosition.y *= surfaceScale;
        eventPosition.z *= surfaceScale;
        return true;
    }
    return false;
}

// OneSheetCapture treats entry, exit, and tangent throat contacts alike.
public bool FindSphericalCaptureEvent(
    Vec4Cart start, Vec4Cart startVelocity, Vec4Cart finish, Vec4Cart finishVelocity,
    Real affineLength, Real captureRadius, out Real fraction,
    out Vec4Cart eventPosition, out Vec4Cart eventTangent) {
    return FindSphericalBoundaryEvent(
        start, startVelocity, finish, finishVelocity, affineLength, captureRadius,
        kSphericalBoundaryAnyContact, fraction, eventPosition, eventTangent);
}
