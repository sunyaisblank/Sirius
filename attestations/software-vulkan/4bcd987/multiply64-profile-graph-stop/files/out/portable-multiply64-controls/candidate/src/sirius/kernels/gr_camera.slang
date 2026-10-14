// Physical camera projection, observer frame, ray and labelled screens shared
// by trace and the compute probes. Camera parameters retain their float ABI.
import gr_types;
import gr_deviation;

// Finite-pupil camera projection shared by the live trace kernel and its
// compute-only parity probe. Host core/camera.h remains the independent
// double-checking implementation.

public struct ThinLensProjectionSample {
    public float pupilRight;
    public float pupilUp;
    public float3 direction;  // camera-rest-frame forward, up, right
};

public ThinLensProjectionSample ProjectThinLensSample(
    float imageRight, float imageUp, float tanHalfFov, float focalLength,
    float aperture, float focusDistance, float sampleU, float sampleV) {
    float focusRight = imageRight * tanHalfFov * focusDistance;
    float focusUp = imageUp * tanHalfFov * focusDistance;
    // 50 mm-equivalent defines one virtual geometric lens unit because no
    // physical mass-to-millimetre scale exists at this boundary.
    float apertureRadius = (focalLength / 50.0f) / (2.0f * aperture);
    float pupilAngle = 2.0f * 3.14159265358979323846f * sampleU;
    float pupilRadius = apertureRadius * sqrt(sampleV);

    ThinLensProjectionSample sample;
    sample.pupilRight = pupilRadius * cos(pupilAngle);
    sample.pupilUp = pupilRadius * sin(pupilAngle);
    sample.direction = normalize(float3(
        focusDistance,
        focusUp - sample.pupilUp,
        focusRight - sample.pupilRight));
    return sample;
}

public struct CameraObserverFrame {
    public Vec4Cart time;
    public Vec4Cart forward;
    public Vec4Cart up;
    public Vec4Cart right;
};

public Real MetricDotCart(Real metric[4][4], Vec4Cart lhs, Vec4Cart rhs) {
    Real a[4] = { lhs.t, lhs.x, lhs.y, lhs.z };
    Real b[4] = { rhs.t, rhs.x, rhs.y, rhs.z };
    Real result = Real(0.0f);
    [unroll] for (int mu = 0; mu < 4; ++mu) {
        [unroll] for (int nu = 0; nu < 4; ++nu) {
            result += metric[mu][nu] * a[mu] * b[nu];
        }
    }
    return result;
}

public Vec4Cart NormaliseSpacelike(Real metric[4][4], Vec4Cart vector) {
    return ScaleVec4Cart(vector, Real(1.0f) / sqrt(MetricDotCart(metric, vector, vector)));
}

// Metric-orthonormal screen frame followed by a local Lorentz boost. This is
// the shader counterpart of core/observer_frame.h. cameraBeta is expressed in
// screen-forward/up/right axes and is validated by the typed host boundary.
public CameraObserverFrame BuildCameraObserverFrame(Real metric[4][4], Real inverseMetric[4][4],
                                              float3 camDir,
                                              float3 camUp, float3 camRight,
                                              float3 cameraBeta) {
    CameraObserverFrame reference;
    Real lapse = Real(1.0f) / sqrt(-inverseMetric[0][0]);
    reference.time = MakeVec4Cart(-lapse * inverseMetric[0][0],
                                  -lapse * inverseMetric[1][0],
                                  -lapse * inverseMetric[2][0],
                                  -lapse * inverseMetric[3][0]);
    reference.forward =
        MakeVec4Cart(Real(0.0f), Real(camDir.x), Real(camDir.y), Real(camDir.z));
    reference.forward = AddVec4Cart(
        reference.forward,
        ScaleVec4Cart(reference.time,
                      MetricDotCart(metric, reference.forward, reference.time)));
    reference.forward = NormaliseSpacelike(metric, reference.forward);

    // The host's camUp vector is the +theta (screen-down) axis.
    reference.up =
        MakeVec4Cart(Real(0.0f), Real(-camUp.x), Real(-camUp.y), Real(-camUp.z));
    reference.up = AddVec4Cart(
        reference.up,
        ScaleVec4Cart(reference.time, MetricDotCart(metric, reference.up, reference.time)));
    reference.up = AddVec4Cart(
        reference.up,
        ScaleVec4Cart(reference.forward,
                      -MetricDotCart(metric, reference.up, reference.forward)));
    reference.up = NormaliseSpacelike(metric, reference.up);

    reference.right =
        MakeVec4Cart(Real(0.0f), Real(camRight.x), Real(camRight.y), Real(camRight.z));
    reference.right = AddVec4Cart(
        reference.right,
        ScaleVec4Cart(reference.time, MetricDotCart(metric, reference.right, reference.time)));
    reference.right = AddVec4Cart(
        reference.right,
        ScaleVec4Cart(reference.forward,
                      -MetricDotCart(metric, reference.right, reference.forward)));
    reference.right = AddVec4Cart(
        reference.right,
        ScaleVec4Cart(reference.up, -MetricDotCart(metric, reference.right, reference.up)));
    reference.right = NormaliseSpacelike(metric, reference.right);

    // Preserve the float camera ABI, but form the boost norm in active precision.
    Real3 representedBeta = Real3(Real(cameraBeta.x), Real(cameraBeta.y), Real(cameraBeta.z));
    Real betaSquared = dot(representedBeta, representedBeta);
    if (betaSquared == Real(0.0f)) return reference;
    Real gamma = Real(1.0f) / sqrt(Real(1.0f) - betaSquared);
    Vec4Cart betaVector = AddVec4Cart(
        AddVec4Cart(ScaleVec4Cart(reference.forward, Real(cameraBeta.x)),
                    ScaleVec4Cart(reference.up, Real(cameraBeta.y))),
        ScaleVec4Cart(reference.right, Real(cameraBeta.z)));
    CameraObserverFrame boosted;
    boosted.time = ScaleVec4Cart(AddVec4Cart(reference.time, betaVector), gamma);
    Real coefficient = (gamma - Real(1.0f)) / betaSquared;
    boosted.forward = AddVec4Cart(
        AddVec4Cart(reference.forward,
                    ScaleVec4Cart(betaVector, coefficient * Real(cameraBeta.x))),
        ScaleVec4Cart(reference.time, gamma * Real(cameraBeta.x)));
    boosted.up = AddVec4Cart(
        AddVec4Cart(reference.up,
                    ScaleVec4Cart(betaVector, coefficient * Real(cameraBeta.y))),
        ScaleVec4Cart(reference.time, gamma * Real(cameraBeta.y)));
    boosted.right = AddVec4Cart(
        AddVec4Cart(reference.right,
                    ScaleVec4Cart(betaVector, coefficient * Real(cameraBeta.z))),
        ScaleVec4Cart(reference.time, gamma * Real(cameraBeta.z)));
    return boosted;
}

public Vec4Cart PastDirectedCameraRay(CameraObserverFrame frame, float3 screenDirection) {
    Real3 direction = Real3(Real(screenDirection.x), Real(screenDirection.y),
                            Real(screenDirection.z));
    direction /= sqrt(dot(direction, direction));
    return AddVec4Cart(
        AddVec4Cart(
            AddVec4Cart(ScaleVec4Cart(frame.time, Real(-1.0f)),
                        ScaleVec4Cart(frame.forward, direction.x)),
            ScaleVec4Cart(frame.up, direction.y)),
        ScaleVec4Cart(frame.right, direction.z));
}

public void SachsScreenBasis(CameraObserverFrame frame, float3 restDirection,
                      out Vec4Cart first, out Vec4Cart second) {
    Real3 firstComponents;
    Real3 secondComponents;
    MakeCelestialTangentBasis(Real3(Real(restDirection.x), Real(restDirection.y),
                                    Real(restDirection.z)), firstComponents, secondComponents);
    first = AddVec4Cart(
        AddVec4Cart(ScaleVec4Cart(frame.forward, Real(firstComponents.x)),
                    ScaleVec4Cart(frame.up, Real(firstComponents.y))),
        ScaleVec4Cart(frame.right, Real(firstComponents.z)));
    second = AddVec4Cart(
        AddVec4Cart(ScaleVec4Cart(frame.forward, Real(secondComponents.x)),
                    ScaleVec4Cart(frame.up, Real(secondComponents.y))),
        ScaleVec4Cart(frame.right, Real(secondComponents.z)));
}
