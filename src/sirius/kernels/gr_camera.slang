// Finite-pupil camera projection shared by the live trace kernel and its
// compute-only parity probe. Host core/camera.h remains the independent
// double-checking implementation.

struct ThinLensProjectionSample {
    float pupilRight;
    float pupilUp;
    float3 direction;  // camera-rest-frame forward, up, right
};

ThinLensProjectionSample ProjectThinLensSample(
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
