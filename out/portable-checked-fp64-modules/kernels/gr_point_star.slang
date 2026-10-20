// Angular point-source filtering shared by the live tracer and direct probe.
// This fixes geometric precision; the catalogue flux/display calibration is
// deliberately supplied by the caller.
module gr_point_star;

public float PointStarAngularWeight(float3 direction, float3 star,
                                    float3 tangentA, float3 tangentB,
                                    float cutoff, float cosOrientation, float sinOrientation,
                                    float invMajorSquared, float invMinorSquared) {
    float normSquared = dot(direction, direction);
    if (normSquared < 0.25f || normSquared > 4.0f) direction /= sqrt(normSquared);
    // acos(dot) loses finite subpixel offsets when dot rounds to one. Cross
    // the difference first to preserve first-order information near alignment.
    float3 normal = cross(direction, star - direction);
    float sine = length(normal);
    float angle = atan2(sine, dot(direction, star));
    if (angle > cutoff) return 0.0f;
    float2 offset = float2(0.0f);
    if (sine > 0.0f) {
        float3 tangent = cross(normal, normalize(direction)) / sine;
        offset = angle * float2(dot(tangent, tangentA), dot(tangent, tangentB));
    } else if (angle > 0.0f) {
        // The anisotropic logarithmic map has no unique tangent at the
        // antipode. Exclude that boundary; retain the circular limit.
        if (invMajorSquared != invMinorSquared) return 0.0f;
        offset.x = angle;
    }
    float alongMajor = cosOrientation * offset.x + sinOrientation * offset.y;
    float alongMinor = -sinOrientation * offset.x + cosOrientation * offset.y;
    return exp(-0.5f * (alongMajor * alongMajor * invMajorSquared +
                        alongMinor * alongMinor * invMinorSquared));
}
