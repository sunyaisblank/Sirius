// Per-pixel geodesic tracer, the Vulkan render path. Each thread renders one
// pixel of a TILE: it builds the camera null ray, integrates it through the
// dispatched metric family in the Cartesian chart, and terminates on horizon
// capture or escape. Captured rays shade black (the shadow); escaping rays shade
// by their asymptotic direction against the starfield; an optional thin disk
// adds Novikov-Thorne blackbody emission with g^4 relativistic beaming where the
// geodesic crosses the equatorial plane. Volumetric mode instead integrates
// emission/absorption through every geodesic segment, with the same optional
// procedural density modulation represented by the CPU tracer.
// The kernel returns LINEAR radiance per pixel; tonemapping and grading are the
// host display pipeline's job.
//
// This is render orchestration, not a physics authority: the metric,
// Christoffel, integrator, blackbody, and disk-redshift routines it calls are
// the ported authorities in the gr_* modules. Termination radii and adaptive
// tolerances arrive through the host/session ABI or are explicit per-precision
// device policy pinned by the live kernel-parity probes.
//
// Metric dispatch (params[2]) matches the registry gpu_supported set the render
// session wires: Kerr-Schild family (Minkowski/Schwarzschild/Kerr), the
// Alcubierre warp family, and the Morris-Thorne wormhole through its Cartesian
// embedding (gr_metrics GetMorrisThorneCartesian*, mirroring the CPU render
// authority; one sheet, throat capture). The spherical-chart Ellis/Morris-
// Thorne routines remain parity-probe references, not render dispatches.
//
// Buffer ABI (bound 0..5 in declaration order):
//   0 radiance   : RWStructuredBuffer<float>, tileWidth*tileHeight*4 (RGBA/px).
//   1 params     : StructuredBuffer<float>, layout documented below.
//   2 starfield  : StructuredBuffer<uint>, packed RGBA8 equirectangular texels
//                  (R | G<<8 | B<<16 | A<<24), size starfieldWidth*starfieldHeight;
//                  bound with a 1-texel dummy when starfieldEnabled is 0.
//   3 pointStars : StructuredBuffer<float>, eight floats per StarEntry.
//   4 pointOffsets/5 pointIndices: the exact latitude/longitude CSR index
//                  built by core::StarfieldSpatialIndex; dummy buffers when
//                  pointStarfieldEnabled is 0.
import gr_types;
import gr_metrics;
import gr_integrator;
import gr_disk;
import gr_deviation;

RWStructuredBuffer<float> radiance;  // binding 0
StructuredBuffer<float> params;      // binding 1
StructuredBuffer<uint> starfield;    // binding 2
StructuredBuffer<float> pointStars;  // binding 3
StructuredBuffer<uint> pointOffsets; // binding 4
StructuredBuffer<uint> pointIndices; // binding 5

// Kernel-side metric dispatch ids (subset of gr_types kMetric*; the ids the
// render session emits). Kept local so this file reads without cross-reference.
static const uint kDispatchKerrSchild = 0;
static const uint kDispatchWarpDrive = 2;
static const uint kDispatchMorrisThorne = 3;

float3 P3(uint base) { return float3(params[base], params[base + 1], params[base + 2]); }

struct CameraObserverFrame {
    Vec4Cart time;
    Vec4Cart forward;
    Vec4Cart up;
    Vec4Cart right;
};

Real MetricDotCart(Real metric[4][4], Vec4Cart lhs, Vec4Cart rhs) {
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

Vec4Cart NormaliseSpacelike(Real metric[4][4], Vec4Cart vector) {
    return ScaleVec4Cart(vector, Real(1.0f) / sqrt(MetricDotCart(metric, vector, vector)));
}

// Metric-orthonormal screen frame followed by a local Lorentz boost. This is
// the shader counterpart of core/observer_frame.h. cameraBeta is expressed in
// screen-forward/up/right axes and is validated by the typed host boundary.
CameraObserverFrame BuildCameraObserverFrame(Real metric[4][4], Real inverseMetric[4][4],
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

    Real betaSquared = Real(dot(cameraBeta, cameraBeta));
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

Vec4Cart PastDirectedCameraRay(CameraObserverFrame frame, float3 screenDirection) {
    float3 direction = normalize(screenDirection);
    return AddVec4Cart(
        AddVec4Cart(
            AddVec4Cart(ScaleVec4Cart(frame.time, Real(-1.0f)),
                        ScaleVec4Cart(frame.forward, Real(direction.x))),
            ScaleVec4Cart(frame.up, Real(direction.y))),
        ScaleVec4Cart(frame.right, Real(direction.z)));
}

void SachsScreenBasis(CameraObserverFrame frame, float3 restDirection,
                      out Vec4Cart first, out Vec4Cart second) {
    float3 direction = normalize(restDirection);
    float3 reference;
    if (abs(direction.x) <= abs(direction.y) && abs(direction.x) <= abs(direction.z)) {
        reference = float3(1.0f, 0.0f, 0.0f);
    } else if (abs(direction.y) <= abs(direction.z)) {
        reference = float3(0.0f, 1.0f, 0.0f);
    } else {
        reference = float3(0.0f, 0.0f, 1.0f);
    }
    float3 firstComponents = normalize(reference - dot(reference, direction) * direction);
    float3 secondComponents = cross(direction, firstComponents);
    first = AddVec4Cart(
        AddVec4Cart(ScaleVec4Cart(frame.forward, Real(firstComponents.x)),
                    ScaleVec4Cart(frame.up, Real(firstComponents.y))),
        ScaleVec4Cart(frame.right, Real(firstComponents.z)));
    second = AddVec4Cart(
        AddVec4Cart(ScaleVec4Cart(frame.forward, Real(secondComponents.x)),
                    ScaleVec4Cart(frame.up, Real(secondComponents.y))),
        ScaleVec4Cart(frame.right, Real(secondComponents.z)));
}

// Smooth escape-direction background: a bounded gradient so the escaped field is
// finite and non-constant across directions (the smoke gate, and the analytic
// background the CPU-vs-Vulkan parity gate compares against). Used when no
// starfield texture is uploaded (starfieldEnabled == 0).
float3 SampleBackgroundGradient(float3 dir) {
    return float3(0.5f, 0.5f, 0.5f) + 0.5f * dir;
}

// Unpack one packed-RGBA8 starfield texel to linear-ish [0,1] RGB, matching the
// CPU path's byte/255 read (RenderSession::SampleStarfield).
float3 StarfieldTexel(int x, int y, int width) {
    uint packed = starfield[y * width + x];
    float r = float(packed & 0xFFu) / 255.0f;
    float g = float((packed >> 8) & 0xFFu) / 255.0f;
    float b = float((packed >> 16) & 0xFFu) / 255.0f;
    return float3(r, g, b);
}

// Equirectangular bilinear starfield lookup along an escape direction. Byte for
// byte the CPU projection (RenderSession::SampleStarfield): theta = acos(z),
// phi = atan2(y, x) wrapped to [0, 2pi), (u, v) -> pixel with clamp-to-edge
// neighbours and (1-fx)(1-fy) bilinear weights.
float3 SampleStarfield(float3 dir, int width, int height) {
    float len = length(dir);
    if (len < 1e-10f) {
        return float3(0.001f, 0.001f, 0.001f);
    }
    float3 d = dir / len;

    float theta = acos(clamp(d.z, -1.0f, 1.0f));
    float phi = atan2(d.y, d.x);
    if (phi < 0.0f) phi += TWO_PI;

    float u = phi / TWO_PI;
    float v = theta / PI;

    float px = u * float(width - 1);
    float py = v * float(height - 1);

    int x0 = int(floor(px));
    int y0 = int(floor(py));
    int x1 = min(x0 + 1, width - 1);
    int y1 = min(y0 + 1, height - 1);

    float fx = px - float(x0);
    float fy = py - float(y0);

    float3 c00 = StarfieldTexel(x0, y0, width);
    float3 c10 = StarfieldTexel(x1, y0, width);
    float3 c01 = StarfieldTexel(x0, y1, width);
    float3 c11 = StarfieldTexel(x1, y1, width);

    float w00 = (1.0f - fx) * (1.0f - fy);
    float w10 = fx * (1.0f - fy);
    float w01 = (1.0f - fx) * fy;
    float w11 = fx * fy;

    return c00 * w00 + c10 * w10 + c01 * w01 + c11 * w11;
}

float3 PointStarColor(uint starIndex) {
    uint base = starIndex * 8u;
    float temperature = pointStars[base + 6u];
    return BlackbodyColor(temperature);
}

void AccumulatePointStarCell(int thetaBin, int phiBin, float3 dir,
                             float3 tangentA, float3 tangentB, float cosCut,
                             float cosOrientation, float sinOrientation,
                             float invMajorSquared, float invMinorSquared,
                             float brightnessScale, inout float3 color) {
    static const int kPhiBins = 512;
    uint cell = uint(thetaBin * kPhiBins + phiBin);
    for (uint at = pointOffsets[cell]; at < pointOffsets[cell + 1u]; ++at) {
        uint starIndex = pointIndices[at];
        uint base = starIndex * 8u;
        float cosine = dot(dir, float3(pointStars[base], pointStars[base + 1u],
                                      pointStars[base + 2u]));
        if (cosine < cosCut) continue;
        float angle = acos(clamp(cosine, -1.0f, 1.0f));
        float tangentX = 0.0f;
        float tangentY = 0.0f;
        float sinAngle = sin(angle);
        if (sinAngle > 1.0e-8f) {
            float3 starDirection =
                float3(pointStars[base], pointStars[base + 1u],
                       pointStars[base + 2u]);
            float3 tangentDirection = (starDirection - cosine * dir) / sinAngle;
            tangentX = angle * dot(tangentDirection, tangentA);
            tangentY = angle * dot(tangentDirection, tangentB);
        }
        float alongMajor =
            cosOrientation * tangentX + sinOrientation * tangentY;
        float alongMinor =
            -sinOrientation * tangentX + cosOrientation * tangentY;
        float weight = exp(-0.5f * (
            alongMajor * alongMajor * invMajorSquared +
            alongMinor * alongMinor * invMinorSquared));
        float intensity = pow(10.0f, -0.4f * pointStars[base + 4u]) *
                          weight * brightnessScale;
        color += PointStarColor(starIndex) * intensity;
    }
}

// Exact GPU counterpart of the elliptical
// StarfieldGenerator::AccumulateThroughBeam. The CSR cells use the major axis
// as a conservative circular candidate set; the tangent-plane anisotropic
// Gaussian is the represented star-field decision.
float3 SamplePointStarfield(float3 direction, float sigmaMajor, float sigmaMinor,
                            float orientation, float brightnessScale) {
    static const int kThetaBins = 256;
    static const int kPhiBins = 512;
    float lengthDirection = length(direction);
    if (lengthDirection < 1e-20f || sigmaMajor <= 0.0f || sigmaMinor <= 0.0f)
        return float3(0.0f);
    float3 dir = direction / lengthDirection;
    float3 reference = abs(dir.z) > 0.9f
        ? float3(1.0f, 0.0f, 0.0f)
        : float3(0.0f, 0.0f, 1.0f);
    float3 tangentA = normalize(reference - dot(reference, dir) * dir);
    float3 tangentB = cross(dir, tangentA);
    float theta = acos(clamp(dir.z, -1.0f, 1.0f));
    float phi = atan2(dir.y, dir.x);
    if (phi < 0.0f) phi += TWO_PI;
    float major = max(sigmaMajor, sigmaMinor);
    float minor = min(sigmaMajor, sigmaMinor);
    float cutoff = min(4.0f * major, PI);
    int firstTheta = clamp(int(floor(max(theta - cutoff, 0.0f) / PI * kThetaBins)),
                           0, kThetaBins - 1);
    int lastTheta = clamp(int(floor(min(theta + cutoff, PI) / PI * kThetaBins)),
                          0, kThetaBins - 1);

    float phiHalfWidth = PI;
    float sinTheta = sin(theta);
    if (theta > cutoff && PI - theta > cutoff && sinTheta > 1.0e-8f) {
        phiHalfWidth = asin(clamp(sin(cutoff) / sinTheta, 0.0f, 1.0f));
    }
    float phiBinWidth = TWO_PI / float(kPhiBins);
    int centrePhi = min(int(floor(phi / TWO_PI * kPhiBins)), kPhiBins - 1);
    int phiRadius = min(int(ceil(phiHalfWidth / phiBinWidth)) + 1, kPhiBins);
    bool allPhi = 2 * phiRadius + 1 >= kPhiBins;
    float cosCut = cos(cutoff);
    float invMajorSquared = 1.0f / (major * major);
    float invMinorSquared = 1.0f / (minor * minor);
    float cosOrientation = cos(orientation);
    float sinOrientation = sin(orientation);
    float3 color = float3(0.0f);

    for (int thetaBin = firstTheta; thetaBin <= lastTheta; ++thetaBin) {
        if (allPhi) {
            for (int phiBin = 0; phiBin < kPhiBins; ++phiBin) {
                AccumulatePointStarCell(
                    thetaBin, phiBin, dir, tangentA, tangentB, cosCut,
                    cosOrientation, sinOrientation, invMajorSquared,
                    invMinorSquared, brightnessScale, color);
            }
        } else {
            for (int delta = -phiRadius; delta <= phiRadius; ++delta) {
                int phiBin = (centrePhi + delta) % kPhiBins;
                if (phiBin < 0) phiBin += kPhiBins;
                AccumulatePointStarCell(
                    thetaBin, phiBin, dir, tangentA, tangentB, cosCut,
                    cosOrientation, sinOrientation, invMajorSquared,
                    invMinorSquared, brightnessScale, color);
            }
        }
    }
    return color;
}

// Deterministic procedural-noise authority mirrored from core/disk/turbulence.h. The
// config surface currently exposes the feature toggle only, so these are the
// exact shared defaults rather than a second operator-visible parameter set.
uint TurbulenceHash(uint x, uint y, uint z, uint seed) {
    uint h = seed;
    h ^= x * 374761393u;
    h ^= y * 668265263u;
    h ^= z * 1013904223u;
    h = ((h >> 16u) ^ h) * 2654435761u;
    h = ((h >> 16u) ^ h) * 2654435761u;
    return h;
}

float3 TurbulenceGradient(uint h) {
    static const float3 gradients[12] = {
        float3(1.0f, 1.0f, 0.0f), float3(-1.0f, 1.0f, 0.0f),
        float3(1.0f, -1.0f, 0.0f), float3(-1.0f, -1.0f, 0.0f),
        float3(1.0f, 0.0f, 1.0f), float3(-1.0f, 0.0f, 1.0f),
        float3(1.0f, 0.0f, -1.0f), float3(-1.0f, 0.0f, -1.0f),
        float3(0.0f, 1.0f, 1.0f), float3(0.0f, -1.0f, 1.0f),
        float3(0.0f, 1.0f, -1.0f), float3(0.0f, -1.0f, -1.0f)
    };
    return gradients[h % 12u];
}

float TurbulenceFade(float t) {
    return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
}

float TurbulencePerlin(float3 p, uint seed) {
    int3 cell = int3(floor(p));
    float3 f = p - float3(cell);
    float3 fade = float3(TurbulenceFade(f.x), TurbulenceFade(f.y), TurbulenceFade(f.z));
    float corner[8];
    for (int z = 0; z <= 1; ++z) {
        for (int y = 0; y <= 1; ++y) {
            for (int x = 0; x <= 1; ++x) {
                int slot = z * 4 + y * 2 + x;
                float3 delta = f - float3(float(x), float(y), float(z));
                corner[slot] = dot(
                    TurbulenceGradient(TurbulenceHash(uint(cell.x + x), uint(cell.y + y),
                                                      uint(cell.z + z), seed)),
                    delta);
            }
        }
    }
    float x00 = lerp(corner[0], corner[1], fade.x);
    float x10 = lerp(corner[2], corner[3], fade.x);
    float x01 = lerp(corner[4], corner[5], fade.x);
    float x11 = lerp(corner[6], corner[7], fade.x);
    return lerp(lerp(x00, x10, fade.y), lerp(x01, x11, fade.y), fade.z);
}

float TurbulenceDensityFactor(float3 p) {
    float value = 0.0f;
    float octaveAmplitude = 1.0f;
    float frequency = 0.2f;  // 1 / default outer_scale_M (5 M).
    float totalAmplitude = 0.0f;
    for (uint octave = 0u; octave < 6u; ++octave) {
        value += octaveAmplitude *
                 TurbulencePerlin(p * frequency, 12345u + octave);
        totalAmplitude += octaveAmplitude;
        octaveAmplitude *= 0.5f;
        frequency *= 2.0f;
    }
    return max(1.0f + 0.3f * value / totalAmplitude, 0.01f);
}

float VolumeScaleHeight(float r, float innerRadius, float hOverR, float hPower) {
    if (r <= 0.0f || innerRadius <= 0.0f) return 0.0f;
    float localHOverR = clamp(hOverR * pow(r / innerRadius, hPower), 0.01f, 0.5f);
    return localHOverR * r;
}

bool IsInsideVolumeDisk(float cylindricalR, float z, float innerRadius,
                        float outerRadius, float hOverR, float hPower) {
    if (cylindricalR < innerRadius || cylindricalR > outerRadius) return false;
    float H = VolumeScaleHeight(cylindricalR, innerRadius, hOverR, hPower);
    return abs(z) <= 3.0f * H;
}

float VolumeDiskTemperature(float cylindricalR, float z, float innerRadius,
                            float innerTemperature, int temperatureModel,
                            float rIsco, float M, float aStar,
                            float hOverR, float hPower) {
    float midplane = temperatureModel == 1
        ? ShakuraSunyaevTemperature(innerTemperature, cylindricalR, innerRadius)
        : FullPageThorneTemperature(innerTemperature, cylindricalR, rIsco, M, aStar);
    // Vertically isothermal closure. z, hOverR, and hPower define membership
    // and opacity only; no synthetic atmosphere-temperature law is implied.
    return max(midplane, 0.0f);
}

// Cubic Hermite sample of one accepted central-ray segment. The position and
// affine tangent match both endpoints, unlike a chord interpolation through a
// curved geodesic.
Vec4Cart AcceptedSegmentPosition(Vec4Cart start, Vec4Cart startVelocity,
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

Vec4Cart AcceptedSegmentTangent(Vec4Cart start, Vec4Cart startVelocity,
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

Real EvaluateCubic(Real coefficientA, Real coefficientB, Real coefficientC,
                   Real coefficientD, Real fraction) {
    return ((coefficientA * fraction + coefficientB) * fraction + coefficientC) *
               fraction +
           coefficientD;
}

// Integrate one geodesic segment using the same midpoint transfer recurrence as
// GeodesicTracer::AccumulateVolumetricEmission.
void AccumulateVolumeSegment(Vec4Cart start, Vec4Cart finish,
                             Vec4Cart startVelocity, Vec4Cart finishVelocity,
                             Real affineLength, int sampleCount,
                             bool volumeEnabled, bool turbulenceEnabled,
                             bool dopplerBeaming, float innerRadius,
                             float outerRadius, float innerTemperature,
                             int temperatureModel, float rIsco, float M,
                             float a, float aStar, float hOverR, float hPower,
                             float tauMidplane, float observerFrequency, inout float opticalDepth,
                             inout float3 emission) {
    if (affineLength <= 0.0f || opticalDepth >= 10.0f) return;
    float deltaLambda = float(affineLength) / float(sampleCount);
    for (int sample = 0; sample < sampleCount; ++sample) {
        Real fraction = Real((float(sample) + 0.5f) / float(sampleCount));
        Vec4Cart position = AcceptedSegmentPosition(
            start, startVelocity, finish, finishVelocity, affineLength, fraction);
        Vec4Cart photonVelocity = AcceptedSegmentTangent(
            start, startVelocity, finish, finishVelocity, affineLength, fraction);
        float3 point = float3(float(position.x), float(position.y), float(position.z));
        float diskR =
            float(ComputeKerrRadiusCart(point.x, point.y, point.z, Real(a)));

        Real localMetric[4][4];
        Real localInverse[4][4];
        GetKerrSchildFamilyMetric(position, M, a, 0.0f, 0.0f, localMetric, localInverse);
        Real lapse = Real(1.0f) / sqrt(-localInverse[0][0]);
        Vec4Cart eulerian = MakeVec4Cart(-lapse * localInverse[0][0],
                                         -lapse * localInverse[1][0],
                                         -lapse * localInverse[2][0],
                                         -lapse * localInverse[3][0]);
        float eulerianFrequency = float(MetricDotCart(localMetric, photonVelocity, eulerian));

        float diskDeltaTau = 0.0f;
        float3 diskSource = float3(0.0f);

        if (volumeEnabled &&
            IsInsideVolumeDisk(diskR, point.z, innerRadius, outerRadius,
                               hOverR, hPower)) {
            float H = VolumeScaleHeight(diskR, innerRadius, hOverR, hPower);
            float zOverH = point.z / H;
            float diskOpacity = tauMidplane / (sqrt(TWO_PI) * H) *
                          pow(diskR / innerRadius, -1.5f) *
                          exp(-0.5f * zOverH * zOverH);
            if (turbulenceEnabled) {
                // Match the CPU authority: its procedural field is defined in
                // (Kerr radial coordinate, azimuth, z), not the oblate
                // Kerr-Schild Cartesian transverse radius.
                float diskPhi = atan2(point.y, point.x);
                float3 turbulencePoint =
                    float3(diskR * cos(diskPhi), diskR * sin(diskPhi), point.z);
                diskOpacity *= TurbulenceDensityFactor(turbulencePoint);
            }
            float temperature =
                VolumeDiskTemperature(diskR, point.z, innerRadius,
                                      innerTemperature, temperatureModel, rIsco,
                                      M, aStar, hOverR, hPower);
            if (isfinite(diskOpacity) && diskOpacity > 0.0f &&
                isfinite(temperature) && temperature > 0.0f) {
                float sqrtM = sqrt(M);
                float omega = sqrtM / (pow(diskR, 1.5f) + a * sqrtM);
                Vec4Cart emitter = MakeVec4Cart(Real(1.0f), Real(-omega * point.y),
                                                Real(omega * point.x), Real(0.0f));
                Real emitterNorm = MetricDotCart(localMetric, emitter, emitter);
                if (!isfinite(emitterNorm) || !(emitterNorm < Real(0.0f))) continue;
                emitter = ScaleVec4Cart(emitter, Real(1.0f) / sqrt(-emitterNorm));
                float emittedFrequency =
                    float(MetricDotCart(localMetric, photonVelocity, emitter));
                if (!isfinite(emittedFrequency) || !(emittedFrequency > 0.0f) ||
                    !isfinite(eulerianFrequency) || !(eulerianFrequency > 0.0f)) {
                    continue;
                }
                diskDeltaTau = diskOpacity * emittedFrequency * deltaLambda;
                float emittedSource = pow(temperature / innerTemperature, 4.0f);
                float sourceFrequency = dopplerBeaming ? emittedFrequency : eulerianFrequency;
                float gFactor = observerFrequency / sourceFrequency;
                float observedSource = ObservedDiskBolometricSource(emittedSource, gFactor);
                diskSource = BlackbodyColor(temperature * gFactor) * observedSource;
            }
        }

        float uncappedDeltaTau = diskDeltaTau;
        if (uncappedDeltaTau <= 0.0f) continue;
        float deltaTau = min(uncappedDeltaTau, 10.0f - opticalDepth);
        float3 source = diskSource;
        float transmission = exp(-deltaTau);
        emission += exp(-opticalDepth) * source * (1.0f - transmission);
        opticalDepth += deltaTau;
        if (opticalDepth >= 10.0f) break;
    }
}

[shader("compute")]
[numthreads(8, 8, 1)]
void ComputeMain(uint3 tid: SV_DispatchThreadID) {
    int imageWidth = int(params[0]);
    int imageHeight = int(params[1]);
    uint metricId = uint(params[2]);

    int tileOriginX = int(params[31]);
    int tileOriginY = int(params[32]);
    int tileWidth = int(params[33]);
    int tileHeight = int(params[34]);

    int tx = int(tid.x);
    int ty = int(tid.y);
    if (tx >= tileWidth || ty >= tileHeight) return;

    int ix = tileOriginX + tx;
    int iy = tileOriginY + ty;
    if (ix >= imageWidth || iy >= imageHeight) return;

    float M = params[3];
    float a = params[4];
    float Q = params[5];
    float Lambda = params[6];

    float3 camPos = P3(7);
    float3 camDir = P3(10);
    float3 camRight = P3(13);
    float3 camUp = P3(16);
    float fov = params[19];
    float aspect = params[20];

    int maxSteps = int(params[21]);
    float stepScale = params[22];
    float minStep = params[23];
    float maxStep = params[24];
    float escapeRadius = params[25];
    float captureFactor = params[26];

    bool diskEnabled = params[27] > 0.5f;
    float innerRadius = params[28];
    float outerRadius = params[29];
    float innerTemperature = params[30];

    bool starfieldEnabled = params[35] > 0.5f;
    int starfieldWidth = int(params[36]);
    int starfieldHeight = int(params[37]);

    float warpVs = params[38];
    float warpSigma = params[39];
    float warpR = params[40];

    // Morris-Thorne (dispatch id 3): throat radius and shape selector
    // (0 Ellis, 1 zero-tidal, 2 absurdly benign); Phi0 = 0 is the zero-tidal
    // subclass the CPU authority carries.
    float throatB0 = params[42];
    int wormShape = int(params[43]);
    float sampleU = params[44];
    float sampleV = params[45];
    float sampleIndex = params[46];
    float3 cameraBeta = P3(47);
    int temperatureModel = int(params[50]);  // 0 Page/Novikov-Thorne, 1 Shakura-Sunyaev.
    bool dopplerBeaming = params[51] > 0.5f;
    int lensModel = int(params[52]);  // 0 pinhole, 1 deterministic thin lens.
    float focalLength = params[53];
    float aperture = params[54];
    float focusDistance = params[55];
    bool pointStarfieldEnabled = params[56] > 0.5f;
    float pixelAngularSize = params[57];
    float pointBrightnessScale = params[58];
    bool volumeEnabled = params[59] > 0.5f;
    float volumeHOverR = params[60];
    float volumeHPower = params[61];
    float volumeTauMidplane = params[62];
    int volumeSamples = int(params[63]);
    bool turbulenceEnabled = params[64] > 0.5f;

    // Ray bundles (P2): propagate two geodesic-deviation vectors through the render
    // loop and report the beam's transverse ellipse. Default off (params[41]==0)
    // leaves the RGBA output untouched, so the parity and pinned renders do not
    // move. Uses gr_deviation's full-Riemann acceleration (DNGR eq A.19) in the
    // same Cartesian chart the integrator runs in.
    bool beamsEnabled = params[41] > 0.5f;

    // Camera ray in world Cartesian. NDC follows the CPU pinhole camera geometry
    // (the renderer builds camDir/camRight/camUp so the two paths trace the same
    // scene): u spans [-aspect, aspect], v spans [-1, 1] top-to-bottom.
    float imageX = (2.0f * (float(ix) + sampleU) / float(imageWidth) - 1.0f) * aspect;
    float imageY = 1.0f - 2.0f * (float(iy) + sampleV) / float(imageHeight);
    float tanHalfFov = tan(fov * 0.5f);
    float3 localDir;
    if (lensModel == 1) {
        // The centre-aperture projection is exactly the pinhole projection.
        // 50 mm-equivalent defines one virtual lens unit because the geometric
        // spacetime has no physical M-to-millimetre scale (core/camera.h).
        float focusX = imageX * tanHalfFov * focusDistance;
        float focusY = imageY * tanHalfFov * focusDistance;
        float apertureRadius = (focalLength / 50.0f) / (2.0f * aperture);
        float lensTheta = TWO_PI * sampleU;
        float lensRadius = apertureRadius * sqrt(sampleV);
        float lensX = lensRadius * cos(lensTheta);
        float lensY = lensRadius * sin(lensTheta);
        localDir = normalize(float3(-focusDistance, -(focusY - lensY), focusX - lensX));
    } else {
        localDir = normalize(float3(-1.0f, -imageY * tanHalfFov, imageX * tanHalfFov));
    }
    float3 screenDirection = float3(-localDir.x, -localDir.y, localDir.z);

    bool isWarp = (metricId == kDispatchWarpDrive);
    bool isWormhole = (metricId == kDispatchMorrisThorne);

    // Initial 4-velocity from the null condition at the camera position, using
    // the dispatched metric.
    Real g[4][4];
    Real g_inv[4][4];
    Vec4Cart pos0 = MakeVec4Cart(0.0f, camPos.x, camPos.y, camPos.z);
    if (isWarp) {
        GetWarpDriveMetric(pos0, warpVs, warpSigma, warpR, 0.0f, 0.0f, 0.0f, g, g_inv);
    } else if (isWormhole) {
        GetMorrisThorneCartesianMetric(pos0, throatB0, 0.0f, wormShape, g, g_inv);
    } else {
        GetKerrSchildFamilyMetric(pos0, M, a, Q, Lambda, g, g_inv);
    }
    // Capture radius: the buffered Kerr-Schild outer horizon r_+ = M +
    // sqrt(M^2 - a^2 - Q^2) for the black-hole family; the throat radius for
    // the one-sheet wormhole (its far mouth is unmodelled, so rays reaching
    // the throat terminate as captured, matching the CPU path). The warp
    // bubble is horizonless: capture is disabled.
    float disc = M * M - a * a - Q * Q;
    float horizonRadius = (disc > 0.0f) ? (M + sqrt(disc)) : M;
    if (isWormhole) horizonRadius = throatB0;
    horizonRadius *= captureFactor;

    GeodesicStateCart state;
    state.x = pos0;
    CameraObserverFrame cameraFrame =
        BuildCameraObserverFrame(g, g_inv, camDir, camUp, camRight, cameraBeta);
    state.u = PastDirectedCameraRay(cameraFrame, screenDirection);
    float observerFrequency = 1.0f;

    float3 color = float3(0.0f, 0.0f, 0.0f);
    float transmittance = 1.0f;
    float volumeOpticalDepth = 0.0f;
    float3 volumeEmission = float3(0.0f);
    bool escaped = false;
    bool captured = false;

    // Two deviation vectors span the ray bundle (P2). A point-source/pupil bundle
    // is used when the point-star field consumes the ellipse; the diagnostic beam
    // path otherwise retains the parallel-bundle initial condition.
    Vec4Cart launchScreenA;
    Vec4Cart launchScreenB;
    SachsScreenBasis(cameraFrame, screenDirection, launchScreenA, launchScreenB);
    const float kBeamSeed = 1.0e-3f;
    bool pupilBundle = pointStarfieldEnabled;
    Real4 xiA = pupilBundle
        ? Real4(0.0f)
        : Real4(launchScreenA.t, launchScreenA.x, launchScreenA.y, launchScreenA.z) *
              Real(kBeamSeed);
    Real4 xiB = pupilBundle
        ? Real4(0.0f)
        : Real4(launchScreenB.t, launchScreenB.x, launchScreenB.y, launchScreenB.z) *
              Real(kBeamSeed);
    Real4 xiDotA = pupilBundle
        ? Real4(launchScreenA.t, launchScreenA.x, launchScreenA.y, launchScreenA.z) *
              Real(kBeamSeed)
        : Real4(0.0f);
    Real4 xiDotB = pupilBundle
        ? Real4(launchScreenB.t, launchScreenB.x, launchScreenB.y, launchScreenB.z) *
              Real(kBeamSeed)
        : Real4(0.0f);

#ifdef SIRIUS_FP32_COMP
    // Compensated rung: running Kahan compensation for the state accumulation
    // across steps, recovering the low-order bits plain fp32 addition drops.
    // The RK4 increment itself is unchanged; only the commit differs.
    Vec4Cart compX = MakeVec4Cart(0.0f, 0.0f, 0.0f, 0.0f);
    Vec4Cart compU = MakeVec4Cart(0.0f, 0.0f, 0.0f, 0.0f);
#endif

    Real adaptiveStep = Real(maxStep);

    for (int step = 0; step < maxSteps; step++) {
        Real px = state.x.x;
        Real py = state.x.y;
        Real pz = state.x.z;
        Real rBig = sqrt(px * px + py * py + pz * pz);

        if (isWormhole) {
            if (rBig <= horizonRadius) {
                captured = true;
                break;
            }  // reached the throat: dark
        } else if (!isWarp) {
            Real rKerr = ComputeKerrRadiusCart(px, py, pz, Real(a));
            if (rKerr <= horizonRadius) {
                captured = true;
                break;
            }  // captured: shadow
        }
        if (rBig >= escapeRadius) { escaped = true; break; }

        // Preserve scale covariance for the exotic metrics. Alcubierre's sigma
        // is an inverse length, so a diffuse wall is governed by 1/sigma while
        // a compact wall is governed by R; the wormhole is governed by b0.
        Real exoticScale = isWarp
            ? max(Real(warpR), Real(1.0f / warpSigma))
            : Real(throatB0);
        Real rStep = (isWarp || isWormhole)
            ? max(rBig, exoticScale)
            : ComputeKerrRadiusCart(px, py, pz, Real(a));
        Real geometricStep = clamp(stepScale * rStep, Real(minStep), Real(maxStep));
        Real dt = min(geometricStep, adaptiveStep);

        AdaptiveGeodesicStep attempted = AttemptAdaptiveGeodesicStep(
            state, dt, isWarp, isWormhole, M, a, Q, Lambda,
            warpVs, warpSigma, warpR, throatB0, wormShape,
            AdaptiveIntegrationTolerance(), AdaptiveNullTolerance());
        if (!attempted.accepted) {
            if (dt <= Real(minStep) * Real(1.0001f)) break;
            adaptiveStep = max(Real(minStep), dt * attempted.stepFactor);
            continue;
        }
        adaptiveStep = clamp(dt * attempted.stepFactor, Real(minStep), Real(maxStep));

        GeodesicStateCart next = attempted.candidate;
        bool terminalDiskHit = false;
#ifdef SIRIUS_FP32_COMP
        Vec4Cart acceptedDX = AddVec4Cart(next.x, ScaleVec4Cart(state.x, Real(-1.0f)));
        Vec4Cart acceptedDU = AddVec4Cart(next.u, ScaleVec4Cart(state.u, Real(-1.0f)));
        next.x = KahanAddVec4Cart(state.x, compX, acceptedDX);
        next.u = KahanAddVec4Cart(state.u, compU, acceptedDU);
#endif

        // Optional volumetric disk (Kerr-Schild black holes only).
        if (volumeEnabled && !isWarp && !isWormhole) {
            float aStar = (M > 1e-6f) ? (a / M) : a;
            float rIsco = ComputeKerrISCO(M, aStar);
            AccumulateVolumeSegment(
                state.x, next.x, state.u, next.u, dt, volumeSamples,
                volumeEnabled, turbulenceEnabled, dopplerBeaming,
                innerRadius, outerRadius, innerTemperature, temperatureModel,
                rIsco, M, a, aStar, volumeHOverR, volumeHPower,
                volumeTauMidplane, observerFrequency,
                volumeOpticalDepth, volumeEmission);
        }

        // Optional thin disk (Kerr-Schild black holes only): emission when the
        // segment crosses the equatorial plane z=0 within [innerRadius,
        // outerRadius], with Novikov-Thorne temperature, blackbody colour, and
        // g^4 relativistic beaming (gr_disk authorities).
        if (diskEnabled && !volumeEnabled && !isWarp && !isWormhole &&
            transmittance > 0.01f) {
            Real dzStart = dt * state.u.z;
            Real dzEnd = dt * next.u.z;
            Real coefficientA = Real(2.0f) * state.x.z - Real(2.0f) * next.x.z +
                                dzStart + dzEnd;
            Real coefficientB = Real(-3.0f) * state.x.z + Real(3.0f) * next.x.z -
                                Real(2.0f) * dzStart - dzEnd;
            Real coefficientC = dzStart;
            Real coefficientD = state.x.z;
            Real coefficientScale = max(
                Real(1.0f), max(max(abs(coefficientA), abs(coefficientB)),
                                max(abs(coefficientC), abs(coefficientD))));
#ifdef SIRIUS_FP64
            Real rootTolerance = Real(1.5e-14) * coefficientScale;
            Real fractionTolerance = Real(2.0e-12);
#else
            Real rootTolerance = Real(8.0e-6f) * coefficientScale;
            Real fractionTolerance = Real(2.0e-5f);
#endif
            Real rootFractions[4];
            rootFractions[0] = Real(0.0f);
            rootFractions[1] = Real(0.0f);
            rootFractions[2] = Real(0.0f);
            rootFractions[3] = Real(0.0f);
            int rootCount = 0;
            bool coplanar = abs(coefficientA) <= rootTolerance &&
                            abs(coefficientB) <= rootTolerance &&
                            abs(coefficientC) <= rootTolerance &&
                            abs(coefficientD) <= rootTolerance;
            if (coplanar) {
                // The governed device step is at most 2 M, below the supported
                // ISCO-to-20 M annulus width. Bracket its first radial entry.
                float startRadius = float(ComputeKerrRadiusCart(
                    state.x.x, state.x.y, Real(0.0f), Real(a)));
                if (startRadius >= innerRadius && startRadius <= outerRadius) {
                    rootFractions[rootCount++] = Real(0.0f);
                } else {
                    Real radialLower = Real(0.0f);
                    for (int radialInterval = 1; radialInterval <= 64 && rootCount == 0;
                         ++radialInterval) {
                        Real radialUpper = Real(float(radialInterval) / 64.0f);
                        Vec4Cart upperPosition = AcceptedSegmentPosition(
                            state.x, state.u, next.x, next.u, dt, radialUpper);
                        float upperRadius = float(ComputeKerrRadiusCart(
                            upperPosition.x, upperPosition.y, Real(0.0f), Real(a)));
                        bool upperInside = upperRadius >= innerRadius && upperRadius <= outerRadius;
                        if (upperInside) {
                            for (int radialIteration = 0; radialIteration < 40;
                                 ++radialIteration) {
                                Real radialMidpoint =
                                    (radialLower + radialUpper) * Real(0.5f);
                                Vec4Cart midpointPosition = AcceptedSegmentPosition(
                                    state.x, state.u, next.x, next.u, dt, radialMidpoint);
                                float midpointRadius = float(ComputeKerrRadiusCart(
                                    midpointPosition.x, midpointPosition.y, Real(0.0f), Real(a)));
                                bool midpointInside = midpointRadius >= innerRadius &&
                                                      midpointRadius <= outerRadius;
                                if (midpointInside) {
                                    radialUpper = radialMidpoint;
                                } else {
                                    radialLower = radialMidpoint;
                                }
                            }
                            rootFractions[rootCount++] = radialUpper;
                        }
                        radialLower = radialUpper;
                    }
                }
            } else {
                // The Hermite z(s) is cubic. Split it at the roots of z'(s),
                // making every interval monotone, then retain all plane roots
                // in observer-to-scene order. This includes tangent contacts
                // and two crossings with equal-sign segment endpoints.
                Real breakpoints[4];
                breakpoints[0] = Real(0.0f);
                breakpoints[1] = Real(1.0f);
                int breakpointCount = 2;
                Real derivativeDiscriminant =
                    coefficientB * coefficientB -
                    Real(3.0f) * coefficientA * coefficientC;
                Real discriminantScale = max(
                    Real(1.0f), max(coefficientB * coefficientB,
                                    abs(Real(3.0f) * coefficientA * coefficientC)));
#ifdef SIRIUS_FP64
                Real discriminantTolerance = Real(1.5e-14) * discriminantScale;
#else
                Real discriminantTolerance = Real(8.0e-6f) * discriminantScale;
#endif
                if (abs(coefficientA) > rootTolerance &&
                    derivativeDiscriminant >= -discriminantTolerance) {
                    Real derivativeRoot = sqrt(max(derivativeDiscriminant, Real(0.0f)));
                    Real firstCritical = (-coefficientB - derivativeRoot) /
                                         (Real(3.0f) * coefficientA);
                    Real secondCritical = (-coefficientB + derivativeRoot) /
                                          (Real(3.0f) * coefficientA);
                    if (firstCritical > Real(0.0f) && firstCritical < Real(1.0f))
                        breakpoints[breakpointCount++] = firstCritical;
                    if (secondCritical > Real(0.0f) && secondCritical < Real(1.0f))
                        breakpoints[breakpointCount++] = secondCritical;
                } else if (abs(coefficientB) > rootTolerance) {
                    Real critical = -coefficientC / (Real(2.0f) * coefficientB);
                    if (critical > Real(0.0f) && critical < Real(1.0f))
                        breakpoints[breakpointCount++] = critical;
                }
                for (int left = 0; left < breakpointCount; ++left) {
                    for (int right = left + 1; right < breakpointCount; ++right) {
                        if (breakpoints[right] < breakpoints[left]) {
                            Real swapValue = breakpoints[left];
                            breakpoints[left] = breakpoints[right];
                            breakpoints[right] = swapValue;
                        }
                    }
                }
                for (int rootInterval = 0; rootInterval + 1 < breakpointCount;
                     ++rootInterval) {
                    Real lower = breakpoints[rootInterval];
                    Real upper = breakpoints[rootInterval + 1];
                    Real lowerZ = EvaluateCubic(coefficientA, coefficientB, coefficientC,
                                                coefficientD, lower);
                    Real upperZ = EvaluateCubic(coefficientA, coefficientB, coefficientC,
                                                coefficientD, upper);
                    if (abs(lowerZ) <= rootTolerance && rootCount < 4 &&
                        (rootCount == 0 ||
                         abs(lower - rootFractions[rootCount - 1]) >
                             fractionTolerance)) {
                        rootFractions[rootCount++] = lower;
                    }
                    if ((lowerZ < Real(0.0f)) != (upperZ < Real(0.0f))) {
                        Real rootLower = lower;
                        Real rootUpper = upper;
                        Real rootLowerZ = lowerZ;
                        for (int rootIteration = 0; rootIteration < 40; ++rootIteration) {
                            Real midpoint = (rootLower + rootUpper) * Real(0.5f);
                            Real midpointZ = EvaluateCubic(
                                coefficientA, coefficientB, coefficientC, coefficientD,
                                midpoint);
                            if ((rootLowerZ < Real(0.0f)) ==
                                (midpointZ < Real(0.0f))) {
                                rootLower = midpoint;
                                rootLowerZ = midpointZ;
                            } else {
                                rootUpper = midpoint;
                            }
                        }
                        Real rootFraction = (rootLower + rootUpper) * Real(0.5f);
                        if (rootCount < 4 &&
                            (rootCount == 0 ||
                             abs(rootFraction - rootFractions[rootCount - 1]) >
                                 fractionTolerance)) {
                            rootFractions[rootCount++] = rootFraction;
                        }
                    }
                }
                Real finalZ = EvaluateCubic(coefficientA, coefficientB, coefficientC,
                                            coefficientD, Real(1.0f));
                if (abs(finalZ) <= rootTolerance && rootCount < 4 &&
                    (rootCount == 0 ||
                     abs(Real(1.0f) - rootFractions[rootCount - 1]) >
                         fractionTolerance)) {
                    rootFractions[rootCount++] = Real(1.0f);
                }
            }

            bool foundDiskCrossing = false;
            Real frac = Real(0.0f);
            Vec4Cart crossingPosition;
            Vec4Cart photonVelocity;
            float cr = 0.0f;
            for (int rootIndex = 0; rootIndex < rootCount && !foundDiskCrossing;
                 ++rootIndex) {
                Real candidateFraction = rootFractions[rootIndex];
                Vec4Cart candidatePosition = AcceptedSegmentPosition(
                    state.x, state.u, next.x, next.u, dt, candidateFraction);
                float candidateRadius = float(ComputeKerrRadiusCart(
                    candidatePosition.x, candidatePosition.y, Real(0.0f), Real(a)));
                if (candidateRadius >= innerRadius && candidateRadius <= outerRadius) {
                    foundDiskCrossing = true;
                    frac = candidateFraction;
                    crossingPosition = candidatePosition;
                    photonVelocity = AcceptedSegmentTangent(
                        state.x, state.u, next.x, next.u, dt, candidateFraction);
                    cr = candidateRadius;
                }
            }
            if (foundDiskCrossing) {
                    // Shading is float (gr_disk carries trig/pow, 32-bit-only in
                    // SPIR-V); the crossing point itself was found at Real.
                    float cx = float(crossingPosition.x);
                    float cy = float(crossingPosition.y);
                    float aStar = (M > 1e-6f) ? (a / M) : a;
                    float rIsco = ComputeKerrISCO(M, aStar);
                    float T;
                    if (temperatureModel == 1) {
                        T = ShakuraSunyaevTemperature(innerTemperature, cr, innerRadius);
                    } else {
                        T = FullPageThorneTemperature(innerTemperature, cr, rIsco, M, aStar);
                    }
                    float cphi = atan2(cy, cx);
                    float2 transfer = ComputeKerrDiskTransferCart(
                        cr, cphi, photonVelocity, observerFrequency, M, a);
                    float gFactor = dopplerBeaming ? transfer.x : transfer.y;
                    float emittedSource = pow(T / innerTemperature, 4.0f);
                    float observedScale =
                        ObservedDiskBolometricSource(emittedSource, gFactor);
                    float3 emit = BlackbodyColor(T * gFactor) * observedScale;
                    // Novikov-Thorne is optically thick: the first surface met
                    // along the past ray hides all farther crossings.
                    color = emit;
                    transmittance = 0.0f;
                    next.x = crossingPosition;
                    next.u = photonVelocity;
                    dt *= frac;
                    terminalDiskHit = true;
            }
        }
        // RK4 on the covariant first-order Jacobi system. Connection and
        // curvature are evaluated at the start, Hermite midpoint, and end of
        // the accepted central-ray step; freezing them at the start lowers the
        // bundle update to first order in a varying tidal field.
        if (beamsEnabled && !isWarp && !isWormhole) {
            Real beamRiemann[4][4][4][4];
            Real beamGamma[4][4][4];
            Real4 k1XiA, k1VA, k2XiA, k2VA, k3XiA, k3VA, k4XiA, k4VA;
            Real4 k1XiB, k1VB, k2XiB, k2VB, k3XiB, k3VB, k4XiB, k4VB;

            GetRiemannTensorCart(state.x, M, a, Q, Lambda, beamRiemann);
            GetKerrSchildFamilyChristoffel(state.x, M, a, Q, Lambda, beamGamma);
            GeodesicDeviationRhs(beamGamma, beamRiemann, state.u, xiA, xiDotA, k1XiA, k1VA);
            GeodesicDeviationRhs(beamGamma, beamRiemann, state.u, xiB, xiDotB, k1XiB, k1VB);

            // Cubic Hermite midpoint of the accepted central ray. It matches
            // both endpoint events and tangents and supplies the same stage
            // geometry used by the CPU bundle path.
            Vec4Cart beamMidPosition = AddVec4Cart(
                ScaleVec4Cart(AddVec4Cart(state.x, next.x), Real(0.5f)),
                ScaleVec4Cart(AddVec4Cart(state.u, ScaleVec4Cart(next.u, Real(-1.0f))),
                              dt / Real(8.0f)));
            Vec4Cart beamMidTangent = AddVec4Cart(
                ScaleVec4Cart(AddVec4Cart(next.x, ScaleVec4Cart(state.x, Real(-1.0f))),
                              Real(1.5f) / dt),
                ScaleVec4Cart(AddVec4Cart(state.u, next.u), Real(-0.25f)));
            GetRiemannTensorCart(beamMidPosition, M, a, Q, Lambda, beamRiemann);
            GetKerrSchildFamilyChristoffel(beamMidPosition, M, a, Q, Lambda, beamGamma);

            GeodesicDeviationRhs(beamGamma, beamRiemann, beamMidTangent,
                                 xiA + k1XiA * (Real(0.5f) * dt),
                                 xiDotA + k1VA * (Real(0.5f) * dt), k2XiA, k2VA);
            GeodesicDeviationRhs(beamGamma, beamRiemann, beamMidTangent,
                                 xiA + k2XiA * (Real(0.5f) * dt),
                                 xiDotA + k2VA * (Real(0.5f) * dt), k3XiA, k3VA);
            GeodesicDeviationRhs(beamGamma, beamRiemann, beamMidTangent,
                                 xiB + k1XiB * (Real(0.5f) * dt),
                                 xiDotB + k1VB * (Real(0.5f) * dt), k2XiB, k2VB);
            GeodesicDeviationRhs(beamGamma, beamRiemann, beamMidTangent,
                                 xiB + k2XiB * (Real(0.5f) * dt),
                                 xiDotB + k2VB * (Real(0.5f) * dt), k3XiB, k3VB);

            GetRiemannTensorCart(next.x, M, a, Q, Lambda, beamRiemann);
            GetKerrSchildFamilyChristoffel(next.x, M, a, Q, Lambda, beamGamma);
            GeodesicDeviationRhs(beamGamma, beamRiemann, next.u,
                                 xiA + k3XiA * dt, xiDotA + k3VA * dt, k4XiA, k4VA);
            GeodesicDeviationRhs(beamGamma, beamRiemann, next.u,
                                 xiB + k3XiB * dt, xiDotB + k3VB * dt, k4XiB, k4VB);

            const Real sixthStep = dt / Real(6.0f);
            xiA += (k1XiA + k2XiA * Real(2.0f) + k3XiA * Real(2.0f) + k4XiA) * sixthStep;
            xiDotA += (k1VA + k2VA * Real(2.0f) + k3VA * Real(2.0f) + k4VA) * sixthStep;
            xiB += (k1XiB + k2XiB * Real(2.0f) + k3XiB * Real(2.0f) + k4XiB) * sixthStep;
            xiDotB += (k1VB + k2VB * Real(2.0f) + k3VB * Real(2.0f) + k4VB) * sixthStep;
        }

        state = next;
        if (terminalDiskHit) break;

        // Classify the accepted endpoint in this iteration. Delaying this until
        // the next loop made a final-attempt capture look like an ordinary
        // max-attempt endpoint and allowed terminal-frame construction inside
        // the capture surface. Segment disk events retain observer-first
        // precedence because they were processed above.
        Real nextRadius = sqrt(state.x.x * state.x.x + state.x.y * state.x.y +
                               state.x.z * state.x.z);
        if (isWormhole) {
            if (nextRadius <= horizonRadius) {
                captured = true;
                break;
            }
        } else if (!isWarp) {
            Real nextKerrRadius =
                ComputeKerrRadiusCart(state.x.x, state.x.y, state.x.z, Real(a));
            if (nextKerrRadius <= horizonRadius) {
                captured = true;
                break;
            }
        }
        if (nextRadius >= escapeRadius) {
            escaped = true;
            break;
        }
    }

    float3 terminalDirection = float3(0.0f, 0.0f, 1.0f);
    float beamOrientation = 0.0f;
    float beamMajorExpansion = 1.0f;
    float beamMinorExpansion = 1.0f;
    float beamExpansion = 1.0f;
    // A captured endpoint has no physical stationary Eulerian observer frame.
    // Escaped rays require the sky direction; other non-captured endpoints need
    // terminal geometry only when beam diagnostics were requested.
    if (escaped || (beamsEnabled && !captured)) {
        Real terminalMetric[4][4];
        Real terminalInverse[4][4];
        if (isWarp) {
            GetWarpDriveMetric(state.x, warpVs, warpSigma, warpR, 0.0f, 0.0f, 0.0f,
                               terminalMetric, terminalInverse);
        } else if (isWormhole) {
            GetMorrisThorneCartesianMetric(state.x, throatB0, 0.0f, wormShape,
                                           terminalMetric, terminalInverse);
        } else {
            GetKerrSchildFamilyMetric(state.x, M, a, Q, Lambda, terminalMetric, terminalInverse);
        }
        CameraObserverFrame terminalFrame = BuildCameraObserverFrame(
            terminalMetric, terminalInverse, float3(1.0f, 0.0f, 0.0f),
            float3(0.0f, -1.0f, 0.0f), float3(0.0f, 0.0f, 1.0f), float3(0.0f));
        Real terminalFrequency = MetricDotCart(terminalMetric, state.u, terminalFrame.time);
        terminalDirection = normalize(float3(
            float(MetricDotCart(terminalMetric, state.u, terminalFrame.forward) /
                  terminalFrequency),
            float(MetricDotCart(terminalMetric, state.u, terminalFrame.up) / terminalFrequency),
            float(MetricDotCart(terminalMetric, state.u, terminalFrame.right) /
                  terminalFrequency)));
        Vec4Cart terminalScreenA;
        Vec4Cart terminalScreenB;
        SachsScreenBasis(terminalFrame, terminalDirection, terminalScreenA, terminalScreenB);
        Vec4Cart xiACart = MakeVec4Cart(xiA.x, xiA.y, xiA.z, xiA.w);
        Vec4Cart xiBCart = MakeVec4Cart(xiB.x, xiB.y, xiB.z, xiB.w);
        float beamA = float(MetricDotCart(terminalMetric, xiACart, terminalScreenA));
        float beamB = float(MetricDotCart(terminalMetric, xiBCart, terminalScreenA));
        float beamC = float(MetricDotCart(terminalMetric, xiACart, terminalScreenB));
        float beamD = float(MetricDotCart(terminalMetric, xiBCart, terminalScreenB));
        float beamSemiMajor;
        float beamSemiMinor;
        float beamDeterminant;
        BeamEllipseFromTransverseMap(beamA, beamB, beamC, beamD,
                                     beamSemiMajor, beamSemiMinor, beamOrientation,
                                     beamDeterminant);
        float terminalRadius =
            length(float3(float(state.x.x), float(state.x.y), float(state.x.z)));
        float beamNormalization = pupilBundle
            ? max(terminalRadius * kBeamSeed, 1.0e-12f)
            : kBeamSeed;
        beamMajorExpansion = beamSemiMajor / beamNormalization;
        beamMinorExpansion = beamSemiMinor / beamNormalization;
        beamExpansion = sqrt(max(beamMajorExpansion * beamMinorExpansion, 0.0f));
    }

    // Captured (shadow) pixels keep only whatever disk emission accumulated
    // (black when none), exactly as the CPU horizon outcome does.
    float3 result = color;
    if (escaped) {
        float3 dir = terminalDirection;
        float sigmaMajor = beamsEnabled
            ? max(pixelAngularSize * beamMajorExpansion, pixelAngularSize)
            : 0.3f * pixelAngularSize;
        float sigmaMinor = beamsEnabled
            ? max(pixelAngularSize * beamMinorExpansion, pixelAngularSize)
            : 0.3f * pixelAngularSize;
        float3 bg = pointStarfieldEnabled
            ? SamplePointStarfield(dir, sigmaMajor, sigmaMinor, beamOrientation,
                                   pointBrightnessScale)
            : (starfieldEnabled
                ? SampleStarfield(dir, starfieldWidth, starfieldHeight)
                : SampleBackgroundGradient(dir));
        result += transmittance * bg;
    }
    if (volumeOpticalDepth > 0.01f) {
        result = result * exp(-volumeOpticalDepth) + volumeEmission;
    }

    // Beam expansion: the transverse deviation magnitude relative to its seed. The
    // alpha channel carries it when ray bundles are on (a parity handle against
    // the CPU bundle); otherwise alpha stays 1, so the default render is unchanged.
    float alpha = 1.0f;
    if (beamsEnabled && !captured) {
        alpha = beamExpansion;
    }

    int pixel = ty * tileWidth + tx;
    if (sampleIndex > 0.5f) {
        float previousWeight = sampleIndex;
        float inverseCount = 1.0f / (previousWeight + 1.0f);
        result = (float3(radiance[pixel * 4 + 0], radiance[pixel * 4 + 1],
                         radiance[pixel * 4 + 2])
                      * previousWeight
                  + result)
                 * inverseCount;
        alpha = (radiance[pixel * 4 + 3] * previousWeight + alpha) * inverseCount;
    }
    radiance[pixel * 4 + 0] = result.x;
    radiance[pixel * 4 + 1] = result.y;
    radiance[pixel * 4 + 2] = result.z;
    radiance[pixel * 4 + 3] = alpha;
}
