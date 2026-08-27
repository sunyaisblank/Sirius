// Per-pixel geodesic tracer, the Vulkan render path. Each thread renders one
// pixel of a TILE: it builds the camera null ray, integrates it through the
// dispatched metric family in the Cartesian chart, and terminates on horizon
// capture or escape. Captured rays shade black (the shadow); escaping rays shade
// by their asymptotic direction against the starfield; an optional thin disk
// adds Novikov-Thorne blackbody emission with g^4 relativistic beaming where the
// geodesic crosses the equatorial plane. Volumetric mode instead integrates
// emission/absorption through every geodesic segment, with the same deterministic
// turbulence and extended inverse-Compton corona represented by the CPU tracer.
// The kernel returns LINEAR radiance per pixel; tonemapping and grading are the
// host display pipeline's job.
//
// This is render orchestration, not a physics authority: the metric,
// Christoffel, integrator, blackbody, and disk-redshift routines it calls are
// the ported authorities in the gr_* modules. Termination radii and the horizon
// buffer are carried from RDOP002A.cu / core/constants.h.
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

// Same local-frame light-cone boost as core/camera.h::AberrateRay. beta and n
// use the camera-local (radial, theta, phi) component order.
float3 AberrateLocal(float3 n, float3 beta) {
    float b2 = dot(beta, beta);
    if (b2 <= 0.0f) return n;
    float b = sqrt(b2);
    float3 axis = beta / b;
    float nPar = dot(n, axis);
    float3 perpendicular = n - nPar * axis;
    float gamma = 1.0f / sqrt(1.0f - b2);
    float denom = 1.0f - b * nPar;
    float nParPrime = (nPar - b) / denom;
    float perpendicularScale = 1.0f / (gamma * denom);
    return normalize(nParPrime * axis + perpendicular * perpendicularScale);
}

// Solve the null condition g_uv u^u u^v = 0 for u^t given the spatial
// direction, choosing the future-directed root (RDOP002A.cu
// cartesianDirToFourVelocity, Cartesian form). Real: on the fp64 rung the
// launch state is as precise as the trajectory it seeds.
Real SolveNullTimeComponent(Real g[4][4], float3 dir) {
    Real A = g[0][0];
    Real B = 2.0f * (g[0][1] * dir.x + g[0][2] * dir.y + g[0][3] * dir.z);
    Real C = g[1][1] * dir.x * dir.x + g[2][2] * dir.y * dir.y + g[3][3] * dir.z * dir.z
            + 2.0f * (g[1][2] * dir.x * dir.y + g[1][3] * dir.x * dir.z + g[2][3] * dir.y * dir.z);
    Real discriminant = B * B - 4.0f * A * C;
    if (discriminant >= 0.0f) {
        Real sqrtD = sqrt(discriminant);
        Real root1 = (-B + sqrtD) / (2.0f * A);
        Real root2 = (-B - sqrtD) / (2.0f * A);
        return (root1 > 0.0f) ? root1 : ((root2 > 0.0f) ? root2 : Real(1.0f));
    }
    return sqrt(max(-C / A, Real(0.0001f)));
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
    if (temperature <= 0.0f) temperature = 5778.0f;
    float normalized = temperature / 6500.0f;
    float3 color;
    if (temperature <= 6500.0f) {
        color = float3(1.0f, pow(normalized, 0.4f), pow(normalized, 0.8f));
    } else {
        color = float3(pow(1.0f / normalized, 0.3f),
                       pow(1.0f / normalized, 0.15f), 1.0f);
    }
    return color / max(color.x, max(color.y, color.z));
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

// Deterministic turbulence authority mirrored from core/disk/turbulence.h. The
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
        ? innerTemperature * pow(innerRadius / cylindricalR, 0.75f)
        : FullPageThorneTemperature(innerTemperature, cylindricalR, rIsco, M, aStar);
    if (midplane <= 0.0f) return 0.0f;
    float H = VolumeScaleHeight(cylindricalR, innerRadius, hOverR, hPower);
    float zOverHSquared = min((z / H) * (z / H), 9.0f);
    static const float kAtmosphereFourthPower = 0.4096f;  // 0.8^4.
    float fourthPower =
        max(1.0f - zOverHSquared * (1.0f - kAtmosphereFourthPower) / 9.0f,
            kAtmosphereFourthPower);
    return midplane * pow(fourthPower, 0.25f);
}

// Integrate one geodesic segment using the same midpoint transfer recurrence as
// GeodesicTracer::AccumulateVolumetricEmission. The corona config surface is a
// boolean today, so the Extended geometry and physical constants are the exact
// CPU defaults.
void AccumulateVolumeSegment(float3 start, float3 finish, int sampleCount,
                             bool volumeEnabled, bool turbulenceEnabled,
                             bool coronaEnabled, bool dopplerBeaming, float innerRadius,
                             float outerRadius, float innerTemperature,
                             int temperatureModel, float rIsco, float M,
                             float a, float aStar, float hOverR, float hPower,
                             float tauMidplane, inout float opticalDepth,
                             inout float3 emission) {
    float3 delta = finish - start;
    float pathLength = length(delta);
    if (pathLength < 1.0e-10f || opticalDepth >= 10.0f) return;
    float ds = pathLength / float(sampleCount);
    float3 direction = delta / pathLength;
    float referenceH = VolumeScaleHeight(innerRadius, innerRadius, hOverR, hPower);

    for (int sample = 0; sample < sampleCount; ++sample) {
        float3 point = start + (float(sample) + 0.5f) * ds * direction;
        float sphericalR = length(point);
        float diskR =
            float(ComputeKerrRadiusCart(point.x, point.y, point.z, Real(a)));
        float diskOpacity = 0.0f;
        float diskSource = 0.0f;

        if (volumeEnabled &&
            IsInsideVolumeDisk(diskR, point.z, innerRadius, outerRadius,
                               hOverR, hPower)) {
            float H = VolumeScaleHeight(diskR, innerRadius, hOverR, hPower);
            float zOverH = point.z / H;
            diskOpacity = tauMidplane / (sqrt(TWO_PI) * H) *
                          pow(diskR / innerRadius, -1.5f) *
                          (referenceH / H) * exp(-0.5f * zOverH * zOverH);
            if (turbulenceEnabled) {
                diskOpacity *= TurbulenceDensityFactor(point);
            }
            float temperature =
                VolumeDiskTemperature(diskR, point.z, innerRadius,
                                      innerTemperature, temperatureModel, rIsco,
                                      M, aStar, hOverR, hPower);
            if (diskOpacity > 0.0f && temperature > 0.0f) {
                float emittedSource = pow(temperature / innerTemperature, 4.0f);
                float phi = atan2(point.y, point.x);
                float gFactor = dopplerBeaming
                    ? ComputeDiskGFactorCart(diskR, phi, direction, M, a)
                    : ComputeDiskGravitationalFactorCart(diskR, M, a);
                diskSource = ObservedDiskBolometricSource(emittedSource, gFactor);
            }
        }

        float coronaOpacity = 0.0f;
        float coronaSource = 0.0f;
        if (coronaEnabled && sphericalR >= innerRadius && sphericalR <= outerRadius) {
            float coronaH = 5.0f * sqrt(sphericalR / innerRadius);
            if (abs(point.z) < coronaH) {
                coronaOpacity = 0.5f / 5.0f;
                float seed = pow(innerRadius / sphericalR, 3.0f) *
                             exp(-0.5f * point.z * point.z / (coronaH * coronaH));
                // Default y=4*(100/511)*max(0.5,0.5^2), and y<1 => boost=1+y.
                coronaSource = seed * (1.0f + 4.0f * (100.0f / 511.0f) * 0.5f);
            }
        }

        float totalOpacity = diskOpacity + coronaOpacity;
        if (totalOpacity <= 0.0f) continue;
        float3 source = float3(
            diskOpacity * diskSource + coronaOpacity * coronaSource * 0.65f,
            diskOpacity * diskSource + coronaOpacity * coronaSource * 0.80f,
            diskOpacity * diskSource + coronaOpacity * coronaSource) / totalOpacity;
        float deltaTau = min(totalOpacity * ds, 10.0f - opticalDepth);
        float transmission = exp(-deltaTau);
        emission = emission * transmission + source * (1.0f - transmission);
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

    int tileOriginX = int(params[32]);
    int tileOriginY = int(params[33]);
    int tileWidth = int(params[34]);
    int tileHeight = int(params[35]);

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
    float emissionScale = params[31];

    bool starfieldEnabled = params[36] > 0.5f;
    int starfieldWidth = int(params[37]);
    int starfieldHeight = int(params[38]);
    // params[39] reserved.

    float warpVs = params[40];
    float warpSigma = params[41];
    float warpR = params[42];

    // Morris-Thorne (dispatch id 3): throat radius and shape selector
    // (0 Ellis, 1 zero-tidal, 2 absurdly benign); Phi0 = 0 is the zero-tidal
    // subclass the CPU authority carries.
    float throatB0 = params[44];
    int wormShape = int(params[45]);
    float sampleU = params[46];
    float sampleV = params[47];
    float sampleIndex = params[48];
    float3 cameraBeta = P3(49);
    int temperatureModel = int(params[52]);  // 0 Page/Novikov-Thorne, 1 Shakura-Sunyaev.
    bool dopplerBeaming = params[53] > 0.5f;
    int lensModel = int(params[54]);  // 0 pinhole, 1 deterministic thin lens.
    float focalLength = params[55];
    float aperture = params[56];
    float focusDistance = params[57];
    bool pointStarfieldEnabled = params[58] > 0.5f;
    float pixelAngularSize = params[59];
    float pointBrightnessScale = params[60];
    bool volumeEnabled = params[61] > 0.5f;
    float volumeHOverR = params[62];
    float volumeHPower = params[63];
    float volumeTauMidplane = params[64];
    int volumeSamples = int(params[65]);
    bool turbulenceEnabled = params[66] > 0.5f;
    bool coronaEnabled = params[67] > 0.5f;

    // Ray bundles (P2): propagate two geodesic-deviation vectors through the render
    // loop and report the beam's transverse ellipse. Default off (params[43]==0)
    // leaves the RGBA output untouched, so the parity and pinned renders do not
    // move. Uses gr_deviation's full-Riemann acceleration (DNGR eq A.19) in the
    // same Cartesian chart the integrator runs in.
    bool beamsEnabled = params[43] > 0.5f;

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
    localDir = AberrateLocal(localDir, cameraBeta);
    float3 rayDir =
        normalize(camDir * (-localDir.x) + camUp * localDir.y + camRight * localDir.z);

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
    Real ut = SolveNullTimeComponent(g, rayDir);

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
    state.u = MakeVec4Cart(ut, rayDir.x, rayDir.y, rayDir.z);

    float3 color = float3(0.0f, 0.0f, 0.0f);
    float transmittance = 1.0f;
    float volumeOpticalDepth = 0.0f;
    float3 volumeEmission = float3(0.0f);
    bool escaped = false;

    // Two deviation vectors span the ray bundle (P2). A point-source/pupil bundle
    // is used when the point-star field consumes the ellipse; the diagnostic beam
    // path otherwise retains the parallel-bundle initial condition.
    float3 seed = (abs(rayDir.z) < 0.9f) ? float3(0.0f, 0.0f, 1.0f) : float3(1.0f, 0.0f, 0.0f);
    float3 xiPerpA = normalize(seed - dot(seed, rayDir) * rayDir);
    float3 xiPerpB = cross(rayDir, xiPerpA);
    const float kBeamSeed = 1.0e-3f;
    bool pupilBundle = pointStarfieldEnabled;
    Real4 xiA = pupilBundle
        ? Real4(0.0f)
        : Real4(0.0f, xiPerpA.x, xiPerpA.y, xiPerpA.z) * Real(kBeamSeed);
    Real4 xiB = pupilBundle
        ? Real4(0.0f)
        : Real4(0.0f, xiPerpB.x, xiPerpB.y, xiPerpB.z) * Real(kBeamSeed);
    Real4 xiDotA = pupilBundle
        ? Real4(0.0f, xiPerpA.x, xiPerpA.y, xiPerpA.z) * Real(kBeamSeed)
        : Real4(0.0f);
    Real4 xiDotB = pupilBundle
        ? Real4(0.0f, xiPerpB.x, xiPerpB.y, xiPerpB.z) * Real(kBeamSeed)
        : Real4(0.0f);
    Real beamAffineLength = 0.0f;

#ifdef SIRIUS_FP32_COMP
    // Compensated rung: running Kahan compensation for the state accumulation
    // across steps, recovering the low-order bits plain fp32 addition drops.
    // The RK4 increment itself is unchanged; only the commit differs.
    Vec4Cart compX = MakeVec4Cart(0.0f, 0.0f, 0.0f, 0.0f);
    Vec4Cart compU = MakeVec4Cart(0.0f, 0.0f, 0.0f, 0.0f);
#endif

    for (int step = 0; step < maxSteps; step++) {
        Real px = state.x.x;
        Real py = state.x.y;
        Real pz = state.x.z;
        Real rBig = sqrt(px * px + py * py + pz * pz);

        if (isWormhole) {
            if (rBig <= horizonRadius) { break; }  // reached the throat: dark
        } else if (!isWarp) {
            Real rKerr = ComputeKerrRadiusCart(px, py, pz, Real(a));
            if (rKerr <= horizonRadius) { break; }  // captured: shadow
        }
        if (rBig >= escapeRadius) { escaped = true; break; }

        Real rStep = (isWarp || isWormhole)
            ? max(rBig, Real(1.0f))
            : ComputeKerrRadiusCart(px, py, pz, Real(a));
        Real dt = clamp(stepScale * rStep, Real(minStep), Real(maxStep));

        Vec4Cart dX, dU;
        if (isWarp) {
            IntegrateGeodesicRK4WarpCartDelta(state, dt, warpVs, warpSigma, warpR,
                                              0.0f, 0.0f, 0.0f, dX, dU);
        } else if (isWormhole) {
            IntegrateGeodesicRK4MorrisThorneCartDelta(state, dt, throatB0, 0.0f, wormShape,
                                                      dX, dU);
        } else {
            IntegrateGeodesicRK4CartDelta(state, dt, M, a, Q, Lambda, dX, dU);
        }
        GeodesicStateCart next;
#ifdef SIRIUS_FP32_COMP
        next.x = KahanAddVec4Cart(state.x, compX, dX);
        next.u = KahanAddVec4Cart(state.u, compU, dU);
#else
        // Plain accumulation: identical arithmetic to the pre-refactor
        // integrator, which performed these adds internally.
        next.x = AddVec4Cart(state.x, dX);
        next.u = AddVec4Cart(state.u, dU);
#endif

        // Optional volumetric disk/corona (Kerr-Schild black holes only).
        if ((volumeEnabled || coronaEnabled) && !isWarp && !isWormhole) {
            float aStar = (M > 1e-6f) ? (a / M) : a;
            float rIsco = ComputeKerrISCO(M, aStar);
            AccumulateVolumeSegment(
                float3(float(state.x.x), float(state.x.y), float(state.x.z)),
                float3(float(next.x.x), float(next.x.y), float(next.x.z)),
                volumeSamples, volumeEnabled, turbulenceEnabled, coronaEnabled,
                dopplerBeaming,
                innerRadius, outerRadius, innerTemperature, temperatureModel,
                rIsco, M, a, aStar, volumeHOverR, volumeHPower,
                volumeTauMidplane, volumeOpticalDepth, volumeEmission);
        }

        // Optional thin disk (Kerr-Schild black holes only): emission when the
        // segment crosses the equatorial plane z=0 within [innerRadius,
        // outerRadius], with Novikov-Thorne temperature, blackbody colour, and
        // g^4 relativistic beaming (gr_disk authorities).
        if (diskEnabled && !volumeEnabled && !isWarp && !isWormhole &&
            transmittance > 0.01f) {
            if (state.x.z * next.x.z < 0.0f) {
                Real frac = state.x.z / (state.x.z - next.x.z);
                Real cxr = state.x.x + frac * (next.x.x - state.x.x);
                Real cyr = state.x.y + frac * (next.x.y - state.x.y);
                // Shading is float (gr_disk carries trig/pow, 32-bit-only in
                // SPIR-V); the crossing point itself was found at Real.
                float cx = float(cxr);
                float cy = float(cyr);
                float cr =
                    float(ComputeKerrRadiusCart(cxr, cyr, Real(0.0f), Real(a)));
                if (cr >= innerRadius && cr <= outerRadius) {
                    float aStar = (M > 1e-6f) ? (a / M) : a;
                    float rIsco = ComputeKerrISCO(M, aStar);
                    float T;
                    if (temperatureModel == 1) {
                        T = innerTemperature * pow(innerRadius / cr, 0.75f);
                    } else {
                        T = FullPageThorneTemperature(innerTemperature, cr, rIsco, M, aStar);
                    }
                    float cphi = atan2(cy, cx);
                    float3 photonVel = float3(float(next.u.x), float(next.u.y), float(next.u.z));
                    float gFactor = dopplerBeaming
                        ? ComputeDiskGFactorCart(cr, cphi, photonVel, M, a)
                        : ComputeDiskGravitationalFactorCart(cr, M, a);
                    float observedScale =
                        ObservedDiskBolometricSource(emissionScale, gFactor);
                    float3 emit = BlackbodyColor(T) * observedScale;
                    color += transmittance * emit;
                    transmittance *= 0.2f;  // partial occlusion per crossing
                }
            }
        }
        // Advance the deviation vector over this step (velocity Verlet on the
        // full-Riemann Jacobi acceleration), coefficients at the step start.
        if (beamsEnabled && !isWarp && !isWormhole) {
            Real beamRiemann[4][4][4][4];
            GetRiemannTensorCart(state.x, M, a, Q, Lambda, beamRiemann);
            Real4 accelA =
                GeodesicDeviationAccelFromRiemann(beamRiemann, state.u, xiA);
            Real4 accelB =
                GeodesicDeviationAccelFromRiemann(beamRiemann, state.u, xiB);
            xiA += xiDotA * dt + accelA * (Real(0.5f) * dt * dt);
            xiB += xiDotB * dt + accelB * (Real(0.5f) * dt * dt);
            xiDotA += accelA * dt;
            xiDotB += accelB * dt;
            beamAffineLength += dt;
        }

        state = next;
    }

    float3 terminalDirection =
        normalize(float3(float(state.u.x), float(state.u.y), float(state.u.z)));
    float3 beamReference = abs(terminalDirection.z) > 0.9f
        ? float3(1.0f, 0.0f, 0.0f)
        : float3(0.0f, 0.0f, 1.0f);
    float3 beamTangentA =
        normalize(beamReference - dot(beamReference, terminalDirection) * terminalDirection);
    float3 beamTangentB = cross(terminalDirection, beamTangentA);
    float3 xiASpatial = float3(float(xiA.y), float(xiA.z), float(xiA.w));
    float3 xiBSpatial = float3(float(xiB.y), float(xiB.z), float(xiB.w));
    float beamA = dot(xiASpatial, beamTangentA);
    float beamB = dot(xiBSpatial, beamTangentA);
    float beamC = dot(xiASpatial, beamTangentB);
    float beamD = dot(xiBSpatial, beamTangentB);
    float beamSemiMajor, beamSemiMinor, beamOrientation, beamDeterminant;
    BeamEllipseFromTransverseMap(beamA, beamB, beamC, beamD,
                                 beamSemiMajor, beamSemiMinor, beamOrientation,
                                 beamDeterminant);
    float beamNormalization = pupilBundle
        ? max(float(beamAffineLength) * kBeamSeed, 1.0e-12f)
        : kBeamSeed;
    float beamMajorExpansion = beamSemiMajor / beamNormalization;
    float beamMinorExpansion = beamSemiMinor / beamNormalization;
    float beamExpansion = sqrt(max(beamMajorExpansion * beamMinorExpansion, 0.0f));

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
        float effectiveTemperature =
            clamp(pow(max(volumeEmission.x, 0.0f), 0.25f) * innerTemperature,
                  1000.0f, 100000.0f);
        float3 volumeColor =
            BlackbodyColor(effectiveTemperature) * volumeEmission;
        result = result * exp(-volumeOpticalDepth) + volumeColor;
    }

    // Beam expansion: the transverse deviation magnitude relative to its seed. The
    // alpha channel carries it when ray bundles are on (a parity handle against
    // the CPU bundle); otherwise alpha stays 1, so the default render is unchanged.
    float alpha = 1.0f;
    if (beamsEnabled) {
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
