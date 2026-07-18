// Per-pixel geodesic tracer, the Vulkan render path. Each thread renders one
// pixel of a TILE: it builds the camera null ray, integrates it through the
// dispatched metric family in the Cartesian chart, and terminates on horizon
// capture or escape. Captured rays shade black (the shadow); escaping rays shade
// by their asymptotic direction against the starfield; an optional thin disk
// adds Novikov-Thorne blackbody emission with g^4 relativistic beaming where the
// geodesic crosses the equatorial plane. The kernel returns LINEAR radiance per
// pixel; tonemapping and grading are the host display pipeline's job.
//
// This is render orchestration, not a physics authority: the metric,
// Christoffel, integrator, blackbody, and disk-redshift routines it calls are
// the ported authorities in the gr_* modules. Termination radii and the horizon
// buffer are carried from RDOP002A.cu / core/constants.h.
//
// Metric dispatch (params[2]) matches the registry gpu_supported set the render
// session wires: Kerr-Schild family (Minkowski/Schwarzschild/Kerr) and the
// Alcubierre warp family. The gr_* modules also carry Ellis/Morris-Thorne
// physics (parity_probe covers them); their spherical-chart render integration
// is a recorded follow-up and the session declines them, so no dispatch id here.
//
// Buffer ABI (bound 0..2 in declaration order):
//   0 radiance   : RWStructuredBuffer<float>, tileWidth*tileHeight*4 (RGBA/px).
//   1 params     : StructuredBuffer<float>, layout documented below.
//   2 starfield  : StructuredBuffer<uint>, packed RGBA8 equirectangular texels
//                  (R | G<<8 | B<<16 | A<<24), size starfieldWidth*starfieldHeight;
//                  bound with a 1-texel dummy when starfieldEnabled is 0.
import gr_types;
import gr_metrics;
import gr_integrator;
import gr_disk;
import gr_deviation;

RWStructuredBuffer<float> radiance;  // binding 0
StructuredBuffer<float> params;      // binding 1
StructuredBuffer<uint> starfield;    // binding 2

// Kernel-side metric dispatch ids (subset of gr_types kMetric*; the ids the
// render session emits). Kept local so this file reads without cross-reference.
static const uint kDispatchKerrSchild = 0;
static const uint kDispatchWarpDrive = 2;

float3 P3(uint base) { return float3(params[base], params[base + 1], params[base + 2]); }

// Solve the null condition g_uv u^u u^v = 0 for u^t given the spatial
// direction, choosing the future-directed root (RDOP002A.cu
// cartesianDirToFourVelocity, Cartesian form).
float SolveNullTimeComponent(float g[4][4], float3 dir) {
    float A = g[0][0];
    float B = 2.0f * (g[0][1] * dir.x + g[0][2] * dir.y + g[0][3] * dir.z);
    float C = g[1][1] * dir.x * dir.x + g[2][2] * dir.y * dir.y + g[3][3] * dir.z * dir.z
            + 2.0f * (g[1][2] * dir.x * dir.y + g[1][3] * dir.x * dir.z + g[2][3] * dir.y * dir.z);
    float discriminant = B * B - 4.0f * A * C;
    if (discriminant >= 0.0f) {
        float sqrtD = sqrt(discriminant);
        float root1 = (-B + sqrtD) / (2.0f * A);
        float root2 = (-B - sqrtD) / (2.0f * A);
        return (root1 > 0.0f) ? root1 : ((root2 > 0.0f) ? root2 : 1.0f);
    }
    return sqrt(max(-C / A, 0.0001f));
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
    // params[39] reserved (background fallback floor); the render path uses the
    // analytic gradient below when no texture is uploaded.

    float warpVs = params[40];
    float warpSigma = params[41];
    float warpR = params[42];

    // Ray bundles (P2): propagate one geodesic-deviation vector through the render
    // loop and report the beam's transverse expansion. Default off (params[43]==0)
    // leaves the RGBA output untouched, so the parity and pinned renders do not
    // move. Uses gr_deviation's full-Riemann acceleration (DNGR eq A.19) in the
    // same Cartesian chart the integrator runs in.
    bool beamsEnabled = params[43] > 0.5f;

    // Camera ray in world Cartesian. NDC follows the CPU pinhole camera geometry
    // (the renderer builds camDir/camRight/camUp so the two paths trace the same
    // scene): u spans [-aspect, aspect], v spans [-1, 1] top-to-bottom.
    float u = (2.0f * (float(ix) + 0.5f) / float(imageWidth) - 1.0f) * aspect;
    float v = 2.0f * (float(iy) + 0.5f) / float(imageHeight) - 1.0f;
    float tanHalfFov = tan(fov * 0.5f);
    u *= tanHalfFov;
    v *= tanHalfFov;
    float3 rayDir = normalize(camDir + camRight * u + camUp * v);

    bool isWarp = (metricId == kDispatchWarpDrive);

    // Initial 4-velocity from the null condition at the camera position, using
    // the dispatched metric.
    float g[4][4];
    float g_inv[4][4];
    Vec4Cart pos0 = MakeVec4Cart(0.0f, camPos.x, camPos.y, camPos.z);
    if (isWarp) {
        GetWarpDriveMetric(pos0, warpVs, warpSigma, warpR, 0.0f, 0.0f, 0.0f, g, g_inv);
    } else {
        GetKerrSchildFamilyMetric(pos0, M, a, Q, Lambda, g, g_inv);
    }
    float ut = SolveNullTimeComponent(g, rayDir);

    // Horizon capture radius (Kerr-Schild family only): r_+ = M + sqrt(M^2 - a^2
    // - Q^2), buffered. The warp bubble is horizonless: capture is disabled.
    float disc = M * M - a * a - Q * Q;
    float horizonRadius = (disc > 0.0f) ? (M + sqrt(disc)) : M;
    horizonRadius *= captureFactor;

    GeodesicStateCart state;
    state.x = pos0;
    state.u = MakeVec4Cart(ut, rayDir.x, rayDir.y, rayDir.z);

    float3 color = float3(0.0f, 0.0f, 0.0f);
    float transmittance = 1.0f;
    bool escaped = false;

    // Deviation vector for the ray bundle (P2), transverse to the initial ray, and
    // its affine derivative; a parallel bundle (xiDot = 0) so the expansion below
    // measures the lensing distortion of the pixel footprint.
    float3 seed = (abs(rayDir.z) < 0.9f) ? float3(0.0f, 0.0f, 1.0f) : float3(1.0f, 0.0f, 0.0f);
    float3 xiPerp = normalize(cross(rayDir, seed));
    const float kBeamSeed = 1.0e-3f;
    float4 xi = float4(0.0f, xiPerp.x, xiPerp.y, xiPerp.z) * kBeamSeed;
    float4 xiDot = float4(0.0f, 0.0f, 0.0f, 0.0f);

    for (int step = 0; step < maxSteps; step++) {
        float px = state.x.x;
        float py = state.x.y;
        float pz = state.x.z;
        float rBig = sqrt(px * px + py * py + pz * pz);

        if (!isWarp) {
            float rKerr = ComputeKerrRadiusCart(px, py, pz, a);
            if (rKerr <= horizonRadius) { break; }  // captured: shadow
        }
        if (rBig >= escapeRadius) { escaped = true; break; }

        float rStep = isWarp ? max(rBig, 1.0f) : ComputeKerrRadiusCart(px, py, pz, a);
        float dt = clamp(stepScale * rStep, minStep, maxStep);

        GeodesicStateCart next = isWarp
            ? IntegrateGeodesicRK4WarpCart(state, dt, warpVs, warpSigma, warpR, 0.0f, 0.0f, 0.0f)
            : IntegrateGeodesicRK4Cart(state, dt, M, a, Q, Lambda);

        // Optional thin disk (Kerr-Schild black holes only): emission when the
        // segment crosses the equatorial plane z=0 within [innerRadius,
        // outerRadius], with Novikov-Thorne temperature, blackbody colour, and
        // g^4 relativistic beaming (gr_disk authorities).
        if (diskEnabled && !isWarp && transmittance > 0.01f) {
            if (state.x.z * next.x.z < 0.0f) {
                float frac = state.x.z / (state.x.z - next.x.z);
                float cx = state.x.x + frac * (next.x.x - state.x.x);
                float cy = state.x.y + frac * (next.x.y - state.x.y);
                float cr = sqrt(cx * cx + cy * cy);
                if (cr >= innerRadius && cr <= outerRadius) {
                    float aStar = (M > 1e-6f) ? (a / M) : a;
                    float rIsco = ComputeKerrISCO(M, aStar);
                    float Qflux = ComputeNovikovThorneQ(cr, rIsco, M, aStar);
                    float T = NovikovThorneTemperature(innerTemperature, Qflux);
                    float cphi = atan2(cy, cx);
                    float3 photonVel = float3(next.u.x, next.u.y, next.u.z);
                    float gFactor = ComputeDiskGFactorCart(cr, cphi, photonVel, M, a);
                    float g4 = gFactor * gFactor * gFactor * gFactor;
                    float3 emit = BlackbodyColor(T) * emissionScale * g4;
                    color += transmittance * emit;
                    transmittance *= 0.2f;  // partial occlusion per crossing
                }
            }
        }
        // Advance the deviation vector over this step (velocity Verlet on the
        // full-Riemann Jacobi acceleration), coefficients at the step start.
        if (beamsEnabled && !isWarp) {
            float4 accel = GeodesicDeviationAccelFullRiemann(state.x, state.u, xi, M, a, Q, Lambda);
            xi += xiDot * dt + accel * (0.5f * dt * dt);
            xiDot += accel * dt;
        }

        state = next;
    }

    // Captured (shadow) pixels keep only whatever disk emission accumulated
    // (black when none), exactly as the CPU horizon outcome does.
    float3 result = color;
    if (escaped) {
        float3 dir = normalize(float3(state.u.x, state.u.y, state.u.z));
        float3 bg = starfieldEnabled
            ? SampleStarfield(dir, starfieldWidth, starfieldHeight)
            : SampleBackgroundGradient(dir);
        result += transmittance * bg;
    }

    // Beam expansion: the transverse deviation magnitude relative to its seed. The
    // alpha channel carries it when ray bundles are on (a parity handle against
    // the CPU bundle); otherwise alpha stays 1, so the default render is unchanged.
    float alpha = 1.0f;
    if (beamsEnabled) {
        float xiMag = length(float3(xi.y, xi.z, xi.w));
        alpha = xiMag / kBeamSeed;
    }

    int pixel = ty * tileWidth + tx;
    radiance[pixel * 4 + 0] = result.x;
    radiance[pixel * 4 + 1] = result.y;
    radiance[pixel * 4 + 2] = result.z;
    radiance[pixel * 4 + 3] = alpha;
}
