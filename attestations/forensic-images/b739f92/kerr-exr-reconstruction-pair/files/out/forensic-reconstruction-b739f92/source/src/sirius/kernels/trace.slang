// Per-pixel geodesic tracer. Each thread renders one pixel of a tile: it builds
// the camera null ray, integrates it through the Kerr-Schild family in the
// Cartesian chart (gr_integrator IntegrateGeodesicRK4Cart), and terminates on
// horizon capture or escape. Captured rays shade black (the shadow); escaping
// rays shade by their asymptotic direction. An optional thin disk adds
// Novikov-Thorne blackbody emission where the geodesic crosses the equatorial
// plane.
//
// This is render orchestration, not a physics authority: the metric,
// Christoffel, integrator, and blackbody routines it calls are the ported
// authorities in the gr_* modules. Termination radii and the horizon buffer are
// carried from RDOP002A.cu / core/constants.h.
//
// Buffer ABI (bound 0..1 in declaration order):
//   0 radiance : RWStructuredBuffer<float>, width*height*4 (RGBA per pixel).
//   1 params   : StructuredBuffer<float>, layout documented below.
import gr_types;
import gr_metrics;
import gr_integrator;
import gr_disk;

RWStructuredBuffer<float> radiance;  // binding 0
StructuredBuffer<float> params;      // binding 1

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

// Smooth escape-direction background: a bounded gradient so the escaped field
// is finite and non-constant across directions (the smoke gate). A production
// starfield replaces this behind the same seam.
float3 SampleBackground(float3 dir) {
    return float3(0.5f, 0.5f, 0.5f) + 0.5f * dir;
}

[shader("compute")]
[numthreads(8, 8, 1)]
void ComputeMain(uint3 tid: SV_DispatchThreadID) {
    int width = int(params[0]);
    int height = int(params[1]);
    int ix = int(tid.x);
    int iy = int(tid.y);
    if (ix >= width || iy >= height) return;

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

    // Camera ray in world Cartesian (RDOP002A.cu raygen NDC construction).
    float u = (2.0f * (float(ix) + 0.5f) / float(width) - 1.0f) * aspect;
    float v = 2.0f * (float(iy) + 0.5f) / float(height) - 1.0f;
    float tanHalfFov = tan(fov * 0.5f);
    u *= tanHalfFov;
    v *= tanHalfFov;
    float3 rayDir = normalize(camDir + camRight * u + camUp * v);

    // Horizon capture radius: r_+ = M + sqrt(M^2 - a^2 - Q^2), buffered.
    float disc = M * M - a * a - Q * Q;
    float horizonRadius = (disc > 0.0f) ? (M + sqrt(disc)) : M;
    horizonRadius *= captureFactor;

    // Initial 4-velocity from the null condition at the camera position.
    float g[4][4];
    float g_inv[4][4];
    GetKerrSchildFamilyMetric(MakeVec4Cart(0.0f, camPos.x, camPos.y, camPos.z), M, a, Q, Lambda, g, g_inv);
    float ut = SolveNullTimeComponent(g, rayDir);

    GeodesicStateCart state;
    state.x = MakeVec4Cart(0.0f, camPos.x, camPos.y, camPos.z);
    state.u = MakeVec4Cart(ut, rayDir.x, rayDir.y, rayDir.z);

    float3 color = float3(0.0f, 0.0f, 0.0f);
    float transmittance = 1.0f;
    bool captured = false;
    bool escaped = false;

    for (int step = 0; step < maxSteps; step++) {
        float px = state.x.x;
        float py = state.x.y;
        float pz = state.x.z;
        float rBig = sqrt(px * px + py * py + pz * pz);
        float rKerr = ComputeKerrRadiusCart(px, py, pz, a);

        if (rKerr <= horizonRadius) { captured = true; break; }
        if (rBig >= escapeRadius) { escaped = true; break; }

        float dt = clamp(stepScale * rKerr, minStep, maxStep);
        GeodesicStateCart next = IntegrateGeodesicRK4Cart(state, dt, M, a, Q, Lambda);

        // Optional thin disk: emission when the segment crosses the equatorial
        // plane z=0 within [innerRadius, outerRadius].
        if (diskEnabled && transmittance > 0.01f) {
            if (state.x.z * next.x.z < 0.0f) {
                float frac = state.x.z / (state.x.z - next.x.z);
                float cx = state.x.x + frac * (next.x.x - state.x.x);
                float cy = state.x.y + frac * (next.x.y - state.x.y);
                float cr = sqrt(cx * cx + cy * cy);
                if (cr >= innerRadius && cr <= outerRadius) {
                    float rIsco = ComputeKerrISCO(M, a);
                    float Qflux = ComputeNovikovThorneQ(cr, rIsco, M, a);
                    float T = NovikovThorneTemperature(innerTemperature, Qflux);
                    float3 emit = BlackbodyColor(T) * emissionScale;
                    color += transmittance * emit;
                    transmittance *= 0.2f;  // partial occlusion per crossing
                }
            }
        }
        state = next;
    }

    float3 result = color;
    if (escaped) {
        float3 dir = normalize(float3(state.u.x, state.u.y, state.u.z));
        result += transmittance * SampleBackground(dir);
    }
    // Captured (shadow) pixels keep only whatever disk emission accumulated.

    int pixel = iy * width + ix;
    radiance[pixel * 4 + 0] = result.x;
    radiance[pixel * 4 + 1] = result.y;
    radiance[pixel * 4 + 2] = result.z;
    radiance[pixel * 4 + 3] = 1.0f;
}
