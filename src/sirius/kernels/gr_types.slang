// State structs and launch parameters shared by every device-physics kernel.
//
// Re-authored from the retired CUDA monolith
// (src/Sirius.Render/Acceleration/OptiX/RDOP002A.cu and RDOP003A.h). These are
// value types only: no global resources are declared here, so descriptor
// bindings stay owned by the entry-point kernels (trace, parity_probe) and
// remain deterministic in declaration order 0..N-1.
//
// Coordinate charts follow the .cu exactly: Vec4 carries Boyer-Lindquist
// (t, r, theta, phi); Vec4Cart carries Kerr-Schild Cartesian (t, x, y, z).
// All arithmetic is single precision, matching the .cu's float physics
// (docs/ARCHITECTURE.md section 6 precision ladder, plain-fp32 rung).
module gr_types;

// pi to float precision, as the .cu spelled it (constexpr float PI).
public static const float PI = 3.14159265358979323846f;
public static const float TWO_PI = 6.28318530717958647692f;

// --- 4-vectors ---------------------------------------------------------------

// Boyer-Lindquist / spherical 4-vector (t, r, theta, phi).
public struct Vec4 {
    public float t;
    public float r;
    public float theta;
    public float phi;
}

// Kerr-Schild Cartesian 4-vector (t, x, y, z). The Kerr-Schild family is
// authored in this chart because its Christoffels are polynomial with no pole
// singularity.
public struct Vec4Cart {
    public float t;
    public float x;
    public float y;
    public float z;
}

public Vec4 MakeVec4(float t, float r, float theta, float phi) {
    Vec4 v;
    v.t = t; v.r = r; v.theta = theta; v.phi = phi;
    return v;
}

public Vec4Cart MakeVec4Cart(float t, float x, float y, float z) {
    Vec4Cart v;
    v.t = t; v.x = x; v.y = y; v.z = z;
    return v;
}

public Vec4 AddVec4(Vec4 a, Vec4 b) {
    return MakeVec4(a.t + b.t, a.r + b.r, a.theta + b.theta, a.phi + b.phi);
}

public Vec4 ScaleVec4(Vec4 a, float s) {
    return MakeVec4(a.t * s, a.r * s, a.theta * s, a.phi * s);
}

public Vec4Cart AddVec4Cart(Vec4Cart a, Vec4Cart b) {
    return MakeVec4Cart(a.t + b.t, a.x + b.x, a.y + b.y, a.z + b.z);
}

public Vec4Cart ScaleVec4Cart(Vec4Cart a, float s) {
    return MakeVec4Cart(a.t * s, a.x * s, a.y * s, a.z * s);
}

// --- Geodesic state ----------------------------------------------------------

// Position and 4-velocity in Boyer-Lindquist coordinates.
public struct GeodesicState {
    public Vec4 x;   // position (t, r, theta, phi)
    public Vec4 u;   // 4-velocity dx^mu/dlambda
}

// Position and 4-velocity in Kerr-Schild Cartesian coordinates. The render
// path integrates here to dodge the pole singularity (RDOP002A.cu
// integrateGeodesicRK4Cart).
public struct GeodesicStateCart {
    public Vec4Cart x;
    public Vec4Cart u;
}

// Canonical (q, p) state for the time-transformed symplectic integrator
// (RDOP002A.cu HamiltonianState). q = (t, r, theta, phi); p are conjugate
// momenta p_mu = g_mu_nu u^nu.
public struct HamiltonianState {
    public float q[4];
    public float p[4];
}

public HamiltonianState ZeroHamiltonian() {
    HamiltonianState h;
    [unroll] for (int i = 0; i < 4; i++) { h.q[i] = 0.0f; h.p[i] = 0.0f; }
    return h;
}

// --- Ray bundle (DNGR geodesic deviation, RDOP003A.h RayBundleState) ----------

// The central ray plus two deviation vectors defining the pixel's celestial
// cross-section. xi1/xi2 are the deviation vectors, dxi1/dxi2 their affine
// derivatives. float4 components map to (t, r, theta, phi).
public struct RayBundleState {
    public float4 x;
    public float4 u;
    public float4 xi1;
    public float4 dxi1;
    public float4 xi2;
    public float4 dxi2;
    public float semiMajor;
    public float semiMinor;
    public float orientation;
    public float affineParam;
    public bool terminated;
}

// --- Metric family parameters ------------------------------------------------

// Metric selector. The values match the render kernel's dispatch id used by
// trace/parity_probe; they are a kernel-local encoding, not the C++ MetricId
// authority (docs/ARCHITECTURE.md section 2 deletes the kernel's private
// MetricType enum in favour of a generated table, a later workstream).
public static const uint kMetricKerrSchild = 0;   // Kerr-Schild family (M, a, Q, Lambda)
public static const uint kMetricEllis      = 1;   // Ellis drainhole wormhole
public static const uint kMetricWarpDrive  = 2;   // Alcubierre / warp-drive family
public static const uint kMetricMorrisThorne = 3; // Morris-Thorne wormhole family

// Family parameters. The Kerr-Schild slot carries (M, a, Q, Lambda); the
// wormhole slot carries (m, n) or (b0, Phi0); the warp slot carries
// (vs, sigma, R). One flat struct keeps the probe ABI simple.
public struct MetricParamsGPU {
    public uint metricId;
    public float M;       // mass, or Ellis m
    public float a;       // spin, or Ellis n
    public float Q;       // charge, or wormhole redshift Phi0
    public float Lambda;  // cosmological constant, or warp wall thickness sigma
    public float p5;      // warp radius R
    public float p6;      // warp velocity vs
    public float p7;      // spare
}
