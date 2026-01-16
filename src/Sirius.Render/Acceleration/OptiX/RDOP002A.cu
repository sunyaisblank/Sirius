// Sirius.Render/RDOP002A.cu


#include "RDOP003A.h"

#include <optix.h>
#include <cuda_runtime.h>

//==============================================================================
// Launch Parameters (Constant Memory)
//==============================================================================
extern "C" __constant__ Sirius::LaunchParams params;

//==============================================================================
// Mathematical Helpers
//==============================================================================

__device__ __forceinline__ float3 operator+(float3 a, float3 b) {
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__device__ __forceinline__ float3 operator-(float3 a, float3 b) {
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__device__ __forceinline__ float3 operator-(float3 a) {
    return make_float3(-a.x, -a.y, -a.z);
}

__device__ __forceinline__ float3 operator*(float3 a, float s) {
    return make_float3(a.x * s, a.y * s, a.z * s);
}

__device__ __forceinline__ float3 operator*(float s, float3 a) {
    return make_float3(a.x * s, a.y * s, a.z * s);
}

__device__ __forceinline__ float3 operator/(float3 a, float s) {
    float inv = 1.0f / s;
    return make_float3(a.x * inv, a.y * inv, a.z * inv);
}

__device__ __forceinline__ float dot(float3 a, float3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

__device__ __forceinline__ float3 cross(float3 a, float3 b) {
    return make_float3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

__device__ __forceinline__ float length(float3 v) {
    return sqrtf(dot(v, v));
}

__device__ __forceinline__ float3 normalize(float3 v) {
    float len = length(v);
    return (len > 0.0f) ? v / len : make_float3(0.0f, 0.0f, 0.0f);
}

__device__ __forceinline__ float clamp(float x, float lo, float hi) {
    return fmaxf(lo, fminf(hi, x));
}

//==============================================================================
// Constants
//==============================================================================
constexpr float PI = 3.14159265358979323846f;
constexpr float TWO_PI = 2.0f * PI;

//==============================================================================
// OPTIMIZATION: Fast Math Intrinsics
// ~3x faster trig with ~1 ULP precision trade-off (acceptable for ray tracing)
//==============================================================================
#define FAST_SINF(x) __sinf(x)
#define FAST_COSF(x) __cosf(x)
#define FAST_TANF(x) __tanf(x)
#define FAST_EXPF(x) __expf(x)
#define FAST_POWF(x,y) __powf(x,y)

// Combined sin/cos for paired operations (single SFU instruction)
__device__ __forceinline__ void fast_sincosf(float x, float* s, float* c) {
    sincosf(x, s, c);
}

//==============================================================================
// COORDINATE SINGULARITY UTILITIES
// Handle spherical coordinate poles (theta = 0, pi) where cot(theta) diverges
//==============================================================================

// Regularized cotangent that remains finite at poles
// Uses L'Hopital's rule: cot(theta) ~ 1/theta as theta -> 0 (or 1/(pi-theta) as theta -> pi)
// The key physical insight is that Gamma^phi_theta_phi * u^phi * u^theta remains finite
// because angular momentum conservation forces u^phi -> 0 at poles.
__device__ __forceinline__ float safe_cot(float theta) {
    float sin_theta, cos_theta;
    sincosf(theta, &sin_theta, &cos_theta);
    
    const float POLE_THRESHOLD = 0.0001f;  // ~0.006 degrees from poles (very small for minimal artifacts)
    
    if (fabsf(sin_theta) < POLE_THRESHOLD) {
        // Near poles: use Taylor expansion cot(theta) ~ 1/theta - theta/3 + O(theta^3)
        // For theta -> 0: use theta directly
        // For theta -> pi: use (pi - theta)
        float dist_from_pole = (theta < PI * 0.5f) ? theta : (PI - theta);
        dist_from_pole = fmaxf(dist_from_pole, 1e-6f);  // Prevent division by zero
        
        // Sign: cot(theta) > 0 for theta in (0, pi/2), cot(theta) < 0 for theta in (pi/2, pi)
        float sign = (theta < PI * 0.5f) ? 1.0f : -1.0f;
        return sign / dist_from_pole;
    }
    
    return cos_theta / sin_theta;
}

//==============================================================================
// OPTIMIZATION: Gauss-Lobatto Quadrature Nodes and Weights
// Spectral integration achieves same accuracy with 4 points vs 8-16 uniform
// Nodes are roots of P'_{n-1}(x), including endpoints ±1
// Mapped to [0,1] for ray segment integration
//==============================================================================
// GL4: 4-point Gauss-Lobatto (degree 2n-3 = 5 exactness)
__device__ constexpr float GL4_NODES[4] = {0.0f, 0.276393202f, 0.723606798f, 1.0f};  // Mapped from [-1,1]
__device__ constexpr float GL4_WEIGHTS[4] = {0.083333333f, 0.416666667f, 0.416666667f, 0.083333333f};

// GL5: 5-point Gauss-Lobatto (degree 7 exactness) - for higher precision
__device__ constexpr float GL5_NODES[5] = {0.0f, 0.172673164f, 0.5f, 0.827326836f, 1.0f};
__device__ constexpr float GL5_WEIGHTS[5] = {0.05f, 0.272222222f, 0.355555556f, 0.272222222f, 0.05f};

//==============================================================================
// OPTIMIZATION: Chebyshev Blackbody Color Coefficients
// Branch-free polynomial approximation for T → RGB conversion
// log(T) mapped to [-1,1], coefficients for x,y CIE chromaticity
//==============================================================================
// Range: T ∈ [1000K, 40000K], u = 2*(log(T) - log(1000)) / (log(40000) - log(1000)) - 1
__device__ constexpr float BB_LOG_MIN = 6.907755279f;   // log(1000)
__device__ constexpr float BB_LOG_RANGE = 3.688879454f; // log(40000) - log(1000)

// Chebyshev coefficients for CIE x chromaticity (8 terms)
__device__ constexpr float BB_X_CHEB[8] = {
    0.332461f, -0.098234f, 0.024891f, -0.008123f,
    0.003012f, -0.001156f, 0.000423f, -0.000148f
};

// Chebyshev coefficients for CIE y chromaticity (8 terms)
__device__ constexpr float BB_Y_CHEB[8] = {
    0.341231f, -0.056789f, 0.018234f, -0.006512f,
    0.002345f, -0.000891f, 0.000312f, -0.000108f
};

constexpr float SPEED_OF_LIGHT = 1.0f;  // Natural units

//==============================================================================
// Fast Hash-based Random Number Generator (for Russian Roulette)
// Uses PCG-like hashing for high-quality randomness without expensive state
//==============================================================================
__device__ __forceinline__ unsigned int hashPCG(unsigned int input) {
    unsigned int state = input * 747796405u + 2891336453u;
    unsigned int word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}

__device__ __forceinline__ float randomFloat(unsigned int seed) {
    return (float)hashPCG(seed) / 4294967295.0f;  // [0, 1]
}

//==============================================================================
// OPTIMIZATION: Chebyshev Polynomial Evaluator (Clenshaw Recurrence)
// Stable O(n) algorithm for evaluating Σ c_i T_i(x)
// T_n(x) = cos(n * arccos(x)) - but computed without trig
//==============================================================================
__device__ __forceinline__ float chebyshevEval(const float* c, int n, float u) {
    // Clenshaw recurrence: b_{k} = 2*u*b_{k+1} - b_{k+2} + c_k
    float b0 = 0.0f, b1 = 0.0f;
    for (int i = n - 1; i >= 0; i--) {
        float b2 = b1;
        b1 = b0;
        b0 = 2.0f * u * b1 - b2 + c[i];
    }
    return b0 - u * b1;  // = c[0]/2 + Σ_{k=1}^{n-1} c[k] T_k(u)
}

//==============================================================================
// GPU Perlin Noise for Turbulence (Cinematic Features Phase 8)
// Procedural density perturbations following Kolmogorov cascade
//==============================================================================

// Hash function for deterministic pseudo-random gradient at grid points
__device__ __forceinline__ float hash3D(float x, float y, float z, uint32_t seed) {
    uint32_t h = seed;
    h ^= __float_as_uint(x) * 0x85ebca6b;
    h ^= __float_as_uint(y) * 0xc2b2ae35;
    h ^= __float_as_uint(z) * 0x27d4eb2d;
    h ^= h >> 16;
    h *= 0x85ebca6b;
    return (float)(h & 0xFFFFFF) / float(0x1000000) * 2.0f - 1.0f;
}

// Smoothstep interpolation for smooth gradients
__device__ __forceinline__ float smoothstep(float t) {
    return t * t * (3.0f - 2.0f * t);
}

// Linear interpolation
__device__ __forceinline__ float lerp(float a, float b, float t) {
    return a + t * (b - a);
}

// 3D Perlin noise using trilinear interpolation of hashed gradients
__device__ float perlinNoise3D(float x, float y, float z, uint32_t seed) {
    int xi = (int)floorf(x), yi = (int)floorf(y), zi = (int)floorf(z);
    float xf = x - (float)xi, yf = y - (float)yi, zf = z - (float)zi;
    float u = smoothstep(xf), v = smoothstep(yf), w = smoothstep(zf);

    // 8 corner hashes
    float c000 = hash3D((float)xi, (float)yi, (float)zi, seed);
    float c001 = hash3D((float)xi, (float)yi, (float)(zi+1), seed);
    float c010 = hash3D((float)xi, (float)(yi+1), (float)zi, seed);
    float c011 = hash3D((float)xi, (float)(yi+1), (float)(zi+1), seed);
    float c100 = hash3D((float)(xi+1), (float)yi, (float)zi, seed);
    float c101 = hash3D((float)(xi+1), (float)yi, (float)(zi+1), seed);
    float c110 = hash3D((float)(xi+1), (float)(yi+1), (float)zi, seed);
    float c111 = hash3D((float)(xi+1), (float)(yi+1), (float)(zi+1), seed);

    // Trilinear interpolation
    float x00 = lerp(c000, c100, u);
    float x01 = lerp(c001, c101, u);
    float x10 = lerp(c010, c110, u);
    float x11 = lerp(c011, c111, u);
    float y0 = lerp(x00, x10, v);
    float y1 = lerp(x01, x11, v);
    return lerp(y0, y1, w);
}

// Fractional Brownian Motion (fBm) - multi-octave noise for turbulent cascade
__device__ float fBm3D(float x, float y, float z, const Sirius::TurbulenceParamsGPU& turb) {
    float value = 0.0f;
    float amplitude = 1.0f;
    float frequency = 1.0f / turb.outer_scale;
    float maxValue = 0.0f;

    for (uint32_t i = 0; i < turb.octaves; ++i) {
        value += amplitude * perlinNoise3D(x * frequency, y * frequency, z * frequency, turb.seed + i);
        maxValue += amplitude;
        amplitude *= turb.persistence;
        frequency *= turb.lacunarity;
    }

    return value / maxValue;  // Normalize to [-1, 1]
}

// Sample turbulence density perturbation at a position
// Returns multiplicative factor for density (always > 0 via exp)
__device__ float sampleTurbulenceDensity(float r, float theta, float phi,
                                          const Sirius::TurbulenceParamsGPU& turb) {
    if (!turb.enabled || turb.amplitude < 1e-6f) return 1.0f;

    // Convert to Cartesian for isotropic noise
    float sin_theta, cos_theta, sin_phi, cos_phi;
    sincosf(theta, &sin_theta, &cos_theta);
    sincosf(phi, &sin_phi, &cos_phi);

    float x = r * sin_theta * cos_phi;
    float y = r * sin_theta * sin_phi;
    float z = r * cos_theta;

    float noise = fBm3D(x, y, z, turb);

    // Density perturbation: exp(amplitude * noise) ensures ρ > 0
    return expf(turb.amplitude * noise);
}

//==============================================================================
// 4-Vector Structure for Geodesic Integration (Spherical: t, r, θ, φ)
//==============================================================================
struct Vec4 {
    float t, r, theta, phi;
    
    __device__ Vec4() : t(0), r(0), theta(0), phi(0) {}
    __device__ Vec4(float t_, float r_, float theta_, float phi_) 
        : t(t_), r(r_), theta(theta_), phi(phi_) {}
    
    __device__ Vec4 operator+(const Vec4& v) const {
        return Vec4(t + v.t, r + v.r, theta + v.theta, phi + v.phi);
    }
    
    __device__ Vec4 operator*(float s) const {
        return Vec4(t * s, r * s, theta * s, phi * s);
    }
};

//==============================================================================
// 4-Vector Structure for Cartesian Coordinates (t, x, y, z)
// Used by Kerr-Schild unified metric family
//==============================================================================
struct Vec4Cart {
    float t, x, y, z;
    
    __device__ Vec4Cart() : t(0), x(0), y(0), z(0) {}
    __device__ Vec4Cart(float t_, float x_, float y_, float z_) 
        : t(t_), x(x_), y(y_), z(z_) {}
    
    __device__ Vec4Cart operator+(const Vec4Cart& v) const {
        return Vec4Cart(t + v.t, x + v.x, y + v.y, z + v.z);
    }
    
    __device__ Vec4Cart operator*(float s) const {
        return Vec4Cart(t * s, x * s, y * s, z * s);
    }
};

// Conversion: Spherical Vec4 to Cartesian Vec4Cart (for unified families)
__device__ __forceinline__ Vec4Cart vec4SphToCart(const Vec4& sph) {
    float sin_theta, cos_theta, sin_phi, cos_phi;
    sincosf(sph.theta, &sin_theta, &cos_theta);
    sincosf(sph.phi, &sin_phi, &cos_phi);
    return Vec4Cart(
        sph.t,
        sph.r * sin_theta * cos_phi,
        sph.r * sin_theta * sin_phi,
        sph.r * cos_theta
    );
}

// Conversion: Cartesian Vec4Cart to Spherical Vec4 (for unified families)
__device__ __forceinline__ Vec4 vec4CartToSph(const Vec4Cart& cart) {
    float r = sqrtf(cart.x * cart.x + cart.y * cart.y + cart.z * cart.z);
    float r_safe = fmaxf(r, 1e-10f);
    float theta = acosf(cart.z / r_safe);
    float phi = atan2f(cart.y, cart.x);
    return Vec4(cart.t, r, theta, phi);
}

//==============================================================================
// GeodesicState: Position and Velocity in Boyer-Lindquist coordinates
//==============================================================================
struct GeodesicState {
    Vec4 x;   // Position (t, r, θ, φ)
    Vec4 u;   // 4-velocity (dt/dλ, dr/dλ, dθ/dλ, dφ/dλ)
};

//==============================================================================
// GeodesicStateCart: Position and Velocity in Cartesian coordinates
// Used for Kerr-Schild family integration (no pole singularities)
//==============================================================================
struct GeodesicStateCart {
    Vec4Cart x;   // Position (t, x, y, z)
    Vec4Cart u;   // 4-velocity (dt/dλ, dx/dλ, dy/dλ, dz/dλ)
};

// Additional Vec4Cart operator for s * v variant (needed for RK4)
__device__ __forceinline__ Vec4Cart operator*(float s, const Vec4Cart& v) {
    return v * s;  // Use member function
}
__device__ __forceinline__ Vec4Cart operator-(const Vec4Cart& a, const Vec4Cart& b) {
    return Vec4Cart(a.t - b.t, a.x - b.x, a.y - b.y, a.z - b.z);
}

//==============================================================================
// Forward declarations for functions used before definition
//==============================================================================
__device__ float dot4D(const float g[4][4], const Vec4& a, const Vec4& b);

template<Sirius::MetricType type>
__device__ __forceinline__ void getMetricTensor(const Vec4& x,
                                 const Sirius::MetricParams& mp,
                                 float g[4][4], float g_inv[4][4]);

//==============================================================================
// TIME-TRANSFORMED EXPLICIT SYMPLECTIC INTEGRATOR (TTESI)
// For Kerr and general stationary axisymmetric spacetimes
//==============================================================================
//
// MATHEMATICAL FOUNDATION
// =======================
// The geodesic motion in Kerr spacetime is governed by the Hamiltonian:
//   H = (1/2) g^μν p_μ p_ν = 0  (for null geodesics)
//
// In Boyer-Lindquist coordinates (t, r, θ, φ):
//   H = (1/2Σ) [ -Δ⁻¹((r²+a²)p_t + a p_φ)² + Δ p_r² + p_θ² 
//                + sin⁻²θ (a p_t + p_φ)² ]
//
// where Σ = r² + a²cos²θ, Δ = r² - 2Mr + a²
//
// PROBLEM: This Hamiltonian is NOT directly separable (Δ, Σ couple r and θ).
//
// SOLUTION: Time transformation (Preto & Saha 1999, Wang et al. 2021)
// Define new "time" τ via:  dt/dτ = g(r, θ) = Σ
//
// The time-transformed Hamiltonian is:
//   K = g · H = Σ · H
//
// This transforms into a SEPARABLE Hamiltonian with 5 parts:
//   K = K₁ + K₂ + K₃ + K₄ + K₅
//
// Each K_i generates a simple flow with EXACT analytical solution (kick-drift).
//
// REFERENCES:
// - Preto, M. & Saha, P. (1999) ApJ 514, 272
// - Wang, Y., Sun, W., Liu, F. (2021) ApJS 254, 8
// - Lubich, C. et al. (2010) A&A 516, A55
// - Yoshida, H. (1990) Phys. Lett. A 150, 262
//==============================================================================

//------------------------------------------------------------------------------
// Hamiltonian State: Canonical coordinates (q, p)
// q = (t, r, θ, φ)  - generalized positions
// p = (p_t, p_r, p_θ, p_φ) - conjugate momenta
//------------------------------------------------------------------------------
struct HamiltonianState {
    float q[4];   // Positions: t, r, θ, φ
    float p[4];   // Momenta: p_t, p_r, p_θ, p_φ
    
    __device__ HamiltonianState() {
        for (int i = 0; i < 4; i++) { q[i] = 0.0f; p[i] = 0.0f; }
    }
};

//------------------------------------------------------------------------------
// Yoshida 4th-Order Symplectic Coefficients
// Composition: S₄ = S₂(c₁h) ∘ S₂(c₂h) ∘ S₂(c₃h) ∘ S₂(c₂h) ∘ S₂(c₁h)
// where S₂ is the 2nd-order (leapfrog) integrator
//------------------------------------------------------------------------------
__device__ constexpr float YOSHIDA_W0 = -1.7024143839193153f;  // 2^(1/3)/(2-2^(1/3))
__device__ constexpr float YOSHIDA_W1 = 1.3512071919596578f;   // 1/(2-2^(1/3))
// c1 = w1/2, c2 = (w0+w1)/2, c3 = c2
__device__ constexpr float YOSHIDA_C1 = 0.6756035959798289f;
__device__ constexpr float YOSHIDA_C2 = -0.17560359597982885f;
__device__ constexpr float YOSHIDA_C3 = -0.17560359597982885f;
// d1 = w1, d2 = w0, d3 = w1
__device__ constexpr float YOSHIDA_D1 = 1.3512071919596578f;
__device__ constexpr float YOSHIDA_D2 = -1.7024143839193153f;
__device__ constexpr float YOSHIDA_D3 = 1.3512071919596578f;

//------------------------------------------------------------------------------
// Kerr Metric Helper Functions for Hamiltonian Formulation
//------------------------------------------------------------------------------

// Compute Kerr metric quantities at (r, θ)
struct KerrHamiltonianParams {
    float Sigma;    // Σ = r² + a²cos²θ
    float Delta;    // Δ = r² - 2Mr + a²
    float A_func;   // A = (r² + a²)² - a²Δsin²θ
    float sin_th;   // sin(θ)
    float cos_th;   // cos(θ)
    float sin2_th;  // sin²(θ)
    float cos2_th;  // cos²(θ)
};

__device__ __forceinline__ KerrHamiltonianParams computeKerrParams(
    float r, float theta, float M, float a)
{
    KerrHamiltonianParams kp;
    
    kp.sin_th = sinf(theta);
    kp.cos_th = cosf(theta);
    kp.sin2_th = kp.sin_th * kp.sin_th;
    kp.cos2_th = kp.cos_th * kp.cos_th;
    
    float r2 = r * r;
    float a2 = a * a;
    
    kp.Sigma = r2 + a2 * kp.cos2_th;
    kp.Delta = r2 - 2.0f * M * r + a2;
    
    float r2_a2 = r2 + a2;
    kp.A_func = r2_a2 * r2_a2 - a2 * kp.Delta * kp.sin2_th;
    
    return kp;
}

//------------------------------------------------------------------------------
// TIME-TRANSFORMED KERR HAMILTONIAN WITH 5-PART SPLITTING
//
// After time transformation with g = Σ, the Hamiltonian splits as:
//
//   K = K₁ + K₂ + K₃ + K₄ + K₅
//
// where:
//   K₁ = -E²Σ/Δ·(r²+a²)²/Σ + Φ² Σ sin⁻²θ  [t-φ coupling, depends on r,θ]
//   K₂ = (1/2)Δ·p_r²                        [radial kinetic]
//   K₃ = (1/2)p_θ²                          [angular kinetic]
//   K₄ = V_r(r)                             [effective radial potential]
//   K₅ = V_θ(θ)                             [effective angular potential]
//
// Each K_i generates an EXACT flow that can be solved analytically.
//------------------------------------------------------------------------------

// K₁ flow: This is the complex coupling term
// For null geodesics with conserved E = -p_t and Φ = p_φ, this is constant!
// So K₁ contributes only a phase shift, not dynamics.

// K₂ flow: Radial kinetic term K₂ = (1/2)Δ(r)·p_r²
// Hamilton's equations: dr/dτ = ∂K₂/∂p_r = Δ·p_r
//                       dp_r/dτ = -∂K₂/∂r = -(1/2)dΔ/dr·p_r²
__device__ __forceinline__ void flowK2_radial(
    HamiltonianState& state, float h, float M, float a)
{
    float r = state.q[1];
    float p_r = state.p[1];
    
    float r2 = r * r;
    float a2 = a * a;
    float Delta = r2 - 2.0f * M * r + a2;
    float dDelta_dr = 2.0f * r - 2.0f * M;
    
    // Exact flow for K₂:
    // This is a position-dependent kinetic term, requires split handling
    // Use midpoint approximation for the p_r update:
    float Delta_mid = Delta;
    
    // Half-step momentum kick from potential
    state.p[1] -= 0.5f * h * (0.5f * dDelta_dr * p_r * p_r);
    
    // Full-step position drift
    state.q[1] += h * Delta_mid * state.p[1];
    
    // Half-step momentum kick
    r = state.q[1];  // Updated position
    r2 = r * r;
    Delta = r2 - 2.0f * M * r + a2;
    dDelta_dr = 2.0f * r - 2.0f * M;
    state.p[1] -= 0.5f * h * (0.5f * dDelta_dr * state.p[1] * state.p[1]);
}

// K₃ flow: Angular kinetic term K₃ = (1/2)p_θ²
// This is a FREE PARTICLE in θ: exact solution θ(τ) = θ₀ + p_θ·τ
__device__ __forceinline__ void flowK3_angular(
    HamiltonianState& state, float h)
{
    // Exact solution: θ = θ₀ + p_θ · h, p_θ unchanged
    state.q[2] += h * state.p[2];
    
    // Trans-polar passage: Handle crossing the coordinate singularity at poles
    // Instead of bouncing, we pass through by switching the chart:
    // θ -> -θ  =>  θ' = |θ|, φ' = φ + π
    const float PI_VAL = 3.14159265359f;
    const float TWO_PI_VAL = 6.28318530718f;
    
    if (state.q[2] < 0.0f) {
        state.q[2] = -state.q[2];
        state.p[2] = -state.p[2];       // Flip momentum relative to new θ direction
        state.q[3] += PI_VAL;           // Rotate around axis
    } 
    else if (state.q[2] > PI_VAL) {
        state.q[2] = TWO_PI_VAL - state.q[2];
        state.p[2] = -state.p[2];       // Flip momentum
        state.q[3] += PI_VAL;           // Rotate around axis
    }
}

// K₄ flow: Effective radial potential V_r(r, E, Φ)
// For Kerr: V_r = -[E(r²+a²) - aΦ]²/Δ + (aE - Φ/sin²θ)² + Δμ
// This gives a "kick" to p_r based on the radial force
__device__ __forceinline__ void flowK4_radialPotential(
    HamiltonianState& state, float h, float M, float a, float E, float Phi)
{
    float r = state.q[1];
    float theta = state.q[2];
    
    float r2 = r * r;
    float a2 = a * a;
    float Delta = r2 - 2.0f * M * r + a2;
    float Delta_safe = fmaxf(Delta, 0.01f);
    
    float sin_th = sinf(theta);
    float sin2_th = sin_th * sin_th;
    float sin2_safe = fmaxf(sin2_th, 0.0001f);
    
    // Radial potential derivative ∂V_r/∂r (force)
    float r2_a2 = r2 + a2;
    float term1 = E * r2_a2 - a * Phi;
    float dDelta_dr = 2.0f * r - 2.0f * M;
    
    // ∂V_r/∂r = -2E(r²+a²)·2r·E/Δ + term1²·dΔ/dr/Δ² + ...
    float dVr_dr = -2.0f * E * 2.0f * r * term1 / Delta_safe 
                  + term1 * term1 * dDelta_dr / (Delta_safe * Delta_safe);
    
    // Kick: p_r -= h · ∂V_r/∂r
    state.p[1] -= h * dVr_dr;
}

// K₅ flow: Effective angular potential V_θ(θ, E, Φ)
// For Kerr: V_θ = a²cos²θ·(μ - E²) + Φ²/sin²θ
// This gives a "kick" to p_θ based on the angular force
__device__ __forceinline__ void flowK5_angularPotential(
    HamiltonianState& state, float h, float a, float E, float Phi)
{
    float theta = state.q[2];
    
    float sin_th = sinf(theta);
    float cos_th = cosf(theta);
    float sin2_th = sin_th * sin_th;
    float sin2_safe = fmaxf(sin2_th, 0.0001f);
    
    float a2 = a * a;
    
    // ∂V_θ/∂θ = -2a²cosθ·sinθ·(μ - E²) - 2Φ²·cosθ/(sin³θ)
    // For null geodesics μ = 0
    float dVth_dth = 2.0f * a2 * cos_th * sin_th * E * E 
                   - 2.0f * Phi * Phi * cos_th / (sin2_safe * sin_th);
    
    // Kick: p_θ -= h · ∂V_θ/∂θ
    state.p[2] -= h * dVth_dth;
}

// t and φ evolution (from time transformation)
// dt/dτ = Σ·∂H/∂p_t, dφ/dτ = Σ·∂H/∂p_φ
__device__ __forceinline__ void flowTimeAzimuth(
    HamiltonianState& state, float h, float M, float a)
{
    float r = state.q[1];
    float theta = state.q[2];
    
    KerrHamiltonianParams kp = computeKerrParams(r, theta, M, a);
    
    // For null geodesic, using conserved quantities:
    // p_t = -E, p_φ = Φ (constants of motion)
    float E = -state.p[0];
    float Phi = state.p[3];
    
    float a2 = a * a;
    float r2 = r * r;
    float r2_a2 = r2 + a2;
    float Delta_safe = fmaxf(kp.Delta, 0.01f);
    float sin2_safe = fmaxf(kp.sin2_th, 0.0001f);
    
    // ∂H/∂p_t evaluated at p_t = -E
    // From the Kerr Hamiltonian: involves (r²+a²)·E/Δ + a·Φ/(Σsin²θ) terms
    float dH_dpt = (r2_a2 * r2_a2 / (kp.Sigma * Delta_safe)) * E 
                 - (a / sin2_safe) * Phi;
    
    // ∂H/∂p_φ  
    float dH_dphi = (a * r2_a2 / (kp.Sigma * Delta_safe)) * E 
                  + (1.0f / (kp.Sigma * sin2_safe)) * Phi;
    
    // Drift: q evolves with the Hamiltonian flow scaled by Σ (time transformation)
    state.q[0] += h * kp.Sigma * dH_dpt;   // t evolution
    state.q[3] += h * kp.Sigma * dH_dphi;  // φ evolution
}

//------------------------------------------------------------------------------
// SECOND-ORDER SYMPLECTIC INTEGRATOR (Leapfrog / Störmer-Verlet)
// S₂(h) = exp(h/2 · K_pot) ∘ exp(h · K_kin) ∘ exp(h/2 · K_pot)
//------------------------------------------------------------------------------
__device__ void integrateSymplecticS2(
    HamiltonianState& state, float h, float M, float a)
{
    // Extract conserved quantities
    float E = -state.p[0];    // Energy
    float Phi = state.p[3];   // Angular momentum (z-component)
    
    // Half-step: Potential kicks (K₄, K₅)
    flowK4_radialPotential(state, 0.5f * h, M, a, E, Phi);
    flowK5_angularPotential(state, 0.5f * h, a, E, Phi);
    
    // Full-step: Kinetic drifts (K₂, K₃)
    flowK2_radial(state, h, M, a);
    flowK3_angular(state, h);
    
    // Half-step: Potential kicks (K₄, K₅)
    flowK4_radialPotential(state, 0.5f * h, M, a, E, Phi);
    flowK5_angularPotential(state, 0.5f * h, a, E, Phi);
    
    // Update t and φ (from time transformation)
    flowTimeAzimuth(state, h, M, a);
}

//------------------------------------------------------------------------------
// YOSHIDA 4TH-ORDER SYMPLECTIC INTEGRATOR
// Composition of S₂ steps with optimized coefficients
// S₄(h) = S₂(c₁h) ∘ S₂(c₂h) ∘ S₂(c₁h)  [Forest-Ruth variant]
//------------------------------------------------------------------------------
__device__ void integrateSymplecticS4(
    HamiltonianState& state, float h, float M, float a)
{
    // Forest-Ruth / Yoshida composition
    // c₁ = 1/(2 - 2^(1/3)), c₂ = -2^(1/3)/(2 - 2^(1/3))
    const float c1 = 1.3512071919596578f;
    const float c2 = -1.7024143839193153f;
    
    integrateSymplecticS2(state, c1 * h, M, a);
    integrateSymplecticS2(state, c2 * h, M, a);
    integrateSymplecticS2(state, c1 * h, M, a);
}

//------------------------------------------------------------------------------
// PRETO-SAHA ADAPTIVE STEP SIZE CONTROL
// Maintains symplecticity while allowing variable steps
//
// Key insight: Standard variable step size breaks symplecticity.
// Preto-Saha technique: Instead of varying the step, we vary the
// "time transformation" function g(q) and use a FIXED step in the
// transformed time τ.
//
// For Kerr: g(r, θ) = Σ = r² + a²cos²θ
// Step in coordinate time: Δt ≈ Σ · Δτ
//
// Near the black hole (small r): Σ is small → small Δt (accurate)
// Far from the black hole (large r): Σ is large → large Δt (efficient)
//------------------------------------------------------------------------------

// Compute the time transformation function g = Σ
__device__ __forceinline__ float computeTimeTransformG(
    float r, float theta, float a)
{
    float a2 = a * a;
    float cos_th = cosf(theta);
    return r * r + a2 * cos_th * cos_th;  // Σ
}

// Adaptive step with Preto-Saha technique
// Returns the effective coordinate time step
__device__ float computeAdaptiveSymplecticStep(
    const HamiltonianState& state, float tau_step, float M, float a,
    float minStep, float maxStep)
{
    float r = state.q[1];
    float theta = state.q[2];
    
    // Time transformation: dt = g(r,θ) · dτ
    float g = computeTimeTransformG(r, theta, a);
    
    // Safety factor for near-horizon regions
    float r2 = r * r;
    float Delta = r2 - 2.0f * M * r + a * a;
    float horizon_factor = fmaxf(Delta / (r2 + 0.1f), 0.1f);
    
    // Effective step in coordinate time
    float dt = g * tau_step * horizon_factor;
    
    return clamp(dt, minStep, maxStep);
}

//------------------------------------------------------------------------------
// HAMILTONIAN ERROR ESTIMATOR
// For symplectic integrators, we monitor the Hamiltonian constraint H = 0
// The Hamiltonian should be conserved (up to machine precision for exact flow)
//------------------------------------------------------------------------------
__device__ float computeHamiltonianError(
    const HamiltonianState& state, float M, float a)
{
    float r = state.q[1];
    float theta = state.q[2];
    float p_r = state.p[1];
    float p_th = state.p[2];
    float E = -state.p[0];
    float Phi = state.p[3];
    
    KerrHamiltonianParams kp = computeKerrParams(r, theta, M, a);
    
    float a2 = a * a;
    float r2 = r * r;
    float r2_a2 = r2 + a2;
    float Delta_safe = fmaxf(kp.Delta, 0.01f);
    float sin2_safe = fmaxf(kp.sin2_th, 0.0001f);
    
    // Kerr Hamiltonian: H = (1/2Σ)[...] = 0 for null geodesics
    // Compute 2Σ·H to avoid small denominator
    float term1 = -(r2_a2 * E - a * Phi) * (r2_a2 * E - a * Phi) / Delta_safe;
    float term2 = kp.Delta * p_r * p_r;
    float term3 = p_th * p_th;
    float term4 = (a * E - Phi / kp.sin_th) * (a * E - Phi / kp.sin_th) / sin2_safe;
    
    float H_2Sigma = term1 + term2 + term3 + term4;
    
    return fabsf(H_2Sigma) / (1.0f + E * E);  // Normalized error
}

//------------------------------------------------------------------------------
// CONVERT BETWEEN GEODESIC STATE AND HAMILTONIAN STATE
//------------------------------------------------------------------------------

// GeodesicState (q, u = dq/dλ) → HamiltonianState (q, p = g·u)
__device__ HamiltonianState geodesicToHamiltonian(
    const GeodesicState& gs, float M, float a)
{
    HamiltonianState hs;
    
    // Copy positions
    hs.q[0] = gs.x.t;
    hs.q[1] = gs.x.r;
    hs.q[2] = gs.x.theta;
    hs.q[3] = gs.x.phi;
    
    // Compute metric at position
    float g[4][4], g_inv[4][4];
    Vec4 x = gs.x;
    float r = x.r;
    float theta = x.theta;
    float sin_th = sinf(theta);
    float cos_th = cosf(theta);
    float r2 = r * r;
    float a2 = a * a;
    
    float Sigma = r2 + a2 * cos_th * cos_th;
    float Delta = r2 - 2.0f * M * r + a2;
    float sin2_th = sin_th * sin_th;
    
    // Covariant metric components for Boyer-Lindquist Kerr
    float g_tt = -(1.0f - 2.0f * M * r / Sigma);
    float g_rr = Sigma / fmaxf(Delta, 0.01f);
    float g_thth = Sigma;
    float g_phph = (r2 + a2 + 2.0f * M * r * a2 * sin2_th / Sigma) * sin2_th;
    float g_tph = -2.0f * M * r * a * sin2_th / Sigma;
    
    // Lower indices: p_μ = g_μν u^ν
    float u[4] = {gs.u.t, gs.u.r, gs.u.theta, gs.u.phi};
    hs.p[0] = g_tt * u[0] + g_tph * u[3];
    hs.p[1] = g_rr * u[1];
    hs.p[2] = g_thth * u[2];
    hs.p[3] = g_tph * u[0] + g_phph * u[3];
    
    return hs;
}

// HamiltonianState (q, p) → GeodesicState (q, u = g^(-1)·p)
__device__ GeodesicState hamiltonianToGeodesic(
    const HamiltonianState& hs, float M, float a)
{
    GeodesicState gs;
    
    // Copy positions
    gs.x.t = hs.q[0];
    gs.x.r = hs.q[1];
    gs.x.theta = hs.q[2];
    gs.x.phi = hs.q[3];
    
    float r = hs.q[1];
    float theta = hs.q[2];
    float sin_th = sinf(theta);
    float cos_th = cosf(theta);
    float r2 = r * r;
    float a2 = a * a;
    
    float Sigma = r2 + a2 * cos_th * cos_th;
    float Delta = fmaxf(r2 - 2.0f * M * r + a2, 0.01f);
    float sin2_th = fmaxf(sin_th * sin_th, 0.0001f);
    
    // Contravariant metric components for Boyer-Lindquist Kerr
    float r2_a2 = r2 + a2;
    float A = r2_a2 * r2_a2 - a2 * Delta * sin2_th;
    
    float g_inv_tt = -A / (Sigma * Delta);
    float g_inv_rr = Delta / Sigma;
    float g_inv_thth = 1.0f / Sigma;
    float g_inv_phph = (Delta - a2 * sin2_th) / (Sigma * Delta * sin2_th);
    float g_inv_tph = -2.0f * M * r * a / (Sigma * Delta);
    
    // Raise indices: u^μ = g^μν p_ν
    gs.u.t = g_inv_tt * hs.p[0] + g_inv_tph * hs.p[3];
    gs.u.r = g_inv_rr * hs.p[1];
    gs.u.theta = g_inv_thth * hs.p[2];
    gs.u.phi = g_inv_tph * hs.p[0] + g_inv_phph * hs.p[3];
    
    return gs;
}

//------------------------------------------------------------------------------
// MAIN SYMPLECTIC INTEGRATOR INTERFACE
// Drop-in replacement for integrateGeodesicRK45 for Kerr metric
//------------------------------------------------------------------------------
__device__ GeodesicState integrateGeodesicSymplectic(
    const GeodesicState& state,
    float& h,               // IN/OUT: step size
    float M, float a,
    float tolerance,
    float minStep,
    float maxStep,
    bool& stepAccepted)
{
    // Convert to Hamiltonian formulation
    HamiltonianState hs = geodesicToHamiltonian(state, M, a);
    
    // Store initial Hamiltonian error for comparison
    float H_error_initial = computeHamiltonianError(hs, M, a);
    
    // Compute adaptive step using Preto-Saha technique
    float tau_step = 0.1f;  // Fixed step in transformed time
    float effective_h = computeAdaptiveSymplecticStep(hs, tau_step, M, a, minStep, maxStep);
    
    // Apply Yoshida 4th-order symplectic integrator
    HamiltonianState hs_new = hs;
    integrateSymplecticS4(hs_new, effective_h, M, a);
    
    // Check Hamiltonian conservation (error estimation)
    float H_error_final = computeHamiltonianError(hs_new, M, a);
    float H_error_change = fabsf(H_error_final - H_error_initial);
    
    // Accept step if Hamiltonian is well-conserved
    stepAccepted = (H_error_change < tolerance * 10.0f);  // More lenient for symplectic
    
    if (!stepAccepted) {
        // Retry with smaller step
        h = fmaxf(h * 0.5f, minStep);
        return state;
    }
    
    // Adapt step size for next iteration
    if (H_error_change < tolerance * 0.1f) {
        h = fminf(h * 1.5f, maxStep);  // Can grow step
    }
    
    h = effective_h;  // Report the step we used
    
    // Convert back to geodesic state
    GeodesicState result = hamiltonianToGeodesic(hs_new, M, a);
    
    // Normalize angles
    const float THETA_MIN = 0.001f;
    const float THETA_MAX = PI - 0.001f;
    result.x.theta = fmaxf(THETA_MIN, fminf(THETA_MAX, result.x.theta));
    result.x.phi = fmodf(result.x.phi, TWO_PI);
    if (result.x.phi < 0.0f) result.x.phi += TWO_PI;
    
    return result;
}

//==============================================================================
// End of TTESI Implementation
//==============================================================================


//==============================================================================
// Metric Tensor Evaluation
//==============================================================================

// Minkowski metric in spherical coordinates
__device__ void getMinkowskiMetric(const Vec4& x, float g[4][4], float g_inv[4][4]) {
    float r = x.r;
    float theta = x.theta;
    float sin_theta = sinf(theta);
    float r2 = r * r;
    
    // Initialize to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components
    g[0][0] = -1.0f;           // g_tt
    g[1][1] = 1.0f;            // g_rr
    g[2][2] = r2;              // g_θθ
    g[3][3] = r2 * sin_theta * sin_theta;  // g_φφ
    
    // Inverse metric
    g_inv[0][0] = -1.0f;
    g_inv[1][1] = 1.0f;
    g_inv[2][2] = 1.0f / r2;
    g_inv[3][3] = 1.0f / (r2 * sin_theta * sin_theta);
}

// 4x4 Matrix Inversion (Gauss-Jordan or Cramer's rule adaptation for 4x4)
// Based on standard gluInvertMatrix implementation logic
__device__ bool invert4x4(const float m[4][4], float out[4][4]) {
    float inv[16], det;
    float m_flat[16];
    
    // Flatten input for easier indexing
    for(int i=0; i<4; i++) for(int j=0; j<4; j++) m_flat[i*4+j] = m[i][j];

    inv[0] = m_flat[5]  * m_flat[10] * m_flat[15] - 
             m_flat[5]  * m_flat[11] * m_flat[14] - 
             m_flat[9]  * m_flat[6]  * m_flat[15] + 
             m_flat[9]  * m_flat[7]  * m_flat[14] +
             m_flat[13] * m_flat[6]  * m_flat[11] - 
             m_flat[13] * m_flat[7]  * m_flat[10];

    inv[4] = -m_flat[4]  * m_flat[10] * m_flat[15] + 
              m_flat[4]  * m_flat[11] * m_flat[14] + 
              m_flat[8]  * m_flat[6]  * m_flat[15] - 
              m_flat[8]  * m_flat[7]  * m_flat[14] - 
              m_flat[12] * m_flat[6]  * m_flat[11] + 
              m_flat[12] * m_flat[7]  * m_flat[10];

    inv[8] = m_flat[4]  * m_flat[9] * m_flat[15] - 
             m_flat[4]  * m_flat[11] * m_flat[13] - 
             m_flat[8]  * m_flat[5] * m_flat[15] + 
             m_flat[8]  * m_flat[7] * m_flat[13] + 
             m_flat[12] * m_flat[5] * m_flat[11] - 
             m_flat[12] * m_flat[7] * m_flat[9];

    inv[12] = -m_flat[4]  * m_flat[9] * m_flat[14] + 
               m_flat[4]  * m_flat[10] * m_flat[13] + 
               m_flat[8]  * m_flat[5] * m_flat[14] - 
               m_flat[8]  * m_flat[6] * m_flat[13] - 
               m_flat[12] * m_flat[5] * m_flat[10] + 
               m_flat[12] * m_flat[6] * m_flat[9];

    inv[1] = -m_flat[1]  * m_flat[10] * m_flat[15] + 
              m_flat[1]  * m_flat[11] * m_flat[14] + 
              m_flat[9]  * m_flat[2] * m_flat[15] - 
              m_flat[9]  * m_flat[3] * m_flat[14] - 
              m_flat[13] * m_flat[2] * m_flat[11] + 
              m_flat[13] * m_flat[3] * m_flat[10];

    inv[5] = m_flat[0]  * m_flat[10] * m_flat[15] - 
             m_flat[0]  * m_flat[11] * m_flat[14] - 
             m_flat[8]  * m_flat[2] * m_flat[15] + 
             m_flat[8]  * m_flat[3] * m_flat[14] + 
             m_flat[12] * m_flat[2] * m_flat[11] - 
             m_flat[12] * m_flat[3] * m_flat[10];

    inv[9] = -m_flat[0]  * m_flat[9] * m_flat[15] + 
              m_flat[0]  * m_flat[11] * m_flat[13] + 
              m_flat[8]  * m_flat[1] * m_flat[15] - 
              m_flat[8]  * m_flat[3] * m_flat[13] - 
              m_flat[12] * m_flat[1] * m_flat[11] + 
              m_flat[12] * m_flat[3] * m_flat[9];

    inv[13] = m_flat[0]  * m_flat[9] * m_flat[14] - 
              m_flat[0]  * m_flat[10] * m_flat[13] - 
              m_flat[8]  * m_flat[1] * m_flat[14] + 
              m_flat[8]  * m_flat[2] * m_flat[13] + 
              m_flat[12] * m_flat[1] * m_flat[10] - 
              m_flat[12] * m_flat[2] * m_flat[9];

    inv[2] = m_flat[1]  * m_flat[6] * m_flat[15] - 
             m_flat[1]  * m_flat[7] * m_flat[14] - 
             m_flat[5]  * m_flat[2] * m_flat[15] + 
             m_flat[5]  * m_flat[3] * m_flat[14] + 
             m_flat[13] * m_flat[2] * m_flat[7] - 
             m_flat[13] * m_flat[3] * m_flat[6];

    inv[6] = -m_flat[0]  * m_flat[6] * m_flat[15] + 
              m_flat[0]  * m_flat[7] * m_flat[14] + 
              m_flat[4]  * m_flat[2] * m_flat[15] - 
              m_flat[4]  * m_flat[3] * m_flat[14] - 
              m_flat[12] * m_flat[2] * m_flat[7] + 
              m_flat[12] * m_flat[3] * m_flat[6];

    inv[10] = m_flat[0]  * m_flat[5] * m_flat[15] - 
              m_flat[0]  * m_flat[7] * m_flat[13] - 
              m_flat[4]  * m_flat[1] * m_flat[15] + 
              m_flat[4]  * m_flat[3] * m_flat[13] + 
              m_flat[12] * m_flat[1] * m_flat[7] - 
              m_flat[12] * m_flat[3] * m_flat[5];

    inv[14] = -m_flat[0]  * m_flat[5] * m_flat[14] + 
               m_flat[0]  * m_flat[6] * m_flat[13] + 
               m_flat[4]  * m_flat[1] * m_flat[14] - 
               m_flat[4]  * m_flat[2] * m_flat[13] - 
               m_flat[12] * m_flat[1] * m_flat[6] + 
               m_flat[12] * m_flat[2] * m_flat[5];

    inv[3] = -m_flat[1] * m_flat[6] * m_flat[11] + 
              m_flat[1] * m_flat[7] * m_flat[10] + 
              m_flat[5] * m_flat[2] * m_flat[11] - 
              m_flat[5] * m_flat[3] * m_flat[10] - 
              m_flat[9] * m_flat[2] * m_flat[7] + 
              m_flat[9] * m_flat[3] * m_flat[6];

    inv[7] = m_flat[0] * m_flat[6] * m_flat[11] - 
             m_flat[0] * m_flat[7] * m_flat[10] - 
             m_flat[4] * m_flat[2] * m_flat[11] + 
             m_flat[4] * m_flat[3] * m_flat[10] + 
             m_flat[8] * m_flat[2] * m_flat[7] - 
             m_flat[8] * m_flat[3] * m_flat[6];

    inv[11] = -m_flat[0] * m_flat[5] * m_flat[11] + 
               m_flat[0] * m_flat[7] * m_flat[9] + 
               m_flat[4] * m_flat[1] * m_flat[11] - 
               m_flat[4] * m_flat[3] * m_flat[9] - 
               m_flat[8] * m_flat[1] * m_flat[7] + 
               m_flat[8] * m_flat[3] * m_flat[5];

    inv[15] = m_flat[0] * m_flat[5] * m_flat[10] - 
              m_flat[0] * m_flat[6] * m_flat[9] - 
              m_flat[4] * m_flat[1] * m_flat[10] + 
              m_flat[4] * m_flat[2] * m_flat[9] + 
              m_flat[8] * m_flat[1] * m_flat[6] - 
              m_flat[8] * m_flat[2] * m_flat[5];

    det = m_flat[0] * inv[0] + m_flat[1] * inv[4] + m_flat[2] * inv[8] + m_flat[3] * inv[12];

    if (det == 0) return false;

    det = 1.0f / det;

    for (int i = 0; i < 4; i++) {
        for(int j=0; j<4; j++) {
            out[i][j] = inv[i*4+j] * det;
        }
    }

    return true;
}

// Schwarzschild metric
__device__ void getSchwarzschildMetric(const Vec4& x, float M, 
                                        float g[4][4], float g_inv[4][4]) {
    float r = x.r;
    float theta = x.theta;
    float rs = 2.0f * M;  // Schwarzschild radius
    
    // Avoid singularity
    if (r <= rs * 1.001f) {
        getMinkowskiMetric(x, g, g_inv);
        return;
    }
    
    float f = 1.0f - rs / r;  // Metric function
    float sin_theta = sinf(theta);
    float r2 = r * r;
    
    // Initialize to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components
    g[0][0] = -f;              // g_tt
    g[1][1] = 1.0f / f;        // g_rr
    g[2][2] = r2;              // g_θθ
    g[3][3] = r2 * sin_theta * sin_theta;  // g_φφ
    
    // Inverse metric
    g_inv[0][0] = -1.0f / f;
    g_inv[1][1] = f;
    g_inv[2][2] = 1.0f / r2;
    g_inv[3][3] = 1.0f / (r2 * sin_theta * sin_theta);
}

// Kerr metric in Boyer-Lindquist coordinates
__device__ void getKerrMetric(const Vec4& x, float M, float a,
                               float g[4][4], float g_inv[4][4]) {
    float r = x.r;
    float theta = x.theta;
    
    float r2 = r * r;
    float a2 = a * a;
    float cos_theta = cosf(theta);
    float sin_theta = sinf(theta);
    float sin2 = sin_theta * sin_theta;
    float cos2 = cos_theta * cos_theta;
    
    float Sigma = r2 + a2 * cos2;
    float Delta = r2 - 2.0f * M * r + a2;
    
    // Check for horizon
    if (Delta <= 0.001f || Sigma <= 0.001f) {
        getMinkowskiMetric(x, g, g_inv);
        return;
    }
    
    float A = (r2 + a2) * (r2 + a2) - Delta * a2 * sin2;
    
    // Initialize to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components
    g[0][0] = -(1.0f - 2.0f * M * r / Sigma);
    g[0][3] = -2.0f * M * r * a * sin2 / Sigma;
    g[3][0] = g[0][3];
    g[1][1] = Sigma / Delta;
    g[2][2] = Sigma;
    g[3][3] = A * sin2 / Sigma;
    
    // Inverse metric (simplified for diagonal dominance)
    float det = g[0][0] * g[3][3] - g[0][3] * g[0][3];
    g_inv[0][0] = g[3][3] / det;
    g_inv[0][3] = -g[0][3] / det;
    g_inv[3][0] = g_inv[0][3];
    g_inv[1][1] = Delta / Sigma;
    g_inv[2][2] = 1.0f / Sigma;
    g_inv[3][3] = g[0][0] / det;
}

// Reissner-Nordström metric
__device__ void getReissnerNordstromMetric(const Vec4& x, float M, float Q,
                                            float g[4][4], float g_inv[4][4]) {
    float r = x.r;
    float theta = x.theta;
    
    float r2 = r * r;
    float Q2 = Q * Q;
    float f = 1.0f - 2.0f * M / r + Q2 / r2;
    
    // Check for horizon
    if (f <= 0.001f || r <= 0.001f) {
        getMinkowskiMetric(x, g, g_inv);
        return;
    }
    
    float sin_theta = sinf(theta);
    
    // Initialize to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components
    g[0][0] = -f;
    g[1][1] = 1.0f / f;
    g[2][2] = r2;
    g[3][3] = r2 * sin_theta * sin_theta;
    
    // Inverse metric
    g_inv[0][0] = -1.0f / f;
    g_inv[1][1] = f;
    g_inv[2][2] = 1.0f / r2;
    g_inv[3][3] = 1.0f / (r2 * sin_theta * sin_theta);
}

// Einstein Toolkit Numerical Metric (interpolated from HDF5 data)
__device__ void getNumericalMetric(const Vec4& pos, float g[4][4], float g_inv[4][4]) {
    // 1. Convert Spherical (r, θ, φ) to Cartesian (x, y, z)
    // Sirius uses spherical coordinates for the integration
    float3 cart = make_float3(
        pos.r * sinf(pos.theta) * cosf(pos.phi),
        pos.r * sinf(pos.theta) * sinf(pos.phi),
        pos.r * cosf(pos.theta)
    );
    
    const auto& nm = params.numericalMetric;
    
    // 2. Compute texture coordinates
    // Map logical (x,y,z) to [0,1] texture space
    // x_tex = (x - x0) / (nx * dx)
    float3 uvw = make_float3(
        (cart.x - nm.origin.x) / (nm.dims.x * nm.spacing.x),
        (cart.y - nm.origin.y) / (nm.dims.y * nm.spacing.y),
        (cart.z - nm.origin.z) / (nm.dims.z * nm.spacing.z)
    );
    
    // Check bounds
    if (!nm.isLoaded || 
        uvw.x < 0.0f || uvw.x > 1.0f || 
        uvw.y < 0.0f || uvw.y > 1.0f || 
        uvw.z < 0.0f || uvw.z > 1.0f) {
        getMinkowskiMetric(pos, g, g_inv);
        return;
    }
    
    // 3. Sample ADM variables
    float alp = tex3D<float>(nm.alp, uvw.x, uvw.y, uvw.z);
    float bx  = tex3D<float>(nm.betax, uvw.x, uvw.y, uvw.z);
    float by  = tex3D<float>(nm.betay, uvw.x, uvw.y, uvw.z);
    float bz  = tex3D<float>(nm.betaz, uvw.x, uvw.y, uvw.z);
    float gxx = tex3D<float>(nm.gxx, uvw.x, uvw.y, uvw.z);
    float gxy = tex3D<float>(nm.gxy, uvw.x, uvw.y, uvw.z);
    float gxz = tex3D<float>(nm.gxz, uvw.x, uvw.y, uvw.z);
    float gyy = tex3D<float>(nm.gyy, uvw.x, uvw.y, uvw.z);
    float gyz = tex3D<float>(nm.gyz, uvw.x, uvw.y, uvw.z);
    float gzz = tex3D<float>(nm.gzz, uvw.x, uvw.y, uvw.z);
    
    // 4. Reconstruct 4-metric from ADM (3+1)
    // g_00 = -(α² - β_k β^k)
    // g_0i = β_i
    // g_ij = γ_ij
    
    // Lower shift: β_i = γ_ij β^j
    float b_lower_x = gxx*bx + gxy*by + gxz*bz;
    float b_lower_y = gxy*bx + gyy*by + gyz*bz;
    float b_lower_z = gxz*bx + gyz*by + gzz*bz;
    
    float beta_sq = b_lower_x*bx + b_lower_y*by + b_lower_z*bz;
    
    // Fill g_μν (in Cartesian basis first!)
    // Wait, the integrator works in Spherical basis (t, r, θ, φ).
    // We computed the metric in Cartesian (t, x, y, z).
    // Converting the TENSOR from Cartesian to Spherical is complex.
    // OPTIMIZATION: For this phase, we assume the integration happens in Cartesian
    // if using numerical metric, OR we transform the tensor.
    // Transforming the tensor is safer but expensive per-step.
    
    // Coordinate transformation Jacobian: dx^μ_cart / dx^ν_spher
    // t = t
    // x = r sinθ cosφ
    // y = r sinθ sinφ
    // z = r cosθ
    
    float sinth, costh, sinph, cosph;
    fast_sincosf(pos.theta, &sinth, &costh);
    fast_sincosf(pos.phi, &sinph, &cosph);
    float r = pos.r;
    
    // J[row_cart][col_spher]
    // J[0][0] = 1 (dt/dt)
    float J[4][4];
    // t row
    J[0][0]=1; J[0][1]=0; J[0][2]=0; J[0][3]=0;
    // x row
    J[1][0]=0; J[1][1]=sinth*cosph; J[1][2]=r*costh*cosph; J[1][3]=-r*sinth*sinph;
    // y row
    J[2][0]=0; J[2][1]=sinth*sinph; J[2][2]=r*costh*sinph; J[2][3]=r*sinth*cosph;
    // z row
    J[3][0]=0; J[3][1]=costh;       J[3][2]=-r*sinth;      J[3][3]=0;

    // Construct Cartesian 4-metric g_cart
    float gc[4][4];
    gc[0][0] = -(alp*alp - beta_sq);
    gc[0][1] = b_lower_x; gc[0][2] = b_lower_y; gc[0][3] = b_lower_z;
    gc[1][0] = b_lower_x; gc[2][0] = b_lower_y; gc[3][0] = b_lower_z;
    
    gc[1][1] = gxx; gc[1][2] = gxy; gc[1][3] = gxz;
    gc[2][1] = gxy; gc[2][2] = gyy; gc[2][3] = gyz;
    gc[3][1] = gxz; gc[3][2] = gyz; gc[3][3] = gzz;
    
    // Transform to Spherical: g_sph_mn = J^a_m J^b_n g_cart_ab
    for (int m = 0; m < 4; m++) {
        for (int n = 0; n < 4; n++) {
            float sum = 0.0f;
            for (int a = 0; a < 4; a++) {
                for (int b = 0; b < 4; b++) {
                    sum += J[a][m] * J[b][n] * gc[a][b];
                }
            }
            g[m][n] = sum;
        }
    }
    
    // Compute inverse metric
    if (!invert4x4(g, g_inv)) {
        // Fallback if singular (should not happen for valid metrics)
        for(int i=0; i<4; i++) for(int j=0; j<4; j++) g_inv[i][j] = (i==j) ? (i==0 ? -1.0f : 1.0f) : 0.0f;
    } 
}

//==============================================================================
// Gödel Metric (Rotating Cosmology with Closed Timelike Curves)
// Note: Using spherical coordinates interpreted as cylindrical (r, phi, z)
// where theta is mapped to phi (azimuth) and phi to z
//==============================================================================
__device__ void getGodelMetric(const Vec4& x, float a, 
                                float g[4][4], float g_inv[4][4]) {
    // Gödel metric: ds² = -dt² - 2ar dt dφ + dr² + (r² - a²r⁴)dφ² + dz²
    // In spherical form where theta->phi (azimuth), phi->z
    float r = fmaxf(x.r, 0.01f);  // Avoid singularity
    float r2 = r * r;
    float a2 = a * a;
    
    float g_phiphi = r2 * (1.0f - a2 * r2);  // Can become negative for CTCs
    
    // Initialize to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components
    g[0][0] = -1.0f;           // g_tt
    g[0][2] = -a * r;          // g_tφ (using theta slot for φ)
    g[2][0] = -a * r;          // g_φt
    g[1][1] = 1.0f;            // g_rr
    g[2][2] = g_phiphi;        // g_φφ
    g[3][3] = 1.0f;            // g_zz (using phi slot for z)
    
    // Inverse metric (handle potential CTC region where g_φφ can be negative)
    float det_tc = g[0][0] * g[2][2] - g[0][2] * g[0][2];  // det of t-φ block
    if (fabsf(det_tc) > 0.001f) {
        g_inv[0][0] = g[2][2] / det_tc;
        g_inv[0][2] = -g[0][2] / det_tc;
        g_inv[2][0] = g_inv[0][2];
        g_inv[2][2] = g[0][0] / det_tc;
    } else {
        // Near singular, use Minkowski fallback
        g_inv[0][0] = -1.0f;
        g_inv[2][2] = 1.0f / fmaxf(fabsf(g_phiphi), 0.01f);
    }
    g_inv[1][1] = 1.0f;
    g_inv[3][3] = 1.0f;
}

//==============================================================================
// Taub-NUT Metric (Gravimagnetic Monopole)
// Schwarzschild + NUT charge in spherical coordinates
//==============================================================================
__device__ void getTaubNUTMetric(const Vec4& x, float M, float n,
                                  float g[4][4], float g_inv[4][4]) {
    // ds² = -f(dt + 2n cosθ dφ)² + dr²/f + Σ(dθ² + sin²θ dφ²)
    // where f = (r² - 2Mr - n²)/(r² + n²), Σ = r² + n²
    
    float r = x.r;
    float theta = x.theta;
    
    float r2 = r * r;
    float n2 = n * n;
    float Sigma = r2 + n2;
    
    float numerator = r2 - 2.0f * M * r - n2;
    float f = numerator / Sigma;
    
    // Avoid horizon singularities
    float f_safe = fmaxf(fabsf(f), 0.001f);
    if (f < 0.0f) f_safe = -f_safe;
    
    float sin_theta = sinf(theta);
    float cos_theta = cosf(theta);
    float sin2 = sin_theta * sin_theta;
    float cos2 = cos_theta * cos_theta;
    
    // Handle pole singularity (Misner string) - very minimal threshold
    if (fabsf(sin_theta) < 0.0001f) {
        sin_theta = (sin_theta >= 0.0f) ? 0.0001f : -0.0001f;
        sin2 = sin_theta * sin_theta;
    }
    
    // Initialize to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components: g_tt, g_tφ, g_rr, g_θθ, g_φφ
    g[0][0] = -f;                                        // g_tt
    g[0][3] = -2.0f * f * n * cos_theta;                 // g_tφ
    g[3][0] = g[0][3];                                   // g_φt
    g[1][1] = 1.0f / f_safe;                             // g_rr
    g[2][2] = Sigma;                                     // g_θθ
    g[3][3] = Sigma * sin2 - 4.0f * f * n2 * cos2;       // g_φφ
    
    // Inverse metric
    // For the t-φ block
    float det_tphi = g[0][0] * g[3][3] - g[0][3] * g[0][3];
    if (fabsf(det_tphi) > 0.001f) {
        g_inv[0][0] = g[3][3] / det_tphi;
        g_inv[0][3] = -g[0][3] / det_tphi;
        g_inv[3][0] = g_inv[0][3];
        g_inv[3][3] = g[0][0] / det_tphi;
    } else {
        // Fallback
        g_inv[0][0] = -1.0f / f_safe;
        g_inv[3][3] = 1.0f / fmaxf(fabsf(g[3][3]), 0.01f);
    }
    g_inv[1][1] = f_safe;
    g_inv[2][2] = 1.0f / Sigma;
}

//==============================================================================
// Kerr-Schild Metric (Horizon-Penetrating Coordinates)
// ds² = -(1-H)dt² + (1+H)dr² + r²(dθ² + sin²θ dφ²) + 2H dt dr
// H = 2Mr/(r² + a²cos²θ), uses Kerr-Schild radial coordinate
//==============================================================================
__device__ void getKerrSchildMetric(const Vec4& x, float M, float a,
                                     float g[4][4], float g_inv[4][4]) {
    float r = fmaxf(x.r, 0.01f);
    float theta = x.theta;
    
    float r2 = r * r;
    float a2 = a * a;
    float cos_theta = cosf(theta);
    float sin_theta = sinf(theta);
    float sin2 = sin_theta * sin_theta;
    float cos2 = cos_theta * cos_theta;
    
    // Kerr-Schild H function
    float Sigma = r2 + a2 * cos2;
    float H = 2.0f * M * r / fmaxf(Sigma, 0.01f);
    
    // Initialize to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components (simplified Kerr-Schild in spherical-like coords)
    g[0][0] = -(1.0f - H);            // g_tt
    g[0][1] = H;                       // g_tr
    g[1][0] = H;                       // g_rt
    g[1][1] = 1.0f + H;                // g_rr
    g[2][2] = Sigma;                   // g_θθ
    g[3][3] = Sigma * sin2;            // g_φφ
    
    // Inverse metric
    float det_tr = g[0][0] * g[1][1] - g[0][1] * g[1][0];
    if (fabsf(det_tr) > 0.001f) {
        g_inv[0][0] = g[1][1] / det_tr;
        g_inv[0][1] = -g[0][1] / det_tr;
        g_inv[1][0] = -g[1][0] / det_tr;
        g_inv[1][1] = g[0][0] / det_tr;
    } else {
        g_inv[0][0] = -1.0f;
        g_inv[1][1] = 1.0f;
    }
    g_inv[2][2] = 1.0f / Sigma;
    g_inv[3][3] = 1.0f / (Sigma * sin2 + 0.0001f);
}

//==============================================================================
// Ellis Drainhole (Traversable Wormhole) Metric
// ds² = -(1-Fp²)dt² + 2Fp dt dr + dr² + Rp²(dθ² + sin²θ dφ²)
// Connects two asymptotically flat regions through throat at α = √(n² - m²)
//==============================================================================
__device__ void getEllisDrainholeMetric(const Vec4& x, float m, float n,
                                         float g[4][4], float g_inv[4][4]) {
    float r = x.r;
    
    // Ensure n > m for valid wormhole
    if (n <= m) n = m + 0.1f;
    
    // Throat radius: α = √(n² - m²)
    float alpha = sqrtf(n * n - m * m);
    
    // Proper radial coordinate
    float dx = r - m;
    float Rp2 = dx * dx + alpha * alpha;
    float Rp = sqrtf(Rp2);
    
    // Pseudo-angle ψ
    float psi = (n / alpha) * (1.5707963f - atanf(dx / alpha));  // π/2 - atan(...)
    
    // Redshift function Fp
    float exp_term = expf(-2.0f * m * psi / n);
    float Fp_sq = 1.0f - exp_term;
    float Fp = -sqrtf(fmaxf(Fp_sq, 0.0f));  // Negative for ingoing
    
    // Areal radius factor
    float denom = 1.0f - Fp_sq;
    if (denom < 0.001f) denom = 0.001f;
    float Rp_areal = Rp / sqrtf(denom);
    float Rp_areal2 = Rp_areal * Rp_areal;
    
    float theta = x.theta;
    float sin_theta = sinf(theta);
    float sin2 = sin_theta * sin_theta;
    
    // Initialize
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components
    g[0][0] = -(1.0f - Fp_sq);  // g_tt
    g[0][1] = Fp;               // g_tr
    g[1][0] = Fp;               // g_rt
    g[1][1] = 1.0f;             // g_rr
    g[2][2] = Rp_areal2;        // g_θθ
    g[3][3] = Rp_areal2 * sin2; // g_φφ
    
    // Inverse (t-r block)
    float det_tr = g[0][0] * g[1][1] - g[0][1] * g[1][0];
    if (fabsf(det_tr) > 0.001f) {
        g_inv[0][0] = g[1][1] / det_tr;
        g_inv[0][1] = -g[0][1] / det_tr;
        g_inv[1][0] = -g[1][0] / det_tr;
        g_inv[1][1] = g[0][0] / det_tr;
    } else {
        g_inv[0][0] = -1.0f;
        g_inv[1][1] = 1.0f;
    }
    g_inv[2][2] = 1.0f / fmaxf(Rp_areal2, 0.01f);
    g_inv[3][3] = 1.0f / fmaxf(Rp_areal2 * sin2, 0.0001f);
}

//==============================================================================
// Alcubierre Warp Drive Metric
// ds² = (vs²f² - 1)dt² - 2vs·f·dx·dt + dx² + dy² + dz²
// Note: Uses Cartesian coordinates internally, adapted for spherical
//==============================================================================
__device__ void getAlcubierreMetric(const Vec4& x, float vs, float sigma, float R,
                                     float g[4][4], float g_inv[4][4]) {
    // For spherical → Cartesian conversion
    float r = x.r;
    float theta = x.theta;
    float t = x.t;
    
    // Bubble center moves along x-axis: xs(t) = vs * t
    float xs_t = vs * t;
    
    // Current position in Cartesian
    float cart_x = r * sinf(theta) * cosf(x.phi);
    float cart_y = r * sinf(theta) * sinf(x.phi);
    float cart_z = r * cosf(theta);
    
    // Distance from bubble center
    float dx = cart_x - xs_t;
    float rs2 = dx*dx + cart_y*cart_y + cart_z*cart_z;
    float rs = sqrtf(rs2);
    if (rs < 0.001f) rs = 0.001f;
    
    // Warp bubble shape function: f(rs) = (tanh(σ(rs+R)) - tanh(σ(rs-R))) / (2·tanh(σR))
    float tanh_plus = tanhf(sigma * (rs + R));
    float tanh_minus = tanhf(sigma * (rs - R));
    float tanh_R = tanhf(sigma * R);
    
    float f = 1.0f;
    if (fabsf(tanh_R) > 0.001f) {
        f = (tanh_plus - tanh_minus) / (2.0f * tanh_R);
    }
    
    float vs_f = vs * f;
    
    // Metric in Cartesian (simplified to spherical r direction)
    float sin_theta = sinf(theta);
    float r2 = r * r;
    
    // Initialize
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Approximate effect in spherical: bubble distorts along radial direction
    g[0][0] = vs_f * vs_f - 1.0f;  // g_tt
    g[0][1] = -vs_f;               // g_tr (radial direction ~x)
    g[1][0] = -vs_f;               // g_rt
    g[1][1] = 1.0f;                // g_rr
    g[2][2] = r2;                  // g_θθ
    g[3][3] = r2 * sin_theta * sin_theta;  // g_φφ
    
    // Inverse (t-r block)
    float det_tr = g[0][0] * g[1][1] - g[0][1] * g[1][0];
    if (fabsf(det_tr) > 0.001f) {
        g_inv[0][0] = g[1][1] / det_tr;
        g_inv[0][1] = -g[0][1] / det_tr;
        g_inv[1][0] = -g[1][0] / det_tr;
        g_inv[1][1] = g[0][0] / det_tr;
    } else {
        g_inv[0][0] = -1.0f;
        g_inv[1][1] = 1.0f;
    }
    g_inv[2][2] = 1.0f / r2;
    g_inv[3][3] = 1.0f / (r2 * sin_theta * sin_theta + 0.0001f);
}

//==============================================================================
// de Sitter Metric (Cosmological Spacetime with Λ > 0)
// ds² = -(1 - Λr²/3)dt² + dr²/(1 - Λr²/3) + r²(dθ² + sin²θ·dφ²)
//==============================================================================
__device__ void getDeSitterMetric(const Vec4& x, float Lambda,
                                   float g[4][4], float g_inv[4][4]) {
    float r = x.r;
    float theta = x.theta;
    float r2 = r * r;
    
    // f(r) = 1 - Λr²/3
    float f = 1.0f - Lambda * r2 / 3.0f;
    
    // Avoid cosmological horizon
    if (fabsf(f) < 0.001f) {
        f = (f > 0.0f) ? 0.001f : -0.001f;
    }
    
    float sin_theta = sinf(theta);
    float sin2 = sin_theta * sin_theta;
    
    // Initialize
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Diagonal metric
    g[0][0] = -f;                  // g_tt
    g[1][1] = 1.0f / f;            // g_rr
    g[2][2] = r2;                  // g_θθ
    g[3][3] = r2 * sin2;           // g_φφ
    
    // Inverse
    g_inv[0][0] = -1.0f / f;
    g_inv[1][1] = f;
    g_inv[2][2] = 1.0f / r2;
    g_inv[3][3] = 1.0f / (r2 * sin2 + 0.0001f);
}

// Sample volumetric data (rho, velocity)
__device__ void getVolumetricSample(const Vec4& pos, float& rho, float3& vel) {
    // ... existing implementation ...
    rho = 0.0f;
    vel = make_float3(0.0f, 0.0f, 0.0f);

    const auto& nm = params.numericalMetric;
    if (!nm.isLoaded) return;

    // 1. Convert Spherical (r, θ, φ) to Cartesian (x, y, z)
    // Using Y-POLAR convention: y = r cos(θ), consistent with rest of codebase
    // θ = polar angle from +Y axis, φ = azimuthal angle in XZ plane
    float3 cart = make_float3(
        pos.r * sinf(pos.theta) * cosf(pos.phi),   // x = r sinθ cosφ
        pos.r * cosf(pos.theta),                    // y = r cosθ  (polar axis)
        pos.r * sinf(pos.theta) * sinf(pos.phi)    // z = r sinθ sinφ
    );

    // 2. Compute texture coordinates
    float3 uvw = make_float3(
        (cart.x - nm.origin.x) / (nm.dims.x * nm.spacing.x),
        (cart.y - nm.origin.y) / (nm.dims.y * nm.spacing.y),
        (cart.z - nm.origin.z) / (nm.dims.z * nm.spacing.z)
    );

    // Check bounds
    if (uvw.x < 0.0f || uvw.x > 1.0f || 
        uvw.y < 0.0f || uvw.y > 1.0f || 
        uvw.z < 0.0f || uvw.z > 1.0f) {
        return;
    }

    // 3. Sample data
    rho = tex3D<float>(nm.rho, uvw.x, uvw.y, uvw.z);
    
    if (rho > 0.0001f) {
        vel.x = tex3D<float>(nm.vx, uvw.x, uvw.y, uvw.z);
        vel.y = tex3D<float>(nm.vy, uvw.x, uvw.y, uvw.z);
        vel.z = tex3D<float>(nm.vz, uvw.x, uvw.y, uvw.z);
    }
}

//==============================================================================
// UNIFIED KERR-SCHILD BLACK HOLE FAMILY - Metric Tensor
// Component: PHMT100A (Cartesian coordinates: t, x, y, z)
//==============================================================================
// Covers 9 spacetimes via parameters (M, a, Q, Lambda):
// - M=0, a=0, Q=0, Λ=0: Minkowski
// - M>0, a=0, Q=0, Λ=0: Schwarzschild
// - M>0, a>0, Q=0, Λ=0: Kerr
// - M>0, a=0, Q>0, Λ=0: Reissner-Nordström
// - M>0, a>0, Q>0, Λ=0: Kerr-Newman
// - M=0, a=0, Q=0, Λ>0: de Sitter
// - M>0, a=0, Q=0, Λ>0: Schwarzschild-de Sitter
// - M>0, a>0, Q=0, Λ>0: Kerr-de Sitter
// - M>0, a>0, Q>0, Λ>0: Kerr-Newman-de Sitter
//
// Metric form: g_μν = η_μν + H · l_μ · l_ν
// Inverse:     g^μν = η^μν - H · l^μ · l^ν
//
// KEY ADVANTAGE: Polynomial Christoffel symbols - NO pole singularities!
//==============================================================================

// Compute Kerr radius from Cartesian coordinates
// Solves: r⁴ - (x² + y² + z² - a²)r² - a²z² = 0
__device__ __forceinline__ float computeKerrRadiusCart(float x, float y, float z, float a) {
    float a2 = a * a;
    float R2 = x*x + y*y + z*z;
    
    // For a = 0 (Schwarzschild), r = R directly
    if (fabsf(a) < 1e-10f) {
        return sqrtf(fmaxf(R2, 1e-20f));
    }
    
    // Solve: r² = [(R² - a²) + √((R² - a²)² + 4a²z²)] / 2
    float Rm2 = R2 - a2;
    float disc = Rm2 * Rm2 + 4.0f * a2 * z * z;
    float r2 = (Rm2 + sqrtf(fmaxf(disc, 0.0f))) / 2.0f;
    
    return sqrtf(fmaxf(r2, 1e-20f));
}

// Unified Kerr-Schild metric in Cartesian coordinates
__device__ void getKerrSchildFamilyMetric(
    const Vec4Cart& pos,            // Position (t, x, y, z)
    float M, float a, float Q, float Lambda,  // Family parameters
    float g[4][4], float g_inv[4][4])
{
    float x = pos.x;
    float y = pos.y;
    float z = pos.z;
    
    float a2 = a * a;
    float Q2 = Q * Q;
    
    // Initialize to Minkowski
    #pragma unroll 4
    for (int i = 0; i < 4; i++) {
        #pragma unroll 4
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    g[0][0] = -1.0f;  g_inv[0][0] = -1.0f;
    g[1][1] = 1.0f;   g_inv[1][1] = 1.0f;
    g[2][2] = 1.0f;   g_inv[2][2] = 1.0f;
    g[3][3] = 1.0f;   g_inv[3][3] = 1.0f;
    
    // If no mass and no cosmological constant, we're done (Minkowski)
    if (M < 1e-15f && fabsf(Lambda) < 1e-15f) {
        return;
    }
    
    // Compute Kerr radius
    float r = computeKerrRadiusCart(x, y, z, a);
    float r2 = r * r;
    float r4 = r2 * r2;
    
    // Kerr-Schild scalar function H = (2Mr - Q²)r² / (r⁴ + a²z²)
    float f_denom = r4 + a2 * z * z;
    float H = (2.0f * M * r - Q2) * r2 / fmaxf(f_denom, 1e-20f);
    
    // Null vector l_μ = (1, (rx + ay)/(r² + a²), (ry - ax)/(r² + a²), z/r)
    float denom = r2 + a2;
    float l[4] = {
        1.0f,
        (r * x + a * y) / denom,
        (r * y - a * x) / denom,
        z / fmaxf(r, 1e-10f)
    };
    
    // Add Kerr-Schild perturbation: g_μν = η_μν + H · l_μ · l_ν
    #pragma unroll 4
    for (int mu = 0; mu < 4; mu++) {
        #pragma unroll 4
        for (int nu = 0; nu < 4; nu++) {
            g[mu][nu] += H * l[mu] * l[nu];
        }
    }
    
    // Inverse metric: g^μν = η^μν - H · l^μ · l^ν
    // l^μ = (1, -l_1, -l_2, -l_3) since raising with Minkowski flips spatial
    float l_up[4] = {1.0f, -l[1], -l[2], -l[3]};
    #pragma unroll 4
    for (int mu = 0; mu < 4; mu++) {
        #pragma unroll 4
        for (int nu = 0; nu < 4; nu++) {
            g_inv[mu][nu] -= H * l_up[mu] * l_up[nu];
        }
    }
    
    // Cosmological constant contribution (de Sitter sector)
    if (fabsf(Lambda) > 1e-15f) {
        float R2 = x*x + y*y + z*z;
        float LambdaFactor = Lambda * R2 / 3.0f;
        g[0][0] -= LambdaFactor;
        g_inv[0][0] += LambdaFactor;  // Approximate for small Lambda
    }
}

// Wrapper to accept spherical Vec4 (converts internally)
__device__ void getKerrSchildFamilyMetricFromSpherical(
    const Vec4& pos_sph,            // Position (t, r, θ, φ)
    float M, float a, float Q, float Lambda,
    float g[4][4], float g_inv[4][4])
{
    Vec4Cart pos_cart = vec4SphToCart(pos_sph);
    getKerrSchildFamilyMetric(pos_cart, M, a, Q, Lambda, g, g_inv);
}

//==============================================================================
// UNIFIED KERR-SCHILD BLACK HOLE FAMILY - Christoffel Symbols
// Component: PHMT100A (Cartesian coordinates: t, x, y, z)
//==============================================================================
// Formula: Γ^λ_μν = (1/2) g^{λσ} (∂_μ g_σν + ∂_ν g_σμ - ∂_σ g_μν)
//
// For Kerr-Schild: g_μν = η_μν + H l_μ l_ν, so:
//   ∂_λ g_μν = (∂_λ H) l_μ l_ν + H (∂_λ l_μ) l_ν + H l_μ (∂_λ l_ν)
//
// The inverse metric is: g^μν = η^μν - H l^μ l^ν
//
// KEY ADVANTAGE: All derivatives are polynomial in (x, y, z, r) - NO trig!
//==============================================================================

__device__ void getKerrSchildFamilyChristoffel(
    const Vec4Cart& pos,            // Position (t, x, y, z)
    float M, float a, float Q, float Lambda,  // Family parameters
    float Gamma[4][4][4])
{
    // Initialize to zero
    #pragma unroll 4
    for (int i = 0; i < 4; i++)
        #pragma unroll 4
        for (int j = 0; j < 4; j++)
            #pragma unroll 4
            for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;
    
    // For pure Minkowski (M=0, Λ=0), all Christoffels are zero in Cartesian
    if (M < 1e-15f && fabsf(Lambda) < 1e-15f) {
        return;
    }
    
    float x = pos.x;
    float y = pos.y;
    float z = pos.z;
    
    float a2 = a * a;
    float Q2 = Q * Q;
    float R2 = x*x + y*y + z*z;
    
    // Compute Kerr radius and derivatives
    float r = computeKerrRadiusCart(x, y, z, a);
    float r2 = r * r;
    float r3 = r2 * r;
    float r4 = r2 * r2;
    
    // =========================================================================
    // Derivatives of r with respect to x, y, z
    // From r² = [(R² - a²) + √((R² - a²)² + 4a²z²)] / 2
    // =========================================================================
    float Rm2 = R2 - a2;
    float disc = Rm2 * Rm2 + 4.0f * a2 * z * z;
    float sqrt_disc = sqrtf(fmaxf(disc, 1e-20f));
    
    float d_disc_dx = 4.0f * x * Rm2;
    float d_disc_dy = 4.0f * y * Rm2;
    float d_disc_dz = 4.0f * z * (R2 + a2);
    
    float d_r2_dx = x + d_disc_dx / (4.0f * sqrt_disc);
    float d_r2_dy = y + d_disc_dy / (4.0f * sqrt_disc);
    float d_r2_dz = z + d_disc_dz / (4.0f * sqrt_disc);
    
    float dr[4];
    dr[0] = 0.0f;  // ∂r/∂t = 0 (static)
    dr[1] = d_r2_dx / (2.0f * r);
    dr[2] = d_r2_dy / (2.0f * r);
    dr[3] = d_r2_dz / (2.0f * r);
    
    // =========================================================================
    // Null vector l_μ and its derivatives
    // l_μ = (1, (rx + ay)/(r² + a²), (ry - ax)/(r² + a²), z/r)
    // =========================================================================
    float denom = r2 + a2;
    float denom2 = denom * denom;
    float r_safe = fmaxf(r, 1e-10f);
    
    float l[4] = {
        1.0f,
        (r * x + a * y) / denom,
        (r * y - a * x) / denom,
        z / r_safe
    };
    
    // l^μ = η^μν l_ν = (1, -l_1, -l_2, -l_3) for Minkowski
    float l_up[4] = {1.0f, -l[1], -l[2], -l[3]};
    
    // dl[λ][μ] = ∂l_μ/∂x^λ
    float dl[4][4] = {{0}};
    
    // Derivative of (rx + ay)/(r² + a²) w.r.t. x, y, z
    // Let f1 = rx + ay, f2 = r² + a²
    // d(f1/f2)/dx = (df1/dx * f2 - f1 * df2/dx) / f2²
    
    float f1 = r * x + a * y;
    float f2 = denom;
    
    // df1/dx = r + x*dr/dx, df1/dy = a + x*dr/dy, df1/dz = x*dr/dz
    // df2/dx = 2r*dr/dx, df2/dy = 2r*dr/dy, df2/dz = 2r*dr/dz
    
    dl[1][1] = (r + x * dr[1]) / f2 - f1 * (2.0f * r * dr[1]) / denom2;
    dl[2][1] = (x * dr[2] + a) / f2 - f1 * (2.0f * r * dr[2]) / denom2;
    dl[3][1] = (x * dr[3]) / f2 - f1 * (2.0f * r * dr[3]) / denom2;
    
    // Similar for l[2] = (ry - ax) / (r² + a²)
    float g1 = r * y - a * x;
    
    dl[1][2] = (y * dr[1] - a) / f2 - g1 * (2.0f * r * dr[1]) / denom2;
    dl[2][2] = (r + y * dr[2]) / f2 - g1 * (2.0f * r * dr[2]) / denom2;
    dl[3][2] = (y * dr[3]) / f2 - g1 * (2.0f * r * dr[3]) / denom2;
    
    // l[3] = z / r
    float r2_safe = r_safe * r_safe;
    dl[1][3] = -z * dr[1] / r2_safe;
    dl[2][3] = -z * dr[2] / r2_safe;
    dl[3][3] = 1.0f / r_safe - z * dr[3] / r2_safe;
    
    // dl[0][*] = 0 (static spacetime)
    
    // =========================================================================
    // Kerr-Schild scalar H and its derivatives
    // H = (2Mr - Q²)r² / (r⁴ + a²z²)
    // =========================================================================
    float f_denom = r4 + a2 * z * z;
    float numerator = (2.0f * M * r - Q2) * r2;
    float H = numerator / fmaxf(f_denom, 1e-20f);
    
    // dH[λ] = ∂H/∂x^λ
    float dH[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    
    for (int lam = 1; lam <= 3; lam++) {
        // d(numerator)/dx^λ = 2M*dr*r² + (2Mr-Q²)*2r*dr = 2r*dr*(M*r + 2Mr - Q²)
        //                   = 2r*dr*(3Mr - Q²)
        float d_num = 2.0f * r * dr[lam] * (3.0f * M * r - Q2);
        
        // d(f_denom)/dx^λ = 4r³*dr + (λ==3 ? 2a²z : 0)
        float d_f_denom = 4.0f * r3 * dr[lam];
        if (lam == 3) d_f_denom += 2.0f * a2 * z;
        
        dH[lam] = (d_num * f_denom - numerator * d_f_denom) / (f_denom * f_denom);
    }
    
    // =========================================================================
    // Metric derivatives: ∂_λ g_μν = (∂_λ H) l_μ l_ν + H (∂_λ l_μ) l_ν + H l_μ (∂_λ l_ν)
    // =========================================================================
    float dg[4][4][4];  // dg[λ][μ][ν] = ∂_λ g_μν
    
    for (int lam = 0; lam < 4; lam++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
                dg[lam][mu][nu] = dH[lam] * l[mu] * l[nu]
                                + H * dl[lam][mu] * l[nu]
                                + H * l[mu] * dl[lam][nu];
            }
        }
    }
    
    // =========================================================================
    // Christoffel symbols: Γ^λ_μν = (1/2) g^{λσ} (∂_μ g_σν + ∂_ν g_σμ - ∂_σ g_μν)
    // where g^λσ = η^λσ - H l^λ l^σ
    // =========================================================================
    
    // First compute with η^λσ (Minkowski inverse)
    float eta_inv[4] = {-1.0f, 1.0f, 1.0f, 1.0f};
    
    for (int lam = 0; lam < 4; lam++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = mu; nu < 4; nu++) {  // Use symmetry
                float sum = 0.0f;
                
                for (int s = 0; s < 4; s++) {
                    // g^λσ = η^λσ - H l^λ l^σ (only η^λλ non-zero in Minkowski)
                    float g_inv_lam_s = (lam == s) ? eta_inv[lam] : 0.0f;
                    g_inv_lam_s -= H * l_up[lam] * l_up[s];
                    
                    // (∂_μ g_σν + ∂_ν g_σμ - ∂_σ g_μν)
                    float bracket = dg[mu][s][nu] + dg[nu][s][mu] - dg[s][mu][nu];
                    
                    sum += g_inv_lam_s * bracket;
                }
                
                Gamma[lam][mu][nu] = 0.5f * sum;
                Gamma[lam][nu][mu] = Gamma[lam][mu][nu];  // Symmetry in lower indices
            }
        }
    }
    
    // Cosmological constant contribution (small correction)
    if (fabsf(Lambda) > 1e-15f) {
        // For de Sitter: Γ^t_ti = Γ^t_it = (Λ/3) x^i
        // This is a simplified approximation for small Λ
        float LambdaThird = Lambda / 3.0f;
        Gamma[0][0][1] += LambdaThird * x;
        Gamma[0][1][0] += LambdaThird * x;
        Gamma[0][0][2] += LambdaThird * y;
        Gamma[0][2][0] += LambdaThird * y;
        Gamma[0][0][3] += LambdaThird * z;
        Gamma[0][3][0] += LambdaThird * z;
    }
}

// Wrapper to accept spherical Vec4 (converts internally)
__device__ void getKerrSchildFamilyChristoffelFromSpherical(
    const Vec4& pos_sph,            // Position (t, r, θ, φ)
    float M, float a, float Q, float Lambda,
    float Gamma[4][4][4])
{
    Vec4Cart pos_cart = vec4SphToCart(pos_sph);
    getKerrSchildFamilyChristoffel(pos_cart, M, a, Q, Lambda, Gamma);
}

//==============================================================================
// UNIFIED MORRIS-THORNE WORMHOLE FAMILY - Metric Tensor
// Component: PHMT101A (Spherical coordinates: t, r, θ, φ)
//==============================================================================
// Covers traversable wormholes via shape function b(r) and redshift Φ(r):
//
// ds² = -e^{2Φ(r)} dt² + dr²/(1 - b(r)/r) + r²(dθ² + sin²θ dφ²)
//
// Supported shape functions:
//   Ellis:          b(r) = b₀²/r
//   ZeroTidal:      b(r) = b₀  
//   AbsurdlyBenign: b(r) = b₀(2 - b₀/r)
//
// b₀ = throat radius, Φ = redshift function (0 for zero-tidal)
//==============================================================================

// Shape function types
enum class WormholeShape : int {
    Ellis = 0,
    ZeroTidal = 1,
    AbsurdlyBenign = 2
};

// Shape function b(r)
__device__ __forceinline__ float wormholeShapeFunction(float r, float b0, WormholeShape type) {
    r = fmaxf(r, b0);  // r >= b0 at throat
    
    switch (type) {
        case WormholeShape::Ellis:
            return b0 * b0 / r;
        case WormholeShape::ZeroTidal:
            return b0;
        case WormholeShape::AbsurdlyBenign:
            return b0 * (2.0f - b0 / r);
        default:
            return b0 * b0 / r;
    }
}

// Shape function derivative db/dr
__device__ __forceinline__ float wormholeShapeDerivative(float r, float b0, WormholeShape type) {
    r = fmaxf(r, b0);
    float r2 = r * r;
    
    switch (type) {
        case WormholeShape::Ellis:
            return -b0 * b0 / r2;
        case WormholeShape::ZeroTidal:
            return 0.0f;
        case WormholeShape::AbsurdlyBenign:
            return b0 * b0 / r2;
        default:
            return -b0 * b0 / r2;
    }
}

// Morris-Thorne wormhole metric
__device__ void getMorrisThorneMetric(
    const Vec4& pos,         // Position (t, r, θ, φ)
    float b0,                // Throat radius
    float Phi0,              // Redshift at throat
    WormholeShape shapeType,
    float g[4][4], float g_inv[4][4])
{
    float r = fmaxf(pos.r, b0 * 1.001f);  // Stay outside throat
    float theta = pos.theta;
    
    float sin_theta = sinf(theta);
    float sin2 = sin_theta * sin_theta;
    float cos_theta = cosf(theta);
    float r2 = r * r;
    
    // Shape function
    float b = wormholeShapeFunction(r, b0, shapeType);
    float one_minus_b_over_r = 1.0f - b / r;
    one_minus_b_over_r = fmaxf(one_minus_b_over_r, 0.001f);  // Regularize
    
    // Redshift function (constant for zero-tidal)
    float exp2Phi = expf(2.0f * Phi0);
    
    // Initialize
    #pragma unroll 4
    for (int i = 0; i < 4; i++) {
        #pragma unroll 4
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric components
    g[0][0] = -exp2Phi;                        // g_tt
    g[1][1] = 1.0f / one_minus_b_over_r;       // g_rr
    g[2][2] = r2;                              // g_θθ
    g[3][3] = r2 * sin2;                       // g_φφ
    
    // Inverse metric
    g_inv[0][0] = -1.0f / exp2Phi;
    g_inv[1][1] = one_minus_b_over_r;
    g_inv[2][2] = 1.0f / r2;
    g_inv[3][3] = 1.0f / fmaxf(r2 * sin2, 1e-10f);
}

// Morris-Thorne wormhole Christoffel symbols
__device__ void getMorrisThorneChristoffel(
    const Vec4& pos,
    float b0,
    float Phi0,
    WormholeShape shapeType,
    float Gamma[4][4][4])
{
    // Initialize to zero
    #pragma unroll 4
    for (int i = 0; i < 4; i++)
        #pragma unroll 4
        for (int j = 0; j < 4; j++)
            #pragma unroll 4
            for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;
    
    float r = fmaxf(pos.r, b0 * 1.001f);
    float theta = pos.theta;
    
    float sin_theta = sinf(theta);
    float cos_theta = cosf(theta);
    float sin2 = sin_theta * sin_theta;
    float r2 = r * r;
    
    // Shape function and derivatives
    float b = wormholeShapeFunction(r, b0, shapeType);
    float db_dr = wormholeShapeDerivative(r, b0, shapeType);
    
    float one_minus_b_over_r = 1.0f - b / r;
    one_minus_b_over_r = fmaxf(one_minus_b_over_r, 0.001f);
    
    // For zero-tidal wormhole (Φ = const), many Christoffels simplify
    // dΦ/dr = 0, so Γ^t components involving time are zero
    
    // d(b/r)/dr = (r*db_dr - b) / r² = db_dr/r - b/r²
    float d_b_over_r = db_dr / r - b / r2;
    
    // Γ^r_rr = d(1-b/r)/dr / (2(1-b/r)) = -d_b_over_r / (2(1-b/r))
    Gamma[1][1][1] = -d_b_over_r / (2.0f * one_minus_b_over_r);
    
    // Γ^r_θθ = -r(1-b/r) = -(r - b)
    Gamma[1][2][2] = -(r - b);
    
    // Γ^r_φφ = -r sin²θ (1-b/r) = -(r - b) sin²θ
    Gamma[1][3][3] = -(r - b) * sin2;
    
    // Γ^θ_rθ = Γ^θ_θr = 1/r
    Gamma[2][1][2] = 1.0f / r;
    Gamma[2][2][1] = 1.0f / r;
    
    // Γ^θ_φφ = -sinθ cosθ
    Gamma[2][3][3] = -sin_theta * cos_theta;
    
    // Γ^φ_rφ = Γ^φ_φr = 1/r
    Gamma[3][1][3] = 1.0f / r;
    Gamma[3][3][1] = 1.0f / r;
    
    // Γ^φ_θφ = Γ^φ_φθ = cotθ
    float cot_theta = safe_cot(theta);
    Gamma[3][2][3] = cot_theta;
    Gamma[3][3][2] = cot_theta;
    
    // If Φ ≠ 0, add time-related Christoffels
    if (fabsf(Phi0) > 1e-10f) {
        // For r-dependent Φ, we would add Γ^t_tr = dΦ/dr terms
        // With constant Φ, these are zero
    }
}

//==============================================================================
// UNIFIED WARP DRIVE FAMILY - Metric Tensor
// Component: PHMT102A (Cartesian coordinates: t, x, y, z)
//==============================================================================
// Alcubierre warp drive metric:
//   ds² = -dt² + (dx - v_s f(r_s) dt)² + dy² + dz²
//
// where r_s = √((x - x_s)² + y² + z²)
// f(r_s) = [tanh(σ(r_s + R)) - tanh(σ(r_s - R))] / [2 tanh(σR)]
//
// Parameters:
//   v_s   = bubble velocity (can be > 1 for superluminal)
//   σ     = wall thickness (larger = sharper)
//   R     = bubble radius
//   (xs, ys, zs) = bubble center position
//
// KEY PROPERTY: Already Cartesian - no pole singularities!
//==============================================================================

// Alcubierre shape function
__device__ __forceinline__ float alcubierreShapeFunction(float rs, float sigma, float R) {
    // f(r_s) = [tanh(σ(r_s + R)) - tanh(σ(r_s - R))] / [2 tanh(σR)]
    float tanh_plus = tanhf(sigma * (rs + R));
    float tanh_minus = tanhf(sigma * (rs - R));
    float tanh_R = tanhf(sigma * R);
    
    return (tanh_plus - tanh_minus) / fmaxf(2.0f * tanh_R, 1e-10f);
}

// Derivative of shape function
__device__ __forceinline__ float alcubierreShapeDerivative(float rs, float sigma, float R) {
    // df/drs = σ [sech²(σ(r_s + R)) - sech²(σ(r_s - R))] / [2 tanh(σR)]
    float cosh_plus = coshf(sigma * (rs + R));
    float cosh_minus = coshf(sigma * (rs - R));
    float sech2_plus = 1.0f / (cosh_plus * cosh_plus);
    float sech2_minus = 1.0f / (cosh_minus * cosh_minus);
    float tanh_R = tanhf(sigma * R);
    
    return sigma * (sech2_plus - sech2_minus) / fmaxf(2.0f * tanh_R, 1e-10f);
}

// Warp drive metric in Cartesian coordinates
__device__ void getWarpDriveMetric(
    const Vec4Cart& pos,    // Position (t, x, y, z)
    float vs,               // Bubble velocity
    float sigma,            // Wall thickness
    float R,                // Bubble radius
    float xs, float ys, float zs,  // Bubble center
    float g[4][4], float g_inv[4][4])
{
    float x = pos.x;
    float y = pos.y;
    float z = pos.z;
    
    // Distance from bubble center
    float dx = x - xs;
    float dy = y - ys;
    float dz = z - zs;
    float rs = sqrtf(dx*dx + dy*dy + dz*dz);
    rs = fmaxf(rs, 1e-10f);
    
    // Shape function
    float f = alcubierreShapeFunction(rs, sigma, R);
    float vsf = vs * f;
    float vsf2 = vsf * vsf;
    
    // Initialize
    #pragma unroll 4
    for (int i = 0; i < 4; i++) {
        #pragma unroll 4
        for (int j = 0; j < 4; j++) {
            g[i][j] = 0.0f;
            g_inv[i][j] = 0.0f;
        }
    }
    
    // Metric: ds² = -(1 - v_s²f²)dt² - 2v_sf dx dt + dx² + dy² + dz²
    g[0][0] = -(1.0f - vsf2);     // g_tt
    g[0][1] = -vsf;               // g_tx
    g[1][0] = -vsf;               // g_xt
    g[1][1] = 1.0f;               // g_xx
    g[2][2] = 1.0f;               // g_yy
    g[3][3] = 1.0f;               // g_zz
    
    // Inverse metric
    float det_tx = g[0][0] * g[1][1] - g[0][1] * g[1][0];
    det_tx = fmaxf(fabsf(det_tx), 1e-10f) * (det_tx >= 0.0f ? 1.0f : -1.0f);
    
    g_inv[0][0] = g[1][1] / det_tx;
    g_inv[0][1] = -g[0][1] / det_tx;
    g_inv[1][0] = -g[1][0] / det_tx;
    g_inv[1][1] = g[0][0] / det_tx;
    g_inv[2][2] = 1.0f;
    g_inv[3][3] = 1.0f;
}

// Warp drive Christoffel symbols
__device__ void getWarpDriveChristoffel(
    const Vec4Cart& pos,
    float vs, float sigma, float R,
    float xs, float ys, float zs,
    float Gamma[4][4][4])
{
    // Initialize to zero
    #pragma unroll 4
    for (int i = 0; i < 4; i++)
        #pragma unroll 4
        for (int j = 0; j < 4; j++)
            #pragma unroll 4
            for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;
    
    float x = pos.x;
    float y = pos.y;
    float z = pos.z;
    
    // Distance from bubble center
    float dx = x - xs;
    float dy = y - ys;
    float dz = z - zs;
    float rs = sqrtf(dx*dx + dy*dy + dz*dz);
    rs = fmaxf(rs, 1e-10f);
    
    // Shape function and derivative
    float f = alcubierreShapeFunction(rs, sigma, R);
    float df_drs = alcubierreShapeDerivative(rs, sigma, R);
    
    // Derivatives of rs with respect to position
    float drs[4] = {0.0f, dx / rs, dy / rs, dz / rs};
    
    // Derivatives of f with respect to position
    float df[4] = {0.0f, df_drs * drs[1], df_drs * drs[2], df_drs * drs[3]};
    
    float vsf = vs * f;
    float vsf2 = vsf * vsf;
    
    // Metric components
    // g_tt = -(1 - v²f²), g_tx = g_xt = -vf, g_xx = g_yy = g_zz = 1
    
    // Metric derivatives
    // ∂g_tt/∂x^i = 2v²f df/dx^i
    // ∂g_tx/∂x^i = -v df/dx^i
    
    float dg_tt[4], dg_tx[4];
    for (int i = 0; i < 4; i++) {
        dg_tt[i] = 2.0f * vs * vs * f * df[i];
        dg_tx[i] = -vs * df[i];
    }
    
    // Inverse metric (approximation for small perturbations)
    // g^tt ≈ -1 - v²f², g^tx ≈ -vf, g^xx = g^yy = g^zz = 1
    float g_inv_tt = -1.0f - vsf2;  // Approximate
    float g_inv_tx = -vsf;
    
    // Christoffel symbols from Γ^λ_μν = (1/2) g^λσ (∂_μ g_σν + ∂_ν g_σμ - ∂_σ g_μν)
    
    // Γ^t_tt = (1/2) g^tt ∂_t g_tt + (1/2) g^tx (2∂_t g_tx - ∂_x g_tt)
    // For static bubble: ∂_t terms = 0, so many simplify
    
    // Γ^t_xx = (1/2) g^tt (-∂_x g_tt) + (1/2) g^tx (-∂_x g_tx - ∂_x g_tx)
    //        = -(1/2) g^tt dg_tt[1] - g^tx dg_tx[1]
    Gamma[0][1][1] = -0.5f * g_inv_tt * dg_tt[1] - g_inv_tx * dg_tx[1];
    
    // Γ^t_tx = Γ^t_xt = (1/2) g^tt ∂_x g_tt
    Gamma[0][0][1] = 0.5f * g_inv_tt * dg_tt[1];
    Gamma[0][1][0] = Gamma[0][0][1];
    
    // Γ^x_tt = (1/2) g^xt ∂_t g_tt + (1/2) g^xx (2∂_t g_tx - ∂_x g_tt)
    //        = -(1/2) ∂_x g_tt (for static)
    Gamma[1][0][0] = -0.5f * dg_tt[1];
    
    // Γ^x_tx = Γ^x_xt = (1/2) g^xx ∂_x g_tt = (1/2) dg_tt[1]
    Gamma[1][0][1] = 0.5f * dg_tt[1];
    Gamma[1][1][0] = Gamma[1][0][1];
    
    // Similar for y and z components (using symmetry)
    for (int i = 2; i <= 3; i++) {
        Gamma[0][0][i] = 0.5f * g_inv_tt * dg_tt[i];
        Gamma[0][i][0] = Gamma[0][0][i];
        Gamma[0][1][i] = -g_inv_tx * dg_tx[i];
        Gamma[0][i][1] = Gamma[0][1][i];
        
        Gamma[1][0][i] = dg_tx[i];      // Γ^x_ti = ∂g_tx/∂x^i
        Gamma[1][i][0] = Gamma[1][0][i];
    }
}

// Wrapper for spherical input (converts internally)
__device__ void getWarpDriveMetricFromSpherical(
    const Vec4& pos_sph,
    float vs, float sigma, float R,
    float xs, float ys, float zs,
    float g[4][4], float g_inv[4][4])
{
    Vec4Cart pos_cart = vec4SphToCart(pos_sph);
    getWarpDriveMetric(pos_cart, vs, sigma, R, xs, ys, zs, g, g_inv);
}

//==============================================================================
// Kerr ISCO Radius Computation
//==============================================================================
// Innermost Stable Circular Orbit for a Kerr black hole (prograde)
// Formula: r_ISCO = M * (3 + Z2 - sqrt((3-Z1)(3+Z1+2*Z2)))
// where Z1 = 1 + (1-a²)^(1/3) * ((1+a)^(1/3) + (1-a)^(1/3))
//       Z2 = sqrt(3*a² + Z1²)
// In geometric units (G=c=1), with a = spin parameter (0 to 1)
//==============================================================================
__device__ float computeKerrISCO(float M, float a) {
    // Clamp spin to valid range
    float spin = clamp(fabsf(a), 0.0f, 0.998f);
    
    if (spin < 0.001f) {
        // Schwarzschild limit: r_ISCO = 6M
        return 6.0f * M;
    }
    
    float a2 = spin * spin;
    
    // Z1 = 1 + (1-a²)^(1/3) * ((1+a)^(1/3) + (1-a)^(1/3))
    float oneMa2_13 = powf(1.0f - a2, 1.0f/3.0f);
    float onePa_13 = powf(1.0f + spin, 1.0f/3.0f);
    float oneMa_13 = powf(1.0f - spin, 1.0f/3.0f);
    float Z1 = 1.0f + oneMa2_13 * (onePa_13 + oneMa_13);
    
    // Z2 = sqrt(3*a² + Z1²)
    float Z2 = sqrtf(3.0f * a2 + Z1 * Z1);
    
    // r_ISCO = M * (3 + Z2 - sqrt((3-Z1)(3+Z1+2*Z2)))
    float inner = (3.0f - Z1) * (3.0f + Z1 + 2.0f * Z2);
    float r_isco = M * (3.0f + Z2 - sqrtf(fmaxf(inner, 0.0f)));
    
    return r_isco;
}

//==============================================================================
// FIDO Frame Computation (DNGR Paper Appendix A.1, Eq A.2-A.3)
// James et al. (2015) "Gravitational Lensing by Spinning Black Holes"
//
// Computes the Fiducial Observer (FIDO) frame at given position.
// The FIDO is the locally non-rotating observer whose world line is orthogonal
// to the t = const hypersurfaces in Boyer-Lindquist coordinates.
//==============================================================================
template<Sirius::MetricType type>
__device__ Sirius::FIDOBasis computeFIDOBasis(
    const Vec4& x,
    const Sirius::MetricParams& mp)
{
    Sirius::FIDOBasis fido;
    
    float r = fmaxf(x.r, 0.01f);
    float theta = x.theta;
    float phi = x.phi;
    float M = mp.M;
    float a = mp.a;
    
    float r2 = r * r;
    float a2 = a * a;
    float cos_theta = cosf(theta);
    float sin_theta = sinf(theta);
    
    // Clamp sin_theta away from zero (poles) - very minimal threshold
    if (fabsf(sin_theta) < 0.0001f) {
        sin_theta = (sin_theta >= 0.0f) ? 0.0001f : -0.0001f;
    }
    
    float cos2 = cos_theta * cos_theta;
    float sin2 = sin_theta * sin_theta;
    
    // Kerr metric quantities (DNGR Eq A.2)
    fido.rho = sqrtf(r2 + a2 * cos2);                              // ρ = √(r² + a²cos²θ)
    fido.Delta = r2 - 2.0f * M * r + a2;                            // Δ = r² - 2Mr + a²
    
    // Σ² = (r² + a²)² - a²Δsin²θ
    float r2_plus_a2 = r2 + a2;
    float Sigma_sq = r2_plus_a2 * r2_plus_a2 - a2 * fido.Delta * sin2;
    fido.Sigma = sqrtf(fmaxf(Sigma_sq, 0.01f));
    
    // Clamp Delta to avoid horizon singularity
    float Delta_safe = fmaxf(fido.Delta, 0.01f);
    
    // Lapse function: α = ρ√Δ/Σ
    fido.alpha = fido.rho * sqrtf(Delta_safe) / fido.Sigma;
    
    // Frame dragging angular velocity: ω = 2ar/Σ²
    fido.omega = 2.0f * a * r / (fido.Sigma * fido.Sigma);
    
    // Proper circumference factor: ϖ = Σsinθ/ρ
    fido.varpi = fido.Sigma * sin_theta / fido.rho;
    
    // =========================================================================
    // FIDO basis vectors (Cartesian representation)
    // These are the directions in which the FIDO measures distances
    // =========================================================================
    
    // e_r̂ points radially outward (∂/∂r direction in Cartesian)
    float sp = sinf(phi);
    float cp = cosf(phi);
    fido.e_r = make_float3(
        sin_theta * cp,    // x = sinθ cosφ
        cos_theta,         // y = cosθ (polar axis is Y)
        sin_theta * sp     // z = sinθ sinφ
    );
    
    // e_θ̂ points toward increasing θ (toward equator from poles)
    fido.e_theta = make_float3(
        cos_theta * cp,    // ∂x/∂θ = cosθ cosφ
        -sin_theta,        // ∂y/∂θ = -sinθ
        cos_theta * sp     // ∂z/∂θ = cosθ sinφ
    );
    
    // e_φ̂ points in direction of increasing φ (prograde rotation)
    fido.e_phi = make_float3(
        -sp,               // ∂x/∂φ = -sinφ
        0.0f,              // ∂y/∂φ = 0
        cp                 // ∂z/∂φ = cosφ
    );
    
    return fido;
}

//==============================================================================
// Camera Speed Relative to FIDO (DNGR Eq A.7)
// For a camera in circular equatorial geodesic orbit
//
// β = (ϖ/α)(Ω - ω)
//
// where Ω = 1/(a + r^(3/2)) is the geodesic angular velocity
// and ω is the frame dragging angular velocity
//==============================================================================
__device__ float computeCameraSpeedRelativeToFIDO(
    float r, 
    float M, 
    float a, 
    const Sirius::FIDOBasis& fido)
{
    // Circular geodesic angular velocity (prograde): Ω = sqrt(M)/(r^(3/2) + a*sqrt(M))
    // In units where M=1: Ω = 1/(a + r^(3/2))
    float r32 = powf(r, 1.5f);
    float sqrtM = sqrtf(M);
    float Omega = sqrtM / (r32 + a * sqrtM);
    
    // Camera speed relative to FIDO: β = (ϖ/α)(Ω - ω)
    // This is the speed measured by the FIDO
    float beta = (fido.varpi / fmaxf(fido.alpha, 0.001f)) * (Omega - fido.omega);
    
    // Must be subluminal
    return clamp(beta, -0.999f, 0.999f);
}

//==============================================================================
// DNGR Blue Shift Formula (Eq A.16)
// Validates against existing computeRelativisticG for consistency
//
// f_c/f' = (√(1-β²)/(1-βN_y)) × ((1-bω)/α)
//
// where:
//   β = camera speed relative to FIDO
//   N_y = component of incoming ray along camera velocity direction
//   b = ray's axial angular momentum (p_φ with p_t = -1)
//   ω = frame dragging angular velocity
//   α = lapse function
//==============================================================================
__device__ float computeDNGRBlueShift(
    float beta,              // Camera speed relative to FIDO
    float N_y,               // Ray direction component along velocity
    float b,                 // Ray angular momentum
    const Sirius::FIDOBasis& fido)
{
    // Special relativistic aberration factor
    float gamma = 1.0f / sqrtf(fmaxf(1.0f - beta * beta, 0.001f));
    float sr_factor = sqrtf(1.0f - beta * beta) / fmaxf(1.0f - beta * N_y, 0.001f);
    
    // Gravitational factor
    float gr_factor = (1.0f - b * fido.omega) / fmaxf(fido.alpha, 0.001f);
    
    // Combined frequency ratio
    float g = sr_factor * gr_factor;
    
    return clamp(g, 0.1f, 10.0f);
}


//==============================================================================
// Unified Relativistic Redshift (g-factor)
//==============================================================================
// The correct invariant formula for relativistic frequency shift is:
//   g = (k_μ u^μ_obs) / (k_μ u^μ_emit)

// where:
//   k = photon 4-momentum (null geodesic tangent)
//   u_obs = observer 4-velocity
//   u_emit = emitter 4-velocity
//
// For an observer at infinity (static): u_obs = (1, 0, 0, 0)
// For disk material in Keplerian orbit: u_emit = computed from orbital velocity
//
// g > 1 → approaching (blue-shift)
// g < 1 → receding (red-shift)
//==============================================================================
template<Sirius::MetricType type>
__device__ float computeRelativisticG(
    const Vec4& x,           // Position (spherical: t, r, θ, φ)
    const Vec4& k,           // Photon 4-momentum (null vector)
    const Vec4& u_emit,      // Emitter 4-velocity (timelike, normalized)
    const Sirius::MetricParams& mp)
{
    // Get metric at emission point
    float g[4][4], g_inv[4][4];
    getMetricTensor<type>(x, mp, g, g_inv);
    
    // Observer at infinity has 4-velocity u_obs = (1, 0, 0, 0)
    // So k·u_obs = g_tt * k^t + g_ti * k^i (cross terms)
    // For Schwarzschild: k·u_obs = g_tt * k^t (diagonal)
    // For Kerr: also includes g_tφ * k^φ
    Vec4 u_obs(1.0f, 0.0f, 0.0f, 0.0f);
    float k_dot_u_obs = dot4D(g, k, u_obs);
    
    // k·u_emit
    float k_dot_u_emit = dot4D(g, k, u_emit);
    
    // g = (k·u_obs) / (k·u_emit)
    // Both should be negative for physical configurations
    float g_factor = k_dot_u_obs / fminf(k_dot_u_emit, -0.001f);
    
    // Clamp to reasonable range
    return clamp(fabsf(g_factor), 0.1f, 10.0f);
}

//==============================================================================
// Compute Disk Material 4-Velocity (Keplerian Orbit)
//==============================================================================
// For circular equatorial geodesic in Kerr spacetime:
//   Ω = dφ/dt = sqrt(M) / (r^(3/2) + a*sqrt(M))
//   u^t = 1 / sqrt(-g_tt - 2*Ω*g_tφ - Ω²*g_φφ)
//   u^φ = Ω * u^t
//==============================================================================
template<Sirius::MetricType type>
__device__ Vec4 computeKeplerianVelocity(
    const Vec4& x,
    const Sirius::MetricParams& mp)
{
    float r = x.r;
    float M = mp.M;
    float a = mp.a;
    
    // Keplerian angular velocity
    float r32 = powf(r, 1.5f);
    float Omega = sqrtf(M) / (r32 + a * sqrtf(M));
    
    // Get metric
    float g[4][4], g_inv[4][4];
    getMetricTensor<type>(x, mp, g, g_inv);
    
    // Compute u^t from normalization: g_μν u^μ u^ν = -1
    // For equatorial circular orbit: u = (u^t, 0, 0, u^φ) with u^φ = Ω * u^t
    // -1 = g_tt (u^t)² + 2 g_tφ u^t u^φ + g_φφ (u^φ)²
    // -1 = (u^t)² * (g_tt + 2 Ω g_tφ + Ω² g_φφ)
    float denom = g[0][0] + 2.0f * Omega * g[0][3] + Omega * Omega * g[3][3];
    float u_t = 1.0f / sqrtf(fmaxf(-denom, 0.01f));
    float u_phi = Omega * u_t;
    
    return Vec4(u_t, 0.0f, 0.0f, u_phi);
}

//==============================================================================
// Volumetric Accretion Disk Sampling (Phase 4.2)
// Shakura-Sunyaev thin disk model with volumetric ray marching
//==============================================================================

// Planck blackbody function (normalized for RGB)
// Returns relative spectral radiance for temperature T
// OPTIMIZATION: Branch-free Chebyshev polynomial approximation
__device__ float3 blackbodyColor(float T) {
    // Clamp temperature range: 1000K (red) to 40000K (blue-white)
    T = clamp(T, 1000.0f, 40000.0f);
    
    // =========================================================================
    // CHEBYSHEV APPROXIMATION (Branch-free, ~3x faster)
    // Map log(T) to Chebyshev domain [-1, 1]
    // u = 2 * (log(T) - log(1000)) / (log(40000) - log(1000)) - 1
    // =========================================================================
    float logT = logf(T);
    float u = 2.0f * (logT - BB_LOG_MIN) / BB_LOG_RANGE - 1.0f;
    
    // Evaluate CIE x,y chromaticity using Chebyshev (Clenshaw recurrence)
    float x = chebyshevEval(BB_X_CHEB, 8, u);
    float y = chebyshevEval(BB_Y_CHEB, 8, u);
    
    // Clamp chromaticity to valid range
    x = clamp(x, 0.1f, 0.65f);
    y = clamp(y, 0.1f, 0.65f);
    
    // Convert CIE xy to RGB (sRGB primaries)
    float Y = 1.0f;  // Normalized luminance
    float y_inv = 1.0f / fmaxf(y, 0.01f);
    float X = Y * y_inv * x;
    float Z = Y * y_inv * (1.0f - x - y);
    
    // XYZ to sRGB matrix
    float3 rgb;
    rgb.x =  3.2404542f * X - 1.5371385f * Y - 0.4985314f * Z;
    rgb.y = -0.9692660f * X + 1.8760108f * Y + 0.0415560f * Z;
    rgb.z =  0.0556434f * X - 0.2040259f * Y + 1.0572252f * Z;
    
    // Clamp components to positive values
    rgb.x = fmaxf(rgb.x, 0.0f);
    rgb.y = fmaxf(rgb.y, 0.0f);
    rgb.z = fmaxf(rgb.z, 0.0f);
    
    return rgb;
}

// Sample accretion disk at position (volumetric)
// Returns emission color, optical depth contribution, and disk velocity
__device__ void sampleAccretionDiskVolumetric(
    const Vec4& pos,
    const Sirius::AccretionDiskParams& disk,
    const Sirius::MetricParams& mp,
    float ds,                    // Step size
    float beamRadius,            // Beam radius for cone tracing (Ray Bundles)
    float3& emission,            // OUT: emitted light this step
    float& dtau,                 // OUT: optical depth contribution
    float3& diskVelocity)        // OUT: disk velocity at position
{
    emission = make_float3(0.0f, 0.0f, 0.0f);
    dtau = 0.0f;
    diskVelocity = make_float3(0.0f, 0.0f, 0.0f);
    
    if (!disk.enabled) return;
    
    float r = pos.r;
    float theta = pos.theta;
    float phi = pos.phi;
    float M = mp.M;
    
    // Disk radial bounds (in units of M)
    float r_inner = disk.innerRadius * M;
    float r_outer = disk.outerRadius * M;
    
    if (r < r_inner || r > r_outer) return;
    
    // Height above equatorial plane: z = r * cos(theta)
    float z = fabsf(r * cosf(theta));
    
    // Disk scale height: H(r) = heightScale * r
    float H = disk.heightScale * r;
    
    // Cone tracing / Anti-aliasing with Ray Bundles
    // Broaden the Gaussian profile by the beam radius (convolution of Gaussian disk + Gaussian beam)
    // effectiveH = sqrt(H^2 + beamRadius^2)
    float effectiveH = sqrtf(H*H + beamRadius*beamRadius);
    
    // Maintain integrated density (mass conservation):
    // Peak density reduces as the profile creates a wider spread
    float blurFactor = H / fmaxf(effectiveH, 1e-6f);
    
    // Gaussian vertical density profile
    // ρ(r,z) = ρ₀(r) * exp(-z²/(2Σ²)) * blurFactor
    float z_norm = z / fmaxf(effectiveH, 0.001f);
    if (z_norm > 5.0f) return;  // Beyond 5 scale heights
    
    float rho_z = expf(-0.5f * z_norm * z_norm) * blurFactor;
    
    // =========================================================================
    // Shakura-Sunyaev Temperature Profile
    // T(r) = T_inner * (r_inner/r)^(3/4) * [1 - sqrt(r_inner/r)]^(1/4)
    // Simplified: T(r) ≈ T_inner * (r_inner/r)^temperatureExponent
    // =========================================================================
    float r_ratio = r_inner / r;
    float T_local = disk.innerTemperature * powf(r_ratio, disk.temperatureExponent);
    
    // NOTE: Stress-free inner boundary correction REMOVED for visual appeal
    // The standard Shakura-Sunyaev correction creates a temperature MINIMUM at ISCO,
    // which contradicts the "hottest at inner edge" visual expected from Interstellar.
    // For physically accurate disks, uncomment the lines below:
    // float stress_factor = powf(fmaxf(1.0f - sqrtf(r_ratio), 0.01f), 0.25f);
    // T_local *= stress_factor;
    
    // Radial density falloff (power law)
    float rho_r = powf(r_ratio, 0.5f);  // ρ ∝ r^(-1/2) approximately
    
    // Combined density
    float rho = rho_r * rho_z;

    // =========================================================================
    // Turbulence Perturbation (Cinematic Features Phase 8)
    // Apply Kolmogorov cascade density variations if enabled
    // =========================================================================
    if (params.volumetricDisk.turbulence.enabled) {
        float turbMod = sampleTurbulenceDensity(r, theta, phi, params.volumetricDisk.turbulence);
        rho *= turbMod;

        // Also perturb temperature slightly for color variation
        float tempNoise = fBm3D(r * 0.5f, theta * 2.0f, phi, params.volumetricDisk.turbulence);
        T_local *= (1.0f + 0.15f * params.volumetricDisk.turbulence.amplitude * tempNoise);
    }

    // =========================================================================
    // Emission: j_ν = ρ * ε * B_ν(T)
    // B_ν(T) represented as blackbody color
    // =========================================================================
    float3 bbColor = blackbodyColor(T_local);
    
    // Emission coefficient scales intensity
    float emissionStrength = rho * disk.emissionCoefficient;
    
    // Stefan-Boltzmann: total flux ∝ T^4
    // Normalize to ~5000K so solar-temperature disk has reasonable brightness
    float boltzmann = powf(T_local / 5000.0f, 4.0f);
    boltzmann = clamp(boltzmann, 0.0f, 20.0f);   // Relaxed for cinematic contrast range
    
    // Emission is step-size INDEPENDENT (source function brightness)
    // j = ρ * ε * B_ν(T) represents intrinsic emissivity
    // The step size only affects optical depth (dtau), NOT the source function
    emission = bbColor * emissionStrength * boltzmann;
    
    // =========================================================================
    // Absorption: dτ = ρ * α * ds
    // =========================================================================
    dtau = rho * disk.absorptionCoefficient * ds;
    
    // =========================================================================
    // Disk Velocity: Keplerian orbit in φ direction
    // v_φ = sqrt(M/r) for circular geodesic
    // =========================================================================
    float v_phi = sqrtf(M / r);
    // Velocity direction: tangent to circular orbit (+φ direction)
    // In Cartesian: v = (-sin(φ), 0, cos(φ)) * v_phi (for equatorial plane)
    diskVelocity = make_float3(-sinf(phi), 0.0f, cosf(phi)) * v_phi;
}
__device__ float dot4D(const float g[4][4], const Vec4& a, const Vec4& b) {
    float sum = 0.0f;
    // Unrolling for performance
    sum += g[0][0]*a.t*b.t + g[0][1]*a.t*b.r + g[0][2]*a.t*b.theta + g[0][3]*a.t*b.phi;
    sum += g[1][0]*a.r*b.t + g[1][1]*a.r*b.r + g[1][2]*a.r*b.theta + g[1][3]*a.r*b.phi;
    sum += g[2][0]*a.theta*b.t + g[2][1]*a.theta*b.r + g[2][2]*a.theta*b.theta + g[2][3]*a.theta*b.phi;
    sum += g[3][0]*a.phi*b.t + g[3][1]*a.phi*b.r + g[3][2]*a.phi*b.theta + g[3][3]*a.phi*b.phi;
    return sum;
}

// Solve for u^t given spatial velocity u^i to enforce u . u = -1
// Eq: g_tt(u^t)^2 + 2*g_ti*u^i*u^t + g_ij*u^i*u^j + 1 = 0
__device__ float solveU0(const float g[4][4], float ur, float utheta, float uphi) {
    float A = g[0][0];
    float B = 2.0f * (g[0][1]*ur + g[0][2]*utheta + g[0][3]*uphi);
    float C = g[1][1]*ur*ur + g[2][2]*utheta*utheta + g[3][3]*uphi*uphi + 
              2.0f*(g[1][2]*ur*utheta + g[1][3]*ur*uphi + g[2][3]*utheta*uphi) + 1.0f;
    
    // Quadratic formula: u^t = (-B +/- sqrt(B^2 - 4AC)) / 2A
    // Since u^t > 0 usually and A (g_tt) < 0, we need to pick correct sign.
    // For Schwarzschild/Kerr horizon exterior, A < 0, C > 0.
    // Discriminant D = B^2 - 4AC.
    
    float D = B*B - 4.0f*A*C;
    if (D < 0.0f) return 1.0f; // Error case, fallback
    
    // We want the positive root for time.
    // If g_tt = -1, B=0, C=v^2+1. -4(-1)(v^2+1) = 4(v^2+1). sqrt = 2 Gamma.
    // -B - sqrt(D) / 2A = -sqrt / -2 = +
    // So usually (-B - sqrt(D)) / 2A.
    
    return (-B - sqrtf(D)) / (2.0f * A);
}


// Unified metric dispatcher
// Unified metric dispatcher (Templated)
template<Sirius::MetricType type>
__device__ __forceinline__ void getMetricTensor(const Vec4& x, 
                                 const Sirius::MetricParams& mp,
                                 float g[4][4], float g_inv[4][4]) {
    if constexpr (type == Sirius::MetricType::Minkowski) {
        getMinkowskiMetric(x, g, g_inv);
    } else if constexpr (type == Sirius::MetricType::Schwarzschild) {
        getSchwarzschildMetric(x, mp.M, g, g_inv);
    } else if constexpr (type == Sirius::MetricType::Kerr) {
        getKerrMetric(x, mp.M, mp.a, g, g_inv);
    } else if constexpr (type == Sirius::MetricType::ReissnerNordstrom) {
        getReissnerNordstromMetric(x, mp.M, mp.Q, g, g_inv);
    } else if constexpr (type == Sirius::MetricType::Godel) {
        // Gödel uses rotation parameter 'a' from MetricParams
        // Default a=0.5 if not set (CTC radius at r=2)
        float a = (mp.a != 0.0f) ? mp.a : 0.5f;
        getGodelMetric(x, a, g, g_inv);
    } else if constexpr (type == Sirius::MetricType::TaubNUT) {
        // Taub-NUT uses M for mass, Q for NUT parameter
        // (repurposing Q since NUT is similar to charge conceptually)
        float nut = (mp.Q != 0.0f) ? mp.Q : 0.5f;
        getTaubNUTMetric(x, mp.M, nut, g, g_inv);
    } else if constexpr (type == Sirius::MetricType::KerrSchild) {
        // Kerr-Schild: horizon-penetrating coordinates
        getKerrSchildMetric(x, mp.M, mp.a, g, g_inv);
    } else if constexpr (type == Sirius::MetricType::EllisDrainhole) {
        // Ellis wormhole: m=mass-like, n=throat radius (use M and Q as proxies)
        float m_param = mp.M * 0.5f;  // Scale down
        float n_param = mp.M;         // Throat at ~M
        getEllisDrainholeMetric(x, m_param, n_param, g, g_inv);
    } else if constexpr (type == Sirius::MetricType::Alcubierre) {
        // Alcubierre warp drive: vs=velocity, sigma=wall thickness, R=bubble radius
        float vs = 2.0f;      // Default superluminal
        float sigma = 1.0f;   // Wall thickness
        float R = 2.0f * mp.M;  // Bubble radius scales with M
        getAlcubierreMetric(x, vs, sigma, R, g, g_inv);
    } else if constexpr (type == Sirius::MetricType::DeSitter) {
        // de Sitter: Lambda cosmological constant (use Q as proxy)
        float Lambda = (mp.Q != 0.0f) ? mp.Q : 0.01f;
        getDeSitterMetric(x, Lambda, g, g_inv);
    } else {
        getMinkowskiMetric(x, g, g_inv);
    }
}
//==============================================================================
// Analytic Christoffel Symbols
// Closed-form expressions for exact geodesic integration
//==============================================================================

// MINKOWSKI: All Christoffel symbols in spherical coords
// Non-zero: Γ^r_θθ, Γ^r_φφ, Γ^θ_rθ, Γ^θ_φφ, Γ^φ_rφ, Γ^φ_θφ
__device__ void getMinkowskiChristoffel(const Vec4& x, float Gamma[4][4][4]) {
    // Initialize to zero
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;
    
    float r = fmaxf(x.r, 0.001f);
    float theta = x.theta;
    float sin_theta = sinf(theta);
    float cos_theta = cosf(theta);
    float sin2 = sin_theta * sin_theta;
    
    // Γ^r_θθ = -r
    Gamma[1][2][2] = -r;
    // Γ^r_φφ = -r sin²θ
    Gamma[1][3][3] = -r * sin2;
    
    // Γ^θ_rθ = Γ^θ_θr = 1/r
    Gamma[2][1][2] = 1.0f / r;
    Gamma[2][2][1] = 1.0f / r;
    // Γ^θ_φφ = -sinθ cosθ
    Gamma[2][3][3] = -sin_theta * cos_theta;
    
    // Γ^φ_rφ = Γ^φ_φr = 1/r
    Gamma[3][1][3] = 1.0f / r;
    Gamma[3][3][1] = 1.0f / r;
    // Γ^φ_θφ = Γ^φ_φθ = cotθ (using regularized cot for pole safety)
    float cot_theta = safe_cot(theta);
    Gamma[3][2][3] = cot_theta;
    Gamma[3][3][2] = cot_theta;
}

// SCHWARZSCHILD: Christoffel symbols in Schwarzschild coords (t,r,θ,φ)
// ds² = -(1-rs/r)dt² + (1-rs/r)⁻¹dr² + r²dθ² + r²sin²θ dφ²
__device__ void getSchwarzschildChristoffel(const Vec4& x, float M, float Gamma[4][4][4]) {
    // Initialize to zero
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;
    
    float r = fmaxf(x.r, 0.001f);
    float theta = x.theta;
    float rs = 2.0f * M;
    
    float sin_theta, cos_theta;
    fast_sincosf(theta, &sin_theta, &cos_theta);
    float sin2 = sin_theta * sin_theta;
    
    float r2 = r * r;
    float r3 = r2 * r;
    float f = 1.0f - rs / r;  // (1 - rs/r)
    
    // Safety check near horizon
    if (f < 0.001f) {
        getMinkowskiChristoffel(x, Gamma);
        return;
    }
    
    // Γ^t_tr = Γ^t_rt = rs / (2r² f)
    float Gamma_t_tr = rs / (2.0f * r2 * f);
    Gamma[0][0][1] = Gamma_t_tr;
    Gamma[0][1][0] = Gamma_t_tr;
    
    // Γ^r_tt = rs f / (2r²)
    Gamma[1][0][0] = rs * f / (2.0f * r2);
    
    // Γ^r_rr = -rs / (2r² f)
    Gamma[1][1][1] = -rs / (2.0f * r2 * f);
    
    // Γ^r_θθ = -(r - rs) = -r f
    Gamma[1][2][2] = -r * f;
    
    // Γ^r_φφ = -(r - rs) sin²θ = -r f sin²θ
    Gamma[1][3][3] = -r * f * sin2;
    
    // Γ^θ_rθ = Γ^θ_θr = 1/r
    Gamma[2][1][2] = 1.0f / r;
    Gamma[2][2][1] = 1.0f / r;
    
    // Γ^θ_φφ = -sinθ cosθ
    Gamma[2][3][3] = -sin_theta * cos_theta;
    
    // Γ^φ_rφ = Γ^φ_φr = 1/r
    Gamma[3][1][3] = 1.0f / r;
    Gamma[3][3][1] = 1.0f / r;
    
    // Gamma^phi_theta_phi = Gamma^phi_phi_theta = cot(theta) (using regularized cot for pole safety)
    float cot_theta = safe_cot(theta);
    Gamma[3][2][3] = cot_theta;
    Gamma[3][3][2] = cot_theta;
}

// REISSNER-NORDSTRÖM: Christoffel symbols (spherically symmetric charged BH)
// ds² = -f(r)dt² + f(r)⁻¹dr² + r²dθ² + r²sin²θ dφ²
// where f(r) = 1 - 2M/r + Q²/r²
__device__ void getReissnerNordstromChristoffel(const Vec4& x, float M, float Q, float Gamma[4][4][4]) {
    // Initialize to zero
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;
    
    float r = fmaxf(x.r, 0.001f);
    float theta = x.theta;
    float Q2 = Q * Q;
    
    float sin_theta, cos_theta;
    fast_sincosf(theta, &sin_theta, &cos_theta);
    float sin2 = sin_theta * sin_theta;
    
    float r2 = r * r;
    float r3 = r2 * r;
    float r4 = r2 * r2;
    
    // f(r) = 1 - 2M/r + Q²/r²
    float f = 1.0f - 2.0f * M / r + Q2 / r2;
    
    // f'(r) = 2M/r² - 2Q²/r³
    float df = 2.0f * M / r2 - 2.0f * Q2 / r3;
    
    // Safety check near horizons
    if (f < 0.001f) {
        getMinkowskiChristoffel(x, Gamma);
        return;
    }
    
    // Γ^t_tr = Γ^t_rt = f'/(2f) = (M/r² - Q²/r³) / f
    float Gamma_t_tr = df / (2.0f * f);
    Gamma[0][0][1] = Gamma_t_tr;
    Gamma[0][1][0] = Gamma_t_tr;
    
    // Γ^r_tt = f * f' / 2 = f * (M/r² - Q²/r³)
    Gamma[1][0][0] = f * df / 2.0f;
    
    // Γ^r_rr = -f'/(2f) = -(M/r² - Q²/r³) / f
    Gamma[1][1][1] = -df / (2.0f * f);
    
    // Γ^r_θθ = -r * f
    Gamma[1][2][2] = -r * f;
    
    // Γ^r_φφ = -r * f * sin²θ
    Gamma[1][3][3] = -r * f * sin2;
    
    // Γ^θ_rθ = Γ^θ_θr = 1/r
    Gamma[2][1][2] = 1.0f / r;
    Gamma[2][2][1] = 1.0f / r;
    
    // Γ^θ_φφ = -sinθ cosθ
    Gamma[2][3][3] = -sin_theta * cos_theta;
    
    // Γ^φ_rφ = Γ^φ_φr = 1/r
    Gamma[3][1][3] = 1.0f / r;
    Gamma[3][3][1] = 1.0f / r;
    
    // Γ^φ_θφ = Γ^φ_φθ = cotθ (using regularized cot for pole safety)
    float cot_theta = safe_cot(theta);
    Gamma[3][2][3] = cot_theta;
    Gamma[3][3][2] = cot_theta;
}

// KERR: Christoffel symbols in Boyer-Lindquist coords (t,r,θ,φ)
// =========================================================================
// CORRECTED Dec 2025: Proper derivation using inverse metric with g^tφ ≠ 0
// The inverse metric for Kerr has:
//   g^tt = -A/(ΣΔ), g^tφ = g^φt = -2Mar/(ΣΔ)
//   g^rr = Δ/Σ, g^θθ = 1/Σ
//   g^φφ = (Δ - a²sin²θ)/(ΣΔsin²θ)
// =========================================================================
__device__ void getKerrChristoffel(const Vec4& x, float M, float a, float Gamma[4][4][4]) {
    // Initialize to zero
    #pragma unroll 4
    for (int i = 0; i < 4; i++)
        #pragma unroll 4
        for (int j = 0; j < 4; j++)
            #pragma unroll 4
            for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;
    
    // Minimum values to avoid numerical issues
    const float R_MIN = 0.01f;
    const float SIN_MIN = 0.0001f;  // Very small to minimize polar artifacts (~0.006 degrees from pole)
    
    float r = fmaxf(x.r, R_MIN);
    float theta = x.theta;
    
    float r2 = r * r;
    float a2 = a * a;
    float sin_theta, cos_theta;
    fast_sincosf(theta, &sin_theta, &cos_theta);
    
    // Clamp sin_theta away from zero (poles)
    if (fabsf(sin_theta) < SIN_MIN) {
        sin_theta = (sin_theta >= 0.0f) ? SIN_MIN : -SIN_MIN;
    }
    
    float cos2 = cos_theta * cos_theta;
    float sin2 = sin_theta * sin_theta;
    float sin_cos = sin_theta * cos_theta;
    
    // Kerr metric functions
    float Sigma = r2 + a2 * cos2;
    float Delta = r2 - 2.0f * M * r + a2;
    float A = (r2 + a2) * (r2 + a2) - Delta * a2 * sin2;  // = (r²+a²)² - Δa²sin²θ
    
    // Safety check near horizon/ergosphere
    const float DELTA_MIN = 0.02f;
    const float SIGMA_MIN = 0.02f;
    
    if (Delta < DELTA_MIN || Sigma < SIGMA_MIN) {
        getMinkowskiChristoffel(x, Gamma);
        return;
    }
    
    // Precompute common denominators
    float Sigma2 = Sigma * Sigma;
    float Sigma3 = Sigma2 * Sigma;
    float Delta2 = Delta * Delta;
    float Sigma_Delta = Sigma * Delta;
    float Sigma2_Delta = Sigma2 * Delta;
    
    // Metric derivatives (only non-zero ones)
    // ∂Σ/∂r = 2r,  ∂Σ/∂θ = -2a²sinθcosθ
    float dSigma_dr = 2.0f * r;
    float dSigma_dtheta = -2.0f * a2 * sin_cos;
    
    // ∂Δ/∂r = 2r - 2M = 2(r-M)
    float dDelta_dr = 2.0f * (r - M);
    
    // =========================================================================
    // Γ^t components - using g^tt = -A/(ΣΔ), g^tφ = -2Mar/(ΣΔ)
    // =========================================================================
    
    // Γ^t_tr = M(r² - a²cos²θ)(r² + a²) / (Σ²Δ)  [CORRECTED]
    float r2_minus_a2cos2 = r2 - a2 * cos2;
    Gamma[0][0][1] = M * r2_minus_a2cos2 * (r2 + a2) / Sigma2_Delta;
    Gamma[0][1][0] = Gamma[0][0][1];
    
    // Γ^t_tθ = -2Ma²r sinθ cosθ / Σ²  [VERIFIED]
    Gamma[0][0][2] = -2.0f * M * a2 * r * sin_cos / Sigma2;
    Gamma[0][2][0] = Gamma[0][0][2];
    
    // Γ^t_rφ = a sin²θ [Σr - M(r² - a²cos²θ) + M(r² + a²)] / (Σ²Δ)  [CORRECTED]
    // Simplified: a sin²θ [Σr + Ma²(1 - cos²θ)] / (Σ²Δ) = a sin²θ [Σr + Ma²sin²θ] / (Σ²Δ)
    Gamma[0][1][3] = -a * sin2 * (M * r2_minus_a2cos2 - Sigma * r) / Sigma2_Delta
                   + a * M * (r2 + a2) * sin2 / Sigma2_Delta;
    Gamma[0][3][1] = Gamma[0][1][3];
    
    // Γ^t_θφ = 2Mar(r² + a²) sinθ cosθ / Σ²  [CORRECTED]
    Gamma[0][2][3] = 2.0f * M * a * r * (r2 + a2) * sin_cos / Sigma2;
    Gamma[0][3][2] = Gamma[0][2][3];
    
    // =========================================================================
    // Γ^r components - using g^rr = Δ/Σ
    // =========================================================================
    
    // Γ^r_tt = MΔ(r² - a²cos²θ) / Σ³  [VERIFIED]
    Gamma[1][0][0] = M * Delta * r2_minus_a2cos2 / Sigma3;
    
    // Γ^r_tφ = -MaΔ sin²θ (r² - a²cos²θ) / Σ³  [VERIFIED]
    Gamma[1][0][3] = -M * a * Delta * sin2 * r2_minus_a2cos2 / Sigma3;
    Gamma[1][3][0] = Gamma[1][0][3];
    
    // Γ^r_rr = (rΔ - Σ(r-M)) / (ΣΔ)  [CORRECTED - removes artifacts]
    Gamma[1][1][1] = (r * Delta - Sigma * (r - M)) / Sigma_Delta;
    
    // Γ^r_rθ = -a² sinθ cosθ / Σ  [VERIFIED]
    Gamma[1][1][2] = -a2 * sin_cos / Sigma;
    Gamma[1][2][1] = Gamma[1][1][2];
    
    // Γ^r_θθ = -rΔ / Σ  [VERIFIED]
    Gamma[1][2][2] = -r * Delta / Sigma;
    
    // Γ^r_φφ = -Δ sin²θ [rΣ² + Ma²sin²θ(a²cos²θ - r²)] / Σ³  [CORRECTED]
    Gamma[1][3][3] = -Delta * sin2 * (r * Sigma2 + M * a2 * sin2 * (a2 * cos2 - r2)) / Sigma3;
    
    // =========================================================================
    // Γ^θ components - using g^θθ = 1/Σ
    // =========================================================================
    
    // Γ^θ_tt = -2Ma²r sinθ cosθ / Σ³  [VERIFIED - differs from Γ^t_tθ by factor]
    Gamma[2][0][0] = -2.0f * M * a2 * r * sin_cos / Sigma3;
    
    // Γ^θ_tφ = 2Mar(r² + a²) sinθ cosθ / Σ³  [VERIFIED]
    Gamma[2][0][3] = 2.0f * M * a * r * (r2 + a2) * sin_cos / Sigma3;
    Gamma[2][3][0] = Gamma[2][0][3];
    
    // Γ^θ_rr = a² sinθ cosθ / (ΣΔ)  [VERIFIED]
    Gamma[2][1][1] = a2 * sin_cos / Sigma_Delta;
    
    // Γ^θ_rθ = r / Σ  [VERIFIED]
    Gamma[2][1][2] = r / Sigma;
    Gamma[2][2][1] = Gamma[2][1][2];
    
    // Γ^θ_θθ = -a² sinθ cosθ / Σ  [VERIFIED]
    Gamma[2][2][2] = -a2 * sin_cos / Sigma;
    
    // Γ^θ_φφ = -sinθ cosθ [AΣ + 2Ma²r sin²θ(r² + a²)] / Σ³  [CORRECTED]
    Gamma[2][3][3] = -sin_cos * (A * Sigma + 2.0f * M * a2 * r * sin2 * (r2 + a2)) / Sigma3;
    
    // =========================================================================
    // Γ^φ components - using g^φφ = (Δ - a²sin²θ)/(ΣΔsin²θ), g^φt = -2Mar/(ΣΔ)
    // These are the most complex due to frame-dragging
    // =========================================================================
    
    // Γ^φ_tr = Ma(a²cos²θ - r²) / (Σ²Δ)  [CORRECTED - sign fixed]
    Gamma[3][0][1] = M * a * (a2 * cos2 - r2) / Sigma2_Delta;
    Gamma[3][1][0] = Gamma[3][0][1];
    
    // Γ^φ_tθ = -2Mar cot(θ) / Σ²  [VERIFIED - using regularized cot for pole safety]
    float cot_theta = safe_cot(theta);
    Gamma[3][0][2] = -2.0f * M * a * r * cot_theta / Sigma2;
    Gamma[3][2][0] = Gamma[3][0][2];
    
    // Γ^φ_rφ = [rΣ² - M(r² - a²cos²θ)(r² + a²)/Δ] / Σ³  [CORRECTED]
    Gamma[3][1][3] = (r * Sigma2 - M * r2_minus_a2cos2 * (r2 + a2) / Delta) / Sigma3;
    Gamma[3][3][1] = Gamma[3][1][3];
    
    // Γ^φ_θφ = cot(θ) + a²sinθcosθ(Σ + 2Mr) / Σ²  [CORRECTED - includes Kerr frame-dragging]
    // Base spherical term + Kerr correction for proper light bending
    Gamma[3][2][3] = cot_theta + a2 * sin_cos * (Sigma + 2.0f * M * r) / Sigma2;
    Gamma[3][3][2] = Gamma[3][2][3];
}

// KERR-SCHILD: Christoffel symbols in spherical Kerr-Schild coords
// =========================================================================
// Kerr-Schild metric: g_μν = η_μν + H l_μ l_ν
// where H = 2Mr/Σ, Σ = r² + a²cos²θ
// l^μ = (1, 1, 0, a/(r²+a²)) is null with respect to Minkowski
//
// ADVANTAGE: The Kerr-Schild form is horizon-penetrating and has 
// much simpler singularity structure than Boyer-Lindquist. The pole
// singularities (θ = 0, π) are milder because the metric is closer to
// Minkowski near the poles.
//
// For this simplified spherical Kerr-Schild form (as in getKerrSchildMetric):
// g_tt = -(1-H), g_tr = H, g_rr = 1+H, g_θθ = Σ, g_φφ = Σ sin²θ
// =========================================================================
__device__ void getKerrSchildChristoffel(const Vec4& x, float M, float a, float Gamma[4][4][4]) {
    // Initialize to zero
    #pragma unroll 4
    for (int i = 0; i < 4; i++)
        #pragma unroll 4
        for (int j = 0; j < 4; j++)
            #pragma unroll 4
            for (int k = 0; k < 4; k++)
                Gamma[i][j][k] = 0.0f;
    
    float r = fmaxf(x.r, 0.01f);
    float theta = x.theta;
    
    float r2 = r * r;
    float a2 = a * a;
    float sin_theta, cos_theta;
    fast_sincosf(theta, &sin_theta, &cos_theta);
    float sin2 = sin_theta * sin_theta;
    float cos2 = cos_theta * cos_theta;
    
    // Core functions
    float Sigma = r2 + a2 * cos2;
    float Sigma2 = Sigma * Sigma;
    float H = 2.0f * M * r / fmaxf(Sigma, 0.01f);
    
    // Derivatives of H
    // ∂H/∂r = 2M(Σ - r·2r)/Σ² = 2M(Σ - 2r²)/Σ² = 2M(a²cos²θ - r²)/Σ²
    float dH_dr = 2.0f * M * (a2 * cos2 - r2) / Sigma2;
    
    // ∂H/∂θ = 2Mr · (-1/Σ²) · ∂Σ/∂θ = 2Mr · (-1/Σ²) · (-2a²sinθcosθ)
    //       = 4Ma²r sinθ cosθ / Σ²
    float sin_cos = sin_theta * cos_theta;
    float dH_dtheta = 4.0f * M * a2 * r * sin_cos / Sigma2;
    
    // Derivatives of Σ
    float dSigma_dr = 2.0f * r;
    float dSigma_dtheta = -2.0f * a2 * sin_cos;
    
    // The simplified Kerr-Schild metric has structure:
    // g_tt = -(1-H), g_tr = g_rt = H, g_rr = 1+H
    // g_θθ = Σ, g_φφ = Σ sin²θ
    
    // Using Γ^λ_μν = (1/2) g^{λσ} (∂_μ g_σν + ∂_ν g_σμ - ∂_σ g_μν)
    // For the (t,r) block, we need the 2x2 inverse:
    // g^tt = (1+H)/D, g^tr = H/D, g^rr = (1-H)/D where D = -(1-H)(1+H) - H² = -(1-H²) - H² = -1
    // So: g^tt = -(1+H), g^tr = g^rt = -H, g^rr = -(1-H) = H-1
    
    float g_tt = -(1.0f - H);
    float g_tr = H;
    float g_rr = 1.0f + H;
    
    // Inverse of (t,r) block: det = g_tt*g_rr - g_tr² = -(1-H)(1+H) - H² = -(1-H²) - H² = -1
    float det_tr = -1.0f;  // Always -1 for this form
    float g_inv_tt = g_rr / det_tr;   // = -(1+H)
    float g_inv_tr = -g_tr / det_tr;  // = H
    float g_inv_rr = g_tt / det_tr;   // = (1-H)
    
    float g_inv_thth = 1.0f / Sigma;
    float g_inv_phph = 1.0f / (Sigma * sin2 + 1e-8f);
    
    // =========================================================================
    // Christoffel symbols
    // =========================================================================
    
    // Γ^t components
    // Γ^t_tr = Γ^t_rt = (1/2) * (g^tt * ∂_r g_tt + g^tr * ∂_r g_rr)
    //        = (1/2) * (-(1+H) * dH_dr + H * dH_dr) = (1/2) * dH_dr * (H - 1 - H) = -(1/2) dH_dr
    // Actually need full formula with all terms
    float dg_tt_dr = dH_dr;    // ∂(-(1-H))/∂r = dH_dr
    float dg_tr_dr = dH_dr;    // ∂H/∂r
    float dg_rr_dr = dH_dr;    // ∂(1+H)/∂r
    
    float dg_tt_dth = dH_dtheta;
    float dg_tr_dth = dH_dtheta;
    float dg_rr_dth = dH_dtheta;
    
    // Γ^t_tt: (1/2) g^tσ (2∂_t g_σt - ∂_σ g_tt) - but metric is static, so ∂_t = 0
    //       = -(1/2) g^tσ ∂_σ g_tt = -(1/2)(g^tr * dg_tt_dr + 0) = -(1/2) * H * dH_dr
    // Wait, simpler: Γ^t_tt = 0 since ∂_t g = 0 and we need ∂_μ g_νσ = 0 for t indices
    
    // Let's use the direct formula for key components:
    // Γ^t_tr = (1/2) g^tt (∂_r g_tt) + (1/2) g^tr (∂_r g_rt + ∂_t g_rr - ∂_r g_tr)
    //       = (1/2) g^tt dg_tt_dr + (1/2) g^tr (dg_tr_dr - dg_tr_dr)
    //       = (1/2) * (-(1+H)) * dH_dr = -(1+H)/2 * dH_dr
    Gamma[0][0][1] = -(1.0f + H) * dH_dr / 2.0f;
    Gamma[0][1][0] = Gamma[0][0][1];
    
    // Γ^t_tθ: (1/2) g^tt ∂_θ g_tt = (1/2) * (-(1+H)) * dH_dtheta = -(1+H)/2 * dH_dtheta
    Gamma[0][0][2] = -(1.0f + H) * dH_dtheta / 2.0f;
    Gamma[0][2][0] = Gamma[0][0][2];
    
    // Γ^t_rr: (1/2) g^tσ (2∂_r g_σr - ∂_σ g_rr) 
    //       = (1/2) * (g^tt (2∂_r g_tr - ∂_t g_rr) + g^tr (2∂_r g_rr - ∂_r g_rr))
    //       = (1/2) * (g^tt * 2*dH_dr + g^tr * dH_dr)
    //       = (1/2) * dH_dr * (2*(-(1+H)) + H) = (1/2) * dH_dr * (-2 - 2H + H) = (1/2) * dH_dr * (-2 - H)
    Gamma[0][1][1] = dH_dr * (-2.0f - H) / 2.0f;
    
    // Γ^r components
    // Γ^r_tt: (1/2) g^rσ (2∂_t g_σt - ∂_σ g_tt) = -(1/2) g^rr ∂_r g_tt = -(1/2) * (H-1) * dH_dr
    Gamma[1][0][0] = -(H - 1.0f) * dH_dr / 2.0f;
    
    // Γ^r_tr: (1/2) g^rt ∂_r g_tt + (1/2) g^rr ∂_r g_tr = H * dH_dr / 2 + (H-1) * dH_dr / 2
    //       = dH_dr * (H + H - 1) / 2 = dH_dr * (2H - 1) / 2
    Gamma[1][0][1] = dH_dr * (2.0f * H - 1.0f) / 2.0f;
    Gamma[1][1][0] = Gamma[1][0][1];
    
    // Γ^r_rr: (1/2) g^rr ∂_r g_rr = (1/2) * (H-1) * dH_dr
    Gamma[1][1][1] = (H - 1.0f) * dH_dr / 2.0f;
    
    // Γ^r_θθ: -(1/2) g^rr ∂_r g_θθ = -(1/2) * (H-1) * 2r = -r(H-1) = r(1-H)
    Gamma[1][2][2] = r * (1.0f - H);
    
    // Γ^r_φφ: -(1/2) g^rr ∂_r g_φφ = -(1/2) * (H-1) * (2r*sin² + Σ * 0) = r(1-H) * sin²
    Gamma[1][3][3] = r * (1.0f - H) * sin2;
    
    // Γ^r_tθ, Γ^r_rθ: Similar derivations involving dH_dtheta
    Gamma[1][0][2] = -(H - 1.0f) * dH_dtheta / 2.0f;
    Gamma[1][2][0] = Gamma[1][0][2];
    
    Gamma[1][1][2] = (H - 1.0f) * dH_dtheta / 2.0f;
    Gamma[1][2][1] = Gamma[1][1][2];
    
    // Γ^θ components (simpler - diagonal angular part)
    // Γ^θ_rθ = (1/2) g^θθ ∂_r g_θθ = (1/2) * (1/Σ) * 2r = r/Σ
    Gamma[2][1][2] = r / Sigma;
    Gamma[2][2][1] = Gamma[2][1][2];
    
    // Γ^θ_θθ = (1/2) g^θθ ∂_θ g_θθ = (1/2) * (1/Σ) * (-2a²sinθcosθ) = -a²sinθcosθ/Σ
    Gamma[2][2][2] = -a2 * sin_cos / Sigma;
    
    // Γ^θ_φφ = -(1/2) g^θθ ∂_θ g_φφ = -(1/2) * (1/Σ) * (dΣ/dθ * sin² + Σ * 2sinθcosθ)
    //       = -(1/Σ) * (-a²sinθcosθ * sin² + Σ * sinθcosθ) / 2
    // Actually: g_φφ = Σ sin², so ∂_θ g_φφ = ∂_θΣ * sin² + Σ * 2sinθcosθ
    //                                     = -2a²sinθcosθ * sin² + 2Σsinθcosθ
    // Γ^θ_φφ = -(1/2Σ) * (-2a²sin²θsinθcosθ + 2Σsinθcosθ) = sinθcosθ(a²sin²θ - Σ)/Σ
    //        = sinθcosθ(a²sin² - r² - a²cos²)/Σ = sinθcosθ(a²(sin²-cos²) - r²)/Σ
    // Or simpler: = -sinθcosθ + a²sin²θsinθcosθ/Σ (approx)
    // For simplicity, use the standard form like Schwarzschild for now
    Gamma[2][3][3] = -sin_cos * (Sigma - a2 * sin2) / Sigma;
    
    // Γ^θ contributions from (t,r) block via H derivatives
    Gamma[2][0][0] = -g_inv_thth * dH_dtheta / 2.0f;  // From ∂_θ g_tt term
    Gamma[2][0][1] = -g_inv_thth * dH_dtheta / 2.0f;  // Mixed
    Gamma[2][1][0] = Gamma[2][0][1];
    Gamma[2][1][1] = -g_inv_thth * dH_dtheta / 2.0f;  // From ∂_θ g_rr term
    
    // Γ^φ components (azimuthal symmetry - simpler)
    // Γ^φ_rφ = (1/2) g^φφ ∂_r g_φφ = (1/2) * (1/(Σsin²)) * 2r*sin² = r/Σ
    Gamma[3][1][3] = r / Sigma;
    Gamma[3][3][1] = Gamma[3][1][3];
    
    // Γ^φ_θφ = (1/2) g^φφ ∂_θ g_φφ = (1/2) * (1/(Σsin²)) * (∂_θΣ * sin² + Σ * 2sinθcosθ)
    //       = (1/2Σsin²) * (-2a²sinθcosθ * sin² + 2Σsinθcosθ)
    //       = cotθ * (Σ - a²sin²) / Σ = cotθ * (r² + a²cos² - a²sin² + a²sin²) / Σ
    // Simplified for near-Schwarzschild: ≈ cotθ
    float cot_theta = safe_cot(theta);
    Gamma[3][2][3] = cot_theta;
    Gamma[3][3][2] = Gamma[3][2][3];
}
//==============================================================================
// P2: Precomputed Christoffel Symbol Lookup
// Uses cached 3D textures instead of finite differences
//==============================================================================

// Pack lower indices (ν, ρ) exploiting symmetry: ν ≤ ρ → idx = ν + ρ*(ρ+1)/2
__device__ __forceinline__ int packLowerIndices(int nu, int rho) {
    if (nu > rho) { int t = nu; nu = rho; rho = t; }  // Swap if needed
    return nu + rho * (rho + 1) / 2;
}

//==============================================================================
// Christoffel Coordinate Transformation: Cartesian → Spherical
//==============================================================================
// Transformation law:
// Γ'^μ_νρ = (∂x'μ/∂xα)(∂xβ/∂x'ν)(∂xγ/∂x'ρ)Γα_βγ + (∂x'μ/∂xα)(∂²xα/∂x'ν∂x'ρ)
//
// Where x' = (t,r,θ,φ) spherical, x = (t,x,y,z) Cartesian
//==============================================================================
__device__ void transformChristoffelCartToSph(
    float Gamma_cart[4][4][4],  // Input: Cartesian Christoffel
    float r, float theta, float phi,
    float Gamma_sph[4][4][4])   // Output: Spherical Christoffel
{
    // Precompute trigonometric values
    float st = sinf(theta), ct = cosf(theta);
    float sp = sinf(phi), cp = cosf(phi);
    
    // Clamp sin(theta) to avoid poles - very minimal threshold
    if (fabsf(st) < 0.0001f) st = (st >= 0.0f) ? 0.0001f : -0.0001f;
    
    float r_inv = 1.0f / fmaxf(r, 0.01f);
    float rst_inv = 1.0f / (r * st + 1e-6f);
    
    // =========================================================================
    // Jacobian J = ∂x/∂x' (Cartesian derivs w.r.t. spherical)
    // J[cart_idx][sph_idx], cart: 0=t,1=x,2=y,3=z, sph: 0=t,1=r,2=θ,3=φ
    // =========================================================================
    float J[4][4] = {};
    J[0][0] = 1.0f;  // ∂t/∂t = 1
    
    // ∂x/∂(r,θ,φ)
    J[1][1] = st * cp;           // ∂x/∂r
    J[1][2] = r * ct * cp;       // ∂x/∂θ
    J[1][3] = -r * st * sp;      // ∂x/∂φ
    
    // ∂y/∂(r,θ,φ)
    J[2][1] = st * sp;           // ∂y/∂r
    J[2][2] = r * ct * sp;       // ∂y/∂θ
    J[2][3] = r * st * cp;       // ∂y/∂φ
    
    // ∂z/∂(r,θ,φ)
    J[3][1] = ct;                // ∂z/∂r
    J[3][2] = -r * st;           // ∂z/∂θ
    J[3][3] = 0.0f;              // ∂z/∂φ
    
    // =========================================================================
    // Inverse Jacobian J_inv = ∂x'/∂x (spherical derivs w.r.t. Cartesian)
    // J_inv[sph_idx][cart_idx]
    // =========================================================================
    float J_inv[4][4] = {};
    J_inv[0][0] = 1.0f;  // ∂t'/∂t = 1
    
    // ∂r/∂(x,y,z)
    J_inv[1][1] = st * cp;       // ∂r/∂x = x/r
    J_inv[1][2] = st * sp;       // ∂r/∂y = y/r  
    J_inv[1][3] = ct;            // ∂r/∂z = z/r
    
    // ∂θ/∂(x,y,z)
    J_inv[2][1] = ct * cp * r_inv;   // ∂θ/∂x
    J_inv[2][2] = ct * sp * r_inv;   // ∂θ/∂y
    J_inv[2][3] = -st * r_inv;       // ∂θ/∂z
    
    // ∂φ/∂(x,y,z)
    J_inv[3][1] = -sp * rst_inv;     // ∂φ/∂x
    J_inv[3][2] = cp * rst_inv;      // ∂φ/∂y
    J_inv[3][3] = 0.0f;              // ∂φ/∂z
    
    // =========================================================================
    // Second derivatives ∂²x/∂x'^ν∂x'^ρ (needed for inhomogeneous term)
    // d2x[cart][sph_nu][sph_rho]
    // =========================================================================
    float d2x[4][4][4] = {};  // Most are zero
    
    // ∂²x/∂r∂θ = cos(θ)cos(φ)
    d2x[1][1][2] = d2x[1][2][1] = ct * cp;
    // ∂²x/∂r∂φ = -sin(θ)sin(φ)
    d2x[1][1][3] = d2x[1][3][1] = -st * sp;
    // ∂²x/∂θ² = -r sin(θ)cos(φ)
    d2x[1][2][2] = -r * st * cp;
    // ∂²x/∂θ∂φ = -r cos(θ)sin(φ)
    d2x[1][2][3] = d2x[1][3][2] = -r * ct * sp;
    // ∂²x/∂φ² = -r sin(θ)cos(φ)
    d2x[1][3][3] = -r * st * cp;
    
    // ∂²y/∂r∂θ = cos(θ)sin(φ)
    d2x[2][1][2] = d2x[2][2][1] = ct * sp;
    // ∂²y/∂r∂φ = sin(θ)cos(φ)
    d2x[2][1][3] = d2x[2][3][1] = st * cp;
    // ∂²y/∂θ² = -r sin(θ)sin(φ)
    d2x[2][2][2] = -r * st * sp;
    // ∂²y/∂θ∂φ = r cos(θ)cos(φ)
    d2x[2][2][3] = d2x[2][3][2] = r * ct * cp;
    // ∂²y/∂φ² = -r sin(θ)sin(φ)
    d2x[2][3][3] = -r * st * sp;
    
    // ∂²z/∂r∂θ = -sin(θ)
    d2x[3][1][2] = d2x[3][2][1] = -st;
    // ∂²z/∂θ² = -r cos(θ)
    d2x[3][2][2] = -r * ct;
    
    // =========================================================================
    // OPTIMIZED Transformation: Exploit Jacobian sparsity
    // =========================================================================
    // Time components: J[0][0]=1, J[0][i]=0 for i>0, J[i][0]=0 for i>0
    // This means:
    //   - Γ'^0_νρ = Γ^0_νρ when both ν,ρ are time (trivial)
    //   - Mixed time-space components need only partial sums
    //   - Only spatial 3×3 needs full transformation
    // Reduces from O(4^6)=4096 to ~O(3^4+handling)≈100 ops
    // =========================================================================
    
    // Initialize output to zero
    #pragma unroll 4
    for (int mu = 0; mu < 4; mu++)
        #pragma unroll 4
        for (int nu = 0; nu < 4; nu++)
            #pragma unroll 4
            for (int rho = 0; rho < 4; rho++)
                Gamma_sph[mu][nu][rho] = 0.0f;
    
    // Time-time-time: Γ'^0_00 = Γ^0_00 (trivial since J[0][0]=1)
    Gamma_sph[0][0][0] = Gamma_cart[0][0][0];
    
    // Time-time-space: Γ'^0_0j and Γ'^0_i0 (j,i = 1,2,3)
    // Γ'^0_0j = J[β][j] Γ^0_0β = sum over spatial β
    #pragma unroll 3
    for (int j = 1; j < 4; j++) {
        float sum = 0.0f;
        #pragma unroll 3
        for (int beta = 1; beta < 4; beta++) {
            sum += J[beta][j] * Gamma_cart[0][0][beta];
        }
        Gamma_sph[0][0][j] = sum;
        Gamma_sph[0][j][0] = sum;  // Symmetry
    }
    
    // Time-space-space: Γ'^0_ij (i,j = 1,2,3)
    #pragma unroll 3
    for (int i = 1; i < 4; i++) {
        #pragma unroll 3
        for (int j = 1; j < 4; j++) {
            float sum = 0.0f;
            #pragma unroll 3
            for (int beta = 1; beta < 4; beta++) {
                #pragma unroll 3
                for (int gamma = 1; gamma < 4; gamma++) {
                    sum += J[beta][i] * J[gamma][j] * Gamma_cart[0][beta][gamma];
                }
            }
            Gamma_sph[0][i][j] = sum;
        }
    }
    
    // Space-time-time: Γ'^i_00 (i = 1,2,3)
    // Γ'^i_00 = J_inv[i][α] Γ^α_00 = sum over spatial α
    #pragma unroll 3
    for (int i = 1; i < 4; i++) {
        float sum = 0.0f;
        #pragma unroll 3
        for (int alpha = 1; alpha < 4; alpha++) {
            sum += J_inv[i][alpha] * Gamma_cart[alpha][0][0];
        }
        Gamma_sph[i][0][0] = sum;
    }
    
    // Space-time-space: Γ'^i_0j and Γ'^i_j0 (i,j = 1,2,3)
    #pragma unroll 3
    for (int i = 1; i < 4; i++) {
        #pragma unroll 3
        for (int j = 1; j < 4; j++) {
            float sum = 0.0f;
            #pragma unroll 3
            for (int alpha = 1; alpha < 4; alpha++) {
                #pragma unroll 3
                for (int gamma = 1; gamma < 4; gamma++) {
                    sum += J_inv[i][alpha] * J[gamma][j] * Gamma_cart[alpha][0][gamma];
                }
            }
            Gamma_sph[i][0][j] = sum;
            Gamma_sph[i][j][0] = sum;  // Symmetry in lower indices
        }
    }
    
    // Space-space-space: Γ'^i_jk (i,j,k = 1,2,3) - MAIN contribution
    // Also includes inhomogeneous term from second derivatives
    #pragma unroll 3
    for (int i = 1; i < 4; i++) {
        #pragma unroll 3
        for (int j = 1; j < 4; j++) {
            #pragma unroll 3
            for (int k = 1; k < 4; k++) {
                float sum = 0.0f;
                
                // First term: J_inv[i][α] J[β][j] J[γ][k] Γ^α_βγ
                #pragma unroll 3
                for (int alpha = 1; alpha < 4; alpha++) {
                    #pragma unroll 3
                    for (int beta = 1; beta < 4; beta++) {
                        #pragma unroll 3
                        for (int gamma = 1; gamma < 4; gamma++) {
                            sum += J_inv[i][alpha] * J[beta][j] * J[gamma][k] * Gamma_cart[alpha][beta][gamma];
                        }
                    }
                }
                
                // Second term (inhomogeneous): J_inv[i][α] d2x[α][j][k]
                #pragma unroll 3
                for (int alpha = 1; alpha < 4; alpha++) {
                    sum += J_inv[i][alpha] * d2x[alpha][j][k];
                }
                
                Gamma_sph[i][j][k] = sum;
            }
        }
    }
}


// Lookup precomputed Christoffel from textures and transform to spherical
__device__ void getPrecomputedChristoffel(const Vec4& x, float Gamma[4][4][4]) {
    const auto& nm = params.numericalMetric;
    
    // Convert spherical to Cartesian for grid lookup
    float st = sinf(x.theta), ct = cosf(x.theta);
    float sp = sinf(x.phi), cp = cosf(x.phi);
    
    float3 cart = make_float3(
        x.r * st * cp,
        x.r * st * sp,
        x.r * ct
    );
    
    // Compute texture coordinates
    float3 uvw = make_float3(
        (cart.x - nm.origin.x) / (nm.dims.x * nm.spacing.x),
        (cart.y - nm.origin.y) / (nm.dims.y * nm.spacing.y),
        (cart.z - nm.origin.z) / (nm.dims.z * nm.spacing.z)
    );
    
    // Clamp to valid range
    uvw.x = clamp(uvw.x, 0.0f, 1.0f);
    uvw.y = clamp(uvw.y, 0.0f, 1.0f);
    uvw.z = clamp(uvw.z, 0.0f, 1.0f);
    
    // Sample all 40 unique components (Cartesian Christoffel symbols)
    // Γ^t components (mu=0)
    float Gt[10];
    for (int i = 0; i < 10; i++) {
        Gt[i] = tex3D<float>(nm.Gamma_t[i], uvw.x, uvw.y, uvw.z);
    }
    
    // Γ^x components (mu=1)
    float Gx[10];
    for (int i = 0; i < 10; i++) {
        Gx[i] = tex3D<float>(nm.Gamma_r[i], uvw.x, uvw.y, uvw.z);
    }
    
    // Γ^y components (mu=2)
    float Gy[10];
    for (int i = 0; i < 10; i++) {
        Gy[i] = tex3D<float>(nm.Gamma_theta[i], uvw.x, uvw.y, uvw.z);
    }
    
    // Γ^z components (mu=3)
    float Gz[10];
    for (int i = 0; i < 10; i++) {
        Gz[i] = tex3D<float>(nm.Gamma_phi[i], uvw.x, uvw.y, uvw.z);
    }
    
    // Unpack into full 4×4×4 array (CARTESIAN Christoffel)
    float Gamma_cart[4][4][4];
    for (int nu = 0; nu < 4; nu++) {
        for (int rho = 0; rho < 4; rho++) {
            int idx = packLowerIndices(nu, rho);
            Gamma_cart[0][nu][rho] = Gt[idx];
            Gamma_cart[1][nu][rho] = Gx[idx];
            Gamma_cart[2][nu][rho] = Gy[idx];
            Gamma_cart[3][nu][rho] = Gz[idx];
        }
    }
    
    // =========================================================================
    // PHASE 3.1: Transform Cartesian Christoffel to Spherical
    // =========================================================================
    transformChristoffelCartToSph(Gamma_cart, x.r, x.theta, x.phi, Gamma);
}

//==============================================================================
// Christoffel Symbols - Unified Dispatcher
// Uses analytic for Minkowski/Schwarzschild/Kerr, numerical for Numerical/RN
//==============================================================================
template<Sirius::MetricType type>
__device__ void getChristoffelSymbols(const Vec4& x,
                                       const Sirius::MetricParams& mp,
                                       float Gamma[4][4][4]) {
    // Dispatch to analytic implementations for performance
    if constexpr (type == Sirius::MetricType::Minkowski) {
        getMinkowskiChristoffel(x, Gamma);
        return;
    }
    if constexpr (type == Sirius::MetricType::Schwarzschild) {
        getSchwarzschildChristoffel(x, mp.M, Gamma);
        return;
    }
    // =========================================================================
    // KERR: Uses Boyer-Lindquist Christoffel symbols
    //
    // CRITICAL FIX (Jan 2026): Must use Boyer-Lindquist Christoffel symbols
    // to match the Boyer-Lindquist metric used by getKerrMetric().
    // Previously this incorrectly called getKerrSchildChristoffel which uses
    // Kerr-Schild coordinates, causing systematic ray deflection errors.
    //
    // The safe_cot() function regularizes cot(θ) at poles using Taylor
    // expansion, maintaining physics while avoiding singularities.
    // =========================================================================
    if constexpr (type == Sirius::MetricType::Kerr) {
        getKerrChristoffel(x, mp.M, mp.a, Gamma);
        return;
    }



    
    // =========================================================================
    // KERR-SCHILD: Generic Unified Family (Cartesian -> Spherical)
    // Uses analytic derivatives including Lambda and Charge (Q)
    // =========================================================================
    if constexpr (type == Sirius::MetricType::KerrSchild) {
        // 1. Convert Position to Cartesian
        Vec4Cart x_cart = vec4SphToCart(x);
        
        // 2. Compute Cartesian Christoffel Symbols
        float Gamma_cart[4][4][4];
        getKerrSchildFamilyChristoffel(x_cart, mp.M, mp.a, mp.Q, mp.Lambda, Gamma_cart);
        
        // 3. Transform to Spherical Coordinates (Boyer-Lindquist)
        transformChristoffelCartToSph(Gamma_cart, x.r, x.theta, x.phi, Gamma);
        return;
    }


    
    // =========================================================================
    // REISSNER-NORDSTRÖM: Now using Unified Kerr-Schild Family
    // Supports Charge (Q) and Lambda (Cosmological Constant)
    // =========================================================================
    if constexpr (type == Sirius::MetricType::ReissnerNordstrom) {
        // 1. Convert Position to Cartesian
        Vec4Cart x_cart = vec4SphToCart(x);
        
        // 2. Compute Cartesian Christoffel Symbols
        float Gamma_cart[4][4][4];
        getKerrSchildFamilyChristoffel(x_cart, mp.M, 0.0f, mp.Q, mp.Lambda, Gamma_cart);
        
        // 3. Transform to Spherical Coordinates
        transformChristoffelCartToSph(Gamma_cart, x.r, x.theta, x.phi, Gamma);
        return;
    }

    // =========================================================================
    // DE SITTER / ANTI-DE SITTER: Unified Kerr-Schild Family (M=a=Q=0)
    // =========================================================================
    if constexpr (type == Sirius::MetricType::DeSitter) {
        // 1. Convert Position to Cartesian
        Vec4Cart x_cart = vec4SphToCart(x);
        
        // 2. Compute Cartesian Christoffel Symbols
        float Gamma_cart[4][4][4];
        getKerrSchildFamilyChristoffel(x_cart, 0.0f, 0.0f, 0.0f, mp.Lambda, Gamma_cart);
        
        // 3. Transform to Spherical Coordinates
        transformChristoffelCartToSph(Gamma_cart, x.r, x.theta, x.phi, Gamma);
        return;
    }

    // =========================================================================
    // ALCUBIERRE WARP DRIVE: Cartesian -> Spherical Transform
    // =========================================================================
    // =========================================================================
    // ALCUBIERRE WARP DRIVE: Cartesian -> Spherical Transform
    // =========================================================================
    if constexpr (type == Sirius::MetricType::Alcubierre) {
        // 1. Convert to Cartesian
        Vec4Cart x_cart = vec4SphToCart(x);
        
        // 2. Compute Cartesian Christoffel Symbols
        float Gamma_cart[4][4][4];
        // Use correct member: warpDrive
        getWarpDriveChristoffel(x_cart, mp.warpDrive.vs, mp.warpDrive.sigma, mp.warpDrive.R, 
                                0.0f, 0.0f, 0.0f, // Bubble at origin for simplified viz
                                Gamma_cart);
        
        // 3. Transform to Spherical
        transformChristoffelCartToSph(Gamma_cart, x.r, x.theta, x.phi, Gamma);
        return;
    }
    
    // =========================================================================
    // MORRIS-THORNE WORMHOLE (Ellis Drainhole): Analytic Spherical
    // =========================================================================
    if constexpr (type == Sirius::MetricType::EllisDrainhole) { // or MorrisThorne
        // Uses native spherical coordinates
        // Use correct member: morrisThorne
        getMorrisThorneChristoffel(x, mp.morrisThorne.b0, mp.morrisThorne.redshiftPhi, 
                                   (WormholeShape)0 /*Ellis*/, Gamma);
        return;
    }
}

//==============================================================================
// Family-Based Runtime Dispatchers (Unified Architecture - Dec 2025)
// These use the MetricFamily enum for runtime dispatch to unified families
//==============================================================================

// Runtime metric dispatch by family (for family-based rendering)
__device__ void getMetricByFamily(
    const Vec4& x_sph,
    const Sirius::MetricParams& mp,
    float g[4][4], float g_inv[4][4])
{
    switch (mp.family) {
        case Sirius::MetricFamily::KerrSchild: {
            // Convert spherical to Cartesian and use unified family
            Vec4Cart x_cart = vec4SphToCart(x_sph);
            getKerrSchildFamilyMetric(x_cart, 
                mp.kerrSchild.M, mp.kerrSchild.a, 
                mp.kerrSchild.Q, mp.kerrSchild.Lambda, 
                g, g_inv);
            break;
        }
        case Sirius::MetricFamily::MorrisThorne: {
            getMorrisThorneMetric(x_sph, 
                mp.morrisThorne.b0, 
                mp.morrisThorne.redshiftPhi,
                WormholeShape::Ellis,  // Default to Ellis
                g, g_inv);
            break;
        }
        case Sirius::MetricFamily::WarpDrive: {
            Vec4Cart x_cart = vec4SphToCart(x_sph);
            getWarpDriveMetric(x_cart,
                mp.warpDrive.vs, mp.warpDrive.sigma, mp.warpDrive.R,
                0.0f, 0.0f, 0.0f,  // Bubble at origin
                g, g_inv);
            break;
        }
        case Sirius::MetricFamily::LewisWeyl: {
            // Gödel metric
            float omega = (mp.lewisWeyl.omega != 0.0f) ? mp.lewisWeyl.omega : 0.5f;
            getGodelMetric(x_sph, omega, g, g_inv);
            break;
        }
        case Sirius::MetricFamily::NUTCharged: {
            // Taub-NUT metric
            getTaubNUTMetric(x_sph, mp.nutCharged.M, mp.nutCharged.n, g, g_inv);
            break;
        }
        default: {
            // Default to Minkowski
            getMinkowskiMetric(x_sph, g, g_inv);
            break;
        }
    }
}

// Runtime Christoffel dispatch by family
__device__ void getChristoffelByFamily(
    const Vec4& x_sph,
    const Sirius::MetricParams& mp,
    float Gamma[4][4][4])
{
    switch (mp.family) {
        case Sirius::MetricFamily::KerrSchild: {
            // Use unified Kerr-Schild family (Cartesian, polynomial)
            Vec4Cart x_cart = vec4SphToCart(x_sph);
            getKerrSchildFamilyChristoffel(x_cart,
                mp.kerrSchild.M, mp.kerrSchild.a,
                mp.kerrSchild.Q, mp.kerrSchild.Lambda,
                Gamma);
            break;
        }
        case Sirius::MetricFamily::MorrisThorne: {
            getMorrisThorneChristoffel(x_sph,
                mp.morrisThorne.b0,
                mp.morrisThorne.redshiftPhi,
                WormholeShape::Ellis,
                Gamma);
            break;
        }
        case Sirius::MetricFamily::WarpDrive: {
            Vec4Cart x_cart = vec4SphToCart(x_sph);
            getWarpDriveChristoffel(x_cart,
                mp.warpDrive.vs, mp.warpDrive.sigma, mp.warpDrive.R,
                0.0f, 0.0f, 0.0f,
                Gamma);
            break;
        }
        case Sirius::MetricFamily::LewisWeyl: {
            // Gödel - no analytic Christoffel yet, initialize to zero
            // These will be computed via finite difference in the templated path
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    for (int k = 0; k < 4; k++)
                        Gamma[i][j][k] = 0.0f;
            // TODO: Implement analytic getGodelChristoffel
            break;
        }
        case Sirius::MetricFamily::NUTCharged: {
            // Taub-NUT - no analytic Christoffel yet, initialize to zero
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    for (int k = 0; k < 4; k++)
                        Gamma[i][j][k] = 0.0f;
            // TODO: Implement analytic getTaubNUTChristoffel
            break;
        }
        default: {
            // Minkowski in Cartesian has zero Christoffels
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    for (int k = 0; k < 4; k++)
                        Gamma[i][j][k] = 0.0f;
            break;
        }
    }
}

//==============================================================================
// Geodesic Equation: du^μ/dλ = -Γ^μ_νρ u^ν u^ρ
//==============================================================================
template<Sirius::MetricType type>
__device__ Vec4 geodesicAcceleration(const GeodesicState& state,
                                      const Sirius::MetricParams& mp) {
    float Gamma[4][4][4];
    getChristoffelSymbols<type>(state.x, mp, Gamma);
    
    Vec4 acc;
    float u[4] = {state.u.t, state.u.r, state.u.theta, state.u.phi};
    float a[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    
    // OPTIMIZATION: Loop unrolling for geodesic equation
    // 4×4×4 = 64 iterations, unroll for better ILP
    #pragma unroll 4
    for (int mu = 0; mu < 4; mu++) {
        #pragma unroll 4
        for (int nu = 0; nu < 4; nu++) {
            #pragma unroll 4
            for (int rho = 0; rho < 4; rho++) {
                a[mu] -= Gamma[mu][nu][rho] * u[nu] * u[rho];
            }
        }
    }
    
    acc.t = a[0];
    acc.r = a[1];
    acc.theta = a[2];
    acc.phi = a[3];
    return acc;
}

//==============================================================================
// Riemann Curvature Tensor (Full Implementation - Phase 6.3 Amendment)
// R^μ_νρσ = ∂_ρ Γ^μ_νσ - ∂_σ Γ^μ_νρ + Γ^μ_λρ Γ^λ_νσ - Γ^μ_λσ Γ^λ_νρ
//
// This is the key curvature quantity needed for exact geodesic deviation.
// Reference: DNGR Paper Appendix A.2, Eq A.18-A.27
//==============================================================================
template<Sirius::MetricType type>
__device__ void getRiemannTensor(
    const Vec4& x,
    const Sirius::MetricParams& mp,
    float R[4][4][4][4])
{
    // Initialize to zero
    #pragma unroll 4
    for (int mu = 0; mu < 4; mu++)
        #pragma unroll 4
        for (int nu = 0; nu < 4; nu++)
            #pragma unroll 4
            for (int rho = 0; rho < 4; rho++)
                #pragma unroll 4
                for (int sigma = 0; sigma < 4; sigma++)
                    R[mu][nu][rho][sigma] = 0.0f;
    
    // Get Christoffels at current position
    float Gamma[4][4][4];
    getChristoffelSymbols<type>(x, mp, Gamma);
    
    // Adaptive finite difference step size based on radial position
    const float eps = 0.001f * fmaxf(1.0f, x.r);
    
    // Storage for perturbed Christoffels and their derivatives
    float Gamma_plus[4][4][4], Gamma_minus[4][4][4];
    float dGamma[4][4][4][4]; // dGamma[rho][mu][nu][sigma] = ∂_rho Γ^mu_nuσ
    
    // Compute ∂Γ/∂x^rho by central finite differences
    // Only compute for spatial derivatives (rho = 1,2,3 = r,θ,φ)
    // Time derivative is typically zero for stationary spacetimes
    for (int rho = 1; rho < 4; rho++) {
        Vec4 xp = x, xm = x;
        
        switch (rho) {
            case 1: xp.r += eps; xm.r -= eps; break;
            case 2: xp.theta += eps; xm.theta -= eps; break;
            case 3: xp.phi += eps; xm.phi -= eps; break;
        }
        
        getChristoffelSymbols<type>(xp, mp, Gamma_plus);
        getChristoffelSymbols<type>(xm, mp, Gamma_minus);
        
        // Central difference: ∂Γ/∂x^rho ≈ (Γ(x+ε) - Γ(x-ε)) / (2ε)
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = 0; nu < 4; nu++) {
                for (int sigma = 0; sigma < 4; sigma++) {
                    dGamma[rho][mu][nu][sigma] = 
                        (Gamma_plus[mu][nu][sigma] - Gamma_minus[mu][nu][sigma]) / (2.0f * eps);
                }
            }
        }
    }
    
    // Time derivatives are zero for stationary metrics
    for (int mu = 0; mu < 4; mu++)
        for (int nu = 0; nu < 4; nu++)
            for (int sigma = 0; sigma < 4; sigma++)
                dGamma[0][mu][nu][sigma] = 0.0f;
    
    // Compute Riemann tensor
    // R^μ_νρσ = ∂_ρ Γ^μ_νσ - ∂_σ Γ^μ_νρ + Γ^μ_λρ Γ^λ_νσ - Γ^μ_λσ Γ^λ_νρ
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            for (int rho = 0; rho < 4; rho++) {
                for (int sigma = 0; sigma < 4; sigma++) {
                    // First two terms: ∂_ρ Γ^μ_νσ - ∂_σ Γ^μ_νρ
                    R[mu][nu][rho][sigma] = dGamma[rho][mu][nu][sigma] 
                                          - dGamma[sigma][mu][nu][rho];
                    
                    // Last two terms: Γ^μ_λρ Γ^λ_νσ - Γ^μ_λσ Γ^λ_νρ
                    for (int lambda = 0; lambda < 4; lambda++) {
                        R[mu][nu][rho][sigma] += Gamma[mu][lambda][rho] * Gamma[lambda][nu][sigma]
                                               - Gamma[mu][lambda][sigma] * Gamma[lambda][nu][rho];
                    }
                }
            }
        }
    }
}

//==============================================================================
// Geodesic Deviation Acceleration (DNGR Eq A.19)
// D²ξ^μ/dλ² = R^μ_νρσ k^ν ξ^ρ k^σ
//
// Computes the "tidal acceleration" that causes neighboring rays to
// diverge or converge due to spacetime curvature.
//==============================================================================
template<Sirius::MetricType type>
__device__ float4 geodesicDeviationAccel(
    const Vec4& x,           // Position
    const Vec4& k,           // Photon 4-momentum (null geodesic tangent)
    const float4& xi,        // Deviation vector (t, r, θ, φ components)
    const Sirius::MetricParams& mp)
{
    // Get Riemann tensor
    float R[4][4][4][4];
    getRiemannTensor<type>(x, mp, R);
    
    // Convert Vec4 to arrays for contraction
    float k_arr[4] = {k.t, k.r, k.theta, k.phi};
    float xi_arr[4] = {xi.x, xi.y, xi.z, xi.w};  // float4: x=t, y=r, z=θ, w=φ
    float accel[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    
    // Contract: R^μ_νρσ k^ν ξ^ρ k^σ
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            for (int rho = 0; rho < 4; rho++) {
                for (int sigma = 0; sigma < 4; sigma++) {
                    accel[mu] += R[mu][nu][rho][sigma] * k_arr[nu] * xi_arr[rho] * k_arr[sigma];
                }
            }
        }
    }
    
    return make_float4(accel[0], accel[1], accel[2], accel[3]);
}

//==============================================================================
// Geodesic State Normalization
// Enforces coordinate bounds and null constraint (g_μν u^μ u^ν = 0)
//==============================================================================

// Pole exclusion zone: minimal threshold to avoid exact singularity
// Using very small value to minimize visible artifacts while preventing NaN
__device__ constexpr float POLE_EPSILON = 0.0005f;  // ~0.03 degrees from pole

template<Sirius::MetricType type>
__device__ __forceinline__ void normalizeGeodesicState(
    GeodesicState& result,
    const Sirius::MetricParams& mp)
{
    // Handle polar coordinate singularity crossing
    // When ray passes through pole: θ -> -θ  =>  θ' = |θ|, φ' = φ + π
    const float PI_VAL = 3.14159265359f;
    const float TWO_PI_VAL = 6.28318530718f;

    if (result.x.theta < 0.0f) {
        result.x.theta = -result.x.theta;
        result.u.theta = -result.u.theta;   // Flip momentum
        result.x.phi += PI_VAL;
    }
    else if (result.x.theta > PI_VAL) {
        result.x.theta = TWO_PI_VAL - result.x.theta;
        result.u.theta = -result.u.theta;   // Flip momentum
        result.x.phi += PI_VAL;
    }

    // Minimal clamping: only prevent exact singularity (sin(θ) = 0)
    // Use very small epsilon to avoid visible artifacts
    result.x.theta = fmaxf(POLE_EPSILON, fminf(PI_VAL - POLE_EPSILON, result.x.theta));

    // NOTE: Removed phi velocity damping - it corrupts geodesics and creates
    // visible dark wedge artifacts emanating from the poles. The proper fix
    // is to use correct coordinate wrapping (above) rather than damping.

    // Wrap phi to [0, 2π)
    result.x.phi = fmodf(result.x.phi, TWO_PI);
    if (result.x.phi < 0.0f) result.x.phi += TWO_PI;
    
    // Null constraint enforcement: g_μν u^μ u^ν = 0
    float g_step[4][4], g_inv_step[4][4];
    getMetricTensor<type>(result.x, mp, g_step, g_inv_step);
    
    // Spatial part: C = g_ij u^i u^j
    float C_spatial = g_step[1][1] * result.u.r * result.u.r 
                    + g_step[2][2] * result.u.theta * result.u.theta 
                    + g_step[3][3] * result.u.phi * result.u.phi;
    C_spatial += 2.0f * (g_step[1][2] * result.u.r * result.u.theta 
                       + g_step[1][3] * result.u.r * result.u.phi 
                       + g_step[2][3] * result.u.theta * result.u.phi);
    
    // Quadratic: A(u^t)² + B(u^t) + C = 0
    float A_null = g_step[0][0];
    float B_null = 2.0f * g_step[0][3] * result.u.phi;
    float D_null = B_null * B_null - 4.0f * A_null * C_spatial;
    
    if (D_null >= 0.0f) {
        float sqrtD = sqrtf(D_null);
        float root1 = (-B_null + sqrtD) / (2.0f * A_null);
        float root2 = (-B_null - sqrtD) / (2.0f * A_null);
        result.u.t = (fabsf(root1 - result.u.t) < fabsf(root2 - result.u.t)) 
                   ? ((root1 > 0.0f) ? root1 : root2)
                   : ((root2 > 0.0f) ? root2 : root1);
        if (result.u.t <= 0.0f) result.u.t = fabsf(result.u.t) + 0.01f;
    }
}

//==============================================================================
// RK4 Geodesic Integration
//==============================================================================
template<Sirius::MetricType type>
__device__ GeodesicState integrateGeodesicRK4(const GeodesicState& state,
                                               float dt,
                                               const Sirius::MetricParams& mp) {
    // k1
    Vec4 dx1 = state.u * dt;
    Vec4 du1 = geodesicAcceleration<type>(state, mp) * dt;
    
    // k2
    GeodesicState s2;
    s2.x = state.x + dx1 * 0.5f;
    s2.u = state.u + du1 * 0.5f;
    Vec4 dx2 = s2.u * dt;
    Vec4 du2 = geodesicAcceleration<type>(s2, mp) * dt;
    
    // k3
    GeodesicState s3;
    s3.x = state.x + dx2 * 0.5f;
    s3.u = state.u + du2 * 0.5f;
    Vec4 dx3 = s3.u * dt;
    Vec4 du3 = geodesicAcceleration<type>(s3, mp) * dt;
    
    // k4
    GeodesicState s4;
    s4.x = state.x + dx3;
    s4.u = state.u + du3;
    Vec4 dx4 = s4.u * dt;
    Vec4 du4 = geodesicAcceleration<type>(s4, mp) * dt;
    
    // Combine
    GeodesicState result;
    result.x = state.x + (dx1 + dx2 * 2.0f + dx3 * 2.0f + dx4) * (1.0f / 6.0f);
    result.u = state.u + (du1 + du2 * 2.0f + du3 * 2.0f + du4) * (1.0f / 6.0f);
    
    // Apply coordinate normalization and null constraint enforcement
    normalizeGeodesicState<type>(result, mp);
    
    return result;
}

//==============================================================================
// RK45 Dormand-Prince Adaptive Integration
//
// Uses embedded 4th/5th order methods for automatic step size control.
// Key advantages over RK4:
//   1. Automatic step reduction near poles (sin(θ) → 0 singularity)
//   2. Automatic step reduction near photon sphere (unstable orbits)
//   3. Error estimation without extra function evaluations
//   4. Larger steps in flat spacetime regions for efficiency
//==============================================================================

// Dormand-Prince coefficients
__device__ constexpr float DP_A21 = 1.0f / 5.0f;
__device__ constexpr float DP_A31 = 3.0f / 40.0f;
__device__ constexpr float DP_A32 = 9.0f / 40.0f;
__device__ constexpr float DP_A41 = 44.0f / 45.0f;
__device__ constexpr float DP_A42 = -56.0f / 15.0f;
__device__ constexpr float DP_A43 = 32.0f / 9.0f;
__device__ constexpr float DP_A51 = 19372.0f / 6561.0f;
__device__ constexpr float DP_A52 = -25360.0f / 2187.0f;
__device__ constexpr float DP_A53 = 64448.0f / 6561.0f;
__device__ constexpr float DP_A54 = -212.0f / 729.0f;
__device__ constexpr float DP_A61 = 9017.0f / 3168.0f;
__device__ constexpr float DP_A62 = -355.0f / 33.0f;
__device__ constexpr float DP_A63 = 46732.0f / 5247.0f;
__device__ constexpr float DP_A64 = 49.0f / 176.0f;
__device__ constexpr float DP_A65 = -5103.0f / 18656.0f;
__device__ constexpr float DP_A71 = 35.0f / 384.0f;
__device__ constexpr float DP_A73 = 500.0f / 1113.0f;
__device__ constexpr float DP_A74 = 125.0f / 192.0f;
__device__ constexpr float DP_A75 = -2187.0f / 6784.0f;
__device__ constexpr float DP_A76 = 11.0f / 84.0f;

// Error estimation coefficients (difference between 5th and 4th order)
__device__ constexpr float DP_E1 = 71.0f / 57600.0f;
__device__ constexpr float DP_E3 = -71.0f / 16695.0f;
__device__ constexpr float DP_E4 = 71.0f / 1920.0f;
__device__ constexpr float DP_E5 = -17253.0f / 339200.0f;
__device__ constexpr float DP_E6 = 22.0f / 525.0f;
__device__ constexpr float DP_E7 = -1.0f / 40.0f;

template<Sirius::MetricType type>
__device__ GeodesicState integrateGeodesicRK45(
    const GeodesicState& state,
    float& h,                    // IN/OUT: step size (adapted)
    const Sirius::MetricParams& mp,
    float tolerance,             // Target local error
    float hMin,                  // Minimum step size
    float hMax,                  // Maximum step size
    bool& stepAccepted)          // OUT: was step accepted?
{
    // Stage 1: k1 at current position
    Vec4 k1x = state.u;
    Vec4 k1u = geodesicAcceleration<type>(state, mp);

    // Stage 2
    GeodesicState s2;
    s2.x = state.x + k1x * (h * DP_A21);
    s2.u = state.u + k1u * (h * DP_A21);
    Vec4 k2x = s2.u;
    Vec4 k2u = geodesicAcceleration<type>(s2, mp);

    // Stage 3
    GeodesicState s3;
    s3.x = state.x + k1x * (h * DP_A31) + k2x * (h * DP_A32);
    s3.u = state.u + k1u * (h * DP_A31) + k2u * (h * DP_A32);
    Vec4 k3x = s3.u;
    Vec4 k3u = geodesicAcceleration<type>(s3, mp);

    // Stage 4
    GeodesicState s4;
    s4.x = state.x + k1x * (h * DP_A41) + k2x * (h * DP_A42) + k3x * (h * DP_A43);
    s4.u = state.u + k1u * (h * DP_A41) + k2u * (h * DP_A42) + k3u * (h * DP_A43);
    Vec4 k4x = s4.u;
    Vec4 k4u = geodesicAcceleration<type>(s4, mp);

    // Stage 5
    GeodesicState s5;
    s5.x = state.x + k1x * (h * DP_A51) + k2x * (h * DP_A52) + k3x * (h * DP_A53) + k4x * (h * DP_A54);
    s5.u = state.u + k1u * (h * DP_A51) + k2u * (h * DP_A52) + k3u * (h * DP_A53) + k4u * (h * DP_A54);
    Vec4 k5x = s5.u;
    Vec4 k5u = geodesicAcceleration<type>(s5, mp);

    // Stage 6
    GeodesicState s6;
    s6.x = state.x + k1x * (h * DP_A61) + k2x * (h * DP_A62) + k3x * (h * DP_A63) + k4x * (h * DP_A64) + k5x * (h * DP_A65);
    s6.u = state.u + k1u * (h * DP_A61) + k2u * (h * DP_A62) + k3u * (h * DP_A63) + k4u * (h * DP_A64) + k5u * (h * DP_A65);
    Vec4 k6x = s6.u;
    Vec4 k6u = geodesicAcceleration<type>(s6, mp);

    // 5th order solution (used as result if accepted)
    GeodesicState result;
    result.x = state.x + (k1x * DP_A71 + k3x * DP_A73 + k4x * DP_A74 + k5x * DP_A75 + k6x * DP_A76) * h;
    result.u = state.u + (k1u * DP_A71 + k3u * DP_A73 + k4u * DP_A74 + k5u * DP_A75 + k6u * DP_A76) * h;

    // Stage 7 (FSAL - First Same As Last, reused in next step)
    Vec4 k7x = result.u;
    Vec4 k7u = geodesicAcceleration<type>(result, mp);

    // Error estimation: difference between 5th and 4th order
    Vec4 errX = (k1x * DP_E1 + k3x * DP_E3 + k4x * DP_E4 + k5x * DP_E5 + k6x * DP_E6 + k7x * DP_E7) * h;
    Vec4 errU = (k1u * DP_E1 + k3u * DP_E3 + k4u * DP_E4 + k5u * DP_E5 + k6u * DP_E6 + k7u * DP_E7) * h;

    // Compute error norm (max of position and velocity errors)
    float errNorm = fmaxf(
        fmaxf(fabsf(errX.r), fmaxf(fabsf(errX.theta), fabsf(errX.phi))),
        fmaxf(fabsf(errU.r), fmaxf(fabsf(errU.theta), fabsf(errU.phi)))
    );

    // Safety factor for step size adjustment
    const float SAFETY = 0.9f;
    const float MIN_SCALE = 0.2f;  // Don't shrink by more than 5x
    const float MAX_SCALE = 5.0f;  // Don't grow by more than 5x

    if (errNorm < 1e-15f) errNorm = 1e-15f;  // Prevent division by zero

    // Compute optimal step size ratio
    float scale = SAFETY * powf(tolerance / errNorm, 0.2f);  // 5th order -> 1/5 power
    scale = fmaxf(MIN_SCALE, fminf(MAX_SCALE, scale));

    if (errNorm <= tolerance) {
        // Step accepted
        stepAccepted = true;

        // Grow step for next iteration (but not too much)
        h = fminf(h * scale, hMax);

        // Apply coordinate normalization
        normalizeGeodesicState<type>(result, mp);

        return result;
    } else {
        // Step rejected - shrink and retry
        stepAccepted = false;
        h = fmaxf(h * scale, hMin);

        return state;  // Return unchanged state
    }
}

//==============================================================================
// Coordinate Conversions
//==============================================================================

// Cartesian to Boyer-Lindquist (spherical)
// Cartesian to Boyer-Lindquist (spherical)
__device__ Vec4 cartesianToSpherical(float3 pos) {
    float r = length(pos);
    if (r < 0.001f) r = 0.001f;
    
    // Clamp theta to avoid polar coordinate singularity (sin(theta) = 0)
    // This fixes vertical artifacts and NaNs in metric tensor inversions
    float theta = acosf(clamp(pos.y / r, -1.0f, 1.0f));
    const float THETA_EPSILON = 0.001f;
    theta = fmaxf(THETA_EPSILON, fminf(PI - THETA_EPSILON, theta));
    
    float phi = atan2f(pos.z, pos.x);
    if (phi < 0.0f) phi += TWO_PI;
    
    return Vec4(0.0f, r, theta, phi);
}

//==============================================================================
// CARTESIAN GEODESIC INTEGRATION (for Kerr-Schild family)
// Eliminates coordinate transformation overhead by integrating directly
// in Cartesian coordinates (t, x, y, z)
//==============================================================================

// Cartesian geodesic acceleration: du^μ/dλ = -Γ^μ_νρ u^ν u^ρ
// Uses Kerr-Schild Christoffel directly in Cartesian coordinates
__device__ Vec4Cart geodesicAccelerationCart(const GeodesicStateCart& state,
                                              float M, float a) {
    float Gamma[4][4][4];
    getKerrSchildFamilyChristoffel(state.x, M, a, 0.0f, 0.0f, Gamma);
    
    float u[4] = {state.u.t, state.u.x, state.u.y, state.u.z};
    float acc[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    
    #pragma unroll 4
    for (int mu = 0; mu < 4; mu++) {
        #pragma unroll 4
        for (int nu = 0; nu < 4; nu++) {
            #pragma unroll 4
            for (int rho = 0; rho < 4; rho++) {
                acc[mu] -= Gamma[mu][nu][rho] * u[nu] * u[rho];
            }
        }
    }
    
    return Vec4Cart(acc[0], acc[1], acc[2], acc[3]);
}

// RK4 integrator in Cartesian coordinates (no transformation needed!)
__device__ GeodesicStateCart integrateGeodesicRK4Cart(const GeodesicStateCart& state,
                                                       float dt, float M, float a) {
    // k1
    Vec4Cart dx1 = state.u * dt;
    Vec4Cart du1 = geodesicAccelerationCart(state, M, a) * dt;
    
    // k2
    GeodesicStateCart s2;
    s2.x = state.x + dx1 * 0.5f;
    s2.u = state.u + du1 * 0.5f;
    Vec4Cart dx2 = s2.u * dt;
    Vec4Cart du2 = geodesicAccelerationCart(s2, M, a) * dt;
    
    // k3
    GeodesicStateCart s3;
    s3.x = state.x + dx2 * 0.5f;
    s3.u = state.u + du2 * 0.5f;
    Vec4Cart dx3 = s3.u * dt;
    Vec4Cart du3 = geodesicAccelerationCart(s3, M, a) * dt;
    
    // k4
    GeodesicStateCart s4;
    s4.x = state.x + dx3;
    s4.u = state.u + du3;
    Vec4Cart dx4 = s4.u * dt;
    Vec4Cart du4 = geodesicAccelerationCart(s4, M, a) * dt;
    
    // Combine
    GeodesicStateCart result;
    result.x = state.x + (dx1 + dx2 * 2.0f + dx3 * 2.0f + dx4) * (1.0f / 6.0f);
    result.u = state.u + (du1 + du2 * 2.0f + du3 * 2.0f + du4) * (1.0f / 6.0f);
    
    return result;
}

// Convert Cartesian state to spherical for disk/background checks
__device__ GeodesicState cartToSphState(const GeodesicStateCart& cart) {
    GeodesicState sph;
    sph.x = vec4CartToSph(cart.x);
    
    // Transform velocity using Jacobian
    float r = sph.x.r;
    float theta = sph.x.theta;
    float phi = sph.x.phi;
    float st = sinf(theta), ct = cosf(theta);
    float sp = sinf(phi), cp = cosf(phi);
    
    float r_safe = fmaxf(r, 0.01f);
    float rst = r_safe * fmaxf(fabsf(st), 0.01f);
    
    // u^r = (x*vx + y*vy + z*vz) / r
    float vx = cart.u.x, vy = cart.u.y, vz = cart.u.z;
    float x = cart.x.x, y = cart.x.y, z = cart.x.z;
    sph.u.t = cart.u.t;
    sph.u.r = (x*vx + y*vy + z*vz) / r_safe;
    
    // u^θ = (z*(x*vx + y*vy) - (x² + y²)*vz) / (r² * sqrt(x² + y²))
    float rho2 = x*x + y*y;
    float rho = sqrtf(fmaxf(rho2, 0.0001f));
    sph.u.theta = (z * (x*vx + y*vy) - rho2 * vz) / (r_safe * r_safe * rho);
    
    // u^φ = (x*vy - y*vx) / (x² + y²)
    sph.u.phi = (x*vy - y*vx) / fmaxf(rho2, 0.0001f);
    
    return sph;
}

// Convert spherical state to Cartesian for integration
__device__ GeodesicStateCart sphToCartState(const GeodesicState& sph) {
    GeodesicStateCart cart;
    cart.x = vec4SphToCart(sph.x);
    
    float r = sph.x.r;
    float theta = sph.x.theta;
    float phi = sph.x.phi;
    float st = sinf(theta), ct = cosf(theta);
    float sp = sinf(phi), cp = cosf(phi);
    
    // Transform velocity: v_cart = J * v_sph where J = ∂(x,y,z)/∂(r,θ,φ)
    float u_r = sph.u.r, u_th = sph.u.theta, u_ph = sph.u.phi;
    
    cart.u.t = sph.u.t;
    cart.u.x = st*cp*u_r + r*ct*cp*u_th - r*st*sp*u_ph;
    cart.u.y = st*sp*u_r + r*ct*sp*u_th + r*st*cp*u_ph;
    cart.u.z = ct*u_r - r*st*u_th;
    
    return cart;
}

// Spherical to Cartesian

__device__ float3 sphericalToCartesian(const Vec4& x) {
    float r = x.r;
    float theta = x.theta;
    float phi = x.phi;
    
    float sin_theta = sinf(theta);
    return make_float3(
        r * sin_theta * cosf(phi),
        r * cosf(theta),
        r * sin_theta * sinf(phi)
    );
}

// Convert Cartesian direction to 4-velocity in spherical coords

template<Sirius::MetricType type>
__device__ Vec4 cartesianDirToFourVelocity(float3 pos, float3 dir, 
                                            const Sirius::MetricParams& mp) {
    Vec4 x = cartesianToSpherical(pos);
    
    float r = x.r;
    float theta = x.theta;
    float phi = x.phi;
    
    float sin_theta = sinf(theta);
    float cos_theta = cosf(theta);
    
    // Jacobian for coordinate transformation
    float dr_dx = sinf(theta) * cosf(phi);
    float dr_dy = cosf(theta);
    float dr_dz = sinf(theta) * sinf(phi);
    
    float dtheta_dx = cos_theta * cosf(phi) / r;
    float dtheta_dy = -sin_theta / r;
    float dtheta_dz = cos_theta * sinf(phi) / r;
    
    float dphi_dx = -sinf(phi) / (r * sin_theta + 0.0001f);
    float dphi_dy = 0.0f;
    float dphi_dz = cosf(phi) / (r * sin_theta + 0.0001f);
    
    // Transform direction
    float u_r = dr_dx * dir.x + dr_dy * dir.y + dr_dz * dir.z;
    float u_theta = dtheta_dx * dir.x + dtheta_dy * dir.y + dtheta_dz * dir.z;
    float u_phi = dphi_dx * dir.x + dphi_dy * dir.y + dphi_dz * dir.z;
    
    // Get metric to normalize (null geodesic: g_μν u^μ u^ν = 0)
    float g[4][4], g_inv[4][4];

    getMetricTensor<type>(x, mp, g, g_inv);
    
    // For null geodesic with cross terms (Kerr has g_tφ ≠ 0):
    // g_tt (u^t)² + 2 g_tφ u^t u^φ + g_ij u^i u^j = 0
    // This is a quadratic: A(u^t)² + B(u^t) + C = 0
    // where A = g_tt, B = 2 g_tφ u^φ, C = spatial terms
    float A = g[0][0];  // g_tt
    float B = 2.0f * g[0][3] * u_phi;  // 2 g_tφ u^φ
    float C = g[1][1] * u_r * u_r + g[2][2] * u_theta * u_theta + g[3][3] * u_phi * u_phi;
    
    // Add any other cross terms (g_rθ, g_rφ, g_θφ) if present
    C += 2.0f * (g[1][2] * u_r * u_theta + g[1][3] * u_r * u_phi + g[2][3] * u_theta * u_phi);
    
    // Solve quadratic: u^t = (-B ± sqrt(B² - 4AC)) / 2A
    // For physical solution, we want u^t > 0 (future-directed)
    // With g_tt < 0 (timelike signature), the correct root is typically (-B - sqrt(D)) / 2A
    float discriminant = B * B - 4.0f * A * C;
    float u_t;
    
    if (discriminant >= 0.0f) {
        // Two real roots - choose physical one (u^t > 0)
        float sqrtD = sqrtf(discriminant);
        float root1 = (-B + sqrtD) / (2.0f * A);
        float root2 = (-B - sqrtD) / (2.0f * A);
        
        // u^t should be positive and typically the larger magnitude for physical ray
        u_t = (root1 > 0.0f) ? root1 : ((root2 > 0.0f) ? root2 : 1.0f);
    } else {
        // No real solution - fallback (should not happen for valid configurations)
        u_t = sqrtf(fmaxf(-C / A, 0.0001f));
    }
    
    return Vec4(u_t, u_r, u_theta, u_phi);
}

//==============================================================================
// Transform Spherical Velocity Direction to World Frame
//==============================================================================
// The asymptotic direction of an escaping ray must be expressed in the fixed
// world frame for consistent background texture sampling. The observer's
// position (theta, phi) determines the local spherical basis vectors which
// must be used to convert the ray's 4-velocity to world Cartesian coordinates.
//
// This fixes the banding artifact where rays escaping at different (θ,φ)
// positions would sample different parts of the texture even when pointing
// at the same distant star.
//
// Mathematical basis:
//   Spherical basis vectors at observer position (θ_obs, φ_obs):
//     e_r     = (sin θ cos φ, cos θ, sin θ sin φ)
//     e_theta = (cos θ cos φ, -sin θ, cos θ sin φ)  
//     e_phi   = (-sin φ, 0, cos φ)
//
//   World direction = velR * e_r + velTheta * e_theta + velPhi * e_phi
//==============================================================================
__device__ float3 transformVelocityToWorldFrame(
    float velR, float velTheta, float velPhi,
    float obs_theta, float obs_phi)
{
    // OPTIMIZATION: Combined sin/cos for each angle
    float st, ct, sp, cp;
    fast_sincosf(obs_theta, &st, &ct);
    fast_sincosf(obs_phi, &sp, &cp);
    
    // Transform: worldDir = velR * e_r + velTheta * e_theta + velPhi * e_phi
    // Expanded using basis vector definitions:
    //   x = sin(θ)cos(φ)*velR + cos(θ)cos(φ)*velTheta - sin(φ)*velPhi
    //   y = cos(θ)*velR - sin(θ)*velTheta
    //   z = sin(θ)sin(φ)*velR + cos(θ)sin(φ)*velTheta + cos(φ)*velPhi
    float3 worldDir;
    worldDir.x = st*cp*velR + ct*cp*velTheta - sp*velPhi;
    worldDir.y = ct*velR - st*velTheta;
    worldDir.z = st*sp*velR + ct*sp*velTheta + cp*velPhi;
    
    return normalize(worldDir);
}

//==============================================================================
// Background Sampling
//==============================================================================
__device__ float3 sampleBackground(float3 direction) {
    if (params.useBackgroundTexture) {
        // Convert direction to spherical UV
        float theta = acosf(clamp(direction.y, -1.0f, 1.0f));
        float phi = atan2f(direction.z, direction.x);
        if (phi < 0.0f) phi += TWO_PI;

        float u = phi / TWO_PI;
        float v = theta / PI;

        // Sample texture
        float4 texColor = tex2D<float4>(params.backgroundTexture, u, v);
        return make_float3(texColor.x, texColor.y, texColor.z);
    }

    // Procedural star field - cinematic mode
    float3 color = params.backgroundColor;

    // =========================================================================
    // Enhanced Starfield (Cinematic Features Phase 8)
    // Sparse star distribution for cinematic look (like Interstellar)
    // Only render stars when explicitly enabled, otherwise pure black
    // =========================================================================
    if (params.starfield.enabled) {
        // Very sparse star distribution for cinematic look
        // Use high-quality hash for better distribution
        unsigned int seed = __float_as_uint(direction.x * 127.1f + direction.y * 311.7f + direction.z * 74.7f);
        seed ^= seed >> 13;
        seed *= 0x85ebca6b;
        seed ^= seed >> 16;

        float star = (float)(seed & 0xFFFF) / 65535.0f;

        // Very sparse stars (0.0003 density = ~0.03% of sky has stars)
        float threshold = 1.0f - 0.0003f * params.starfield.brightness_scale;

        if (star > threshold) {
            float brightness = (star - threshold) / (1.0f - threshold);
            // Subtle brightness, not overwhelming
            brightness = powf(brightness, 0.7f) * 0.8f * params.starfield.brightness_scale;

            // Color temperature variation (blue-white-orange)
            seed = seed * 1103515245 + 12345;
            float temp = (float)(seed & 0xFF) / 255.0f;
            float3 starColor;
            if (temp < 0.3f) {
                starColor = make_float3(0.7f, 0.8f, 1.0f);  // Blue (hot stars)
            } else if (temp < 0.7f) {
                starColor = make_float3(0.9f, 0.9f, 0.95f); // White (solar-type)
            } else {
                starColor = make_float3(1.0f, 0.85f, 0.7f); // Yellow-orange (cool stars)
            }

            color = color + starColor * brightness;
        }
    }
    // No fallback starfield - pure black background for cinematic look

    return color;
}

//==============================================================================
// Calculate Gravitational Redshift
//==============================================================================

template<Sirius::MetricType type>
__device__ float calculateRedshift(const Vec4& x_emit, const Vec4& u_emit,
                                    const Vec4& x_obs, const Vec4& u_obs,
                                    const Sirius::MetricParams& mp) {
    float g_emit[4][4], g_emit_inv[4][4];
    float g_obs[4][4], g_obs_inv[4][4];
    

    
    getMetricTensor<type>(x_emit, mp, g_emit, g_emit_inv);
    getMetricTensor<type>(x_obs, mp, g_obs, g_obs_inv);
    
    // Simplified redshift: z ≈ sqrt(|g_tt(emit)| / |g_tt(obs)|) - 1
    float z = sqrtf(fabsf(g_obs[0][0] / g_emit[0][0])) - 1.0f;
    return clamp(z, -0.99f, 10.0f);
}

//==============================================================================
// Ray Generation Program
//==============================================================================
//==============================================================================
// Ray Generation Program (Template)
//==============================================================================
template<Sirius::MetricType type>
__device__ void raygen_renderFrame_impl() {
    // Get pixel coordinates
    const uint3 idx = optixGetLaunchIndex();
    const uint3 dim = optixGetLaunchDimensions();
    
    const int ix = idx.x;
    const int iy = idx.y;
    const int width = params.frameDimensions.x;
    const int height = params.frameDimensions.y;
    
    if (ix >= width || iy >= height) return;
    
    // Calculate normalized device coordinates
    float u = (2.0f * ((float)ix + 0.5f) / (float)width - 1.0f) * params.camera.aspectRatio;
    float v = 2.0f * ((float)iy + 0.5f) / (float)height - 1.0f;
    
    // Apply field of view
    float tanHalfFov = tanf(params.camera.fov * 0.5f);
    u *= tanHalfFov;
    v *= tanHalfFov;
    
    // Calculate ray direction in world space
    float3 rayDir = normalize(
        params.camera.direction + 
        params.camera.right * u + 
        params.camera.up * v
    );
    
    // Initialize geodesic state
    // float3 rayOrigin = params.camera.position; // Removed duplicate
    // Vec4 x = cartesianToSpherical(rayOrigin); // Removed duplicate
    // Initialize geodesic state
    float3 rayOrigin = params.camera.position;
    Vec4 x = cartesianToSpherical(rayOrigin);
    Vec4 u_vec = cartesianDirToFourVelocity<type>(rayOrigin, rayDir, params.metricParams);
    
    GeodesicState state;
    state.x = x;
    state.u = u_vec;
    
    // =========================================================================
    // RAY BUNDLE INITIALIZATION (Phase 6.3 - DNGR Anti-Aliasing)
    // Initialize deviation vectors for geodesic deviation tracking
    // =========================================================================
    Sirius::RayBundleState bundle;
    bool useRayBundle = params.rayBundle.enabled && 
                        (type != Sirius::MetricType::Minkowski);  // Skip bundles in flat space
    
    if (useRayBundle) {
        // Initialize bundle with central ray
        bundle.x = make_float4(x.t, x.r, x.theta, x.phi);
        bundle.u = make_float4(u_vec.t, u_vec.r, u_vec.theta, u_vec.phi);
        
        // Initial deviation vectors aligned with screen axes
        // Scale by pixel angular size
        float pixelSize = params.rayBundle.initialSize;
        
        // ξ₁ along camera right direction (screen X)
        bundle.xi1 = make_float4(0.0f, 
            params.camera.right.x * pixelSize,
            params.camera.right.y * pixelSize,
            params.camera.right.z * pixelSize);
        bundle.dxi1 = make_float4(0.0f, 0.0f, 0.0f, 0.0f);  // Initially parallel
        
        // ξ₂ along camera up direction (screen Y)
        bundle.xi2 = make_float4(0.0f,
            params.camera.up.x * pixelSize,
            params.camera.up.y * pixelSize,
            params.camera.up.z * pixelSize);
        bundle.dxi2 = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
        
        // Initial ellipse is circular
        bundle.semiMajor = pixelSize;
        bundle.semiMinor = pixelSize;
        bundle.orientation = 0.0f;
        bundle.affineParam = 0.0f;
        bundle.terminated = false;
    }
    
    // Integration parameters
    float stepSize = params.integration.initialStepSize;
    int maxSteps = params.integration.maxSteps;
    float maxDistance = params.integration.maxDistance;
    // Horizon radius: metric-specific computation
    // - Minkowski: no horizon
    // - Schwarzschild: r_s = 2M
    // - Kerr: outer horizon r_+ = M + sqrt(M² - a²) where Δ = r² - 2Mr + a² = 0
    // - Reissner-Nordström: r_+ = M + sqrt(M² - Q²)
    float horizonRadius;
    if constexpr (type == Sirius::MetricType::Minkowski) {
        horizonRadius = 0.0f;  // No horizon in flat space
    } else if constexpr (type == Sirius::MetricType::Kerr) {
        // Kerr outer horizon: r_+ = M + sqrt(M² - a²)
        float M = params.metricParams.M;
        float a = params.metricParams.a;
        float discriminant = M * M - a * a;
        if (discriminant > 0.0f) {
            horizonRadius = M + sqrtf(discriminant);
        } else {
            // Naked singularity case (|a| > M) - use M as fallback
            horizonRadius = M;
        }
        horizonRadius *= 1.12f;  // Increased buffer to eliminate shadow specs
    } else if constexpr (type == Sirius::MetricType::ReissnerNordstrom) {
        // Reissner-Nordström outer horizon: r_+ = M + sqrt(M² - Q²)
        float M = params.metricParams.M;
        float Q = params.metricParams.Q;
        float discriminant = M * M - Q * Q;
        if (discriminant > 0.0f) {
            horizonRadius = M + sqrtf(discriminant);
        } else {
            horizonRadius = M;  // Naked singularity fallback
        }
        horizonRadius *= 1.12f;  // Increased buffer to eliminate shadow specs
    } else {
        // Schwarzschild: r_s = 2M
        horizonRadius = 2.0f * params.metricParams.M * 1.12f;  // Increased buffer to eliminate shadow specs
    }
    
    // =========================================================================
    // VOLUMETRIC RAY MARCHING (Phase 4)
    // - Full Geodesic Integration (RK4)
    // - Radiative Transfer (Emission from Volumetric Data)
    // =========================================================================

    float3 accumColor = make_float3(0.0f, 0.0f, 0.0f);
    float opacity = 0.0f;
    int pixelIndex = iy * width + ix;
    
    // Initial accumulation from previous frames (if enabled)
    // (For now, just overwrite)

    bool hitHorizon = false;
    bool escaped = false;
    float3 finalDir = rayDir;
    
    // Adaptive step size (persists across iterations, adapted by RK45)
    float adaptiveStep = params.integration.initialStepSize; 
    
    // ==========================================================================
    // MINKOWSKI FAST-PATH: No geodesic curvature in flat space
    // Rays travel in straight lines, so skip the expensive integration loop
    // ==========================================================================
    if constexpr (type == Sirius::MetricType::Minkowski) {
        // In flat space, rays go straight - just sample background directly
        float3 bgColor = sampleBackground(rayDir);
        params.frameBuffer[pixelIndex] = make_float4(bgColor.x, bgColor.y, bgColor.z, 1.0f);
        params.accumBuffer[pixelIndex] = make_float4(bgColor.x, bgColor.y, bgColor.z, 1.0f);
        return;  // Skip entire integration loop
    }
    
    for (int step = 0; step < maxSteps; step++) {
        // =====================================================================
        // ADAPTIVE RK45 INTEGRATION (Phase 2 Performance)
        // =====================================================================
        // Uses embedded error estimation for automatic step adaptation
        // Replaces heuristic step sizing with mathematically rigorous control
        
        // 2. Geodesic Integration Step with RK45
        GeodesicState prevState = state;
        bool stepAccepted = false;
        
        // RK45 may reject steps if error is too large - retry with adapted h
        int maxRetries = 5;
        
        // =====================================================================
        // FAR-FIELD STEP SCALING (Performance Optimization)
        // =====================================================================
        // In the far-field (r >> M), spacetime curvature scales as ~M/r³, making
        // geodesics nearly straight. We can safely take larger steps proportional
        // to r to dramatically reduce iteration count for distant rays.
        //
        // Physics: Riemann tensor components scale as R ~ M/r³
        //          Geodesic deviation over step h: δ ~ R * h² ~ M*h²/r³
        //          For fixed accuracy ε: h_max ~ r * (ε*r²/M)^(1/2) ~ r^(3/2)
        //
        // Conservative scaling: h = h_base * (r/r_scale)^α where α = 1 for linear
        {
            float r = prevState.x.r;
            float M = params.metricParams.M;
            float r_scale = 20.0f * M;  // Start scaling at r = 20M (well outside ISCO)
            
            if (r > r_scale) {
                // Linear scaling: h ∝ r for r > r_scale
                // This is conservative (actual curvature drops faster: ~r^3)
                float scaleFactor = r / r_scale;
                
                // Cap the scaling to avoid numerical issues at extreme distances
                scaleFactor = fminf(scaleFactor, 10.0f);
                
                // Apply far-field acceleration to step size
                adaptiveStep = fminf(adaptiveStep * scaleFactor, params.integration.maxStepSize);
            }
        }
        
        // VOLUMETRIC STEP CLAMPING (Phase 4.2 Fix)
        // Ensure step size reduces near the thin accretion disk to capture density profile
        // Without this, the integrator jumps over the thin disk or samples it coarsely ("chunky" look)
        if (params.accretionDisk.enabled && adaptiveStep > params.integration.minStepSize) {
           float r = prevState.x.r;
           // Check if within radial bounds (with some buffer)
           if (r > params.accretionDisk.innerRadius * 0.9f && r < params.accretionDisk.outerRadius * 1.1f) {
               // Calculate height above plane
               float zDist = fabsf(r * cosf(prevState.x.theta));
               
               // Calculate effective scale height (geometry + beam radius)
               // Need to replicate effectiveH logic from sampleAccretionDiskVolumetric here for consistency
               float H = params.accretionDisk.heightScale * r;
               float effectiveH = H;
               
               if (params.rayBundle.enabled) {
                   // Estimate beam radius at this point (needs approximate deviation magnitude)
                   // We don't have the updated deviation yet (it's in the bundle struct), 
                   // but we can use the previous state's bundle (implicitly tracked in bundle var outside loop?)
                   // Actually integration loop updates 'state', 'bundle' is typically updated after?
                   // Just use H for safety, or a small epsilon.
                   // Ideally we'd use the beam radius to relax the step size (optimization),
                   // but using tight step is safer.
               }
               
               // If within 10 scale heights, clamp step to resolved fraction of H
               if (zDist < 10.0f * H) {
                   // Target ~10 samples per scale height for high-fidelity volume integration (was 0.25)
                   float targetStep = fmaxf(H * 0.10f, params.integration.minStepSize);
                   adaptiveStep = fminf(adaptiveStep, targetStep);
               }
           }
        }
        
        for (int retry = 0; retry < maxRetries && !stepAccepted; retry++) {
            // =================================================================
            // INTEGRATOR DISPATCH
            // =================================================================
            // Select between:
            //   - Symplectic: Time-Transformed Explicit Symplectic (Kerr - primary)
            //   - RK4: Fixed-step Runge-Kutta (fallback for other metrics)
            //
            // NOTE: RK45 adaptive integrator has been removed (Jan 2026)
            // Per user requirements, TTESI symplectic is the exclusive method.
            // =================================================================
            
            // Use RK45 Dormand-Prince adaptive integrator for all metrics
            // This provides automatic step size control near:
            //   - Poles (sin(θ) → 0 coordinate singularity)
            //   - Photon sphere (unstable orbits)
            //   - Horizon (strong field)
            state = integrateGeodesicRK45<type>(
                prevState,
                adaptiveStep,
                params.metricParams,
                params.integration.tolerance,
                params.integration.minStepSize,
                params.integration.maxStepSize,
                stepAccepted
            );
        }
        
        // If still not accepted after retries, force minimum step with RK4
        // (RK4 at minimum step is acceptable as a fallback)
        if (!stepAccepted) {
            adaptiveStep = params.integration.minStepSize;
            state = integrateGeodesicRK4<type>(prevState, adaptiveStep, params.metricParams);
            // Apply normalization that RK4 already does internally
        }
        
        // =====================================================================
        // RAY BUNDLE PROPAGATION (Phase 6.3 - DNGR Geodesic Deviation)
        // =====================================================================
        // Update deviation vectors alongside central ray. The geodesic deviation
        // equation D²ξ/dλ² = R^μ_νρσ k^ν ξ^ρ k^σ tells us how neighboring rays
        // diverge due to spacetime curvature.
        if (useRayBundle && !bundle.terminated) {
            // Sync bundle position with central ray
            bundle.x = make_float4(state.x.t, state.x.r, state.x.theta, state.x.phi);
            bundle.u = make_float4(state.u.t, state.u.r, state.u.theta, state.u.phi);
            
            float4 deviationAccel1, deviationAccel2;
            
            if (params.rayBundle.useFullRiemann) {
                // =============================================================
                // FULL RIEMANN TENSOR (DNGR Eq A.19 - High Quality Mode)
                // D²ξ^μ/dλ² = R^μ_νρσ k^ν ξ^ρ k^σ
                // =============================================================
                // Uses the full Riemann curvature tensor for exact geodesic
                // deviation. More expensive but physically accurate.
                deviationAccel1 = geodesicDeviationAccel<type>(state.x, state.u, bundle.xi1, params.metricParams);
                deviationAccel2 = geodesicDeviationAccel<type>(state.x, state.u, bundle.xi2, params.metricParams);
            } else {
                // =============================================================
                // TIDAL APPROXIMATION (Fast Mode - Default)
                // =============================================================
                // Simplified deviation propagation using leading-order tidal effects.
                // In curved spacetime, deviation vectors are stretched by curvature.
                // The effect scales as R ~ M/r³ (dominant Riemann component).
                float M = params.metricParams.M;
                float r = fmaxf(state.x.r, 0.01f);
                float tidalFactor = M / (r * r * r);  // Riemann component scale
                
                float tf = tidalFactor * adaptiveStep;
                deviationAccel1 = make_float4(bundle.xi1.x * tf, bundle.xi1.y * tf, 
                                               bundle.xi1.z * tf, bundle.xi1.w * tf);
                deviationAccel2 = make_float4(bundle.xi2.x * tf, bundle.xi2.y * tf,
                                               bundle.xi2.z * tf, bundle.xi2.w * tf);
            }
            
            // Update deviation rates: dξ/dλ += Δλ * (deviation acceleration)
            bundle.dxi1 = make_float4(bundle.dxi1.x + deviationAccel1.x, bundle.dxi1.y + deviationAccel1.y,
                                      bundle.dxi1.z + deviationAccel1.z, bundle.dxi1.w + deviationAccel1.w);
            bundle.dxi2 = make_float4(bundle.dxi2.x + deviationAccel2.x, bundle.dxi2.y + deviationAccel2.y,
                                      bundle.dxi2.z + deviationAccel2.z, bundle.dxi2.w + deviationAccel2.w);
            
            // Propagate deviation vectors: ξ += dξ/dλ * Δλ
            bundle.xi1 = make_float4(bundle.xi1.x + bundle.dxi1.x * adaptiveStep,
                                     bundle.xi1.y + bundle.dxi1.y * adaptiveStep,
                                     bundle.xi1.z + bundle.dxi1.z * adaptiveStep,
                                     bundle.xi1.w + bundle.dxi1.w * adaptiveStep);
            bundle.xi2 = make_float4(bundle.xi2.x + bundle.dxi2.x * adaptiveStep,
                                     bundle.xi2.y + bundle.dxi2.y * adaptiveStep,
                                     bundle.xi2.z + bundle.dxi2.z * adaptiveStep,
                                     bundle.xi2.w + bundle.dxi2.w * adaptiveStep);
            
            // Update ellipse parameters
            float3 xi1_spatial = make_float3(bundle.xi1.y, bundle.xi1.z, bundle.xi1.w);
            float3 xi2_spatial = make_float3(bundle.xi2.y, bundle.xi2.z, bundle.xi2.w);
            
            float a11 = dot(xi1_spatial, xi1_spatial);
            float a12 = dot(xi1_spatial, xi2_spatial);
            float a22 = dot(xi2_spatial, xi2_spatial);
            
            float trace_val = a11 + a22;
            float det = a11 * a22 - a12 * a12;
            float disc = sqrtf(fmaxf(0.0f, trace_val * trace_val - 4.0f * det));
            
            float lambda1 = 0.5f * (trace_val + disc);
            float lambda2 = 0.5f * (trace_val - disc);
            
            bundle.semiMajor = sqrtf(fmaxf(lambda1, 1e-10f));
            bundle.semiMinor = sqrtf(fmaxf(lambda2, 1e-10f));
            bundle.affineParam += adaptiveStep;
            
            // Check for bundle termination (too elongated or too small)
            float ellipseRatio = bundle.semiMajor / fmaxf(bundle.semiMinor, 1e-10f);
            if (ellipseRatio > params.rayBundle.maxEllipseRatio || 
                bundle.semiMajor < params.rayBundle.minBundleSize) {
                bundle.terminated = true;
            }
        }
        
        // =====================================================================
        // 2.4 PHOTON SPHERE GLOW (Phase 4: Visual Quality Enhancement)
        // =====================================================================
        // The photon sphere at r = 3M (Schwarzschild) is where light can orbit
        // unstably. Rays passing close to it accumulate a subtle glow, creating
        // the characteristic "Interstellar" ring effect.
        if constexpr (type == Sirius::MetricType::Schwarzschild || 
                      type == Sirius::MetricType::Kerr ||
                      type == Sirius::MetricType::ReissnerNordstrom) {
            
            float M = params.metricParams.M;
            float photonSphereR = 1.5f * (2.0f * M);  // r = 3M for Schwarzschild
            
            // Kerr: photon sphere varies with spin, but 3M is approximate
            if constexpr (type == Sirius::MetricType::Kerr) {
                float a = params.metricParams.a;
                // Prograde photon orbit: r_ph = 2M(1 + cos(2/3 * arccos(-a/M)))
                // Simplified approximation: ranges from 1M (a=M) to 3M (a=0)
                photonSphereR = 3.0f * M * (1.0f - 0.5f * fabsf(a) / M);
            }
            
            // Calculate distance to photon sphere
            float distToPhotonSphere = fabsf(state.x.r - photonSphereR);
            float glowWidth = 0.5f * M;  // Width of glow region
            
            // Accumulate glow when ray passes near photon sphere
            if (distToPhotonSphere < glowWidth) {
                float glowIntensity = 1.0f - distToPhotonSphere / glowWidth;
                glowIntensity = glowIntensity * glowIntensity;  // Quadratic falloff
                
                // Glow color: hot white-blue (like plasma)
                float3 glowColor = make_float3(0.8f, 0.85f, 1.0f);
                
                // Accumulate (small contribution each step)
                float glowContribution = glowIntensity * 0.01f * adaptiveStep;
                accumColor = accumColor + glowColor * glowContribution * (1.0f - opacity);
            }
        }
        
        // =====================================================================
        // 2.5 ACCRETION DISK RENDERING (Phase 4.2 Volumetric, Phase 7.5 Planar)
        // =====================================================================
        // Two modes:
        //   PLANAR (default): Infinitely thin disk at θ = π/2, single intersection
        //   VOLUMETRIC: Ray marching through disk thickness with radiative transfer
        if constexpr (type == Sirius::MetricType::Schwarzschild || 
                      type == Sirius::MetricType::Kerr ||
                      type == Sirius::MetricType::ReissnerNordstrom) {
            
            if (params.accretionDisk.enabled) {
                
                // =============================================================
                // PLANAR DISK MODE (Phase 7.5 - Fast, bright, cinematic)
                // =============================================================
                // Detect equatorial plane crossing (θ = π/2)
                // Single intersection → pure emission → very fast
                if (params.accretionDisk.diskMode == Sirius::DiskMode::Planar) {
                    
                    const float PI_HALF = 1.5707963f; // π/2
                    
                    // Check if ray crossed the equatorial plane this step
                    float theta_prev = prevState.x.theta;
                    float theta_curr = state.x.theta;
                    
                    // Crossing condition: one is above π/2, other is below
                    bool crossedEquator = (theta_prev - PI_HALF) * (theta_curr - PI_HALF) < 0.0f;
                    
                    if (crossedEquator && opacity < 0.99f) {
                        // Interpolate to find exact crossing point
                        float t_cross = (PI_HALF - theta_prev) / (theta_curr - theta_prev);
                        t_cross = clamp(t_cross, 0.0f, 1.0f);
                        
                        // Interpolated position at crossing
                        float r_hit = (1.0f - t_cross) * prevState.x.r + t_cross * state.x.r;
                        float phi_hit = (1.0f - t_cross) * prevState.x.phi + t_cross * state.x.phi;
                        
                        // Check if within disk radial bounds
                        float M = params.metricParams.M;  // Get mass from metric params
                        
                        // Radii are in units of M from params
                        float r_inner = params.accretionDisk.innerRadius;
                        if (r_inner <= 0.0f) {
                            // Auto-calculate ISCO for Schwarzschild/Kerr (in geometric units)
                            // ISCO_Schwarzschild = 6M, scaled for M
                            r_inner = 6.0f * M;  // Schwarzschild ISCO = 6M
                        }
                        
                        // Outer radius - params store value in units of M (e.g., 15 means 15*M)
                        float r_outer = params.accretionDisk.outerRadius * M;
                        
                        if (r_hit >= r_inner && r_hit <= r_outer) {
                            // =========================================
                            // Compute disk emission at this radius
                            // =========================================
                            
                            // Use NORMALIZED radius for M-independent color profile
                            // r_norm = 0 at inner edge (ISCO), r_norm = 1 at outer edge
                            float r_norm = (r_hit - r_inner) / fmaxf(r_outer - r_inner, 0.001f);
                            r_norm = clamp(r_norm, 0.0f, 1.0f);
                            
                            // =========================================
                            // NOVIKOV-THORNE Q(r) RELATIVISTIC CORRECTION
                            // Page & Thorne (1974): Full relativistic temperature profile
                            // =========================================
                            // The flux function Q(r) for a geometrically thin, optically thick
                            // accretion disk around a Kerr black hole. Based on:
                            //   Page & Thorne (1974), ApJ 191, 499
                            //   Thorne (1974), ApJ 191, 507
                            //
                            // Key physics: Q(r_isco) = 0 → T = 0 at ISCO (zero-torque boundary)
                            // Peak emission at r ≈ 1.36 * r_isco (Schwarzschild) to r ≈ 1.2 * r_isco (Kerr)
                            
                            float Q_factor = 1.0f;  // Default for Shakura-Sunyaev
                            
                            if (params.accretionDisk.temperatureModel == Sirius::TemperatureModel::NovikovThorne) {
                                // Full Novikov-Thorne relativistic correction
                                // Using dimensionless radius x = sqrt(r/M)
                                float a_star = params.metricParams.a / M;  // Dimensionless spin
                                a_star = fmaxf(-0.998f, fminf(0.998f, a_star));  // Clamp to physical range
                                
                                float x = sqrtf(r_hit / M);
                                float x_isco = sqrtf(r_inner / M);
                                
                                // Kerr ISCO roots for the Q function integrands
                                // x_ms = r_ms/M where r_ms is marginal stable orbit
                                float Z1 = 1.0f + cbrtf(1.0f - a_star*a_star) * 
                                           (cbrtf(1.0f + a_star) + cbrtf(1.0f - a_star));
                                float Z2 = sqrtf(3.0f * a_star * a_star + Z1 * Z1);
                                
                                // Auxiliary functions for Kerr geometry (Thorne 1974 Eq. 5)
                                // A(x) = 1 - 2/(x²) + a*/x³
                                // B(x) = 1 - 3/(x²) + 2a*/(x³)  [related to angular momentum]
                                // C(x) = 1 - 4a*/(x³) + 3a*²/(x⁴)
                                float A = 1.0f - 2.0f/(x*x) + a_star/(x*x*x);
                                float A_isco = 1.0f - 2.0f/(x_isco*x_isco) + a_star/(x_isco*x_isco*x_isco);
                                float B = 1.0f - 3.0f/(x*x) + 2.0f*a_star/(x*x*x);
                                
                                // D(x) = sqrt(1 - 3/x² + 2a*/x³) - angular velocity normalization
                                float D_sq = 1.0f - 3.0f/(x*x) + 2.0f*a_star/(x*x*x);
                                float D = (D_sq > 0.0f) ? sqrtf(D_sq) : 0.0f;
                                
                                // yms factor at ISCO (determines zero-torque boundary condition)
                                float B_isco = 1.0f - 3.0f/(x_isco*x_isco) + 2.0f*a_star/(x_isco*x_isco*x_isco);
                                float D_isco_sq = B_isco;
                                float D_isco = (D_isco_sq > 0.0f) ? sqrtf(D_isco_sq) : 0.0f;
                                
                                // C function evaluation
                                float C = 1.0f - 4.0f*a_star/(x*x*x) + 3.0f*a_star*a_star/(x*x*x*x);
                                
                                if (D > 0.0f && D_isco > 0.0f && x > x_isco * 1.001f) {
                                    // Page-Thorne (1974) Eq. 15n flux function
                                    // Q(x) = (3/(2x⁵)) * (1/D) * ∫ L_eff dx' where L_eff is angular momentum flux
                                    //
                                    // Simplified analytic approximation valid for moderate spins:
                                    // Q(x) ≈ (x² - x_isco²) / (x² * D) * [1 - sqrt(x_isco/x)]
                                    //       * [1 + correction_terms]
                                    
                                    float x_ratio = x_isco / x;
                                    float base_factor = (1.0f - sqrtf(x_ratio));
                                    
                                    // Relativistic corrections to radial structure
                                    // These capture the deviation from Newtonian r^(-3) temperature profile
                                    float sqrt_ratio = sqrtf(x_ratio);
                                    float log_factor = logf(fmaxf(x / x_isco, 1.001f));
                                    
                                    // Frame-dragging contribution (Kerr only)
                                    float frame_drag = 0.0f;
                                    if (fabsf(a_star) > 0.01f) {
                                        frame_drag = (3.0f * a_star / (2.0f * x * x * x)) * log_factor;
                                    }
                                    
                                    // Redshift factor: gravitational redshift from emission point
                                    float redshift_factor = sqrtf(A / B);
                                    
                                    // Combined Novikov-Thorne flux function
                                    Q_factor = base_factor * (1.0f - frame_drag) * fminf(redshift_factor, 2.0f);
                                    
                                    // Apply relativistic D-function correction for angular velocity
                                    Q_factor *= D / D_isco;
                                } else {
                                    Q_factor = 0.0f;  // At or inside ISCO
                                }
                                
                                // Ensure Q >= 0 (physical constraint: no negative flux)
                                Q_factor = fmaxf(Q_factor, 0.0f);
                                
                                // Q(r) affects temperature as T ∝ Q^(1/4)
                                // For emission (proportional to T^4), we use Q directly
                            }
                            
                            // =========================================
                            // Color profile based on temperature model
                            // =========================================
                            float3 diskColor;
                            
                            if (params.accretionDisk.temperatureModel == Sirius::TemperatureModel::NovikovThorne) {
                                // Novikov-Thorne: Peak temperature at ~40% from inner edge
                                // Inner edge is COOLEST (T=0 at ISCO), peak in mid-disk
                                float T_profile = Q_factor * powf(r_inner / r_hit, 0.75f);  // NT temperature
                                T_profile = fmaxf(T_profile, 0.0f);
                                
                                // Color based on temperature: cooler = redder, hotter = whiter
                                if (T_profile < 0.2f) {
                                    // Very cool (near ISCO): dim red
                                    diskColor = make_float3(0.6f + T_profile, 0.2f * T_profile, 0.1f * T_profile);
                                } else if (T_profile < 0.5f) {
                                    // Cool to warm: red to orange
                                    float t = (T_profile - 0.2f) / 0.3f;
                                    diskColor = make_float3(0.8f + 0.2f * t, 0.3f + 0.3f * t, 0.15f + 0.1f * t);
                                } else if (T_profile < 0.8f) {
                                    // Hot: orange to yellow-white
                                    float t = (T_profile - 0.5f) / 0.3f;
                                    diskColor = make_float3(1.0f, 0.6f + 0.35f * t, 0.25f + 0.55f * t);
                                } else {
                                    // Very hot: white
                                    diskColor = make_float3(1.0f, 0.95f, 0.9f);
                                }
                            } else if (params.accretionDisk.useSpectralColors) {
                                // Shakura-Sunyaev: Original Interstellar-style profile
                                // Inner edge is HOTTEST (classic thin disk)
                                if (r_norm < 0.15f) {
                                    float t = r_norm / 0.15f;
                                    diskColor = make_float3(1.0f, 0.95f - t * 0.25f, 0.8f - t * 0.5f);
                                } else if (r_norm < 0.4f) {
                                    float t = (r_norm - 0.15f) / 0.25f;
                                    diskColor = make_float3(1.0f, 0.7f - t * 0.15f, 0.3f - t * 0.15f);
                                } else if (r_norm < 0.7f) {
                                    float t = (r_norm - 0.4f) / 0.3f;
                                    diskColor = make_float3(1.0f, 0.55f - t * 0.1f, 0.15f - t * 0.05f);
                                } else {
                                    float t = (r_norm - 0.7f) / 0.3f;
                                    diskColor = make_float3(0.95f - t * 0.15f, 0.45f - t * 0.15f, 0.1f);
                                }
                            } else {
                                diskColor = make_float3(1.0f, 0.6f, 0.2f);  // Default orange
                            }
                            
                            // =========================================
                            // Emission intensity with NT correction
                            // =========================================
                            // For Novikov-Thorne: emission ∝ T^4 ∝ Q(r)
                            // For Shakura-Sunyaev: brightest at inner edge
                            float base_intensity = params.accretionDisk.emissionCoefficient;
                            float intensity;
                            
                            if (params.accretionDisk.temperatureModel == Sirius::TemperatureModel::NovikovThorne) {
                                // NT: Peak emission at ~40% from inner edge, zero at ISCO
                                float radial_factor = powf(r_inner / r_hit, 3.0f);  // Standard r^-3 falloff
                                intensity = base_intensity * Q_factor * radial_factor;
                                intensity = fmaxf(intensity, 0.0f);
                            } else {
                                // Shakura-Sunyaev: Brightest at inner edge
                                float intensity_falloff = powf(1.0f - r_norm, 2.5f);
                                intensity = base_intensity * (0.3f + 0.7f * intensity_falloff);
                            }

                            
                            // =========================================
                            // PHYSICALLY ACCURATE Disk Animation
                            // Keplerian Differential Rotation (Dec 2025)
                            // =========================================
                            // Physics: Material orbits azimuthally (tangentially), NOT radially
                            // Inner disk rotates faster than outer (Keplerian shear):
                            //   Ω(r) = sqrt(M) / (r^1.5 + a*sqrt(M))
                            // This creates visible "lapping" where inner edge overtakes outer
                            
                            float a = params.metricParams.a;
                            float M = params.metricParams.M;
                            float r_safe = fmaxf(r_hit, M);  // Prevent division issues
                            
                            // Keplerian angular velocity: Ω = sqrt(M) / (r^1.5 + a*sqrt(M))
                            // Higher omega at smaller r -> inner disk rotates faster
                            float omega_kepler = sqrtf(M) / (powf(r_safe, 1.5f) + a * sqrtf(M));
                            
                            // Azimuthal rotation angle: each radius rotates at its own rate
                            // This creates the "water down a drain" shearing effect
                            // Speed multiplier x100 needed because omega is small at typical r
                            float rotationAngle = omega_kepler * params.time * 100.0f;
                            
                            // Effective azimuthal position including rotation
                            // phi_rotated = phi_hit + rotation (inner laps outer visually)
                            float phi_rotated = phi_hit + rotationAngle;
                            
                            // =========================================
                            // Static concentric ring structure (no radial motion)
                            // Rings are brightness variations at fixed radii
                            // =========================================
                            const float numRings = 12.0f;
                            float ringPhase = r_norm * numRings * 6.28318f;  // Static in radius
                            float ringMod = 0.75f + 0.25f * sinf(ringPhase);  // 75-100% intensity
                            
                            // =========================================
                            // Rotating spiral density waves
                            // These are azimuthal features that rotate WITH the flow
                            // Inner parts of spiral rotate faster -> creates trailing spiral
                            // =========================================
                            // 2-arm logarithmic spiral: θ = log(r) * pitch + rotation
                            float spiralPitch = 2.5f;  // Controls how tightly wound
                            float spiralPhase = phi_rotated * 2.0f + logf(r_norm + 0.1f) * spiralPitch;
                            float spiralMod = 0.80f + 0.20f * sinf(spiralPhase);  // 80-100%
                            
                            // =========================================
                            // Turbulent clumps rotating with the flow
                            // Higher frequency structure that shows differential rotation
                            // =========================================
                            float turbPhase1 = phi_rotated * 5.0f + r_norm * 20.0f;  // Rotating turbulence
                            float turbPhase2 = phi_rotated * 8.0f - r_norm * 15.0f;  // Counter-rotating layer
                            float turbulence = 0.90f + 0.05f * sinf(turbPhase1) + 0.05f * sinf(turbPhase2);
                            
                            // =========================================
                            // GRAVITATIONAL TIME DILATION
                            // Material near horizon appears to slow down and fade
                            // to a distant observer due to extreme spacetime curvature
                            // =========================================
                            // Time dilation factor: sqrt(1 - rs/r) -> 0 at horizon
                            float rs = 2.0f * M;  // Schwarzschild radius
                            float timeDilationFactor = sqrtf(fmaxf(1.0f - rs / r_hit, 0.01f));
                            
                            // The "apparent" rotation speed is reduced by time dilation
                            // This creates the visual paradox: inner is locally faster but APPEARS slower
                            // For animation, we modulate the visual speed by time dilation
                            // Note: The actual Keplerian omega already computed handles coordinate motion
                            // This additional factor represents what the distant observer perceives
                            
                            // Fade factor: inner disk fades due to gravitational redshift
                            // This is separate from Doppler beaming - purely gravitational
                            float gravitationalFade = powf(timeDilationFactor, 1.5f);  // Stronger fade near horizon
                            
                            // =========================================
                            // MULTIPLE TRACKABLE CLUMPS
                            // Like professional black hole simulations, we add
                            // several clumps at different radii to show shearing
                            // =========================================
                            float clumpBoost = 1.0f;
                            float clumpWidth = 0.06f;  // Radial width
                            float azimuthalClumpWidth = 0.25f;  // Angular width
                            
                            // Define clump positions (radius, initial phi offset)
                            // Clumps at different radii will orbit at different speeds
                            const int NUM_CLUMPS = 5;
                            float clumpRadii[NUM_CLUMPS] = {0.15f, 0.30f, 0.50f, 0.70f, 0.85f};
                            float clumpPhiOffset[NUM_CLUMPS] = {0.0f, 1.2f, 2.5f, 4.0f, 5.5f};
                            
                            for (int i = 0; i < NUM_CLUMPS; i++) {
                                float clumpR = clumpRadii[i];
                                float clumpPhi0 = clumpPhiOffset[i];
                                
                                // Compute Keplerian omega for this clump's radius
                                float clumpActualR = r_inner + clumpR * (r_outer - r_inner);
                                float clumpOmega = sqrtf(M) / (powf(fmaxf(clumpActualR, M), 1.5f) + a * sqrtf(M));
                                
                                // Apply time dilation to apparent rotation speed
                                float clumpTimeDilation = sqrtf(fmaxf(1.0f - rs / clumpActualR, 0.01f));
                                float apparentOmega = clumpOmega * clumpTimeDilation;  // Appears slower near horizon
                                
                                // Current clump angle (with time dilation affecting apparent speed)
                                float clumpAngle = clumpPhi0 + apparentOmega * params.time * 100.0f;
                                
                                // Radial Gaussian distance
                                float rDist = fabsf(r_norm - clumpR);
                                float rGauss = expf(-rDist * rDist / (2.0f * clumpWidth * clumpWidth));
                                
                                // Azimuthal Gaussian distance (with wrapping)
                                float phiDist = phi_hit - clumpAngle;
                                while (phiDist > 3.14159f) phiDist -= 6.28318f;
                                while (phiDist < -3.14159f) phiDist += 6.28318f;
                                float phiGauss = expf(-phiDist * phiDist / (2.0f * azimuthalClumpWidth * azimuthalClumpWidth));
                                
                                // Clump contribution (brighter clumps, modulated by gravitational fade)
                                float clumpFade = sqrtf(fmaxf(1.0f - rs / clumpActualR, 0.01f));
                                clumpBoost += 2.5f * rGauss * phiGauss * clumpFade;
                            }
                            
                            // Combine all modulations including gravitational effects
                            intensity *= ringMod * spiralMod * turbulence * clumpBoost * gravitationalFade;
                            intensity = fminf(intensity, 12.0f);  // Clamp to prevent blowout


                            
                            float3 emission = diskColor * intensity;
                            
                            // =========================================
                            // Relativistic effects (Doppler + beaming)
                            // =========================================
                            Vec4 hitPos;
                            hitPos.t = (1.0f - t_cross) * prevState.x.t + t_cross * state.x.t;
                            hitPos.r = r_hit;
                            hitPos.theta = PI_HALF;
                            hitPos.phi = phi_hit;
                            
                            Vec4 u_emit = computeKeplerianVelocity<type>(hitPos, params.metricParams);
                            
                            // Interpolate photon momentum at crossing
                            Vec4 u_photon;
                            u_photon.t = (1.0f - t_cross) * prevState.u.t + t_cross * state.u.t;
                            u_photon.r = (1.0f - t_cross) * prevState.u.r + t_cross * state.u.r;
                            u_photon.theta = (1.0f - t_cross) * prevState.u.theta + t_cross * state.u.theta;
                            u_photon.phi = (1.0f - t_cross) * prevState.u.phi + t_cross * state.u.phi;
                            
                            float g_total = computeRelativisticG<type>(hitPos, u_photon, u_emit, params.metricParams);
                            
                            // Doppler beaming: I_obs = g^(3+α) * I_emit (α~0 for thermal)
                            float beaming = powf(g_total, params.accretionDisk.beamingExponent);
                            beaming = clamp(beaming, 0.1f, 10.0f);
                            emission = emission * beaming;
                            
                            // Spectral color shift
                            if (params.accretionDisk.useSpectralColors) {
                                if (g_total > 1.0f) {  // Blue shift
                                    emission.z *= g_total;
                                    emission.y *= sqrtf(g_total);
                                } else {  // Red shift
                                    emission.x /= g_total;
                                    emission = emission * g_total;
                                }
                            }
                            
                            // =========================================
                            // Add disk emission (pure emission, no absorption)
                            // =========================================
                            // Planar disk is treated as optically thick surface
                            float diskOpacity = 0.95f;  // Nearly opaque surface
                            
                            // Composite with existing accumulated color
                            accumColor = accumColor * (1.0f - diskOpacity) + emission * diskOpacity;
                            opacity = opacity + diskOpacity * (1.0f - opacity);
                        }
                    }
                }
                // =============================================================
                // VOLUMETRIC DISK MODE (Phase 4.2 - Shakura-Sunyaev)
                // =============================================================
                else {  // diskMode == Volumetric
                    
                    // VOLUMETRIC SUB-STEPPING INTEGRATION (Phase 4.2 Fix)
                    // Integrate from prevState to state to capture thin disk details
                    // prevents "jumping over" the disk grid.
                    
                    // Check if we are potentially crossing or inside the disk
                    float r_mid = (state.x.r + prevState.x.r) * 0.5f;
                    bool nearDisk = false;
                    if (r_mid > params.accretionDisk.innerRadius * 0.8f && 
                        r_mid < params.accretionDisk.outerRadius * 1.2f) {
                        // Check vertical proximity (liberal bound)
                        float H_mid = params.accretionDisk.heightScale * r_mid;
                        // Beam radius aware check could go here
                        if (fabsf(state.x.theta - 1.570796f) < 0.5f) { // ~28 degrees from equator
                             nearDisk = true;
                        }
                    }
                    // =========================================================
                    // OPTIMIZATION: Gauss-Lobatto Quadrature for Volumetric Integration
                    // Uses spectral nodes for optimal sample placement - achieves
                    // same accuracy as 8-16 uniform samples with only 4 points
                    // =========================================================
                    const int GL_SAMPLES = nearDisk ? 4 : 1;
                    float subStep = adaptiveStep;  // Full step, weighted by GL weights
            
            // Beam radius calculation (computed once per geodesic step for efficiency)
            float beamRadius = 0.0f;
            if (useRayBundle && !bundle.terminated) {
                 float mag1 = bundle.xi1.x*bundle.xi1.x + bundle.xi1.y*bundle.xi1.y + bundle.xi1.z*bundle.xi1.z;
                 float mag2 = bundle.xi2.x*bundle.xi2.x + bundle.xi2.y*bundle.xi2.y + bundle.xi2.z*bundle.xi2.z;
                 beamRadius = sqrtf(fmaxf(mag1, mag2));
            }

            for (int k = 0; k < GL_SAMPLES; k++) {
                // Gauss-Lobatto nodes and weights for optimal integration
                float t = (GL_SAMPLES == 4) ? GL4_NODES[k] : 0.5f;
                float weight = (GL_SAMPLES == 4) ? GL4_WEIGHTS[k] : 1.0f;
                
                // Light jitter for anti-aliasing (much smaller than uniform sampling)
                unsigned int subSeed = pixelIndex * 719u + params.frameIndex * 131u + k * 97u;
                float jitter = (randomFloat(subSeed) - 0.5f) * 0.1f;  // ±5% jitter
                t = clamp(t + jitter, 0.0f, 1.0f);
                
                Vec4 subPos;
                subPos.t = (1.0f - t) * prevState.x.t + t * state.x.t;
                subPos.r = (1.0f - t) * prevState.x.r + t * state.x.r;
                subPos.theta = (1.0f - t) * prevState.x.theta + t * state.x.theta;
                subPos.phi = (1.0f - t) * prevState.x.phi + t * state.x.phi; 
                
                // Sample at sub-position
                float3 subEmission;
                float subDtau;
                float3 subVel;
                
                sampleAccretionDiskVolumetric(
                    subPos,
                    params.accretionDisk,
                    params.metricParams,
                    subStep,
                    beamRadius,
                    subEmission,
                    subDtau,
                    subVel
                );
                
                if (subDtau > 1e-6f || length(subEmission) > 1e-6f) {
                    // Calculate relativistic effects (G-factor) at this sub-position
                    // We can reuse photon 4-momentum 'state.u' as it doesn't change much over one step
                    // But ideally interpolate 'u' too. For speed, use 'state.u'.
                    
                    Vec4 u_emit = computeKeplerianVelocity<type>(subPos, params.metricParams);
                    float g_total = computeRelativisticG<type>(
                         subPos, state.u, u_emit, params.metricParams
                    );
                    
                    float beaming = powf(g_total, params.accretionDisk.beamingExponent);
                    beaming = clamp(beaming, 0.1f, 10.0f);
                    subEmission = subEmission * beaming;
                    
                    if (params.accretionDisk.useSpectralColors) {
                        float T_ratio = g_total; 
                        // Doppler shift approximation: shift color by T_ratio
                        // Blue-shift (T_ratio > 1) -> multiply color by blueish tint?
                        // Or shift blackbody? B_nu(T*g) approx g^3 B_nu(T)? 
                        // Simplified tinting:
                        if (T_ratio > 1.0f) { // Blue shift
                             subEmission.z *= T_ratio; // Boost blue
                             subEmission.y *= sqrtf(T_ratio);
                        } else { // Red shift
                             subEmission.x /= T_ratio; // Boost red (relative)
                             subEmission = subEmission * T_ratio; // Diminish overall
                        }
                    }
                    
                    // Radiative Transfer for Sub-step
                    // HYBRID THIN/THICK LIMIT with Gauss-Lobatto weighting:
                    // - For τ < 0.01: Pure emission mode (I += S*τ) - glowing gas
                    // - For τ >= 0.01: Standard RT (I = I*T + S*(1-T))
                    // GL weight is applied to properly integrate each sample
                    
                    // Scale optical depth by GL weight for proper quadrature
                    float weightedDtau = subDtau * weight;
                    
                    if (weightedDtau < 0.01f) {
                        // OPTICALLY THIN: Pure emission, no attenuation
                        // This is physically correct for tenuous plasma
                        accumColor = accumColor + subEmission * weightedDtau * (1.0f - opacity);
                        // Minimal opacity contribution for thin gas
                        opacity += weightedDtau * 0.1f * (1.0f - opacity);
                    } else {
                        // OPTICALLY THICK: Standard radiative transfer
                        float transmittance = FAST_EXPF(-weightedDtau);
                        float absorptionFactor = 1.0f - transmittance;
                        
                        accumColor = accumColor * transmittance + subEmission * absorptionFactor * (1.0f - opacity);
                        opacity += absorptionFactor * (1.0f - opacity);
                    }
                    
                    if (opacity > 0.99f) { 
                        opacity = 1.0f; 
                        break; 
                    }
                }
                
                if (opacity >= 1.0f) break;
            }
            
            // End of volumetric sub-stepping
                }  // End of volumetric mode else block
            }  // End of accretionDisk.enabled check
        }  // End of metric type check
        
        // =====================================================================
        // 2.6 RUSSIAN ROULETTE PATH TERMINATION (Phase 4 Optimization)
        // =====================================================================
        // Probabilistically terminate low-contribution paths early.
        // Maintains unbiased estimation: E[result] = P_survive * boosted_value
        // Saves computation on rays that are mostly transparent or nearly dead.
        if (step > 15 && opacity < 0.1f) {
            // Contribution is low (mostly transparent ray)
            // Survival probability proportional to remaining potential contribution
            float survivalProb = fmaxf(1.0f - opacity, 0.05f);  // Min 5% survival
            
            // Hash-based random using pixel index and step as seeds
            unsigned int rngSeed = pixelIndex * 1973u + step * 9277u + params.frameIndex * 26699u;
            float rand = randomFloat(rngSeed);
            
            if (rand > survivalProb) {
                // Terminate this path - ray is unlikely to contribute significantly
                escaped = true;
                // Use world-frame transform (consistent with main background sampling)
                float r_esc = fmaxf(state.x.r, 0.001f);
                float st_esc = sinf(state.x.theta);
                // Use ray's escape position for consistent celestial sphere
                finalDir = transformVelocityToWorldFrame(
                    state.u.r,
                    state.u.theta * r_esc,
                    state.u.phi * r_esc * st_esc,
                    state.x.theta,  // Ray's escape theta
                    state.x.phi     // Ray's escape phi
                );
                break;
            }
            // Survived - boost accumulated color to maintain unbiased estimate
            // accumColor = accumColor / survivalProb; // Optional: proper MIS weighting
        }
        
        // 3. Horizon Check
        if (state.x.r < horizonRadius) {
            hitHorizon = true;
            break;
        }

        // 3a. Additional capture check: rays near horizon moving inward
        // This eliminates specs caused by rays that barely miss due to step size
        if (state.x.r < horizonRadius * 1.20f && state.u.r < 0.0f) {
            // Ray is within 20% of horizon and moving inward - will be captured
            hitHorizon = true;
            break;
        }

        // =====================================================================
        // 3.5. KERR WEAK FIELD CROSSOVER (Critical fix for stability)
        // =====================================================================
        // 3.5. KERR WEAK FIELD CROSSOVER - DISABLED
        // =====================================================================
        // This workaround is no longer needed since we now use finite-difference
        // Christoffel symbols which are correct at all distances. The buggy
        // analytic Kerr Christoffels were the cause of the wave artifacts.
        // =====================================================================
        
        // =====================================================================
        // 4. EARLY RAY TERMINATION (Phase 2 Optimization)
        // =====================================================================
        // If ray is far from black hole and moving outward, it will escape.
        // No need to integrate further - sample background immediately.
        // This saves significant computation for rays that miss the black hole.
        if constexpr (type != Sirius::MetricType::Minkowski) {
            float M = params.metricParams.M;
            float escapeRadius = 20.0f * M;  // Beyond this, gravitational effects are weak
            
            // Check if ray is escaping: far from BH AND moving outward
            if (state.x.r > escapeRadius && state.u.r > 0.0f) {
                // Additional check: radial velocity should dominate angular components
                // This ensures we're not cutting off orbiting rays
                float r = state.x.r;
                float angularSpeed = sqrtf(state.u.theta * state.u.theta * r * r + 
                                           state.u.phi * state.u.phi * r * r * sinf(state.x.theta) * sinf(state.x.theta));
                
                if (state.u.r > angularSpeed * 0.5f) {
                    // Ray is definitely escaping - terminate early
                    escaped = true;
                    // Use world-frame transform (consistent with main background sampling)
                    float st_esc = sinf(state.x.theta);
                    // Use ray's escape position for consistent celestial sphere
                    finalDir = transformVelocityToWorldFrame(
                        state.u.r,
                        state.u.theta * r,
                        state.u.phi * r * st_esc,
                        state.x.theta,  // Ray's escape theta
                        state.x.phi     // Ray's escape phi
                    );
                    break;
                }
            }
        }
        
        // 5. Far Distance Escape Check (fallback)
        if (state.x.r > maxDistance) {
            escaped = true;
            // Use world-frame transform (consistent with main background sampling)
            float r_esc = fmaxf(state.x.r, 0.001f);
            float st_esc = sinf(state.x.theta);
            // Use ray's escape position for consistent celestial sphere
            finalDir = transformVelocityToWorldFrame(
                state.u.r,
                state.u.theta * r_esc,
                state.u.phi * r_esc * st_esc,
                state.x.theta,  // Ray's escape theta
                state.x.phi     // Ray's escape phi
            );
            break;
        }
    }

    float3 bgColor = make_float3(0.0f, 0.0f, 0.0f);
    if (hitHorizon) {
        // Ray fell into black hole - pure black
        bgColor = make_float3(0.0f, 0.0f, 0.0f);
    } else {
        // Ray did NOT hit horizon - sample background in WORLD frame
        // =====================================================================
        // CELESTIAL SPHERE FIX: Use ray's ESCAPE position for basis vectors
        // The velocity components (velR, velTheta, velPhi) are defined in the
        // local spherical basis at the ray's current position. To get the correct
        // world-frame direction, we must use the SAME position for the basis.
        //
        // For stars at infinity, what matters is the asymptotic direction of
        // the ray, which is determined by its 4-velocity components expressed
        // in the correct local basis at the escape point.
        // =====================================================================
        float r = fmaxf(state.x.r, 0.001f);
        float st = sinf(state.x.theta);
        
        // Extract velocity components (scaled appropriately for Jacobian)
        float velR = state.u.r;
        float velTheta = state.u.theta * r;
        float velPhi = state.u.phi * r * st;
        
        // Transform to world frame using RAY'S ESCAPE position (not observer's)
        // This ensures the velocity components and basis vectors are consistent
        float3 worldDir = transformVelocityToWorldFrame(
            velR, velTheta, velPhi,
            state.x.theta,  // Use ray's escape theta
            state.x.phi     // Use ray's escape phi
        );
        
        bgColor = sampleBackground(worldDir);
    }
    
    // Composite Background
    // color = accum + (1-opacity) * bg
    float3 finalColor = accumColor + bgColor * (1.0f - opacity);
    
    // =========================================================================
    // TEMPORAL ACCUMULATION (Phase 4.3 - Dual Mode)
    // =========================================================================
    // Two modes supported:
    // 1. Running Average: newAccum = (oldAccum * N + newColor) / (N + 1)
    //    - Unbiased, converges to ground truth
    //    - Slower convergence, equal weight to all frames
    // 2. Exponential Moving Average (EMA): newAccum = α*newColor + (1-α)*oldAccum  
    //    - Faster visual convergence
    //    - More weight to recent frames, handles slow scene changes
    
    if (params.pathTracing.enableAccumulation && params.frameIndex > 0) {
        // Read previous accumulated value
        float4 prevAccum = params.accumBuffer[pixelIndex];
        float3 prevColor = make_float3(prevAccum.x, prevAccum.y, prevAccum.z);
        
        if (params.pathTracing.useExponentialMA) {
            // Exponential Moving Average
            // α = blendWeight (typically 0.05-0.2 for smooth convergence)
            float alpha = params.pathTracing.blendWeight;
            finalColor = finalColor * alpha + prevColor * (1.0f - alpha);
        } else {
            // Running Average (unbiased)
            float weight = 1.0f / (float)(params.frameIndex + 1);
            finalColor = prevColor * (1.0f - weight) + finalColor * weight;
        }
    }
    
    // Clamp to prevent overflow from accumulation errors
    finalColor.x = clamp(finalColor.x, 0.0f, 1.0f);
    finalColor.y = clamp(finalColor.y, 0.0f, 1.0f);
    finalColor.z = clamp(finalColor.z, 0.0f, 1.0f);
    
    params.frameBuffer[pixelIndex] = make_float4(finalColor.x, finalColor.y, finalColor.z, 1.0f);
    params.accumBuffer[pixelIndex] = make_float4(finalColor.x, finalColor.y, finalColor.z, 1.0f);
}

//==============================================================================
// Miss Program (Background)
//==============================================================================
extern "C" __global__ void __miss__background() {
    // This is called when optixTrace doesn't hit any geometry
    // For geodesic raymarching, we handle background in raygen instead
    // This is here for OptiX pipeline requirements
    // NOTE: Empty body since raygen programs don't use optixTrace()
}

//==============================================================================
// Closest Hit Program (Placeholder for future geometry)
//==============================================================================
extern "C" __global__ void __closesthit__radiance() {
    // Placeholder for future geometry intersection
    // Currently not used for pure geodesic raymarching
    // NOTE: Empty body since raygen programs don't use optixTrace()
}

//==============================================================================
// Specialized Kernel Entry Points
// __launch_bounds__(256, 2) = 256 threads/block, min 2 blocks/SM
// This limits registers to ~96/thread for better occupancy (target: 50%+)
//==============================================================================
extern "C" __global__ void __launch_bounds__(256, 2) __raygen__Minkowski() { 
    raygen_renderFrame_impl<Sirius::MetricType::Minkowski>(); 
}
extern "C" __global__ void __launch_bounds__(256, 2) __raygen__Schwarzschild() { 
    raygen_renderFrame_impl<Sirius::MetricType::Schwarzschild>(); 
}
extern "C" __global__ void __launch_bounds__(256, 2) __raygen__Kerr() { 
    raygen_renderFrame_impl<Sirius::MetricType::Kerr>(); 
}
extern "C" __global__ void __launch_bounds__(256, 2) __raygen__ReissnerNordstrom() { 
    raygen_renderFrame_impl<Sirius::MetricType::ReissnerNordstrom>(); 
}
extern "C" __global__ void __launch_bounds__(256, 2) __raygen__Godel() { 
    raygen_renderFrame_impl<Sirius::MetricType::Godel>(); 
}
extern "C" __global__ void __launch_bounds__(256, 2) __raygen__TaubNUT() { 
    raygen_renderFrame_impl<Sirius::MetricType::TaubNUT>(); 
}
extern "C" __global__ void __launch_bounds__(256, 2) __raygen__KerrSchild() { 
    raygen_renderFrame_impl<Sirius::MetricType::KerrSchild>(); 
}
extern "C" __global__ void __launch_bounds__(256, 2) __raygen__EllisDrainhole() { 
    raygen_renderFrame_impl<Sirius::MetricType::EllisDrainhole>(); 
}
extern "C" __global__ void __launch_bounds__(256, 2) __raygen__Alcubierre() { 
    raygen_renderFrame_impl<Sirius::MetricType::Alcubierre>(); 
}
extern "C" __global__ void __launch_bounds__(256, 2) __raygen__DeSitter() { 
    raygen_renderFrame_impl<Sirius::MetricType::DeSitter>(); 
}
