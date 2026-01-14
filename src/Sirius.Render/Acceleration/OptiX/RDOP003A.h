// RDOP003A.h - OptiX Launch Parameters
//
// Shared Host/Device structures for GPU kernels.
// Defines metrics, camera, integration, accretion disk, and post-processing params.

#pragma once

// OptiX SBT record alignment constants (must be defined early)
#ifndef OPTIX_SBT_RECORD_ALIGNMENT
#define OPTIX_SBT_RECORD_ALIGNMENT 16
#endif

#ifndef OPTIX_SBT_RECORD_HEADER_SIZE
#define OPTIX_SBT_RECORD_HEADER_SIZE 32
#endif

// Handle CUDA and OptiX types for different compilation scenarios
#ifdef __CUDACC__
// Compiling with NVCC - full CUDA+OptiX support
#include <cuda_runtime.h>
#include <optix_types.h>
#elif defined(__CUDA_RUNTIME_H__)
// CUDA runtime was already included (but not NVCC) - define OptiX types manually
using OptixTraversableHandle = unsigned long long;
#else
// Pure C++ build without CUDA
struct float3 { float x, y, z; };
struct float4 { float x, y, z, w; };
struct int2 { int x, y; };
struct int3 { int x, y, z; };
using cudaTextureObject_t = unsigned long long;
using OptixTraversableHandle = unsigned long long;
using CUdeviceptr = unsigned long long;

// Alignment macro for non-CUDA builds
#define __align__(x)
#endif

namespace Sirius {

//==============================================================================
// Metric Family Enumeration (Unified Architecture - Dec 2025)
// Analytic solutions only - no numerical/finite difference support
//==============================================================================
enum class MetricFamily : uint32_t {
    // Unified Families (cover multiple spacetimes via parameters)
    KerrSchild = 0,      // Covers: Minkowski, Schwarzschild, Kerr, RN, KN, dS, SdS, KdS
    MorrisThorne = 1,    // Covers: Ellis drainhole, traversable wormholes
    WarpDrive = 2,       // Covers: Alcubierre, Natário, Van den Broeck
    
    // Special Cases (require individual analytic treatment)
    LewisWeyl = 10,      // Covers: Gödel (closed timelike curves)
    NUTCharged = 11,     // Covers: Taub-NUT (Misner string topology)
    
    // Reserved for future analytic families
    Custom = 99
};

// Legacy enum (deprecated - use MetricFamily instead)
// Retained only for compile-time template dispatch during migration
enum class MetricType : int {
    Minkowski = 0,       // → KerrSchild(M=0)
    Schwarzschild = 1,   // → KerrSchild(a=Q=0)
    Kerr = 2,            // → KerrSchild(Q=0)
    ReissnerNordstrom = 3, // → KerrSchild(a=0)
    Godel = 4,           // → LewisWeyl
    TaubNUT = 5,         // → NUTCharged
    KerrSchild = 6,      // → KerrSchild (native)
    EllisDrainhole = 7,  // → MorrisThorne
    Alcubierre = 8,      // → WarpDrive
    DeSitter = 9,        // → KerrSchild(M=a=Q=0, Λ>0)
    Custom = 99
};

//==============================================================================
// Geodesic Integrator Type
// NOTE: RK45/RK4 removed (Jan 2026) - symplectic TTESI is now exclusive
//==============================================================================
enum class IntegratorType : int {
    Verlet = 0,        // Velocity Verlet (symplectic, fast) - legacy
    Symplectic = 1     // Time-Transformed Explicit Symplectic (TTESI) - primary
};

//==============================================================================
// Ray Types for OptiX
//==============================================================================
enum class RayType : int {
    Geodesic = 0,      // Primary rays following geodesics
    Shadow = 1,        // Shadow rays (future path tracing)
    NUM_RAY_TYPES = 2
};

//==============================================================================
// Unified Metric Parameters Structure
// Uses union for family-specific parameters (GPU-compatible)
//==============================================================================

// Kerr-Schild family parameters: covers 9 spacetimes
struct KerrSchildFamilyParams {
    float M;           // Mass (0 for Minkowski/dS)
    float a;           // Spin (0 for static, |a| ≤ M for BH)
    float Q;           // Charge (0 for uncharged)
    float Lambda;      // Cosmological constant (0 for asymp. flat)
};

// Morris-Thorne family parameters: traversable wormholes
struct MorrisThorneParams {
    float b0;          // Throat radius
    float redshiftPhi; // Redshift function parameter (0 = Ellis)
    float shapeK;      // Shape function parameter
    float padding;
};

// Warp Drive family parameters: Alcubierre class
struct WarpDriveParams {
    float vs;          // Bubble velocity
    float sigma;       // Wall thickness
    float R;           // Bubble radius
    float padding;
};

// Lewis-Weyl family parameters: Gödel
struct LewisWeylParams {
    float omega;       // Rotation parameter
    float padding[3];
};

// NUT-Charged family parameters: Taub-NUT
struct NUTChargedParams {
    float M;           // Mass
    float n;           // NUT charge (gravitomagnetic monopole)
    float padding[2];
};

// Unified parameter structure with family dispatch
struct MetricParams {
    MetricFamily family;
    uint32_t padding_family;
    
    // Union for family-specific parameters
    union {
        KerrSchildFamilyParams kerrSchild;
        MorrisThorneParams morrisThorne;
        WarpDriveParams warpDrive;
        LewisWeylParams lewisWeyl;
        NUTChargedParams nutCharged;
        float raw[4];  // Raw access for legacy compatibility
    };
    
    // Derived quantities (computed from parameters)
    float rs;          // Schwarzschild radius (2M)
    float rplus;       // Outer horizon radius
    float rminus;      // Inner horizon radius
    float ergosphere;  // Ergosphere outer boundary
    
    // Legacy fields maintained for backward compatibility
    // These are synced with kerrSchild union member during initialization
    float M;           // Mass (synced with kerrSchild.M)
    float a;           // Spin (synced with kerrSchild.a)
    float Q;           // Charge (synced with kerrSchild.Q)
    float Lambda;      // Cosmological (synced with kerrSchild.Lambda)
};

//==============================================================================
// Numerical Metric Data (GPU-accessible)
// Grid-based metric from Einstein Toolkit / HDF5
//==============================================================================

//==============================================================================
// FIDO Frame (Fiducial Observer) - DNGR Paper Eq A.2-A.3
// James et al. (2015) "Gravitational Lensing by Spinning Black Holes"
//
// The FIDO is a locally non-rotating observer whose 4-velocity is orthogonal
// to the constant-t hypersurfaces. This provides the reference frame for
// measuring velocities and computing relativistic aberration.
//
// Basis vectors (orthonormal in FIDO frame):
//   e_r̂ = (√Δ/ρ)(∂/∂r)
//   e_θ̂ = (1/ρ)(∂/∂θ)
//   e_φ̂ = (1/ϖ)(∂/∂φ)
//
// where ρ² = r² + a²cos²θ, Δ = r² - 2Mr + a²,
//       Σ² = (r² + a²)² - a²Δsin²θ, ϖ = Σsinθ/ρ
//==============================================================================
struct FIDOBasis {
    // Orthonormal basis vectors in Cartesian representation
    float3 e_r;       // Radial unit vector (∂/∂r direction)
    float3 e_theta;   // Polar unit vector (∂/∂θ direction)
    float3 e_phi;     // Azimuthal unit vector (∂/∂φ direction)
    
    // Kerr metric quantities at FIDO location (Eq A.2)
    float alpha;      // Lapse function α = ρ√Δ/Σ
    float omega;      // Frame dragging angular velocity ω = 2ar/Σ²
    float varpi;      // ϖ = Σsinθ/ρ (proper circumference factor)
    
    // Derived quantities for convenience
    float rho;        // ρ = √(r² + a²cos²θ)
    float Delta;      // Δ = r² - 2Mr + a²
    float Sigma;      // Σ = √((r² + a²)² - a²Δsin²θ)
};

//==============================================================================
// Numerical Metric Data (GPU-accessible)
// Grid-based metric from Einstein Toolkit / HDF5
//==============================================================================

struct NumericalMetricData {
    // Grid dimensions
    int3 dims;              // Grid size (nx, ny, nz)
    float3 origin;          // Grid origin (x0, y0, z0)
    float3 spacing;         // Grid spacing (dx, dy, dz)
    
    // 3D texture objects for metric components (trilinear interpolation)
    cudaTextureObject_t gxx;
    cudaTextureObject_t gxy;
    cudaTextureObject_t gxz;
    cudaTextureObject_t gyy;
    cudaTextureObject_t gyz;
    cudaTextureObject_t gzz;
    cudaTextureObject_t alp;    // Lapse
    cudaTextureObject_t betax;  // Shift x
    cudaTextureObject_t betay;  // Shift y
    cudaTextureObject_t betaz;  // Shift z
    
    // Hydrodynamics (Phase 4)
    cudaTextureObject_t rho;    // Density
    cudaTextureObject_t vx;     // Velocity x
    cudaTextureObject_t vy;     // Velocity y
    cudaTextureObject_t vz;     // Velocity z
    
    // =========================================================================
    // P2 Geodesic Lookup Tables: Precomputed Christoffel Symbols
    // =========================================================================
    // Christoffel symbols Γ^μ_νρ have 4×4×4 = 64 components
    // With lower-index symmetry (ν↔ρ), only 4×10 = 40 unique components
    // Stored as 3D textures for O(1) lookup during geodesic integration
    //
    // Indexing: Γ[mu][packed_idx] where packed_idx = ν + ρ*(ρ+1)/2 for ν ≤ ρ
    //
    // Memory: 40 × sizeof(float) × nx × ny × nz
    //         For 128³ grid: 40 × 4 × 2M = 320 MB
    // =========================================================================
    
    // Christoffel textures: Gamma_mu_nu_rho (packed by lower symmetry)
    // mu=0 (t): Γ^t_00, Γ^t_01, Γ^t_02, Γ^t_03, Γ^t_11, Γ^t_12, Γ^t_13, Γ^t_22, Γ^t_23, Γ^t_33
    cudaTextureObject_t Gamma_t[10];
    // mu=1 (r/x): Γ^r_00, ...
    cudaTextureObject_t Gamma_r[10];
    // mu=2 (θ/y): Γ^θ_00, ...
    cudaTextureObject_t Gamma_theta[10];
    // mu=3 (φ/z): Γ^φ_00, ...
    cudaTextureObject_t Gamma_phi[10];
    
    bool christoffelLoaded;  // True if Christoffel textures are precomputed
    
    bool isLoaded;
};

//==============================================================================
// Camera State
//==============================================================================
struct CameraState {
    float3 position;   // Observer 4-position (spatial projection)
    float3 direction;  // Look direction (normalized)
    float3 up;         // Up vector (normalized)
    float3 right;      // Right vector (computed)
    
    float fov;         // Field of view in radians
    float yaw;         // Horizontal rotation
    float pitch;       // Vertical rotation
    float aspectRatio; // Width / Height
    
    // Observer's spherical coordinates for world-frame transforms
    // These are used to transform ray directions to world frame for correct
    // background texture sampling (fixes banding artifact)
    float observerTheta; // Observer's θ position (polar angle from +Y axis)
    float observerPhi;   // Observer's φ position (azimuthal angle in XZ plane)
    
    // Relativistic observer parameters
    float3 velocity;   // 3-velocity (coordinate)
    float gamma;       // Lorentz factor
};

//==============================================================================
// Numerical Metric Host Data (for upload)
//==============================================================================
struct NumericalMetricHostData {
    int3 dims;
    float3 origin;
    float3 spacing;
    
    // Arrays of size dims.x * dims.y * dims.z (Linear layout: x + y*nx + z*nx*ny)
    const float* gxx; const float* gxy; const float* gxz;
    const float* gyy; const float* gyz; const float* gzz;
    const float* alp;
    const float* betax; const float* betay; const float* betaz;
    
    // Hydrodynamics
    const float* rho;
    const float* vx; const float* vy; const float* vz;
};

//==============================================================================
// Integration Control
//==============================================================================
struct IntegrationParams {
    IntegratorType type;
    int maxSteps;          // Maximum integration steps
    float maxDistance;     // Maximum affine parameter
    float initialStepSize; // Initial step size
    float minStepSize;     // Minimum step (adaptive)
    float maxStepSize;     // Maximum step (adaptive)
    float tolerance;       // Error tolerance (adaptive)
    
    // Event horizon detection
    float horizonBuffer;   // Buffer distance from horizon
    bool detectHorizon;    // Enable horizon detection
};

//==============================================================================
// Accretion Disk Mode (Phase 7.5 - Planar vs Volumetric)
//==============================================================================
enum class DiskMode : int {
    Planar = 0,      // Infinitely thin disk (fast, bright, cinematic)
    Volumetric = 1   // Volumetric ray marching (slower, physically thicker)
};

//==============================================================================
// Temperature Model (Phase 8.6 - Novikov-Thorne)
//==============================================================================
enum class TemperatureModel : int {
    ShakuraSunyaev = 0,  // T(r) ∝ r^(-3/4), valid far from BH
    NovikovThorne = 1     // Relativistic: T→0 at ISCO (Page & Thorne 1974)
};

//==============================================================================
// Accretion Disk Parameters (Phase 4.2 Volumetric, Phase 7.5 Planar)
// Shakura-Sunyaev or Novikov-Thorne thin disk model
//==============================================================================
struct AccretionDiskParams {
    bool enabled;              // Enable disk rendering
    DiskMode diskMode;         // Planar (default) or Volumetric
    
    // Disk geometry
    float innerRadius;         // ISCO radius (units of M)
    float outerRadius;         // Outer edge of disk
    float heightScale;         // H/r ratio for disk thickness (volumetric only)
    
    // Temperature model selection (Phase 8.6)
    TemperatureModel temperatureModel;  // ShakuraSunyaev or NovikovThorne
    
    // Temperature parameters: T(r) = T_inner * (r_inner/r)^(3/4) * Q(r)
    // Q(r) = 1 for Shakura-Sunyaev, Q(r) → 0 at ISCO for Novikov-Thorne
    float innerTemperature;    // Temperature at inner edge (Kelvin)
    float temperatureExponent; // Typically 0.75 for thin disk
    
    // Emission/absorption (volumetric mode only)
    float emissionCoefficient; // j_ν scaling
    float absorptionCoefficient; // α_ν scaling
    float scatteringAlbedo;    // ω = σ_s / (σ_s + σ_a)
    
    // Ray marching (volumetric mode only)
    int volumetricSamples;     // Samples per ray through disk
    float opticalDepthMax;     // Stop if τ exceeds this
    
    // Visual enhancements
    float beamingExponent;     // Doppler beaming power (typically 3-4)
    float limbDarkening;       // Edge darkening factor
    bool useSpectralColors;    // Use PHSP001A spectral colors
};


//==============================================================================
// Lens Flare Parameters (Phase 6.5 - DNGR Cinematic Mode)
// Post-process PSF convolution for Interstellar-grade visuals
//==============================================================================
struct LensFlareParams {
    bool enabled;           // Toggle cinematic lens flare
    float intensity;        // Flare strength (0.0 - 1.0, default 0.3)
    float threshold;        // Brightness threshold for flare (default 0.8)
    int starPoints;         // Number of star points (default: 6)
};

//==============================================================================
// Bloom/Glow Parameters (Phase 7 - Cinematic Visual Quality)
// Screen-space bloom for soft, atmospheric disk edges like Interstellar
//==============================================================================
struct BloomParams {
    bool enabled;           // Toggle bloom effect
    float intensity;        // Bloom strength (0.0 - 1.0, default 0.5)
    float threshold;        // Brightness threshold (default 0.7)
    int blurPasses;         // Number of blur iterations (default 3)
    float blurRadius;       // Blur kernel radius in pixels (default 4.0)
};

//==============================================================================
// Ray Bundle State (Phase 6.3 - DNGR Eq A.18-A.27)
// Geodesic deviation for ray bundle anti-aliasing
//
// Tracks the central ray and deviation vectors that define how a bundle
// of initially parallel rays spreads due to spacetime curvature.
//==============================================================================
struct RayBundleState {
    // Central ray state (position and 4-velocity)
    float4 x;                 // Position (t, r, θ, φ)
    float4 u;                 // 4-velocity (dt/dλ, dr/dλ, dθ/dλ, dφ/dλ)
    
    // Deviation vectors (2 independent directions in screen plane)
    // These track how neighboring rays diverge from the central ray
    float4 xi1;               // ξ₁ - deviation in screen-x direction
    float4 dxi1;              // dξ₁/dλ - rate of change
    float4 xi2;               // ξ₂ - deviation in screen-y direction
    float4 dxi2;              // dξ₂/dλ - rate of change
    
    // Ellipse parameters in celestial sphere (computed from ξ₁, ξ₂)
    float semiMajor;          // a - semi-major axis (radians)
    float semiMinor;          // b - semi-minor axis (radians)
    float orientation;        // ψ - orientation angle (radians)
    
    // Integration state
    float affineParam;        // λ - affine parameter
    bool terminated;          // Ray has hit horizon or escaped
};

//==============================================================================
// Ray Bundle Integration Parameters (Phase 6.3)
//==============================================================================
struct RayBundleParams {
    bool enabled;             // Use ray bundles (vs single rays)
    bool useFullRiemann;      // Use full Riemann tensor (vs tidal approximation)
    float initialSize;        // Initial bundle size (pixel angular scale)
    float filterWidth;        // Gaussian filter sigma (in bundle radii)
    int filterSamples;        // Spatial filter samples for ellipse
    float minBundleSize;      // Stop subdivision below this size
    float maxEllipseRatio;    // Max a/b ratio before subdivision
};

//==============================================================================
// Path Tracing Parameters (with Temporal Accumulation - Phase 4.3)
//==============================================================================
struct PathTracingParams {
    int samplesPerPixel;   // SPP for Monte Carlo
    int maxBounces;        // Maximum ray bounces
    unsigned int seed;     // Random seed for this frame
    bool enableAccumulation; // Progressive rendering
    
    // === Temporal Accumulation (Phase 4.3) ===
    unsigned int accumulatedFrames;  // Number of frames accumulated
    float blendWeight;               // Weight for new frame (EMA: 0.05-0.2 typical)
    bool useExponentialMA;           // Use EMA (true) or simple average (false)
    bool resetAccumulation;          // Trigger to reset accumulation buffer
    float motionThreshold;           // Reset if camera motion exceeds this
};

//==============================================================================
// Prepass Configuration (Phase 2.3 Optimization)
//==============================================================================
enum class RayClassification : unsigned char {
    Unknown = 0,      // Not yet classified
    HitHorizon = 1,   // Ray fell into black hole
    Escaped = 2,      // Ray escaped to infinity
    HitDisk = 3,      // Ray hit accretion disk
    Complex = 4       // Requires full integration
};

struct PrepassParams {
    bool enabled;                 // Enable prepass optimization
    int downsampleFactor;         // 2 = half res, 4 = quarter res
    RayClassification* classificationBuffer; // Low-res classification results
    float3* directionBuffer;      // Final ray directions from prepass
    int2 prepassDimensions;       // Prepass buffer size
};

//==============================================================================
// Geodesic Lookup Table (Phase 3.2 Optimization)
// Precomputed ray trajectories for O(1) lookup instead of O(N) integration
//==============================================================================

/// Entry in geodesic lookup table
struct GeodesicLUTEntry {
    float3 finalDirection;   // Final ray direction when escaping
    float deflectionAngle;   // Total bending angle
    float timeDelay;         // Shapiro time delay
    RayClassification result; // Hit horizon, escaped, etc.
    unsigned char nOrbits;   // Number of half-orbits around BH
    unsigned char padding[2];
};

/// Configuration for geodesic LUT generation
struct GeodesicLUTConfig {
    int impactParameterBins;     // Resolution in b dimension (default 512)
    int phiBins;                 // Resolution in initial angle (default 256)
    float impactParamMin;        // Minimum impact parameter (units of M)
    float impactParamMax;        // Maximum impact parameter (typically 10M)
    float observerDistance;      // r_observer for LUT computation
    bool includeHigherOrderImages; // Include n>0 images
};

/// Geodesic lookup table params for GPU
struct GeodesicLUTParams {
    bool enabled;                 // Use LUT instead of integration
    cudaTextureObject_t lutTexture; // 2D texture of GeodesicLUTEntry
    int2 lutDimensions;           // (impactBins, phiBins)
    float impactParamMin;
    float impactParamMax;
    float observerDistance;
};

//==============================================================================
// OptiX Launch Parameters
// This structure is uploaded to GPU constant memory before each launch
//==============================================================================
struct LaunchParams {
    // === Output ===
    float4* frameBuffer;          // RGBA framebuffer (linear)
    float4* accumBuffer;          // Accumulation buffer (progressive)
    int2 frameDimensions;         // (width, height)

    unsigned int frameIndex;      // Current frame number
    float time;                   // Simulation time (for animation)
    
    // === Camera ===
    CameraState camera;
    
    // === Metric Configuration ===
    MetricType metricType;
    MetricParams metricParams;
    NumericalMetricData numericalMetric;
    
    // === Integration ===
    IntegrationParams integration;
    
    // === Path Tracing ===
    PathTracingParams pathTracing;
    
    // === Background ===
    cudaTextureObject_t backgroundTexture;
    bool useBackgroundTexture;
    float3 backgroundColor;       // Fallback solid color
    
    // === Prepass (Phase 2.3) ===
    PrepassParams prepass;
    
    // === Geodesic LUT (Phase 3.2) ===
    GeodesicLUTParams geodesicLUT;
    
    // === Accretion Disk (Phase 4.2) ===
    AccretionDiskParams accretionDisk;
    
    // === Lens Flare (Phase 6.5 - DNGR Cinematic Mode) ===
    LensFlareParams lensFlare;
    
    // === Bloom/Glow (Phase 7 - Cinematic Visual Quality) ===
    BloomParams bloom;
    
    // === Ray Bundles (Phase 6.3 - DNGR Anti-Aliasing) ===
    RayBundleParams rayBundle;
    
    // === OptiX Handles (not typically used for geodesic raymarching) ===
    OptixTraversableHandle traversable;
    
    // === Debug ===
    bool debugMode;
    int debugPixelX;
    int debugPixelY;
};

//==============================================================================
// SBT Record Types
// These wrap per-program data passed through OptiX Shader Binding Table
//==============================================================================

// Header for all SBT records (required by OptiX)
template<typename T>
struct SbtRecord {
    __align__(OPTIX_SBT_RECORD_ALIGNMENT) char header[OPTIX_SBT_RECORD_HEADER_SIZE];
    T data;
};

// Ray generation program data
struct RayGenData {
    // Currently empty - camera data comes from LaunchParams
};

// Miss program data
struct MissData {
    float3 backgroundColor;
};

// Hit group program data (for future geometry intersection)
struct HitGroupData {
    float3 albedo;           // Surface color
    float roughness;         // PBR roughness
    float metallic;          // PBR metallic
};

// Empty record for programs with no per-instance data
struct EmptyData {};

//==============================================================================
// Type aliases for SBT records
//==============================================================================
using RayGenRecord = SbtRecord<RayGenData>;
using MissRecord = SbtRecord<MissData>;
using HitGroupRecord = SbtRecord<HitGroupData>;
using EmptyRecord = SbtRecord<EmptyData>;

//==============================================================================
// Helper Functions (Host-side)
//==============================================================================
#ifndef __CUDACC__

inline float3 make_float3(float x, float y, float z) {
    return float3{x, y, z};
}

inline float4 make_float4(float x, float y, float z, float w) {
    return float4{x, y, z, w};
}

inline int2 make_int2(int x, int y) {
    return int2{x, y};
}

inline int3 make_int3(int x, int y, int z) {
    return int3{x, y, z};
}

// Initialize default launch parameters
inline LaunchParams createDefaultLaunchParams() {
    LaunchParams params = {};
    
    // Camera defaults
    params.camera.position = make_float3(0.0f, 0.0f, -10.0f);
    params.camera.direction = make_float3(0.0f, 0.0f, 1.0f);
    params.camera.up = make_float3(0.0f, 1.0f, 0.0f);
    params.camera.right = make_float3(1.0f, 0.0f, 0.0f);
    params.camera.fov = 1.0472f;  // 60 degrees in radians
    params.camera.aspectRatio = 16.0f / 9.0f;
    params.camera.gamma = 1.0f;
    
    // Metric defaults (Schwarzschild via Kerr-Schild family)
    params.metricType = MetricType::Schwarzschild;  // Legacy compatibility
    params.metricParams.family = MetricFamily::KerrSchild;
    params.metricParams.kerrSchild.M = 1.0f;
    params.metricParams.kerrSchild.a = 0.0f;
    params.metricParams.kerrSchild.Q = 0.0f;
    params.metricParams.kerrSchild.Lambda = 0.0f;
    // Sync legacy fields
    params.metricParams.M = 1.0f;
    params.metricParams.a = 0.0f;
    params.metricParams.Q = 0.0f;
    params.metricParams.Lambda = 0.0f;
    params.metricParams.rs = 2.0f;  // 2M
    
    // Integration defaults
    params.integration.type = IntegratorType::Symplectic;
    params.integration.maxSteps = 10000;         // Increased for fine volumetric sampling
    params.integration.maxDistance = 1000.0f;    // Large celestial sphere for camera orbit
    params.integration.initialStepSize = 0.05f;
    params.integration.minStepSize = 0.0001f;    // Allow very fine steps (1e-4)
    params.integration.maxStepSize = 0.5f;
    params.integration.tolerance = 1e-6f;
    params.integration.horizonBuffer = 0.01f;
    params.integration.detectHorizon = true;
    
    // Path tracing defaults
    params.pathTracing.samplesPerPixel = 1;
    params.pathTracing.maxBounces = 1;
    params.pathTracing.enableAccumulation = false;
    params.pathTracing.accumulatedFrames = 0;
    params.pathTracing.blendWeight = 0.1f;  // EMA weight
    params.pathTracing.useExponentialMA = true;
    params.pathTracing.resetAccumulation = false;
    params.pathTracing.motionThreshold = 0.01f;
    
    // Background
    params.useBackgroundTexture = false;
    params.backgroundColor = make_float3(0.02f, 0.02f, 0.08f);  // Visible dark blue background
    
    // Prepass defaults (Phase 2.3)
    params.prepass.enabled = false;  // Disabled by default
    params.prepass.downsampleFactor = 4;
    params.prepass.classificationBuffer = nullptr;
    params.prepass.directionBuffer = nullptr;
    params.prepass.prepassDimensions = make_int2(0, 0);
    
    // Geodesic LUT defaults (Phase 3.2)
    params.geodesicLUT.enabled = false;  // Disabled by default
    params.geodesicLUT.lutTexture = 0;
    params.geodesicLUT.lutDimensions = make_int2(0, 0);
    params.geodesicLUT.impactParamMin = 0.0f;
    params.geodesicLUT.impactParamMax = 10.0f;
    params.geodesicLUT.observerDistance = 100.0f;
    
    // Accretion disk defaults (Phase 4.2)
    // For visible disk: emission >> absorption gives bright colors
    // heightScale controls thickness: H/r ratio (thin disk ~0.01)
    // Accretion Disk defaults - Tuned for "Interstellar" (Gargantua) look
    params.accretionDisk.enabled = true;
    params.accretionDisk.diskMode = Sirius::DiskMode::Planar;  // Planar = fast, bright; Volumetric = slower
    params.accretionDisk.innerRadius = 0.0f;     // 0 = use ISCO (auto-calculated as 6M)
    params.accretionDisk.outerRadius = 15.0f;    // In units of M, so 15M (typical ~10-30M)
    params.accretionDisk.heightScale = 0.005f;   // Razor thin
    params.accretionDisk.temperatureModel = Sirius::TemperatureModel::ShakuraSunyaev;  // Cinematic: bright inner edge
    params.accretionDisk.innerTemperature = 15000.0f; // Hot inner edge → white/blue-white
    params.accretionDisk.temperatureExponent = 0.75f;  // Standard T ∝ r^(-3/4)
    params.accretionDisk.emissionCoefficient = 10.0f;  // High vibrancy for outer glow
    params.accretionDisk.absorptionCoefficient = 2.0f;  // Translucent glow
    params.accretionDisk.scatteringAlbedo = 0.0f;
    params.accretionDisk.volumetricSamples = 8;  // Reduced for performance; bloom compensates visually
    params.accretionDisk.opticalDepthMax = 2.0f;   // Lower to prevent pure black occlusion
    params.accretionDisk.beamingExponent = 2.0f; // Softer Doppler transitions
    params.accretionDisk.limbDarkening = 0.6f;
    params.accretionDisk.useSpectralColors = true;
    
    // Lens flare defaults (Phase 6.5 - DNGR Cinematic Mode)
    // Disabled by default to preserve physics accuracy
    params.lensFlare.enabled = false;
    params.lensFlare.intensity = 0.3f;
    params.lensFlare.threshold = 0.8f;
    params.lensFlare.starPoints = 6;
    
    // Bloom defaults (Phase 7 - Cinematic Visual Quality)
    // Enabled by default for soft, atmospheric disk edges
    params.bloom.enabled = true;
    params.bloom.intensity = 0.5f;
    params.bloom.threshold = 0.7f;
    params.bloom.blurPasses = 3;
    params.bloom.blurRadius = 4.0f;
    
    // Ray bundle defaults (Phase 6.3 - DNGR Anti-Aliasing)
    // Enabled by default - DNGR methodology for anti-aliasing in curved spacetime
    params.rayBundle.enabled = true;
    params.rayBundle.useFullRiemann = false;  // Start with tidal approx (faster), toggle for quality
    params.rayBundle.initialSize = 0.01f;     // Larger for more visible AA effect
    params.rayBundle.filterWidth = 1.5f;      // Gaussian sigma in bundle radii
    params.rayBundle.filterSamples = 4;       // Samples for ellipse filtering
    params.rayBundle.minBundleSize = 1e-6f;   // Minimum size before termination
    params.rayBundle.maxEllipseRatio = 100.0f; // Max elongation
    
    return params;
}

#endif // !__CUDACC__

} // namespace Sirius
