// RDOP003A.h - OptiX Launch Parameters
//
// Shared Host/Device structures for GPU kernels.
// Defines metrics, camera, integration, accretion disk, and post-processing params.

#pragma once

#include <cstdint>

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
//
// PRECISION NOTE (specification.md §5.2-5.3):
// GPU (single precision): Null constraint tolerance 10^-5
//   - ~100 ops/step × 500 steps = 50,000 ops total
//   - Error bound: ~10^-2 accumulated, 10^-5 per step acceptable
//
// CPU (double precision): Null constraint tolerance 10^-10
//   - Same operations but ~10^-14 per step
//   - Error bound: ~10^-11 accumulated
//
// The tolerance difference is intentional and acceptable for visualization.
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
// Conservation Diagnostics (Optional Debug Mode)
//
// Tracks conservation law violations during geodesic integration.
// Specification requirements (specification.md):
//   - Null constraint: |g_μν k^μ k^ν| < 10^-5 (GPU) / 10^-10 (CPU)
//   - Energy E drift: |ΔE/E| < 10^-4
//   - Angular momentum L_z drift: |ΔL_z/L_z| < 10^-4
//
// Enable via LaunchParams for debugging; disabled by default for performance.
//==============================================================================
struct ConservationDiagnostics {
    float max_null_violation;   // Maximum |g_μν k^μ k^ν| observed
    float max_E_drift;          // Maximum |ΔE/E| relative energy drift
    float max_L_drift;          // Maximum |ΔL_z/L_z| relative angular momentum drift
    int rays_with_violations;   // Count of rays exceeding tolerance
    bool enabled;               // Enable diagnostic tracking (false by default)
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
// SMBH Parameters (GPU-compatible)
// Supermassive black hole astrophysical scaling
//==============================================================================
struct SMBHParamsGPU {
    float mass_M;                  // Mass in geometric units (M = GM/c² = 1)
    float spin;                    // Dimensionless spin a* = a/M
    float inclination_rad;         // Observer inclination [radians]
    float r_isco;                  // ISCO radius in M
    float r_horizon;               // Outer horizon radius in M
    float r_inner_horizon;         // Inner horizon radius in M
    float angular_size_rg;         // Angular size of r_g [radians]
    float padding;                 // Alignment padding
};

//==============================================================================
// Turbulence Parameters (GPU-compatible)
// Kolmogorov cascade density perturbations
//==============================================================================
struct TurbulenceParamsGPU {
    float amplitude;
    float outer_scale;
    float inner_scale;
    float lacunarity;
    float persistence;
    uint32_t octaves;
    uint32_t seed;
    uint32_t enabled;
};

//==============================================================================
// Photon Ring Parameters (GPU-compatible)
// Enhancement for light orbiting black hole before escaping
//==============================================================================
struct PhotonRingParams {
    uint32_t enabled;              // Enable photon ring enhancement
    uint32_t minOrbits;            // Minimum orbits to qualify (typically 1)
    float brightnessBoost;         // Brightness multiplier per orbit (e.g., 2.0)
    float falloffPerOrbit;         // Brightness decay per additional orbit (e.g., 0.5)
    float innerSharpness;          // Edge sharpness near shadow (higher = sharper)
    float colorShift;              // Blue shift for highly deflected photons
    float ringWidth;               // Apparent width factor (0.01 = thin, 0.1 = thick)
    float padding;
};

//==============================================================================
// Corona Parameters (GPU-compatible)
// Inverse-Compton scattering corona model
//==============================================================================
struct CoronaParamsGPU {
    float temperature_keV;
    float optical_depth;
    float scale_height;
    float inner_radius;
    float outer_radius;
    float emissivity_index;
    float intensity_scale;
    uint32_t geometry;
    float lamppost_height;
    float comptonization_y;
    float spectral_index;
    uint32_t enabled;
};

//==============================================================================
// Volumetric Disk Parameters (GPU-compatible)
// Enhanced Shakura-Sunyaev + turbulence + corona
//==============================================================================
struct VolumetricDiskParamsGPU {
    // Structure
    float inner_radius;
    float outer_radius;
    float reference_radius;
    float H_over_r;
    float H_power;
    float density_power;
    float padding1[2];

    // Temperature
    float inner_temperature;
    float temperature_power;
    uint32_t use_novikov_thorne;
    float padding2;

    // Emission
    float emission_coeff;
    float absorption_coeff;
    float scattering_albedo;
    float beaming_exponent;

    // Ray marching
    int volumetric_samples;
    float optical_depth_max;
    uint32_t enabled;
    float padding3;

    // Sub-models
    TurbulenceParamsGPU turbulence;
    CoronaParamsGPU corona;
};

//==============================================================================
// MHD Jet Parameters (GPU-compatible)
// Relativistic jet with magnetic field structure
//==============================================================================
struct JetMHDParamsGPU {
    // Magnetic field
    float B_base;
    float power_law_index;
    float B_field_order;
    float padding1;

    // Geometry
    float opening_angle;
    float z_launch;
    float z_max;
    float collimation;

    // Kinematics
    float lorentz_factor;
    float beta;
    float velocity_profile;
    float padding2;

    // Electron distribution
    float electron_index;
    float spectral_index;
    float gamma_min;
    float gamma_max;
    float n_e_0;
    float n_e_decay;
    float intensity_scale;
    float max_polarisation;

    // Flags
    uint32_t enabled;
    uint32_t enable_polarisation;
    uint32_t padding3[2];
};

//==============================================================================
// Starfield Parameters (GPU-compatible)
// Depth-resolved stellar catalog with parallax
//==============================================================================
struct StarfieldParamsGPU {
    uint64_t star_buffer;          // CUdeviceptr to StarEntry array
    uint32_t star_count;
    float magnitude_limit;
    float brightness_scale;
    float aperture_mm;
    float focus_distance_pc;
    uint32_t enabled;
    uint32_t enable_parallax;
    uint32_t enable_dof;
};

//==============================================================================
// Film Simulation Parameters (GPU-compatible)
// IMAX 70mm post-processing
//==============================================================================
struct FilmParamsGPU {
    // Grain
    float grain_intensity;
    float grain_size;
    float grain_uniformity;
    uint32_t grain_seed;

    // Halation
    float halation_radius;
    float halation_strength;
    float halation_threshold;
    float halation_color_r;
    float halation_color_g;
    float halation_color_b;
    float padding1;
    float padding2;

    // Color
    float saturation;
    float contrast;
    float exposure;
    float toe_strength;
    float shoulder_strength;
    float midtone_point;
    float padding3;
    float padding4;

    // Vignette
    float vignette_strength;
    float vignette_radius;
    float vignette_softness;
    float padding5;

    // Chromatic Aberration (Phase 9 - Cinematic Lens Effect)
    float chromatic_strength;      // RGB channel separation (0.0 = off, 0.01 = subtle)
    float chromatic_radial_power;  // How much effect increases toward edges (2.0 typical)
    float padding6;
    float padding7;

    // Feature flags (packed as bits)
    // bit 0: grain, bit 1: halation, bit 2: vignette, bit 3: enabled, bit 4: chromatic
    uint32_t features;
};

//==============================================================================
// LaunchParamsCore - Essential parameters (~2KB)
// This minimal structure fits comfortably in constant memory and contains
// everything needed for basic geodesic ray tracing without advanced features.
//==============================================================================
struct LaunchParamsCore {
    // === Output ===
    float4* frameBuffer;          // RGBA framebuffer (linear)
    int2 frameDimensions;         // (width, height)
    unsigned int frameIndex;      // Current frame number

    // === Camera ===
    CameraState camera;

    // === Metric Configuration ===
    MetricType metricType;
    MetricParams metricParams;

    // === Integration ===
    IntegrationParams integration;

    // === Background ===
    cudaTextureObject_t backgroundTexture;
    bool useBackgroundTexture;
    float3 backgroundColor;

    // === Debug ===
    bool debugMode;
    int debugPixelX;
    int debugPixelY;
};

//==============================================================================
// LaunchParamsExtended - Optional/advanced features (global memory)
// Access via global memory pointer from LaunchParams.extended
// Contains cinematic features, volumetric rendering, and optimization caches.
//==============================================================================
struct LaunchParamsExtended {
    // === Path Tracing ===
    PathTracingParams pathTracing;

    // === Prepass (Phase 2.3) ===
    PrepassParams prepass;

    // === Geodesic LUT (Phase 3.2) ===
    GeodesicLUTParams geodesicLUT;

    // === Accretion Disk (Phase 4.2) ===
    AccretionDiskParams accretionDisk;

    // === Lens Flare (Phase 6.5) ===
    LensFlareParams lensFlare;

    // === Bloom/Glow (Phase 7) ===
    BloomParams bloom;

    // === Ray Bundles (Phase 6.3) ===
    RayBundleParams rayBundle;

    // === Numerical Metric (if used) ===
    NumericalMetricData numericalMetric;

    // === Cinematic Expansion (2026) ===
    SMBHParamsGPU smbh;
    VolumetricDiskParamsGPU volumetricDisk;
    JetMHDParamsGPU jetMHD;
    StarfieldParamsGPU starfield;
    FilmParamsGPU film;

    // === Accumulation Buffer ===
    float4* accumBuffer;
};

//==============================================================================
// OptiX Launch Parameters (Full)
// This structure is uploaded to GPU constant memory before each launch.
// For performance-critical paths, use LaunchParamsCore only.
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

    // === Photon Ring Enhancement (Phase 9 - Cinematic Realism) ===
    PhotonRingParams photonRing;

    // === Cinematic Expansion (2026) ===
    // SMBH astrophysical parameters
    SMBHParamsGPU smbh;

    // Volumetric disk with turbulence and corona
    VolumetricDiskParamsGPU volumetricDisk;

    // MHD relativistic jet
    JetMHDParamsGPU jetMHD;

    // Depth-resolved starfield
    StarfieldParamsGPU starfield;

    // IMAX 70mm film simulation
    FilmParamsGPU film;

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
    params.accretionDisk.outerRadius = 30.0f;    // Extended for cinematic wide composition
    params.accretionDisk.heightScale = 0.005f;   // Razor thin
    params.accretionDisk.temperatureModel = Sirius::TemperatureModel::ShakuraSunyaev;  // Cinematic: bright inner edge
    params.accretionDisk.innerTemperature = 6500.0f; // Warm orange inner edge (Interstellar look)
    params.accretionDisk.temperatureExponent = 0.75f;  // Standard T ∝ r^(-3/4)
    params.accretionDisk.emissionCoefficient = 2.5f;  // Balanced emission
    params.accretionDisk.absorptionCoefficient = 2.0f;  // Translucent glow
    params.accretionDisk.scatteringAlbedo = 0.0f;
    params.accretionDisk.volumetricSamples = 8;  // Reduced for performance; bloom compensates visually
    params.accretionDisk.opticalDepthMax = 2.0f;   // Lower to prevent pure black occlusion
    params.accretionDisk.beamingExponent = 3.5f; // Strong Doppler beaming for asymmetry
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

    // =========================================================================
    // Cinematic Expansion Defaults (2026)
    // =========================================================================

    // SMBH defaults (normalized geometric units)
    params.smbh.mass_M = 1.0f;
    params.smbh.spin = 0.0f;
    params.smbh.inclination_rad = 1.5708f;  // 90 degrees
    params.smbh.r_isco = 6.0f;              // Schwarzschild ISCO
    params.smbh.r_horizon = 2.0f;           // Schwarzschild horizon
    params.smbh.r_inner_horizon = 0.0f;
    params.smbh.angular_size_rg = 1e-6f;
    params.smbh.padding = 0.0f;

    // Volumetric disk defaults
    params.volumetricDisk.inner_radius = 0.0f;  // Use ISCO
    params.volumetricDisk.outer_radius = 500.0f;
    params.volumetricDisk.reference_radius = 10.0f;
    params.volumetricDisk.H_over_r = 0.1f;
    params.volumetricDisk.H_power = 0.25f;
    params.volumetricDisk.density_power = 1.5f;
    params.volumetricDisk.padding1[0] = params.volumetricDisk.padding1[1] = 0.0f;
    params.volumetricDisk.inner_temperature = 1e7f;
    params.volumetricDisk.temperature_power = 0.75f;
    params.volumetricDisk.use_novikov_thorne = 1;
    params.volumetricDisk.padding2 = 0.0f;
    params.volumetricDisk.emission_coeff = 1.0f;
    params.volumetricDisk.absorption_coeff = 0.1f;
    params.volumetricDisk.scattering_albedo = 0.3f;
    params.volumetricDisk.beaming_exponent = 3.0f;
    params.volumetricDisk.volumetric_samples = 64;
    params.volumetricDisk.optical_depth_max = 10.0f;
    params.volumetricDisk.enabled = 0;  // Disabled by default (use planar disk)
    params.volumetricDisk.padding3 = 0.0f;

    // Turbulence sub-model
    params.volumetricDisk.turbulence.amplitude = 0.3f;
    params.volumetricDisk.turbulence.outer_scale = 5.0f;
    params.volumetricDisk.turbulence.inner_scale = 0.1f;
    params.volumetricDisk.turbulence.lacunarity = 2.0f;
    params.volumetricDisk.turbulence.persistence = 0.5f;
    params.volumetricDisk.turbulence.octaves = 6;
    params.volumetricDisk.turbulence.seed = 12345;
    params.volumetricDisk.turbulence.enabled = 1;

    // Corona sub-model
    params.volumetricDisk.corona.temperature_keV = 100.0f;
    params.volumetricDisk.corona.optical_depth = 0.5f;
    params.volumetricDisk.corona.scale_height = 5.0f;
    params.volumetricDisk.corona.inner_radius = 0.0f;
    params.volumetricDisk.corona.outer_radius = 20.0f;
    params.volumetricDisk.corona.emissivity_index = 3.0f;
    params.volumetricDisk.corona.intensity_scale = 1.0f;
    params.volumetricDisk.corona.geometry = 3;  // Extended
    params.volumetricDisk.corona.lamppost_height = 10.0f;
    params.volumetricDisk.corona.comptonization_y = 0.8f;
    params.volumetricDisk.corona.spectral_index = 1.7f;
    params.volumetricDisk.corona.enabled = 0;

    // MHD Jet defaults
    params.jetMHD.B_base = 1e4f;
    params.jetMHD.power_law_index = 1.0f;
    params.jetMHD.B_field_order = 0.5f;
    params.jetMHD.padding1 = 0.0f;
    params.jetMHD.opening_angle = 0.1f;
    params.jetMHD.z_launch = 3.0f;
    params.jetMHD.z_max = 200.0f;
    params.jetMHD.collimation = 0.5f;
    params.jetMHD.lorentz_factor = 5.0f;
    params.jetMHD.beta = 0.98f;
    params.jetMHD.velocity_profile = 0.0f;
    params.jetMHD.padding2 = 0.0f;
    params.jetMHD.electron_index = 2.2f;
    params.jetMHD.spectral_index = 0.6f;
    params.jetMHD.gamma_min = 10.0f;
    params.jetMHD.gamma_max = 1e6f;
    params.jetMHD.n_e_0 = 1e5f;
    params.jetMHD.n_e_decay = 2.0f;
    params.jetMHD.intensity_scale = 1.0f;
    params.jetMHD.max_polarisation = 0.7f;
    params.jetMHD.enabled = 0;  // Disabled by default
    params.jetMHD.enable_polarisation = 1;
    params.jetMHD.padding3[0] = params.jetMHD.padding3[1] = 0;

    // Starfield defaults
    params.starfield.star_buffer = 0;  // No buffer allocated
    params.starfield.star_count = 0;
    params.starfield.magnitude_limit = 12.0f;
    params.starfield.brightness_scale = 1.0f;
    params.starfield.aperture_mm = 50.0f;
    params.starfield.focus_distance_pc = 100.0f;
    params.starfield.enabled = 0;  // Disabled by default
    params.starfield.enable_parallax = 1;
    params.starfield.enable_dof = 1;

    // Film simulation defaults (Interstellar-style IMAX)
    params.film.grain_intensity = 0.025f;
    params.film.grain_size = 1.5f;
    params.film.grain_uniformity = 0.7f;
    params.film.grain_seed = 0;
    params.film.halation_radius = 8.0f;
    params.film.halation_strength = 0.15f;
    params.film.halation_threshold = 0.8f;
    params.film.halation_color_r = 1.0f;
    params.film.halation_color_g = 0.5f;
    params.film.halation_color_b = 0.2f;
    params.film.padding1 = params.film.padding2 = 0.0f;
    params.film.saturation = 0.95f;
    params.film.contrast = 1.05f;
    params.film.exposure = 1.0f;  // Neutral exposure (1.0 = no change)
    params.film.toe_strength = 0.5f;
    params.film.shoulder_strength = 0.5f;
    params.film.midtone_point = 0.18f;
    params.film.padding3 = params.film.padding4 = 0.0f;
    params.film.vignette_strength = 0.3f;
    params.film.vignette_radius = 1.2f;
    params.film.vignette_softness = 0.5f;
    params.film.padding5 = 0.0f;
    params.film.features = 0;  // All disabled by default for GPU rendering

    return params;
}

#endif // !__CUDACC__

} // namespace Sirius
