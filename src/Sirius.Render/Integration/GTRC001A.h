// =============================================================================
// GTRC001A.h - Geodesic Ray Tracer Component
// Component ID: GTRC001A (Render/Integration/GeodesicTracer)
// =============================================================================
//
// PURPOSE
// =======
// Wraps Sirius.Core geodesic integrator for batch rendering applications.
// Bridges the gap between camera ray generation and physics-based integration.
//
// MATHEMATICAL BASIS
// ==================
// Traces null geodesics: d²x^μ/dλ² + Γ^μ_αβ (dx^α/dλ)(dx^β/dλ) = 0
// Constraint: g_μν k^μ k^ν = 0 (null condition)
//
// COORDINATE SYSTEMS
// ==================
// - Camera generates rays in Boyer-Lindquist coordinates (t, r, θ, φ)
// - Integration uses Kerr-Schild Cartesian (t, x, y, z)
// - This class handles the coordinate transformation
//
// TESTS: TSIN005A.cpp
// =============================================================================

#ifndef GTRC001A_H
#define GTRC001A_H

#include "PHGD001A.h"
#include "PHMT100A.h"
#include "CMBS001A.h"
#include "MTTN001A.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Sirius {

// =============================================================================
// Trace Result Structure
// =============================================================================
struct TraceResult {
    /// @brief Outcome of ray tracing
    enum class Outcome {
        ESCAPED,     ///< Ray escaped to infinity (sample background)
        HORIZON,     ///< Ray captured by black hole
        DISK_HIT,    ///< Ray intersected accretion disk
        MAX_STEPS,   ///< Integration terminated at step limit
        SPIRALING    ///< Ray detected as spiraling near photon sphere (early termination)
    };

    Outcome outcome = Outcome::MAX_STEPS;

    // Final state (valid for ESCAPED outcome)
    Vec4 final_position;    ///< Position when ray terminated
    Vec4 final_direction;   ///< Direction when ray terminated

    // Disk hit data (valid for DISK_HIT outcome)
    float disk_radius = 0.0f;       ///< Radius of disk intersection (in M)
    float disk_temperature = 0.0f;  ///< Temperature at intersection point
    float disk_phi = 0.0f;          ///< Azimuthal angle of intersection

    // Motion blur data (valid for DISK_HIT outcome)
    // Enables analytical computation of g-factor at different azimuthal angles
    // g(φ + δφ) = grav_factor / (gamma × (1 - v_orb × (A·cos(δφ) + B·sin(δφ))))
    // where A = n_x·sin(φ) - n_y·cos(φ), B = n_x·cos(φ) + n_y·sin(φ)
    float gfactor_grav = 1.0f;      ///< Gravitational redshift: sqrt(1 - 2Mr/(r² + a²))
    float gfactor_gamma = 1.0f;     ///< Lorentz factor: 1/sqrt(1 - v²)
    float gfactor_v_orb = 0.0f;     ///< Orbital velocity magnitude: r × Ω
    float gfactor_A = 0.0f;         ///< Doppler coefficient A = n_x·sin(φ) - n_y·cos(φ)
    float gfactor_B = 0.0f;         ///< Doppler coefficient B = n_x·cos(φ) + n_y·sin(φ)

    // Diagnostics
    int steps_taken = 0;            ///< Number of integration steps
    float redshift = 1.0f;          ///< Combined g-factor at original φ
    bool numerical_failure = false; ///< True if NaN/Inf detected

    // Phase 1: Magnification and photon ring
    float magnification = 1.0f;     ///< |det(J)| beam magnification from gravitational lensing
    float min_radius = 1e10f;       ///< Minimum radius reached during integration
    bool photon_ring = false;       ///< True if ray passed through photon ring region

    // ==========================================================================
    // Higher-Order Imaging (Photon Ring Structure)
    // ==========================================================================
    // Tracks winding number and multiple disk crossings for n=0,1,2,3... images
    // See: Gralla, Lupsasca & Marrone (2020) "The shape of the black hole photon ring"
    //
    // Image order classification:
    //   n=0: Direct (weak deflection, b >> b_crit)
    //   n=1: Primary (half-orbit, b slightly > b_crit)
    //   n=2: Secondary (full orbit, one winding)
    //   n=3+: Higher orders (multiple windings)
    //
    // Critical impact parameter: b_crit = 3√3 M ≈ 5.196 M (Schwarzschild)
    // ==========================================================================

    /// @brief Number of equatorial plane crossings (proxy for winding number)
    /// Each crossing adds ~0.5 to the effective image order
    int equatorial_crossings = 0;

    /// @brief Total azimuthal angle change |Δφ| (radians)
    /// Full orbit = 2π, half orbit = π
    float total_phi_change = 0.0f;

    /// @brief Image order n = floor(equatorial_crossings / 2)
    /// n=0: direct, n=1: primary, n=2: secondary, etc.
    int image_order = 0;

    /// @brief Impact parameter b = L/E (angular momentum / energy)
    /// Compare to b_crit to classify ray trajectory
    float impact_parameter = 0.0f;

    // ==========================================================================
    // Multiple Disk Crossings (for accumulating all image orders)
    // ==========================================================================
    static constexpr int MAX_DISK_CROSSINGS = 4;  ///< Track up to n=3 images

    /// @brief Disk crossing data for each intersection
    struct DiskCrossing {
        float r = 0.0f;         ///< Radius at crossing
        float phi = 0.0f;       ///< Azimuthal angle at crossing
        float temperature = 0.0f;  ///< Temperature at crossing
        float redshift = 1.0f;  ///< g-factor at crossing
        int crossing_index = 0; ///< Which crossing (0 = first, 1 = second, etc.)
        bool valid = false;     ///< True if this crossing occurred
    };

    /// @brief All disk crossings (primary, secondary, tertiary, etc.)
    DiskCrossing disk_crossings[MAX_DISK_CROSSINGS] = {};

    /// @brief Number of valid disk crossings
    int num_disk_crossings = 0;

    // ==========================================================================
    // Volumetric Disk Data (Phase 6)
    // ==========================================================================
    // For 3D disks with finite vertical structure, ray marching accumulates
    // emission along the path through the disk volume.
    //
    // Radiative transfer: I_out = I_in × e^(-τ) + ∫ S × e^(-(τ-τ')) dτ'
    // ==========================================================================

    /// @brief Accumulated volumetric emission [R, G, B] from ray marching
    float volumetric_emission[3] = {0.0f, 0.0f, 0.0f};

    /// @brief Total optical depth along ray path through disk
    float optical_depth = 0.0f;

    /// @brief Path length through disk volume [GM/c²]
    float volumetric_path_length = 0.0f;

    /// @brief True if ray passed through volumetric disk
    bool volumetric_hit = false;
};

// =============================================================================
// Tracer Configuration
// =============================================================================
struct TracerConfig {
    // Termination conditions
    float escape_radius = 100.0f;   ///< Radius beyond which ray is considered escaped
    float horizon_factor = 1.05f;   ///< Safety factor above Schwarzschild radius
    int max_steps = 5000;           ///< Maximum integration steps per ray

    // Accretion disk parameters
    bool enable_disk = true;        ///< Whether to check for disk intersection
    float disk_inner = 6.0f;        ///< Inner disk radius (ISCO for Schwarzschild = 6M)
    float disk_outer = 20.0f;       ///< Outer disk radius
    float disk_thickness = 0.01f;   ///< Angular half-thickness (radians from equator)
    float disk_temperature_inner = 1.0f;  ///< Temperature at inner edge (arbitrary units)

    // Volumetric disk parameters (Phase 6)
    bool enable_volumetric = false;     ///< Enable 3D volumetric disk (vs thin disk)
    float volumetric_H_over_r = 0.1f;   ///< Scale height ratio H/r at reference radius
    float volumetric_H_power = 0.25f;   ///< Flaring index: H/r ∝ r^H_power
    float volumetric_tau_midplane = 10.0f;  ///< Midplane optical depth at r_ref
    int volumetric_samples = 32;        ///< Number of ray marching samples
    float volumetric_tau_max = 10.0f;   ///< Early termination optical depth threshold

    // Integration parameters
    IntegratorConfig integrator;    ///< RK45 configuration

    // Default constructor with sensible defaults
    TracerConfig() {
        integrator = Geodesic::getDefaultConfig();
        integrator.abs_tolerance = 1e-6f;
        integrator.rel_tolerance = 1e-6f;
        integrator.min_step = 1e-6f;
        integrator.max_step = 0.1f;
        integrator.initial_step = 0.01f;
        integrator.use_rk45 = true;
    }
};

// =============================================================================
// GeodesicTracer Class
// =============================================================================
class GeodesicTracer {
public:
    /// @brief Construct tracer with metric and configuration
    /// @param metric Pointer to metric (must remain valid for tracer lifetime)
    /// @param config Tracer configuration
    GeodesicTracer(IMetric* metric, const TracerConfig& config);

    /// @brief Default destructor
    ~GeodesicTracer() = default;

    // =========================================================================
    // Main Interface
    // =========================================================================

    /// @brief Trace a camera ray through spacetime
    /// @param camera_ray Ray from camera (origin in BL coords, direction in BL basis)
    /// @return TraceResult with outcome and relevant data
    /// PRECONDITION: camera_ray.origin is in valid domain
    /// POSTCONDITION: result.outcome describes termination condition
    TraceResult trace(const CameraRay& camera_ray);

    /// @brief Trace with custom metric parameters (for batch with varying spin, etc.)
    /// @param camera_ray Ray from camera
    /// @param mass Black hole mass M
    /// @param spin Black hole spin parameter a
    /// @return TraceResult with outcome and relevant data
    TraceResult trace(const CameraRay& camera_ray, double mass, double spin);

    // =========================================================================
    // Configuration
    // =========================================================================

    /// @brief Update configuration
    void setConfig(const TracerConfig& config) { m_Config = config; }

    /// @brief Get current configuration
    const TracerConfig& getConfig() const { return m_Config; }

    /// @brief Update metric pointer
    void setMetric(IMetric* metric) { m_Metric = metric; }

private:
    IMetric* m_Metric;
    TracerConfig m_Config;

    // =========================================================================
    // Internal Methods
    // =========================================================================

    /// @brief Initialize Lightray from CameraRay
    /// Converts from Boyer-Lindquist to Kerr-Schild Cartesian coordinates
    /// @param camera_ray Input ray in BL coordinates
    /// @return Initialized Lightray in Cartesian coordinates
    Lightray initializeLightray(const CameraRay& camera_ray);

    /// @brief Check for disk intersection between two positions
    /// @param pos_old Previous position (Cartesian)
    /// @param pos_new Current position (Cartesian)
    /// @param intersection_r [out] Radius at intersection if found
    /// @param intersection_phi [out] Azimuthal angle at intersection if found
    /// @return true if ray crossed equatorial plane within disk bounds
    bool checkDiskIntersection(const Vec4& pos_old, const Vec4& pos_new,
                               float& intersection_r, float& intersection_phi);

    /// @brief Compute disk temperature at given radius
    /// Uses Novikov-Thorne profile: T ~ r^(-3/4) * correction
    /// @param r Radius in M units
    /// @return Temperature in arbitrary units
    float computeDiskTemperature(float r);

    /// @brief Compute horizon radius from current metric parameters
    /// @return Event horizon radius in M units
    float computeHorizonRadius();

    /// @brief Check if ray state contains NaN or Inf
    /// @param ray Current ray state
    /// @return true if any component is invalid
    bool hasInvalidState(const Lightray& ray);

    /// @brief Compute g-factor (gravitational redshift + Doppler) for disk emission
    /// g = grav_factor / (gamma × (1 - v·n))
    /// @param r Radius at disk (in M units)
    /// @param phi Azimuthal angle at disk
    /// @param ray_vel Ray 4-velocity at intersection
    /// @return g-factor (intensity scales as g^4)
    float computeGFactor(float r, float phi, const Vec4& ray_vel);

    /// @brief Compute g-factor with decomposed components for motion blur
    /// Populates TraceResult with gfactor_grav, gfactor_gamma, gfactor_v_orb, gfactor_A, gfactor_B
    /// enabling analytical computation of g(φ + δφ) = grav / (gamma × (1 - v_orb × (A·cos(δφ) + B·sin(δφ))))
    /// @param r Radius at disk (in M units)
    /// @param phi Azimuthal angle at disk
    /// @param ray_vel Ray 4-velocity at intersection
    /// @param result TraceResult to populate with g-factor components
    void computeGFactorWithComponents(float r, float phi, const Vec4& ray_vel, TraceResult& result);

    /// @brief Compute Keplerian orbital velocity at radius r
    /// For Kerr: Omega = (M^(1/2)) / (r^(3/2) + a*M^(1/2))
    /// @param r Radius in M units
    /// @return Angular velocity Omega
    float computeOrbitalVelocity(float r);

    // =========================================================================
    // Volumetric Disk Methods (Phase 6)
    // =========================================================================

    /// @brief Check if position is inside volumetric disk volume
    /// @param r Cylindrical radius
    /// @param z Height above midplane
    /// @return true if (r, z) is inside the disk
    bool isInVolumetricDisk(float r, float z);

    /// @brief Compute scale height at given radius
    /// H(r) = H_over_r × r × (r/r_ref)^H_power
    /// @param r Radius in M units
    /// @return Scale height H in M units
    float computeScaleHeight(float r);

    /// @brief Compute volumetric disk density × opacity at position
    /// Returns κρ for optical depth integration
    /// @param r Cylindrical radius
    /// @param z Height above midplane
    /// @return κρ in geometric units
    float computeVolumetricOpacityDensity(float r, float z);

    /// @brief Compute volumetric disk temperature at position
    /// Includes vertical temperature gradient
    /// @param r Cylindrical radius
    /// @param z Height above midplane
    /// @return Temperature in arbitrary units
    float computeVolumetricTemperature(float r, float z);

    /// @brief Accumulate emission along ray path through volumetric disk
    /// Performs ray marching with radiative transfer
    /// @param ray Current ray state at disk entry
    /// @param entry_pos Position where ray enters disk
    /// @param exit_pos Position where ray exits disk (or current pos if still inside)
    /// @param result TraceResult to populate with volumetric emission data
    void accumulateVolumetricEmission(const Lightray& ray,
                                       const Vec4& entry_pos,
                                       const Vec4& exit_pos,
                                       TraceResult& result);
};

// =============================================================================
// Inline Implementations
// =============================================================================

inline GeodesicTracer::GeodesicTracer(IMetric* metric, const TracerConfig& config)
    : m_Metric(metric), m_Config(config) {}

inline float GeodesicTracer::computeHorizonRadius() {
    // For Kerr-Schild family, extract M and a
    // Horizon radius: r_+ = M + sqrt(M² - a²)
    // We use horizon_factor * 2M as conservative estimate for Schwarzschild
    auto params = m_Metric->getParameters();
    double M = params.count("mass") ? params.at("mass").value : 1.0;
    double a_over_M = params.count("spin") ? params.at("spin").value : 0.0;
    double a = a_over_M * M;

    double a2 = a * a;
    double M2 = M * M;

    if (a2 >= M2) {
        // Naked singularity or extremal - use 2M as fallback
        return static_cast<float>(2.0 * M * m_Config.horizon_factor);
    }

    double r_plus = M + std::sqrt(M2 - a2);
    return static_cast<float>(r_plus * m_Config.horizon_factor);
}

inline float GeodesicTracer::computeDiskTemperature(float r) {
    // Novikov-Thorne temperature profile (simplified)
    // T(r) ~ r^(-3/4) for thin disk
    // Normalized so T(disk_inner) = disk_temperature_inner

    float r_in = m_Config.disk_inner;
    float T_in = m_Config.disk_temperature_inner;

    // Temperature falls as r^(-3/4)
    float T = T_in * std::pow(r_in / r, 0.75f);

    return T;
}

inline bool GeodesicTracer::hasInvalidState(const Lightray& ray) {
    for (int i = 0; i < 4; ++i) {
        if (std::isnan(ray.position(i)) || std::isinf(ray.position(i))) return true;
        if (std::isnan(ray.velocity(i)) || std::isinf(ray.velocity(i))) return true;
    }
    return false;
}

inline bool GeodesicTracer::checkDiskIntersection(const Vec4& pos_old, const Vec4& pos_new,
                                                   float& intersection_r, float& intersection_phi) {
    // Check if ray crossed equatorial plane (z = 0 in Cartesian)
    double z_old = pos_old(3);
    double z_new = pos_new(3);

    // No crossing if same sign
    if (z_old * z_new >= 0) return false;

    // Linear interpolation to find crossing point
    double alpha = -z_old / (z_new - z_old);

    double x_cross = pos_old(1) + alpha * (pos_new(1) - pos_old(1));
    double y_cross = pos_old(2) + alpha * (pos_new(2) - pos_old(2));

    // Compute cylindrical radius at crossing
    double r_cross = std::sqrt(x_cross * x_cross + y_cross * y_cross);

    // Check if within disk bounds
    if (r_cross < m_Config.disk_inner || r_cross > m_Config.disk_outer) {
        return false;
    }

    intersection_r = static_cast<float>(r_cross);
    intersection_phi = static_cast<float>(std::atan2(y_cross, x_cross));

    return true;
}

} // namespace Sirius

#endif // GTRC001A_H
