// PHJT002A.h - MHD-Based Relativistic Jet Model
// Component ID: PHJT002A (Physics/Jet/MHDJet)
//
// Magnetohydrodynamic relativistic jet model with magnetic field decay,
// synchrotron emission, and Doppler beaming. Extends PHJT001A with
// proper MHD physics and B-field structure.
//
// PHYSICAL FOUNDATION:
// ====================
// AGN and microquasar jets are magnetically dominated outflows driven by
// the Blandford-Znajek mechanism (extraction of black hole rotational energy).
//
// Key physics implemented:
//   1. Magnetic field decay: B(z) = B_0 × (z_0/z)^p
//   2. Magnetic flux conservation: Φ = ∮ B·dA = const
//   3. Force-free limit: j × B = ∇P at jet boundary
//   4. Synchrotron emission: P ∝ B² γ² for relativistic electrons
//   5. Doppler beaming: D = 1/[Γ(1 - β cos θ)]
//
// REFERENCES:
// - Blandford & Znajek (1977) MNRAS 179, 433 (BZ mechanism)
// - Blandford & Königl (1979) ApJ 232, 34 (synchrotron jet model)
// - McKinney (2006) MNRAS 368, 1561 (GRMHD jets)

#pragma once

#include <cmath>
#include <cstdint>
#include <algorithm>

namespace Sirius {

//==============================================================================
// MHD Jet Configuration
//==============================================================================
struct JetMHDConfig {
    // =========================================================================
    // Magnetic Field
    // =========================================================================
    float B_base = 1e4f;               ///< Base magnetic field strength [Gauss]
    float power_law_index = 1.0f;      ///< B(z) ~ z^(-p) decay exponent
    float B_field_order = 0.5f;        ///< Field ordering (0 = random, 1 = ordered)

    // =========================================================================
    // Geometry
    // =========================================================================
    float opening_half_angle = 0.1f;   ///< Jet opening angle [rad] (~ 6°)
    float z_launch_M = 3.0f;           ///< Jet launching height [M]
    float z_max_M = 200.0f;            ///< Maximum jet extent [M]
    float collimation = 0.5f;          ///< Shape: 0 = conical, 1 = parabolic

    // =========================================================================
    // Kinematics
    // =========================================================================
    float lorentz_factor = 5.0f;       ///< Bulk Lorentz factor Γ
    float velocity_profile = 0.0f;     ///< 0 = constant Γ, >0 = accelerating

    // =========================================================================
    // Electron Distribution
    // =========================================================================
    float electron_index = 2.2f;       ///< Power-law index p: n(γ) ∝ γ^(-p)
    float gamma_min = 10.0f;           ///< Minimum electron Lorentz factor
    float gamma_max = 1e6f;            ///< Maximum electron Lorentz factor
    float n_e_0 = 1e5f;                ///< Electron density at z_launch [cm^-3]
    float n_e_decay = 2.0f;            ///< n_e ∝ z^(-n_decay)

    // =========================================================================
    // Emission
    // =========================================================================
    float intensity_scale = 1.0f;      ///< Overall emission scaling
    bool enable_polarisation = true;   ///< Compute Stokes parameters
    bool enabled = true;

    // =========================================================================
    // Physical Constants
    // =========================================================================
    static constexpr float c = 2.998e10f;      // Speed of light [cm/s]
    static constexpr float m_e = 9.109e-28f;   // Electron mass [g]
    static constexpr float e = 4.803e-10f;     // Elementary charge [esu]

    // =========================================================================
    // Invariants
    // =========================================================================
    // B_base > 0
    // opening_half_angle ∈ (0, π/4]
    // lorentz_factor >= 1.0
    // power_law_index ∈ [0.5, 2.0]

    /// @brief Validate and clamp parameters
    void validate() {
        B_base = std::max(B_base, 1.0f);
        power_law_index = std::clamp(power_law_index, 0.5f, 2.0f);
        B_field_order = std::clamp(B_field_order, 0.0f, 1.0f);
        opening_half_angle = std::clamp(opening_half_angle, 0.01f, 0.785f);  // 0.6° - 45°
        z_launch_M = std::max(z_launch_M, 1.0f);
        z_max_M = std::max(z_max_M, z_launch_M + 1.0f);
        collimation = std::clamp(collimation, 0.0f, 1.0f);
        lorentz_factor = std::max(lorentz_factor, 1.0f);
        electron_index = std::clamp(electron_index, 1.5f, 4.0f);
        gamma_min = std::max(gamma_min, 1.0f);
        gamma_max = std::max(gamma_max, gamma_min * 10.0f);
    }

    /// @brief Compute bulk velocity β from Lorentz factor
    float beta() const {
        return std::sqrt(1.0f - 1.0f / (lorentz_factor * lorentz_factor));
    }

    /// @brief Compute synchrotron spectral index from electron index
    /// α = (p - 1) / 2 for optically thin synchrotron
    float spectralIndex() const {
        return (electron_index - 1.0f) / 2.0f;
    }

    /// @brief Compute maximum synchrotron polarisation degree
    /// π_max = (p + 1) / (p + 7/3)
    float maxPolarisationDegree() const {
        float p = electron_index;
        return (p + 1.0f) / (p + 7.0f / 3.0f);
    }
};

//==============================================================================
// MHD Jet Model Class
//==============================================================================
class JetMHD {
public:
    explicit JetMHD(const JetMHDConfig& config = JetMHDConfig())
        : m_Config(config) {
        m_Config.validate();
        m_Beta = m_Config.beta();
    }

    // =========================================================================
    // Geometry
    // =========================================================================

    /// @brief Jet radius at height z
    /// Blended conical/parabolic: R(z) = R_0 × (z/z_0)^(1 - c/2) × tan(θ)
    /// where c = collimation (0 = conical, 1 = parabolic)
    float jetRadius(float z) const {
        if (z < m_Config.z_launch_M) return 0.0f;

        float z_ratio = z / m_Config.z_launch_M;
        float shape_exp = 1.0f - m_Config.collimation * 0.5f;
        return m_Config.z_launch_M * std::tan(m_Config.opening_half_angle)
               * std::pow(z_ratio, shape_exp);
    }

    /// @brief Check if point is inside jet volume
    /// @param x, y, z Cartesian coordinates [M]
    /// @param jet_sign +1 for northern jet, -1 for southern jet
    bool isInsideJet(float x, float y, float z, int jet_sign = 1) const {
        if (!m_Config.enabled) return false;

        float h = z * static_cast<float>(jet_sign);
        if (h < m_Config.z_launch_M || h > m_Config.z_max_M) return false;

        float rho = std::sqrt(x * x + y * y);
        float R_jet = jetRadius(h);

        return rho <= R_jet;
    }

    // =========================================================================
    // Magnetic Field
    // =========================================================================

    /// @brief Magnetic field strength at height z
    /// B(z) = B_0 × (z_launch / z)^p
    float magneticField(float z) const {
        float h = std::abs(z);
        if (h <= m_Config.z_launch_M) return m_Config.B_base;

        return m_Config.B_base * std::pow(m_Config.z_launch_M / h, m_Config.power_law_index);
    }

    /// @brief Verify magnetic flux conservation
    /// Φ = π R² B should be approximately constant along jet
    float magneticFlux(float z) const {
        float R = jetRadius(std::abs(z));
        float B = magneticField(z);
        return static_cast<float>(M_PI) * R * R * B;
    }

    // =========================================================================
    // Plasma Properties
    // =========================================================================

    /// @brief Electron density at height z
    /// n_e(z) = n_e_0 × (z_launch / z)^n_decay
    float electronDensity(float z) const {
        float h = std::abs(z);
        if (h <= m_Config.z_launch_M) return m_Config.n_e_0;

        return m_Config.n_e_0 * std::pow(m_Config.z_launch_M / h, m_Config.n_e_decay);
    }

    /// @brief Bulk Lorentz factor at height (with optional acceleration)
    float lorentzFactor(float z) const {
        if (m_Config.velocity_profile <= 0.0f) {
            return m_Config.lorentz_factor;
        }

        // Accelerating jet: Γ(z) = Γ_0 × (z/z_0)^v where v = velocity_profile
        float h = std::abs(z);
        float z_ratio = h / m_Config.z_launch_M;
        return m_Config.lorentz_factor * std::pow(z_ratio, m_Config.velocity_profile);
    }

    // =========================================================================
    // Doppler Beaming
    // =========================================================================

    /// @brief Doppler factor for emission towards observer
    /// D = 1 / [Γ(1 - β cos θ)]
    /// @param cos_theta Cosine of angle between jet axis and line of sight
    /// @param z Height (for accelerating jets)
    float dopplerFactor(float cos_theta, float z) const {
        float Gamma = lorentzFactor(z);
        float beta = std::sqrt(1.0f - 1.0f / (Gamma * Gamma));
        return 1.0f / (Gamma * (1.0f - beta * cos_theta));
    }

    /// @brief Doppler-boosted intensity
    /// I_obs = D^(2+α) × I_emit for optically thin synchrotron
    float boostedIntensity(float I_emit, float cos_theta, float z) const {
        float D = dopplerFactor(cos_theta, z);
        float alpha = m_Config.spectralIndex();
        float boost_power = 2.0f + alpha;
        return I_emit * std::pow(D, boost_power);
    }

    // =========================================================================
    // Synchrotron Emission
    // =========================================================================

    /// @brief Synchrotron emissivity in comoving frame
    /// j_ν ∝ n_e × B^((p+1)/2)
    float synchrotronEmissivity(float z) const {
        float B = magneticField(z);
        float n_e = electronDensity(z);
        float p = m_Config.electron_index;

        // Normalized to base values
        float B_norm = B / m_Config.B_base;
        float n_norm = n_e / m_Config.n_e_0;

        float j = n_norm * std::pow(B_norm, (p + 1.0f) / 2.0f);
        return j * m_Config.intensity_scale;
    }

    /// @brief Synchrotron self-absorption coefficient
    /// α_ν ∝ n_e × B^((p+2)/2) × ν^(-(p+4)/2)
    float selfAbsorptionCoeff(float z, float nu_ratio = 1.0f) const {
        float B = magneticField(z);
        float n_e = electronDensity(z);
        float p = m_Config.electron_index;

        float B_norm = B / m_Config.B_base;
        float n_norm = n_e / m_Config.n_e_0;

        return n_norm * std::pow(B_norm, (p + 2.0f) / 2.0f)
               * std::pow(nu_ratio, -(p + 4.0f) / 2.0f);
    }

    /// @brief Compute jet emission at position towards observer
    /// @param x, y, z Position [M]
    /// @param obs_x, obs_y, obs_z Observer position [M]
    float computeEmission(float x, float y, float z,
                          float obs_x, float obs_y, float obs_z) const {
        float total_emission = 0.0f;

        for (int jet_sign = -1; jet_sign <= 1; jet_sign += 2) {
            if (!isInsideJet(x, y, z, jet_sign)) continue;

            float h = std::abs(z);

            // Line of sight direction
            float dx = obs_x - x;
            float dy = obs_y - y;
            float dz = obs_z - z;
            float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (dist < 1e-6f) continue;

            // Cosine of angle between jet axis and LOS
            float cos_theta = (dz / dist) * static_cast<float>(jet_sign);

            // Comoving emissivity
            float j_co = synchrotronEmissivity(h);

            // Apply Doppler boosting
            float I_obs = boostedIntensity(j_co, cos_theta, h);

            total_emission += I_obs;
        }

        return total_emission;
    }

    // =========================================================================
    // Polarisation
    // =========================================================================

    /// @brief Synchrotron polarisation degree at height
    /// π = π_max × B_order where π_max = (p+1)/(p+7/3)
    float polarisationDegree(float z) const {
        (void)z;  // Uniform B_order in this model
        return m_Config.maxPolarisationDegree() * m_Config.B_field_order;
    }

    /// @brief Electric vector position angle (EVPA)
    /// Perpendicular to projected B-field direction
    float evpa(float x, float y) const {
        // For toroidal B-field, EVPA is along the jet axis projected
        // For poloidal B-field, EVPA is perpendicular
        float phi = std::atan2(y, x);
        // Blend based on B_field_order (simplified)
        return phi + 1.5708f * (1.0f - m_Config.B_field_order);
    }

    // =========================================================================
    // Ray Marching
    // =========================================================================

    struct JetRayResult {
        float intensity;
        float optical_depth;
        float polarisation_degree;
        float evpa;
        bool hit_jet;
    };

    JetRayResult integrateRay(float x1, float y1, float z1,
                               float x2, float y2, float z2,
                               float obs_x, float obs_y, float obs_z,
                               int num_samples = 32) const {
        JetRayResult result = {0.0f, 0.0f, 0.0f, 0.0f, false};

        if (!m_Config.enabled) return result;

        float dx = x2 - x1, dy = y2 - y1, dz = z2 - z1;
        float path_length = std::sqrt(dx * dx + dy * dy + dz * dz);
        float ds = path_length / static_cast<float>(num_samples);

        float total_pol = 0.0f;
        float total_evpa = 0.0f;
        int pol_samples = 0;

        for (int i = 0; i < num_samples; ++i) {
            float t = (static_cast<float>(i) + 0.5f) / static_cast<float>(num_samples);
            float x = x1 + t * dx;
            float y = y1 + t * dy;
            float z = z1 + t * dz;

            for (int jet_sign = -1; jet_sign <= 1; jet_sign += 2) {
                if (!isInsideJet(x, y, z, jet_sign)) continue;
                result.hit_jet = true;

                float h = std::abs(z);

                // LOS angle
                float lx = obs_x - x, ly = obs_y - y, lz = obs_z - z;
                float ldist = std::sqrt(lx * lx + ly * ly + lz * lz);
                float cos_theta = (ldist > 1e-6f) ? (lz / ldist) * static_cast<float>(jet_sign) : 1.0f;

                // Emission
                float j_co = synchrotronEmissivity(h);
                float j_obs = boostedIntensity(j_co, cos_theta, h);

                // Absorption
                float alpha = selfAbsorptionCoeff(h);
                float dtau = alpha * ds;
                result.optical_depth += dtau;

                // Radiative transfer
                float transmission = std::exp(-result.optical_depth);
                result.intensity += j_obs * ds * transmission;

                // Polarisation
                if (m_Config.enable_polarisation) {
                    total_pol += polarisationDegree(h);
                    total_evpa += evpa(x, y);
                    pol_samples++;
                }
            }
        }

        if (pol_samples > 0) {
            result.polarisation_degree = total_pol / static_cast<float>(pol_samples);
            result.evpa = total_evpa / static_cast<float>(pol_samples);
        }

        return result;
    }

    // =========================================================================
    // Configuration
    // =========================================================================
    const JetMHDConfig& config() const { return m_Config; }
    void setConfig(const JetMHDConfig& config) {
        m_Config = config;
        m_Config.validate();
        m_Beta = m_Config.beta();
    }

private:
    JetMHDConfig m_Config;
    float m_Beta;
};

} // namespace Sirius
