// PHJT001A.h - Relativistic Jet Emission Model
// Component ID: PHJT001A (Physics/Jet/JetModel)
//
// Implements bipolar relativistic jet emission for AGN and X-ray binaries.
// Based on synchrotron radiation from shock-accelerated electrons.
//
// PHYSICAL FOUNDATION:
// ====================
// Jets are collimated outflows of plasma moving at relativistic speeds
// along the rotation axis of the black hole. Key physics:
//
// 1. Geometry: Parabolic or conical expansion from launching radius
// 2. Velocity: Bulk Lorentz factor Γ ~ 2-10 (microquasars) or 10-50 (AGN)
// 3. Emission: Synchrotron from power-law electron distribution
// 4. Beaming: Doppler boosting D = 1 / [Γ(1 - β cos θ)]
//
// SYNCHROTRON SPECTRUM:
// =====================
// For power-law electrons n(γ) ∝ γ^(-p), the synchrotron spectrum is:
//   j_ν ∝ B^((p+1)/2) × ν^(-(p-1)/2)
//
// Typical values: p = 2.2 gives spectral index α = 0.6
//
// Reference: Blandford & Königl (1979), ApJ 232, 34
//            Rybicki & Lightman (1979) Chapter 6

#pragma once

#include <cmath>
#include <algorithm>

namespace Sirius {

//==============================================================================
// Jet Geometry Configuration
//==============================================================================
struct JetConfig {
    // Geometry
    float r_launch = 3.0f;      ///< Jet launching radius [M]
    float r_max = 200.0f;       ///< Maximum jet extent [M]
    float opening_angle = 0.1f;  ///< Half-opening angle at r_launch [radians]
    float collimation = 0.5f;    ///< Shape parameter (0 = conical, 1 = parabolic)

    // Kinematics
    float lorentz_factor = 5.0f; ///< Bulk Lorentz factor Γ
    float velocity_profile = 0.0f; ///< 0 = constant, >0 = accelerating jet

    // Emission
    float spectral_index = 2.2f;  ///< Electron power-law index p
    float B_field_0 = 1e4f;       ///< Magnetic field at r_launch [Gauss]
    float B_field_decay = 1.0f;   ///< B ∝ r^(-B_decay)
    float n_e_0 = 1e5f;           ///< Electron density at r_launch [cm^-3]
    float n_e_decay = 2.0f;       ///< n_e ∝ r^(-n_decay)
    float gamma_min = 10.0f;      ///< Minimum electron Lorentz factor
    float gamma_max = 1e6f;       ///< Maximum electron Lorentz factor

    // Polarisation
    bool enable_polarisation = true;
    float B_field_order = 0.5f;   ///< Magnetic field ordering (0 = random, 1 = ordered)
};

//==============================================================================
// Relativistic Jet Model
//==============================================================================
class RelativisticJet {
public:
    explicit RelativisticJet(const JetConfig& config = JetConfig()) : m_Config(config) {
        // Precompute velocity from Lorentz factor
        m_Beta = std::sqrt(1.0f - 1.0f / (m_Config.lorentz_factor * m_Config.lorentz_factor));
    }

    //--------------------------------------------------------------------------
    // Geometry Methods
    //--------------------------------------------------------------------------

    /// @brief Check if point is inside jet volume
    /// @param x, y, z Cartesian coordinates [M]
    /// @param jet_sign +1 for northern jet, -1 for southern jet
    /// @return true if point is inside jet
    bool isInsideJet(float x, float y, float z, int jet_sign = 1) const {
        // Height along jet axis
        float h = z * jet_sign;
        if (h < m_Config.r_launch || h > m_Config.r_max) return false;

        // Cylindrical radius
        float rho = std::sqrt(x*x + y*y);

        // Jet radius at this height (parabolic/conical blend)
        float r_jet = jetRadius(h);

        return rho <= r_jet;
    }

    /// @brief Compute jet radius at height h
    /// Shape: r = r_0 × (h/h_0)^collimation × tan(opening_angle)
    float jetRadius(float h) const {
        if (h <= m_Config.r_launch) return 0.0f;

        float h_ratio = h / m_Config.r_launch;
        float shape_factor = std::pow(h_ratio, m_Config.collimation);
        return m_Config.r_launch * std::tan(m_Config.opening_angle) * shape_factor;
    }

    //--------------------------------------------------------------------------
    // Velocity and Beaming
    //--------------------------------------------------------------------------

    /// @brief Get jet velocity vector (normalised)
    /// @param jet_sign +1 for northern jet, -1 for southern jet
    void getVelocity(float& vx, float& vy, float& vz, int jet_sign = 1) const {
        vx = 0.0f;
        vy = 0.0f;
        vz = m_Beta * static_cast<float>(jet_sign);
    }

    /// @brief Compute Doppler factor for observer at angle θ to jet axis
    /// D = 1 / [Γ(1 - β cos θ)]
    /// @param cos_theta Cosine of angle between jet axis and line of sight
    /// @return Doppler factor D (approaching side: D > 1, receding: D < 1)
    float dopplerFactor(float cos_theta) const {
        float Gamma = m_Config.lorentz_factor;
        return 1.0f / (Gamma * (1.0f - m_Beta * cos_theta));
    }

    /// @brief Doppler boosted intensity
    /// I_obs = D^n × I_emit, where n depends on emission mechanism
    /// For optically thin synchrotron: n = 2 + α (α = spectral index)
    /// For optically thick: n = 2 - α
    /// @param I_emit Emitted intensity (comoving frame)
    /// @param cos_theta Angle to line of sight
    /// @param optically_thin true for thin synchrotron
    float boostedIntensity(float I_emit, float cos_theta, bool optically_thin = true) const {
        float D = dopplerFactor(cos_theta);
        float alpha = (m_Config.spectral_index - 1.0f) / 2.0f;  // Synchrotron spectral index

        float n = optically_thin ? (2.0f + alpha) : (2.0f - alpha);
        return I_emit * std::pow(D, n);
    }

    //--------------------------------------------------------------------------
    // Emission
    //--------------------------------------------------------------------------

    /// @brief Magnetic field strength at position
    /// B(r) = B_0 × (r_launch / r)^B_decay
    float magneticField(float r) const {
        if (r <= m_Config.r_launch) return m_Config.B_field_0;
        float ratio = m_Config.r_launch / r;
        return m_Config.B_field_0 * std::pow(ratio, m_Config.B_field_decay);
    }

    /// @brief Electron density at position
    /// n_e(r) = n_e_0 × (r_launch / r)^n_decay
    float electronDensity(float r) const {
        if (r <= m_Config.r_launch) return m_Config.n_e_0;
        float ratio = m_Config.r_launch / r;
        return m_Config.n_e_0 * std::pow(ratio, m_Config.n_e_decay);
    }

    /// @brief Synchrotron emissivity (comoving frame)
    /// j_ν ∝ n_e × B^((p+1)/2)
    /// @param r Distance from black hole [M]
    /// @return Emissivity in arbitrary units
    float synchrotronEmissivity(float r) const {
        float B = magneticField(r);
        float n_e = electronDensity(r);
        float p = m_Config.spectral_index;

        // j_ν ∝ n_e × B^((p+1)/2)
        float B_power = std::pow(B / m_Config.B_field_0, (p + 1.0f) / 2.0f);
        return n_e / m_Config.n_e_0 * B_power;
    }

    /// @brief Compute jet emission at position
    /// @param x, y, z Position [M]
    /// @param observer_x, observer_y, observer_z Observer position [M]
    /// @return Intensity (including Doppler boosting)
    float computeEmission(float x, float y, float z,
                          float observer_x, float observer_y, float observer_z) const {
        float total_emission = 0.0f;

        // Check both jets
        for (int jet_sign = -1; jet_sign <= 1; jet_sign += 2) {
            if (!isInsideJet(x, y, z, jet_sign)) continue;

            float h = std::abs(z);
            float r = std::sqrt(x*x + y*y + z*z);

            // Line of sight direction (from emission point to observer)
            float dx = observer_x - x;
            float dy = observer_y - y;
            float dz = observer_z - z;
            float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (dist < 1e-6f) continue;

            // Cosine of angle between jet axis and line of sight
            float cos_theta = (dz / dist) * static_cast<float>(jet_sign);

            // Comoving frame emissivity
            float j_comoving = synchrotronEmissivity(r);

            // Apply Doppler boosting
            float I_obs = boostedIntensity(j_comoving, cos_theta, true);

            total_emission += I_obs;
        }

        return total_emission;
    }

    /// @brief Compute polarisation degree of jet emission
    /// Synchrotron polarisation: π = (p+1)/(p+7/3) × B_order
    float polarisationDegree() const {
        float p = m_Config.spectral_index;
        float pi_max = (p + 1.0f) / (p + 7.0f/3.0f);
        return pi_max * m_Config.B_field_order;
    }

    //--------------------------------------------------------------------------
    // Configuration
    //--------------------------------------------------------------------------
    const JetConfig& getConfig() const { return m_Config; }
    void setConfig(const JetConfig& config) {
        m_Config = config;
        m_Beta = std::sqrt(1.0f - 1.0f / (m_Config.lorentz_factor * m_Config.lorentz_factor));
    }

private:
    JetConfig m_Config;
    float m_Beta;  ///< v/c
};

//==============================================================================
// Jet Ray Marching
//==============================================================================
namespace JetRayMarching {

/// @brief Sample jet emission along ray segment
/// @param jet Jet model
/// @param start_x, start_y, start_z Ray start position
/// @param end_x, end_y, end_z Ray end position
/// @param observer_x, observer_y, observer_z Observer position
/// @param num_samples Number of samples along ray
/// @return Integrated emission along ray
inline float integrateJetEmission(const RelativisticJet& jet,
                                   float start_x, float start_y, float start_z,
                                   float end_x, float end_y, float end_z,
                                   float observer_x, float observer_y, float observer_z,
                                   int num_samples = 32) {
    float total = 0.0f;

    float dx = (end_x - start_x) / num_samples;
    float dy = (end_y - start_y) / num_samples;
    float dz = (end_z - start_z) / num_samples;
    float ds = std::sqrt(dx*dx + dy*dy + dz*dz);

    for (int i = 0; i < num_samples; i++) {
        float t = (static_cast<float>(i) + 0.5f) / num_samples;
        float x = start_x + t * (end_x - start_x);
        float y = start_y + t * (end_y - start_y);
        float z = start_z + t * (end_z - start_z);

        float emission = jet.computeEmission(x, y, z, observer_x, observer_y, observer_z);
        total += emission * ds;
    }

    return total;
}

} // namespace JetRayMarching

} // namespace Sirius
