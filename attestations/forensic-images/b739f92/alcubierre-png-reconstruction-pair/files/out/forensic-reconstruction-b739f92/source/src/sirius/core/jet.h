#pragma once

// Bipolar relativistic jet emission for AGN and X-ray binaries: synchrotron
// from a power-law electron population, with Doppler beaming along the axis.
// Reference: Blandford & Koenigl (1979), ApJ 232, 34; Rybicki & Lightman
// (1979) ch. 6. Ported from PHJT001A.h.

#include <cmath>

namespace sirius::core {

// Jet geometry, kinematics, and emission parameters.
struct JetConfig {
    // Geometry
    float r_launch = 3.0f;       // Jet launching radius [M]
    float r_max = 200.0f;        // Maximum jet extent [M]
    float opening_angle = 0.1f;  // Half-opening angle at r_launch [radians]
    float collimation = 0.5f;    // Shape parameter (0 = conical, 1 = parabolic)

    // Kinematics
    float lorentz_factor = 5.0f;    // Bulk Lorentz factor Gamma
    float velocity_profile = 0.0f;  // 0 = constant, >0 = accelerating jet

    // Emission
    float spectral_index = 2.2f;  // Electron power-law index p
    float B_field_0 = 1e4f;       // Magnetic field at r_launch [Gauss]
    float B_field_decay = 1.0f;   // B propto r^(-B_decay)
    float n_e_0 = 1e5f;           // Electron density at r_launch [cm^-3]
    float n_e_decay = 2.0f;       // n_e propto r^(-n_decay)
    float gamma_min = 10.0f;      // Minimum electron Lorentz factor
    float gamma_max = 1e6f;       // Maximum electron Lorentz factor

    // Polarisation
    bool enable_polarisation = true;
    float B_field_order = 0.5f;  // Magnetic field ordering (0 = random, 1 = ordered)
};

// Relativistic jet emission model.
class RelativisticJet {
  public:
    explicit RelativisticJet(const JetConfig& config = JetConfig()) : config_(config) {
        // Precompute velocity from Lorentz factor
        beta_ = std::sqrt(1.0f - 1.0f / (config_.lorentz_factor * config_.lorentz_factor));
    }

    // Whether the Cartesian point (x, y, z) [M] lies inside the jet volume;
    // jet_sign = +1 for the northern jet, -1 for the southern.
    bool IsInsideJet(float x, float y, float z, int jet_sign = 1) const {
        // Height along jet axis
        float h = z * jet_sign;
        if (h < config_.r_launch || h > config_.r_max) return false;

        // Cylindrical radius
        float rho = std::sqrt(x * x + y * y);

        // Jet radius at this height (parabolic/conical blend)
        float r_jet = JetRadius(h);

        return rho <= r_jet;
    }

    // Jet radius at height h: r = r_0 * (h/h_0)^collimation * tan(opening_angle).
    float JetRadius(float h) const {
        if (h <= config_.r_launch) return 0.0f;

        float h_ratio = h / config_.r_launch;
        float shape_factor = std::pow(h_ratio, config_.collimation);
        return config_.r_launch * std::tan(config_.opening_angle) * shape_factor;
    }

    // Normalised jet velocity vector; jet_sign selects north (+1) or south (-1).
    void GetVelocity(float& vx, float& vy, float& vz, int jet_sign = 1) const {
        vx = 0.0f;
        vy = 0.0f;
        vz = beta_ * static_cast<float>(jet_sign);
    }

    // Doppler factor D = 1 / [Gamma(1 - beta cos theta)] for line of sight at
    // cos_theta to the jet axis (approaching side D > 1, receding D < 1).
    float DopplerFactor(float cos_theta) const {
        float Gamma = config_.lorentz_factor;
        return 1.0f / (Gamma * (1.0f - beta_ * cos_theta));
    }

    // Doppler-boosted intensity I_obs = D^n I_emit; n = 2 + alpha (optically
    // thin) or 2 - alpha (optically thick), with alpha the synchrotron index.
    float BoostedIntensity(float I_emit, float cos_theta, bool optically_thin = true) const {
        float D = DopplerFactor(cos_theta);
        float alpha = (config_.spectral_index - 1.0f) / 2.0f;  // Synchrotron spectral index

        float n = optically_thin ? (2.0f + alpha) : (2.0f - alpha);
        return I_emit * std::pow(D, n);
    }

    // Magnetic field strength B(r) = B_0 * (r_launch / r)^B_decay.
    float MagneticField(float r) const {
        if (r <= config_.r_launch) return config_.B_field_0;
        float ratio = config_.r_launch / r;
        return config_.B_field_0 * std::pow(ratio, config_.B_field_decay);
    }

    // Electron density n_e(r) = n_e_0 * (r_launch / r)^n_decay.
    float ElectronDensity(float r) const {
        if (r <= config_.r_launch) return config_.n_e_0;
        float ratio = config_.r_launch / r;
        return config_.n_e_0 * std::pow(ratio, config_.n_e_decay);
    }

    // Synchrotron emissivity (comoving frame): j_nu propto n_e * B^((p+1)/2).
    float SynchrotronEmissivity(float r) const {
        float B = MagneticField(r);
        float n_e = ElectronDensity(r);
        float p = config_.spectral_index;

        // j_nu propto n_e * B^((p+1)/2)
        float B_power = std::pow(B / config_.B_field_0, (p + 1.0f) / 2.0f);
        return n_e / config_.n_e_0 * B_power;
    }

    // Jet intensity at position (x, y, z) [M] seen by an observer at
    // (observer_x, observer_y, observer_z) [M], summed over both jets.
    float ComputeEmission(float x, float y, float z, float observer_x, float observer_y,
                          float observer_z) const {
        float total_emission = 0.0f;

        // Check both jets
        for (int jet_sign = -1; jet_sign <= 1; jet_sign += 2) {
            if (!IsInsideJet(x, y, z, jet_sign)) continue;

            [[maybe_unused]] float h = std::abs(z);
            float r = std::sqrt(x * x + y * y + z * z);

            // Line of sight direction (from emission point to observer)
            float dx = observer_x - x;
            float dy = observer_y - y;
            float dz = observer_z - z;
            float dist = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (dist < 1e-6f) continue;

            // Cosine of angle between jet axis and line of sight
            float cos_theta = (dz / dist) * static_cast<float>(jet_sign);

            // Comoving frame emissivity
            float j_comoving = SynchrotronEmissivity(r);

            // Apply Doppler boosting
            float I_obs = BoostedIntensity(j_comoving, cos_theta, true);

            total_emission += I_obs;
        }

        return total_emission;
    }

    // Synchrotron polarisation degree: pi = (p+1)/(p+7/3) * B_order.
    float PolarisationDegree() const {
        float p = config_.spectral_index;
        float pi_max = (p + 1.0f) / (p + 7.0f / 3.0f);
        return pi_max * config_.B_field_order;
    }

    const JetConfig& GetConfig() const { return config_; }
    void SetConfig(const JetConfig& config) {
        config_ = config;
        beta_ = std::sqrt(1.0f - 1.0f / (config_.lorentz_factor * config_.lorentz_factor));
    }

  private:
    JetConfig config_;
    float beta_;  // v/c
};

// Ray marching over jet emission.
namespace jet_ray_marching {

// Integrate jet emission along the segment from (start_*) to (end_*), sampling
// num_samples points and accumulating emission * ds toward the observer.
inline float IntegrateJetEmission(const RelativisticJet& jet, float start_x, float start_y,
                                  float start_z, float end_x, float end_y, float end_z,
                                  float observer_x, float observer_y, float observer_z,
                                  int num_samples = 32) {
    float total = 0.0f;

    float dx = (end_x - start_x) / num_samples;
    float dy = (end_y - start_y) / num_samples;
    float dz = (end_z - start_z) / num_samples;
    float ds = std::sqrt(dx * dx + dy * dy + dz * dz);

    for (int i = 0; i < num_samples; i++) {
        float t = (static_cast<float>(i) + 0.5f) / num_samples;
        float x = start_x + t * (end_x - start_x);
        float y = start_y + t * (end_y - start_y);
        float z = start_z + t * (end_z - start_z);

        float emission = jet.ComputeEmission(x, y, z, observer_x, observer_y, observer_z);
        total += emission * ds;
    }

    return total;
}

}  // namespace jet_ray_marching

}  // namespace sirius::core
