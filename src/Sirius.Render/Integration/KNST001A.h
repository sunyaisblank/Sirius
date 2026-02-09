// KNST001A.h - Spectral Transport
// Component ID: KNST001A
// Integrates spectral radiance along geodesics with redshift
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_KNST001A_H
#define SIRIUS_RENDER_KNST001A_H

#include "KNAA001A.h"
#include <PHAD001A.h>
#include <MTSB001A.h>  // For SpectralRadiance
#include <MTTP001A.h>  // For Vec4d
#include <cmath>
#include <vector>

// SpectralRadiance is defined in MTSB001A.h in Sirius namespace

namespace sirius::kernel {

//==============================================================================
// SpectralTransport
// Integrates radiance along beam path with redshift
//==============================================================================

class SpectralTransport {
public:
    struct IntegrationResult {
        Sirius::SpectralRadiance totalRadiance;
        int diskHits = 0;
        bool hitHorizon = false;
    };

    explicit SpectralTransport(const Sirius::AccretionDiskD* disk)
        : m_disk(disk) {}

    //--------------------------------------------------------------------------
    // G-Factor Computation
    //--------------------------------------------------------------------------

    /// Compute gravitational redshift factor for disk emission at equator
    ///
    /// FORMULA (exact for equatorial disk):
    /// g_grav = √(1 - 2Mr/(r² + a²))
    ///
    /// Derivation: From Kerr metric, g_tt = -(1 - 2Mr/(r² + a²cos²θ))
    /// At equator (θ = π/2): g_tt = -(1 - 2Mr/(r² + a²))
    /// Gravitational redshift: g = √(-g_tt)
    ///
    /// For Schwarzschild (a=0): reduces to g = √(1 - 2M/r)
    ///
    /// @param r Radius in units of M (assuming M=1 normalization)
    /// @param a Spin parameter a/M
    /// @return Gravitational redshift factor (0 < g < 1)
    static double computeDiskGFactor(double r, double a) {
        // Full Kerr formula at equator: g = √(1 - 2Mr/(r² + a²))
        // With M=1 normalization: g = √(1 - 2r/(r² + a²))
        double r2 = r * r;
        double a2 = a * a;
        double arg = 1.0 - 2.0 * r / (r2 + a2);

        // Clamp for numerical safety (arg should be positive outside horizon)
        arg = std::max(arg, 1e-10);

        return std::sqrt(arg);
    }

    //--------------------------------------------------------------------------
    // Emission Evaluation
    //--------------------------------------------------------------------------

    /// Evaluate emission at a single point
    Sirius::SpectralRadiance evaluateEmission(const BeamStateD& beam,
                                                         double g) const {
        if (!m_disk) return Sirius::SpectralRadiance();

        double r = beam.x.r;
        double theta = beam.x.theta;

        // Check if in disk plane
        if (std::abs(theta - M_PI / 2.0) > 0.01) {
            return Sirius::SpectralRadiance();
        }

        // Get disk temperature at this radius
        double T = m_disk->effectiveTemperature(r);
        if (T <= 0) return Sirius::SpectralRadiance();

        // Create blackbody spectrum and apply redshift
        auto emission = Sirius::SpectralRadiance::blackbody(T);
        return emission.applyRedshift(g);
    }

    //--------------------------------------------------------------------------
    // Path Integration
    //--------------------------------------------------------------------------

    /// Integrate radiance along beam path
    IntegrationResult integrate(const std::vector<BeamStateD>& path,
                                const std::vector<double>& gFactors) const {
        IntegrationResult result;

        if (path.empty()) return result;

        for (size_t i = 0; i < path.size(); ++i) {
            const auto& beam = path[i];
            double g = (i < gFactors.size()) ? gFactors[i] : 1.0;

            // Check for horizon
            double M = m_disk ? m_disk->config().M : 1.0;
            if (beam.x.r < 2.0 * M * 1.001) {
                result.hitHorizon = true;
                break;
            }

            // Check for disk intersection
            if (std::abs(beam.x.theta - M_PI / 2.0) < 0.01) {
                auto emission = evaluateEmission(beam, g);
                if (emission.totalEnergy() > 0) {
                    result.totalRadiance += emission;
                    result.diskHits++;
                }
            }
        }

        return result;
    }

private:
    const Sirius::AccretionDiskD* m_disk;
};

} // namespace sirius::kernel

#endif // SIRIUS_RENDER_KNST001A_H
