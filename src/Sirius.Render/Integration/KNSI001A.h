// KNSI001A.h - Scene Intersection
// Component ID: KNSI001A
// Determines beam termination condition (disk, horizon, celestial sphere)
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_KNSI001A_H
#define SIRIUS_RENDER_KNSI001A_H

#include "KNAA001A.h"
#include "KNST001A.h"  // For SpectralRadiance
#include <PHAD001A.h>
#include <cmath>

namespace sirius::kernel {

//==============================================================================
// DiskJacobian
// 2D Jacobian for disk surface coordinates (r, phi)
//==============================================================================

struct DiskJacobian {
    double J[2][2] = {};  // ∂(r,φ)/∂(θ₀,φ₀)
    double determinant = 0.0;
};

//==============================================================================
// SceneIntersector
// Determines beam termination and extracts intersection data
//==============================================================================

class SceneIntersector {
public:
    struct Config {
        double escapeRadius;
        double horizonFactor;

        Config() : escapeRadius(100.0), horizonFactor(1.001) {}
    };

    struct IntersectionResult {
        enum Type {
            NONE,
            DISK,
            HORIZON,
            CELESTIAL_SPHERE
        };

        Type type = NONE;
        double r = 0.0;       // Radial coordinate at intersection
        double theta = 0.0;   // Polar angle
        double phi = 0.0;     // Azimuthal angle
        double magnification = 1.0;  // |det(J)|
        double g_factor = 1.0;       // Gravitational redshift factor
        Sirius::SpectralRadiance emission;  // Disk emission at intersection
    };

    explicit SceneIntersector(const Sirius::AccretionDiskD* disk)
        : m_disk(disk), m_config() {}

    SceneIntersector(const Sirius::AccretionDiskD* disk, const Config& config)
        : m_disk(disk), m_config(config) {}

    //--------------------------------------------------------------------------
    // Scene Intersection
    //--------------------------------------------------------------------------

    /// Determine intersection type and extract data
    IntersectionResult intersectScene(const BeamStateD& beam) const {
        IntersectionResult result;

        double r = beam.x.r;
        double theta = beam.x.theta;
        double phi = beam.x.phi;

        // Check if beam was terminated (fell into horizon)
        if (beam.terminated) {
            // Check if inside event horizon
            double M = m_disk ? m_disk->config().M : 1.0;
            double rH = 2.0 * M * m_config.horizonFactor;
            if (r < rH) {
                result.type = IntersectionResult::HORIZON;
                result.r = r;
                result.theta = theta;
                result.phi = phi;
                return result;
            }
        }

        // Check if escaped to celestial sphere
        if (r > m_config.escapeRadius) {
            result.type = IntersectionResult::CELESTIAL_SPHERE;
            result.r = r;
            result.theta = theta;
            result.phi = phi;
            return result;
        }

        // Check if in disk plane (theta ≈ π/2)
        if (m_disk && std::abs(theta - M_PI / 2.0) < 0.01) {
            double rISCO = Sirius::AccretionDiskD::computeISCO(
                m_disk->config().a_star);
            double rOuter = m_disk->config().r_outer * m_disk->config().M;

            if (r >= rISCO && r <= rOuter) {
                result.type = IntersectionResult::DISK;
                result.r = r;
                result.theta = theta;
                result.phi = phi;
                return result;
            }
        }

        result.type = IntersectionResult::NONE;
        return result;
    }

    //--------------------------------------------------------------------------
    // Jacobian Transformation
    //--------------------------------------------------------------------------

    /// Extract 2D disk Jacobian from 4D beam Jacobian
    DiskJacobian transformJacobianToDisk(const BeamStateD& beam) const {
        DiskJacobian dj;

        // Extract ∂(r,φ)/∂(θ₀,φ₀) from 4D Jacobian
        // beam.J[i][j] = ∂x^i/∂x₀^j
        // We want: ∂r/∂θ₀, ∂r/∂φ₀, ∂φ/∂θ₀, ∂φ/∂φ₀

        dj.J[0][0] = beam.J[1][2];  // ∂r/∂θ₀
        dj.J[0][1] = beam.J[1][3];  // ∂r/∂φ₀
        dj.J[1][0] = beam.J[3][2];  // ∂φ/∂θ₀
        dj.J[1][1] = beam.J[3][3];  // ∂φ/∂φ₀

        // Determinant
        dj.determinant = dj.J[0][0] * dj.J[1][1] - dj.J[0][1] * dj.J[1][0];

        return dj;
    }

private:
    const Sirius::AccretionDiskD* m_disk;
    Config m_config;
};

} // namespace sirius::kernel

#endif // SIRIUS_RENDER_KNSI001A_H
