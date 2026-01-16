// KNAA001A.h - Beam Sampling / Anti-Aliasing
// Component ID: KNAA001A
// Beam-based starfield sampling with Gaussian PSF weighting
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_KNAA001A_H
#define SIRIUS_RENDER_KNAA001A_H

#include <MTTP001A.h>  // For Vec4d
#include <array>
#include <cmath>
#include <vector>

namespace sirius::kernel {

// Use Vec4d from MTTP001A.h (sirius::math namespace)
using sirius::math::Vec4d;
using sirius::math::GeodesicStateD;
using sirius::math::Mat4d;

//==============================================================================
// BeamStateD
// Full beam propagation state including position, Jacobian, and ellipse
//==============================================================================

struct BeamStateD {
    // Position in Boyer-Lindquist coordinates
    Vec4d x;

    // 4x4 Jacobian matrix for beam propagation
    double J[4][4] = {};

    // Beam ellipse parameters
    double majorAxis = 0.0;       // Semi-major axis in radians
    double minorAxis = 0.0;       // Semi-minor axis in radians
    double orientation = 0.0;     // Angle of major axis from reference direction

    // Solid angles
    double solidAngle = 1e-6;
    double initialPixelSolidAngle = 1e-6;

    // State flags
    bool terminated = false;      // Has beam reached termination condition?
};

//==============================================================================
// StarEntry
// A catalog star for beam sampling
//==============================================================================

struct StarEntry {
    double alpha = 0.0;      // Right ascension offset (radians)
    double delta = 0.0;      // Declination offset (radians)
    double magnitude = 5.0;  // Apparent magnitude
    double temperature = 6000.0;  // Effective temperature (K)
};

//==============================================================================
// StarSample
// Result of beam-weighted star sampling
//==============================================================================

struct StarSample {
    const StarEntry* star = nullptr;
    double weight = 0.0;     // Gaussian beam weight [0, 1]
    double flux = 0.0;       // Weighted flux contribution
};

//==============================================================================
// BeamSampler
// Samples starfield with beam-shaped PSF
//==============================================================================

class BeamSampler {
public:
    /// Compute beam solid angle (π × a × b for ellipse)
    static double beamSolidAngle(const BeamStateD& beam) {
        return M_PI * beam.majorAxis * beam.minorAxis;
    }

    /// Sample stars within beam footprint
    /// @param beam Current beam state (ellipse parameters)
    /// @param catalog Star catalog entries
    /// @param beamAlpha Beam centre right ascension
    /// @param beamDelta Beam centre declination
    /// @return Weighted star samples
    std::vector<StarSample> sampleStarfield(const BeamStateD& beam,
                                            const std::vector<StarEntry>& catalog,
                                            double beamAlpha,
                                            double beamDelta) const {
        std::vector<StarSample> samples;

        // Cutoff at 3σ
        double cutoffMajor = 3.0 * beam.majorAxis;
        double cutoffMinor = 3.0 * beam.minorAxis;
        double maxCutoff = std::max(cutoffMajor, cutoffMinor);

        double cosOrient = std::cos(-beam.orientation);
        double sinOrient = std::sin(-beam.orientation);

        for (const auto& star : catalog) {
            // Offset from beam centre
            double dalpha = star.alpha - beamAlpha;
            double ddelta = star.delta - beamDelta;

            // Quick distance check
            double distSq = dalpha * dalpha + ddelta * ddelta;
            if (distSq > maxCutoff * maxCutoff) continue;

            // Transform to beam ellipse coordinates
            double u = dalpha * cosOrient - ddelta * sinOrient;
            double v = dalpha * sinOrient + ddelta * cosOrient;

            // Normalized ellipse distance
            double normDistSq = 0.0;
            if (beam.majorAxis > 0.0 && beam.minorAxis > 0.0) {
                double uNorm = u / beam.majorAxis;
                double vNorm = v / beam.minorAxis;
                normDistSq = uNorm * uNorm + vNorm * vNorm;
            }

            // 3σ cutoff
            if (normDistSq > 9.0) continue;

            // Gaussian weight
            double weight = std::exp(-0.5 * normDistSq);

            StarSample sample;
            sample.star = &star;
            sample.weight = weight;
            sample.flux = weight * std::pow(10.0, -0.4 * star.magnitude);
            samples.push_back(sample);
        }

        return samples;
    }
};

} // namespace sirius::kernel

#endif // SIRIUS_RENDER_KNAA001A_H
