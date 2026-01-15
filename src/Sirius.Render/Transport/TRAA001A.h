// KNAA001A.h - Anti-Aliasing via Beam Footprint
// Component ID: KNAA001A
// Purpose: Starfield sampling with Gaussian beam weighting for anti-aliasing
//
// MATHEMATICAL BASIS:
// A beam with semi-major axis σ_a and semi-minor axis σ_b defines an elliptical
// Gaussian filter. Stars within the filter footprint contribute weighted by:
//
//   w(Δα, Δδ) = exp(-0.5 × [(Δα/σ_a)² + (Δδ/σ_b)²])
//
// where (Δα, Δδ) are angular offsets from beam centre in aligned coordinates.
//
// REFERENCE: James et al. (2015) "DNGR" Section 4.2

#pragma once

#include "KNBI001A.h"
#include "../Sirius.Core/Symplectic/MTSB001A.h"
#include <cmath>
#include <vector>

namespace sirius::kernel {

using namespace sirius::math;
using namespace sirius::spectral;

//==============================================================================
// Star Catalog Entry
//==============================================================================

struct StarEntry {
    double alpha;       // Right ascension [rad]
    double delta;       // Declination [rad]
    double magnitude;   // Apparent magnitude
    double temperature; // Effective temperature [K]
    
    // Compute spectral radiance (flux per solid angle)
    SpectralRadiance spectrum() const {
        // Flux from magnitude: F = F₀ × 10^(-0.4 × m)
        // where F₀ is flux for m=0 star (Vega)
        double flux_ratio = std::pow(10.0, -0.4 * magnitude);
        
        // Blackbody spectrum scaled by flux
        SpectralRadiance bb = SpectralRadiance::blackbody(temperature);
        double scale = flux_ratio / bb.totalEnergy();
        return bb * scale;
    }
};

//==============================================================================
// BeamSampler: Sample starfield with beam filtering
//==============================================================================

class BeamSampler {
public:
    struct Config {
        double truncationRadius = 3.0;  // Truncate filter at 3σ
        double minWeight = 0.01;        // Minimum weight to include
    };
    
    explicit BeamSampler(const Config& config = Config()) : m_config(config) {}
    
    //--------------------------------------------------------------------------
    // Sample stars within beam footprint
    //--------------------------------------------------------------------------
    
    struct SampledStar {
        const StarEntry* star;
        double weight;
        SpectralRadiance contribution;
    };
    
    std::vector<SampledStar> sampleStarfield(
        const BeamStateD& beam,
        const std::vector<StarEntry>& catalog,
        double beamAlpha, double beamDelta  // Beam centre in sky coordinates
    ) const {
        std::vector<SampledStar> result;
        
        // Get beam ellipse parameters
        double sigma_a = beam.majorAxis;
        double sigma_b = beam.minorAxis;
        double pa = beam.orientation;  // Position angle
        
        // Truncation radius in angular space
        double truncation = m_config.truncationRadius;
        double max_radius = truncation * std::max(sigma_a, sigma_b);
        
        // Pre-compute rotation for aligned coordinates
        double cos_pa = std::cos(pa);
        double sin_pa = std::sin(pa);
        
        for (const auto& star : catalog) {
            // Angular offset from beam centre
            double dalpha = star.alpha - beamAlpha;
            double ddelta = star.delta - beamDelta;
            
            // Quick reject: bounding box
            double dist2 = dalpha*dalpha + ddelta*ddelta;
            if (dist2 > max_radius * max_radius) continue;
            
            // Rotate to ellipse-aligned coordinates
            double dx =  cos_pa * dalpha + sin_pa * ddelta;
            double dy = -sin_pa * dalpha + cos_pa * ddelta;
            
            // Gaussian weight
            double u = (sigma_a > 0) ? (dx / sigma_a) : 0;
            double v = (sigma_b > 0) ? (dy / sigma_b) : 0;
            double weight = std::exp(-0.5 * (u*u + v*v));
            
            if (weight < m_config.minWeight) continue;
            
            // Weighted contribution
            SampledStar ss;
            ss.star = &star;
            ss.weight = weight;
            ss.contribution = star.spectrum() * weight;
            
            result.push_back(ss);
        }
        
        return result;
    }
    
    //--------------------------------------------------------------------------
    // Integrate total contribution from sampled stars
    //--------------------------------------------------------------------------
    
    SpectralRadiance integrateContributions(const std::vector<SampledStar>& samples) const {
        SpectralRadiance total = SpectralRadiance::zero();
        double totalWeight = 0;
        
        for (const auto& s : samples) {
            total = total + s.contribution;
            totalWeight += s.weight;
        }
        
        // Normalise by total weight (so filter integrates to 1)
        if (totalWeight > 0) {
            double normFactor = 1.0 / totalWeight;
            total = total * normFactor;
        }
        
        return total;
    }
    
    //--------------------------------------------------------------------------
    // Compute effective solid angle of beam
    // Ω = π × σ_a × σ_b (area of ellipse in steradian)
    //--------------------------------------------------------------------------
    
    static double beamSolidAngle(const BeamStateD& beam) {
        return M_PI * beam.majorAxis * beam.minorAxis;
    }
    
private:
    Config m_config;
};

} // namespace sirius::kernel
