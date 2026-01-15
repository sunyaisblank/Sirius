// KNST001A.h - Spectral Transport Module
// Component ID: KNST001A
// Purpose: Integrate spectral radiance along beam path with redshift
//
// MATHEMATICAL BASIS:
// The specific intensity I_ν transforms under Lorentz boost as:
//   I_ν/ν³ = invariant
//
// For gravitational redshift with factor g = ν_obs/ν_emit:
//   I_obs = g⁴ × I_emit (already in SpectralRadiance::applyRedshift)
//
// Total flux along beam path:
//   F = ∫ I(s) × Ω(s) ds
// where Ω(s) is the solid angle subtended by the beam at s.
//
// REFERENCE: James et al. (2015) "DNGR" Section 4

#pragma once

#include "KNBI001A.h"
#include "KNSI001A.h"
#include "KNAA001A.h"
#include "../Sirius.Core/Symplectic/MTSB001A.h"
#include <vector>

namespace sirius::kernel {

using namespace sirius::physics;
using namespace sirius::spectral;

//==============================================================================
// SpectralTransport: Integrate spectral radiance along beam
//==============================================================================

class SpectralTransport {
public:
    struct Config {
        double minContribution = 1e-6;  // Skip contributions below this
        int maxPathSegments = 1000;     // Maximum integration segments
    };
    
    SpectralTransport(
        const AccretionDiskD* disk,
        const std::vector<StarEntry>* starfield = nullptr,
        const Config& config = Config()
    ) : m_disk(disk), m_starfield(starfield), m_config(config) {}
    
    //--------------------------------------------------------------------------
    // Accumulated Spectral Result
    //--------------------------------------------------------------------------
    
    struct TransportResult {
        SpectralRadiance totalRadiance;  // Final accumulated radiance
        SpectralRadiance diskContrib;    // Disk contribution
        SpectralRadiance starContrib;    // Starfield contribution
        double minG = 1.0;               // Minimum g-factor encountered
        double maxG = 1.0;               // Maximum g-factor encountered
        int diskHits = 0;                // Number of disk intersections
        int starSamples = 0;             // Stars contributing
        bool hitHorizon = false;
    };
    
    //--------------------------------------------------------------------------
    // Integrate Spectral Radiance Along Pre-computed Beam Path
    //--------------------------------------------------------------------------
    
    TransportResult integrate(
        const std::vector<BeamStateD>& beamPath,
        const std::vector<double>& gFactors  // g-factor at each point
    ) const {
        TransportResult result;
        result.totalRadiance = SpectralRadiance::zero();
        result.diskContrib = SpectralRadiance::zero();
        result.starContrib = SpectralRadiance::zero();
        
        if (beamPath.empty()) return result;
        
        SceneIntersector intersector(m_disk);
        BeamSampler sampler;
        
        for (size_t i = 0; i < beamPath.size(); ++i) {
            const BeamStateD& beam = beamPath[i];
            double g = (i < gFactors.size()) ? gFactors[i] : 1.0;
            
            result.minG = std::min(result.minG, g);
            result.maxG = std::max(result.maxG, g);
            
            if (beam.terminated) {
                // Check if terminated at horizon
                if (beam.x.r < m_disk->iscoRadius()) {
                    result.hitHorizon = true;
                }
                break;
            }
            
            // Check for disk intersection
            auto intersection = intersector.intersectScene(beam);
            
            if (intersection.type == SceneIntersector::IntersectionResult::DISK) {
                // Disk emission with redshift
                SpectralRadiance diskEmission = m_disk->emissionSpectrum(intersection.r);
                
                // Apply limb darkening
                double cos_incl = std::abs(std::cos(beam.x.theta - M_PI/2));
                double limb = AccretionDiskD::limbDarkening(cos_incl);
                diskEmission = diskEmission * limb;
                
                // Apply redshift (g⁴ scaling + wavelength shift)
                diskEmission = diskEmission.applyRedshift(g);
                
                // Weight by beam solid angle
                double weight = beam.solidAngle / beam.initialPixelSolidAngle;
                if (weight > 0 && weight < 1e10) {
                    diskEmission = diskEmission * weight;
                }
                
                result.diskContrib = result.diskContrib + diskEmission;
                result.diskHits++;
            }
            
            // Check for celestial sphere (starfield)
            if (intersection.type == SceneIntersector::IntersectionResult::CELESTIAL_SPHERE) {
                if (m_starfield != nullptr) {
                    // Sample starfield with beam footprint
                    auto samples = sampler.sampleStarfield(
                        beam, *m_starfield, 
                        beam.x.phi, beam.x.theta - M_PI/2  // RA, Dec from coords
                    );
                    
                    auto starLight = sampler.integrateContributions(samples);
                    
                    // Redshift star light (usually g ≈ 1 at celestial sphere)
                    starLight = starLight.applyRedshift(g);
                    
                    result.starContrib = result.starContrib + starLight;
                    result.starSamples += static_cast<int>(samples.size());
                }
            }
        }
        
        // Combine contributions
        result.totalRadiance = result.diskContrib + result.starContrib;
        
        return result;
    }
    
    //--------------------------------------------------------------------------
    // Single-Point Emission (for simple integration)
    //--------------------------------------------------------------------------
    
    SpectralRadiance evaluateEmission(const BeamStateD& beam, double g) const {
        SceneIntersector::Config siConfig;
        SceneIntersector intersector(m_disk, siConfig);
        
        auto intersection = intersector.intersectScene(beam);
        
        if (intersection.type == SceneIntersector::IntersectionResult::DISK) {
            SpectralRadiance emission = intersection.emission;
            return emission.applyRedshift(g);
        }
        
        return SpectralRadiance::zero();
    }
    
    //--------------------------------------------------------------------------
    // Compute g-factor for Disk Point
    // g = √(1 - 2M/r) for Schwarzschild (ignores disk rotation)
    //--------------------------------------------------------------------------
    
    static double computeDiskGFactor(double r, double a = 0) {
        double M = 1.0;  // Geometric units
        
        // Simplified: ignoring disk angular velocity
        // Full formula: g = 1 / (u^t × (1 - bΩ))
        // where Ω is disk angular velocity
        
        if (a == 0) {
            // Schwarzschild
            double g2 = 1.0 - 2*M/r;
            return (g2 > 0) ? std::sqrt(g2) : 0;
        } else {
            // Kerr - approximate
            double r2 = r * r;
            double a2 = a * a;
            double Delta = r2 - 2*M*r + a2;
            double g2 = Delta / (r2 + a2 + 2*M*a2/r);
            return (g2 > 0) ? std::sqrt(g2) : 0;
        }
    }
    
private:
    const AccretionDiskD* m_disk;
    const std::vector<StarEntry>* m_starfield;
    Config m_config;
};

} // namespace sirius::kernel
