// KNSI001A.h - Scene Intersection and Spectral Integration
// Component ID: KNSI001A
// Purpose: CPU implementation of scene intersection (disk, celestial sphere)
//
// MATHEMATICAL BASIS:
// The beam intersects with scene elements (accretion disk, stars) and accumulates
// spectral radiance. The Jacobian determines the sampled area on each surface.
//
// For disk intersection:
//   J_disk = |det(∂(r,φ)/∂(θ₀,φ₀))| - Jacobian projected onto disk plane
//
// REFERENCE: James et al. (2015) "DNGR" Section 4

#pragma once

#include "KNBI001A.h"
#include "KNAA001A.h"
#include "../Sirius.Physics/Disk/PHAD001A.h"
#include "../Sirius.Math/MTSB001A.h"
#include <cmath>

namespace sirius::kernel {

using namespace sirius::physics;
using namespace sirius::spectral;

//==============================================================================
// SceneIntersector: Handle disk and celestial sphere intersections
//==============================================================================

class SceneIntersector {
public:
    struct Config {
        double escapeRadius = 1000;      // Celestial sphere radius [M]
        double diskThickness = 0.01;     // Disk half-thickness [rad from equator]
    };
    
    SceneIntersector(const AccretionDiskD* disk, const Config& config = Config())
        : m_disk(disk), m_config(config) {}
    
    //--------------------------------------------------------------------------
    // Disk Plane Jacobian Transform
    // Transforms 4×4 spacetime Jacobian to 2×2 disk tangent plane Jacobian
    //--------------------------------------------------------------------------
    
    struct DiskJacobian {
        double J[2][2];     // ∂(r,φ)/∂(θ₀,φ₀)
        double determinant;
        double area;        // Sampled area on disk
    };
    
    DiskJacobian transformJacobianToDisk(const BeamStateD& beam) const {
        DiskJacobian result;
        
        // The full 4×4 Jacobian maps initial conditions to current position:
        // J^i_j = ∂x^i/∂x^j_0 where x = (t, r, θ, φ)
        //
        // For disk intersection (θ = π/2), we need the 2×2 submatrix
        // that maps initial angular directions to disk coordinates (r, φ):
        //
        // J_disk[a][b] = ∂(r,φ)^a / ∂(θ₀,φ₀)^b
        //
        // This is extracted from rows 1,3 (r,φ) and columns 2,3 (θ₀,φ₀)
        
        result.J[0][0] = beam.J[1][2];  // ∂r/∂θ₀
        result.J[0][1] = beam.J[1][3];  // ∂r/∂φ₀
        result.J[1][0] = beam.J[3][2];  // ∂φ/∂θ₀
        result.J[1][1] = beam.J[3][3];  // ∂φ/∂φ₀
        
        // Determinant gives area scaling
        result.determinant = result.J[0][0] * result.J[1][1] - result.J[0][1] * result.J[1][0];
        
        // Sampled area = |det| × initial pixel solid angle × r² (for disk metric)
        double r = beam.x.r;
        result.area = std::abs(result.determinant) * beam.initialPixelSolidAngle * r * r;
        
        return result;
    }
    
    //--------------------------------------------------------------------------
    // Beam Area on Celestial Sphere
    //--------------------------------------------------------------------------
    
    double celestialSphereArea(const BeamStateD& beam) const {
        // At celestial sphere, angular Jacobian gives solid angle directly
        // Area = solid_angle × R² where R is sphere radius
        double R = m_config.escapeRadius;
        return beam.solidAngle * R * R;
    }
    
    //--------------------------------------------------------------------------
    // Full Scene Intersection
    //--------------------------------------------------------------------------
    
    struct IntersectionResult {
        enum Type { NONE, DISK, CELESTIAL_SPHERE, HORIZON };
        
        Type type = NONE;
        double r = 0;           // Radius at intersection
        double g_factor = 1.0;  // Gravitational redshift factor
        SpectralRadiance emission;
        double area = 0;        // Sampled area
    };
    
    IntersectionResult intersectScene(const BeamStateD& beam) const {
        IntersectionResult result;
        
        // Check termination conditions
        if (beam.terminated) {
            if (beam.x.r < m_disk->iscoRadius()) {
                result.type = IntersectionResult::HORIZON;
            }
            return result;
        }
        
        // Check celestial sphere (escaped)
        if (beam.x.r > m_config.escapeRadius) {
            result.type = IntersectionResult::CELESTIAL_SPHERE;
            result.r = m_config.escapeRadius;
            result.area = celestialSphereArea(beam);
            result.g_factor = 1.0;  // Far from BH
            // Celestial sphere emission would come from starfield
            return result;
        }
        
        // Check disk intersection
        double theta_equator = std::abs(beam.x.theta - M_PI/2);
        if (theta_equator < m_config.diskThickness) {
            double r = beam.x.r;
            if (r >= m_disk->innerRadius() && r <= m_disk->outerRadius()) {
                result.type = IntersectionResult::DISK;
                result.r = r;
                
                // Compute disk plane Jacobian
                auto diskJ = transformJacobianToDisk(beam);
                result.area = diskJ.area;
                
                // Gravitational redshift factor at disk
                // g = √(1 - 2M/r) for Schwarzschild (simplified)
                double M = 1.0;  // Mass in geometric units
                double g2 = 1.0 - 2*M/r;
                result.g_factor = (g2 > 0) ? std::sqrt(g2) : 0;
                
                // Disk emission with limb darkening
                SpectralRadiance local_emission = m_disk->emissionSpectrum(r);
                
                // Viewing angle (approximate: use θ deviation)
                double cos_incl = std::cos(theta_equator);
                double limb = AccretionDiskD::limbDarkening(cos_incl);
                
                // Apply limb darkening and redshift
                result.emission = local_emission * limb;
                result.emission = result.emission.applyRedshift(result.g_factor);
                
                return result;
            }
        }
        
        result.type = IntersectionResult::NONE;
        return result;
    }
    
    //--------------------------------------------------------------------------
    // Accumulate Spectral Contribution from Intersection
    //--------------------------------------------------------------------------
    
    SpectralRadiance accumulateContribution(const IntersectionResult& intersection) const {
        if (intersection.type == IntersectionResult::DISK) {
            // Disk emission already includes redshift and limb darkening
            return intersection.emission;
        }
        
        if (intersection.type == IntersectionResult::CELESTIAL_SPHERE) {
            // Would sample starfield here - return placeholder
            return SpectralRadiance::zero();
        }
        
        // Horizon or no intersection
        return SpectralRadiance::zero();
    }
    
private:
    const AccretionDiskD* m_disk;
    Config m_config;
};

} // namespace sirius::kernel
