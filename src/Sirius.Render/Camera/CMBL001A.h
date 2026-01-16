// CMBL001A.h - Camera Model (Boyer-Lindquist)
// Component ID: CMBL001A
// Camera configuration for ray generation in Boyer-Lindquist coordinates
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_CMBL001A_H
#define SIRIUS_RENDER_CMBL001A_H

#include <array>
#include <cmath>

namespace sirius::render {

//==============================================================================
// CameraModel
// Configures virtual camera for ray generation
//==============================================================================

struct CameraModel {
    // Position in Boyer-Lindquist coordinates (r, theta, phi)
    double r = 50.0;       // Radial distance (in M units)
    double theta = M_PI_2; // Polar angle (pi/2 = equatorial plane)
    double phi = 0.0;      // Azimuthal angle

    // Orientation (angles relative to local FIDO frame)
    double inclination = 0.0; // Camera tilt (0 = radially inward)
    double roll = 0.0;        // Camera roll

    // Lens parameters
    double fov = 90.0;        // Field of view (degrees)
    double focalLength = 1.0; // Normalized focal length

    // Compute ray direction for pixel (u, v) normalized to [-1, 1]
    std::array<double, 3> rayDirection(double u, double v) const {
        double tanHalfFov = std::tan(fov * M_PI / 360.0);
        double x = u * tanHalfFov;
        double y = v * tanHalfFov;
        double z = -1.0;  // Looking radially inward

        // Normalize
        double len = std::sqrt(x * x + y * y + z * z);
        return {x / len, y / len, z / len};
    }
};

} // namespace sirius::render

#endif // SIRIUS_RENDER_CMBL001A_H
