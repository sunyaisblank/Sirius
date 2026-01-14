// CRGC001A.h - Geodesic Camera Controller
//
// Camera motion along geodesic paths in curved spacetime.
// Follows timelike geodesics with proper time, dilation, and frame-dragging.
// Orbit types: Circular, ISCO, Plunging, Freefall

#pragma once

#include <cmath>

namespace Sirius {

//==============================================================================
// Local vector types for geodesic calculations
//==============================================================================
struct Vec3 {
    float x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}
};

struct Vec4 {
    float t, r, theta, phi;
    Vec4() : t(0), r(0), theta(0), phi(0) {}
    Vec4(float t_, float r_, float th_, float ph_) : t(t_), r(r_), theta(th_), phi(ph_) {}
    
    Vec4 operator+(const Vec4& o) const { return Vec4(t+o.t, r+o.r, theta+o.theta, phi+o.phi); }
    Vec4 operator*(float s) const { return Vec4(t*s, r*s, theta*s, phi*s); }
};


//==============================================================================
// Geodesic Camera Mode
//==============================================================================
enum class CameraMode : int {
    FreeWASD = 0,         // Standard non-geodesic movement (default)
    GeodesicFreefall = 1, // Follow geodesic from current state
    CircularOrbit = 2,    // Stable circular orbit at given r
    ISCO = 3,             // Innermost stable circular orbit
    PlungingOrbit = 4     // Radial infall from given r
};

//==============================================================================
// Geodesic Camera State
// Represents the camera's position and velocity in Boyer-Lindquist coordinates
//==============================================================================
struct GeodesicCameraState {
    Vec4 position;       // (t, r, θ, φ) in Boyer-Lindquist
    Vec4 velocity;       // 4-velocity u^μ = dx^μ/dτ
    float properTime;    // τ accumulated since orbit start
    float coordinateTime;// t accumulated
    float gamma;         // Lorentz factor relative to FIDO
    
    GeodesicCameraState() : properTime(0.0f), coordinateTime(0.0f), gamma(1.0f) {
        position = Vec4(0.0f, 50.0f, 1.5708f, 0.0f);  // Default: r=50M, θ=π/2
        velocity = Vec4(1.0f, 0.0f, 0.0f, 0.0f);      // At rest initially
    }
};

//==============================================================================
// GeodesicCamera Class
// Manages camera motion along geodesic paths in curved spacetime
//==============================================================================
class GeodesicCamera {
public:
    GeodesicCamera();
    ~GeodesicCamera() = default;
    
    // =========================================================================
    // Initialization
    // =========================================================================
    
    /// Initialize camera at given position
    void initialize(const Vec4& position, float M, float a);
    
    /// Set metric parameters (mass and spin)
    void setMetricParams(float M, float a);
    
    // =========================================================================
    // Orbit Presets
    // =========================================================================
    
    /// Set camera to circular equatorial orbit at radius r
    /// @param r Orbital radius in units of M
    /// @param prograde True for prograde orbit, false for retrograde
    void setCircularOrbit(float r, bool prograde = true);
    
    /// Set camera to ISCO for current metric parameters
    void setISCO(bool prograde = true);
    
    /// Set camera to plunging orbit from given radius
    /// @param r_start Starting radius in units of M
    void setPlungingOrbit(float r_start);
    
    /// Enable/disable geodesic motion
    void setMode(CameraMode mode);
    CameraMode getMode() const { return m_Mode; }
    
    // =========================================================================
    // Time Evolution
    // =========================================================================
    
    /// Advance camera along geodesic by coordinate time dt
    /// @param dt Coordinate time step
    void advanceByCoordinateTime(float dt);
    
    /// Advance camera along geodesic by proper time dτ
    /// @param dtau Proper time step
    void advanceByProperTime(float dtau);
    
    /// Reset proper time counter
    void resetProperTime();
    
    // =========================================================================
    // State Access
    // =========================================================================
    
    /// Get current camera state
    const GeodesicCameraState& getState() const { return m_State; }
    
    /// Get current position in Boyer-Lindquist
    Vec4 getPosition() const { return m_State.position; }
    
    /// Get current 4-velocity
    Vec4 getVelocity() const { return m_State.velocity; }
    
    /// Get accumulated proper time
    float getProperTime() const { return m_State.properTime; }
    
    /// Get accumulated coordinate time
    float getCoordinateTime() const { return m_State.coordinateTime; }
    
    /// Get time dilation factor (dτ/dt)
    float getTimeDilationFactor() const;
    
    /// Get orbital angular velocity for circular orbits
    float getOrbitalAngularVelocity() const;
    
    // =========================================================================
    // Conversion
    // =========================================================================
    
    /// Convert current state to Cartesian position for rendering
    Vec3 getCartesianPosition() const;
    
    /// Get look direction in local frame
    Vec3 getLookDirection() const;
    
private:
    // =========================================================================
    // Internal Methods
    // =========================================================================
    
    /// RK45/DOPRI adaptive step for geodesic integration
    /// @param dt Step size (modified adaptively)
    void rk45Step(float& dt);
    
    /// Compute geodesic acceleration: a^μ = -Γ^μ_αβ u^α u^β
    Vec4 computeAcceleration(const Vec4& pos, const Vec4& vel);
    
    /// Normalize 4-velocity to satisfy g_μν u^μ u^ν = -1
    void normalizeVelocity();
    
    /// Compute ISCO radius for current metric
    float computeISCO(bool prograde) const;
    
    /// Compute circular orbit velocity at radius r
    Vec4 computeCircularVelocity(float r, bool prograde) const;
    
    // =========================================================================
    // Member Variables
    // =========================================================================
    GeodesicCameraState m_State;
    CameraMode m_Mode;
    
    // Metric parameters
    float m_M;  // Black hole mass
    float m_a;  // Spin parameter (|a| <= M)
    
    // Orbit parameters
    float m_OrbitRadius;  // For circular orbits
    bool m_Prograde;      // Orbit direction
    
    // Adaptive integration
    float m_StepSize;     // Current RK45 step size (adapts automatically)
};

} // namespace Sirius
