// CMBS001A.h - Camera Interface and Lens Models
// Component ID: CMBS001A (Camera/Base Interface)
//
// Abstract interface for camera systems. All camera implementations
// (Pinhole, ThinLens, Fisheye, etc.) derive from ICamera.
//
// Architecture: Strategy pattern for interchangeable lens models.
// Each camera generates rays for given pixel coordinates.

#pragma once

#include <MTTN001A.h>
#include <cmath>
#include <memory>
#include <string>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Sirius {

//==============================================================================
// Ray Structure (for camera output)
//==============================================================================
struct CameraRay {
    Vec4 origin;       ///< Ray origin (observer position)
    Vec4 direction;    ///< Ray direction (unit 4-vector)
    float weight;      ///< Ray weight (for importance sampling)
    
    CameraRay() : weight(1.0f) {}
};

//==============================================================================
// Camera Configuration
//==============================================================================
struct CameraConfig {
    // Position (Boyer-Lindquist coordinates: t, r, θ, φ)
    double t = 0.0;
    double r = 50.0;           ///< Distance from centre (in M)
    double theta = M_PI / 2.0; ///< Polar angle (π/2 = equatorial)
    double phi = 0.0;          ///< Azimuthal angle
    
    // Orientation
    float yaw = 0.0f;          ///< Yaw rotation (radians)
    float pitch = 0.0f;        ///< Pitch rotation (radians)
    float roll = 0.0f;         ///< Roll rotation (radians)
    
    // Lens properties
    float fov = 60.0f;         ///< Field of view (degrees)
    float focalLength = 50.0f; ///< Focal length (mm, for ThinLens)
    float aperture = 2.8f;     ///< f-number (for depth of field)
    float focusDistance = 50.0f; ///< Focus distance (M)
    
    // Image properties
    int width = 1920;
    int height = 1080;
};

//==============================================================================
// Lens Type Enumeration
//==============================================================================
enum class LensType {
    Pinhole,       ///< Ideal pinhole (infinite depth of field)
    ThinLens,      ///< Thin lens model with depth of field
    Fisheye,       ///< Equidistant fisheye projection
    Orthographic,  ///< Orthographic projection
    Panoramic      ///< 360° equirectangular
};

//==============================================================================
// ICamera Interface
//==============================================================================
class ICamera {
public:
    virtual ~ICamera() = default;
    
    /// @brief Generate a ray for the given pixel coordinates
    /// @param x Pixel x-coordinate (0 to width-1)
    /// @param y Pixel y-coordinate (0 to height-1)
    /// @param u Sample offset within pixel [0, 1)
    /// @param v Sample offset within pixel [0, 1)
    /// @return CameraRay with origin and direction
    virtual CameraRay generateRay(int x, int y, float u = 0.5f, float v = 0.5f) const = 0;
    
    /// @brief Get camera type
    virtual LensType getLensType() const = 0;
    
    /// @brief Get camera name
    virtual const char* getName() const = 0;
    
    /// @brief Get current configuration
    virtual const CameraConfig& getConfig() const = 0;
    
    /// @brief Update configuration
    virtual void setConfig(const CameraConfig& config) = 0;
    
    /// @brief Get position as Vec4
    Vec4 getPosition() const {
        const auto& cfg = getConfig();
        Vec4 pos;
        pos(0) = cfg.t;
        pos(1) = cfg.r;
        pos(2) = cfg.theta;
        pos(3) = cfg.phi;
        return pos;
    }
};

//==============================================================================
// PinholeCamera - Ideal perspective projection
//==============================================================================
class PinholeCamera : public ICamera {
public:
    explicit PinholeCamera(const CameraConfig& config = CameraConfig())
        : m_Config(config) {
        updateInternals();
    }
    
    CameraRay generateRay(int x, int y, float u = 0.5f, float v = 0.5f) const override {
        CameraRay ray;
        
        // Set origin
        ray.origin(0) = m_Config.t;
        ray.origin(1) = m_Config.r;
        ray.origin(2) = m_Config.theta;
        ray.origin(3) = m_Config.phi;
        
        // Normalised device coordinates (-1 to 1)
        float px = (2.0f * (x + u) / m_Config.width - 1.0f) * m_AspectRatio;
        float py = 1.0f - 2.0f * (y + v) / m_Config.height;
        
        // Direction in camera space (looking along -Z)
        float dx = px * m_TanHalfFov;
        float dy = py * m_TanHalfFov;
        float dz = -1.0f;
        
        // Apply rotation (yaw around Y, pitch around X, roll around Z)
        float cosY = std::cos(m_Config.yaw), sinY = std::sin(m_Config.yaw);
        float cosP = std::cos(m_Config.pitch), sinP = std::sin(m_Config.pitch);
        float cosR = std::cos(m_Config.roll), sinR = std::sin(m_Config.roll);
        
        // Roll
        float rx = dx * cosR - dy * sinR;
        float ry = dx * sinR + dy * cosR;
        
        // Pitch
        float ry2 = ry * cosP - dz * sinP;
        float rz = ry * sinP + dz * cosP;
        
        // Yaw
        float rx2 = rx * cosY + rz * sinY;
        float rz2 = -rx * sinY + rz * cosY;
        
        // Normalise
        float len = std::sqrt(rx2 * rx2 + ry2 * ry2 + rz2 * rz2);

        // =======================================================================
        // Map camera-local direction to spherical velocity components:
        //
        // Camera convention at (r, θ, φ):
        //   -Z (forward) → radial inward (-r direction)
        //   +X (right)   → +φ direction (increasing azimuth)
        //   +Y (up)      → -θ direction (toward North pole, decreasing θ)
        //
        // At camera position, global Cartesian frame:
        //   vx_global = dz (camera -Z maps to -x for observer at +x axis)
        //   vy_global = dx (camera +X maps to +y)
        //   vz_global = dy (camera +Y maps to +z toward North)
        //
        // Spherical velocity from Cartesian (at θ=π/2, φ=0):
        //   vr = vx_global
        //   vθ = -vz_global / r  →  r·vθ = -vz_global = -dy
        //   vφ = vy_global / (r·sinθ)  →  r·sinθ·vφ = vy_global = dx
        //
        // The tracer expects: dir(1)=vr, dir(2)=r·vθ, dir(3)=r·sinθ·vφ
        // =======================================================================
        ray.direction(0) = 0.0f;   // dt/dλ (set by geodesic normalisation)
        ray.direction(1) = rz2 / len;   // vr = vx_global (toward -r)
        ray.direction(2) = -ry2 / len;  // r·vθ = -vz_global (minus sign: +Y → -θ)
        ray.direction(3) = rx2 / len;   // r·sinθ·vφ = vy_global
        
        ray.weight = 1.0f;
        return ray;
    }
    
    LensType getLensType() const override { return LensType::Pinhole; }
    const char* getName() const override { return "Pinhole Camera"; }
    const CameraConfig& getConfig() const override { return m_Config; }
    
    void setConfig(const CameraConfig& config) override {
        m_Config = config;
        updateInternals();
    }
    
private:
    void updateInternals() {
        m_AspectRatio = static_cast<float>(m_Config.width) / m_Config.height;
        m_TanHalfFov = std::tan(m_Config.fov * static_cast<float>(M_PI) / 360.0f);
    }
    
    CameraConfig m_Config;
    float m_AspectRatio = 1.0f;
    float m_TanHalfFov = 1.0f;
};

//==============================================================================
// ThinLensCamera - Depth of field simulation
//==============================================================================
class ThinLensCamera : public ICamera {
public:
    explicit ThinLensCamera(const CameraConfig& config = CameraConfig())
        : m_Config(config) {
        updateInternals();
    }
    
    CameraRay generateRay(int x, int y, float u = 0.5f, float v = 0.5f) const override {
        CameraRay ray;
        
        ray.origin(0) = m_Config.t;
        ray.origin(1) = m_Config.r;
        ray.origin(2) = m_Config.theta;
        ray.origin(3) = m_Config.phi;
        
        // Point on image plane
        float px = (2.0f * (x + u) / m_Config.width - 1.0f) * m_AspectRatio;
        float py = 1.0f - 2.0f * (y + v) / m_Config.height;
        
        // Point on focus plane
        float focusX = px * m_Config.focusDistance / m_Config.focalLength;
        float focusY = py * m_Config.focusDistance / m_Config.focalLength;
        float focusZ = -m_Config.focusDistance;
        
        // Random point on aperture (for depth of field)
        float apertureRadius = m_Config.focalLength / (2.0f * m_Config.aperture);
        // Use deterministic jitter based on u, v
        float theta = 2.0f * static_cast<float>(M_PI) * u;
        float r_lens = apertureRadius * std::sqrt(v);
        float lensX = r_lens * std::cos(theta);
        float lensY = r_lens * std::sin(theta);
        
        // Direction from lens point to focus point
        float dx = focusX - lensX;
        float dy = focusY - lensY;
        float dz = focusZ;
        
        float len = std::sqrt(dx * dx + dy * dy + dz * dz);

        // Same coordinate mapping as PinholeCamera (see detailed comments there)
        ray.direction(0) = 0.0f;
        ray.direction(1) = dz / len;   // vr
        ray.direction(2) = -dy / len;  // r·vθ (minus sign: +Y → -θ)
        ray.direction(3) = dx / len;   // r·sinθ·vφ

        ray.weight = 1.0f;
        return ray;
    }

    LensType getLensType() const override { return LensType::ThinLens; }
    const char* getName() const override { return "Thin Lens Camera"; }
    const CameraConfig& getConfig() const override { return m_Config; }
    
    void setConfig(const CameraConfig& config) override {
        m_Config = config;
        updateInternals();
    }
    
private:
    void updateInternals() {
        m_AspectRatio = static_cast<float>(m_Config.width) / m_Config.height;
    }
    
    CameraConfig m_Config;
    float m_AspectRatio = 1.0f;
};

//==============================================================================
// FisheyeCamera - Equidistant fisheye projection
//==============================================================================
class FisheyeCamera : public ICamera {
public:
    explicit FisheyeCamera(const CameraConfig& config = CameraConfig())
        : m_Config(config) {
        updateInternals();
    }
    
    CameraRay generateRay(int x, int y, float u = 0.5f, float v = 0.5f) const override {
        CameraRay ray;
        
        ray.origin(0) = m_Config.t;
        ray.origin(1) = m_Config.r;
        ray.origin(2) = m_Config.theta;
        ray.origin(3) = m_Config.phi;
        
        // Normalised coordinates from centre
        float px = (2.0f * (x + u) / m_Config.width - 1.0f) * m_AspectRatio;
        float py = 1.0f - 2.0f * (y + v) / m_Config.height;
        
        // Radial distance from centre
        float r_img = std::sqrt(px * px + py * py);
        float phi_img = std::atan2(py, px);
        
        // Equidistant projection: θ = r * (FOV/2)
        float theta_ray = r_img * m_Config.fov * static_cast<float>(M_PI) / 360.0f;
        
        // Clamp to hemisphere
        if (theta_ray > static_cast<float>(M_PI)) {
            ray.weight = 0.0f;
            return ray;
        }
        
        // Direction
        float sinT = std::sin(theta_ray);
        float cosT = std::cos(theta_ray);
        
        float dx = sinT * std::cos(phi_img);
        float dy = sinT * std::sin(phi_img);
        float dz = -cosT;

        // Same coordinate mapping as PinholeCamera (see detailed comments there)
        ray.direction(0) = 0.0f;
        ray.direction(1) = dz;    // vr
        ray.direction(2) = -dy;   // r·vθ (minus sign: +Y → -θ)
        ray.direction(3) = dx;    // r·sinθ·vφ
        
        ray.weight = 1.0f;
        return ray;
    }
    
    LensType getLensType() const override { return LensType::Fisheye; }
    const char* getName() const override { return "Fisheye Camera"; }
    const CameraConfig& getConfig() const override { return m_Config; }
    
    void setConfig(const CameraConfig& config) override {
        m_Config = config;
        updateInternals();
    }
    
private:
    void updateInternals() {
        m_AspectRatio = static_cast<float>(m_Config.width) / m_Config.height;
    }
    
    CameraConfig m_Config;
    float m_AspectRatio = 1.0f;
};

//==============================================================================
// Camera Factory
//==============================================================================
inline std::unique_ptr<ICamera> createCamera(LensType type, const CameraConfig& config = CameraConfig()) {
    switch (type) {
        case LensType::Pinhole:
            return std::make_unique<PinholeCamera>(config);
        case LensType::ThinLens:
            return std::make_unique<ThinLensCamera>(config);
        case LensType::Fisheye:
            return std::make_unique<FisheyeCamera>(config);
        default:
            return std::make_unique<PinholeCamera>(config);
    }
}

} // namespace Sirius
