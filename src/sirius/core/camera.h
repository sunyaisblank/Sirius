#pragma once

// Camera interface and lens models (pinhole, thin-lens, fisheye) that generate
// rays in Boyer-Lindquist coordinates for the geodesic tracer. Strategy pattern
// over ICamera for interchangeable lens models. Ported from CMBS001A.h.

#include "sirius/base/contracts.h"
#include "sirius/core/tensor.h"

#include <cmath>
#include <memory>
#include <numbers>

namespace sirius::core {

// Ray emitted by a camera: observer position and unit 4-direction.
struct CameraRay {
    Vec4 origin;     // Ray origin (observer position)
    Vec4 direction;  // Ray direction (unit 4-vector)
    float weight;    // Ray weight (for importance sampling)

    CameraRay() : weight(1.0f) {}
};

// Observer placement, orientation, and lens/image properties.
struct CameraConfig {
    // Position (Boyer-Lindquist coordinates: t, r, theta, phi)
    double t = 0.0;
    double r = 50.0;                        // Distance from centre (in M)
    double theta = std::numbers::pi / 2.0;  // Polar angle (pi/2 = equatorial)
    double phi = 0.0;                       // Azimuthal angle

    // Orientation
    float yaw = 0.0f;    // Yaw rotation (radians)
    float pitch = 0.0f;  // Pitch rotation (radians)
    float roll = 0.0f;   // Roll rotation (radians)

    // Camera four-velocity: the spatial beta = v/c of the observer worldline,
    // expressed in the local ray-direction component frame (the same orthonormal
    // basis as CameraRay::direction indices 1,2,3, with index 1 the forward/view
    // axis). Special-relativistic aberration of the generated ray is applied in
    // that frame (P5); finite |beta| < 1 is a fail-closed precondition. Default
    // zero leaves every generated ray untouched (byte-pin). Reference: DNGR
    // (James et al. 2015, CQG 32 065001), arbitrary-worldline camera; MTW eq 2.29.
    double beta_x = 0.0;
    double beta_y = 0.0;
    double beta_z = 0.0;

    // Lens properties
    float fov = 60.0f;             // Field of view (degrees)
    float focal_length = 50.0f;    // Focal length (mm-equivalent, for ThinLens)
    float aperture = 2.8f;         // f-number (for depth of field)
    float focus_distance = 50.0f;  // Focus distance (M)

    // Image properties
    int width = 1920;
    int height = 1080;
};

// Lens projection models.
enum class LensType {
    Pinhole,   // Ideal pinhole (infinite depth of field)
    ThinLens,  // Thin lens model with depth of field
    Fisheye    // Equidistant fisheye projection
};

// Special-relativistic aberration of a generated ray in the camera's local
// frame. Given a camera moving with three-velocity beta (in the ray-component
// orthonormal basis), a photon direction n transforms so that, along the boost
// axis, cos theta' = (cos theta - beta) / (1 - beta cos theta) (MTW eq 2.29;
// the light-cone form of the Lorentz boost of the null direction). The
// perpendicular components scale by 1 / (gamma (1 - beta.n)), which keeps the
// result a unit direction. beta = 0 returns the ray untouched with no arithmetic,
// so the default camera path is byte-identical.
inline void AberrateRay(CameraRay& ray, double beta_x, double beta_y, double beta_z) {
    double b2 = beta_x * beta_x + beta_y * beta_y + beta_z * beta_z;
    if (b2 <= 0.0) return;  // Exact no-op: preserves the pinned render.

    // Operator input is rejected by ConfigLoader before session construction.
    // Keep the lower-level typed API fail-closed as well: silently projecting
    // an invalid worldline onto a different, near-luminal one changes the
    // requested camera.
    SIRIUS_PRE(std::isfinite(b2) && b2 < 1.0);

    double b = std::sqrt(b2);
    double ux = beta_x / b, uy = beta_y / b, uz = beta_z / b;  // Boost unit axis.

    double dx = ray.direction(1), dy = ray.direction(2), dz = ray.direction(3);
    double dlen = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (dlen < 1e-30) return;
    dx /= dlen;
    dy /= dlen;
    dz /= dlen;

    double n_par = dx * ux + dy * uy + dz * uz;  // cos theta about the boost axis.
    double px = dx - n_par * ux, py = dy - n_par * uy, pz = dz - n_par * uz;

    double gamma = 1.0 / std::sqrt(1.0 - b2);
    double denom = 1.0 - b * n_par;
    double n_par_prime = (n_par - b) / denom;
    double perp_scale = 1.0 / (gamma * denom);

    double ax = n_par_prime * ux + px * perp_scale;
    double ay = n_par_prime * uy + py * perp_scale;
    double az = n_par_prime * uz + pz * perp_scale;
    double alen = std::sqrt(ax * ax + ay * ay + az * az);
    if (alen > 1e-30) {
        ax /= alen;
        ay /= alen;
        az /= alen;
    }

    ray.direction(1) = ax;
    ray.direction(2) = ay;
    ray.direction(3) = az;
}

// Abstract camera: generates a ray per pixel sample.
class ICamera {
  public:
    virtual ~ICamera() = default;

    // Generate a ray for pixel (x, y) with sample offset (u, v) in [0, 1).
    virtual CameraRay GenerateRay(int x, int y, float u = 0.5f, float v = 0.5f) const = 0;

    // Generate a ray and apply the camera-worldline aberration from the config
    // (P5). Composes over every lens model without each model re-implementing the
    // boost. At beta = 0 this is exactly GenerateRay, so the default path and the
    // pinned render are unchanged.
    CameraRay GenerateRayAberrated(int x, int y, float u = 0.5f, float v = 0.5f) const {
        CameraRay ray = GenerateRay(x, y, u, v);
        const auto& cfg = GetConfig();
        AberrateRay(ray, cfg.beta_x, cfg.beta_y, cfg.beta_z);
        return ray;
    }

    virtual LensType GetLensType() const = 0;
    virtual const char* GetName() const = 0;
    virtual const CameraConfig& GetConfig() const = 0;
    virtual void SetConfig(const CameraConfig& config) = 0;

    // Observer position as a 4-vector (t, r, theta, phi).
    Vec4 GetPosition() const {
        const auto& cfg = GetConfig();
        Vec4 pos;
        pos(0) = cfg.t;
        pos(1) = cfg.r;
        pos(2) = cfg.theta;
        pos(3) = cfg.phi;
        return pos;
    }
};

// Ideal perspective projection (infinite depth of field).
class PinholeCamera : public ICamera {
  public:
    explicit PinholeCamera(const CameraConfig& config = CameraConfig()) : config_(config) {
        UpdateInternals();
    }

    CameraRay GenerateRay(int x, int y, float u = 0.5f, float v = 0.5f) const override {
        CameraRay ray;

        // Set origin
        ray.origin(0) = config_.t;
        ray.origin(1) = config_.r;
        ray.origin(2) = config_.theta;
        ray.origin(3) = config_.phi;

        // Normalised device coordinates (-1 to 1)
        float px = (2.0f * (x + u) / config_.width - 1.0f) * aspect_ratio_;
        float py = 1.0f - 2.0f * (y + v) / config_.height;

        // Direction in camera space (looking along -Z)
        float dx = px * tan_half_fov_;
        float dy = py * tan_half_fov_;
        float dz = -1.0f;

        // Apply rotation (yaw around Y, pitch around X, roll around Z)
        float cos_y = std::cos(config_.yaw), sin_y = std::sin(config_.yaw);
        float cos_p = std::cos(config_.pitch), sin_p = std::sin(config_.pitch);
        float cos_r = std::cos(config_.roll), sin_r = std::sin(config_.roll);

        // Roll
        float rx = dx * cos_r - dy * sin_r;
        float ry = dx * sin_r + dy * cos_r;

        // Pitch
        float ry2 = ry * cos_p - dz * sin_p;
        float rz = ry * sin_p + dz * cos_p;

        // Yaw
        float rx2 = rx * cos_y + rz * sin_y;
        float rz2 = -rx * sin_y + rz * cos_y;

        // Normalise
        float len = std::sqrt(rx2 * rx2 + ry2 * ry2 + rz2 * rz2);

        // =======================================================================
        // Map camera-local direction to spherical velocity components:
        //
        // Camera convention at (r, theta, phi):
        //   -Z (forward) -> radial inward (-r direction)
        //   +X (right)   -> +phi direction (increasing azimuth)
        //   +Y (up)      -> -theta direction (toward North pole, decreasing theta)
        //
        // At camera position, global Cartesian frame:
        //   vx_global = dz (camera -Z maps to -x for observer at +x axis)
        //   vy_global = dx (camera +X maps to +y)
        //   vz_global = dy (camera +Y maps to +z toward North)
        //
        // Spherical velocity from Cartesian (at theta=pi/2, phi=0):
        //   vr = vx_global
        //   vtheta = -vz_global / r  ->  r*vtheta = -vz_global = -dy
        //   vphi = vy_global / (r*sin theta)  ->  r*sin theta*vphi = vy_global = dx
        //
        // The tracer expects: dir(1)=vr, dir(2)=r*vtheta, dir(3)=r*sin theta*vphi
        // =======================================================================
        ray.direction(0) = 0.0f;        // dt/dlambda (set by geodesic normalisation)
        ray.direction(1) = rz2 / len;   // vr = vx_global (toward -r)
        ray.direction(2) = -ry2 / len;  // r*vtheta = -vz_global (minus sign: +Y -> -theta)
        ray.direction(3) = rx2 / len;   // r*sin theta*vphi = vy_global

        ray.weight = 1.0f;
        return ray;
    }

    LensType GetLensType() const override { return LensType::Pinhole; }
    const char* GetName() const override { return "Pinhole Camera"; }
    const CameraConfig& GetConfig() const override { return config_; }

    void SetConfig(const CameraConfig& config) override {
        config_ = config;
        UpdateInternals();
    }

  private:
    void UpdateInternals() {
        aspect_ratio_ = static_cast<float>(config_.width) / config_.height;
        tan_half_fov_ = std::tan(config_.fov * static_cast<float>(std::numbers::pi) / 360.0f);
    }

    CameraConfig config_;
    float aspect_ratio_ = 1.0f;
    float tan_half_fov_ = 1.0f;
};

// Thin-lens model with depth-of-field defocus.
class ThinLensCamera : public ICamera {
  public:
    explicit ThinLensCamera(const CameraConfig& config = CameraConfig()) : config_(config) {
        UpdateInternals();
    }

    CameraRay GenerateRay(int x, int y, float u = 0.5f, float v = 0.5f) const override {
        CameraRay ray;

        ray.origin(0) = config_.t;
        ray.origin(1) = config_.r;
        ray.origin(2) = config_.theta;
        ray.origin(3) = config_.phi;

        // Point on image plane
        float px = (2.0f * (x + u) / config_.width - 1.0f) * aspect_ratio_;
        float py = 1.0f - 2.0f * (y + v) / config_.height;

        // Point on the focus plane. A ray through the centre of the aperture
        // must have exactly the same perspective projection as the pinhole
        // camera; the previous focal-length division silently collapsed a 60
        // degree request to roughly two degrees at the default 50 mm setting.
        float focus_x = px * tan_half_fov_ * config_.focus_distance;
        float focus_y = py * tan_half_fov_ * config_.focus_distance;
        float focus_z = -config_.focus_distance;

        // Random point on the aperture (for depth of field). The spacetime uses
        // geometric M units and has no physical mass-to-millimetre scale, so a
        // 50 mm-equivalent lens defines one virtual lens unit. This preserves
        // the f-number and focal-length controls without mixing raw millimetres
        // directly with a focus distance measured in M.
        constexpr float kReferenceFocalLengthMillimetres = 50.0f;
        float aperture_radius =
            (config_.focal_length / kReferenceFocalLengthMillimetres) / (2.0f * config_.aperture);
        // Use deterministic jitter based on u, v
        float theta = 2.0f * static_cast<float>(std::numbers::pi) * u;
        float r_lens = aperture_radius * std::sqrt(v);
        float lens_x = r_lens * std::cos(theta);
        float lens_y = r_lens * std::sin(theta);

        // Direction from lens point to focus point
        float dx = focus_x - lens_x;
        float dy = focus_y - lens_y;
        float dz = focus_z;

        float len = std::sqrt(dx * dx + dy * dy + dz * dz);

        // Same coordinate mapping as PinholeCamera (see detailed comments there)
        ray.direction(0) = 0.0f;
        ray.direction(1) = dz / len;   // vr
        ray.direction(2) = -dy / len;  // r*vtheta (minus sign: +Y -> -theta)
        ray.direction(3) = dx / len;   // r*sin theta*vphi

        ray.weight = 1.0f;
        return ray;
    }

    LensType GetLensType() const override { return LensType::ThinLens; }
    const char* GetName() const override { return "Thin Lens Camera"; }
    const CameraConfig& GetConfig() const override { return config_; }

    void SetConfig(const CameraConfig& config) override {
        config_ = config;
        UpdateInternals();
    }

  private:
    void UpdateInternals() {
        aspect_ratio_ = static_cast<float>(config_.width) / config_.height;
        tan_half_fov_ = std::tan(config_.fov * static_cast<float>(std::numbers::pi) / 360.0f);
    }

    CameraConfig config_;
    float aspect_ratio_ = 1.0f;
    float tan_half_fov_ = 1.0f;
};

// Equidistant fisheye projection (angle proportional to image radius).
class FisheyeCamera : public ICamera {
  public:
    explicit FisheyeCamera(const CameraConfig& config = CameraConfig()) : config_(config) {
        UpdateInternals();
    }

    CameraRay GenerateRay(int x, int y, float u = 0.5f, float v = 0.5f) const override {
        CameraRay ray;

        ray.origin(0) = config_.t;
        ray.origin(1) = config_.r;
        ray.origin(2) = config_.theta;
        ray.origin(3) = config_.phi;

        // Normalised coordinates from centre
        float px = (2.0f * (x + u) / config_.width - 1.0f) * aspect_ratio_;
        float py = 1.0f - 2.0f * (y + v) / config_.height;

        // Radial distance from centre
        float r_img = std::sqrt(px * px + py * py);
        float phi_img = std::atan2(py, px);

        // Equidistant projection: theta = r * (FOV/2)
        float theta_ray = r_img * config_.fov * static_cast<float>(std::numbers::pi) / 360.0f;

        // Clamp to hemisphere
        if (theta_ray > static_cast<float>(std::numbers::pi)) {
            ray.weight = 0.0f;
            return ray;
        }

        // Direction
        float sin_t = std::sin(theta_ray);
        float cos_t = std::cos(theta_ray);

        float dx = sin_t * std::cos(phi_img);
        float dy = sin_t * std::sin(phi_img);
        float dz = -cos_t;

        // Same coordinate mapping as PinholeCamera (see detailed comments there)
        ray.direction(0) = 0.0f;
        ray.direction(1) = dz;   // vr
        ray.direction(2) = -dy;  // r*vtheta (minus sign: +Y -> -theta)
        ray.direction(3) = dx;   // r*sin theta*vphi

        ray.weight = 1.0f;
        return ray;
    }

    LensType GetLensType() const override { return LensType::Fisheye; }
    const char* GetName() const override { return "Fisheye Camera"; }
    const CameraConfig& GetConfig() const override { return config_; }

    void SetConfig(const CameraConfig& config) override {
        config_ = config;
        UpdateInternals();
    }

  private:
    void UpdateInternals() { aspect_ratio_ = static_cast<float>(config_.width) / config_.height; }

    CameraConfig config_;
    float aspect_ratio_ = 1.0f;
};

// Construct a camera for the requested lens. LensType contains only represented
// models, so an ordinary typed call cannot request a silent approximation.
inline std::unique_ptr<ICamera> CreateCamera(LensType type,
                                             const CameraConfig& config = CameraConfig()) {
    switch (type) {
        case LensType::Pinhole:
            return std::make_unique<PinholeCamera>(config);
        case LensType::ThinLens:
            return std::make_unique<ThinLensCamera>(config);
        case LensType::Fisheye:
            return std::make_unique<FisheyeCamera>(config);
    }
    SIRIUS_ASSERT(false);  // Malformed enum value from an unsafe cast.
    return nullptr;
}

}  // namespace sirius::core
