#pragma once

// Camera interface and lens models (pinhole, thin-lens, fisheye) that generate
// rays in Boyer-Lindquist coordinates for the geodesic tracer. Strategy pattern
// over ICamera for interchangeable lens models. Ported from CMBS001A.h.

#include "sirius/base/contracts.h"
#include "sirius/core/tensor.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>
#include <optional>
#include <string_view>

namespace sirius::core {

// Lens projection models and their operator-visible defaults. Configuration,
// typed-session validation, and camera construction share this authority so a
// focal parameter cannot be accepted by a projection that never consumes it.
enum class LensType {
    Pinhole,   // Ideal pinhole (infinite depth of field)
    ThinLens,  // Thin lens model with depth of field
    Fisheye    // Equidistant fisheye projection
};

inline constexpr float kDefaultCameraFocalLength = 50.0f;
inline constexpr float kDefaultCameraAperture = 2.8f;
inline constexpr float kDefaultCameraFocusDistance = 50.0f;
inline constexpr float kReferenceFocalLengthMillimetres = 50.0f;
inline constexpr double kMaximumThinLensPupilFraction = 0.1;

struct ThinLensProjectionSample {
    float pupil_right = 0.0f;
    float pupil_up = 0.0f;
    float direction_forward = 1.0f;
    float direction_up = 0.0f;
    float direction_right = 0.0f;
};

[[nodiscard]] inline float ThinLensApertureRadius(float focal_length, float aperture) noexcept {
    return (focal_length / kReferenceFocalLengthMillimetres) / (2.0f * aperture);
}

// Deterministic finite-pupil sample in the camera rest frame. image_right and
// image_up are normalised image-plane coordinates before the FOV scale. The
// returned ray starts at (0, pupil_up, pupil_right) on the aperture plane and
// points to the requested point on the focus plane.
[[nodiscard]] inline ThinLensProjectionSample ProjectThinLensSample(
    float image_right, float image_up, float tan_half_fov, float focal_length, float aperture,
    float focus_distance, float sample_u, float sample_v) noexcept {
    const float focus_right = image_right * tan_half_fov * focus_distance;
    const float focus_up = image_up * tan_half_fov * focus_distance;
    const float aperture_radius = ThinLensApertureRadius(focal_length, aperture);
    const float pupil_angle = 2.0f * static_cast<float>(std::numbers::pi) * sample_u;
    const float pupil_radius = aperture_radius * std::sqrt(sample_v);

    ThinLensProjectionSample sample;
    sample.pupil_right = pupil_radius * std::cos(pupil_angle);
    sample.pupil_up = pupil_radius * std::sin(pupil_angle);
    const float forward = focus_distance;
    const float up = focus_up - sample.pupil_up;
    const float right = focus_right - sample.pupil_right;
    const float inverse_length = 1.0f / std::sqrt(forward * forward + up * up + right * right);
    sample.direction_forward = forward * inverse_length;
    sample.direction_up = up * inverse_length;
    sample.direction_right = right * inverse_length;
    return sample;
}

[[nodiscard]] constexpr std::optional<LensType> ParseLensType(std::string_view name) noexcept {
    if (name == "Pinhole") return LensType::Pinhole;
    if (name == "ThinLens") return LensType::ThinLens;
    if (name == "Fisheye") return LensType::Fisheye;
    return std::nullopt;
}

[[nodiscard]] constexpr std::optional<std::string_view> LensSpecificParameterIssue(
    LensType lens, float focal_length, float aperture, float focus_distance) noexcept {
    switch (lens) {
        case LensType::Pinhole:
        case LensType::ThinLens:
        case LensType::Fisheye:
            break;
        default:
            return "unknown lens identity";
    }
    if (lens != LensType::ThinLens &&
        (focal_length != kDefaultCameraFocalLength || aperture != kDefaultCameraAperture ||
         focus_distance != kDefaultCameraFocusDistance)) {
        return "focal length, aperture, and focus distance apply only to ThinLens";
    }
    return std::nullopt;
}

[[nodiscard]] inline std::optional<std::string_view> ThinLensGeometryIssue(
    LensType lens, double observer_radius, float focal_length, float aperture,
    float focus_distance) noexcept {
    if (lens != LensType::ThinLens) return std::nullopt;
    const double pupil_radius = ThinLensApertureRadius(focal_length, aperture);
    const double local_scale = std::min(observer_radius, static_cast<double>(focus_distance));
    if (!std::isfinite(pupil_radius) || !(pupil_radius > 0.0) || !std::isfinite(local_scale) ||
        !(local_scale > 0.0) || pupil_radius > kMaximumThinLensPupilFraction * local_scale) {
        return "thin-lens pupil exceeds the represented local tangent-plane domain";
    }
    return std::nullopt;
}

// Ray emitted by a camera: observer position and unit 4-direction.
struct CameraRay {
    Vec4 origin;                // Ray origin (observer position)
    Vec4 direction;             // Unit spatial direction in the camera rest frame
    double beta_forward = 0.0;  // Observer v/c in the screen-forward direction
    double beta_up = 0.0;       // Observer v/c in the screen-up direction
    double beta_right = 0.0;    // Observer v/c in the screen-right direction
    // Finite-aperture displacement in the camera's instantaneous rest screen.
    // The metric-aware CPU/device launch boundary applies it to the ray event;
    // changing direction alone would still be a pinhole.
    double aperture_up = 0.0;
    double aperture_right = 0.0;
    // False identifies a projection-masked sample (for example, outside the
    // circular fisheye image). Such a sample contributes black and must not be
    // passed to the geodesic tracer with its deliberately zero direction.
    bool active = true;
};

[[nodiscard]] inline bool IsRepresentedCameraRay(const CameraRay& ray) noexcept {
    for (int component = 0; component < 4; ++component) {
        if (!std::isfinite(ray.origin(component)) || !std::isfinite(ray.direction(component))) {
            return false;
        }
    }
    if (!(ray.origin(1) > 0.0) || !(ray.origin(2) > 0.0) || !(ray.origin(2) < std::numbers::pi) ||
        ray.direction(0) != 0.0) {
        return false;
    }
    const double beta_squared = ray.beta_forward * ray.beta_forward + ray.beta_up * ray.beta_up +
                                ray.beta_right * ray.beta_right;
    if (!std::isfinite(beta_squared) || beta_squared >= 1.0 || !std::isfinite(ray.aperture_up) ||
        !std::isfinite(ray.aperture_right)) {
        return false;
    }
    const double direction_norm_squared = ray.direction(1) * ray.direction(1) +
                                          ray.direction(2) * ray.direction(2) +
                                          ray.direction(3) * ray.direction(3);
    if (!std::isfinite(direction_norm_squared)) return false;
    return ray.active ? std::abs(direction_norm_squared - 1.0) <= 2.0e-5
                      : direction_norm_squared == 0.0;
}

// Observer placement, orientation, and lens/image properties.
struct CameraConfig {
    // Position (Boyer-Lindquist coordinates: t, r, theta, phi)
    double t = 0.0;
    double r = 50.0;                        // Geometric coordinate radius from the centre.
    double theta = std::numbers::pi / 2.0;  // Polar angle (pi/2 = equatorial)
    double phi = 0.0;                       // Azimuthal angle

    // Orientation
    float yaw = 0.0f;    // Yaw rotation (radians)
    float pitch = 0.0f;  // Pitch rotation (radians)
    float roll = 0.0f;   // Roll rotation (radians)

    // Camera four-velocity: spatial beta = v/c in the screen's forward, up,
    // right orthonormal axes. The tracer combines this worldline with the rest
    // screen ray in a metric-orthonormal tetrad; it is not an Euclidean edit of
    // the ray direction. Finite |beta| < 1 is a fail-closed precondition.
    double beta_x = 0.0;
    double beta_y = 0.0;
    double beta_z = 0.0;

    // Lens properties
    float fov = 60.0f;                               // Field of view (degrees)
    float focal_length = kDefaultCameraFocalLength;  // mm-equivalent, for ThinLens.
    float aperture = kDefaultCameraAperture;         // f-number, for ThinLens.
    float focus_distance =
        kDefaultCameraFocusDistance;  // Geometric coordinate length, for ThinLens.

    // Image properties
    int width = 1920;
    int height = 1080;
};

// Intrinsic camera-domain authority. Metric-specific observer placement remains
// owned by the session boundary; these are the minimum conditions under which
// a concrete lens represents the supplied typed request without undefined or
// ignored arithmetic.
[[nodiscard]] inline std::optional<std::string_view> CameraConfigIssue(
    LensType lens, const CameraConfig& config) noexcept {
    switch (lens) {
        case LensType::Pinhole:
        case LensType::ThinLens:
        case LensType::Fisheye:
            break;
        default:
            return "unknown lens identity";
    }
    if (!std::isfinite(config.t) || !std::isfinite(config.r) || config.r <= 0.0 ||
        !std::isfinite(config.theta) || config.theta <= 0.0 || config.theta >= std::numbers::pi ||
        !std::isfinite(config.phi) || !std::isfinite(config.yaw) || !std::isfinite(config.pitch) ||
        !std::isfinite(config.roll)) {
        return "camera placement or orientation is outside the intrinsic domain";
    }
    const double beta_squared = config.beta_x * config.beta_x + config.beta_y * config.beta_y +
                                config.beta_z * config.beta_z;
    if (!std::isfinite(beta_squared) || beta_squared >= 1.0) {
        return "camera beta magnitude must be finite and below one";
    }
    const float maximum_fov = lens == LensType::Fisheye ? 360.0f : 179.0f;
    if (!std::isfinite(config.fov) || config.fov <= 0.0f || config.fov > maximum_fov ||
        config.width <= 0 || config.height <= 0) {
        return "field of view and image dimensions are outside the lens domain";
    }
    if (!std::isfinite(config.focal_length) || config.focal_length <= 0.0f ||
        !std::isfinite(config.aperture) || config.aperture <= 0.0f ||
        !std::isfinite(config.focus_distance) || config.focus_distance <= 0.0f) {
        return "thin-lens parameters must be finite and positive";
    }
    if (const auto issue = LensSpecificParameterIssue(lens, config.focal_length, config.aperture,
                                                      config.focus_distance);
        issue.has_value()) {
        return issue;
    }
    return ThinLensGeometryIssue(lens, config.r, config.focal_length, config.aperture,
                                 config.focus_distance);
}

// Abstract camera: generates a ray per pixel sample.
class ICamera {
  public:
    virtual ~ICamera() = default;

    // Generate a ray for pixel (x, y) with sample offset (u, v) in [0, 1).
    virtual CameraRay GenerateRay(int x, int y, float u = 0.5f, float v = 0.5f) const = 0;

    // Generate the rest-frame screen ray and bind its observer worldline. The
    // metric-aware tracer performs the tetrad boost at the launch event.
    CameraRay GenerateRayForObserver(int x, int y, float u = 0.5f, float v = 0.5f) const {
        CameraRay ray = GenerateRay(x, y, u, v);
        const auto& cfg = GetConfig();
        const double beta_squared =
            cfg.beta_x * cfg.beta_x + cfg.beta_y * cfg.beta_y + cfg.beta_z * cfg.beta_z;
        SIRIUS_PRE(std::isfinite(beta_squared) && beta_squared < 1.0);
        ray.beta_forward = cfg.beta_x;
        ray.beta_up = cfg.beta_y;
        ray.beta_right = cfg.beta_z;
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

  protected:
    void RequirePixelSample(int x, int y, float u, float v) const {
        const auto& config = GetConfig();
        SIRIUS_PRE(x >= 0 && x < config.width && y >= 0 && y < config.height);
        SIRIUS_PRE(std::isfinite(u) && u >= 0.0f && u < 1.0f);
        SIRIUS_PRE(std::isfinite(v) && v >= 0.0f && v < 1.0f);
    }
};

// Ideal perspective projection (infinite depth of field).
class PinholeCamera : public ICamera {
  public:
    explicit PinholeCamera(const CameraConfig& config = CameraConfig()) : config_(config) {
        SIRIUS_PRE(!CameraConfigIssue(LensType::Pinhole, config).has_value());
        UpdateInternals();
    }

    CameraRay GenerateRay(int x, int y, float u = 0.5f, float v = 0.5f) const override {
        RequirePixelSample(x, y, u, v);
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

        return ray;
    }

    LensType GetLensType() const override { return LensType::Pinhole; }
    const char* GetName() const override { return "Pinhole Camera"; }
    const CameraConfig& GetConfig() const override { return config_; }

    void SetConfig(const CameraConfig& config) override {
        SIRIUS_PRE(!CameraConfigIssue(LensType::Pinhole, config).has_value());
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
        SIRIUS_PRE(!CameraConfigIssue(LensType::ThinLens, config).has_value());
        UpdateInternals();
    }

    CameraRay GenerateRay(int x, int y, float u = 0.5f, float v = 0.5f) const override {
        RequirePixelSample(x, y, u, v);
        CameraRay ray;

        ray.origin(0) = config_.t;
        ray.origin(1) = config_.r;
        ray.origin(2) = config_.theta;
        ray.origin(3) = config_.phi;

        // Point on image plane
        float px = (2.0f * (x + u) / config_.width - 1.0f) * aspect_ratio_;
        float py = 1.0f - 2.0f * (y + v) / config_.height;

        // The spacetime uses geometric coordinate units and has no physical
        // mass-to-millimetre scale, so 50 mm-equivalent defines one virtual
        // lens unit. The launch boundary applies the returned finite pupil
        // displacement in its metric-orthonormal camera frame.
        const ThinLensProjectionSample sample =
            ProjectThinLensSample(px, py, tan_half_fov_, config_.focal_length, config_.aperture,
                                  config_.focus_distance, u, v);
        ray.aperture_up = sample.pupil_up;
        ray.aperture_right = sample.pupil_right;

        // Same coordinate mapping as PinholeCamera (see detailed comments there)
        ray.direction(0) = 0.0f;
        ray.direction(1) = -sample.direction_forward;  // screen forward = radial inward
        ray.direction(2) = -sample.direction_up;       // screen up = decreasing theta
        ray.direction(3) = sample.direction_right;     // screen right = increasing phi

        return ray;
    }

    LensType GetLensType() const override { return LensType::ThinLens; }
    const char* GetName() const override { return "Thin Lens Camera"; }
    const CameraConfig& GetConfig() const override { return config_; }

    void SetConfig(const CameraConfig& config) override {
        SIRIUS_PRE(!CameraConfigIssue(LensType::ThinLens, config).has_value());
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
        SIRIUS_PRE(!CameraConfigIssue(LensType::Fisheye, config).has_value());
        UpdateInternals();
    }

    CameraRay GenerateRay(int x, int y, float u = 0.5f, float v = 0.5f) const override {
        RequirePixelSample(x, y, u, v);
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
            ray.active = false;
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

        return ray;
    }

    LensType GetLensType() const override { return LensType::Fisheye; }
    const char* GetName() const override { return "Fisheye Camera"; }
    const CameraConfig& GetConfig() const override { return config_; }

    void SetConfig(const CameraConfig& config) override {
        SIRIUS_PRE(!CameraConfigIssue(LensType::Fisheye, config).has_value());
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
