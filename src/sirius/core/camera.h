#pragma once

// Camera interface and lens models (pinhole, thin-lens, fisheye) that generate
// rays in Boyer-Lindquist coordinates for the geodesic tracer. Strategy pattern
// over ICamera for interchangeable lens models. Ported from CMBS001A.h.

#include "sirius/base/contracts.h"
#include "sirius/core/celestial_tangent_basis.h"
#include "sirius/core/tensor.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <expected>
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

// Local differential of normalized launch direction with respect to continuous
// film coordinates measured in pixels (+x right, +y down), at fixed pupil.
// Rows of angular_jacobian use the same rest-frame celestial tangent basis as
// ObserverScreenBasis. The absolute determinant is a solid-angle density
// (steradians per square film pixel), not a pixel-integrated solid angle.
// Observer beta is already represented by the tracer's boosted screen: applying
// aberration here would count that transformation twice.
struct CameraFilmDifferential {
    std::array<std::array<double, 2>, 3> direction_derivative;
    std::array<std::array<double, 2>, 2> angular_jacobian;
    double signed_solid_angle_density;
    double solid_angle_density;
};

// A smooth geometric extension of the lens, in double film-pixel coordinates.
// It does not differentiate or reproduce the legacy float arithmetic staircase.
// A successful inactive ray denotes optical masking; it is never a trace input.
// The antipodal fisheye rim has a ray but no regular angular differential.
enum class CameraProjectionFailure { InvalidInput, Unsupported, Arithmetic, Unrepresentable };
struct CameraFilmProjection {
    CameraRay ray;
    std::optional<CameraFilmDifferential> differential;
};

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
    // Omitting the independent pupil pair selects the aperture centre.
    virtual CameraRay GenerateRay(int x, int y, float image_u = 0.5f, float image_v = 0.5f,
                                  float pupil_u = 0.5f, float pupil_v = 0.0f) const = 0;

    // Generate the rest-frame screen ray and bind its observer worldline. The
    // metric-aware tracer performs the tetrad boost at the launch event.
    CameraRay GenerateRayForObserver(int x, int y, float image_u = 0.5f, float image_v = 0.5f,
                                     float pupil_u = 0.5f, float pupil_v = 0.0f) const {
        CameraRay ray = GenerateRay(x, y, image_u, image_v, pupil_u, pupil_v);
        const auto& cfg = GetConfig();
        const double beta_squared =
            cfg.beta_x * cfg.beta_x + cfg.beta_y * cfg.beta_y + cfg.beta_z * cfg.beta_z;
        SIRIUS_PRE(std::isfinite(beta_squared) && beta_squared < 1.0);
        ray.beta_forward = cfg.beta_x;
        ray.beta_up = cfg.beta_y;
        ray.beta_right = cfg.beta_z;
        return ray;
    }

    // Optional for custom projections. The built-in lenses differentiate their
    // smooth geometric projection using the stored lens coefficients; float
    // quantization itself is not differentiated. GenerateRay remains unchanged.
    // Masked samples and the collapsed fisheye antipodal rim are unavailable.
    [[nodiscard]] virtual std::optional<CameraFilmDifferential> FilmDifferentialForObserver(
        int, int, float = 0.5f, float = 0.5f, float = 0.5f, float = 0.0f) const {
        return std::nullopt;
    }

    // Output-crop-independent film coordinates (+x right, +y down). Built-in
    // lenses preserve stored lens coefficients and the original float pupil
    // sample. This leaves GenerateRay's nominal float results unchanged.
    [[nodiscard]] virtual std::expected<CameraFilmProjection, CameraProjectionFailure>
    ProjectFilmForObserver(double, double, float = 0.5f, float = 0.0f) const {
        return std::unexpected(CameraProjectionFailure::Unsupported);
    }

    // Refinement must not silently reuse a rounded-away coordinate or active
    // direction. Zero offset remains a valid request for the original geometry.
    [[nodiscard]] std::expected<CameraFilmProjection, CameraProjectionFailure>
    ProjectFilmOffsetForObserver(double x, double y, double dx, double dy, float pupil_u = 0.5f,
                                 float pupil_v = 0.0f) const {
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(dx) || !std::isfinite(dy))
            return std::unexpected(CameraProjectionFailure::InvalidInput);
        const double shifted_x = x + dx, shifted_y = y + dy;
        if (!std::isfinite(shifted_x) || !std::isfinite(shifted_y))
            return std::unexpected(CameraProjectionFailure::Arithmetic);
        if ((dx != 0.0 && shifted_x == x) || (dy != 0.0 && shifted_y == y))
            return std::unexpected(CameraProjectionFailure::Unrepresentable);
        auto result = ProjectFilmForObserver(shifted_x, shifted_y, pupil_u, pupil_v);
        if (!result || !result->ray.active || (dx == 0.0 && dy == 0.0)) return result;
        const auto original = ProjectFilmForObserver(x, y, pupil_u, pupil_v);
        if (!original) return std::unexpected(original.error());
        if (original->ray.active && original->ray.direction(1) == result->ray.direction(1) &&
            original->ray.direction(2) == result->ray.direction(2) &&
            original->ray.direction(3) == result->ray.direction(3))
            return std::unexpected(CameraProjectionFailure::Unrepresentable);
        return result;
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
    void RequirePixelSample(int x, int y, float image_u, float image_v, float pupil_u,
                            float pupil_v) const {
        const auto& config = GetConfig();
        SIRIUS_PRE(x >= 0 && x < config.width && y >= 0 && y < config.height);
        SIRIUS_PRE(std::isfinite(image_u) && image_u >= 0.0f && image_u < 1.0f);
        SIRIUS_PRE(std::isfinite(image_v) && image_v >= 0.0f && image_v < 1.0f);
        SIRIUS_PRE(std::isfinite(pupil_u) && pupil_u >= 0.0f && pupil_u < 1.0f);
        SIRIUS_PRE(std::isfinite(pupil_v) && pupil_v >= 0.0f && pupil_v < 1.0f);
    }
};

// Ideal perspective projection (infinite depth of field).
class PinholeCamera : public ICamera {
  public:
    explicit PinholeCamera(const CameraConfig& config = CameraConfig()) : config_(config) {
        SIRIUS_PRE(!CameraConfigIssue(LensType::Pinhole, config).has_value());
        UpdateInternals();
    }

    CameraRay GenerateRay(int x, int y, float image_u = 0.5f, float image_v = 0.5f,
                          [[maybe_unused]] float pupil_u = 0.5f,
                          [[maybe_unused]] float pupil_v = 0.0f) const override {
        RequirePixelSample(x, y, image_u, image_v, pupil_u, pupil_v);
        CameraRay ray;

        // Set origin
        ray.origin(0) = config_.t;
        ray.origin(1) = config_.r;
        ray.origin(2) = config_.theta;
        ray.origin(3) = config_.phi;

        // Normalised device coordinates (-1 to 1)
        float px = (2.0f * (x + image_u) / config_.width - 1.0f) * aspect_ratio_;
        float py = 1.0f - 2.0f * (y + image_v) / config_.height;

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

    [[nodiscard]] std::optional<CameraFilmDifferential> FilmDifferentialForObserver(
        int x, int y, float image_u = 0.5f, float image_v = 0.5f, float pupil_u = 0.5f,
        float pupil_v = 0.0f) const override;

    [[nodiscard]] std::expected<CameraFilmProjection, CameraProjectionFailure>
    ProjectFilmForObserver(double film_x, double film_y, float pupil_u = 0.5f,
                           float pupil_v = 0.0f) const override;

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

    CameraRay GenerateRay(int x, int y, float image_u = 0.5f, float image_v = 0.5f,
                          float pupil_u = 0.5f, float pupil_v = 0.0f) const override {
        RequirePixelSample(x, y, image_u, image_v, pupil_u, pupil_v);
        CameraRay ray;

        ray.origin(0) = config_.t;
        ray.origin(1) = config_.r;
        ray.origin(2) = config_.theta;
        ray.origin(3) = config_.phi;

        // Point on image plane
        float px = (2.0f * (x + image_u) / config_.width - 1.0f) * aspect_ratio_;
        float py = 1.0f - 2.0f * (y + image_v) / config_.height;

        // The spacetime uses geometric coordinate units and has no physical
        // mass-to-millimetre scale, so 50 mm-equivalent defines one virtual
        // lens unit. The launch boundary applies the returned finite pupil
        // displacement in its metric-orthonormal camera frame.
        const ThinLensProjectionSample sample =
            ProjectThinLensSample(px, py, tan_half_fov_, config_.focal_length, config_.aperture,
                                  config_.focus_distance, pupil_u, pupil_v);
        ray.aperture_up = sample.pupil_up;
        ray.aperture_right = sample.pupil_right;

        // Same coordinate mapping as PinholeCamera (see detailed comments there)
        ray.direction(0) = 0.0f;
        ray.direction(1) = -sample.direction_forward;  // screen forward = radial inward
        ray.direction(2) = -sample.direction_up;       // screen up = decreasing theta
        ray.direction(3) = sample.direction_right;     // screen right = increasing phi

        return ray;
    }

    [[nodiscard]] std::optional<CameraFilmDifferential> FilmDifferentialForObserver(
        int x, int y, float image_u = 0.5f, float image_v = 0.5f, float pupil_u = 0.5f,
        float pupil_v = 0.0f) const override;

    [[nodiscard]] std::expected<CameraFilmProjection, CameraProjectionFailure>
    ProjectFilmForObserver(double film_x, double film_y, float pupil_u = 0.5f,
                           float pupil_v = 0.0f) const override;

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

    CameraRay GenerateRay(int x, int y, float image_u = 0.5f, float image_v = 0.5f,
                          [[maybe_unused]] float pupil_u = 0.5f,
                          [[maybe_unused]] float pupil_v = 0.0f) const override {
        RequirePixelSample(x, y, image_u, image_v, pupil_u, pupil_v);
        CameraRay ray;

        ray.origin(0) = config_.t;
        ray.origin(1) = config_.r;
        ray.origin(2) = config_.theta;
        ray.origin(3) = config_.phi;

        // Normalised coordinates from centre
        float px = (2.0f * (x + image_u) / config_.width - 1.0f) * aspect_ratio_;
        float py = 1.0f - 2.0f * (y + image_v) / config_.height;

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

    [[nodiscard]] std::optional<CameraFilmDifferential> FilmDifferentialForObserver(
        int x, int y, float image_u = 0.5f, float image_v = 0.5f, float pupil_u = 0.5f,
        float pupil_v = 0.0f) const override;

    [[nodiscard]] std::expected<CameraFilmProjection, CameraProjectionFailure>
    ProjectFilmForObserver(double film_x, double film_y, float pupil_u = 0.5f,
                           float pupil_v = 0.0f) const override;

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

namespace camera_detail {

// Continuous geometry uses the same stored float coefficients and original
// pupil arithmetic as the nominal projection, but double film evaluation. All
// derivatives below belong to this same smooth function, including outside the
// output crop; no isolated nominal-point rounding branch is introduced.
[[nodiscard]] inline std::expected<CameraFilmProjection, CameraProjectionFailure>
ProjectContinuousFilm(const CameraConfig& config, LensType lens, double film_x, double film_y,
                      float pupil_u, float pupil_v) {
    using Failure = CameraProjectionFailure;
    using Vector = std::array<double, 3>;
    if (CameraConfigIssue(lens, config) || !std::isfinite(film_x) || !std::isfinite(film_y) ||
        !std::isfinite(pupil_u) || pupil_u < 0 || pupil_u >= 1 || !std::isfinite(pupil_v) ||
        pupil_v < 0 || pupil_v >= 1)
        return std::unexpected(Failure::InvalidInput);
    CameraFilmProjection result{};
    auto& ray = result.ray;
    ray.origin(0) = config.t;
    ray.origin(1) = config.r;
    ray.origin(2) = config.theta;
    ray.origin(3) = config.phi;
    ray.beta_forward = config.beta_x;
    ray.beta_up = config.beta_y;
    ray.beta_right = config.beta_z;
    const double aspect = static_cast<float>(config.width) / config.height;
    const double px = (2.0 * film_x / config.width - 1.0) * aspect;
    const double py = 1.0 - 2.0 * film_y / config.height;
    const double sx = 2.0 * aspect / config.width, sy = -2.0 / config.height;
    if (!std::isfinite(px) || !std::isfinite(py)) return std::unexpected(Failure::Arithmetic);
    Vector q{};
    std::array<Vector, 2> dq{};
    bool regular = true;
    if (lens == LensType::Pinhole || lens == LensType::ThinLens) {
        const double t = std::tan(config.fov * static_cast<float>(std::numbers::pi) / 360.0f);
        if (lens == LensType::ThinLens) {
            // Calling at the image centre retains the existing pupil arithmetic
            // without narrowing the continuous film coordinates.
            const auto pupil =
                ProjectThinLensSample(0, 0, static_cast<float>(t), config.focal_length,
                                      config.aperture, config.focus_distance, pupil_u, pupil_v);
            ray.aperture_up = pupil.pupil_up;
            ray.aperture_right = pupil.pupil_right;
            q = {-config.focus_distance, ray.aperture_up - py * t * config.focus_distance,
                 px * t * config.focus_distance - ray.aperture_right};
            dq = {Vector{0, 0, sx * t * config.focus_distance},
                  Vector{0, -sy * t * config.focus_distance, 0}};
        } else {
            const double cy = std::cos(config.yaw), syaw = std::sin(config.yaw);
            const double cp = std::cos(config.pitch), sp = std::sin(config.pitch);
            const double cr = std::cos(config.roll), sr = std::sin(config.roll);
            const auto rotate = [&](Vector v) {
                const double xr = v[0] * cr - v[1] * sr, yr = v[0] * sr + v[1] * cr;
                const double yp = yr * cp - v[2] * sp, zp = yr * sp + v[2] * cp;
                return Vector{-xr * syaw + zp * cy, -yp, xr * cy + zp * syaw};
            };
            q = rotate({px * t, py * t, -1});
            dq = {rotate({sx * t, 0, 0}), rotate({0, sy * t, 0})};
        }
    } else {
        const double a =
            config.fov * static_cast<double>(static_cast<float>(std::numbers::pi)) / 360.0;
        const double radius = std::hypot(px, py), theta = a * radius;
        if (!std::isfinite(theta)) return std::unexpected(Failure::Arithmetic);
        // The nominal camera uses the represented float pi as its mask edge.
        if (theta > static_cast<float>(std::numbers::pi)) {
            ray.active = false;
            return result;
        }
        regular = theta < std::numbers::pi;
        double sinc, radial;
        if (std::abs(theta) < 1e-3) {
            const double t2 = theta * theta;
            sinc = a * (1 - t2 / 6 + t2 * t2 / 120 - t2 * t2 * t2 / 5040);
            radial = a * a * a * (-1.0 / 3 + t2 / 30 - t2 * t2 / 840);
        } else {
            sinc = std::sin(theta) / radius;
            radial = (theta * std::cos(theta) - std::sin(theta)) / (radius * radius * radius);
        }
        q = {-std::cos(theta), -sinc * py, sinc * px};
        dq = {Vector{a * sinc * px * sx, -py * radial * px * sx, (sinc + px * px * radial) * sx},
              Vector{a * sinc * py * sy, -(sinc + py * py * radial) * sy, px * radial * py * sy}};
    }
    for (double value : q)
        if (!std::isfinite(value)) return std::unexpected(Failure::Arithmetic);
    const double length = std::hypot(q[0], q[1], q[2]);
    if (!std::isfinite(length) || !(length > 0)) return std::unexpected(Failure::Arithmetic);
    Vector n{};
    for (unsigned i = 0; i < 3; ++i) ray.direction(i + 1) = n[i] = q[i] / length;
    if (!IsRepresentedCameraRay(ray)) return std::unexpected(Failure::Arithmetic);
    if (!regular) return result;
    const auto basis = relativity::MakeCelestialTangentBasis(n);
    if (!basis) return std::unexpected(Failure::Arithmetic);
    CameraFilmDifferential map{};
    for (unsigned column = 0; column < 2; ++column) {
        double longitudinal = 0;
        for (unsigned i = 0; i < 3; ++i) longitudinal += n[i] * dq[column][i];
        for (unsigned i = 0; i < 3; ++i) {
            const double d = (dq[column][i] - n[i] * longitudinal) / length;
            if (!std::isfinite(d)) return std::unexpected(Failure::Arithmetic);
            map.direction_derivative[i][column] = d;
            map.angular_jacobian[0][column] += basis->first[i] * d;
            map.angular_jacobian[1][column] += basis->second[i] * d;
        }
    }
    map.signed_solid_angle_density =
        std::fma(map.angular_jacobian[0][0], map.angular_jacobian[1][1],
                 -map.angular_jacobian[0][1] * map.angular_jacobian[1][0]);
    map.solid_angle_density = std::abs(map.signed_solid_angle_density);
    if (!std::isfinite(map.solid_angle_density)) return std::unexpected(Failure::Arithmetic);
    if (!(map.solid_angle_density > 0)) return std::unexpected(Failure::Unrepresentable);
    result.differential = map;
    return result;
}

[[nodiscard]] inline std::optional<CameraFilmDifferential> MeasureFilmDifferential(
    const ICamera& camera, int x, int y, float image_u, float image_v, float pupil_u,
    float pupil_v) {
    const CameraRay ray = camera.GenerateRayForObserver(x, y, image_u, image_v, pupil_u, pupil_v);
    if (!ray.active || !IsRepresentedCameraRay(ray)) return std::nullopt;
    const auto& config = camera.GetConfig();
    using Vector = std::array<double, 3>;
    const Vector raw{ray.direction(1), ray.direction(2), ray.direction(3)};
    const double raw_norm = std::hypot(raw[0], raw[1], raw[2]);
    const Vector n{raw[0] / raw_norm, raw[1] / raw_norm, raw[2] / raw_norm};
    const auto basis = relativity::MakeCelestialTangentBasis(raw);
    if (!basis) return std::nullopt;
    const float aspect = static_cast<float>(config.width) / config.height;
    const float px = (2.0f * (x + image_u) / config.width - 1.0f) * aspect;
    const float py = 1.0f - 2.0f * (y + image_v) / config.height;
    const double film_x = 2.0 * static_cast<double>(aspect) / config.width;
    const double film_y = -2.0 / config.height;
    Vector q{};
    std::array<Vector, 2> dq{};
    if (camera.GetLensType() == LensType::Pinhole) {
        const float t = std::tan(config.fov * static_cast<float>(std::numbers::pi) / 360.0f);
        const float cy = std::cos(config.yaw), sy = std::sin(config.yaw);
        const float cp = std::cos(config.pitch), sp = std::sin(config.pitch);
        const float cr = std::cos(config.roll), sr = std::sin(config.roll);
        // Preserve the nominal projection's rounded intermediates without
        // changing GenerateRay. Tangents follow its same roll/pitch/yaw order.
        const float dx = px * t, dy = py * t, dz = -1.0f;
        const float rx = dx * cr - dy * sr, ry = dx * sr + dy * cr;
        const float ry2 = ry * cp - dz * sp, rz = ry * sp + dz * cp;
        const float rx2 = rx * cy + rz * sy, rz2 = -rx * sy + rz * cy;
        q = {rz2, -ry2, rx2};
        const auto rotate = [&](Vector v) {
            const double xr = v[0] * cr - v[1] * sr, yr = v[0] * sr + v[1] * cr;
            const double yp = yr * cp - v[2] * sp, zp = yr * sp + v[2] * cp;
            return Vector{-xr * sy + zp * cy, -yp, xr * cy + zp * sy};
        };
        dq = {rotate({film_x * t, 0, 0}), rotate({0, film_y * t, 0})};
    } else if (camera.GetLensType() == LensType::ThinLens) {
        const float t = std::tan(config.fov * static_cast<float>(std::numbers::pi) / 360.0f);
        const float focus_right = px * t * config.focus_distance;
        const float focus_up = py * t * config.focus_distance;
        const float up = focus_up - static_cast<float>(ray.aperture_up);
        const float right = focus_right - static_cast<float>(ray.aperture_right);
        q = {-config.focus_distance, -up, right};
        dq = {Vector{0, 0, film_x * t * config.focus_distance},
              Vector{0, -film_y * t * config.focus_distance, 0}};
        // ThinLens currently has no orientation rotation in GenerateRay;
        // applying one only to the differential would be inconsistent.
    } else if (camera.GetLensType() == LensType::Fisheye) {
        const float image_radius = std::sqrt(px * px + py * py);
        const float theta_ray =
            image_radius * config.fov * static_cast<float>(std::numbers::pi) / 360.0f;
        if (theta_ray >= static_cast<float>(std::numbers::pi)) return std::nullopt;
        const double a =
            static_cast<double>(config.fov) * static_cast<float>(std::numbers::pi) / 360.0;
        const double radius = std::hypot(static_cast<double>(px), static_cast<double>(py));
        const double theta = a * radius;
        double sinc, radial;
        if (std::abs(theta) < 1e-3) {
            const double t2 = theta * theta;
            sinc = a * (1.0 - t2 / 6.0 + t2 * t2 / 120.0 - t2 * t2 * t2 / 5040.0);
            radial = a * a * a * (-1.0 / 3.0 + t2 / 30.0 - t2 * t2 / 840.0);
        } else {
            sinc = std::sin(theta) / radius;
            radial = (theta * std::cos(theta) - std::sin(theta)) / (radius * radius * radius);
        }
        q = raw;
        dq = {Vector{a * sinc * px * film_x, -py * radial * px * film_x,
                     (sinc + px * px * radial) * film_x},
              Vector{a * sinc * py * film_y, -(sinc + py * py * radial) * film_y,
                     px * radial * py * film_y}};
        // Fisheye likewise currently ignores yaw/pitch/roll. The analytic
        // centre limit avoids differentiating atan2(0,0).
    } else {
        return std::nullopt;
    }
    const double length = std::hypot(q[0], q[1], q[2]);
    if (!std::isfinite(length) || !(length > 0.0)) return std::nullopt;
    CameraFilmDifferential result{};
    for (unsigned column = 0; column < 2; ++column) {
        double longitudinal = 0.0;
        for (unsigned i = 0; i < 3; ++i) longitudinal += n[i] * dq[column][i];
        for (unsigned i = 0; i < 3; ++i) {
            const double derivative = (dq[column][i] - n[i] * longitudinal) / length;
            if (!std::isfinite(derivative)) return std::nullopt;
            result.direction_derivative[i][column] = derivative;
            result.angular_jacobian[0][column] += basis->first[i] * derivative;
            result.angular_jacobian[1][column] += basis->second[i] * derivative;
        }
    }
    result.signed_solid_angle_density =
        result.angular_jacobian[0][0] * result.angular_jacobian[1][1] -
        result.angular_jacobian[0][1] * result.angular_jacobian[1][0];
    result.solid_angle_density = std::abs(result.signed_solid_angle_density);
    if (!std::isfinite(result.solid_angle_density)) return std::nullopt;
    return result;
}
}  // namespace camera_detail

inline std::optional<CameraFilmDifferential> PinholeCamera::FilmDifferentialForObserver(
    int x, int y, float image_u, float image_v, float pupil_u, float pupil_v) const {
    return camera_detail::MeasureFilmDifferential(*this, x, y, image_u, image_v, pupil_u, pupil_v);
}
inline std::optional<CameraFilmDifferential> ThinLensCamera::FilmDifferentialForObserver(
    int x, int y, float image_u, float image_v, float pupil_u, float pupil_v) const {
    return camera_detail::MeasureFilmDifferential(*this, x, y, image_u, image_v, pupil_u, pupil_v);
}
inline std::optional<CameraFilmDifferential> FisheyeCamera::FilmDifferentialForObserver(
    int x, int y, float image_u, float image_v, float pupil_u, float pupil_v) const {
    return camera_detail::MeasureFilmDifferential(*this, x, y, image_u, image_v, pupil_u, pupil_v);
}

inline std::expected<CameraFilmProjection, CameraProjectionFailure>
PinholeCamera::ProjectFilmForObserver(double x, double y, float pupil_u, float pupil_v) const {
    return camera_detail::ProjectContinuousFilm(config_, GetLensType(), x, y, pupil_u, pupil_v);
}
inline std::expected<CameraFilmProjection, CameraProjectionFailure>
ThinLensCamera::ProjectFilmForObserver(double x, double y, float pupil_u, float pupil_v) const {
    return camera_detail::ProjectContinuousFilm(config_, GetLensType(), x, y, pupil_u, pupil_v);
}
inline std::expected<CameraFilmProjection, CameraProjectionFailure>
FisheyeCamera::ProjectFilmForObserver(double x, double y, float pupil_u, float pupil_v) const {
    return camera_detail::ProjectContinuousFilm(config_, GetLensType(), x, y, pupil_u, pupil_v);
}

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
