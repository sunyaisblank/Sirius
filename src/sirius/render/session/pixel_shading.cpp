// Physical pixel radiance for both serial and parallel session workers.
// This file owns source/disk sampling; session lifecycle and publication stay
// in render_session.cpp.

#include "sirius/core/constants.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/spectral/blackbody.h"
#include "sirius/render/pixel_sampling.h"
#include "sirius/render/session/point_source_detector.h"
#include "sirius/render/session/render_session.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <optional>
#include <string>

namespace sirius::render {

using backend::GeodesicTracer;
using backend::TraceResult;
using core::CameraRay;
using core::Vec4;
namespace math = core::constants::math;

// =============================================================================
// Pixel shading helpers (unified for single- and multi-threaded paths).
// =============================================================================
RenderSession::PixelResult RenderSession::ShadeDiskHit(const TraceResult& result) const {
    PixelResult px;

    // Volumetric disk: use the pre-integrated emission from ray marching.
    if (result.volumetric_hit) {
        // Samples were coloured at their own temperature and g-factor before
        // invariant transfer. Reconstructing one effective blackbody here would
        // apply the spectral mapping twice.
        px.r = result.volumetric_emission[0];
        px.g = result.volumetric_emission[1];
        px.b = result.volumetric_emission[2];
        return px;
    }

    // Thin disk: accumulate emission from all disk crossings. Relativistic
    // beaming is applied exactly once: emitted T^4 becomes observed g^4 T^4.
    float total_r = 0.0f, total_g = 0.0f, total_b = 0.0f;
    core::StokesVector total_stokes;
    const bool polarisation_mode = config_.color_mode == core::color_modes::Mode::Polarisation;

    for (int crossing_idx = 0; crossing_idx < result.num_disk_crossings; crossing_idx++) {
        const auto& crossing = result.disk_crossings[crossing_idx];
        if (!crossing.valid) continue;

        float T_emit = crossing.temperature;
        float g = crossing.redshift;

        // The stationary axisymmetric disk has one covariant transfer sample at
        // a crossing. Temporal requests fail at validation until an evolving
        // emissivity (rather than an azimuth-shifted steady flow) is represented.
        const std::array<float, 1> temporal_redshifts{g};

        const float emitted_intensity = std::pow(T_emit, 4.0f);
        const float observed_intensity =
            core::color_modes::ObservedBolometricIntensity(emitted_intensity, g);

        if (polarisation_mode) {
            SIRIUS_ASSERT(crossing.polarisation_valid);
            if (!crossing.polarisation_valid) continue;

            const float chi = crossing.polarisation_evpa;
            const float degree = crossing.polarisation_degree;
            const float atmosphere_intensity =
                observed_intensity * crossing.polarisation_intensity_scale;
            core::StokesVector crossing_stokes{
                atmosphere_intensity, atmosphere_intensity * degree * std::cos(2.0f * chi),
                atmosphere_intensity * degree * std::sin(2.0f * chi), 0.0f};
            total_stokes += crossing_stokes;
            continue;
        }

        const core::spectral::Rgb disk_color = core::color_modes::AverageTemporalColorMode(
            config_.color_mode, T_emit, temporal_redshifts, emitted_intensity,
            config_.disk_temperature_scale);

        total_r += disk_color.r;
        total_g += disk_color.g;
        total_b += disk_color.b;
    }

    // Lensing changes the ray-to-solid-angle map, not radiance along one ray.
    constexpr float output_scale = 1.0f;

    if (polarisation_mode) {
        total_stokes *= output_scale;
        total_stokes.Normalise();
        const core::spectral::Rgb visualised =
            core::color_modes::polarisation_vis::StokesToRgbHsv(total_stokes);
        px.r = visualised.r;
        px.g = visualised.g;
        px.b = visualised.b;
    } else {
        px.r = total_r * output_scale;
        px.g = total_g * output_scale;
        px.b = total_b * output_scale;
    }

    return px;
}

RenderSession::PixelResult RenderSession::ShadeEscaped(const TraceResult& result) const {
    PixelResult px;
    if (config_.point_starfield && star_index_ && star_index_->Size() > 0) {
        SampleStarfieldPoints(result, px.r, px.g, px.b);
    } else {
        SampleStarfield(result.final_direction, px.r, px.g, px.b);
    }

    return px;
}

namespace {
std::expected<PointDetectorProbe, PointDetectorFailure> MeasurePointDetector(
    const core::CameraFilmProjection& film, const TraceResult& traced,
    const core::AngularMatrix2& film_from_standard) {
    if (traced.cancelled) return std::unexpected(PointDetectorFailure::Cancelled);
    if (traced.numerical_failure || traced.outcome == TraceResult::Outcome::MaxSteps)
        return std::unexpected(PointDetectorFailure::TraceFailed);
    PointDetectorProbe value;
    value.inner_attempts = static_cast<std::size_t>(traced.steps_taken);
    if (traced.outcome != TraceResult::Outcome::Escaped) return value;
    if (!film.differential || !traced.beam.infinity_source_map)
        return std::unexpected(PointDetectorFailure::ProjectionUnavailable);
    const auto& sky = *traced.beam.infinity_source_map;
    value.tail_attempts = sky.attempted_steps;
    value.visible = true;
    value.direction = sky.map.direction;
    value.camera_over_source_frequency = 1.0 / sky.frequency;
    value.transmission = traced.volumetric_hit ? std::exp(-double(traced.optical_depth)) : 1;
    for (int row = 0; row < 2; ++row)
        for (int column = 0; column < 2; ++column)
            for (int angular = 0; angular < 2; ++angular)
                for (int axis = 0; axis < 2; ++axis)
                    value.source_derivative[row][column] +=
                        sky.map.jacobian[row][angular] *
                        film.differential->angular_jacobian[angular][axis] *
                        film_from_standard[axis][column];
    return value;
}

std::expected<PointDetectorBatchResult, PointDetectorError> SharedPointSources(
    const core::ICamera& camera, GeodesicTracer& tracer,
    const core::StarfieldSpatialIndex& catalogue, double brightness, double sigma, int x, int y,
    int width, int height, const CameraSample& sample, const std::function<bool()>& cancelled) {
    std::array<PointDetectorFootprint, 16> footprints;
    for (int i = 0; i < width * height; ++i) {
        const int dx = i % width, dy = i / width;
        const auto film = camera.ProjectFilmForObserver(x + dx + double(sample.image_u),
                                                        y + dy + double(sample.image_v),
                                                        sample.pupil_u, sample.pupil_v);
        if (!film || !film->ray.active || !film->differential)
            return std::unexpected(
                PointDetectorError{PointDetectorFailure::ProjectionUnavailable, {}});
        const auto& p = film->differential->angular_jacobian;
        const double determinant = std::fma(p[0][0], p[1][1], -p[0][1] * p[1][0]);
        if (!std::isfinite(determinant) || determinant == 0)
            return std::unexpected(PointDetectorError{PointDetectorFailure::Unresolved, {}});
        footprints[i] = {{double(dx), double(dy)},
                         {{{sigma * p[1][1] / determinant, -sigma * p[0][1] / determinant},
                           {-sigma * p[1][0] / determinant, sigma * p[0][0] / determinant}}}};
    }
    const PointDetectorSampler probe = [&](const DetectorCoordinate& q)
        -> std::expected<PointDetectorProbe, PointDetectorFailure> {
        const auto film = camera.ProjectFilmOffsetForObserver(x + double(sample.image_u),
                                                              y + double(sample.image_v), q[0],
                                                              q[1], sample.pupil_u, sample.pupil_v);
        if (!film) return std::unexpected(PointDetectorFailure::ProjectionUnavailable);
        if (!film->ray.active) return PointDetectorProbe{};
        // Derivative per common film coordinate. The batch engine composes
        // the discovery chart and each original Gaussian separately.
        return MeasurePointDetector(*film, tracer.TracePointSource(film->ray), {{{1, 0}, {0, 1}}});
    };
    return EvaluatePointDetectorBatch(
        catalogue, brightness,
        std::span(footprints).first(static_cast<std::size_t>(width * height)), probe, cancelled);
}
}  // namespace

bool RenderSession::UsesPhysicalPointDetector() const {
    const auto* family = dynamic_cast<const core::KerrSchildFamily*>(metric_.get());
    return config_.point_starfield && star_index_ && family && family->GetParams().Q == 0 &&
           family->GetParams().Lambda == 0;
}

base::Expected<RenderSession::PixelBlock> RenderSession::ShadeBlock(int x, int y, int width,
                                                                    int height,
                                                                    GeodesicTracer* tracer) const {
    PixelBlock result{x, y, width, height};
    std::optional<base::Error> sample_error;
    int sample_index = 0;
    const bool physical_point_detector = UsesPhysicalPointDetector();
    const int samples_taken =
        ForEachCameraSample(config_.samples_per_pixel, [&](const CameraSample& sample) {
            const int current_sample = sample_index++;
            if (sample_error) return;
            // All pixels use this sample's original pupil and image offset. A
            // declined envelope has no partial result: retry the individual
            // footprints, whose smaller support may still be represented.
            std::optional<PointDetectorBatchResult> shared;
            if (physical_point_detector && width * height > 1 && !IsStopping()) {
                auto batch =
                    SharedPointSources(*camera_, *tracer, *star_index_,
                                       config_.point_starfield_config.brightness_scale,
                                       pixel_angular_size_ * (config_.ray_bundles ? 1.0 : .3), x, y,
                                       width, height, sample, [&] { return IsStopping(); });
                if (batch) {
                    shared = std::move(*batch);
                } else if (batch.error().reason == PointDetectorFailure::Cancelled) {
                    sample_error.emplace(base::ErrorDomain::kPhysics, "shade block",
                                         "render cancelled during shared point discovery");
                    return;
                }
            }
            for (int pixel_index = 0; pixel_index < width * height; ++pixel_index) {
                const int px_coord = x + pixel_index % width;
                const int py_coord = y + pixel_index / width;
                // Keep a sample's early-return semantics local to its pixel. Every
                // pixel accumulates float radiance in the original SPP order.
                const auto shade_sample = [&] {
                    const auto fail_sample = [&](const std::string& reason) {
                        sample_error.emplace(base::ErrorDomain::kPhysics, "shade pixel",
                                             std::format("pixel ({}, {}), sample {}: {}", px_coord,
                                                         py_coord, current_sample, reason));
                    };
                    if (IsStopping()) {
                        fail_sample("render cancelled");
                        return;
                    }
                    const auto projection = camera_->ProjectFilmForObserver(
                        static_cast<double>(px_coord) + sample.image_u,
                        static_cast<double>(py_coord) + sample.image_v, sample.pupil_u,
                        sample.pupil_v);
                    if (!projection) {
                        fail_sample("camera projection is not represented");
                        return;
                    }
                    const CameraRay& camera_ray = projection->ray;
                    SIRIUS_ASSERT(core::IsRepresentedCameraRay(camera_ray));
                    if (!camera_ray.active) return;
                    // The packet centre and all its offset probes must consume the same
                    // infinity map. Their surface/volume contributions still finish before
                    // an outward vacuum handoff can succeed.
                    TraceResult trace_result = physical_point_detector
                                                   ? tracer->TracePointSource(camera_ray)
                                                   : tracer->Trace(camera_ray);
                    if (trace_result.cancelled) {
                        fail_sample("ray cancelled");
                        return;
                    }
                    if (trace_result.numerical_failure) {
                        const char* reason =
                            trace_result.coupled_failure == core::CoupledStepFailure::WorkLimit
                                ? "ray work limit exhausted"
                                : "numerical ray failure";
                        fail_sample(std::format(
                            "{} (coupled failure: {}, integrator termination: {}, attempts: {}, "
                            "accepted affine distance: {})",
                            reason, core::CoupledStepFailureName(trace_result.coupled_failure),
                            trace_result.integrator_termination, trace_result.steps_taken,
                            trace_result.affine_length));
                        return;
                    }
                    if (trace_result.outcome == TraceResult::Outcome::Escaped) {
                        for (int component = 0; component < 4; ++component) {
                            if (!std::isfinite(trace_result.final_direction(component))) {
                                fail_sample("non-finite escaped direction");
                                return;
                            }
                        }
                    }
                    if (trace_result.volumetric_hit && !std::isfinite(trace_result.optical_depth)) {
                        fail_sample("non-finite optical depth");
                        return;
                    }

                    float sr = 0.0f, sg = 0.0f, sb = 0.0f;

                    switch (trace_result.outcome) {
                        case TraceResult::Outcome::Horizon:
                        case TraceResult::Outcome::Throat:
                            break;

                        case TraceResult::Outcome::DiskHit: {
                            PixelResult disk = ShadeDiskHit(trace_result);
                            sr = disk.r;
                            sg = disk.g;
                            sb = disk.b;
                            break;
                        }

                        case TraceResult::Outcome::Escaped: {
                            if (!physical_point_detector) {
                                PixelResult esc = ShadeEscaped(trace_result);
                                sr = esc.r;
                                sg = esc.g;
                                sb = esc.b;
                            }
                            break;
                        }

                        case TraceResult::Outcome::MaxSteps:
                            // An unfinished ray has no terminal background to compose
                            // with accumulated volume emission, even if its producer
                            // omitted the numerical-failure flag.
                            fail_sample("ray work limit exhausted");
                            return;
                        default:
                            SIRIUS_ASSERT(false);
                            sr = 1.0f;
                            sg = 0.0f;
                            sb = 1.0f;
                            break;
                    }

                    // Volumetric transfer composes with the terminal surface/background;
                    // it is not a terminal ray outcome. Apply I = I_bg exp(-tau) + I_vol
                    // after shading the actual fate of the central ray.
                    if (trace_result.volumetric_hit) {
                        PixelResult volume = ShadeDiskHit(trace_result);
                        const float transmission =
                            std::exp(-std::max(trace_result.optical_depth, 0.0f));
                        sr = sr * transmission + volume.r;
                        sg = sg * transmission + volume.g;
                        sb = sb * transmission + volume.b;
                    }

                    if (physical_point_detector && shared) {
                        sr += static_cast<float>(shared->samples[pixel_index].rgb[0]);
                        sg += static_cast<float>(shared->samples[pixel_index].rgb[1]);
                        sb += static_cast<float>(shared->samples[pixel_index].rgb[2]);
                    } else if (physical_point_detector) {
                        if (!projection->differential) {
                            fail_sample("point detector camera differential is unavailable");
                            return;
                        }
                        // Freeze the original angular packet in its smooth film chart.
                        // All refinement uses this same L, pupil and Gaussian support.
                        const auto& p = projection->differential->angular_jacobian;
                        const double determinant = std::fma(p[0][0], p[1][1], -p[0][1] * p[1][0]);
                        const double sigma = pixel_angular_size_ * (config_.ray_bundles ? 1.0 : .3);
                        if (!std::isfinite(determinant) || determinant == 0) {
                            fail_sample("point detector camera map is singular");
                            return;
                        }
                        const core::AngularMatrix2 film_from_standard{
                            {{sigma * p[1][1] / determinant, -sigma * p[0][1] / determinant},
                             {-sigma * p[1][0] / determinant, sigma * p[0][0] / determinant}}};
                        const auto measure = [&](const core::CameraFilmProjection& film,
                                                 const TraceResult& traced) {
                            return MeasurePointDetector(film, traced, film_from_standard);
                        };
                        const PointDetectorSampler probe = [&](const DetectorCoordinate& z)
                            -> std::expected<PointDetectorProbe, PointDetectorFailure> {
                            if (z == DetectorCoordinate{})
                                return measure(*projection, trace_result);
                            const double dx =
                                film_from_standard[0][0] * z[0] + film_from_standard[0][1] * z[1];
                            const double dy =
                                film_from_standard[1][0] * z[0] + film_from_standard[1][1] * z[1];
                            const auto film = camera_->ProjectFilmOffsetForObserver(
                                static_cast<double>(px_coord) + sample.image_u,
                                static_cast<double>(py_coord) + sample.image_v, dx, dy,
                                sample.pupil_u, sample.pupil_v);
                            if (!film)
                                return std::unexpected(PointDetectorFailure::ProjectionUnavailable);
                            if (!film->ray.active) return PointDetectorProbe{};
                            return measure(*film, tracer->TracePointSource(film->ray));
                        };
                        const auto detector = EvaluatePointDetector(
                            *star_index_, config_.point_starfield_config.brightness_scale, probe,
                            [&] { return IsStopping(); });
                        if (!detector) {
                            fail_sample(
                                std::format("point detector failure {} after {} probes, {} cells "
                                            "and {} candidate visits",
                                            static_cast<int>(detector.error().reason),
                                            detector.error().statistics.probes,
                                            detector.error().statistics.cells,
                                            detector.error().statistics.candidate_visits));
                            return;
                        }
                        // Image-specific transmission is already included. The central
                        // volume attenuation above must not be applied a second time.
                        sr += static_cast<float>(detector->rgb[0]);
                        sg += static_cast<float>(detector->rgb[1]);
                        sb += static_cast<float>(detector->rgb[2]);
                    }

                    if (!std::isfinite(sr) || !std::isfinite(sg) || !std::isfinite(sb)) {
                        fail_sample("non-finite sample radiance");
                        return;
                    }
                    auto& accumulated = result.pixels[pixel_index];
                    accumulated.r += sr;
                    accumulated.g += sg;
                    accumulated.b += sb;
                    if (!std::isfinite(accumulated.r) || !std::isfinite(accumulated.g) ||
                        !std::isfinite(accumulated.b)) {
                        fail_sample("non-finite accumulated radiance");
                    }
                };
                shade_sample();
                if (sample_error) return;
            }
        });
    if (sample_error) return std::unexpected(std::move(*sample_error));
    const float inv_samples = 1.0f / static_cast<float>(samples_taken);
    for (int i = 0; i < width * height; ++i) {
        result.pixels[i].r *= inv_samples;
        result.pixels[i].g *= inv_samples;
        result.pixels[i].b *= inv_samples;
    }
    return result;
}

base::Expected<std::vector<float>> RenderSession::ShadeTile(const Tile& tile,
                                                            GeodesicTracer* tracer,
                                                            PixelBlock& cache) const {
    std::vector<float> pixels(static_cast<std::size_t>(tile.width) * tile.height * 4, 0);
    const int edge = UsesPhysicalPointDetector() ? 4 : 1;
    // Walk block intersections rather than image rows so the bounded cache is
    // sufficient even at high SPP. Anchors always refer to the complete frame.
    for (int y = (tile.y / edge) * edge; y < tile.y + tile.height; y += edge) {
        for (int x = (tile.x / edge) * edge; x < tile.x + tile.width; x += edge) {
            if (IsStopping())
                return base::Fail(base::ErrorDomain::kPhysics, "shade tile", "render cancelled");
            if (cache.x != x || cache.y != y) {
                auto block = ShadeBlock(x, y, std::min(edge, config_.width - x),
                                        std::min(edge, config_.height - y), tracer);
                if (!block) return std::unexpected(block.error());
                cache = std::move(*block);
            }
            for (int py = std::max(y, tile.y); py < std::min(y + edge, tile.y + tile.height);
                 ++py) {
                for (int px = std::max(x, tile.x); px < std::min(x + edge, tile.x + tile.width);
                     ++px) {
                    const auto& pixel = cache.pixels[(py - y) * cache.width + px - x];
                    const auto index =
                        (static_cast<std::size_t>(py - tile.y) * tile.width + px - tile.x) * 4;
                    pixels[index] = pixel.r;
                    pixels[index + 1] = pixel.g;
                    pixels[index + 2] = pixel.b;
                    pixels[index + 3] = 1;
                }
            }
        }
    }
    return pixels;
}

// =============================================================================
// Starfield background sampling (equirectangular projection).
// =============================================================================
void RenderSession::SampleStarfield(const Vec4& direction, float& r, float& g, float& b) const {
    if (!starfield_loaded_ || starfield_data_.empty()) {
        SIRIUS_ASSERT(starfield_loaded_ && !starfield_data_.empty());
        r = g = b = 0.0f;
        return;
    }

    // Direction is Cartesian (x, y, z) from the geodesic tracer.
    double dx = direction(1);
    double dy = direction(2);
    double dz = direction(3);

    double len = std::sqrt(dx * dx + dy * dy + dz * dz);
    if (len < 1e-10) {
        SIRIUS_ASSERT(len >= 1e-10);
        r = g = b = 0.0f;
        return;
    }

    dx /= len;
    dy /= len;
    dz /= len;

    // Spherical coordinates: theta 0 at +Z to pi at -Z; phi 0 at +X toward +Y.
    double theta = std::acos(std::clamp(dz, -1.0, 1.0));
    double phi = std::atan2(dy, dx);
    if (phi < 0) phi += math::kTwoPi;

    double u = phi / math::kTwoPi;
    double v = theta / math::kPi;

    double px = u * (starfield_width_ - 1);
    double py = v * (starfield_height_ - 1);

    int x0 = static_cast<int>(std::floor(px));
    int y0 = static_cast<int>(std::floor(py));
    int x1 = std::min(x0 + 1, starfield_width_ - 1);
    int y1 = std::min(y0 + 1, starfield_height_ - 1);

    double fx = px - x0;
    double fy = py - y0;

    auto sample = [this](int x, int y) -> std::array<float, 3> {
        int idx = (y * starfield_width_ + x) * 4;  // 4 bytes per pixel (RGBA).
        return {starfield_data_[idx + 0] / 255.0f, starfield_data_[idx + 1] / 255.0f,
                starfield_data_[idx + 2] / 255.0f};
    };

    auto c00 = sample(x0, y0);
    auto c10 = sample(x1, y0);
    auto c01 = sample(x0, y1);
    auto c11 = sample(x1, y1);

    float w00 = static_cast<float>((1.0 - fx) * (1.0 - fy));
    float w10 = static_cast<float>(fx * (1.0 - fy));
    float w01 = static_cast<float>((1.0 - fx) * fy);
    float w11 = static_cast<float>(fx * fy);

    r = c00[0] * w00 + c10[0] * w10 + c01[0] * w01 + c11[0] * w11;
    g = c00[1] * w00 + c10[1] * w10 + c01[1] * w01 + c11[1] * w11;
    b = c00[2] * w00 + c10[2] * w10 + c01[2] * w01 + c11[2] * w11;
}

// =============================================================================
// Filtered point-source star field sampling (P3).
// =============================================================================
void RenderSession::SampleStarfieldPoints(const TraceResult& result, float& r, float& g,
                                          float& b) const {
    // Beam footprint on the sky. With ray bundles the tracer supplies the full
    // lensed ellipse; a pinhole (bundles off) samples at a fraction of the pixel
    // angular size, so a star pops in and out as the camera rotates (the flicker
    // the beam filter removes). Both ellipse axes are floored at the pixel size so
    // a pixel always integrates at least its own solid angle.
    constexpr float kPinholeFraction = 0.3f;  // Pinhole sigma as a fraction of a pixel.
    float pixel = static_cast<float>(pixel_angular_size_);
    float sigma_major;
    float sigma_minor;
    float orientation = 0.0f;
    if (config_.ray_bundles && result.beam.valid) {
        sigma_major = std::max(result.beam.footprint_major, pixel);
        sigma_minor = std::max(result.beam.footprint_minor, pixel);
        orientation = result.beam.orientation;
    } else {
        sigma_major = kPinholeFraction * pixel;
        sigma_minor = sigma_major;
    }

    const auto& d = result.final_direction;
    star_generator_->AccumulateThroughBeam(static_cast<float>(d(1)), static_cast<float>(d(2)),
                                           static_cast<float>(d(3)), sigma_major, sigma_minor,
                                           orientation, *star_index_, r, g, b);
}

}  // namespace sirius::render
