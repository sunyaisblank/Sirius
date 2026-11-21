#include "sirius/render/session/kerr_point_detector.h"

#include "sirius/backend/trace_step_executor.h"
#include "sirius/core/camera.h"
#include "sirius/core/kerr_detector/finite_screen.h"

#include <mutex>
#include <new>

namespace sirius::render {
namespace {
using namespace core::relativity::finite_tail_detail;

struct SourceBand {
    Range weighted_temperature;
    double clearance = kInfinity;
};

std::optional<std::array<SourceBand, 3>> CollectReachableSources(
    const AngularCatalogue& catalogue, Range amplitude, std::size_t& remaining,
    const std::function<bool()>& cancelled) {
    std::array<SourceBand, 3> result{};
    const auto add = [&](unsigned group, Range weight, double clearance) {
        auto& band = result[group];
        band.weighted_temperature = band.weighted_temperature + weight;
        band.clearance = std::min(band.clearance, clearance);
    };
    for (bool north : {false, true}) {
        for (const auto& row : catalogue.Rows(north)) {
            if (cancelled && cancelled()) return {};
            if (row.entries.empty() || row.latitude.lo > amplitude.hi) continue;
            const Range gap = amplitude - row.latitude;
            if (gap.lo > 0.0 && LatitudeBand(gap.lo) == LatitudeBand(gap.hi)) {
                add(LatitudeBand(gap.lo), row.weighted_temperature, gap.lo);
                continue;
            }
            // Only rows touching a band edge or a polar turn require individual
            // visits. A source which may touch a turn needs a different proof.
            for (const auto& entry : row.entries) {
                if (remaining == 0) return {};
                --remaining;
                if (entry.latitude.lo > amplitude.hi) continue;
                const Range star_gap = amplitude - entry.latitude;
                if (!(star_gap.lo > 0.0)) return {};
                add(LatitudeBand(star_gap.lo), entry.weighted_temperature, star_gap.lo);
            }
        }
    }
    return result;
}

std::optional<StationaryChart> EvaluateFilmChart(const core::CameraConfig& camera,
                                                 const core::KerrSchildParams& metric,
                                                 const backend::DetectorLaunchFrame& launch,
                                                 const DetectorCoordinate& origin,
                                                 const PointDetectorFootprint& footprint,
                                                 const FilmVector& box) {
    const FilmJet z0(box[0], {Range(1.0), Range(0.0)});
    const FilmJet z1(box[1], {Range(0.0), Range(1.0)});
    const auto film = [&](unsigned axis) {
        const auto& row = footprint.chart_from_standard[axis];
        return FilmJet(origin[axis]) + FilmJet(footprint.centre[axis]) +
               (FilmJet(row[0]) * z0 + FilmJet(row[1]) * z1);
    };
    const auto& pupil = launch.pupil;
    const auto q = core::camera_detail::ThinLensFilmVector(
        core::camera_detail::FilmPlaneX(film(0), camera),
        core::camera_detail::FilmPlaneY(film(1), camera), core::camera_detail::FilmTangent(camera),
        camera.focus_distance, pupil.aperture_up, pupil.aperture_right);
    const FilmJet norm = sqrt(Square(q[0]) + Square(q[1]) + Square(q[2]));
    const std::array direction{q[0] / norm, q[1] / norm, q[2] / norm};
    std::array<FilmJet, 4> event, tangent;
    for (unsigned i = 0; i < 4; ++i) {
        event[i] = FilmJet(launch.position[i]);
        tangent[i] = FilmJet(0.0) - FilmJet(launch.tetrad[0][i]);
        for (unsigned axis = 0; axis < 3; ++axis)
            tangent[i] = tangent[i] + FilmJet(launch.tetrad[axis + 1][i]) * direction[axis];
    }
    const auto state = core::relativity::kerr_infinity_detail::EncloseDetectorRadialState<true>(
        metric, event, tangent);
    return state ? MakeStationaryChart(*state, metric.a) : std::nullopt;
}

template <typename Evaluate>
FilmVector CoverEscapingChart(const FilmVector& cell, const StationaryChart& whole,
                              const Evaluate& evaluate) {
    auto target = whole.coordinate;
    const FilmVector middle{Range(std::midpoint(cell[0].lo, cell[0].hi)),
                            Range(std::midpoint(cell[1].lo, cell[1].hi))};
    if (const auto anchor = evaluate(middle)) {
        // Use a common film point so the radius and minimum enclosures
        // retain the correlations lost in direct interval evaluation.
        for (unsigned i = 0; i < 2; ++i) {
            const Range anchored = anchor->coordinate[i] +
                                   whole.jacobian[i][0] * (cell[0] - middle[0]) +
                                   whole.jacobian[i][1] * (cell[1] - middle[1]);
            if (anchored.Finite())
                target[i] =
                    Range(std::max(target[i].lo, anchored.lo), std::min(target[i].hi, anchored.hi));
        }
    }
    const auto& [radial, potential] = whole.jacobian;
    std::array<double, 2> shallow, deep;
    std::array<bool, 2> increasing;
    for (unsigned axis = 0; axis < 2; ++axis) {
        if (!(potential[axis].lo > 0.0 || potential[axis].hi < 0.0) ||
            !(radial[axis].lo > 0.0 || radial[axis].hi < 0.0))
            return target;
        const bool reversed = potential[axis].lo > 0.0;
        shallow[axis] = reversed ? cell[axis].hi : cell[axis].lo;
        deep[axis] = reversed ? cell[axis].lo : cell[axis].hi;
        increasing[axis] = (radial[axis].lo > 0.0) != reversed;
    }
    const FilmVector deepest{Range(deep[0]), Range(deep[1])};
    if (const auto corner = evaluate(deepest))
        target[1].lo = std::max(target[1].lo, corner->coordinate[1].lo);
    if (increasing[0] == increasing[1]) return target;

    // Along either axis the chosen deep direction decreases the potential.
    // Extrema on the escaping side occur on an edge with the other axis
    // moved fully toward the extremum, then at the edge's shadow crossing.
    const auto crossing = [&](unsigned axis) -> std::optional<Range> {
        auto edge = deepest;
        const auto potential_at = [&](double value) {
            edge[axis] = Range(value);
            const auto point = evaluate(edge);
            return point ? point->coordinate[1] : InvalidRange();
        };
        if (!(potential_at(shallow[axis]).lo > 0.0 && potential_at(deep[axis]).hi < 0.0)) return {};
        const auto tighten = [&](double known, double other, bool positive) {
            for (unsigned i = 0; i < 45; ++i) {
                const double mid = std::midpoint(known, other);
                if (mid == known || mid == other) break;
                const auto value = potential_at(mid);
                if (positive ? value.lo > 0.0 : value.hi < 0.0)
                    known = mid;
                else
                    other = mid;
            }
            return known;
        };
        const double left = tighten(shallow[axis], deep[axis], true);
        const double right = tighten(deep[axis], left, false);
        edge[axis] = Range(std::min(left, right), std::max(left, right));
        const auto enclosure = evaluate(edge);
        return enclosure ? std::optional(enclosure->coordinate[0]) : std::nullopt;
    };
    const unsigned increasing_axis = increasing[0] ? 0 : 1;
    if (const auto edge = crossing(increasing_axis))
        target[0].lo = std::max(target[0].lo, edge->lo);
    if (const auto edge = crossing(1 - increasing_axis))
        target[0].hi = std::min(target[0].hi, edge->hi);
    return target;
}

Range BoundGaussian(const FilmVector& box) {
    const Range radius2 = Square(box[0]) + Square(box[1]);
    const Range weight = ExpPoint((Range(-0.5) * Range(radius2.lo)).hi);
    // The original Gaussian integral is 2*pi*(1 - exp(-8)) > 6.
    return Range(weight.hi) / Range(6.0);
}

std::optional<Range> ColourEnvelope(Range camera_frequency, Range gaussian, double brightness) {
    // Rayleigh-Jeans bounds Planck at every represented wavelength. Use
    // the absolute CIE envelope and the D65 matrix row norm, with margin
    // for the represented spectral arithmetic.
    constexpr double first = core::constants::spectral::kLambdaMin;
    constexpr double last = core::constants::spectral::kLambdaMax;
    constexpr double step = (last - first) / 32.0;
    for (unsigned i = 0; i < 32; ++i) {
        const auto weight = core::spectral::WavelengthToXyz((first + (i + .5) * step) * 1e9);
        if (!(std::abs(weight.X) < 2.0 && std::abs(weight.Y) < 2.0 && std::abs(weight.Z) < 2.0))
            return {};
    }
    const Range result =
        gaussian * Range(64.0) * Range(32.0) * Range(step) *
        Range(core::constants::physical::kPlanckC1) * camera_frequency * Range(brightness) /
        (Range(core::constants::physical::kPlanckC2) * Square(Square(Range(first))));
    return result.Finite() ? std::optional(result) : std::nullopt;
}

std::optional<Range> BoundAzimuthImages(double mass, double spin, Range centre, Range scale,
                                        const StationaryChart& chart, const CatalogueScreen& screen,
                                        const PolarIntegralTable& polar, FiniteCandidate candidate,
                                        const std::function<bool()>& cancelled,
                                        std::size_t& remaining,
                                        const NearCriticalPhases* shared_phases = nullptr) {
    if (cancelled && cancelled()) return {};
    const auto& star = *candidate.source;
    std::optional<RefinedCandidate> refined;
    if (shared_phases && candidate.log_scale.hi <= shared_phases->maximum_log_scale)
        refined = shared_phases->Refine(candidate);
    const bool used_shared_phases = refined.has_value();
    if (!refined)
        refined = RefineFiniteAzimuth(mass, spin, centre, scale, chart, candidate, star.latitude,
                                      star.azimuth);
    if (!refined) return {};
    if (!refined->possible) return Range{};
    // This verifies the latitude root exists continuously over radius.
    // Thus a sign-definite azimuth slope permits at most one root per wrap.
    if (!EntireLatitudeBranchExists(screen, polar, candidate)) return {};
    const Range wraps = (refined->azimuth - star.azimuth) / (Range(2.0) * kPiBounds);
    const Range separation = ExpBounds(refined->log_scale);
    if (!wraps.Finite() || !separation.Finite()) return {};
    const double first = std::ceil(wraps.lo), last = std::floor(wraps.hi);
    if (first > last) return Range{};
    const auto derivative = DifferentiateFiniteAzimuth(
        mass, spin, centre, Range(separation.lo, std::min(separation.hi, scale.hi)),
        refined->log_scale, chart, candidate.winding, candidate.approaching_equator, star.latitude);
    if (derivative) {
        const Range slope = Abs(derivative->conditioned_radial_slope);
        if (slope.lo > 0.0) {
            const Range roots = Range(last) - Range(first) + Range(1.0);
            const Range density =
                roots * Range(2.0) * Square(separation) /
                (slope * derivative->latitude_log_slope_magnitude * derivative->mu_speed);
            if (density.Finite()) return density;
        }
    }

    candidate.log_scale = refined->log_scale;
    if (used_shared_phases) {
        if (remaining == 0) return {};
        --remaining;
        return BoundAzimuthImages(mass, spin, centre, scale, chart, screen, polar, candidate,
                                  cancelled, remaining);
    }

    // An inconclusive derivative need not imply a caustic image. Retest
    // longitude on each radius half before requesting a derivative there.
    // Closed halves may count a shared endpoint twice, conservatively.
    const double middle = std::midpoint(centre.lo, centre.hi);
    if (!(centre.lo < middle && middle < centre.hi)) return {};
    Range total;
    for (Range half : {Range(centre.lo, middle), Range(middle, centre.hi)}) {
        if (remaining == 0) return {};
        --remaining;
        const auto bound = BoundAzimuthImages(mass, spin, half, scale, chart, screen, polar,
                                              candidate, cancelled, remaining);
        if (!bound) return {};
        total = total + *bound;
    }
    return total.Finite() ? std::optional(total) : std::nullopt;
}
}  // namespace

class KerrDetectorCatalogue {
  public:
    explicit KerrDetectorCatalogue(std::span<const core::StarEntry> entries) : stars_(entries) {}

    const AngularCatalogue* Get(const std::function<bool()>& cancelled) const {
        std::call_once(prepare_, [&] {
            try {
                angular_ = AngularCatalogue::Build(stars_, cancelled);
            } catch (const std::bad_alloc&) {
                // The optional enclosure cannot crowd out the original detector.
            }
        });
        return angular_ ? &*angular_ : nullptr;
    }

  private:
    std::span<const core::StarEntry> stars_;
    mutable std::once_flag prepare_;
    mutable std::optional<AngularCatalogue> angular_;
};

std::shared_ptr<const KerrDetectorCatalogue> PrepareKerrDetectorCatalogue(
    const core::StarfieldSpatialIndex& stars) {
    return std::make_shared<KerrDetectorCatalogue>(stars.Stars());
}

std::optional<std::array<double, 3>> BoundKerrDetectorCell(
    const KerrDetectorCatalogue& catalogue, double brightness, const core::CameraConfig& camera,
    const backend::DetectorLaunchFrame& launch, const DetectorCoordinate& film_origin,
    const PointDetectorFootprint& footprint, const DetectorCoordinate& lower,
    const DetectorCoordinate& upper, const std::function<bool()>& cancelled,
    std::size_t& remaining) {
    using core::KerrSchildParams;
    const auto& [mass, spin, charge, lambda] = launch.metric;
    if (!(mass > spin && spin > 0.0) || charge != 0.0 || lambda != 0.0 || !RoundingEnvironment() ||
        (cancelled && cancelled()))
        return {};
    const KerrSchildParams metric{mass, spin, charge, lambda};
    FilmVector box{Range(lower[0], upper[0]), Range(lower[1], upper[1])};
    if (!box[0].Finite() || !box[1].Finite()) return {};
    // Tighten the rectangle to the axis-aligned hull of its intersection
    // with the original Gaussian disk. At its rim the active part can be
    // much narrower than the cell, and only that part needs a sky bound.
    const std::array nearest{std::clamp(0.0, box[0].lo, box[0].hi),
                             std::clamp(0.0, box[1].lo, box[1].hi)};
    for (unsigned axis = 0; axis < 2; ++axis) {
        const Range extent2 = Range(16.0) - Square(Range(nearest[1 - axis]));
        if (extent2.hi < 0.0) return std::array<double, 3>{};
        const double extent = sqrt(Range(extent2.hi)).hi;
        box[axis] = Range(std::max(box[axis].lo, -extent), std::min(box[axis].hi, extent));
        if (box[axis].hi < box[axis].lo) return std::array<double, 3>{};
    }
    const auto evaluate = [&](const FilmVector& region) {
        return EvaluateFilmChart(camera, metric, launch, film_origin, footprint, region);
    };
    const auto whole = evaluate(box);
    if (!whole || !(whole->initial_mu.lo > 0.0 && whole->mu_rate.hi < 0.0) ||
        !ProveInjectiveChart(whole->jacobian))
        return {};
    const auto& [radial_gradient, potential_gradient] = whole->jacobian;
    const Range determinant = Abs(radial_gradient[0] * potential_gradient[1] -
                                  radial_gradient[1] * potential_gradient[0]);
    const auto target = CoverEscapingChart(box, *whole, evaluate);
    if (!target[0].Finite() || !target[1].Finite()) return {};
    // MakeStationaryChart admits the only positive radial minimum and an
    // inward launch. Positive bounds there and at both ends rule out a
    // turning point throughout the supported film rectangle.
    if (target[1].lo > 0.0 && whole->positive_radial_endpoints) return std::array<double, 3>{};
    if (!(determinant.lo > 0.0 && target[1].lo < 0.0)) return {};
    const Range scale(0.0, sqrt(Range(-target[1].lo)).hi);

    const auto phase =
        EncloseRegularizedRadialPhase(mass, spin, target[0], scale, whole->inverse_start, 128);
    const auto shape = MakeFinitePolarShape(mass, spin, target[0], scale);
    if (!phase || !shape) return {};
    const Range amplitude = sqrt(shape->amplitude2);
    const auto* angular = catalogue.Get(cancelled);
    if (!angular) return {};
    const auto bands = CollectReachableSources(*angular, amplitude, remaining, cancelled);
    if (!bands) return {};
    double clearance = kInfinity;
    for (const auto& band : *bands) clearance = std::min(clearance, band.clearance);
    if (clearance == kInfinity) return std::array<double, 3>{};
    const auto latitude = EncloseLatitudeRootMonotonicity(mass, spin, target[0], scale,
                                                          whole->initial_mu, clearance, *phase);
    const auto colour =
        ColourEnvelope(Range(1.0) / (Range(0.0) - whole->energy), BoundGaussian(box), brightness);
    if (!latitude || !colour) return {};

    Range upper_bound;
    std::array<unsigned, 3> cutoffs{};
    for (std::size_t i = 0; i < bands->size(); ++i) {
        const auto& band = (*bands)[i];
        if (band.clearance == kInfinity) continue;
        const auto tail = EncloseWindingTail(mass, spin, target[0], scale, whole->initial_mu,
                                             band.clearance, determinant, *phase, *latitude, 8);
        if (!tail) return {};
        const double factor =
            (*colour * band.weighted_temperature * Range(tail->tail_prefactor)).hi;
        const double bound = std::nextafter(std::scalbn(factor, tail->tail_exponent), kInfinity);
        if (!std::isfinite(bound)) return {};
        upper_bound = upper_bound + Range(0.0, bound);
        cutoffs[i] = tail->first_winding;
    }

    const auto polar = MakePolarIntegralTable(*shape);
    if (!polar) return {};
    const auto phases = PrepareNearCriticalPhases(mass, spin, target[0], scale, *whole);
    const auto screen = ScreenFiniteWindings(mass, spin, target[0], scale, *whole, *polar, *phase,
                                             *angular, cutoffs, amplitude, clearance, remaining,
                                             cancelled, phases ? &*phases : nullptr);
    if (!screen) return {};
    for (const auto& candidate : screen->candidates) {
        const auto density =
            BoundAzimuthImages(mass, spin, target[0], scale, *whole, *screen, *polar, candidate,
                               cancelled, remaining, phases ? &*phases : nullptr);
        if (!density) return {};
        if (density->hi == 0.0) continue;
        const Range contribution =
            *colour * candidate.source->weighted_temperature * *density / determinant;
        if (!contribution.Finite() || contribution.hi < 0.0) return {};
        upper_bound = upper_bound + Range(0.0, contribution.hi);
    }
    if (!upper_bound.Finite() || (cancelled && cancelled())) return {};
    return std::array{upper_bound.hi, upper_bound.hi, upper_bound.hi};
}

}  // namespace sirius::render
