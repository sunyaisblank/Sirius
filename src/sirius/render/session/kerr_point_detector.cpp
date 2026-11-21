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

template <typename ExcludeTurn>
std::optional<std::array<SourceBand, 3>> CollectReachableSources(
    const AngularCatalogue& catalogue, Range amplitude, std::size_t& remaining,
    const std::function<bool()>& cancelled, const ExcludeTurn& exclude_turn) {
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
            // Only rows touching a band edge or a polar turn require individual visits.
            for (const auto& entry : row.entries) {
                if (remaining == 0) return {};
                --remaining;
                if (entry.latitude.lo > amplitude.hi) continue;
                const Range star_gap = amplitude - entry.latitude;
                if (!(star_gap.lo > 0.0)) {
                    if (!exclude_turn(entry, north)) return {};
                    continue;
                }
                add(LatitudeBand(star_gap.lo), entry.weighted_temperature, star_gap.lo);
            }
        }
    }
    return result;
}

struct FiniteTurnScreen {
    FiniteRadialPhase radial;
    Range radial_time, fixed_azimuth;
    PhasePair half, initial;
};

std::optional<FiniteTurnScreen> PrepareFiniteTurnScreen(double mass, double spin, Range radius,
                                                        Range scale, const StationaryChart& chart,
                                                        const FinitePolarShape& shape,
                                                        Range amplitude) {
    const auto radial = EvaluateFiniteRadialPhase(mass, spin, radius, scale, chart.inverse_start);
    if (!radial) return {};
    const Range time = Range(-2.0) * radial->leading * LogBounds(scale) + radial->time_residual;
    const Range azimuth = FixedPupilAzimuth(mass, spin, chart);
    const auto half = IntegratePolarPartial(shape, kPiBounds * Range(.5));
    const auto initial = IntegratePolarPartial(shape, AsinBounds(chart.initial_mu / amplitude));
    if (!time.Finite() || !azimuth.Finite() || !std::ranges::all_of(half, &Range::Finite) ||
        !std::ranges::all_of(initial, &Range::Finite))
        return {};
    return FiniteTurnScreen{*radial, time, azimuth, half, initial};
}

bool ExcludeFinitePolarTurn(const FiniteTurnScreen& screen, const FinitePolarShape& shape,
                            Range amplitude, const AngularCatalogue::Entry& source, bool north,
                            unsigned cutoff, std::size_t& remaining,
                            const std::function<bool()>& cancelled) {
    // The latitude restricts the endpoint on either side of the polar turn.
    // Reflection about pi/2 covers both legs without differentiating asin.
    const Range ratio = Range(source.latitude.lo) / Range(amplitude.hi);
    const Range lower_angle = AsinBounds(Range(ratio.lo));
    if (!lower_angle.Finite()) return false;
    const auto lower = IntegratePolarPartial(shape, Range(lower_angle.lo));
    const PhasePair endpoint{
        Range(lower[0].lo, (Range(2.0) * screen.half[0] - lower[0]).hi),
        Range(lower[1].lo, (Range(2.0) * screen.half[1] - lower[1]).hi),
    };
    if (!std::ranges::all_of(endpoint, &Range::Finite)) return false;
    for (unsigned m = north ? 1U : 0U; m < cutoff; m += 2) {
        if (remaining == 0 || (cancelled && cancelled())) return false;
        --remaining;
        const Range time = Range(2.0 * m) * screen.half[0] + screen.initial[0] + endpoint[0];
        if (!time.Finite()) return false;
        if (time.hi < screen.radial_time.lo || time.lo > screen.radial_time.hi) continue;
        const Range angle = Range(2.0 * m) * screen.half[1] + screen.initial[1] + endpoint[1];
        const Range longitude = screen.fixed_azimuth - (angle + screen.radial.drag_at_turn * time +
                                                        screen.radial.spin_residual);
        if (!longitude.Finite() || AzimuthMayMatch(longitude, source.azimuth)) return false;
    }
    return true;
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
        if (!(Abs(potential[axis]).lo > 0.0)) return target;
        const bool reversed = potential[axis].lo > 0.0;
        shallow[axis] = reversed ? cell[axis].hi : cell[axis].lo;
        deep[axis] = reversed ? cell[axis].lo : cell[axis].hi;
        increasing[axis] = (radial[axis].lo > 0.0) != reversed;
    }
    const FilmVector deepest{Range(deep[0]), Range(deep[1])};
    const auto lowest = evaluate(deepest);
    if (lowest) target[1].lo = std::max(target[1].lo, lowest->coordinate[1].lo);
    if (target[1].lo >= 0.0 || target[1].hi < 0.0) return target;
    const FilmVector shallowest{Range(shallow[0]), Range(shallow[1])};
    if (const auto highest = evaluate(shallowest))
        target[1].hi = std::min(target[1].hi, highest->coordinate[1].hi);
    if (!lowest || !(lowest->coordinate[1].hi < 0.0) || target[1].hi < 0.0 ||
        !(Abs(radial[0]).lo > 0.0 && Abs(radial[1]).lo > 0.0) || increasing[0] == increasing[1])
        return target;

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
        if (!(potential_at(shallow[axis]).lo > 0.0)) return {};
        Range root = cell[axis];
        const bool sign_is_positive = potential[axis].lo > 0.0;
        const Range slope = Abs(potential[axis]);
        // The bracket contains a crossing and the whole-edge derivative
        // excludes zero. Intersect it with the mean-value enclosure of the
        // root, retaining both signs of the pointwise rounding uncertainty.
        for (unsigned i = 0; i < 8; ++i) {
            const double mid = std::midpoint(root.lo, root.hi);
            const Range value = potential_at(mid);
            const Range oriented = sign_is_positive ? value : Range(0.0) - value;
            const Range image = Range(mid) - oriented / slope;
            if (!image.Finite()) break;
            const Range narrowed(std::max(root.lo, image.lo), std::min(root.hi, image.hi));
            if (!narrowed.Finite()) return {};
            if (narrowed.lo == root.lo && narrowed.hi == root.hi) break;
            root = narrowed;
        }
        edge[axis] = root;
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
    // Bound the Rayleigh-Jeans spectrum at the same 32 wavelengths as
    // point-source transfer. Absolute XYZ and D65 coefficients bound each
    // linear-light channel even when the transform involves cancellation.
    constexpr double first = core::constants::spectral::kLambdaMin;
    constexpr double last = core::constants::spectral::kLambdaMax;
    constexpr double step = (last - first) / 32.0;
    std::array<Range, 3> xyz{};
    for (unsigned i = 0; i < 32; ++i) {
        const double wavelength = first + (i + .5) * step;
        const auto weight = core::spectral::WavelengthToXyz(wavelength * 1e9);
        const Range response = Range(step) / Square(Square(Range(wavelength)));
        xyz[0] = xyz[0] + response * Abs(Range(weight.X));
        xyz[1] = xyz[1] + response * Abs(Range(weight.Y));
        xyz[2] = xyz[2] + response * Abs(Range(weight.Z));
    }
    constexpr auto x = core::colour::XyzD65ToLinearSrgb(1.0, 0.0, 0.0);
    constexpr auto y = core::colour::XyzD65ToLinearSrgb(0.0, 1.0, 0.0);
    constexpr auto z = core::colour::XyzD65ToLinearSrgb(0.0, 0.0, 1.0);
    const auto channel = [&](double a, double b, double c) {
        return (xyz[0] * Abs(Range(a)) + xyz[1] * Abs(Range(b)) + xyz[2] * Abs(Range(c))).hi;
    };
    const double band =
        std::max({channel(x.r, y.r, z.r), channel(x.g, y.g, z.g), channel(x.b, y.b, z.b)});
    // Retain a factor-four margin over the absolute-coefficient envelope
    // for the represented spectral arithmetic.
    const Range result = gaussian * Range(4.0) * Range(band) *
                         Range(core::constants::physical::kPlanckC1) * camera_frequency *
                         Range(brightness) / Range(core::constants::physical::kPlanckC2);
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

std::optional<std::array<double, 3>> BoundKerrDetectorCells(
    const KerrDetectorCatalogue& catalogue, double brightness, const core::CameraConfig& camera,
    const backend::DetectorLaunchFrame& launch, const DetectorCoordinate& film_origin,
    const PointDetectorFootprint& footprint, std::span<const PointDetectorCell> cells,
    const std::function<bool()>& cancelled, std::size_t& remaining) {
    using core::KerrSchildParams;
    const auto& [mass, spin, charge, lambda] = launch.metric;
    if (!(mass > spin && spin > 0.0) || charge != 0.0 || lambda != 0.0 || !RoundingEnvironment())
        return {};
    const KerrSchildParams metric{mass, spin, charge, lambda};
    const auto evaluate = [&](const FilmVector& region) {
        return EvaluateFilmChart(camera, metric, launch, film_origin, footprint, region);
    };
    const auto hull = [](Range a, Range b) {
        return Range(std::min(a.lo, b.lo), std::max(a.hi, b.hi));
    };
    std::optional<FilmVector> film_hull;
    Range radius;
    double minimum_potential = 0.0;
    double maximum_potential = -kInfinity;
    double minimum_determinant = kInfinity;
    double maximum_gaussian = 0.0;
    for (const auto& cell : cells) {
        if (cancelled && cancelled()) return {};
        FilmVector box{Range(cell.lower[0], cell.upper[0]), Range(cell.lower[1], cell.upper[1])};
        if (!box[0].Finite() || !box[1].Finite()) return {};
        // Retain only the axis-aligned hull of the original Gaussian disk
        // in this cell. The omitted portion cannot reach a catalogue source.
        const std::array nearest{std::clamp(0.0, box[0].lo, box[0].hi),
                                 std::clamp(0.0, box[1].lo, box[1].hi)};
        for (unsigned axis = 0; axis < 2; ++axis) {
            const Range extent2 = Range(16.0) - Square(Range(nearest[1 - axis]));
            if (extent2.hi < 0.0) {
                box[axis] = InvalidRange();
                break;
            }
            const double extent = sqrt(Range(extent2.hi)).hi;
            box[axis] = Range(std::max(box[axis].lo, -extent), std::min(box[axis].hi, extent));
        }
        if (!box[0].Finite() || !box[1].Finite()) continue;
        const auto chart = evaluate(box);
        if (!chart) return {};
        const auto target = CoverEscapingChart(box, *chart, evaluate);
        if (!target[0].Finite() || !target[1].Finite()) return {};
        // A positive sole radial minimum and positive endpoints exclude
        // any turning point on this inward trajectory.
        if (target[1].lo > 0.0 && chart->positive_radial_endpoints) continue;
        const auto& [radial, potential] = chart->jacobian;
        const Range determinant = Abs(radial[0] * potential[1] - radial[1] * potential[0]);
        const Range gaussian = BoundGaussian(box);
        if (!(determinant.lo > 0.0 && target[1].lo < 0.0 && gaussian.Finite())) return {};
        radius = film_hull ? hull(radius, target[0]) : target[0];
        minimum_potential = std::min(minimum_potential, target[1].lo);
        maximum_potential = std::max(maximum_potential, target[1].hi);
        minimum_determinant = std::min(minimum_determinant, determinant.lo);
        maximum_gaussian = std::max(maximum_gaussian, gaussian.hi);
        film_hull = film_hull
                        ? FilmVector{hull((*film_hull)[0], box[0]), hull((*film_hull)[1], box[1])}
                        : box;
    }
    if (!film_hull) return std::array<double, 3>{};
    // Prove injectivity across the group so a sky root in the common
    // invariant rectangle can contribute at most once across these leaves.
    // The detailed per-cell charts supply the tighter radius, depth and gain.
    const auto whole = evaluate(*film_hull);
    if (!whole || !(whole->initial_mu.lo > 0.0 && whole->mu_rate.hi < 0.0) ||
        !ProveInjectiveChart(whole->jacobian))
        return {};
    const Range determinant(minimum_determinant);
    const Range scale(0.0, sqrt(Range(-minimum_potential)).hi);

    const auto phase =
        EncloseRegularizedRadialPhase(mass, spin, radius, scale, whole->inverse_start, 128);
    const auto shape = MakeFinitePolarShape(mass, spin, radius, scale);
    if (!phase || !shape) return {};
    const Range amplitude = sqrt(shape->amplitude2);

    // Each complete polar half-libration consumes at least the lower
    // period bound, even when a source can lie on a turning latitude.
    constexpr unsigned first_tail_winding = 8;
    Range finite_scale;
    std::optional<unsigned> finite_cutoff;
    if (maximum_potential < 0.0) {
        finite_scale = Range(sqrt(Range(-maximum_potential)).lo, scale.hi);
        const Range maximum_time = Range(-2.0) * phase->leading_factor * LogPoint(finite_scale.lo) +
                                   phase->time_remainder.value;
        const Range half_period =
            kPiBounds / sqrt(shape->omega2 + shape->modulation * Range(0.0, 1.0));
        const double maximum_winding = (maximum_time / half_period).hi;
        if (maximum_winding >= 0.0 && maximum_winding < first_tail_winding)
            finite_cutoff = static_cast<unsigned>(maximum_winding) + 1;
    }
    const auto* angular = catalogue.Get(cancelled);
    if (!angular) return {};
    std::optional<FiniteTurnScreen> turn_screen;
    std::vector<const AngularCatalogue::Entry*> excluded_turns;
    const auto exclude_turn = [&](const AngularCatalogue::Entry& source, bool north) {
        if (!finite_cutoff) return false;
        if (!turn_screen)
            turn_screen = PrepareFiniteTurnScreen(mass, spin, radius, finite_scale, *whole, *shape,
                                                  amplitude);
        if (!turn_screen || !ExcludeFinitePolarTurn(*turn_screen, *shape, amplitude, source, north,
                                                    *finite_cutoff, remaining, cancelled))
            return false;
        excluded_turns.push_back(&source);
        return true;
    };
    const auto bands =
        CollectReachableSources(*angular, amplitude, remaining, cancelled, exclude_turn);
    if (!bands) return {};
    double clearance = kInfinity;
    for (const auto& band : *bands) clearance = std::min(clearance, band.clearance);
    if (clearance == kInfinity) return std::array<double, 3>{};
    const auto latitude = EncloseLatitudeRootMonotonicity(mass, spin, radius, scale,
                                                          whole->initial_mu, clearance, *phase);
    const auto colour = ColourEnvelope(Range(1.0) / (Range(0.0) - whole->energy),
                                       Range(maximum_gaussian), brightness);
    if (!latitude || !colour) return {};

    Range upper_bound;
    std::array<unsigned, 3> cutoffs{};
    for (std::size_t i = 0; i < bands->size(); ++i) {
        const auto& band = (*bands)[i];
        if (band.clearance == kInfinity) continue;
        if (finite_cutoff) {
            cutoffs[i] = *finite_cutoff;
            continue;
        }
        const auto tail =
            EncloseWindingTail(mass, spin, radius, scale, whole->initial_mu, band.clearance,
                               determinant, *phase, *latitude, first_tail_winding);
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
    const auto phases = PrepareNearCriticalPhases(mass, spin, radius, scale, *whole);
    const auto screen = ScreenFiniteWindings(
        mass, spin, radius, scale, *whole, *polar, *phase, *angular, cutoffs, amplitude, clearance,
        remaining, cancelled, phases ? &*phases : nullptr, excluded_turns,
        finite_cutoff ? std::optional(finite_scale) : std::nullopt, *colour / determinant,
        PointDetectorPolicy{}.absolute_rgb_error * .01);
    if (!screen) return {};
    upper_bound = upper_bound + screen->grouped_rgb;
    for (const auto& candidate : screen->candidates) {
        const auto density =
            BoundAzimuthImages(mass, spin, radius, scale, *whole, *screen, *polar, candidate,
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
