#pragma once

#include "sirius/core/kerr_detector/angular_catalogue.h"
#include "sirius/core/kerr_detector/finite_phase.h"

#include <functional>
#include <span>
#include <vector>

namespace sirius::core::relativity::finite_tail_detail {

inline Range AsinPoint(double y) {
    if (!(y >= 0.0 && y < 1.0)) return InvalidRange();
    if (y == 0.0) return Range(0.0);
    const Range t(y), complement = sqrt(Range(1.0) - Square(t));
    if (y < 0x1.6a09e667f3bcdp-1) return AtanUnit(t / complement);
    return kPiBounds * Range(.5) - AtanUnit(complement / t);
}

inline Range AsinBounds(Range t) {
    if (!t.Finite() || t.lo < 0.0 || t.hi >= 1.0) return InvalidRange();
    return Range(AsinPoint(t.lo).lo, AsinPoint(t.hi).hi);
}

inline Range CosPoint(double angle) {
    Range term(1.0), sum(1.0);
    const Range square = Square(Range(angle));
    for (unsigned n = 1; n <= 13; ++n) {
        term = term * square / Range((2 * n - 1) * (2 * n));
        sum = n % 2 == 0 ? sum + term : sum - term;
    }
    const Range omitted = term * square / Range(27 * 28);
    return Range(sum.lo, (sum + omitted).hi);
}

inline Range PolarSineSquared(Range angle) { return Square(SineLobe(angle)); }

template <unsigned Order>
TaylorSeries<Order> PolarSineSquared(const TaylorSeries<Order>& angle) {
    // Here angle is the independent integration variable.
    TaylorSeries<Order> sine, cosine;
    const Range x = angle.coefficient[0];
    sine.coefficient[0] = SineLobe(x);
    cosine.coefficient[0] = Range(CosPoint(x.hi).lo, CosPoint(x.lo).hi);
    for (unsigned n = 1; n <= Order; ++n) {
        sine.coefficient[n] = cosine.coefficient[n - 1] / Range(n);
        cosine.coefficient[n] = (Range(0.0) - sine.coefficient[n - 1]) / Range(n);
    }
    return sine * sine;
}

struct FinitePolarShape {
    Range omega2, modulation, amplitude2, angular;
};

inline std::optional<FinitePolarShape> MakeFinitePolarShape(double mass, double spin, Range u,
                                                            Range scale) {
    const auto constants = InvertStationaryScalars(
        mass, spin,
        std::array{PupilAffine::Coordinate(u, true),
                   PupilAffine::Coordinate(Range(0.0) - Square(scale), false)});
    if (!constants) return {};
    const auto a2 = Square(PupilAffine(spin));
    const auto b = a2 - constants->eta - Square(constants->lambda);
    const auto omega2 =
        (sqrt(Square(b) + PupilAffine(4.0) * a2 * constants->eta) - b) * PupilAffine(.5);
    const auto amp2 = constants->eta / omega2;
    const auto modulation = a2 * amp2;
    const FinitePolarShape result{omega2.Bounds(), modulation.Bounds(), amp2.Bounds(),
                                  constants->lambda.Bounds()};
    if (!(result.omega2.lo > 0.0 && result.modulation.lo >= 0.0 && result.amplitude2.lo > 0.0 &&
          result.amplitude2.hi < 1.0 && result.angular.lo > 0.0))
        return {};
    return result;
}

struct PolarIntegralTable {
    static constexpr unsigned count = 1024;
    static constexpr double end = kPiBounds.hi * .5;
    FinitePolarShape shape;
    std::array<double, count + 1> nodes;
    std::array<PhasePair, count + 1> cumulative;
    PhasePair interpolation_error;

    static double node(unsigned i) { return std::scalbn(double(i) * end, -10); }

    PhasePair At(double angle) const {
        const auto upper = std::upper_bound(nodes.begin(), nodes.end(), angle);
        const auto interval =
            std::min(std::size_t(upper - nodes.begin() - 1), std::size_t(count - 1));
        const Range ratio = (Range(angle) - Range(nodes[interval])) /
                            (Range(nodes[interval + 1]) - Range(nodes[interval]));
        PhasePair answer;
        for (unsigned j = 0; j < 2; ++j) {
            answer[j] = (Range(1.0) - ratio) * cumulative[interval][j] +
                        ratio * cumulative[interval + 1][j] + interpolation_error[j];
        }
        return answer;
    }

    PhasePair Partial(Range angle) const {
        if (!angle.Finite() || angle.lo < 0.0 || angle.hi > end)
            return {InvalidRange(), InvalidRange()};
        const auto lo = At(angle.lo), hi = At(angle.hi);
        return {Range(lo[0].lo, hi[0].hi), Range(lo[1].lo, hi[1].hi)};
    }
};

inline std::optional<PolarIntegralTable> MakePolarIntegralTable(FinitePolarShape shape) {
    PolarIntegralTable table{};
    table.shape = shape;
    const auto in_sine2 = [&](auto s2) {
        using N = decltype(s2);
        const N time = N(1.0) / sqrt(N(shape.omega2) + N(shape.modulation) * s2);
        return std::array{time, time * N(shape.angular) / (N(1.0) - N(shape.amplitude2) * s2)};
    };
    const auto integrand = [&](auto angle) { return in_sine2(PolarSineSquared(angle)); };
    for (unsigned i = 0; i <= table.count; ++i) table.nodes[i] = table.node(i);
    double widest = 0.0;
    for (unsigned i = 0; i < table.count; ++i) {
        const Range lo(table.node(i)), hi(table.node(i + 1));
        widest = std::max(widest, (hi - lo).hi);
        const auto sample = PhaseGaussPanel(lo, hi, integrand);
        for (unsigned j = 0; j < 2; ++j) {
            table.cumulative[i + 1][j] = table.cumulative[i][j] + sample[j];
            if (!table.cumulative[i + 1][j].Finite()) return {};
        }
    }
    const auto derivative = in_sine2(TaylorSeries<1>::Variable(0.0, 1.0));
    for (unsigned j = 0; j < 2; ++j) {
        // The second derivative of the antiderivative is the derivative
        // of h(sin² psi). Its absolute chain factor |sin(2psi)| <= 1.
        const double magnitude =
            (Square(Range(widest)) * Abs(derivative[j].coefficient[1]) / Range(8.0)).hi;
        table.interpolation_error[j] = Range(-magnitude, magnitude);
    }
    return table;
}

inline Range ExpPoint(double value) {
    if (!std::isfinite(value) || std::abs(value) > 700.0) return InvalidRange();
    const int n = int(std::floor(value / 0x1.62e42fefa39efp-1));
    const Range r = Range(value) - Range(n) * LogPoint(2.0);
    if (Magnitude(r) > .75) return InvalidRange();
    Range term(1.0), sum(1.0);
    for (unsigned i = 1; i <= 24; ++i) {
        term = term * r / Range(i);
        sum = sum + term;
    }
    // Taylor remainder using exp(x) < 3 for |x| <= 3/4.
    const double error = (Range(3.0) * Abs(term * r) / Range(25.0)).hi;
    const Range reduced = sum + Range(-error, error);
    return Widen(std::scalbn(reduced.lo, n), std::scalbn(reduced.hi, n));
}

inline Range ExpBounds(Range value) {
    if (!value.Finite()) return InvalidRange();
    return Range(ExpPoint(value.lo).lo, ExpPoint(value.hi).hi);
}

struct FiniteCandidate {
    unsigned winding;
    const AngularCatalogue::Entry* source;
    bool approaching_equator;
    Range log_scale;
};

inline bool AzimuthMayMatch(Range unwrapped, Range sky) {
    if (!unwrapped.Finite() || !sky.Finite()) return true;
    const Range difference = (unwrapped - sky) / (Range(2.0) * kPiBounds);
    return std::ceil(difference.lo) <= std::floor(difference.hi);
}

inline Range FixedPupilAzimuth(double mass, double spin, const StationaryChart& chart) {
    const Range m(mass), a(spin), distance = sqrt(Square(m) - Square(a));
    const Range high = m + distance, low = m - distance;
    const Range u = chart.inverse_start;
    if (!(chart.nx.lo > 0.0)) return InvalidRange();
    return AtanUnit(chart.ny / chart.nx) +
           a / (Range(2.0) * distance) *
               LogBounds((Range(1.0) - low * u) / (Range(1.0) - high * u));
}

inline PhasePair IntegratePolarPartial(const FinitePolarShape& polar, Range phase) {
    const auto integrand = [&](auto angle) {
        using N = decltype(angle);
        const N s2 = PolarSineSquared(angle);
        const N time = N(1.0) / sqrt(N(polar.omega2) + N(polar.modulation) * s2);
        return std::array{time, time * N(polar.angular) / (N(1.0) - N(polar.amplitude2) * s2)};
    };
    if (!phase.Finite() || phase.lo < 0.0 || phase.hi > kPiBounds.hi) {
        return {InvalidRange(), InvalidRange()};
    }
    PhasePair result{};
    constexpr unsigned panels = 24;
    const Range step = Range(phase.lo) / Range(panels);
    if (phase.lo > 0.0) {
        for (unsigned i = 0; i < panels; ++i) {
            const auto evaluated = PhaseGaussPanel(Range(i) * step, Range(i + 1) * step, integrand);
            for (unsigned axis = 0; axis < 2; ++axis) result[axis] = result[axis] + evaluated[axis];
        }
    }
    const auto extra = integrand(phase);
    const auto remaining = Range(phase.hi) - Range(phase.lo);
    for (unsigned axis = 0; axis < 2; ++axis)
        result[axis] = result[axis] + Range(0.0, std::max(0.0, remaining.hi)) * extra[axis];
    return result;
}

inline PhasePair PolarTargetAt(double mass, double spin, Range centre, Range span, Range launch_mu,
                               Range source_mu, unsigned winding, bool approach) {
    const auto shape = MakeFinitePolarShape(mass, spin, centre, span);
    if (!shape) return {InvalidRange(), InvalidRange()};
    const auto amplitude = sqrt(shape->amplitude2);
    const auto half = IntegratePolarPartial(*shape, Range(kPiBounds.lo * .5, kPiBounds.hi * .5));
    const auto initial = IntegratePolarPartial(*shape, AsinBounds(launch_mu / amplitude));
    const auto final = IntegratePolarPartial(*shape, AsinBounds(Abs(source_mu) / amplitude));
    PhasePair target;
    for (unsigned axis = 0; axis < 2; ++axis) {
        const Range period = Range(2.0) * half[axis];
        target[axis] = Range(winding) * period + initial[axis] +
                       (approach ? period - final[axis] : final[axis]);
    }
    return target;
}

struct RefinedCandidate {
    bool possible;
    Range log_scale, azimuth;
};

// High windings reach the same near-critical strip of the radial chart.
// Share enclosures of the radial integrals and polar antiderivatives over
// the strip instead of reintegrating them for each catalogue candidate.
struct NearCriticalPhases {
    FiniteRadialPhase radial;
    PolarIntegralTable polar;
    PhasePair period, initial;
    Range absolute_amplitude, fixed_azimuth;
    double maximum_log_scale;

    std::optional<RefinedCandidate> Refine(const FiniteCandidate& candidate) const {
        const auto& star = *candidate.source;
        const auto final = polar.Partial(AsinBounds(star.latitude / absolute_amplitude));
        PhasePair target;
        for (unsigned axis = 0; axis < 2; ++axis) {
            const Range endpoint =
                candidate.approaching_equator ? period[axis] - final[axis] : final[axis];
            target[axis] = Range(candidate.winding) * period[axis] + initial[axis] + endpoint;
        }
        const Range azimuth =
            fixed_azimuth - (target[1] + radial.drag_at_turn * target[0] + radial.spin_residual);
        const Range timing = (radial.time_residual - target[0]) / (Range(2.0) * radial.leading);
        if (!azimuth.Finite() || !timing.Finite()) return {};
        const Range log_scale(std::max(timing.lo, candidate.log_scale.lo),
                              std::min(timing.hi, candidate.log_scale.hi));
        return RefinedCandidate{log_scale.Finite() && AzimuthMayMatch(azimuth, star.azimuth),
                                log_scale, azimuth};
    }
};

inline std::optional<NearCriticalPhases> PrepareNearCriticalPhases(double mass, double spin,
                                                                   Range centre, Range scale,
                                                                   const StationaryChart& chart) {
    // The strip lies inside the admitted radial domain. This cutoff selects
    // a shared integration domain and does not discard any source.
    const double edge = std::min(0x1p-40, scale.hi);
    const Range limit(0.0, edge);
    const auto radial = EvaluateFiniteRadialPhase(mass, spin, centre, limit, chart.inverse_start);
    const auto shape = MakeFinitePolarShape(mass, spin, centre, limit);
    const Range azimuth = FixedPupilAzimuth(mass, spin, chart);
    const Range log_edge = LogPoint(edge);
    if (!radial || !shape || !azimuth.Finite() || !log_edge.Finite()) return {};
    const Range amplitude = sqrt(shape->amplitude2);
    const auto half = IntegratePolarPartial(*shape, Range(kPiBounds.lo * .5, kPiBounds.hi * .5));
    const auto initial = IntegratePolarPartial(*shape, AsinBounds(chart.initial_mu / amplitude));
    auto table = MakePolarIntegralTable(*shape);
    if (!table || !half[0].Finite() || !half[1].Finite() || !initial[0].Finite() ||
        !initial[1].Finite())
        return {};
    return NearCriticalPhases{
        *radial,    std::move(*table), {Range(2.0) * half[0], Range(2.0) * half[1]},
        initial,    amplitude,         azimuth,
        log_edge.lo};
}

// Share the differentiated fixed-angle integrals across a strip of small
// radial turning separations. Each partial then needs only one local panel
// and the Leibniz boundary term for its moving polar endpoint.
class PolarDerivativeStrip {
  public:
    static std::optional<PolarDerivativeStrip> Prepare(double mass, double spin, Range centre,
                                                       Range scale, const StationaryChart& chart) {
        PolarDerivativeStrip result;
        const Range domain(0.0, std::min(scale.hi, 0x1p-10));
        const auto radial =
            EncloseRegularizedRadialPhase(mass, spin, centre, domain, chart.inverse_start, 128);
        const FilmJet u(centre, {Range(1.0), Range(0.0)}), s(domain, {Range(0.0), Range(1.0)});
        const auto constants =
            InvertStationaryScalars(mass, spin, std::array{u, FilmJet(0.0) - Square(s)});
        if (!radial || !constants) return {};
        result.radial_ = *radial;
        result.a2_ = Square(FilmJet(spin));
        result.angular_ = constants->lambda;
        const auto p = result.a2_ - constants->eta - Square(result.angular_);
        result.omega2_ =
            (sqrt(Square(p) + FilmJet(4.0) * result.a2_ * constants->eta) - p) * FilmJet(.5);
        result.amplitude2_ = constants->eta / result.omega2_;
        result.amplitude_ = sqrt(result.amplitude2_);
        result.maximum_log_scale_ = LogPoint(domain.hi).lo;
        if (!result.amplitude_.value.Finite() || !std::isfinite(result.maximum_log_scale_))
            return {};
        for (unsigned i = 0; i < kPanels; ++i) {
            const auto sample = result.At({Node(i), Node(i + 1)});
            const FilmJet width(Range(Node(i + 1)) - Range(Node(i)));
            for (unsigned j = 0; j < sample.size(); ++j)
                result.prefix_[i + 1][j] = result.prefix_[i][j] + sample[j] * width;
        }
        const auto half = result.Integrate(FilmJet(Range(kPiBounds.lo * .5, kPiBounds.hi * .5)));
        const auto initial = result.ToLatitude(chart.initial_mu);
        if (!half || !initial) return {};
        result.initial_ = *initial;
        for (unsigned j = 0; j < half->size(); ++j) result.period_[j] = (*half)[j] * FilmJet(2.0);
        return result;
    }

    using PolarJet = std::array<FilmJet, 2>;

    std::optional<PolarJet> ToLatitude(Range latitude) const {
        const FilmJet ratio = FilmJet(latitude) / amplitude_;
        const Range value = AsinBounds(ratio.value),
                    divisor = sqrt(Range(1.0) - Square(ratio.value));
        return Integrate({value, {ratio.derivative[0] / divisor, ratio.derivative[1] / divisor}});
    }

    bool Covers(Range log_separation) const { return log_separation.hi <= maximum_log_scale_; }

    std::optional<Range> Density(Range latitude, const PolarJet& final, unsigned winding,
                                 bool approach, Range log_separation, Range separation,
                                 Range sky) const {
        if (!sky.Finite()) return {};
        PolarJet target;
        for (unsigned j = 0; j < target.size(); ++j)
            target[j] = period_[j] * FilmJet(winding) + initial_[j] +
                        (approach ? period_[j] - final[j] : final[j]);
        const Range radial_u = Range(-2.0) * radial_.leading_gradient[0] * log_separation +
                               radial_.time_remainder.radial;
        const Range transverse = separation * target[0].derivative[1] - radial_.time_log_derivative;
        if (!(transverse.lo > 0.0)) return {};
        const Range log_s_u = (radial_u - target[0].derivative[0]) / transverse;
        const Range slope = Abs(
            target[1].derivative[0] + radial_.spin_remainder.radial +
            (separation * target[1].derivative[1] + radial_.spin_remainder.log_scale) * log_s_u);
        if (!(slope.lo > 0.0)) return {};
        const Range z2 = Square(latitude);
        const Range speed = sqrt((amplitude2_.value - z2) * (omega2_.value + a2_.value * z2));
        // Each fixed latitude has one monotone azimuth branch across the
        // admitted chart. One source can occur at most once per lift in sky.
        const Range span = (Range(sky.hi) - Range(sky.lo)) / (Range(2.0) * kPiBounds);
        const Range roots = Range(1.0) + Range(std::floor(span.hi));
        const Range density =
            roots * Range(2.0) * Square(separation) / (slope * transverse * speed);
        return density.Finite() ? std::optional(density) : std::nullopt;
    }

  private:
    static constexpr unsigned kPanels = 128;
    static constexpr double kEnd = kPiBounds.hi * .5;
    static double Node(unsigned i) { return std::scalbn(double(i) * kEnd, -7); }

    PolarJet At(Range angle) const {
        const auto sin2 = FilmJet(Square(SineLobe(angle)));
        const FilmJet time = FilmJet(1.0) / sqrt(omega2_ + a2_ * amplitude2_ * sin2);
        return {time, time * (angular_ / (FilmJet(1.0) - amplitude2_ * sin2) +
                              FilmJet(radial_.limiting_drag, radial_.drag_gradient))};
    }

    std::optional<PolarJet> Integrate(const FilmJet& phase) const {
        const Range angle = phase.value;
        if (!angle.Finite() || angle.lo < 0.0 || angle.hi > kEnd) return {};
        unsigned i = std::min(kPanels, unsigned(angle.lo / kEnd * kPanels));
        if (Node(i) > angle.lo) --i;
        const auto sample = At({Node(i), angle.hi});
        const auto boundary = At(angle);
        const FilmJet length(angle - Range(Node(i)));
        PolarJet result;
        for (unsigned j = 0; j < sample.size(); ++j) {
            result[j] = prefix_[i][j] + sample[j] * length;
            for (unsigned axis = 0; axis < 2; ++axis)
                result[j].derivative[axis] =
                    result[j].derivative[axis] + boundary[j].value * phase.derivative[axis];
            if (!result[j].value.Finite() || !result[j].derivative[0].Finite() ||
                !result[j].derivative[1].Finite())
                return {};
        }
        return result;
    }

    FilmJet a2_, angular_, omega2_, amplitude2_, amplitude_;
    RegularizedRadialPhase radial_;
    double maximum_log_scale_;
    PolarJet period_, initial_;
    std::array<PolarJet, kPanels + 1> prefix_{};
};

struct ReachableCatalogueRow {
    const AngularCatalogue::Row* row;
    Range latitude;
    unsigned band;
    std::vector<const AngularCatalogue::Entry*> subset;
    Range weighted_temperature;
};

using ReachableCatalogueRows = std::array<std::vector<ReachableCatalogueRow>, 2>;

struct CatalogueScreen {
    std::vector<FiniteCandidate> candidates;
    Range grouped_rgb;
    Range edge_time;
    PhasePair period, initial;
};

inline std::optional<CatalogueScreen> ScreenFiniteWindings(
    double mass, double spin, Range centre, Range scale, const StationaryChart& chart,
    const PolarIntegralTable& polar, const RegularizedRadialPhase& prior,
    const ReachableCatalogueRows& reachable, const std::array<unsigned, 3>& cutoffs,
    std::size_t& remaining, const std::function<bool()>& cancelled,
    const NearCriticalPhases* shared_phases, std::optional<Range> finite_scale,
    Range density_to_rgb, double grouped_allowance) {
    const Range fixed_phi = FixedPupilAzimuth(mass, spin, chart);
    const Range amplitude = sqrt(polar.shape.amplitude2);
    if (!fixed_phi.Finite() ||
        !(chart.energy.hi < 0.0 && chart.initial_mu.lo > 0.0 && chart.mu_rate.hi < 0.0))
        return {};
    const auto half = polar.Partial(Range(kPiBounds.lo * .5, kPiBounds.hi * .5));
    const PhasePair period{half[0] * Range(2.0), half[1] * Range(2.0)};
    const auto initial = polar.Partial(AsinBounds(chart.initial_mu / amplitude));
    const Range largest_log = LogPoint(scale.hi);
    const auto edge =
        EvaluateFiniteRadialPhase(mass, spin, centre, Range(scale.hi), chart.inverse_start);
    const auto critical =
        EvaluateFiniteRadialPhase(mass, spin, centre, Range(0.0), chart.inverse_start);
    if (!edge || !critical || !initial[0].Finite() || !period[0].Finite()) return {};
    const Range edge_time = Range(-2.0) * edge->leading * largest_log + edge->time_residual;
    std::optional<FiniteRadialPhase> physical_radial;
    double earliest = edge_time.lo;
    double latest = kInfinity;
    if (finite_scale) {
        physical_radial =
            EvaluateFiniteRadialPhase(mass, spin, centre, *finite_scale, chart.inverse_start);
        if (physical_radial) {
            const Range time = Range(-2.0) * physical_radial->leading * LogBounds(*finite_scale) +
                               physical_radial->time_residual;
            if (time.Finite()) {
                earliest = std::max(earliest, time.lo);
                latest = time.hi;
            }
        }
    }
    // Screen against the physical escape interval. Root counting below
    // still needs the complete latitude branch over the wider radial chart.
    const auto can_arrive = [&](Range time) { return !(time.hi < earliest || time.lo > latest); };
    const auto sky_window = [&](Range angle, Range time, const FiniteRadialPhase& radial) {
        const Range arrival = std::isfinite(latest) && time.Finite()
                                  ? Range(std::max(time.lo, earliest), std::min(time.hi, latest))
                                  : time;
        const Range original =
            fixed_phi - (angle + radial.drag_at_turn * arrival + radial.spin_residual);
        if (!physical_radial) return original;
        const Range physical = fixed_phi - (angle + physical_radial->drag_at_turn * arrival +
                                            physical_radial->spin_residual);
        if (!physical.Finite()) return original;
        if (!original.Finite()) return physical;
        return Range(std::max(original.lo, physical.lo), std::min(original.hi, physical.hi));
    };
    const Range twice_leading = Range(2.0) * prior.leading_factor;
    const auto radial_bound = [&](Range log_scale) -> std::optional<FiniteRadialPhase> {
        const Range span = ExpBounds(log_scale);
        if (!span.Finite() || span.lo > scale.hi) return {};
        if (span.hi > 1e-6) {
            return EvaluateFiniteRadialPhase(mass, spin, centre,
                                             Range(span.lo, std::min(span.hi, scale.hi)),
                                             chart.inverse_start);
        }
        const Range decay = Range(span.hi) * (Range(prior.spin_derivative_constant) +
                                              Range(prior.spin_derivative_log_coefficient) *
                                                  (Range(1.0) - LogPoint(span.hi)));
        auto phase = *critical;
        phase.drag_at_turn = phase.drag_at_turn + Range(0.0, span.hi) * prior.drag_gradient[1];
        phase.spin_residual = phase.spin_residual + Range(-decay.hi, decay.hi);
        return phase;
    };
    const unsigned last = *std::ranges::max_element(cutoffs);
    // At small radial separations an entire row can contribute very
    // little even when its longitude window contains many sources.
    constexpr unsigned first_grouped_winding = 8;
    constexpr unsigned first_whole_row_winding = 12;
    const auto strip = finite_scale || last <= first_grouped_winding
                           ? std::nullopt
                           : PolarDerivativeStrip::Prepare(mass, spin, centre, scale, chart);
    struct ScreenRow {
        const ReachableCatalogueRow* sources;
        std::optional<PolarDerivativeStrip::PolarJet> endpoint;
        PhasePair equator, critical_equator;
    };
    std::array<std::vector<ScreenRow>, 2> screen_rows;
    for (unsigned hemisphere = 0; hemisphere < reachable.size(); ++hemisphere) {
        auto& output = screen_rows[hemisphere];
        output.reserve(reachable[hemisphere].size());
        for (const auto& row : reachable[hemisphere]) {
            if (cancelled && cancelled()) return {};
            const auto equator = polar.Partial(AsinBounds(row.latitude / amplitude));
            if (!equator[0].Finite() || !equator[1].Finite()) return {};
            PhasePair critical_equator = equator;
            if (shared_phases) {
                const auto partial = shared_phases->polar.Partial(
                    AsinBounds(row.latitude / shared_phases->absolute_amplitude));
                if (partial[0].Finite() && partial[1].Finite()) critical_equator = partial;
            }
            const auto endpoint = strip ? strip->ToLatitude(row.latitude) : std::nullopt;
            output.push_back({&row, endpoint, equator, critical_equator});
        }
    }
    CatalogueScreen result{{}, {}, edge_time, period, initial};
    for (unsigned m = 0; m < last; ++m) {
        if (cancelled && cancelled()) return {};
        const Range m_value(m);
        const Range full = m_value * period[0] + initial[0] + Range(0.0, period[0].hi);
        if (!can_arrive(full)) continue;
        const Range log_span = (prior.time_remainder.value - full) / twice_leading;
        const Range log_branch(log_span.lo, std::min(largest_log.hi, log_span.hi));
        if (!log_span.Finite()) return {};
        if (log_branch.lo > log_branch.hi) continue;
        const auto grouped_separation =
            strip && strip->Covers(log_branch) ? ExpBounds(log_branch) : InvalidRange();
        const bool near_critical =
            shared_phases && log_branch.hi <= shared_phases->maximum_log_scale;
        const auto coarse =
            near_critical ? std::optional(shared_phases->radial) : radial_bound(log_branch);
        if (!coarse) return {};
        const bool full_radial = near_critical || ExpPoint(log_branch.hi).lo > 1e-6;
        const auto& active_period = near_critical ? shared_phases->period : period;
        const auto& active_initial = near_critical ? shared_phases->initial : initial;
        const auto EquatorToEnd = [&](PhasePair equator, bool approach) {
            return approach
                       ? PhasePair{active_period[0] - equator[0], active_period[1] - equator[1]}
                       : equator;
        };
        const auto inspect = [&](const AngularCatalogue::Entry& star, bool approach) -> bool {
            auto equator = near_critical ? shared_phases->polar.Partial(AsinBounds(
                                               star.latitude / shared_phases->absolute_amplitude))
                                         : polar.Partial(AsinBounds(star.latitude / amplitude));
            if (near_critical && (!equator[0].Finite() || !equator[1].Finite()))
                equator = polar.Partial(AsinBounds(star.latitude / amplitude));
            if (!equator[0].Finite() || !equator[1].Finite()) return false;
            const auto part = EquatorToEnd(equator, approach);
            const Range tau = m_value * active_period[0] + active_initial[0] + part[0];
            const Range phi = m_value * active_period[1] + active_initial[1] + part[1];
            if (!can_arrive(tau)) return true;
            const Range sky = sky_window(phi, tau, *coarse);
            if (sky.lo > sky.hi || !AzimuthMayMatch(sky, star.azimuth)) return true;

            const Range local_log = (coarse->time_residual - tau) / (Range(2.0) * coarse->leading);
            // Without the shared radial integrals the critical-limit
            // estimate adjusts only the spin terms. Retain the full
            // logarithmic timing envelope in that case.
            const Range preliminary = (prior.time_remainder.value - tau) / twice_leading;
            Range valid(std::max(log_branch.lo, preliminary.lo),
                        std::min(log_branch.hi, preliminary.hi));
            if (!preliminary.Finite() || (full_radial && !local_log.Finite())) return false;
            if (full_radial) {
                valid = Range(std::max(valid.lo, local_log.lo), std::min(valid.hi, local_log.hi));
            }
            if (valid.lo > valid.hi) return true;
            if (!valid.Finite()) return false;
            if (!near_critical) {
                const auto refined = radial_bound(valid);
                if (!refined) return false;
                const Range tighter =
                    fixed_phi - (phi + refined->drag_at_turn * tau + refined->spin_residual);
                if (!AzimuthMayMatch(tighter, star.azimuth)) return true;
            }
            result.candidates.push_back({m, &star, approach, valid});
            return true;
        };
        // The represented polar motion begins north of the equator and
        // toward it. Each complete half-libration switches hemispheres.
        for (const auto& prepared : screen_rows[m % 2]) {
            const auto& sources = *prepared.sources;
            if (m >= cutoffs[sources.band]) continue;
            for (const bool approach : {false, true}) {
                const auto part = EquatorToEnd(
                    near_critical ? prepared.critical_equator : prepared.equator, approach);
                const Range tau = m_value * active_period[0] + active_initial[0] + part[0];
                if (!can_arrive(tau)) continue;
                const Range phi = m_value * active_period[1] + active_initial[1] + part[1];
                const Range sky = sky_window(phi, tau, *coarse);
                if (sky.lo > sky.hi) continue;
                if (m >= first_grouped_winding && grouped_separation.Finite() &&
                    prepared.endpoint && sky.Finite()) {
                    const Range full_part =
                        approach ? period[0] - prepared.equator[0] : prepared.equator[0];
                    const Range full_target = m_value * period[0] + initial[0] + full_part;
                    if (full_target.Finite() && full_target.lo > edge_time.hi) {
                        // Reduce the matching longitude slices via prefix
                        // sums, or use the row's outward aggregate for a
                        // sufficiently faint winding.
                        const auto density =
                            strip->Density(sources.latitude, *prepared.endpoint, m, approach,
                                           log_branch, grouped_separation, sky);
                        if (density) {
                            const bool whole_row =
                                m >= first_whole_row_winding || !sources.subset.empty();
                            const auto weight = whole_row
                                                    ? std::optional(sources.weighted_temperature)
                                                    : sources.row->QueryWeight(sky);
                            const Range contribution =
                                weight ? density_to_rgb * *weight * *density : InvalidRange();
                            const Range total = result.grouped_rgb + Range(0.0, contribution.hi);
                            // Keep most of the grouped allowance available to
                            // deeper windings even when early windows match.
                            const double allowance = m >= first_whole_row_winding
                                                         ? grouped_allowance
                                                         : grouped_allowance * .25;
                            if (contribution.Finite() && contribution.hi >= 0.0 &&
                                total.hi <= allowance) {
                                result.grouped_rgb = total;
                                continue;
                            }
                        }
                    }
                }
                const auto visit = [&](const AngularCatalogue::Entry& entry) {
                    if (remaining == 0 || (cancelled && cancelled())) return false;
                    --remaining;
                    return inspect(entry, approach);
                };
                const auto dereference =
                    [](const AngularCatalogue::Entry* entry) -> const AngularCatalogue::Entry& {
                    return *entry;
                };
                const bool completed = sources.subset.empty()
                                           ? sources.row->Query(sky, visit)
                                           : AngularCatalogue::Row::QueryEntries(
                                                 sources.subset, sources.row->longitude_error, sky,
                                                 dereference, visit);
                if (!completed) return {};
            }
        }
    }
    return result;
}

inline std::optional<RefinedCandidate> RefineFiniteAzimuth(double mass, double spin, Range centre,
                                                           Range outer_scale,
                                                           const StationaryChart& chart,
                                                           const FiniteCandidate& candidate,
                                                           Range source_mu, Range star_azimuth) {
    Range log_scale = candidate.log_scale;
    const Range chart_azimuth = FixedPupilAzimuth(mass, spin, chart);
    Range azimuth;
    for (unsigned i = 0; i < 5; ++i) {
        const Range expanded = ExpBounds(log_scale);
        if (!expanded.Finite()) return {};
        if (expanded.lo > outer_scale.hi) return RefinedCandidate{};
        const Range s(expanded.lo, std::min(expanded.hi, outer_scale.hi));
        const auto radial = EvaluateFiniteRadialPhase(mass, spin, centre, s, chart.inverse_start);
        const auto target = PolarTargetAt(mass, spin, centre, s, chart.initial_mu, source_mu,
                                          candidate.winding, candidate.approaching_equator);
        if (!radial || !target[0].Finite() || !target[1].Finite()) return {};
        azimuth =
            chart_azimuth - (target[1] + radial->drag_at_turn * target[0] + radial->spin_residual);
        if (!azimuth.Finite()) return {};
        if (!AzimuthMayMatch(azimuth, star_azimuth))
            return RefinedCandidate{false, log_scale, azimuth};
        const Range next = (radial->time_residual - target[0]) / (Range(2.0) * radial->leading);
        if (!next.Finite()) return {};
        const Range intersect(std::max(log_scale.lo, next.lo), std::min(log_scale.hi, next.hi));
        if (!intersect.Finite()) return RefinedCandidate{};
        if ((intersect.hi - intersect.lo) >= .8 * (log_scale.hi - log_scale.lo)) break;
        log_scale = intersect;
    }
    return RefinedCandidate{true, log_scale, azimuth};
}

inline bool EntireLatitudeBranchExists(const CatalogueScreen& screen,
                                       const PolarIntegralTable& polar,
                                       const FiniteCandidate& candidate) {
    const Range amplitude = sqrt(polar.shape.amplitude2);
    const Range at_source = polar.Partial(AsinBounds(candidate.source->latitude / amplitude))[0];
    const Range part = candidate.approaching_equator ? screen.period[0] - at_source : at_source;
    const Range target = Range(candidate.winding) * screen.period[0] + screen.initial[0] + part;
    return target.Finite() && target.lo > screen.edge_time.hi;
}

struct AzimuthDerivative {
    Range conditioned_radial_slope;
    Range latitude_log_slope_magnitude;
    Range mu_speed;
};

// Integrate the fixed-angle parameter derivatives with the same enclosed
// Gauss rule used for the phase values. Leibniz's boundary term carries the
// changing polar endpoint.
inline std::array<FilmJet, 2> IntegratePolarDerivatives(const FilmJet& phase, const FilmJet& omega2,
                                                        const FilmJet& modulation,
                                                        const FilmJet& amplitude2,
                                                        const FilmJet& angular,
                                                        const FilmJet& drag) {
    std::array<FilmJet, 2> result{};
    if (!phase.value.Finite() || phase.value.lo < 0.0 || phase.value.hi > kPiBounds.hi * .5) {
        for (auto& term : result) term.derivative.fill(InvalidRange());
        return result;
    }
    const auto density = [&](auto sin2, auto frequency, auto variation, auto amplitude,
                             auto momentum, auto radial_drag) {
        using N = decltype(sin2);
        const N time = N(1.0) / sqrt(frequency + variation * sin2);
        return std::array{time, time * (momentum / (N(1.0) - amplitude * sin2) + radial_drag)};
    };
    const auto boundary = density(PolarSineSquared(phase.value), omega2.value, modulation.value,
                                  amplitude2.value, angular.value, drag.value);
    constexpr unsigned panels = 32;
    const Range step = Range(phase.value.lo) / Range(panels);
    const double tail = (Range(phase.value.hi) - Range(phase.value.lo)).hi;
    for (unsigned axis = 0; axis < phase.derivative.size(); ++axis) {
        const auto differential = [&](auto angle) {
            using N = decltype(angle);
            using Jet = RadialJet<N>;
            const auto jet = [&](const FilmJet& p) {
                return Jet(N(p.value), N(p.derivative[axis]));
            };
            const auto q = density(Jet(PolarSineSquared(angle)), jet(omega2), jet(modulation),
                                   jet(amplitude2), jet(angular), jet(drag));
            return std::array{q[0].dual, q[1].dual};
        };
        PhasePair integral{};
        if (phase.value.lo > 0.0) {
            for (unsigned i = 0; i < panels; ++i) {
                const auto part =
                    PhaseGaussPanel(Range(i) * step, Range(i + 1) * step, differential);
                for (unsigned j = 0; j < part.size(); ++j) integral[j] = integral[j] + part[j];
            }
        }
        const auto end = differential(phase.value);
        for (unsigned j = 0; j < result.size(); ++j)
            result[j].derivative[axis] =
                integral[j] + Range(0.0, tail) * end[j] + boundary[j] * phase.derivative[axis];
    }
    return result;
}

inline std::optional<AzimuthDerivative> DifferentiateFiniteAzimuth(
    double mass, double spin, Range centre, Range s, Range log_s, const StationaryChart& chart,
    unsigned winding, bool approach, Range source_mu) {
    const auto radial =
        EncloseRegularizedRadialPhase(mass, spin, centre, s, chart.inverse_start, 128);
    const FilmJet u(centre, {Range(1.0), Range(0.0)}), scale(s, {Range(0.0), Range(1.0)});
    const auto constants =
        InvertStationaryScalars(mass, spin, std::array{u, FilmJet(0.0) - Square(scale)});
    if (!radial || !constants) return {};
    const FilmJet a2 = Square(FilmJet(spin));
    const FilmJet p = a2 - constants->eta - Square(constants->lambda);
    const FilmJet omega2 = (sqrt(Square(p) + FilmJet(4.0) * a2 * constants->eta) - p) * FilmJet(.5);
    const FilmJet amplitude2 = constants->eta / omega2;
    const FilmJet drag(radial->limiting_drag, radial->drag_gradient);
    const auto asin_mu = [&](Range mu) {
        const FilmJet ratio = FilmJet(Abs(mu)) / sqrt(amplitude2);
        const Range value = AsinBounds(ratio.value),
                    divisor = sqrt(Range(1.0) - Square(ratio.value));
        return FilmJet(value, {ratio.derivative[0] / divisor, ratio.derivative[1] / divisor});
    };
    const auto integrate = [&](const FilmJet& phase) {
        std::array<FilmJet, 2> result{};
        constexpr unsigned panels = 128;
        for (unsigned i = 0; i < panels; ++i) {
            // Differentiate the integral after changing its variable to
            // [0,1]; this carries the moving partial polar endpoint.
            const FilmJet angle = phase * FilmJet(Range(i, i + 1) / Range(panels));
            const Range v = SineLobe(angle.value);
            const Range w(CosPoint(angle.value.hi).lo, CosPoint(std::max(0.0, angle.value.lo)).hi);
            const FilmJet sin(v, {w * angle.derivative[0], w * angle.derivative[1]});
            const FilmJet sin2 = Square(sin);
            const FilmJet endpoint_time = phase / sqrt(omega2 + a2 * amplitude2 * sin2);
            const FilmJet endpoint_azimuth =
                endpoint_time * (constants->lambda / (FilmJet(1.0) - amplitude2 * sin2) + drag);
            result[0] = result[0] + endpoint_time / FilmJet(panels);
            result[1] = result[1] + endpoint_azimuth / FilmJet(panels);
        }
        return result;
    };
    const auto derivative = [&](const auto& quadrature) -> std::optional<AzimuthDerivative> {
        const auto half = quadrature(FilmJet(Range(kPiBounds.lo * .5, kPiBounds.hi * .5)));
        const auto initial = quadrature(asin_mu(chart.initial_mu));
        const auto final = quadrature(asin_mu(source_mu));
        std::array<FilmJet, 2> target;
        for (unsigned i = 0; i < 2; ++i) {
            const FilmJet period = half[i] * FilmJet(2.0);
            target[i] =
                period * FilmJet(winding) + initial[i] + (approach ? period - final[i] : final[i]);
        }
        const Range radial_u =
            Range(-2.0) * radial->leading_gradient[0] * log_s + radial->time_remainder.radial;
        const Range transverse = s * target[0].derivative[1] - radial->time_log_derivative;
        if (!(transverse.lo > 0.0)) return {};
        const Range log_s_u = (radial_u - target[0].derivative[0]) / transverse;
        const Range conditioned =
            Range(0.0) -
            (target[1].derivative[0] + radial->spin_remainder.radial +
             (s * target[1].derivative[1] + radial->spin_remainder.log_scale) * log_s_u);
        if (!conditioned.Finite() || !transverse.Finite()) return {};
        const Range z2 = Square(Abs(source_mu));
        const Range speed = sqrt((amplitude2.value - z2) * (omega2.value + a2.value * z2));
        if (!speed.Finite() || !(speed.lo > 0.0)) return {};
        return AzimuthDerivative{conditioned, transverse, speed};
    };
    const auto coarse = derivative(integrate);
    if (s.hi >= 0x1p-20 || (coarse && Abs(coarse->conditioned_radial_slope).lo > 0.0))
        return coarse;
    const auto enclosed = [&](const FilmJet& endpoint) {
        return IntegratePolarDerivatives(endpoint, omega2, a2 * amplitude2, amplitude2,
                                         constants->lambda, drag);
    };
    const auto sharp = derivative(enclosed);
    return sharp ? sharp : coarse;
}

}  // namespace sirius::core::relativity::finite_tail_detail
