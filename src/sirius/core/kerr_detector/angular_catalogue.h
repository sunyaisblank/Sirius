#pragma once

#include "sirius/core/kerr_detector/winding_tail.h"
#include "sirius/core/spectral/blackbody.h"
#include "sirius/core/starfield.h"

#include <functional>

namespace sirius::core::relativity::finite_tail_detail {

inline Range AtanSmall(Range x) {
    using finite_tail_detail::Abs;
    using finite_tail_detail::Square;
    const Range x2 = Square(x);
    Range term = x, sum = x;
    for (unsigned n = 1; n < 20; ++n) {
        term = term * x2;
        sum = n % 2 != 0 ? sum - term / Range(2 * n + 1) : sum + term / Range(2 * n + 1);
    }
    const double remainder = (Abs(term) * x2 / Range(41.0)).hi;
    return sum + Range(-remainder, remainder);
}

inline Range AtanPoint(double value) {
    const double magnitude = std::abs(value);
    const Range x(magnitude);
    const Range result =
        magnitude <= 0.5 ? AtanSmall(x)
                         : kPiBounds * Range(0.25) + AtanSmall((x - Range(1.0)) / (x + Range(1.0)));
    return value < 0.0 ? Range(-result.hi, -result.lo) : result;
}

inline Range AtanUnit(Range value) {
    if (!value.Finite() || finite_tail_detail::Magnitude(value) > 1.00000000000001)
        return finite_tail_detail::InvalidRange();
    return Range(AtanPoint(value.lo).lo, AtanPoint(value.hi).hi);
}

inline Range CatalogueAzimuth(double x, double y) {
    if (x == 0.0 && y == 0.0) return finite_tail_detail::InvalidRange();
    if (std::abs(x) >= std::abs(y)) {
        const auto theta = AtanUnit(Range(x < 0 ? -y : y) / Range(std::abs(x)));
        if (x < 0.0) return theta + (y >= 0.0 ? kPiBounds : Range(-kPiBounds.hi, -kPiBounds.lo));
        return theta;
    }
    const auto theta = AtanUnit(Range(x) / Range(std::abs(y)));
    return y > 0.0 ? kPiBounds * Range(0.5) - theta : theta - kPiBounds * Range(0.5);
}

// The discrete groups share one high-winding bound at their smallest
// certified distance from a polar turn.
inline unsigned LatitudeBand(double gap) { return gap >= .01 ? 0 : gap >= .0001 ? 1 : 2; }

// Catalogue coordinates are independent of the pupil and film geometry.
// The angular index can be built once and shared across winding queries.
class AngularCatalogue {
  public:
    static constexpr unsigned kLatitudeRows = 512;

    struct Entry {
        Range latitude, azimuth;
        Range weighted_temperature;
        double longitude_key;
    };

    struct Row {
        std::vector<Entry> entries;
        Range latitude, weighted_temperature;
        double longitude_error = 0.0;

        template <class Visit>
        bool Query(Range unwrapped, Visit&& visit) const {
            if (entries.empty()) return true;
            if (!unwrapped.Finite()) return false;
            const Range full(-kPiBounds.hi, kPiBounds.hi);
            const Range two_pi = Range(2.0) * kPiBounds;
            const Range wrap_range = (unwrapped - full) / two_pi;
            if (!wrap_range.Finite() || Magnitude(wrap_range) > 1000000.0) return false;
            const int first = int(std::ceil(wrap_range.lo)), last = int(std::floor(wrap_range.hi));
            // Multiple lifts that cover a full revolution need only one
            // visit per catalogue entry.
            if (last - first > 1 ||
                (Range(unwrapped.hi) - Range(unwrapped.lo)).hi + 2 * longitude_error >= two_pi.lo) {
                return std::all_of(entries.begin(), entries.end(), visit);
            }
            using Iterator = std::vector<Entry>::const_iterator;
            std::array<std::pair<Iterator, Iterator>, 2> slices{};
            unsigned count = 0;
            for (int winding = first; winding <= last; ++winding) {
                const Range arc =
                    unwrapped - two_pi * Range(winding) + Range(-longitude_error, longitude_error);
                const auto lo =
                    std::ranges::lower_bound(entries, arc.lo, {}, &Entry::longitude_key);
                const auto hi =
                    std::ranges::upper_bound(entries, arc.hi, {}, &Entry::longitude_key);
                if (lo != hi) slices[count++] = {lo, hi};
            }
            if (count == 2 && slices[1].first < slices[0].first) std::swap(slices[0], slices[1]);
            Iterator previous_end = entries.begin();
            for (unsigned i = 0; i < count; ++i) {
                const auto [lo, hi] = slices[i];
                for (auto at = std::max(lo, previous_end); at < hi; ++at)
                    if (!visit(*at)) return false;
                previous_end = std::max(previous_end, hi);
            }
            return true;
        }
    };

    static std::optional<AngularCatalogue> Build(std::span<const StarEntry> stars,
                                                 const std::function<bool()>& cancelled) {
        if (stars.size() > std::numeric_limits<std::uint32_t>::max()) return {};
        AngularCatalogue index;
        for (std::size_t i = 0; i < stars.size(); ++i) {
            if (i % 256 == 0 && cancelled && cancelled()) return {};
            const auto& star = stars[i];
            const auto colour = spectral::BlackbodyBandRgb(star.temperature_K);
            if (!colour || !std::isfinite(colour->r) || !std::isfinite(colour->g) ||
                !std::isfinite(colour->b))
                return {};
            const Range normalizer(std::max({colour->r, colour->g, colour->b, .001f}));
            const Range weighted = Range(star.Intensity()) * Range(star.temperature_K) / normalizer;
            if (!weighted.Finite()) return {};
            const Range x(star.direction_x), y(star.direction_y), z(star.direction_z);
            const Range normalized = z / sqrt(Square(x) + Square(y) + Square(z));
            // An equatorial source belongs to both adjacent winding
            // charts. This continuous fixed-pupil proof declines it.
            if (!normalized.Finite() || !(normalized.lo > 0.0 || normalized.hi < 0.0)) return {};
            const Range latitude = Abs(normalized);
            Range azimuth = CatalogueAzimuth(x.lo, y.lo);
            if (!azimuth.Finite()) azimuth = Range(-kPiBounds.hi, kPiBounds.hi);
            const unsigned hemisphere = normalized.lo > 0.0 ? 1 : 0;
            const unsigned band =
                std::min(kLatitudeRows - 1,
                         unsigned(std::midpoint(latitude.lo, latitude.hi) * kLatitudeRows));
            Row& row = index.rows_[hemisphere][band];
            if (row.entries.empty())
                row.latitude = latitude;
            else
                row.latitude = Range(std::min(row.latitude.lo, latitude.lo),
                                     std::max(row.latitude.hi, latitude.hi));
            row.longitude_error = std::max(row.longitude_error, MidpointError(azimuth));
            row.weighted_temperature = row.weighted_temperature + weighted;
            row.entries.push_back(
                {latitude, azimuth, weighted, std::midpoint(azimuth.lo, azimuth.hi)});
        }
        for (auto& hemisphere : index.rows_)
            for (auto& row : hemisphere) std::ranges::sort(row.entries, {}, &Entry::longitude_key);
        return index;
    }

    std::span<const Row> Rows(bool north) const { return rows_[north ? 1 : 0]; }

  private:
    std::array<std::array<Row, kLatitudeRows>, 2> rows_;
};

}  // namespace sirius::core::relativity::finite_tail_detail
