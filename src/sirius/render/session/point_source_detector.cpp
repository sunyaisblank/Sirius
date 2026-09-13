#include "sirius/render/session/point_source_detector.h"

#include "sirius/core/observer_frame.h"
#include "sirius/core/spectral/point_source_transfer.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <numbers>
#include <optional>
#include <utility>
#include <vector>

namespace sirius::render {
namespace {
using Coordinate = DetectorCoordinate;
using Matrix = core::AngularMatrix2;
using Rgb = std::array<double, 3>;
using Direction = std::array<double, 3>;

struct Cell {
    Coordinate lower, upper;
    unsigned depth;
    Coordinate Centre() const { return {(lower[0] + upper[0]) * .5, (lower[1] + upper[1]) * .5}; }
    double Width() const { return upper[0] - lower[0]; }
    bool IntersectsSupport() const {
        return std::hypot(std::clamp(0.0, lower[0], upper[0]),
                          std::clamp(0.0, lower[1], upper[1])) <= 4;
    }
    bool Contains(Coordinate z, double margin = 0) const {
        return z[0] >= lower[0] - margin && z[0] <= upper[0] + margin &&
               z[1] >= lower[1] - margin && z[1] <= upper[1] + margin;
    }
    std::array<Cell, 4> Children() const {
        const auto c = Centre();
        return {{{lower, c, depth + 1},
                 {{c[0], lower[1]}, {upper[0], c[1]}, depth + 1},
                 {{lower[0], c[1]}, {c[0], upper[1]}, depth + 1},
                 {c, upper, depth + 1}}};
    }
};

std::optional<Coordinate> Solve(const Matrix& a, Coordinate b) {
    const double scale =
        std::max({std::abs(a[0][0]), std::abs(a[0][1]), std::abs(a[1][0]), std::abs(a[1][1])});
    if (!(scale > 0) || !std::isfinite(scale)) return std::nullopt;
    const double aa = a[0][0] / scale, ab = a[0][1] / scale;
    const double ba = a[1][0] / scale, bb = a[1][1] / scale;
    const double determinant = std::fma(aa, bb, -ab * ba);
    if (std::abs(determinant) <=
        32 * std::numeric_limits<double>::epsilon() * (std::abs(aa * bb) + std::abs(ab * ba)))
        return std::nullopt;
    b[0] /= scale;
    b[1] /= scale;
    Coordinate result{std::fma(bb, b[0], -ab * b[1]) / determinant,
                      std::fma(aa, b[1], -ba * b[0]) / determinant};
    if (!std::isfinite(result[0]) || !std::isfinite(result[1])) return std::nullopt;
    return result;
}

std::optional<Coordinate> SkyOffset(const Direction& from, const Direction& to) {
    const auto basis = core::relativity::MakeCelestialTangentBasis(from);
    if (!basis) return std::nullopt;
    const auto separation = core::relativity::MeasureCelestialSeparation(from, to);
    if (!std::isfinite(separation.angle)) return std::nullopt;
    if (separation.sine == 0) {
        if (separation.angle == 0) return Coordinate{};
        return std::nullopt;
    }
    const auto& n = separation.normal;
    const Direction tangent{(n[1] * from[2] - n[2] * from[1]) / separation.sine,
                            (n[2] * from[0] - n[0] * from[2]) / separation.sine,
                            (n[0] * from[1] - n[1] * from[0]) / separation.sine};
    Coordinate result{};
    for (int axis = 0; axis < 3; ++axis) {
        result[0] += separation.angle * tangent[axis] * basis->first[axis];
        result[1] += separation.angle * tangent[axis] * basis->second[axis];
    }
    return result;
}

struct CachedProbe {
    Coordinate z;
    PointDetectorProbe value;
};
struct ImageRoot {
    std::uint32_t star;
    Coordinate z;
    std::size_t probe;
    double uncertainty;
    std::size_t previous_probe;
};
struct Estimate {
    Rgb rgb{};
    Rgb error{};
    std::vector<std::size_t> images;
    bool regular = true;
};

class Detector {
  public:
    Detector(const core::StarfieldSpatialIndex& catalogue, double brightness,
             const PointDetectorSampler& sampler, const std::function<bool()>& cancelled,
             const PointDetectorPolicy& policy)
        : catalogue_(catalogue),
          brightness_(brightness),
          sampler_(sampler),
          cancelled_(cancelled),
          policy_(policy) {}

    std::expected<PointDetectorResult, PointDetectorError> Run() {
        try {
            if (!sampler_ || !std::isfinite(brightness_) || brightness_ < 0 ||
                !Positive(policy_.absolute_rgb_error) || !Positive(policy_.relative_rgb_error) ||
                !Positive(policy_.geometry_error) || !Positive(policy_.root_error) ||
                !Positive(policy_.maximum_linearization_residual) ||
                policy_.maximum_linearization_residual >= 1 || policy_.maximum_probes == 0 ||
                policy_.maximum_probes > 65536 || policy_.maximum_cells == 0 ||
                policy_.maximum_candidate_visits == 0 ||
                policy_.maximum_candidate_visits > 1048576 || policy_.maximum_newton_steps == 0 ||
                policy_.maximum_depth > 20 || policy_.minimum_depth > policy_.maximum_depth)
                Fail(PointDetectorFailure::InvalidInput);
            probes_.reserve(policy_.maximum_probes);
            roots_.reserve(policy_.maximum_candidate_visits);
            statistics_.reserved_cache_bytes =
                probes_.capacity() * sizeof(CachedProbe) + roots_.capacity() * sizeof(ImageRoot);
            CheckCancellation();
            const Cell support{{-4, -4}, {4, 4}, 0};
            auto result = Refine(support, Measure(support));
            statistics_.roots = roots_.size();
            return PointDetectorResult{result.rgb, result.error, statistics_};
        } catch (PointDetectorFailure reason) {
            statistics_.roots = roots_.size();
            return std::unexpected(PointDetectorError{reason, statistics_});
        } catch (const std::bad_alloc&) {
            statistics_.roots = roots_.size();
            return std::unexpected(
                PointDetectorError{PointDetectorFailure::WorkLimit, statistics_});
        }
    }

  private:
    static bool Positive(double value) { return std::isfinite(value) && value > 0; }
    [[noreturn]] static void Fail(PointDetectorFailure reason) { throw reason; }
    void CheckCancellation() const {
        if (cancelled_ && cancelled_()) Fail(PointDetectorFailure::Cancelled);
    }
    std::size_t Probe(Coordinate z) {
        CheckCancellation();
        for (std::size_t i = 0; i < probes_.size(); ++i)
            if (probes_[i].z == z) return i;
        if (probes_.size() == policy_.maximum_probes) Fail(PointDetectorFailure::WorkLimit);
        auto value = sampler_(z);
        ++statistics_.probes;
        CheckCancellation();
        if (!value) Fail(value.error());
        statistics_.inner_attempts += value->inner_attempts;
        statistics_.tail_attempts += value->tail_attempts;
        if (value->visible) {
            const double norm =
                std::hypot(value->direction[0], value->direction[1], value->direction[2]);
            if (!Positive(norm) || !Positive(value->camera_over_source_frequency) ||
                !std::isfinite(value->transmission) || value->transmission < 0 ||
                value->transmission > 1)
                Fail(PointDetectorFailure::Arithmetic);
            for (auto& component : value->direction) component /= norm;
            for (const auto& row : value->source_derivative)
                for (double component : row)
                    if (!std::isfinite(component)) Fail(PointDetectorFailure::Arithmetic);
        }
        probes_.push_back({z, *value});
        return probes_.size() - 1;
    }

    std::optional<std::size_t> FindImage(std::uint32_t star, const Cell& cell, std::size_t seed,
                                         bool regular, double margin, bool& unresolved) {
        const auto& entry = catalogue_.Stars()[star];
        const Direction target{entry.direction_x, entry.direction_y, entry.direction_z};
        // A regular validated cell has a single local branch. Reuse its existing
        // original-coordinate root at a shared edge; never merge by star ID.
        if (regular) {
            for (std::size_t i = 0; i < roots_.size(); ++i)
                if (roots_[i].star == star && cell.Contains(roots_[i].z, roots_[i].uncertainty))
                    return i;
        }
        auto index = seed;
        Coordinate z = probes_[index].z;
        for (unsigned iteration = 0; iteration < policy_.maximum_newton_steps; ++iteration) {
            ++statistics_.newton_steps;
            const auto& current = probes_[index].value;
            if (!current.visible) return std::nullopt;
            const auto residual = SkyOffset(current.direction, target);
            const auto delta =
                residual ? Solve(current.source_derivative, *residual) : std::nullopt;
            if (!delta) {
                unresolved = true;
                return std::nullopt;
            }
            const double error = std::hypot((*delta)[0], (*delta)[1]);
            if (error <= policy_.root_error) {
                // Validate with an independently evaluated corrected image ray.
                const Coordinate corrected{z[0] + (*delta)[0], z[1] + (*delta)[1]};
                if (!cell.Contains(corrected, margin)) return std::nullopt;
                const auto validated = Probe(corrected);
                const auto& actual = probes_[validated].value;
                if (!actual.visible) return std::nullopt;
                const auto remainder = SkyOffset(actual.direction, target);
                const auto remaining =
                    remainder ? Solve(actual.source_derivative, *remainder) : std::nullopt;
                if (!remaining ||
                    std::hypot((*remaining)[0], (*remaining)[1]) > policy_.root_error) {
                    unresolved = true;
                    return std::nullopt;
                }
                const double uncertainty = std::hypot((*remaining)[0], (*remaining)[1]) +
                                           32 * std::numeric_limits<double>::epsilon() *
                                               (1 + std::hypot(corrected[0], corrected[1]));
                for (std::size_t i = 0; i < roots_.size(); ++i) {
                    if (roots_[i].star != star) continue;
                    if (roots_[i].z == corrected) return i;
                    if (std::hypot(roots_[i].z[0] - corrected[0], roots_[i].z[1] - corrected[1]) <=
                        roots_[i].uncertainty + uncertainty) {
                        if (regular) return i;
                        unresolved = true;
                        return std::nullopt;
                    }
                }
                if (roots_.size() == policy_.maximum_candidate_visits)
                    Fail(PointDetectorFailure::WorkLimit);
                statistics_.maximum_root_coordinate_error =
                    std::max(statistics_.maximum_root_coordinate_error, uncertainty);
                roots_.push_back({star, corrected, validated, uncertainty, index});
                return roots_.size() - 1;
            }
            Coordinate next{z[0] + (*delta)[0], z[1] + (*delta)[1]};
            if (regular && !cell.Contains(next, margin)) return std::nullopt;
            // A bounded trust region keeps nonlinear trials on this cell's
            // branch. Failed visibility trials backtrack toward the seed.
            double weight = 1;
            bool advanced = false;
            for (unsigned backtrack = 0; backtrack < 12; ++backtrack) {
                next = {z[0] + weight * (*delta)[0], z[1] + weight * (*delta)[1]};
                if (cell.Contains(next, margin) && next != z) {
                    const auto candidate = Probe(next);
                    if (probes_[candidate].value.visible) {
                        const auto candidate_residual =
                            SkyOffset(probes_[candidate].value.direction, target);
                        if (candidate_residual &&
                            std::hypot((*candidate_residual)[0], (*candidate_residual)[1]) <
                                std::hypot((*residual)[0], (*residual)[1])) {
                            z = next;
                            index = candidate;
                            advanced = true;
                            break;
                        }
                    }
                }
                weight *= .5;
            }
            if (!advanced) {
                unresolved = true;
                return std::nullopt;
            }
        }
        unresolved = true;
        return std::nullopt;
    }

    Estimate Measure(const Cell& cell) {
        Estimate result;
        if (!cell.IntersectsSupport()) return result;
        if (++statistics_.cells > policy_.maximum_cells) Fail(PointDetectorFailure::WorkLimit);
        CheckCancellation();
        const auto centre = cell.Centre();
        const double half = cell.Width() * .5;
        std::array<std::size_t, 13> nodes{};
        std::size_t count = 0;
        nodes[count++] = Probe(centre);
        for (int y = -1; y <= 1; ++y)
            for (int x = -1; x <= 1; ++x)
                if (x != 0 || y != 0)
                    nodes[count++] = Probe({centre[0] + x * half, centre[1] + y * half});
        // Staggered positions are not the quadtree's next-level corners.
        for (const Coordinate shift : {Coordinate{-.37, -.61}, Coordinate{.61, -.37},
                                       Coordinate{.37, .61}, Coordinate{-.61, .37}})
            nodes[count++] = Probe({centre[0] + shift[0] * half, centre[1] + shift[1] * half});
        const auto seed = std::find_if(nodes.begin(), nodes.end(), [&](std::size_t i) {
            return probes_[i].value.visible &&
                   Solve(probes_[i].value.source_derivative, {1, 0}).has_value();
        });
        if (seed == nodes.end()) {
            result.regular = std::none_of(nodes.begin(), nodes.end(),
                                          [&](std::size_t i) { return probes_[i].value.visible; });
            return result;
        }
        const auto& anchor = probes_[*seed];
        double residual_bound = 0;
        Coordinate residual_by_axis{};
        Coordinate source_residual{};
        Coordinate source_reach{};
        bool mixed = false;
        const auto& matrix = anchor.value.source_derivative;
        const double parity = std::fma(matrix[0][0], matrix[1][1], -matrix[0][1] * matrix[1][0]);
        for (auto i : nodes) {
            const auto& other = probes_[i];
            if (!other.value.visible) {
                mixed = true;
                continue;
            }
            const auto offset = SkyOffset(anchor.value.direction, other.value.direction);
            if (offset) {
                for (int axis = 0; axis < 2; ++axis) {
                    const double predicted = matrix[axis][0] * (other.z[0] - anchor.z[0]) +
                                             matrix[axis][1] * (other.z[1] - anchor.z[1]);
                    source_residual[axis] =
                        std::max(source_residual[axis], std::abs((*offset)[axis] - predicted));
                }
            }
            const auto inferred = offset ? Solve(matrix, *offset) : std::nullopt;
            if (!inferred) {
                result.regular = false;
                continue;
            }
            const double error = std::hypot((*inferred)[0] - (other.z[0] - anchor.z[0]),
                                            (*inferred)[1] - (other.z[1] - anchor.z[1]));
            for (int axis = 0; axis < 2; ++axis)
                residual_by_axis[axis] =
                    std::max(residual_by_axis[axis],
                             std::abs((*inferred)[axis] - (other.z[axis] - anchor.z[axis])));
            residual_bound = std::max(residual_bound, error);
            const auto& b = other.value.source_derivative;
            const double other_parity = std::fma(b[0][0], b[1][1], -b[0][1] * b[1][0]);
            // Curvature controls branch discovery and its query padding. Flux
            // uses the retraced image and its local derivative, so it does not
            // inherit this affine approximation's residual as a radiance error.
            if (!(parity * other_parity > 0) ||
                error > policy_.maximum_linearization_residual * cell.Width())
                result.regular = false;
        }
        // Visibility discontinuities require extra spatial sampling. Actual
        // image rays decide transmission; the centre cannot mask the support.
        if (mixed && cell.depth < policy_.minimum_depth + 2) {
            // This cell cannot be admitted regardless of its catalogue sum.
            // Refine visibility before spending visits on its broad query;
            // only represented descendant estimates can own radiance.
            result.regular = false;
            return result;
        }
        double reach = 0;
        for (auto i : nodes) {
            if (!probes_[i].value.visible) continue;
            const auto& b = probes_[i].value.source_derivative;
            reach = std::max(reach, std::hypot(b[0][0], b[0][1]) + std::hypot(b[1][0], b[1][1]));
        }
        // Keep the two source directions separate near a fold. Inverting an
        // almost singular map can make a distant catalogue point appear to
        // need arbitrarily deep film refinement. The forward map still has a
        // small, measured range in its critical direction. Express every node
        // in the anchor's sky chart so tangent-basis switches cannot mix axes.
        for (auto i : nodes) {
            if (!probes_[i].value.visible) continue;
            const auto offset = SkyOffset(anchor.value.direction, probes_[i].value.direction);
            if (!offset) continue;
            for (int axis = 0; axis < 2; ++axis)
                source_reach[axis] = std::max(source_reach[axis], std::abs((*offset)[axis]));
        }
        for (int axis = 0; axis < 2; ++axis)
            source_reach[axis] += 2 * source_residual[axis] +
                                  policy_.root_error * std::hypot(matrix[axis][0], matrix[axis][1]);
        double query_radius = reach * (2 * cell.Width() + residual_bound);
        if (!mixed) {
            // The exact same rectangular forward-range predicate below
            // already rejects everything outside this circumscribed disk.
            // Apply that bound to the index query as well, without dropping
            // a candidate which could pass the existing predicate.
            query_radius =
                std::min(query_radius, std::nextafter(std::hypot(source_reach[0], source_reach[1]),
                                                      std::numeric_limits<double>::infinity()));
        }
        const float query_sigma = std::nextafter(
            static_cast<float>(query_radius * .25 + 8 * std::numeric_limits<float>::epsilon()),
            std::numeric_limits<float>::infinity());
        bool unresolved = false;
        bool has_candidate = false;
        const auto& n = anchor.value.direction;
        catalogue_.ForEachCandidateWhile(
            static_cast<float>(n[0]), static_cast<float>(n[1]), static_cast<float>(n[2]),
            query_sigma, [&](std::uint32_t star) {
                if (++statistics_.candidate_visits > policy_.maximum_candidate_visits)
                    Fail(PointDetectorFailure::WorkLimit);
                CheckCancellation();
                const auto& candidate = catalogue_.Stars()[star];
                const auto separation = core::relativity::MeasureCelestialSeparation(
                    n,
                    Direction{candidate.direction_x, candidate.direction_y, candidate.direction_z});
                if (separation.angle > query_radius) return true;
                const auto offset = SkyOffset(
                    n,
                    Direction{candidate.direction_x, candidate.direction_y, candidate.direction_z});
                const auto predicted = offset ? Solve(matrix, *offset) : std::nullopt;
                if (!mixed && offset &&
                    (std::abs((*offset)[0]) > source_reach[0] ||
                     std::abs((*offset)[1]) > source_reach[1]))
                    return true;
                if (predicted) {
                    for (int axis = 0; axis < 2; ++axis) {
                        const double value = anchor.z[axis] + (*predicted)[axis];
                        const double margin = 2 * residual_by_axis[axis] + policy_.root_error;
                        if (value < cell.lower[axis] - margin || value > cell.upper[axis] + margin)
                            return true;
                    }
                }
                has_candidate = true;
                // An irregular cell only needs an existence witness before
                // subdivision; further catalogue visits cannot change its estimate.
                if (!result.regular) return false;
                const auto image = FindImage(star, cell, *seed, true,
                                             2 * residual_bound + policy_.root_error, unresolved);
                if (!image) return true;
                const auto& root = roots_[*image];
                const auto& point = probes_[root.probe].value;
                const auto response = core::MakeRestrictedAffinePointResponse(
                    point.source_derivative, cell.lower, cell.upper, policy_.geometry_error);
                if (!response) {
                    unresolved = true;
                    return true;
                }
                const auto density = response->DensityAtOriginalRoot(root.z);
                if (!density) Fail(PointDetectorFailure::Arithmetic);
                if (*density == 0) return true;
                const auto& entry = catalogue_.Stars()[star];
                const auto rgb = core::spectral::TransferPointSourceBand(
                    entry.temperature_K, point.camera_over_source_frequency,
                    static_cast<double>(entry.Intensity()) * brightness_, *density);
                if (!rgb) Fail(PointDetectorFailure::Arithmetic);
                const auto& previous = probes_[root.previous_probe];
                const double correction =
                    std::hypot(previous.z[0] - root.z[0], previous.z[1] - root.z[1]);
                Rgb prior_rgb{};
                if (correction > 0) {
                    const auto prior_response = core::MakeRestrictedAffinePointResponse(
                        previous.value.source_derivative, cell.lower, cell.upper,
                        policy_.geometry_error);
                    if (!prior_response) {
                        unresolved = true;
                        return true;
                    }
                    const auto prior_density = prior_response->DensityAtOriginalRoot(root.z);
                    if (!prior_density) Fail(PointDetectorFailure::Arithmetic);
                    const auto transferred = core::spectral::TransferPointSourceBand(
                        entry.temperature_K, previous.value.camera_over_source_frequency,
                        static_cast<double>(entry.Intensity()) * brightness_, *prior_density);
                    if (!transferred) Fail(PointDetectorFailure::Arithmetic);
                    for (int channel = 0; channel < 3; ++channel)
                        prior_rgb[channel] = (*transferred)[channel] * previous.value.transmission;
                }
                result.images.push_back(*image);
                for (int channel = 0; channel < 3; ++channel) {
                    const double contribution = (*rgb)[channel] * point.transmission;
                    result.rgb[channel] += contribution;
                    // Estimate root-location, local transfer and represented
                    // response errors separately from level differences.
                    const double smooth_error =
                        correction > 0 ? 2 * std::abs(contribution - prior_rgb[channel]) *
                                             root.uncertainty / correction
                                       : 0;
                    const double weight_error =
                        std::expm1(4 * root.uncertainty + .5 * root.uncertainty * root.uncertainty);
                    result.error[channel] +=
                        smooth_error +
                        contribution * (weight_error + response->query.arithmetic_area_bound +
                                        128 * std::numeric_limits<double>::epsilon());
                    if (!std::isfinite(result.rgb[channel]) ||
                        !std::isfinite(result.error[channel]))
                        Fail(PointDetectorFailure::Arithmetic);
                }
                return true;
            });
        result.regular = !unresolved && (result.regular || !has_candidate);
        std::sort(result.images.begin(), result.images.end());
        return result;
    }

    static Estimate Sum(const std::array<Estimate, 4>& children) {
        Estimate result;
        for (const auto& child : children) {
            result.regular = result.regular && child.regular;
            result.images.insert(result.images.end(), child.images.begin(), child.images.end());
            for (int channel = 0; channel < 3; ++channel) {
                result.rgb[channel] += child.rgb[channel];
                result.error[channel] += child.error[channel];
            }
        }
        std::sort(result.images.begin(), result.images.end());
        return result;
    }
    bool Agrees(const Cell& cell, const Estimate& coarse, Estimate& fine) const {
        if (!coarse.regular || !fine.regular || coarse.images != fine.images) return false;
        bool accepted = true;
        // A LOWER Gaussian-mass bound allocates no more than the original
        // packet allowance across disjoint leaves. Boundary cells retain the
        // relative allowance; they do not gain a fresh absolute allowance.
        const double far_radius =
            std::hypot(std::max(std::abs(cell.lower[0]), std::abs(cell.upper[0])),
                       std::max(std::abs(cell.lower[1]), std::abs(cell.upper[1])));
        const double mass_bound = far_radius <= 4 ? cell.Width() * cell.Width() *
                                                        std::exp(-.5 * far_radius * far_radius) /
                                                        (2 * std::numbers::pi * -std::expm1(-8.0))
                                                  : 0;
        for (int channel = 0; channel < 3; ++channel) {
            fine.error[channel] =
                std::max(fine.error[channel], std::abs(fine.rgb[channel] - coarse.rgb[channel]));
            if (fine.error[channel] > policy_.absolute_rgb_error * mass_bound +
                                          policy_.relative_rgb_error * fine.rgb[channel])
                accepted = false;
        }
        return accepted;
    }
    Estimate Refine(const Cell& cell, const Estimate& coarse) {
        if (!cell.IntersectsSupport()) return coarse;
        if (cell.depth >= policy_.maximum_depth) Fail(PointDetectorFailure::Unresolved);
        const auto children = cell.Children();
        std::array<Estimate, 4> estimates;
        for (int i = 0; i < 4; ++i) estimates[i] = Measure(children[i]);
        auto fine = Sum(estimates);
        if (cell.depth >= policy_.minimum_depth && cell.depth + 2 <= policy_.maximum_depth &&
            Agrees(cell, coarse, fine)) {
            std::array<Estimate, 4> validation;
            for (int i = 0; i < 4; ++i) {
                const auto grandchildren = children[i].Children();
                std::array<Estimate, 4> leaves;
                for (int j = 0; j < 4; ++j) leaves[j] = Measure(grandchildren[j]);
                validation[i] = Sum(leaves);
            }
            auto checked = Sum(validation);
            if (Agrees(cell, fine, checked)) return checked;
        }
        for (int i = 0; i < 4; ++i) estimates[i] = Refine(children[i], estimates[i]);
        return Sum(estimates);
    }

    const core::StarfieldSpatialIndex& catalogue_;
    double brightness_;
    const PointDetectorSampler& sampler_;
    const std::function<bool()>& cancelled_;
    const PointDetectorPolicy& policy_;
    PointDetectorStatistics statistics_;
    std::vector<CachedProbe> probes_;
    std::vector<ImageRoot> roots_;
};
}  // namespace

std::expected<PointDetectorResult, PointDetectorError> EvaluatePointDetector(
    const core::StarfieldSpatialIndex& catalogue, double brightness_scale,
    const PointDetectorSampler& sample, const std::function<bool()>& cancelled,
    const PointDetectorPolicy& policy) {
    return Detector(catalogue, brightness_scale, sample, cancelled, policy).Run();
}

}  // namespace sirius::render
