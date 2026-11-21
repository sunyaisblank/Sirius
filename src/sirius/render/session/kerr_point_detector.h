#pragma once

#include "sirius/render/session/point_source_detector.h"

#include <memory>

namespace sirius::backend {
struct DetectorLaunchFrame;
}

namespace sirius::core {
struct CameraConfig;
}

namespace sirius::render {

class KerrDetectorCatalogue;

struct PointDetectorCell {
    DetectorCoordinate lower, upper;
};

// Prepare the immutable sky data once for the fixed render catalogue.
[[nodiscard]] std::shared_ptr<const KerrDetectorCatalogue> PrepareKerrDetectorCatalogue(
    const core::StarfieldSpatialIndex& stars);

// Bound one or more cells of the same original Gaussian using a single
// inverse chart and catalogue search. Cells must be interior-disjoint.
[[nodiscard]] std::optional<std::array<double, 3>> BoundKerrDetectorCells(
    const KerrDetectorCatalogue& catalogue, double brightness, const core::CameraConfig& camera,
    const backend::DetectorLaunchFrame& launch, const DetectorCoordinate& film_origin,
    const PointDetectorFootprint& footprint, std::span<const PointDetectorCell> cells,
    const std::function<bool()>& cancelled, std::size_t& remaining);

}  // namespace sirius::render
