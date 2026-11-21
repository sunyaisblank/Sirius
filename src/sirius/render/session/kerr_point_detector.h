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

// Prepare the immutable sky data once for the fixed render catalogue.
[[nodiscard]] std::shared_ptr<const KerrDetectorCatalogue> PrepareKerrDetectorCatalogue(
    const core::StarfieldSpatialIndex& stars);

[[nodiscard]] std::optional<std::array<double, 3>> BoundKerrDetectorCell(
    const KerrDetectorCatalogue& catalogue, double brightness, const core::CameraConfig& camera,
    const backend::DetectorLaunchFrame& launch, const DetectorCoordinate& film_origin,
    const PointDetectorFootprint& footprint, const DetectorCoordinate& lower,
    const DetectorCoordinate& upper, const std::function<bool()>& cancelled,
    std::size_t& remaining);

}  // namespace sirius::render
