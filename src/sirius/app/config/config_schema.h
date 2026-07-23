#pragma once

// App-layer view of the render configuration tree. The structs themselves live
// in sirius/render/render_config.h (the render library consumes them and must
// not depend on the JSON codec); this header re-homes them into sirius::app as
// aliases and re-attaches the nlohmann serialisation the legacy CRCF003A.h
// carried intrusively. Ported from CRCF003A.h and CRFM001A.h.
//
// Why non-intrusive: the JSON dependency is an app/config concern. Defining
// to_json/from_json here rather than on the struct keeps sirius_render free of
// nlohmann while leaving config files byte-compatible with the legacy schema.
// Every field spelling below is copied from the legacy
// NLOHMANN_DEFINE_TYPE_INTRUSIVE_WITH_DEFAULT lists so existing config JSON
// parses identically.

#include "sirius/render/render_config.h"

#include <nlohmann/json.hpp>

// The serialisation free functions must sit in the namespace of the types they
// act on so argument-dependent lookup finds them; the types are declared in
// sirius::render. Leaf types are defined before the aggregates that contain
// them so each nested from_json/to_json is visible where it is used.
namespace sirius::render {

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(RenderConfig, width, height, samplesPerPixel,
                                                tileSize, threadCount, outputPath)

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(MetricConfig, name, mass, spin, charge, lambda,
                                                temperatureModel, diskTemperature, throatRadius,
                                                warpVelocity, bubbleRadius, bubbleSigma)

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(ObserverConfig, distance, inclination, azimuth, fov)

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(PostProcessConfig, enableBloom, bloomIntensity,
                                                bloomThreshold, exposure, contrast, saturation,
                                                tonemapper)

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(VolumetricConfig, enabled, hOverR, hPower,
                                                tauMidplane, samples, enableTurbulence,
                                                enableCorona)

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(FilmSimConfig, enabled, preset, grainIntensity,
                                                halationStrength, vignetteStrength)

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(BackendConfig, preferred, enableDenoiser,
                                                cudaDevice)

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(SiriusConfig, render, metric, observer, postprocess,
                                                backend, volumetric, film)

}  // namespace sirius::render

namespace sirius::app {

// Aliases: the app layer names the render-layer config types as its own so CLI
// and loader code reads in sirius::app terms without flattening the render
// namespace at every use site.
using RenderConfig = render::RenderConfig;
using MetricConfig = render::MetricConfig;
using ObserverConfig = render::ObserverConfig;
using PostProcessConfig = render::PostProcessConfig;
using VolumetricConfig = render::VolumetricConfig;
using FilmSimConfig = render::FilmSimConfig;
using BackendConfig = render::BackendConfig;
using SiriusConfig = render::SiriusConfig;

// Global CLI options (not persisted to a config file). Ported from CRCF003A.h;
// these never cross the JSON boundary, so they take new-tree snake_case names.
struct GlobalOptions {
    bool verbose = false;
    bool json_output = false;
    bool no_color = false;
    std::string config_path;  // Override config file path.
    bool show_help = false;
    bool show_version = false;
};

}  // namespace sirius::app
