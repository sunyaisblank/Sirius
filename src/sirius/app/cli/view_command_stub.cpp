// Fail-closed implementation for a build host without the OpenGL development
// surface. The rest of the CLI and both render backends remain available; the
// missing optional capability is represented explicitly at its command seam.

#include "sirius/app/cli/cli_output.h"
#include "sirius/app/cli/view_command.h"
#include "sirius/app/config/config_loader.h"

#include <cmath>
#include <stdexcept>

namespace sirius::app {

namespace cli = cli_output;

namespace {

const std::string& RequireValue(const std::vector<std::string>& args, std::size_t& index,
                                const std::string& option) {
    if (++index >= args.size()) {
        throw std::invalid_argument(option + " requires a value");
    }
    return args[index];
}

int ParseInteger(const std::string& text) {
    std::size_t consumed = 0;
    const int value = std::stoi(text, &consumed);
    if (consumed != text.size()) throw std::invalid_argument("trailing characters");
    return value;
}

double ParseDouble(const std::string& text) {
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != text.size() || !std::isfinite(value)) {
        throw std::invalid_argument("expected one finite number");
    }
    return value;
}

}  // namespace

std::string ViewCommand::Usage() const {
    return R"(Usage: sirius view [options]

The interactive viewer is unavailable in this build. Reconfigure with
SIRIUS_BUILD_VIEWER=ON on a host with OpenGL development files.
)";
}

bool ViewCommand::ParseArgs(const std::vector<std::string>& args, const GlobalOptions& /*globals*/,
                            SiriusConfig& config) {
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];
        try {
            if (arg == "--width") {
                config.render.width = ParseInteger(RequireValue(args, i, arg));
            } else if (arg == "--height") {
                config.render.height = ParseInteger(RequireValue(args, i, arg));
            } else if (arg == "--spin") {
                config.metric.spin = ParseDouble(RequireValue(args, i, arg));
                config.metric.name = config.metric.spin == 0.0 ? "Schwarzschild" : "Kerr";
            } else if (arg == "--distance") {
                config.observer.distance = ParseDouble(RequireValue(args, i, arg));
            } else if (arg == "--inclination") {
                config.observer.inclination = ParseDouble(RequireValue(args, i, arg));
            } else if (arg == "--fov") {
                config.observer.fov = ParseDouble(RequireValue(args, i, arg));
            } else if (arg == "--no-disk") {
                config.disk_enabled = false;
            } else if (arg == "--jets") {
                jets_enabled_ = true;
            } else if (arg == "--cpu") {
                config.backend.preferred = "cpu";
            } else if (arg == "--gpu") {
                config.backend.preferred = "vulkan";
            } else if (arg == "--backend") {
                config.backend.preferred = RequireValue(args, i, arg);
            } else {
                cli::Error("Unknown view option: " + arg);
                return false;
            }
        } catch (const std::exception& error) {
            cli::Error("Invalid value for " + arg + ": " + error.what());
            return false;
        }
    }
    return true;
}

int ViewCommand::Execute(const std::vector<std::string>& args, const GlobalOptions& globals,
                         SiriusConfig& config) {
    jets_enabled_ = false;
    if (!ParseArgs(args, globals, config)) return 1;

    const auto errors = ConfigLoader::Validate(config);
    if (!errors.empty()) {
        for (const auto& error : errors) cli::Error(error);
        return 1;
    }
    if (jets_enabled_) {
        cli::Error(
            "Relativistic jets require covariant geodesic radiative transfer, which is not "
            "represented");
        return 1;
    }

    cli::Error(
        "Interactive viewer was not compiled; CPU and Vulkan render commands remain "
        "available");
    return 1;
}

}  // namespace sirius::app
