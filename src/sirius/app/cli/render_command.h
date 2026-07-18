#pragma once

// The `render` command: parses render flags into the config tree and drives the
// FSM RenderSession to written output. Ported from CRCL002A.h.

#include "sirius/app/cli/command_router.h"

namespace sirius::app {

// Renders an image via the CPU geodesic tracer.
class RenderCommand : public Command {
  public:
    std::string Name() const override { return "render"; }
    std::string Description() const override { return "Render an image"; }
    std::string Usage() const override;

    int Execute(const std::vector<std::string>& args, const GlobalOptions& globals,
                SiriusConfig& config) override;

  private:
    // Parse render-specific flags into config. Sets gpu_backend_requested_ when
    // an explicit GPU backend is selected.
    bool ParseArgs(const std::vector<std::string>& args, const GlobalOptions& globals,
                   SiriusConfig& config);

    // Print a human-readable configuration summary.
    void PrintConfig(const SiriusConfig& config, bool verbose);

    // Configure and run the render session; returns the process exit code.
    // `use_vulkan` routes the session onto the Vulkan render path.
    int ExecuteSession(const SiriusConfig& config, const GlobalOptions& globals, bool use_vulkan);

    // Set when --gpu or --backend vulkan selects the Vulkan backend. The Vulkan
    // render path runs when a device is present and declines cleanly when not.
    bool gpu_backend_requested_ = false;
};

}  // namespace sirius::app
