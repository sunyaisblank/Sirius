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
    int ExecuteSession(const SiriusConfig& config, const GlobalOptions& globals);

    // Set when --gpu or --backend selects a GPU backend; the CPU path is the
    // only wired one, so such a request declines cleanly.
    bool gpu_backend_requested_ = false;
};

}  // namespace sirius::app
