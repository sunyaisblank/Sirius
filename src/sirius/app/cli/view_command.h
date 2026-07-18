#pragma once

// The `view` command: launches the OpenGL interactive viewer with progressive
// refinement and orbital camera controls. Ported from CRCL006A.h.

#include "sirius/app/cli/command_router.h"

namespace sirius::app {

// Launches the real-time viewer window.
class ViewCommand : public Command {
  public:
    std::string Name() const override { return "view"; }
    std::string Description() const override { return "Launch interactive viewer"; }
    std::string Usage() const override;

    int Execute(const std::vector<std::string>& args, const GlobalOptions& globals,
                SiriusConfig& config) override;

  private:
    bool ParseArgs(const std::vector<std::string>& args, const GlobalOptions& globals,
                   SiriusConfig& config);
};

}  // namespace sirius::app
