#pragma once

// The `info` command: system capabilities, the metric catalogue, the effective
// configuration, and the render backends. Ported from CRCL003A.h. OptiX is
// retired; backend reporting is honest about what is compiled in.

#include "sirius/app/cli/command_router.h"

namespace sirius::app {

// Displays system, metric, configuration, and backend information.
class InfoCommand : public Command {
  public:
    std::string Name() const override { return "info"; }
    std::string Description() const override { return "Display system information"; }
    std::string Usage() const override;

    int Execute(const std::vector<std::string>& args, const GlobalOptions& globals,
                SiriusConfig& config) override;

  private:
    int ShowSystem(const GlobalOptions& globals);
    int ShowMetrics(const GlobalOptions& globals);
    int ShowConfig(const GlobalOptions& globals, const SiriusConfig& config);
    int ShowBackends(const GlobalOptions& globals);
};

}  // namespace sirius::app
