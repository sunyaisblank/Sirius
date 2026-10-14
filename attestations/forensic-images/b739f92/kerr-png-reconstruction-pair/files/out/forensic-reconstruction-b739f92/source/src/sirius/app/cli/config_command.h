#pragma once

// The `config` command: show effective configuration, validate a file, write a
// default file, and list search paths. Ported from CRCL004A.h.

#include "sirius/app/cli/command_router.h"

namespace sirius::app {

// Configuration management subcommands.
class ConfigCommand : public Command {
  public:
    std::string Name() const override { return "config"; }
    std::string Description() const override { return "Configuration management"; }
    std::string Usage() const override;

    int Execute(const std::vector<std::string>& args, const GlobalOptions& globals,
                SiriusConfig& config) override;

  private:
    int ShowConfig(const GlobalOptions& globals, const SiriusConfig& config);
    int ValidateConfig(const std::vector<std::string>& args, const GlobalOptions& globals);
    int InitConfig(const std::vector<std::string>& args, const GlobalOptions& globals);
    int ShowPaths(const GlobalOptions& globals);
};

}  // namespace sirius::app
