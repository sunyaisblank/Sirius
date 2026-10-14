// CRCL003A.h - Info Command
// Component ID: CRCL003A (Cli/InfoCommand)
//
// CLI command for displaying system and configuration information.

#pragma once

#include "CRCL001A.h"

namespace Sirius::Cli {

/// @brief Info command implementation
class InfoCommand : public Command {
public:
    std::string name() const override { return "info"; }
    std::string description() const override { return "Display system information"; }
    std::string usage() const override;

    int execute(const std::vector<std::string>& args,
               const Configuration::GlobalOptions& globals,
               Configuration::SiriusConfig& config) override;

private:
    /// @brief Display system information (GPU, CUDA, OptiX)
    int showSystem(const Configuration::GlobalOptions& globals);

    /// @brief Display available metrics
    int showMetrics(const Configuration::GlobalOptions& globals);

    /// @brief Display current configuration
    int showConfig(const Configuration::GlobalOptions& globals,
                   const Configuration::SiriusConfig& config);

    /// @brief Display available backends
    int showBackends(const Configuration::GlobalOptions& globals);
};

} // namespace Sirius::Cli
