// CRCL006A.h - View Command (Interactive Viewer)
// Component ID: CRCL006A (Cli/ViewCommand)
//
// CLI command for launching the interactive viewer with real-time
// progressive refinement and camera controls.

#pragma once

#include "CRCL001A.h"

namespace Sirius::Cli {

/// @brief View command implementation (interactive viewer)
class ViewCommand : public Command {
public:
    std::string name() const override { return "view"; }
    std::string description() const override { return "Launch interactive viewer"; }
    std::string usage() const override;

    int execute(const std::vector<std::string>& args,
               const Configuration::GlobalOptions& globals,
               Configuration::SiriusConfig& config) override;

private:
    /// @brief Parse view-specific arguments
    bool parseArgs(const std::vector<std::string>& args,
                   const Configuration::GlobalOptions& globals,
                   Configuration::SiriusConfig& config);
};

} // namespace Sirius::Cli
