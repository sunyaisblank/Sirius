// CRCL002A.h - Render Command
// Component ID: CRCL002A (Cli/RenderCommand)
//
// CLI command for rendering images using FSM RenderSession.

#pragma once

#include "CRCL001A.h"

namespace Sirius::Cli {

/// @brief Render command implementation
class RenderCommand : public Command {
public:
    std::string name() const override { return "render"; }
    std::string description() const override { return "Render an image"; }
    std::string usage() const override;

    int execute(const std::vector<std::string>& args,
               const Configuration::GlobalOptions& globals,
               Configuration::SiriusConfig& config) override;

private:
    /// @brief Parse render-specific arguments into config
    bool parseArgs(const std::vector<std::string>& args,
                   const Configuration::GlobalOptions& globals,
                   Configuration::SiriusConfig& config);

    /// @brief Print render configuration summary
    void printConfig(const Configuration::SiriusConfig& config, bool verbose);

    /// @brief Execute render using FSM RenderSession
    int executeSession(const Configuration::SiriusConfig& config,
                       const Configuration::GlobalOptions& globals);
};

} // namespace Sirius::Cli
