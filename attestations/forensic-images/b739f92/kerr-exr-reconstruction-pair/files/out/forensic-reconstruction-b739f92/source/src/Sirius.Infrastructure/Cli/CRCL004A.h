// CRCL004A.h - Config Command
// Component ID: CRCL004A (Cli/ConfigCommand)
//
// CLI command for configuration management: show, validate, init.

#pragma once

#include "CRCL001A.h"

namespace Sirius::Cli {

/// @brief Config command implementation
class ConfigCommand : public Command {
public:
    std::string name() const override { return "config"; }
    std::string description() const override { return "Configuration management"; }
    std::string usage() const override;

    int execute(const std::vector<std::string>& args,
               const Configuration::GlobalOptions& globals,
               Configuration::SiriusConfig& config) override;

private:
    /// @brief Show effective configuration
    int showConfig(const Configuration::GlobalOptions& globals,
                   const Configuration::SiriusConfig& config);

    /// @brief Validate a configuration file
    int validateConfig(const std::vector<std::string>& args,
                       const Configuration::GlobalOptions& globals);

    /// @brief Create default configuration file
    int initConfig(const std::vector<std::string>& args,
                   const Configuration::GlobalOptions& globals);

    /// @brief Show config search paths
    int showPaths(const Configuration::GlobalOptions& globals);
};

} // namespace Sirius::Cli
