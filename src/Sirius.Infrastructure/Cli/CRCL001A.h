// CRCL001A.h - CLI Command Router
// Component ID: CRCL001A (Cli/Router)
//
// Main CLI entry point. Parses global options and routes to subcommands.
// Maintains backward compatibility with legacy CLI syntax.

#pragma once

#include "CRCF003A.h"

#include <string>
#include <vector>
#include <functional>

namespace Sirius::Cli {

/// @brief CLI command interface
class Command {
public:
    virtual ~Command() = default;

    /// @brief Get command name (e.g., "render", "info", "config")
    virtual std::string name() const = 0;

    /// @brief Get command description
    virtual std::string description() const = 0;

    /// @brief Get command usage/help text
    virtual std::string usage() const = 0;

    /// @brief Execute the command
    /// @param args Command arguments (after command name)
    /// @param globals Global options from parent
    /// @param config Loaded configuration
    /// @return Exit code (0 = success)
    virtual int execute(const std::vector<std::string>& args,
                       const Configuration::GlobalOptions& globals,
                       Configuration::SiriusConfig& config) = 0;
};

/// @brief Main CLI router
class CommandRouter {
public:
    CommandRouter();

    /// @brief Run the CLI with command-line arguments
    /// @param argc Argument count
    /// @param argv Argument values
    /// @return Exit code
    int run(int argc, char* argv[]);

    /// @brief Print global help
    void printHelp();

    /// @brief Print version information
    void printVersion();

private:
    /// @brief Parse global options from arguments
    /// @param args All arguments
    /// @param globals Output global options
    /// @return Remaining arguments after global options
    std::vector<std::string> parseGlobalOptions(
        const std::vector<std::string>& args,
        Configuration::GlobalOptions& globals);

    /// @brief Extract command name from arguments
    /// @param args Arguments (with global options removed)
    /// @return Command name, or empty for implicit render command
    std::string extractCommand(const std::vector<std::string>& args);

    /// @brief Check if argument looks like legacy CLI syntax
    /// @param args Arguments to check
    /// @return true if using legacy --option style without command
    bool isLegacySyntax(const std::vector<std::string>& args);

    /// @brief Route to appropriate command
    /// @param commandName Command to run (empty = render)
    /// @param args Command arguments
    /// @param globals Global options
    /// @param config Configuration
    /// @return Exit code
    int routeCommand(const std::string& commandName,
                    const std::vector<std::string>& args,
                    const Configuration::GlobalOptions& globals,
                    Configuration::SiriusConfig& config);

    /// @brief Registered commands
    std::vector<std::unique_ptr<Command>> m_commands;
};

} // namespace Sirius::Cli
