#pragma once

// CLI entry point: parses global options anywhere on the command line and
// routes to a subcommand, with legacy implicit-render syntax preserved. Ported
// from CRCL001A.h.

#include "sirius/app/config/config_schema.h"

#include <memory>
#include <string>
#include <vector>

namespace sirius::app {

// One CLI subcommand (render, info, config, view).
class Command {
  public:
    virtual ~Command() = default;

    // Command name as typed on the command line.
    virtual std::string Name() const = 0;

    // One-line description for the help listing.
    virtual std::string Description() const = 0;

    // Full usage text for `<command> --help`.
    virtual std::string Usage() const = 0;

    // Run the command. args are the tokens after the command name; returns the
    // process exit code (0 = success).
    virtual int Execute(const std::vector<std::string>& args, const GlobalOptions& globals,
                        SiriusConfig& config) = 0;
};

// Registers the commands and dispatches a command line to one of them.
class CommandRouter {
  public:
    CommandRouter();

    // Parse argc/argv, load configuration, and route. Returns the exit code.
    int Run(int argc, char* argv[]);

    void PrintHelp();
    void PrintVersion();

  private:
    // Extract recognised global flags; returns the remaining tokens.
    std::vector<std::string> ParseGlobalOptions(const std::vector<std::string>& args,
                                                GlobalOptions& globals);

    // First token if it names a command (or "help"); empty otherwise.
    std::string ExtractCommand(const std::vector<std::string>& args);

    // Whether the arguments use the legacy implicit-render `--flag` form.
    bool IsLegacySyntax(const std::vector<std::string>& args);

    // Dispatch to the named command (empty routes to help).
    int RouteCommand(const std::string& command_name, const std::vector<std::string>& args,
                     const GlobalOptions& globals, SiriusConfig& config);

    std::vector<std::unique_ptr<Command>> commands_;
};

}  // namespace sirius::app
