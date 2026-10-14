// CRCL001A.cpp - CLI Command Router
// Component ID: CRCL001A (Cli/Router)

#include "CRCL001A.h"
#include "CRCL002A.h"  // Render command
#include "CRCL003A.h"  // Info command
#include "CRCL004A.h"  // Config command
#include "CRCL005A.h"  // Output utilities
#include "CRCL006A.h"  // View command (interactive viewer)
#include "CRCF002A.h"  // Config loader

#include <iostream>
#include <algorithm>

namespace Sirius::Cli {

// Version info
static const char* SIRIUS_VERSION = "1.0.0";

CommandRouter::CommandRouter() {
    // Register commands
    m_commands.push_back(std::make_unique<RenderCommand>());
    m_commands.push_back(std::make_unique<InfoCommand>());
    m_commands.push_back(std::make_unique<ConfigCommand>());
    m_commands.push_back(std::make_unique<ViewCommand>());  // Interactive viewer
}

int CommandRouter::run(int argc, char* argv[]) {
    // Convert to string vector
    std::vector<std::string> args;
    for (int i = 1; i < argc; ++i) {
        args.push_back(argv[i]);
    }

    // Parse global options
    Configuration::GlobalOptions globals;
    std::vector<std::string> remaining = parseGlobalOptions(args, globals);

    // Apply global settings
    Output::setColorEnabled(!globals.noColor);

    // Handle help/version early. '--help' beside a command shows that
    // command's usage; on its own it shows the global help.
    if (globals.showHelp) {
        std::string helpTarget = extractCommand(remaining);
        for (const auto& cmd : m_commands) {
            if (cmd->name() == helpTarget) {
                std::cout << cmd->usage() << std::endl;
                return 0;
            }
        }
        printHelp();
        return 0;
    }
    if (globals.showVersion) {
        printVersion();
        return 0;
    }

    // Load configuration
    std::optional<std::string> configPath;
    if (!globals.configPath.empty()) {
        configPath = globals.configPath;
    }
    Configuration::SiriusConfig config = Configuration::ConfigLoader::load(configPath);

    // Determine command
    std::string commandName;
    std::vector<std::string> commandArgs;

    if (remaining.empty()) {
        // No arguments - show help
        printHelp();
        return 0;
    } else if (isLegacySyntax(remaining)) {
        // Legacy syntax: sirius --width 1920 --output foo.ppm
        // Treat as implicit render command
        commandName = "render";
        commandArgs = remaining;

        if (globals.verbose) {
            Output::warning("Using legacy CLI syntax. Consider: sirius render [options]");
        }
    } else {
        // Modern syntax: sirius render --width 1920
        commandName = extractCommand(remaining);
        if (!commandName.empty() && remaining.size() > 0 && remaining[0] == commandName) {
            commandArgs = std::vector<std::string>(remaining.begin() + 1, remaining.end());
        } else {
            commandArgs = remaining;
        }
    }

    // Route to command
    return routeCommand(commandName, commandArgs, globals, config);
}

void CommandRouter::printHelp() {
    Output::banner();

    std::cout << "Usage: sirius [global-options] <command> [options]\n\n";

    std::cout << "Global Options (recognised anywhere on the command line):\n";
    std::cout << "  -v, --verbose       Enable verbose output\n";
    std::cout << "  --json              Output in JSON format\n";
    std::cout << "  --no-color          Disable colored output\n";
    std::cout << "  --config <path>     Specify config file path\n";
    std::cout << "  --help              Show this help message\n";
    std::cout << "  --version           Show version information\n";
    std::cout << "\n";

    std::cout << "Configuration layering: defaults, then config file, then\n";
    std::cout << "SIRIUS_* environment variables, then command-line flags.\n";
    std::cout << "Environment overrides: SIRIUS_WIDTH, SIRIUS_HEIGHT, SIRIUS_SAMPLES,\n";
    std::cout << "  SIRIUS_TILE_SIZE, SIRIUS_THREADS, SIRIUS_OUTPUT, SIRIUS_METRIC,\n";
    std::cout << "  SIRIUS_MASS, SIRIUS_SPIN, SIRIUS_CHARGE, SIRIUS_DISTANCE,\n";
    std::cout << "  SIRIUS_INCLINATION, SIRIUS_AZIMUTH, SIRIUS_FOV, SIRIUS_BLOOM,\n";
    std::cout << "  SIRIUS_EXPOSURE, SIRIUS_BACKEND, SIRIUS_CUDA_DEVICE\n";
    std::cout << "\n";

    std::cout << "Commands:\n";
    for (const auto& cmd : m_commands) {
        std::cout << "  " << cmd->name();
        // Pad to align descriptions
        size_t pad = 16 - cmd->name().length();
        std::cout << std::string(pad, ' ') << cmd->description() << "\n";
    }
    std::cout << "  help            Show this help message\n";
    std::cout << "\n";

    std::cout << "Run 'sirius <command> --help' for command-specific help.\n";
    std::cout << "\n";

    std::cout << "Examples:\n";
    std::cout << "  sirius render -o output.ppm -s 256\n";
    std::cout << "  sirius render -m Kerr -a 0.9 --fov 90\n";
    std::cout << "  sirius info system\n";
    std::cout << "  sirius config show\n";
}

void CommandRouter::printVersion() {
    std::cout << "Sirius v" << SIRIUS_VERSION << "\n";
    std::cout << "General Relativistic Ray Tracer\n";
    std::cout << "\n";
    std::cout << "Build info:\n";
#ifdef SIRIUS_HAS_OPTIX
    std::cout << "  OptiX backend: enabled\n";
#else
    std::cout << "  OptiX backend: disabled\n";
#endif
    std::cout << "  C++ standard: C++17\n";
}

std::vector<std::string> CommandRouter::parseGlobalOptions(
    const std::vector<std::string>& args,
    Configuration::GlobalOptions& globals)
{
    // Global flags are recognised anywhere on the command line. The previous
    // version stopped parsing them at the first command name (and its command
    // list omitted 'view'), so documented invocations such as
    // 'sirius info system --json' reached the command parser as unknown
    // options and either errored or were silently ignored. No command defines
    // a flag that collides with the global set, so position independence is
    // unambiguous.
    std::vector<std::string> remaining;

    for (size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];

        if (arg == "-v" || arg == "--verbose") {
            globals.verbose = true;
        } else if (arg == "--json") {
            globals.jsonOutput = true;
        } else if (arg == "--no-color") {
            globals.noColor = true;
        } else if (arg == "--config" && i + 1 < args.size()) {
            globals.configPath = args[++i];
        } else if (arg == "--help") {
            // Only --help; -h belongs to commands (render uses it for height)
            globals.showHelp = true;
        } else if (arg == "--version") {
            globals.showVersion = true;
        } else {
            remaining.push_back(arg);
        }
    }

    return remaining;
}

std::string CommandRouter::extractCommand(const std::vector<std::string>& args) {
    if (args.empty()) {
        return "";
    }

    const std::string& first = args[0];

    // Check if first arg is a known command
    for (const auto& cmd : m_commands) {
        if (first == cmd->name()) {
            return first;
        }
    }

    // Check special commands
    if (first == "help") {
        return "help";
    }

    // Not a command - might be legacy syntax or invalid
    return "";
}

bool CommandRouter::isLegacySyntax(const std::vector<std::string>& args) {
    if (args.empty()) {
        return false;
    }

    const std::string& first = args[0];

    // Legacy syntax starts with -- options
    if (first.substr(0, 2) == "--" || first.substr(0, 1) == "-") {
        // Check it's not a command
        for (const auto& cmd : m_commands) {
            if (first == cmd->name()) {
                return false;
            }
        }
        return true;
    }

    return false;
}

int CommandRouter::routeCommand(
    const std::string& commandName,
    const std::vector<std::string>& args,
    const Configuration::GlobalOptions& globals,
    Configuration::SiriusConfig& config)
{
    // Handle special commands
    if (commandName == "help" || commandName.empty()) {
        printHelp();
        return 0;
    }

    // Find and execute command
    for (auto& cmd : m_commands) {
        if (cmd->name() == commandName) {
            // Check for --help in command args (only --help, not -h which may be used by commands)
            for (const auto& arg : args) {
                if (arg == "--help") {
                    std::cout << cmd->usage() << std::endl;
                    return 0;
                }
            }

            return cmd->execute(args, globals, config);
        }
    }

    // Unknown command
    Output::error("Unknown command: " + commandName);
    std::cout << "\nRun 'sirius help' for available commands.\n";
    return 1;
}

} // namespace Sirius::Cli
