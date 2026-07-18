// CLI command router. Ported from CRCL001A.cpp.

#include "sirius/app/cli/command_router.h"

#include "sirius/app/cli/cli_output.h"
#include "sirius/app/cli/config_command.h"
#include "sirius/app/cli/info_command.h"
#include "sirius/app/cli/render_command.h"
#include "sirius/app/cli/view_command.h"
#include "sirius/app/config/config_loader.h"

#include <iostream>
#include <optional>

namespace sirius::app {

namespace cli = cli_output;

// Version string for `--version`.
static const char* kSiriusVersion = "2.0.0";

CommandRouter::CommandRouter() {
    commands_.push_back(std::make_unique<RenderCommand>());
    commands_.push_back(std::make_unique<InfoCommand>());
    commands_.push_back(std::make_unique<ConfigCommand>());
    commands_.push_back(std::make_unique<ViewCommand>());
}

int CommandRouter::Run(int argc, char* argv[]) {
    std::vector<std::string> args;
    for (int i = 1; i < argc; ++i) {
        args.push_back(argv[i]);
    }

    GlobalOptions globals;
    std::vector<std::string> remaining = ParseGlobalOptions(args, globals);

    cli::SetColorEnabled(!globals.no_color);

    // '--help' beside a command shows that command's usage; on its own it shows
    // the global help.
    if (globals.show_help) {
        std::string help_target = ExtractCommand(remaining);
        for (const auto& cmd : commands_) {
            if (cmd->Name() == help_target) {
                std::cout << cmd->Usage() << std::endl;
                return 0;
            }
        }
        PrintHelp();
        return 0;
    }
    if (globals.show_version) {
        PrintVersion();
        return 0;
    }

    std::optional<std::string> config_path;
    if (!globals.config_path.empty()) {
        config_path = globals.config_path;
    }
    SiriusConfig config = ConfigLoader::Load(config_path);

    std::string command_name;
    std::vector<std::string> command_args;

    if (remaining.empty()) {
        PrintHelp();
        return 0;
    } else if (IsLegacySyntax(remaining)) {
        // Legacy: sirius --width 1920 --output foo.ppm, treated as render.
        command_name = "render";
        command_args = remaining;

        if (globals.verbose) {
            cli::Warning("Using legacy CLI syntax. Consider: sirius render [options]");
        }
    } else {
        command_name = ExtractCommand(remaining);
        if (!command_name.empty() && !remaining.empty() && remaining[0] == command_name) {
            command_args = std::vector<std::string>(remaining.begin() + 1, remaining.end());
        } else {
            command_args = remaining;
        }
    }

    return RouteCommand(command_name, command_args, globals, config);
}

void CommandRouter::PrintHelp() {
    cli::Banner();

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
    for (const auto& cmd : commands_) {
        std::cout << "  " << cmd->Name();
        size_t pad = 16 - cmd->Name().length();
        std::cout << std::string(pad, ' ') << cmd->Description() << "\n";
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

void CommandRouter::PrintVersion() {
    std::cout << "Sirius v" << kSiriusVersion << "\n";
    std::cout << "General Relativistic Ray Tracer\n";
    std::cout << "\n";
    std::cout << "Build info:\n";
#ifdef SIRIUS_HAS_VULKAN_BACKEND
    std::cout << "  Vulkan backend: compiled in\n";
#else
    std::cout << "  Vulkan backend: not compiled\n";
#endif
    std::cout << "  C++ standard: C++26\n";
}

std::vector<std::string> CommandRouter::ParseGlobalOptions(const std::vector<std::string>& args,
                                                           GlobalOptions& globals) {
    // Global flags are position-independent: no command defines a flag that
    // collides with this set, so a global recognised after a command name is
    // unambiguous.
    std::vector<std::string> remaining;

    for (size_t i = 0; i < args.size(); ++i) {
        const std::string& arg = args[i];

        if (arg == "-v" || arg == "--verbose") {
            globals.verbose = true;
        } else if (arg == "--json") {
            globals.json_output = true;
        } else if (arg == "--no-color") {
            globals.no_color = true;
        } else if (arg == "--config" && i + 1 < args.size()) {
            globals.config_path = args[++i];
        } else if (arg == "--help") {
            // Only --help; -h belongs to commands (render uses it for height).
            globals.show_help = true;
        } else if (arg == "--version") {
            globals.show_version = true;
        } else {
            remaining.push_back(arg);
        }
    }

    return remaining;
}

std::string CommandRouter::ExtractCommand(const std::vector<std::string>& args) {
    if (args.empty()) {
        return "";
    }

    const std::string& first = args[0];

    for (const auto& cmd : commands_) {
        if (first == cmd->Name()) {
            return first;
        }
    }

    if (first == "help") {
        return "help";
    }

    return "";
}

bool CommandRouter::IsLegacySyntax(const std::vector<std::string>& args) {
    if (args.empty()) {
        return false;
    }

    const std::string& first = args[0];

    if (first.substr(0, 2) == "--" || first.substr(0, 1) == "-") {
        for (const auto& cmd : commands_) {
            if (first == cmd->Name()) {
                return false;
            }
        }
        return true;
    }

    return false;
}

int CommandRouter::RouteCommand(const std::string& command_name,
                                const std::vector<std::string>& args, const GlobalOptions& globals,
                                SiriusConfig& config) {
    if (command_name == "help" || command_name.empty()) {
        PrintHelp();
        return 0;
    }

    for (auto& cmd : commands_) {
        if (cmd->Name() == command_name) {
            // --help (but not -h, which some commands use for height).
            for (const auto& arg : args) {
                if (arg == "--help") {
                    std::cout << cmd->Usage() << std::endl;
                    return 0;
                }
            }

            return cmd->Execute(args, globals, config);
        }
    }

    cli::Error("Unknown command: " + command_name);
    std::cout << "\nRun 'sirius help' for available commands.\n";
    return 1;
}

}  // namespace sirius::app
