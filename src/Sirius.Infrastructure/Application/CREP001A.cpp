// CREP001A.cpp - Program Entry Point
// Component ID: CREP001A (Application/Entry Point)
//
// Main entry point for Sirius CLI. Routes to CommandRouter for all commands.
// Maintains backward compatibility with --interactive mode (deprecated).

#include "CRAP001A.h"
#include "CRCL001A.h"

#include <iostream>
#include <exception>
#include <string>
#include <cstring>

/// @brief Run deprecated interactive mode
/// Kept for backward compatibility only.
int runInteractiveMode() {
    std::cout << "================================================================================\n";
    std::cout << "                      SIRIUS INTERACTIVE MODE (DEPRECATED)                      \n";
    std::cout << "================================================================================\n";
    std::cout << "\nWARNING: Interactive mode is deprecated and will be removed.\n";
    std::cout << "         Please use batch rendering: sirius render [options]\n\n";

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    Application app;
    app.run();
#pragma GCC diagnostic pop
    return 0;
}

int main(int argc, char* argv[]) {
    try {
        // Check for deprecated --interactive flag
        for (int i = 1; i < argc; ++i) {
            if (std::strcmp(argv[i], "--interactive") == 0) {
                return runInteractiveMode();
            }
        }

        // Route to CLI command router
        Sirius::Cli::CommandRouter router;
        return router.run(argc, argv);
    }
    catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Unknown fatal error occurred." << std::endl;
        return 1;
    }
}
