// Program entry point. Routes every command line through the CLI command
// router; uncaught exceptions become a non-zero exit with a diagnostic. Ported
// from CREP001A.cpp.

#include "sirius/app/cli/command_router.h"

#include <exception>
#include <iostream>

int main(int argc, char* argv[]) {
    try {
        sirius::app::CommandRouter router;
        return router.Run(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown fatal error occurred." << std::endl;
        return 1;
    }
}
