// CREP001A.cpp - Program Entry Point
// Component ID: CREP001A (Application/Entry Point)
//
// Main entry point for Sirius CLI. Routes to CommandRouter for all commands.

#include "CRCL001A.h"

#include <iostream>
#include <exception>

int main(int argc, char* argv[]) {
    try {
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
