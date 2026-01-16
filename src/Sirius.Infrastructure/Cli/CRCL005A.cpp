// CRCL005A.cpp - CLI Output Formatters
// Component ID: CRCL005A (Cli/Output)

#include "CRCL005A.h"
#include "CRCF003A.h"
#include "CRPF001A.h"

#include <ftxui/dom/elements.hpp>
#include <ftxui/screen/screen.hpp>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

#ifdef SIRIUS_HAS_OPTIX
#include <cuda_runtime.h>
#endif

namespace Sirius::Cli::Output {

using namespace ftxui;

// Global color state
bool g_colorEnabled = true;

void setColorEnabled(bool enabled) {
    g_colorEnabled = enabled;
}

bool isColorEnabled() {
    return g_colorEnabled;
}

// =============================================================================
// Colored Output Utilities
// =============================================================================

void success(const std::string& message) {
    if (g_colorEnabled) {
        std::cout << "\033[32m✓\033[0m " << message << std::endl;
    } else {
        std::cout << "[OK] " << message << std::endl;
    }
}

void error(const std::string& message) {
    if (g_colorEnabled) {
        std::cerr << "\033[31m✗\033[0m " << message << std::endl;
    } else {
        std::cerr << "[ERROR] " << message << std::endl;
    }
}

void warning(const std::string& message) {
    if (g_colorEnabled) {
        std::cout << "\033[33m!\033[0m " << message << std::endl;
    } else {
        std::cout << "[WARN] " << message << std::endl;
    }
}

void info(const std::string& message) {
    if (g_colorEnabled) {
        std::cout << "\033[34mℹ\033[0m " << message << std::endl;
    } else {
        std::cout << "[INFO] " << message << std::endl;
    }
}

void rule(const std::string& title) {
    int width = 60;
    std::string dash = "-";
    if (title.empty()) {
        for (int i = 0; i < width; ++i) std::cout << dash;
        std::cout << std::endl;
    } else {
        int padding = (width - static_cast<int>(title.length()) - 2) / 2;
        for (int i = 0; i < padding; ++i) std::cout << dash;
        std::cout << " " << title << " ";
        for (int i = 0; i < padding; ++i) std::cout << dash;
        std::cout << std::endl;
    }
}

void banner() {
    if (g_colorEnabled) {
        std::cout << "\033[36m";
    }
    std::cout << R"(
   _____ _      _
  / ____(_)    (_)
 | (___  _ _ __ _ _   _ ___
  \___ \| | '__| | | | / __|
  ____) | | |  | | |_| \__ \
 |_____/|_|_|  |_|\__,_|___/
)" << std::endl;
    if (g_colorEnabled) {
        std::cout << "\033[0m";
    }
    std::cout << "  General Relativistic Ray Tracer" << std::endl;
    std::cout << std::endl;
}

// =============================================================================
// Progress Display
// =============================================================================

ftxui::Element renderProgressBar(const ProgressState& state) {
    // Build progress bar string
    const int barWidth = 30;
    int filled = static_cast<int>(state.progress * barWidth);

    std::string bar;
    for (int i = 0; i < barWidth; ++i) {
        if (i < filled) {
            bar += "█";
        } else {
            bar += "░";
        }
    }

    // Format percentage
    std::ostringstream pct;
    pct << std::fixed << std::setprecision(0) << (state.progress * 100) << "%";

    // Format tiles
    std::ostringstream tiles;
    tiles << "Tile " << state.tilesCompleted << "/" << state.tilesTotal;

    // Format ETA
    std::ostringstream eta;
    if (state.etaSeconds > 0 && !state.isComplete) {
        int mins = static_cast<int>(state.etaSeconds) / 60;
        int secs = static_cast<int>(state.etaSeconds) % 60;
        eta << "ETA: " << mins << "m " << secs << "s";
    } else if (state.isComplete) {
        eta << "Complete!";
    }

    return hbox({
        text("["),
        text(bar) | color(Color::Green),
        text("] "),
        text(pct.str()) | bold,
        text(" | "),
        text(tiles.str()),
        text(" | "),
        text(eta.str()) | dim,
    });
}

void printProgress(const ProgressState& state) {
    Element element = renderProgressBar(state);
    Screen screen = Screen::Create(Dimension::Full(), Dimension::Fixed(1));
    Render(screen, element);

    // Print with carriage return to overwrite line
    std::cout << "\r" << screen.ToString() << std::flush;
}

void clearProgress() {
    std::cout << "\r" << std::string(80, ' ') << "\r" << std::flush;
}

// =============================================================================
// Table Formatting
// =============================================================================

ftxui::Element renderTable(const std::string& title,
                           const std::vector<TableRow>& rows) {
    Elements tableRows;

    // Find max key width for alignment
    size_t maxKeyWidth = title.length();
    for (const auto& row : rows) {
        if (!row.isHeader) {
            maxKeyWidth = std::max(maxKeyWidth, row.key.length());
        }
    }

    for (const auto& row : rows) {
        if (row.isHeader) {
            tableRows.push_back(
                text(row.key) | bold | color(Color::Blue) | center
            );
        } else {
            std::string paddedKey = row.key;
            paddedKey.resize(maxKeyWidth, ' ');

            Elements rowElements;
            rowElements.push_back(text("| "));
            rowElements.push_back(text(paddedKey) | dim);
            rowElements.push_back(text(" | "));
            rowElements.push_back(text(row.value));
            rowElements.push_back(text(" |"));

            tableRows.push_back(hbox(std::move(rowElements)));
        }
    }

    Elements result;
    result.push_back(text(title) | bold | center);
    result.push_back(separator());
    for (auto& r : tableRows) {
        result.push_back(std::move(r));
    }

    return vbox(std::move(result)) | border;
}

void printTable(const std::string& title, const std::vector<TableRow>& rows) {
    // Find max widths
    size_t maxKeyWidth = title.length();
    size_t maxValueWidth = 0;
    for (const auto& row : rows) {
        if (!row.isHeader) {
            maxKeyWidth = std::max(maxKeyWidth, row.key.length());
            maxValueWidth = std::max(maxValueWidth, row.value.length());
        }
    }

    size_t totalWidth = maxKeyWidth + maxValueWidth + 7;

    // Print header
    std::cout << "+" << std::string(totalWidth - 2, '-') << "+" << std::endl;
    size_t titlePad = (totalWidth - 2 - title.length()) / 2;
    std::cout << "|" << std::string(titlePad, ' ') << title
              << std::string(totalWidth - 2 - titlePad - title.length(), ' ')
              << "|" << std::endl;
    std::cout << "+" << std::string(totalWidth - 2, '-') << "+" << std::endl;

    // Print rows
    for (const auto& row : rows) {
        if (row.isHeader) {
            std::cout << "+" << std::string(totalWidth - 2, '-') << "+" << std::endl;
            size_t pad = (totalWidth - 2 - row.key.length()) / 2;
            std::cout << "|" << std::string(pad, ' ') << row.key
                      << std::string(totalWidth - 2 - pad - row.key.length(), ' ')
                      << "|" << std::endl;
            std::cout << "+" << std::string(totalWidth - 2, '-') << "+" << std::endl;
        } else {
            std::string paddedKey = row.key;
            paddedKey.resize(maxKeyWidth, ' ');
            std::string paddedValue = row.value;
            paddedValue.resize(maxValueWidth, ' ');
            std::cout << "| " << paddedKey << " | " << paddedValue << " |" << std::endl;
        }
    }

    std::cout << "+" << std::string(totalWidth - 2, '-') << "+" << std::endl;
}

void printSystemInfo() {
    std::vector<TableRow> rows;

    // Platform info
    rows.push_back({"Platform", Platform::PathResolver::platformName()});

    if (Platform::PathResolver::isWSL2()) {
        rows.push_back({"Environment", "WSL2"});
    }

#ifdef SIRIUS_HAS_OPTIX
    // CUDA info
    int cudaVersion = 0;
    cudaRuntimeGetVersion(&cudaVersion);
    int major = cudaVersion / 1000;
    int minor = (cudaVersion % 1000) / 10;
    std::ostringstream cudaVer;
    cudaVer << major << "." << minor;
    rows.push_back({"CUDA", cudaVer.str()});

    // GPU info
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, 0);
        rows.push_back({"GPU", prop.name});

        std::ostringstream vram;
        vram << (prop.totalGlobalMem / (1024 * 1024)) << " MB";
        rows.push_back({"VRAM", vram.str()});

        std::ostringstream cc;
        cc << prop.major << "." << prop.minor;
        rows.push_back({"Compute Cap.", cc.str()});
    }

    rows.push_back({"", "", true});  // Section header
    rows.push_back({"Available Backends", "", true});
    rows.push_back({"OptiX", "Available (preferred)"});
    rows.push_back({"CPU", "Fallback"});
#else
    rows.push_back({"CUDA", "Not available"});
    rows.push_back({"", "", true});
    rows.push_back({"Available Backends", "", true});
    rows.push_back({"CPU", "Available"});
#endif

    printTable("System Information", rows);
}

// =============================================================================
// Configuration Display
// =============================================================================

void printConfig(const Sirius::Configuration::SiriusConfig& config) {
    std::vector<TableRow> rows;

    // Render section
    rows.push_back({"Render Settings", "", true});
    rows.push_back({"Resolution", std::to_string(config.render.width) + " x " + std::to_string(config.render.height)});
    rows.push_back({"Samples/Pixel", std::to_string(config.render.samplesPerPixel)});
    rows.push_back({"Tile Size", std::to_string(config.render.tileSize) + " px"});
    rows.push_back({"Output", config.render.outputPath});

    // Metric section
    rows.push_back({"Metric Settings", "", true});
    rows.push_back({"Metric", config.metric.name});
    rows.push_back({"Mass (M)", std::to_string(config.metric.mass)});
    rows.push_back({"Spin (a)", std::to_string(config.metric.spin)});

    // Observer section
    rows.push_back({"Observer Settings", "", true});
    std::ostringstream dist;
    dist << std::fixed << std::setprecision(1) << config.observer.distance << " M";
    rows.push_back({"Distance", dist.str()});

    std::ostringstream incl;
    incl << std::fixed << std::setprecision(1) << config.observer.inclination << "°";
    rows.push_back({"Inclination", incl.str()});

    std::ostringstream fov;
    fov << std::fixed << std::setprecision(1) << config.observer.fov << "°";
    rows.push_back({"FOV", fov.str()});

    // Post-process section
    rows.push_back({"Post-Processing", "", true});
    rows.push_back({"Bloom", config.postprocess.enableBloom ? "Enabled" : "Disabled"});

    std::ostringstream exp;
    exp << std::fixed << std::setprecision(2) << config.postprocess.exposure;
    rows.push_back({"Exposure", exp.str()});

    printTable("Configuration", rows);
}

// =============================================================================
// JSON Output
// =============================================================================

void printJson(const std::string& json) {
    std::cout << json << std::endl;
}

// =============================================================================
// Spinner
// =============================================================================

Spinner::Spinner(const std::string& message)
    : m_message(message) {}

Spinner::~Spinner() {
    stop(true);
}

void Spinner::start() {
    m_running = true;
    // Simple inline spinner - just print the message
    std::cout << m_message << "..." << std::flush;
}

void Spinner::stop(bool success) {
    if (m_running.exchange(false)) {
        std::cout << " ";
        if (success) {
            if (g_colorEnabled) {
                std::cout << "\033[32m✓\033[0m";
            } else {
                std::cout << "[OK]";
            }
        } else {
            if (g_colorEnabled) {
                std::cout << "\033[31m✗\033[0m";
            } else {
                std::cout << "[FAIL]";
            }
        }
        std::cout << std::endl;
    }
}

void Spinner::setMessage(const std::string& message) {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_message = message;
}

} // namespace Sirius::Cli::Output
