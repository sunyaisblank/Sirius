// Terminal output formatters. Ported from CRCL005A.cpp.

#include "sirius/app/cli/cli_output.h"

#include "sirius/app/config/config_schema.h"

#include <ftxui/dom/elements.hpp>
#include <ftxui/screen/screen.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace sirius::app::cli_output {

using namespace ftxui;

namespace {

// Global colour state; the --no-color flag clears it at startup.
bool g_color_enabled = true;

// Single-frame FTXUI progress bar element. File-local: PrintProgress is the
// only consumer.
ftxui::Element RenderProgressBar(const ProgressState& state) {
    const int bar_width = 30;
    int filled = static_cast<int>(state.progress * bar_width);

    std::string bar;
    for (int i = 0; i < bar_width; ++i) {
        if (i < filled) {
            bar += "█";
        } else {
            bar += "░";
        }
    }

    std::ostringstream pct;
    pct << std::fixed << std::setprecision(0) << (state.progress * 100) << "%";

    std::ostringstream tiles;
    tiles << "Tile " << state.tiles_completed << "/" << state.tiles_total;

    std::ostringstream eta;
    if (state.eta_seconds > 0 && !state.is_complete) {
        int mins = static_cast<int>(state.eta_seconds) / 60;
        int secs = static_cast<int>(state.eta_seconds) % 60;
        eta << "ETA: " << mins << "m " << secs << "s";
    } else if (state.is_complete) {
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

}  // namespace

void SetColorEnabled(bool enabled) { g_color_enabled = enabled; }

// --- Coloured status lines --------------------------------------------------

void Success(const std::string& message) {
    if (g_color_enabled) {
        std::cout << "\033[32m✓\033[0m " << message << std::endl;
    } else {
        std::cout << "[OK] " << message << std::endl;
    }
}

void Error(const std::string& message) {
    if (g_color_enabled) {
        std::cerr << "\033[31m✗\033[0m " << message << std::endl;
    } else {
        std::cerr << "[ERROR] " << message << std::endl;
    }
}

void Warning(const std::string& message) {
    if (g_color_enabled) {
        std::cout << "\033[33m!\033[0m " << message << std::endl;
    } else {
        std::cout << "[WARN] " << message << std::endl;
    }
}

void Info(const std::string& message) {
    if (g_color_enabled) {
        std::cout << "\033[34mℹ\033[0m " << message << std::endl;
    } else {
        std::cout << "[INFO] " << message << std::endl;
    }
}

void Rule(const std::string& title) {
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

void Banner() {
    if (g_color_enabled) {
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
    if (g_color_enabled) {
        std::cout << "\033[0m";
    }
    std::cout << "  General Relativistic Ray Tracer" << std::endl;
    std::cout << std::endl;
}

// --- Progress display -------------------------------------------------------

void PrintProgress(const ProgressState& state) {
    Element element = RenderProgressBar(state);
    Screen screen = Screen::Create(Dimension::Full(), Dimension::Fixed(1));
    Render(screen, element);

    std::cout << "\r" << screen.ToString() << std::flush;
}

void ClearProgress() { std::cout << "\r" << std::string(80, ' ') << "\r" << std::flush; }

// --- Key-value table --------------------------------------------------------

void PrintTable(const std::string& title, const std::vector<TableRow>& rows) {
    std::size_t max_key_width = title.length();
    std::size_t max_value_width = 0;
    for (const auto& row : rows) {
        if (!row.is_header) {
            max_key_width = std::max(max_key_width, row.key.length());
            max_value_width = std::max(max_value_width, row.value.length());
        }
    }

    std::size_t total_width = max_key_width + max_value_width + 7;

    std::cout << "+" << std::string(total_width - 2, '-') << "+" << std::endl;
    std::size_t title_pad = (total_width - 2 - title.length()) / 2;
    std::cout << "|" << std::string(title_pad, ' ') << title
              << std::string(total_width - 2 - title_pad - title.length(), ' ') << "|" << std::endl;
    std::cout << "+" << std::string(total_width - 2, '-') << "+" << std::endl;

    for (const auto& row : rows) {
        if (row.is_header) {
            std::cout << "+" << std::string(total_width - 2, '-') << "+" << std::endl;
            std::size_t pad = (total_width - 2 - row.key.length()) / 2;
            std::cout << "|" << std::string(pad, ' ') << row.key
                      << std::string(total_width - 2 - pad - row.key.length(), ' ') << "|"
                      << std::endl;
            std::cout << "+" << std::string(total_width - 2, '-') << "+" << std::endl;
        } else {
            std::string padded_key = row.key;
            padded_key.resize(max_key_width, ' ');
            std::string padded_value = row.value;
            padded_value.resize(max_value_width, ' ');
            std::cout << "| " << padded_key << " | " << padded_value << " |" << std::endl;
        }
    }

    std::cout << "+" << std::string(total_width - 2, '-') << "+" << std::endl;
}

// --- Configuration view -----------------------------------------------------

void PrintConfig(const SiriusConfig& config) {
    std::vector<TableRow> rows;

    rows.push_back({"Render Settings", "", true});
    rows.push_back(
        {"Resolution",
         std::to_string(config.render.width) + " x " + std::to_string(config.render.height),
         false});
    rows.push_back({"Samples/Pixel", std::to_string(config.render.samples_per_pixel), false});
    rows.push_back({"Tile Size", std::to_string(config.render.tile_size) + " px", false});
    rows.push_back({"Output", config.render.output_path, false});

    rows.push_back({"Metric Settings", "", true});
    rows.push_back({"Metric", config.metric.name, false});
    rows.push_back({"Mass (M)", std::to_string(config.metric.mass), false});
    {
        std::ostringstream spin_str;
        spin_str << std::fixed << std::setprecision(3) << config.metric.spin;
        rows.push_back({"Spin (a)", spin_str.str(), false});
    }

    rows.push_back({"Observer Settings", "", true});
    {
        std::ostringstream dist;
        dist << std::fixed << std::setprecision(1) << config.observer.distance << " M";
        rows.push_back({"Distance", dist.str(), false});
    }
    {
        std::ostringstream incl;
        incl << std::fixed << std::setprecision(1) << config.observer.inclination << "°";
        rows.push_back({"Inclination", incl.str(), false});
    }
    {
        std::ostringstream fov;
        fov << std::fixed << std::setprecision(1) << config.observer.fov << "°";
        rows.push_back({"FOV", fov.str(), false});
    }

    rows.push_back({"Post-Processing", "", true});
    rows.push_back({"Bloom", config.postprocess.enable_bloom ? "Enabled" : "Disabled", false});
    {
        std::ostringstream exp;
        exp << std::fixed << std::setprecision(2) << config.postprocess.exposure;
        rows.push_back({"Exposure", exp.str(), false});
    }

    PrintTable("Configuration", rows);
}

// --- JSON passthrough -------------------------------------------------------

void PrintJson(const std::string& json) { std::cout << json << std::endl; }

}  // namespace sirius::app::cli_output
