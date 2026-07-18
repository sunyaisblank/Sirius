#pragma once

// Terminal output helpers: coloured status lines, an ASCII key-value table, a
// single-frame FTXUI progress bar, and the config view. Ported from CRCL005A.h.
// System/backend enumeration moved to info_command (it needs the backend seam);
// this header stays presentation-only.

#include <string>
#include <vector>

// The config view takes the render-layer config tree by const reference; a
// forward declaration keeps the JSON codec out of this header.
namespace sirius::render {
struct SiriusConfig;
}

namespace sirius::app::cli_output {

// Enable or disable ANSI colour (the --no-color global flag drives this).
void SetColorEnabled(bool enabled);

// --- Coloured status lines --------------------------------------------------
void Success(const std::string& message);
void Error(const std::string& message);
void Warning(const std::string& message);
void Info(const std::string& message);

// Horizontal rule, optionally captioned.
void Rule(const std::string& title = "");

// The Sirius banner.
void Banner();

// --- Progress display -------------------------------------------------------
struct ProgressState {
    float progress = 0.0f;  // 0.0 to 1.0.
    int tiles_completed = 0;
    int tiles_total = 0;
    double elapsed_seconds = 0.0;
    double eta_seconds = 0.0;
    bool is_complete = false;
};

// Print a progress bar, overwriting the current terminal line.
void PrintProgress(const ProgressState& state);

// Clear the progress line.
void ClearProgress();

// --- Key-value table --------------------------------------------------------
struct TableRow {
    std::string key;
    std::string value;
    bool is_header = false;
};

// Print a bordered key-value table.
void PrintTable(const std::string& title, const std::vector<TableRow>& rows);

// --- Configuration view -----------------------------------------------------
void PrintConfig(const render::SiriusConfig& config);

// --- JSON passthrough (for the --json flag) ---------------------------------
void PrintJson(const std::string& json);

}  // namespace sirius::app::cli_output
