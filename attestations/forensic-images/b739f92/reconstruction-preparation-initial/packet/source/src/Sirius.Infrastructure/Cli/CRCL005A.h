// CRCL005A.h - CLI Output Formatters
// Component ID: CRCL005A (Cli/Output)
//
// FTXUI-based rich terminal output: progress bars, tables, colored text.
// Provides consistent formatting across all CLI commands.

#pragma once

#include <ftxui/dom/elements.hpp>
#include <ftxui/screen/screen.hpp>
#include <ftxui/screen/terminal.hpp>

#include <chrono>
#include <functional>
#include <string>
#include <vector>
#include <atomic>
#include <mutex>

namespace Sirius::Cli::Output {

/// @brief Global color state (can be disabled with --no-color)
extern bool g_colorEnabled;

/// @brief Set whether colors are enabled
void setColorEnabled(bool enabled);

/// @brief Check if colors are enabled
bool isColorEnabled();

// =============================================================================
// Colored Output Utilities
// =============================================================================

/// @brief Print success message (green checkmark)
void success(const std::string& message);

/// @brief Print error message (red X)
void error(const std::string& message);

/// @brief Print warning message (yellow !)
void warning(const std::string& message);

/// @brief Print info message (blue i)
void info(const std::string& message);

/// @brief Print a horizontal rule/divider
void rule(const std::string& title = "");

/// @brief Print the Sirius banner
void banner();

// =============================================================================
// Progress Display
// =============================================================================

/// @brief Render progress state (snapshot for display)
struct ProgressState {
    float progress = 0.0f;        ///< 0.0 to 1.0
    int tilesCompleted = 0;
    int tilesTotal = 0;
    double elapsedSeconds = 0.0;
    double etaSeconds = 0.0;
    double raysPerSecond = 0.0;
    bool isComplete = false;
    bool hasFailed = false;
    std::string errorMessage;
};

/// @brief Render a progress bar (non-interactive, single frame)
/// @param state Current progress state
/// @return Rendered FTXUI element
ftxui::Element renderProgressBar(const ProgressState& state);

/// @brief Print progress bar to terminal (overwrites current line)
void printProgress(const ProgressState& state);

/// @brief Clear progress line
void clearProgress();

// =============================================================================
// Table Formatting
// =============================================================================

/// @brief Table row (key-value pairs)
struct TableRow {
    std::string key;
    std::string value;
    bool isHeader = false;
};

/// @brief Render a key-value table
/// @param title Table title
/// @param rows Table rows
/// @return Rendered FTXUI element
ftxui::Element renderTable(const std::string& title,
                           const std::vector<TableRow>& rows);

/// @brief Print a key-value table to terminal
void printTable(const std::string& title,
                const std::vector<TableRow>& rows);

/// @brief Render system info table
ftxui::Element renderSystemInfo();

/// @brief Print system info table
void printSystemInfo();

// =============================================================================
// Configuration Display
// =============================================================================

} // namespace Sirius::Cli::Output

// Forward declaration
namespace Sirius::Configuration {
    struct SiriusConfig;
}

namespace Sirius::Cli::Output {

/// @brief Render configuration as table
ftxui::Element renderConfig(const Sirius::Configuration::SiriusConfig& config);

/// @brief Print configuration to terminal
void printConfig(const Sirius::Configuration::SiriusConfig& config);

// =============================================================================
// JSON Output (for --json flag)
// =============================================================================

/// @brief Print JSON to stdout
void printJson(const std::string& json);

// =============================================================================
// Spinner (for long operations)
// =============================================================================

/// @brief Simple spinner for blocking operations
class Spinner {
public:
    Spinner(const std::string& message);
    ~Spinner();

    /// @brief Start the spinner (non-blocking, runs in background)
    void start();

    /// @brief Stop the spinner with success
    void stop(bool success = true);

    /// @brief Update spinner message
    void setMessage(const std::string& message);

private:
    std::string m_message;
    std::atomic<bool> m_running{false};
    std::mutex m_mutex;
    int m_frameIndex{0};
};

} // namespace Sirius::Cli::Output
