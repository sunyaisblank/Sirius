// SRER001A.h - Error Accumulator
// Component ID: SRER001A
// Thread-safe accumulation of render errors, warnings, and fatal conditions
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_SRER001A_H
#define SIRIUS_RENDER_SRER001A_H

#include <atomic>
#include <exception>
#include <mutex>
#include <string>
#include <vector>

namespace sirius::render {

//==============================================================================
// ErrorSeverity
// Severity levels for render errors
//==============================================================================

enum class ErrorSeverity {
    WARNING,  // Non-blocking issue (e.g., NaN pixel clamped)
    ERROR,    // Tile-level failure (e.g., integration diverged)
    FATAL     // Session-level failure (e.g., out of memory)
};

//==============================================================================
// RenderError
// Represents a single error/warning that occurred during rendering
//==============================================================================

struct RenderError {
    ErrorSeverity severity;
    std::string message;
    std::string context;

    // Factory methods
    static RenderError warning(const std::string& msg) {
        return RenderError{ErrorSeverity::WARNING, msg, ""};
    }

    static RenderError error(const std::string& msg, const std::string& ctx = "") {
        return RenderError{ErrorSeverity::ERROR, msg, ctx};
    }

    static RenderError fatal(const std::string& msg) {
        return RenderError{ErrorSeverity::FATAL, msg, ""};
    }

    // Create from exception
    static RenderError fromException(std::exception_ptr eptr, const std::string& ctx = "") {
        try {
            if (eptr) std::rethrow_exception(eptr);
        } catch (const std::exception& e) {
            return RenderError{ErrorSeverity::ERROR, e.what(), ctx};
        } catch (...) {
            return RenderError{ErrorSeverity::ERROR, "Non-standard exception", ctx};
        }
        return RenderError{ErrorSeverity::ERROR, "Unknown exception", ctx};
    }

    // Severity checks
    bool isWarning() const { return severity == ErrorSeverity::WARNING; }
    bool isError() const { return severity == ErrorSeverity::ERROR; }
    bool isFatal() const { return severity == ErrorSeverity::FATAL; }
};

//==============================================================================
// ErrorAccumulator
// Thread-safe collection of render errors
//==============================================================================

class ErrorAccumulator {
public:
    ErrorAccumulator() : m_hasFatal(false) {}

    // Non-copyable
    ErrorAccumulator(const ErrorAccumulator&) = delete;
    ErrorAccumulator& operator=(const ErrorAccumulator&) = delete;

    //--------------------------------------------------------------------------
    // Add Errors (Thread-Safe)
    //--------------------------------------------------------------------------

    void add(const RenderError& error) {
        if (error.isFatal()) {
            m_hasFatal.store(true, std::memory_order_release);
        }

        std::lock_guard<std::mutex> lock(m_mutex);
        m_errors.push_back(error);
    }

    //--------------------------------------------------------------------------
    // Query State
    //--------------------------------------------------------------------------

    /// Check if any errors or warnings have been recorded
    bool hasErrors() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return !m_errors.empty();
    }

    /// Check if a fatal error has occurred (lock-free check)
    bool hasFatalError() const {
        return m_hasFatal.load(std::memory_order_acquire);
    }

    /// Count warnings
    size_t warningCount() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return countBySeverity(ErrorSeverity::WARNING);
    }

    /// Count non-fatal errors
    size_t errorCount() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return countBySeverity(ErrorSeverity::ERROR);
    }

    /// Get all errors (copies)
    std::vector<RenderError> all() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_errors;
    }

    //--------------------------------------------------------------------------
    // Control
    //--------------------------------------------------------------------------

    /// Clear all accumulated errors
    void clear() {
        m_hasFatal.store(false, std::memory_order_release);
        std::lock_guard<std::mutex> lock(m_mutex);
        m_errors.clear();
    }

private:
    size_t countBySeverity(ErrorSeverity sev) const {
        size_t count = 0;
        for (const auto& e : m_errors) {
            if (e.severity == sev) ++count;
        }
        return count;
    }

    mutable std::mutex m_mutex;
    std::vector<RenderError> m_errors;
    std::atomic<bool> m_hasFatal;
};

} // namespace sirius::render

#endif // SIRIUS_RENDER_SRER001A_H
