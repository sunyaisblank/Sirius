// SRCT001A.h - Cancellation Token
// Component ID: SRCT001A
// Provides thread-safe cancellation signaling for long-running render operations
//
// Governance: docs/specification.md

#ifndef SIRIUS_RENDER_SRCT001A_H
#define SIRIUS_RENDER_SRCT001A_H

#include <atomic>
#include <stdexcept>
#include <string>

namespace sirius::render {

//==============================================================================
// CancelledException
// Thrown when a cancelled operation is checked via throwIfCancelled()
//==============================================================================

class CancelledException : public std::runtime_error {
public:
    CancelledException()
        : std::runtime_error("Render operation was cancelled by user request") {}

    explicit CancelledException(const std::string& message)
        : std::runtime_error(message) {}
};

//==============================================================================
// CancellationToken
// Thread-safe cancellation signaling mechanism
//==============================================================================

class CancellationToken {
public:
    CancellationToken() : m_cancelled(false) {}

    // Non-copyable, moveable
    CancellationToken(const CancellationToken&) = delete;
    CancellationToken& operator=(const CancellationToken&) = delete;
    CancellationToken(CancellationToken&&) = default;
    CancellationToken& operator=(CancellationToken&&) = default;

    //--------------------------------------------------------------------------
    // Query State
    //--------------------------------------------------------------------------

    /// Check if cancellation has been requested
    /// Thread-safe: can be called from any thread
    bool isCancelled() const {
        return m_cancelled.load(std::memory_order_acquire);
    }

    /// Throw CancelledException if cancelled
    /// Use in hot paths that should abort on cancellation
    void throwIfCancelled() const {
        if (isCancelled()) {
            throw CancelledException();
        }
    }

    //--------------------------------------------------------------------------
    // Control
    //--------------------------------------------------------------------------

    /// Request cancellation
    /// Thread-safe: can be called from any thread (e.g., UI thread)
    void request() {
        m_cancelled.store(true, std::memory_order_release);
    }

    /// Reset cancellation state
    /// Should only be called when no render operation is active
    void reset() {
        m_cancelled.store(false, std::memory_order_release);
    }

private:
    std::atomic<bool> m_cancelled;
};

} // namespace sirius::render

#endif // SIRIUS_RENDER_SRCT001A_H
