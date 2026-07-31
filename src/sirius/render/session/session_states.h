#pragma once

// Render-session FSM states and their queries. Ported from SNST001A.h.

#include <cstdint>
#include <string_view>

namespace sirius::render {

// All valid states of the render session.
enum class SessionState : uint8_t {
    Idle,          // Awaiting configuration.
    Initialising,  // Loading resources, creating tiles.
    Scheduling,    // Assigning the next tile to a worker.
    Rendering,     // Active tile processing.
    Completing,    // All tiles done, writing output.
    Complete,      // Success (terminal).
    Failed,        // Error occurred (terminal).
    Cancelled      // User cancelled (terminal).
};

// Whether the session has reached a terminal state.
constexpr bool IsTerminal(SessionState state) {
    return state == SessionState::Complete || state == SessionState::Failed ||
           state == SessionState::Cancelled;
}

// Whether the session is actively running (schedulable or paused).
constexpr bool IsActive(SessionState state) {
    return state == SessionState::Initialising || state == SessionState::Scheduling ||
           state == SessionState::Rendering || state == SessionState::Completing;
}

// Human-readable state name (for logging).
constexpr std::string_view StateName(SessionState state) {
    switch (state) {
        case SessionState::Idle:
            return "Idle";
        case SessionState::Initialising:
            return "Initialising";
        case SessionState::Scheduling:
            return "Scheduling";
        case SessionState::Rendering:
            return "Rendering";
        case SessionState::Completing:
            return "Completing";
        case SessionState::Complete:
            return "Complete";
        case SessionState::Failed:
            return "Failed";
        case SessionState::Cancelled:
            return "Cancelled";
    }
    return "Unknown";
}

}  // namespace sirius::render
