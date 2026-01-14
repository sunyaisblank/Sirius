// SNST001A.h - Session State Enumeration
// Component ID: SNST001A (Session/State)
//
// Defines all valid states for the render session FSM.

#pragma once

#include <cstdint>
#include <string_view>

namespace Sirius {

//==============================================================================
// Session States
//==============================================================================
enum class SessionState : uint8_t {
    Idle,           ///< Awaiting configuration
    Initialising,   ///< Loading resources, creating tiles
    Scheduling,     ///< Assigning next tile to worker
    Rendering,      ///< Active tile processing
    Paused,         ///< User-requested pause
    Completing,     ///< All tiles done, writing output
    Complete,       ///< Success (terminal)
    Failed,         ///< Error occurred (terminal)
    Cancelled       ///< User cancelled (terminal)
};

//==============================================================================
// State Queries
//==============================================================================
constexpr bool isTerminal(SessionState state) {
    return state == SessionState::Complete ||
           state == SessionState::Failed ||
           state == SessionState::Cancelled;
}

constexpr bool isActive(SessionState state) {
    return state == SessionState::Scheduling ||
           state == SessionState::Rendering ||
           state == SessionState::Paused;
}

//==============================================================================
// State Name (for logging)
//==============================================================================
constexpr std::string_view stateName(SessionState state) {
    switch (state) {
        case SessionState::Idle:         return "Idle";
        case SessionState::Initialising: return "Initialising";
        case SessionState::Scheduling:   return "Scheduling";
        case SessionState::Rendering:    return "Rendering";
        case SessionState::Paused:       return "Paused";
        case SessionState::Completing:   return "Completing";
        case SessionState::Complete:     return "Complete";
        case SessionState::Failed:       return "Failed";
        case SessionState::Cancelled:    return "Cancelled";
    }
    return "Unknown";
}

} // namespace Sirius
