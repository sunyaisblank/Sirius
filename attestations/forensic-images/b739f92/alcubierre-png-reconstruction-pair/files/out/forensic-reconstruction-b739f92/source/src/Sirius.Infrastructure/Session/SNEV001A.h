// SNEV001A.h - Session Event Enumeration
// Component ID: SNEV001A (Session/Event)
//
// Defines all valid events/inputs for the render session FSM.

#pragma once

#include <cstdint>
#include <string_view>

namespace Sirius {

//==============================================================================
// Session Events
//==============================================================================
enum class SessionEvent : uint8_t {
    Start,              ///< User initiates render
    Ready,              ///< Initialisation complete
    TileAvailable,      ///< Tile ready for processing
    TileComplete,       ///< Tile finished rendering
    AllTilesComplete,   ///< All tiles processed
    OutputWritten,      ///< Output file written
    Pause,              ///< User requests pause
    Resume,             ///< User resumes from pause
    Cancel,             ///< User cancels render
    Error               ///< Error occurred
};

//==============================================================================
// Event Name (for logging)
//==============================================================================
constexpr std::string_view eventName(SessionEvent event) {
    switch (event) {
        case SessionEvent::Start:            return "Start";
        case SessionEvent::Ready:            return "Ready";
        case SessionEvent::TileAvailable:    return "TileAvailable";
        case SessionEvent::TileComplete:     return "TileComplete";
        case SessionEvent::AllTilesComplete: return "AllTilesComplete";
        case SessionEvent::OutputWritten:    return "OutputWritten";
        case SessionEvent::Pause:            return "Pause";
        case SessionEvent::Resume:           return "Resume";
        case SessionEvent::Cancel:           return "Cancel";
        case SessionEvent::Error:            return "Error";
    }
    return "Unknown";
}

} // namespace Sirius
