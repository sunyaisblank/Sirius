#pragma once

// Render-session FSM events (inputs). Ported from SNEV001A.h.

#include <cstdint>
#include <string_view>

namespace sirius::render {

// All valid events driving the render session.
enum class SessionEvent : uint8_t {
    Start,             // User initiates a render.
    Ready,             // Initialisation complete.
    TileAvailable,     // Tile ready for processing.
    TileComplete,      // Tile finished rendering.
    AllTilesComplete,  // All tiles processed.
    OutputWritten,     // Output file written.
    Pause,             // User requests a pause.
    Resume,            // User resumes from pause.
    Cancel,            // User cancels the render.
    Error              // Error occurred.
};

// Human-readable event name (for logging).
constexpr std::string_view EventName(SessionEvent event) {
    switch (event) {
        case SessionEvent::Start:
            return "Start";
        case SessionEvent::Ready:
            return "Ready";
        case SessionEvent::TileAvailable:
            return "TileAvailable";
        case SessionEvent::TileComplete:
            return "TileComplete";
        case SessionEvent::AllTilesComplete:
            return "AllTilesComplete";
        case SessionEvent::OutputWritten:
            return "OutputWritten";
        case SessionEvent::Pause:
            return "Pause";
        case SessionEvent::Resume:
            return "Resume";
        case SessionEvent::Cancel:
            return "Cancel";
        case SessionEvent::Error:
            return "Error";
    }
    return "Unknown";
}

}  // namespace sirius::render
