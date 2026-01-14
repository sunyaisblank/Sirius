// SNRS001A.h - Render Session
// Component ID: SNRS001A (Session/RenderSession)
//
// FSM-based render session orchestrator. Manages tile scheduling,
// progress tracking, and display buffer updates.

#pragma once

#include "SMSM001A.h"
#include "SNST001A.h"
#include "SNEV001A.h"
#include "SNTL001A.h"
#include "SNPR001A.h"
#include "SNDP001A.h"
#include "IRRJ001A.h"
#include <thread>
#include <vector>
#include <functional>
#include <iostream>

namespace Sirius {

//==============================================================================
// Session Transition Table (Compile-time DFA)
//==============================================================================
constexpr std::array<Transition<SessionState, SessionEvent>, 14> SESSION_TRANSITIONS = {{
    // From Idle
    {SessionState::Idle, SessionEvent::Start, SessionState::Initialising},
    
    // From Initialising
    {SessionState::Initialising, SessionEvent::Ready, SessionState::Scheduling},
    {SessionState::Initialising, SessionEvent::Error, SessionState::Failed},
    
    // From Scheduling
    {SessionState::Scheduling, SessionEvent::TileAvailable, SessionState::Rendering},
    {SessionState::Scheduling, SessionEvent::AllTilesComplete, SessionState::Completing},
    
    // From Rendering
    {SessionState::Rendering, SessionEvent::TileComplete, SessionState::Scheduling},
    {SessionState::Rendering, SessionEvent::Pause, SessionState::Paused},
    {SessionState::Rendering, SessionEvent::Cancel, SessionState::Cancelled},
    {SessionState::Rendering, SessionEvent::Error, SessionState::Failed},
    
    // From Paused
    {SessionState::Paused, SessionEvent::Resume, SessionState::Rendering},
    {SessionState::Paused, SessionEvent::Cancel, SessionState::Cancelled},
    
    // From Completing
    {SessionState::Completing, SessionEvent::OutputWritten, SessionState::Complete},
    {SessionState::Completing, SessionEvent::Error, SessionState::Failed},
    
    // Cancel available from Scheduling
    {SessionState::Scheduling, SessionEvent::Cancel, SessionState::Cancelled}
}};

constexpr StateMachineConfig<SessionState, SessionEvent, 14> SESSION_CONFIG = {
    SessionState::Idle,
    SESSION_TRANSITIONS
};

//==============================================================================
// Render Session Configuration
//==============================================================================
struct SessionConfig {
    int width = 1920;
    int height = 1080;
    int tileSize = 64;
    int samplesPerPixel = 64;
    int threadCount = 0;    ///< 0 = auto-detect
    std::string outputPath = "render.ppm";
    
    // Black hole parameters
    std::string metricName = "Schwarzschild";
    double blackHoleMass = 1.0;
    double blackHoleSpin = 0.0;
    double observerDistance = 50.0;
    double observerInclination = 1.5708;  // 90 degrees
    float cameraFOV = 60.0f;
    
    // Post-processing
    bool enableBloom = true;
    float bloomIntensity = 0.3f;
    float exposure = 1.0f;
};

//==============================================================================
// Render Session
//==============================================================================
class RenderSession {
public:
    using FSM = StateMachine<SessionState, SessionEvent, 14>;
    using CompletionCallback = std::function<void(SessionState finalState, const std::string& message)>;
    
    RenderSession() : m_FSM(SESSION_CONFIG) {
        setupActions();
    }
    
    /// @brief Configure session (must be in Idle state)
    bool configure(const SessionConfig& config) {
        if (m_FSM.getState() != SessionState::Idle) {
            return false;
        }
        m_Config = config;
        return true;
    }
    
    /// @brief Start the render session
    bool start() {
        return m_FSM.process(SessionEvent::Start);
    }
    
    /// @brief Pause rendering
    bool pause() {
        return m_FSM.process(SessionEvent::Pause);
    }
    
    /// @brief Resume from pause
    bool resume() {
        return m_FSM.process(SessionEvent::Resume);
    }
    
    /// @brief Cancel rendering
    bool cancel() {
        m_Progress.getCancellationToken().cancel();
        return m_FSM.process(SessionEvent::Cancel);
    }
    
    /// @brief Wait for session to complete (blocking)
    void waitForCompletion() {
        if (m_RenderThread.joinable()) {
            m_RenderThread.join();
        }
    }
    
    /// @brief Execute synchronously (blocking)
    SessionState execute() {
        if (!start()) {
            return SessionState::Failed;
        }
        waitForCompletion();
        return m_FSM.getState();
    }
    
    // State queries
    SessionState getState() const { return m_FSM.getState(); }
    bool isComplete() const { return isTerminal(m_FSM.getState()); }
    bool isRunning() const { return isActive(m_FSM.getState()); }
    
    // Progress access
    float getProgress() const { return m_Progress.getProgress(); }
    double getETA() const { return m_Progress.getETA(); }
    const TileScheduler& getTileScheduler() const { return m_Tiles; }
    DisplayBuffer& getDisplayBuffer() { return m_Display; }
    
    // Callbacks
    void setCompletionCallback(CompletionCallback cb) { m_CompletionCallback = cb; }
    void setProgressCallback(ProgressCallback cb) { m_Progress.setCallback(cb); }
    
private:
    void setupActions() {
        m_FSM.setEntryAction([this](SessionState state) {
            onEnterState(state);
        });
        
        m_FSM.setTransitionAction([this](SessionState from, SessionEvent event, SessionState to) {
            std::cout << "[Session] " << stateName(from) << " --" 
                      << eventName(event) << "--> " << stateName(to) << std::endl;
        });
    }
    
    void onEnterState(SessionState state) {
        switch (state) {
            case SessionState::Initialising:
                initialise();
                break;
            case SessionState::Scheduling:
                scheduleNextTile();
                break;
            case SessionState::Rendering:
                // Tile is being processed
                break;
            case SessionState::Completing:
                writeOutput();
                break;
            case SessionState::Complete:
            case SessionState::Failed:
            case SessionState::Cancelled:
                onSessionEnd(state);
                break;
            default:
                break;
        }
    }
    
    void initialise();
    void scheduleNextTile();
    void renderTile(Tile* tile);
    void writeOutput();
    void onSessionEnd(SessionState state);
    
    FSM m_FSM;
    SessionConfig m_Config;
    TileScheduler m_Tiles;
    ProgressTracker m_Progress;
    DisplayBuffer m_Display;
    std::thread m_RenderThread;
    std::string m_ErrorMessage;
    CompletionCallback m_CompletionCallback;
};

} // namespace Sirius
