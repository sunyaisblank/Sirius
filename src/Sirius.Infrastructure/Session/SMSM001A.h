// SMSM001A.h - Generic Finite State Machine
// Component ID: SMSM001A (StateMachine/Machine)
//
// Deterministic Mealy machine template with compile-time transition table.
// Supports entry/exit actions per state and transition actions.
//
// Reference: Finite-state automaton (DFA) with Mealy output model

#pragma once

#include <array>
#include <functional>
#include <optional>
#include <mutex>
#include <atomic>

namespace Sirius {

//==============================================================================
// Transition Definition
//==============================================================================
template<typename State, typename Event>
struct Transition {
    State from;
    Event event;
    State to;
    
    constexpr Transition(State f, Event e, State t) : from(f), event(e), to(t) {}
};

//==============================================================================
// State Machine Configuration
//==============================================================================
template<typename State, typename Event, size_t TransitionCount>
struct StateMachineConfig {
    State initialState;
    std::array<Transition<State, Event>, TransitionCount> transitions;
    
    /// @brief Find target state for given transition
    constexpr std::optional<State> findTransition(State current, Event event) const {
        for (const auto& t : transitions) {
            if (t.from == current && t.event == event) {
                return t.to;
            }
        }
        return std::nullopt;
    }
};

//==============================================================================
// Action Types
//==============================================================================
template<typename State>
using StateAction = std::function<void(State)>;

template<typename State, typename Event>
using TransitionAction = std::function<void(State from, Event event, State to)>;

//==============================================================================
// State Machine
//==============================================================================
template<typename State, typename Event, size_t TransitionCount>
class StateMachine {
public:
    using Config = StateMachineConfig<State, Event, TransitionCount>;
    
    explicit StateMachine(const Config& config)
        : m_Config(config)
        , m_CurrentState(config.initialState)
    {}
    
    /// @brief Get current state (thread-safe)
    State getState() const {
        std::lock_guard<std::recursive_mutex> lock(m_Mutex);
        return m_CurrentState;
    }
    
    /// @brief Check if transition is valid without executing
    bool canProcess(Event event) const {
        std::lock_guard<std::recursive_mutex> lock(m_Mutex);
        return m_Config.findTransition(m_CurrentState, event).has_value();
    }
    
    /// @brief Process event and execute transition if valid
    /// @return true if transition occurred, false if invalid
    bool process(Event event) {
        std::lock_guard<std::recursive_mutex> lock(m_Mutex);
        
        auto nextState = m_Config.findTransition(m_CurrentState, event);
        if (!nextState.has_value()) {
            return false;
        }
        
        State prevState = m_CurrentState;
        State newState = *nextState;
        
        // Exit action
        if (m_ExitAction) {
            m_ExitAction(prevState);
        }
        
        // Transition action
        if (m_TransitionAction) {
            m_TransitionAction(prevState, event, newState);
        }
        
        // State change
        m_CurrentState = newState;
        
        // Entry action
        if (m_EntryAction) {
            m_EntryAction(newState);
        }
        
        return true;
    }
    
    /// @brief Force state (for error recovery only)
    void forceState(State state) {
        std::lock_guard<std::recursive_mutex> lock(m_Mutex);
        m_CurrentState = state;
    }
    
    /// @brief Check if in terminal state
    template<typename... States>
    bool isTerminal(States... terminalStates) const {
        State current = getState();
        return ((current == terminalStates) || ...);
    }
    
    // Action setters
    void setEntryAction(StateAction<State> action) { m_EntryAction = action; }
    void setExitAction(StateAction<State> action) { m_ExitAction = action; }
    void setTransitionAction(TransitionAction<State, Event> action) { m_TransitionAction = action; }
    
private:
    const Config& m_Config;
    State m_CurrentState;
    mutable std::recursive_mutex m_Mutex;  // Recursive to allow actions to call process()
    
    StateAction<State> m_EntryAction;
    StateAction<State> m_ExitAction;
    TransitionAction<State, Event> m_TransitionAction;
};

//==============================================================================
// Helper macros for defining transition tables
//==============================================================================
#define SIRIUS_TRANSITION(from, event, to) \
    Transition<decltype(from), decltype(event)>(from, event, to)

} // namespace Sirius
