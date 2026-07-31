#pragma once

// Generic deterministic finite-state machine with a compile-time transition
// table and entry/exit/transition actions. Ported from SMSM001A.h.
//
// The transition table is a std::array searched linearly; for the render
// session's handful of states this is smaller and faster than a hashed map, and
// it is constexpr-constructible so the table is validated at compile time.

#include <array>
#include <functional>
#include <mutex>
#include <optional>

namespace sirius::render {

// A single (from, event) -> to edge.
template <typename State, typename Event>
struct Transition {
    State from;
    Event event;
    State to;

    constexpr Transition(State f, Event e, State t) : from(f), event(e), to(t) {}
};

// Initial state plus the full transition table.
template <typename State, typename Event, size_t TransitionCount>
struct StateMachineConfig {
    State initial_state;
    std::array<Transition<State, Event>, TransitionCount> transitions;

    // Target state for (current, event), or nullopt when no edge exists.
    constexpr std::optional<State> FindTransition(State current, Event event) const {
        for (const auto& t : transitions) {
            if (t.from == current && t.event == event) {
                return t.to;
            }
        }
        return std::nullopt;
    }
};

template <typename State>
using StateAction = std::function<void(State)>;

template <typename State, typename Event>
using TransitionAction = std::function<void(State from, Event event, State to)>;

// Mealy machine: transitions fire exit, transition, then entry actions.
template <typename State, typename Event, size_t TransitionCount>
class StateMachine {
  public:
    using Config = StateMachineConfig<State, Event, TransitionCount>;

    explicit StateMachine(const Config& config)
        : config_(config), current_state_(config.initial_state) {}

    // Current state (thread-safe).
    State GetState() const {
        std::lock_guard<std::recursive_mutex> lock(mutex_);
        return current_state_;
    }

    // Whether an event has a valid transition from the current state.
    bool CanProcess(Event event) const {
        std::lock_guard<std::recursive_mutex> lock(mutex_);
        return config_.FindTransition(current_state_, event).has_value();
    }

    // Process an event; returns true when a transition occurred. State changes
    // are serialised, but actions run after releasing the state mutex so a
    // long-running entry action cannot make GetState() observationally block
    // for the duration of the operation it launched.
    bool Process(Event event) {
        State prev_state;
        State new_state;
        StateAction<State> exit_action;
        StateAction<State> entry_action;
        TransitionAction<State, Event> transition_action;
        {
            std::lock_guard<std::recursive_mutex> lock(mutex_);
            auto next_state = config_.FindTransition(current_state_, event);
            if (!next_state.has_value()) {
                return false;
            }
            prev_state = current_state_;
            new_state = *next_state;
            current_state_ = new_state;
            exit_action = exit_action_;
            entry_action = entry_action_;
            transition_action = transition_action_;
        }

        if (exit_action) {
            exit_action(prev_state);
        }
        if (transition_action) {
            transition_action(prev_state, event, new_state);
        }
        if (entry_action) {
            entry_action(new_state);
        }

        return true;
    }

    // Force a state (error recovery only).
    void ForceState(State state) {
        std::lock_guard<std::recursive_mutex> lock(mutex_);
        current_state_ = state;
    }

    // Whether the current state matches any of the given terminal states.
    template <typename... States>
    bool IsTerminal(States... terminal_states) const {
        State current = GetState();
        return ((current == terminal_states) || ...);
    }

    void SetEntryAction(StateAction<State> action) {
        std::lock_guard<std::recursive_mutex> lock(mutex_);
        entry_action_ = std::move(action);
    }
    void SetExitAction(StateAction<State> action) {
        std::lock_guard<std::recursive_mutex> lock(mutex_);
        exit_action_ = std::move(action);
    }
    void SetTransitionAction(TransitionAction<State, Event> action) {
        std::lock_guard<std::recursive_mutex> lock(mutex_);
        transition_action_ = std::move(action);
    }

  private:
    const Config& config_;
    State current_state_;
    // Recursive so an entry/transition action may itself call Process().
    mutable std::recursive_mutex mutex_;

    StateAction<State> entry_action_;
    StateAction<State> exit_action_;
    TransitionAction<State, Event> transition_action_;
};

}  // namespace sirius::render
