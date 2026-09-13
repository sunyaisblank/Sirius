#pragma once

#include "sirius/backend/device.h"

#include <array>
#include <memory>

namespace sirius::backend {

// Three-term binary32 expansion with a separately retained arithmetic radius.
// Camera outputs use two terms and set tail to zero.
// This describes the finite arithmetic operation, not the ODE truncation error
// or a global trajectory enclosure. Both low parts and radii cross readback.
struct RetainedValue {
    float high = 0;
    float low = 0;
    float tail = 0;
    float radius = 0;
    std::uint32_t valid = 0;

    [[nodiscard]] static RetainedValue FromDouble(double value);
    [[nodiscard]] bool IsRepresented() const;
    // Rounded diagnostic value; consumers needing the complete expansion must
    // use all three stored components directly.
    [[nodiscard]] double Center() const { return (double(high) + double(low)) + double(tail); }
};
static_assert(sizeof(RetainedValue) == 20);

struct RetainedCameraInput {
    // Includes continuous film pixels and physical Cartesian pupil offsets;
    // these must not be rounded to one float before the device launch.
    std::array<RetainedValue, 32> values{};
};
struct RetainedCameraOutput {
    // x, k, frame, X, coordinate K, covariant V, observer du, repeated frame.
    std::array<RetainedValue, 104> values{};
    bool valid = false;
};

struct RetainedStepInput {
    // M, a, Q, Lambda; x[4], p[4]; four (X[4], delta-p[4]) columns;
    // chart reflection (+1 ingoing, -1 outgoing); positive affine interval.
    std::array<RetainedValue, 46> values{};
};
struct RetainedStepOutput {
    // All phase records use x[4], p[4], four (X[4], delta-p[4]) columns.
    std::array<RetainedValue, 40> fifth{}, fourth{}, increment{}, error{};
    // Attempted RHS evaluations, including an evaluation that declined.
    std::uint32_t stages = 0;
    bool valid = false;
};

struct RetainedEndpointInput {
    // Metric parameters, complete phase record, and exact chart reflection.
    std::array<RetainedValue, 45> values{};
};
struct RetainedEndpointOutput {
    // Projected phase (x,p,X,P) stays retained for subsequent device stages.
    // Physical (x,k,X,V) is the coupled acceptance/dense-event representation.
    std::array<RetainedValue, 40> phase{}, physical{};
    std::uint32_t component = 4;
    bool valid = false;
};

struct RetainedDenseInput {
    // M,a,Q,L; two physical (x,k,X,V) endpoints; retained x and four X
    // increments; affine interval, fraction, normal[4], moving flag, chart.
    std::array<RetainedValue, 112> values{};
};
struct RetainedDenseOutput {
    std::array<RetainedValue, 40> physical{};
    bool valid = false;
};

// Fixed-capacity device stages. The caller owns admission, continuation and
// source publication; a private RK candidate is never an accepted trajectory.
// No device allocation occurs during Camera or Step. Device must outlive this
// object, and callers must serialize its synchronous submissions.
class RetainedCompute {
  public:
    [[nodiscard]] static base::Expected<std::unique_ptr<RetainedCompute>> Create(
        ComputeDevice& device, std::size_t capacity);
    [[nodiscard]] base::Expected<std::vector<RetainedCameraOutput>> Camera(
        std::span<const RetainedCameraInput> inputs, DispatchTiming* timing = nullptr);
    [[nodiscard]] base::Expected<std::vector<RetainedStepOutput>> Step(
        std::span<const RetainedStepInput> inputs, DispatchTiming* timing = nullptr);
    [[nodiscard]] base::Expected<std::vector<RetainedEndpointOutput>> Endpoint(
        std::span<const RetainedEndpointInput> inputs, DispatchTiming* timing = nullptr);
    [[nodiscard]] base::Expected<std::vector<RetainedDenseOutput>> Dense(
        std::span<const RetainedDenseInput> inputs, DispatchTiming* timing = nullptr);
    [[nodiscard]] std::size_t Capacity() const { return capacity_; }

  private:
    explicit RetainedCompute(ComputeDevice& device, std::size_t capacity)
        : device_(device), capacity_(capacity) {}
    struct Stage {
        KernelHandle kernel;
        std::array<BufferHandle, 2> buffers;
        std::vector<std::uint32_t> input, output;
    };
    [[nodiscard]] base::Expected<void> Dispatch(Stage& stage, DispatchTiming* timing);
    ComputeDevice& device_;
    std::size_t capacity_;
    Stage camera_, transport_, endpoint_, dense_;
};

}  // namespace sirius::backend
