#include "sirius/render/session/render_session.h"
#include "sirius/render/vulkan_renderer.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <numbers>
#include <vector>

namespace {
using namespace sirius::render;
constexpr int kEdge = 16;
constexpr int kPixels = kEdge * kEdge;
constexpr char kCertifiedMode[] = "SIRIUS_CERTIFIED_SPHERICAL_ESCAPE";

void Check(bool passed, const char* message) {
    if (passed) return;
    std::fprintf(stderr, "Failed: %s\n", message);
    std::exit(EXIT_FAILURE);
}

struct Frame {
    VulkanRenderStats stats;
    std::vector<float> rgba;
};

Frame Render(const SessionConfig& config, bool certified) {
    const int selected = certified ? setenv(kCertifiedMode, "1", 1) : unsetenv(kCertifiedMode);
    Check(selected == 0, "set development mode");
    DisplayBuffer display;
    display.Initialise(kEdge, kEdge);
    std::array<unsigned, kPixels> coverage{};
    std::mutex mutex;
    auto rendered = RenderVulkanToDisplay(
        config, display, {}, {},
        [&](int x, int y, int width, int height, std::span<const float> values) {
            Check(x >= 0 && y >= 0 && width > 0 && height > 0 && x + width <= kEdge &&
                      y + height <= kEdge &&
                      values.size() == static_cast<unsigned>(4 * width * height),
                  "represented tile extent");
            Check(std::all_of(values.begin(), values.end(),
                              [](float value) { return std::isfinite(value) && value >= 0.0f; }),
                  "complete finite tile");
            std::lock_guard lock(mutex);
            for (int row = y; row < y + height; ++row)
                for (int column = x; column < x + width; ++column) ++coverage[row * kEdge + column];
        });
    if (!rendered) {
        std::fprintf(stderr, "%s\n", rendered.error().Description().c_str());
        std::exit(EXIT_FAILURE);
    }
    Check(std::all_of(coverage.begin(), coverage.end(), [](unsigned count) { return count == 1; }),
          "every pixel is published once");
    auto rgba = display.SnapshotFloatData();
    Check(rgba.size() == kPixels * 4, "committed RGBA extent");
    bool nonzero_rgb = false;
    for (unsigned i = 0; i < rgba.size(); ++i) {
        nonzero_rgb |= i % 4 != 3 && rgba[i] > 0.0f;
        Check(std::isfinite(rgba[i]) && rgba[i] >= 0.0f && (i % 4 != 3 || rgba[i] == 1.0f),
              "represented committed radiance");
    }
    Check(nonzero_rgb, "nonzero committed RGB");
    return {std::move(*rendered), std::move(rgba)};
}

void CheckCancelled(const SessionConfig& config) {
    DisplayBuffer display;
    display.Initialise(kEdge, kEdge);
    display.MutateFloatData([](auto& values) { std::fill(values.begin(), values.end(), .125f); });
    const auto before = display.SnapshotFloatData();
    std::atomic<bool> cancel{false};
    const auto result = RenderVulkanToDisplay(
        config, display, {}, [&] { return cancel.load(); },
        [&](int, int, int, int, std::span<const float>) { cancel.store(true); });
    Check(cancel.load() && !result && display.SnapshotFloatData() == before,
          "cancelled render preserves the committed framebuffer");
}

double Error(const Frame& a, const Frame& b) {
    double maximum = 0.0;
    for (unsigned i = 0; i < a.rgba.size(); ++i) {
        const double x = a.rgba[i], y = b.rgba[i];
        maximum = std::max(maximum, std::abs(x - y) / std::max({1.0, std::abs(x), std::abs(y)}));
    }
    return maximum;
}
}  // namespace

int main() {
    SessionConfig config;
    config.metric_id = sirius::core::MetricId::Kerr;
    config.backend = RenderBackend::Vulkan;
    config.observer_inclination = std::numbers::pi / 2;
    config.width = config.height = kEdge;
    config.samples_per_pixel = 1;
    config.write_output = false;
    const auto ordinary = Render(config, false);
    const auto certified = Render(config, true);
    CheckCancelled(config);
    const auto retried = Render(config, true);
    const auto warm_ordinary = Render(config, false);
    const auto& a = ordinary.stats;
    const auto& b = certified.stats;
    Check(a.retained_intervals && b.retained_intervals && a.camera_batches > 0 &&
              b.camera_batches > 0,
          "same retained renderer route");
    Check(a.certified_spherical_escapes == 0 && b.certified_spherical_escapes > 0 &&
              b.certified_spherical_turns > 0 && b.certified_spherical_panels > 0,
          "both analytic continuations publish through the renderer");
    Check(a.accepted_intervals > b.accepted_intervals &&
              a.tile_plan.usable_bytes == b.tile_plan.usable_bytes &&
              a.explicit_buffer_allocation_bytes == b.explicit_buffer_allocation_bytes &&
              a.continuation_capacity == b.continuation_capacity &&
              a.work_tile_edge == b.work_tile_edge && a.precision == b.precision,
          "work reduction under the same residency and route");
    const auto error = Error(ordinary, certified);
    Check(
        error <= 1e-5 && Error(certified, retried) == 0.0 && Error(ordinary, warm_ordinary) == 0.0,
        "paired radiance and deterministic retry");
    std::printf(
        "pixels=%d ordinary_retained=%llu certified_retained=%llu terminal_events=%llu turns=%llu "
        "quadrature_panels=%llu scaled_radiance_error=%.9g retry_error=%.9g "
        "ordinary_s=%.6f warm_ordinary_s=%.6f certified_s=%.6f retried_s=%.6f "
        "allocation_bytes=%llu capacity=%llu\n",
        kPixels, static_cast<unsigned long long>(a.accepted_intervals),
        static_cast<unsigned long long>(b.accepted_intervals),
        static_cast<unsigned long long>(b.certified_spherical_escapes),
        static_cast<unsigned long long>(b.certified_spherical_turns),
        static_cast<unsigned long long>(b.certified_spherical_panels), error,
        Error(certified, retried), a.seconds, warm_ordinary.stats.seconds, b.seconds,
        retried.stats.seconds, static_cast<unsigned long long>(b.explicit_buffer_allocation_bytes),
        static_cast<unsigned long long>(b.continuation_capacity));
}
