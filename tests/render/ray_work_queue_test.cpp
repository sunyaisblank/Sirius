#include "sirius/render/session/ray_work_queue.h"

#include <gtest/gtest.h>

#include <array>
#include <atomic>
#include <chrono>
#include <future>
#include <stdexcept>

namespace sirius::test {
namespace {
using backend::TraceResult;
using core::CameraRay;
using render::RayWorkQueue;

TEST(RayWorkQueue, IndependentRaysOverlapAndReturnInRequestOrder) {
    std::mutex mutex;
    std::condition_variable completed;
    bool second_finished = false;
    bool overlapped = false;
    std::array<std::atomic<int>, 2> active{};
    std::atomic<bool> worker_reentered{false};
    RayWorkQueue queue(2, [&](std::size_t worker, const CameraRay& ray) {
        if (active[worker].fetch_add(1) != 0) worker_reentered = true;
        const int index = static_cast<int>(ray.aperture_up);
        if (index == 0) {
            std::unique_lock lock(mutex);
            overlapped =
                completed.wait_for(lock, std::chrono::seconds(5), [&] { return second_finished; });
        } else if (index == 1) {
            {
                std::lock_guard lock(mutex);
                second_finished = true;
            }
            completed.notify_all();
        }
        TraceResult result;
        result.steps_taken = index;
        result.outcome = TraceResult::Outcome::Escaped;
        --active[worker];
        return result;
    });
    // Exceeds queue residency, exercising producer backpressure and reuse.
    std::array<CameraRay, RayWorkQueue::kMaximumBatch> rays;
    for (std::size_t i = 0; i < rays.size(); ++i) rays[i].aperture_up = static_cast<double>(i);
    const auto result = queue.Execute(rays);
    ASSERT_TRUE(result);
    EXPECT_TRUE(overlapped);
    EXPECT_FALSE(worker_reentered.load());
    ASSERT_EQ(result->size(), rays.size());
    for (std::size_t i = 0; i < rays.size(); ++i)
        EXPECT_EQ((*result)[i].steps_taken, static_cast<int>(i));
}

TEST(RayWorkQueue, ConcurrentCallersKeepTheirOwnResults) {
    RayWorkQueue queue(3, [](std::size_t, const CameraRay& ray) {
        TraceResult result;
        result.affine_length = ray.aperture_up;
        return result;
    });
    std::array<std::future<bool>, 4> callers;
    for (std::size_t caller = 0; caller < callers.size(); ++caller)
        callers[caller] = std::async(std::launch::async, [&, caller] {
            std::array<CameraRay, 32> rays;
            for (std::size_t i = 0; i < rays.size(); ++i)
                rays[i].aperture_up = static_cast<double>(100 * caller + i);
            const auto result = queue.Execute(rays);
            if (!result || result->size() != rays.size()) return false;
            for (std::size_t i = 0; i < rays.size(); ++i)
                if ((*result)[i].affine_length != rays[i].aperture_up) return false;
            return true;
        });
    for (auto& caller : callers) EXPECT_TRUE(caller.get());
}

TEST(RayWorkQueue, WorkerFailureDrainsBatchAndQueueRemainsUsable) {
    std::atomic<int> traced{0};
    RayWorkQueue queue(2, [&](std::size_t, const CameraRay& ray) {
        ++traced;
        if (ray.aperture_up == 3) throw std::runtime_error("injected trace failure");
        TraceResult result;
        result.cancelled = ray.aperture_up < 0;
        result.numerical_failure = ray.aperture_up > 0;
        return result;
    });
    std::array<CameraRay, 13> rays;
    rays[5].aperture_up = 3;
    const auto failed = queue.Execute(rays);
    ASSERT_FALSE(failed);
    EXPECT_EQ(failed.error().domain(), base::ErrorDomain::kInternal);
    EXPECT_EQ(traced.load(), 13);
    rays[5].aperture_up = -1;
    rays[6].aperture_up = 1;
    const auto next = queue.Execute(rays);
    ASSERT_TRUE(next);
    EXPECT_TRUE((*next)[5].cancelled);
    EXPECT_TRUE((*next)[6].numerical_failure);
    EXPECT_EQ(traced.load(), 26);
}

TEST(RayWorkQueue, RejectsUnboundedWorkBeforeCallingTracer) {
    std::size_t calls = 0;
    RayWorkQueue queue(1, [&](std::size_t, const CameraRay&) {
        ++calls;
        return TraceResult{};
    });
    std::array<CameraRay, RayWorkQueue::kMaximumBatch + 1> rays;
    EXPECT_FALSE(queue.Execute(rays));
    const auto empty = queue.Execute({});
    ASSERT_TRUE(empty);
    EXPECT_TRUE(empty->empty());
    EXPECT_EQ(calls, 0u);
    EXPECT_THROW((RayWorkQueue(0, {})), std::invalid_argument);
    EXPECT_THROW((RayWorkQueue(1025, {})), std::invalid_argument);
    EXPECT_THROW((RayWorkQueue(1, {})), std::invalid_argument);
}
}  // namespace
}  // namespace sirius::test
