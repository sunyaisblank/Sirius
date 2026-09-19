#include "sirius/render/session/progress_tracker.h"

#include <gtest/gtest.h>

#include <chrono>

using sirius::render::ProgressTracker;
using namespace std::chrono_literals;

TEST(ProgressTrackerTest, WorkerCompletionBurstRetainsItsTracingTime) {
    ProgressTracker::Clock::time_point now{};
    ProgressTracker progress([&] { return now; });
    progress.Start();
    progress.SetTotals(24, 3);
    EXPECT_LT(progress.GetEta(), 0);
    double callback_eta = -1;
    progress.SetCallback([&](float, int, int, double eta) { callback_eta = eta; });

    // Fifteen workers spend 100 seconds tracing, then publish their first wave
    // over three seconds. Publication speed cannot stand in for tracing speed.
    now += 100s;
    for (int tile = 0; tile < 15; ++tile) {
        now += 200ms;
        progress.CompleteTile(3);
    }
    EXPECT_EQ(progress.GetTilesComplete(), 15);
    EXPECT_FLOAT_EQ(progress.GetProgress(), .625f);
    EXPECT_GT(callback_eta, 60);
    EXPECT_LT(callback_eta, 70);
    EXPECT_GT(progress.GetEta(), 60);
    EXPECT_LT(progress.GetTilesPerSecond(), .15);
}

TEST(ProgressTrackerTest, QuietWorkerWaveUpdatesEtaAndRestartClearsItsHistory) {
    ProgressTracker::Clock::time_point now{};
    ProgressTracker progress([&] { return now; });
    progress.Start();
    progress.SetTotals(4, 1);
    now += 10s;
    progress.CompleteTile();
    progress.CompleteTile();
    const double previous_eta = progress.GetEta();
    now += 30s;
    EXPECT_GT(progress.GetEta(), previous_eta);
    EXPECT_DOUBLE_EQ(progress.GetElapsedSeconds(), 40);
    progress.CompleteTile();
    progress.CompleteTile();
    EXPECT_DOUBLE_EQ(progress.GetEta(), 0);
    EXPECT_EQ(progress.GetEtaString(), "0s");

    progress.Start();
    progress.SetTotals(24, 3);
    EXPECT_EQ(progress.GetTilesComplete(), 0);
    EXPECT_FLOAT_EQ(progress.GetProgress(), 0);
    EXPECT_LT(progress.GetEta(), 0);
    EXPECT_DOUBLE_EQ(progress.GetTilesPerSecond(), 0);
}
