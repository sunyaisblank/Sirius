// TSOF002A.cpp - Cancellation Token Tests
// Component ID: TSOF002A
// Tests for: OFCT001A.h

#include <gtest/gtest.h>
#include "OFCT001A.h"
#include <thread>
#include <chrono>
#include <atomic>

using namespace sirius::offline;

//==============================================================================
// Basic Functionality
//==============================================================================

TEST(CancellationTokenTest, InitiallyNotCancelled) {
    CancellationToken token;
    EXPECT_FALSE(token.isCancelled());
}

TEST(CancellationTokenTest, RequestSetsCancelled) {
    CancellationToken token;
    token.request();
    EXPECT_TRUE(token.isCancelled());
}

TEST(CancellationTokenTest, ResetClearsCancelled) {
    CancellationToken token;
    token.request();
    EXPECT_TRUE(token.isCancelled());
    token.reset();
    EXPECT_FALSE(token.isCancelled());
}

TEST(CancellationTokenTest, ThrowIfCancelledDoesNotThrowWhenNotCancelled) {
    CancellationToken token;
    EXPECT_NO_THROW(token.throwIfCancelled());
}

TEST(CancellationTokenTest, ThrowIfCancelledThrowsWhenCancelled) {
    CancellationToken token;
    token.request();
    EXPECT_THROW(token.throwIfCancelled(), CancelledException);
}

TEST(CancellationTokenTest, CancelledExceptionMessage) {
    CancelledException ex;
    EXPECT_NE(std::string(ex.what()).find("cancelled"), std::string::npos);
}

//==============================================================================
// Thread Safety
//==============================================================================

TEST(CancellationTokenTest, AtomicVisibilityAcrossThreads) {
    CancellationToken token;
    std::atomic<bool> seen{false};
    
    // Worker thread waits for cancellation
    std::thread worker([&]() {
        while (!token.isCancelled()) {
            std::this_thread::sleep_for(std::chrono::microseconds(100));
        }
        seen = true;
    });
    
    // Give worker time to start
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    
    // Cancel from main thread
    token.request();
    
    // Worker should see it quickly
    worker.join();
    EXPECT_TRUE(seen);
}

TEST(CancellationTokenTest, ConcurrentRequests) {
    CancellationToken token;
    
    // Multiple threads calling request() simultaneously
    // Should not cause data races or undefined behaviour
    std::vector<std::thread> threads;
    for (int i = 0; i < 10; ++i) {
        threads.emplace_back([&]() {
            token.request();
        });
    }
    
    for (auto& t : threads) {
        t.join();
    }
    
    EXPECT_TRUE(token.isCancelled());
}

TEST(CancellationTokenTest, ConcurrentReadsAndWrite) {
    CancellationToken token;
    std::atomic<int> readCount{0};
    std::atomic<bool> stopReading{false};
    
    // Multiple reader threads
    std::vector<std::thread> readers;
    for (int i = 0; i < 4; ++i) {
        readers.emplace_back([&]() {
            while (!stopReading) {
                [[maybe_unused]] bool cancelled = token.isCancelled();
                readCount++;
                std::this_thread::yield();
            }
        });
    }
    
    // Writer thread
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    token.request();
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    
    stopReading = true;
    for (auto& t : readers) {
        t.join();
    }
    
    EXPECT_TRUE(token.isCancelled());
    EXPECT_GT(readCount.load(), 0);
}

//==============================================================================
// Performance (Smoke Test)
//==============================================================================

TEST(CancellationTokenTest, CheckPerformance) {
    CancellationToken token;
    
    // Should be able to check many times per millisecond
    auto start = std::chrono::high_resolution_clock::now();
    constexpr int iterations = 100000;
    
    for (int i = 0; i < iterations; ++i) {
        [[maybe_unused]] bool c = token.isCancelled();
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    
    // Should complete in < 100ms (very conservative)
    EXPECT_LT(ms, 100) << "100k checks took " << ms << "ms";
}
