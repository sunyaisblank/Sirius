// TSOF003A.cpp - Error Accumulator Tests
// Component ID: TSOF003A
// Tests for: OFER001A.h

#include <gtest/gtest.h>
#include "OFER001A.h"
#include <thread>
#include <vector>

using namespace sirius::offline;

//==============================================================================
// Basic Functionality
//==============================================================================

TEST(ErrorAccumulatorTest, InitiallyEmpty) {
    ErrorAccumulator acc;
    EXPECT_FALSE(acc.hasErrors());
    EXPECT_FALSE(acc.hasFatalError());
    EXPECT_EQ(acc.warningCount(), 0);
    EXPECT_EQ(acc.errorCount(), 0);
}

TEST(ErrorAccumulatorTest, AddWarning) {
    ErrorAccumulator acc;
    acc.add(RenderError::warning("Test warning"));
    
    EXPECT_TRUE(acc.hasErrors());
    EXPECT_FALSE(acc.hasFatalError());
    EXPECT_EQ(acc.warningCount(), 1);
    EXPECT_EQ(acc.errorCount(), 0);
}

TEST(ErrorAccumulatorTest, AddError) {
    ErrorAccumulator acc;
    acc.add(RenderError::error("Test error", "Tile 42"));
    
    EXPECT_TRUE(acc.hasErrors());
    EXPECT_FALSE(acc.hasFatalError());
    EXPECT_EQ(acc.warningCount(), 0);
    EXPECT_EQ(acc.errorCount(), 1);
}

TEST(ErrorAccumulatorTest, AddFatalError) {
    ErrorAccumulator acc;
    acc.add(RenderError::fatal("Out of memory"));
    
    EXPECT_TRUE(acc.hasErrors());
    EXPECT_TRUE(acc.hasFatalError());
}

TEST(ErrorAccumulatorTest, FatalFlagImmediatelyVisible) {
    ErrorAccumulator acc;
    
    // hasFatalError should be visible immediately after add
    acc.add(RenderError::fatal("GPU error"));
    EXPECT_TRUE(acc.hasFatalError());
}

TEST(ErrorAccumulatorTest, Clear) {
    ErrorAccumulator acc;
    acc.add(RenderError::warning("w1"));
    acc.add(RenderError::error("e1"));
    acc.add(RenderError::fatal("f1"));
    
    acc.clear();
    
    EXPECT_FALSE(acc.hasErrors());
    EXPECT_FALSE(acc.hasFatalError());
    EXPECT_EQ(acc.warningCount(), 0);
    EXPECT_EQ(acc.errorCount(), 0);
}

TEST(ErrorAccumulatorTest, GetAllErrors) {
    ErrorAccumulator acc;
    acc.add(RenderError::warning("w1"));
    acc.add(RenderError::error("e1", "ctx1"));
    
    auto errors = acc.all();
    ASSERT_EQ(errors.size(), 2);
    EXPECT_EQ(errors[0].message, "w1");
    EXPECT_EQ(errors[1].message, "e1");
    EXPECT_EQ(errors[1].context, "ctx1");
}

//==============================================================================
// Exception Capture
//==============================================================================

TEST(ErrorAccumulatorTest, FromException) {
    try {
        throw std::runtime_error("Test exception");
    } catch (...) {
        auto error = RenderError::fromException(std::current_exception(), "Test context");
        EXPECT_EQ(error.severity, ErrorSeverity::ERROR);
        EXPECT_NE(error.message.find("Test exception"), std::string::npos);
        EXPECT_EQ(error.context, "Test context");
    }
}

TEST(ErrorAccumulatorTest, FromNonStdException) {
    try {
        throw 42;  // Non-std::exception
    } catch (...) {
        auto error = RenderError::fromException(std::current_exception());
        EXPECT_NE(error.message.find("Non-standard"), std::string::npos);
    }
}

//==============================================================================
// Thread Safety
//==============================================================================

TEST(ErrorAccumulatorTest, ConcurrentAdds) {
    ErrorAccumulator acc;
    constexpr int threadsCount = 10;
    constexpr int errorsPerThread = 100;
    
    std::vector<std::thread> threads;
    for (int t = 0; t < threadsCount; ++t) {
        threads.emplace_back([&acc, t, errorsPerThread]() {
            for (int i = 0; i < errorsPerThread; ++i) {
                acc.add(RenderError::warning("Thread " + std::to_string(t) + " error " + std::to_string(i)));
            }
        });
    }
    
    for (auto& th : threads) {
        th.join();
    }
    
    EXPECT_EQ(acc.warningCount(), threadsCount * errorsPerThread);
}

TEST(ErrorAccumulatorTest, ConcurrentFatalCheck) {
    ErrorAccumulator acc;
    std::atomic<bool> stopReading{false};
    std::atomic<int> checkCount{0};
    
    // Multiple reader threads checking hasFatalError
    std::vector<std::thread> readers;
    for (int i = 0; i < 4; ++i) {
        readers.emplace_back([&]() {
            while (!stopReading) {
                [[maybe_unused]] bool fatal = acc.hasFatalError();
                checkCount++;
                std::this_thread::yield();
            }
        });
    }
    
    // Add fatal error from another thread
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    acc.add(RenderError::fatal("Fatal"));
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    
    stopReading = true;
    for (auto& t : readers) {
        t.join();
    }
    
    EXPECT_TRUE(acc.hasFatalError());
    EXPECT_GT(checkCount.load(), 0);
}

//==============================================================================
// Severity Helpers
//==============================================================================

TEST(RenderErrorTest, SeverityFactories) {
    auto w = RenderError::warning("warn");
    EXPECT_EQ(w.severity, ErrorSeverity::WARNING);
    
    auto e = RenderError::error("err");
    EXPECT_EQ(e.severity, ErrorSeverity::ERROR);
    
    auto f = RenderError::fatal("fatal");
    EXPECT_EQ(f.severity, ErrorSeverity::FATAL);
}
