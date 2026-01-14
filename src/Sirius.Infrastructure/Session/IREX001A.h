// IREX001A.h - Job Executor
//
// Manages render job queue, GPU resource allocation, and optional preview window.
// Central coordinator for batch rendering operations.

#pragma once

#include "IRRJ001A.h"
#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>

namespace Sirius {

//==============================================================================
// Executor Configuration
//==============================================================================
struct ExecutorConfig {
    int maxConcurrentJobs = 1;      // GPU jobs are serialized by default
    bool enablePreview = false;      // Optional preview window
    int previewWidth = 640;
    int previewHeight = 360;
    int previewUpdateInterval = 100; // ms between preview updates
    bool verboseLogging = true;
};

//==============================================================================
// Job Executor
//==============================================================================
class JobExecutor {
public:
    explicit JobExecutor(const ExecutorConfig& config = ExecutorConfig());
    ~JobExecutor();
    
    // Non-copyable
    JobExecutor(const JobExecutor&) = delete;
    JobExecutor& operator=(const JobExecutor&) = delete;
    
    // Job management
    void submit(std::shared_ptr<RenderJob> job);
    void submitBatch(const std::vector<RenderConfig>& configs);
    
    // Execution control
    void start();
    void stop();
    void waitForCompletion();
    
    // Status
    bool isRunning() const { return m_Running.load(); }
    int getPendingJobCount() const;
    int getCompletedJobCount() const { return m_CompletedJobs.load(); }
    
    // Current job access
    std::shared_ptr<RenderJob> getCurrentJob() const;
    
    // Global callbacks
    using JobStartCallback = std::function<void(const RenderJob& job)>;
    using JobCompleteCallback = std::function<void(const RenderJob& job, JobStatus status)>;
    using AllCompleteCallback = std::function<void(int totalJobs, int successful)>;
    
    void setJobStartCallback(JobStartCallback callback);
    void setJobCompleteCallback(JobCompleteCallback callback);
    void setAllCompleteCallback(AllCompleteCallback callback);
    
private:
    void workerThread();
    void previewThread();
    
    ExecutorConfig m_Config;
    
    // Job queue
    std::queue<std::shared_ptr<RenderJob>> m_JobQueue;
    mutable std::mutex m_QueueMutex;
    std::condition_variable m_QueueCV;
    
    // Current execution
    std::shared_ptr<RenderJob> m_CurrentJob;
    mutable std::mutex m_CurrentJobMutex;
    
    // Worker thread
    std::thread m_WorkerThread;
    std::thread m_PreviewThread;
    std::atomic<bool> m_Running{false};
    std::atomic<bool> m_StopRequested{false};
    
    // Statistics
    std::atomic<int> m_CompletedJobs{0};
    std::atomic<int> m_SuccessfulJobs{0};
    
    // Callbacks
    JobStartCallback m_JobStartCallback;
    JobCompleteCallback m_JobCompleteCallback;
    AllCompleteCallback m_AllCompleteCallback;
};

//==============================================================================
// Convenience: Single Job Execution
//==============================================================================
inline void executeRenderJob(const RenderConfig& config) {
    auto job = std::make_shared<RenderJob>(config);
    job->execute();
}

//==============================================================================
// Convenience: Batch Execution
//==============================================================================
inline void executeRenderBatch(const std::vector<RenderConfig>& configs) {
    ExecutorConfig execConfig;
    execConfig.verboseLogging = true;
    
    JobExecutor executor(execConfig);
    executor.submitBatch(configs);
    executor.start();
    executor.waitForCompletion();
}

} // namespace Sirius
