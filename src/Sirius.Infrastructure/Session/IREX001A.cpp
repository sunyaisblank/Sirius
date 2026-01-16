// IREX001A.cpp - Job Executor Implementation
//
// Manages job queue and worker thread for batch rendering.

#include "IREX001A.h"
#include <iostream>
#include <chrono>

namespace Sirius {

//==============================================================================
// Constructor / Destructor
//==============================================================================
JobExecutor::JobExecutor(const ExecutorConfig& config)
    : m_Config(config)
{
    if (m_Config.verboseLogging) {
        std::cout << "[JobExecutor] Initialized with config:" << std::endl;
        std::cout << "  Max concurrent jobs: " << m_Config.maxConcurrentJobs << std::endl;
        std::cout << "  Preview enabled: " << (m_Config.enablePreview ? "yes" : "no") << std::endl;
    }
}

JobExecutor::~JobExecutor() {
    stop();
}

//==============================================================================
// Job Management
//==============================================================================
void JobExecutor::submit(std::shared_ptr<RenderJob> job) {
    std::lock_guard<std::mutex> lock(m_QueueMutex);
    m_JobQueue.push(std::move(job));
    m_QueueCV.notify_one();
    
    if (m_Config.verboseLogging) {
        std::cout << "[JobExecutor] Job submitted. Queue size: " << m_JobQueue.size() << std::endl;
    }
}

void JobExecutor::submitBatch(const std::vector<RenderConfig>& configs) {
    std::lock_guard<std::mutex> lock(m_QueueMutex);
    
    for (const auto& config : configs) {
        auto job = std::make_shared<RenderJob>(config);
        m_JobQueue.push(std::move(job));
    }
    
    m_QueueCV.notify_one();
    
    if (m_Config.verboseLogging) {
        std::cout << "[JobExecutor] Batch submitted: " << configs.size() 
                  << " jobs. Queue size: " << m_JobQueue.size() << std::endl;
    }
}

//==============================================================================
// Status
//==============================================================================
int JobExecutor::getPendingJobCount() const {
    std::lock_guard<std::mutex> lock(m_QueueMutex);
    return static_cast<int>(m_JobQueue.size());
}

std::shared_ptr<RenderJob> JobExecutor::getCurrentJob() const {
    std::lock_guard<std::mutex> lock(m_CurrentJobMutex);
    return m_CurrentJob;
}

//==============================================================================
// Callbacks
//==============================================================================
void JobExecutor::setJobStartCallback(JobStartCallback callback) {
    m_JobStartCallback = std::move(callback);
}

void JobExecutor::setJobCompleteCallback(JobCompleteCallback callback) {
    m_JobCompleteCallback = std::move(callback);
}

void JobExecutor::setAllCompleteCallback(AllCompleteCallback callback) {
    m_AllCompleteCallback = std::move(callback);
}

//==============================================================================
// Execution Control
//==============================================================================
void JobExecutor::start() {
    if (m_Running.load()) {
        std::cerr << "[JobExecutor] Already running" << std::endl;
        return;
    }
    
    m_Running.store(true);
    m_StopRequested.store(false);
    m_CompletedJobs.store(0);
    m_SuccessfulJobs.store(0);
    
    m_WorkerThread = std::thread(&JobExecutor::workerThread, this);
    
    if (m_Config.enablePreview) {
        m_PreviewThread = std::thread(&JobExecutor::previewThread, this);
    }
    
    if (m_Config.verboseLogging) {
        std::cout << "[JobExecutor] Started" << std::endl;
    }
}

void JobExecutor::stop() {
    if (!m_Running.load()) return;
    
    m_StopRequested.store(true);
    m_QueueCV.notify_all();
    
    // Cancel current job if running
    {
        std::lock_guard<std::mutex> lock(m_CurrentJobMutex);
        if (m_CurrentJob && m_CurrentJob->isRunning()) {
            m_CurrentJob->cancel();
        }
    }
    
    if (m_WorkerThread.joinable()) {
        m_WorkerThread.join();
    }
    
    if (m_PreviewThread.joinable()) {
        m_PreviewThread.join();
    }
    
    m_Running.store(false);
    
    if (m_Config.verboseLogging) {
        std::cout << "[JobExecutor] Stopped" << std::endl;
    }
}

void JobExecutor::waitForCompletion() {
    if (m_WorkerThread.joinable()) {
        m_WorkerThread.join();
    }
    
    if (m_PreviewThread.joinable()) {
        m_PreviewThread.join();
    }
    
    m_Running.store(false);
}

//==============================================================================
// Worker Thread
//==============================================================================
void JobExecutor::workerThread() {
    if (m_Config.verboseLogging) {
        std::cout << "[JobExecutor] Worker thread started" << std::endl;
    }
    
    while (!m_StopRequested.load()) {
        std::shared_ptr<RenderJob> job;
        
        // Get next job from queue
        {
            std::unique_lock<std::mutex> lock(m_QueueMutex);
            
            // Wait for job or stop signal
            m_QueueCV.wait(lock, [this] {
                return !m_JobQueue.empty() || m_StopRequested.load();
            });
            
            if (m_StopRequested.load() && m_JobQueue.empty()) {
                break;
            }
            
            if (!m_JobQueue.empty()) {
                job = m_JobQueue.front();
                m_JobQueue.pop();
            }
        }
        
        if (!job) continue;
        
        // Set as current job
        {
            std::lock_guard<std::mutex> lock(m_CurrentJobMutex);
            m_CurrentJob = job;
        }
        
        // Notify job start
        if (m_JobStartCallback) {
            m_JobStartCallback(*job);
        }
        
        if (m_Config.verboseLogging) {
            std::cout << "[JobExecutor] Starting job: " << job->getConfig().outputPath << std::endl;
        }
        
        // Execute the job
        job->execute();
        
        // Update statistics
        m_CompletedJobs.fetch_add(1);
        if (job->getStatus() == JobStatus::Completed) {
            m_SuccessfulJobs.fetch_add(1);
        }
        
        // Notify job complete
        if (m_JobCompleteCallback) {
            m_JobCompleteCallback(*job, job->getStatus());
        }
        
        // Clear current job
        {
            std::lock_guard<std::mutex> lock(m_CurrentJobMutex);
            m_CurrentJob = nullptr;
        }
        
        // Check if all jobs are done
        {
            std::lock_guard<std::mutex> lock(m_QueueMutex);
            if (m_JobQueue.empty()) {
                if (m_AllCompleteCallback) {
                    m_AllCompleteCallback(m_CompletedJobs.load(), m_SuccessfulJobs.load());
                }
                break;  // Exit worker thread
            }
        }
    }
    
    if (m_Config.verboseLogging) {
        std::cout << "[JobExecutor] Worker thread finished. Completed: " 
                  << m_CompletedJobs.load() << " jobs (" 
                  << m_SuccessfulJobs.load() << " successful)" << std::endl;
    }
}

//==============================================================================
// Preview Thread
//==============================================================================
void JobExecutor::previewThread() {
    // TODO: Implement preview window using GLFW
    // For now, just poll and log progress
    
    if (m_Config.verboseLogging) {
        std::cout << "[JobExecutor] Preview thread started (placeholder)" << std::endl;
    }
    
    while (!m_StopRequested.load() && m_Running.load()) {
        auto job = getCurrentJob();
        
        if (job && job->isRunning()) {
            [[maybe_unused]] float progress = job->getProgress();
            // TODO: In a real implementation, update preview window here
        }
        
        std::this_thread::sleep_for(
            std::chrono::milliseconds(m_Config.previewUpdateInterval)
        );
    }
    
    if (m_Config.verboseLogging) {
        std::cout << "[JobExecutor] Preview thread finished" << std::endl;
    }
}

} // namespace Sirius
