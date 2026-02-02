// IRLG001A.h - Logger
// Component ID: IRLG001A (Infrastructure/Render/Logger)
//
// Thread-safe logging with timestamps, log levels, and optional file output.

#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <mutex>
#include <memory>

namespace Sirius {

//==============================================================================
// Log Levels
//==============================================================================
enum class LogLevel {
    Debug,
    Info,
    Warning,
    Error,
    Fatal
};

//==============================================================================
// Logger - Thread-safe logging
//==============================================================================
class Logger {
public:
    static Logger& instance() {
        static Logger logger;
        return logger;
    }
    
    /// @brief Set minimum log level
    void setLevel(LogLevel level) { m_Level = level; }
    
    /// @brief Enable file logging
    void setLogFile(const std::string& path) {
        std::lock_guard<std::mutex> lock(m_Mutex);
        m_LogFile = std::make_unique<std::ofstream>(path, std::ios::app);
    }
    
    /// @brief Log a message
    void log(LogLevel level, const std::string& message, 
             const char* file = nullptr, int line = 0) {
        if (level < m_Level) return;
        
        std::lock_guard<std::mutex> lock(m_Mutex);
        
        std::string formatted = formatMessage(level, message, file, line);
        
        // Console output (with colour)
        std::ostream& out = (level >= LogLevel::Warning) ? std::cerr : std::cout;
        out << colourForLevel(level) << formatted << "\033[0m" << std::endl;
        
        // File output
        if (m_LogFile && m_LogFile->is_open()) {
            *m_LogFile << formatted << std::endl;
        }
    }
    
    // Convenience methods
    void debug(const std::string& msg) { log(LogLevel::Debug, msg); }
    void info(const std::string& msg) { log(LogLevel::Info, msg); }
    void warn(const std::string& msg) { log(LogLevel::Warning, msg); }
    void error(const std::string& msg) { log(LogLevel::Error, msg); }
    void fatal(const std::string& msg) { log(LogLevel::Fatal, msg); }
    
private:
    Logger() = default;
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;
    
    std::string formatMessage(LogLevel level, const std::string& message,
                               const char* file, int line) {
        std::ostringstream ss;

        // Timestamp (thread-safe version)
        auto now = std::chrono::system_clock::now();
        auto time = std::chrono::system_clock::to_time_t(now);

        // Use localtime_r on POSIX systems for thread-safety
        std::tm timeinfo{};
#if defined(_WIN32) || defined(_WIN64)
        localtime_s(&timeinfo, &time);  // Windows thread-safe version
#else
        localtime_r(&time, &timeinfo);  // POSIX thread-safe version
#endif
        ss << std::put_time(&timeinfo, "%H:%M:%S");
        
        // Level
        ss << " [" << levelString(level) << "] ";
        
        // Message
        ss << message;
        
        // File/line (for debug)
        if (file && level == LogLevel::Debug) {
            ss << " (" << file << ":" << line << ")";
        }
        
        return ss.str();
    }
    
    const char* levelString(LogLevel level) {
        switch (level) {
            case LogLevel::Debug:   return "DEBUG";
            case LogLevel::Info:    return "INFO ";
            case LogLevel::Warning: return "WARN ";
            case LogLevel::Error:   return "ERROR";
            case LogLevel::Fatal:   return "FATAL";
        }
        return "?????";
    }
    
    const char* colourForLevel(LogLevel level) {
        switch (level) {
            case LogLevel::Debug:   return "\033[90m";  // Grey
            case LogLevel::Info:    return "\033[37m";  // White
            case LogLevel::Warning: return "\033[33m";  // Yellow
            case LogLevel::Error:   return "\033[31m";  // Red
            case LogLevel::Fatal:   return "\033[91m";  // Bright red
        }
        return "";
    }
    
    LogLevel m_Level = LogLevel::Info;
    std::unique_ptr<std::ofstream> m_LogFile;
    std::mutex m_Mutex;
};

// Convenience macros
#define LOG_DEBUG(msg) Sirius::Logger::instance().log(Sirius::LogLevel::Debug, msg, __FILE__, __LINE__)
#define LOG_INFO(msg)  Sirius::Logger::instance().info(msg)
#define LOG_WARN(msg)  Sirius::Logger::instance().warn(msg)
#define LOG_ERROR(msg) Sirius::Logger::instance().error(msg)
#define LOG_FATAL(msg) Sirius::Logger::instance().fatal(msg)

} // namespace Sirius
