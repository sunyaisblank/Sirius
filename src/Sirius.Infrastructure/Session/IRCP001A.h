// IRCP001A.h - CLI Parser
// Component ID: IRCP001A (Infrastructure/Render/CLI Parser)
//
// Command-line argument parsing with validation and help generation.

#pragma once

#include <string>
#include <vector>
#include <map>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace Sirius {

//==============================================================================
// Argument Types
//==============================================================================
enum class ArgType {
    Flag,       ///< Boolean flag (--verbose)
    String,     ///< String value (--output path.exr)
    Integer,    ///< Integer value (--width 1920)
    Float,      ///< Float value (--exposure 1.5)
    Double      ///< Double value (--spin 0.9)
};

//==============================================================================
// Argument Definition
//==============================================================================
struct ArgDef {
    std::string name;           ///< Long name (--output)
    char shortName = 0;         ///< Short name (-o), 0 for none
    ArgType type = ArgType::Flag;
    std::string description;
    std::string defaultValue;
    bool required = false;
};

//==============================================================================
// CLI Parser
//==============================================================================
class CLIParser {
public:
    CLIParser(const std::string& programName, const std::string& description = "")
        : m_ProgramName(programName), m_Description(description) {}
    
    /// @brief Add an argument definition
    void addArgument(const ArgDef& arg) {
        m_Arguments.push_back(arg);
        if (arg.shortName != 0) {
            m_ShortMap[arg.shortName] = arg.name;
        }
    }
    
    /// @brief Add flag (--flag)
    void addFlag(const std::string& name, char shortName, const std::string& desc) {
        addArgument({name, shortName, ArgType::Flag, desc, "false", false});
    }
    
    /// @brief Add string option (--option value)
    void addString(const std::string& name, char shortName, const std::string& desc,
                   const std::string& defaultVal = "", bool required = false) {
        addArgument({name, shortName, ArgType::String, desc, defaultVal, required});
    }
    
    /// @brief Add integer option (--width 1920)
    void addInt(const std::string& name, char shortName, const std::string& desc,
                int defaultVal = 0, bool required = false) {
        addArgument({name, shortName, ArgType::Integer, desc, std::to_string(defaultVal), required});
    }
    
    /// @brief Add float option (--exposure 1.5)
    void addFloat(const std::string& name, char shortName, const std::string& desc,
                  float defaultVal = 0.0f, bool required = false) {
        addArgument({name, shortName, ArgType::Float, desc, std::to_string(defaultVal), required});
    }
    
    /// @brief Add double option (--spin 0.9)
    void addDouble(const std::string& name, char shortName, const std::string& desc,
                   double defaultVal = 0.0, bool required = false) {
        addArgument({name, shortName, ArgType::Double, desc, std::to_string(defaultVal), required});
    }
    
    /// @brief Parse command line arguments
    /// @return true if parsing succeeded, false if help requested or error
    bool parse(int argc, char* argv[]) {
        // Set defaults
        for (const auto& arg : m_Arguments) {
            m_Values[arg.name] = arg.defaultValue;
        }
        
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            
            if (arg == "--help" || arg == "-h") {
                printHelp();
                return false;
            }
            
            // Long option
            if (arg.substr(0, 2) == "--") {
                std::string name = arg.substr(2);
                const ArgDef* def = findArg(name);
                if (!def) {
                    std::cerr << "Unknown option: " << arg << std::endl;
                    return false;
                }
                
                if (def->type == ArgType::Flag) {
                    m_Values[name] = "true";
                } else if (i + 1 < argc) {
                    m_Values[name] = argv[++i];
                } else {
                    std::cerr << "Missing value for: " << arg << std::endl;
                    return false;
                }
            }
            // Short option
            else if (arg[0] == '-' && arg.length() == 2) {
                char c = arg[1];
                auto it = m_ShortMap.find(c);
                if (it == m_ShortMap.end()) {
                    std::cerr << "Unknown option: " << arg << std::endl;
                    return false;
                }
                
                const ArgDef* def = findArg(it->second);
                if (def->type == ArgType::Flag) {
                    m_Values[it->second] = "true";
                } else if (i + 1 < argc) {
                    m_Values[it->second] = argv[++i];
                } else {
                    std::cerr << "Missing value for: " << arg << std::endl;
                    return false;
                }
            }
            else {
                m_Positional.push_back(arg);
            }
        }
        
        // Check required arguments
        for (const auto& arg : m_Arguments) {
            if (arg.required && m_Values[arg.name].empty()) {
                std::cerr << "Missing required argument: --" << arg.name << std::endl;
                return false;
            }
        }
        
        return true;
    }
    
    // Value getters
    bool getFlag(const std::string& name) const {
        auto it = m_Values.find(name);
        return it != m_Values.end() && it->second == "true";
    }
    
    std::string getString(const std::string& name) const {
        auto it = m_Values.find(name);
        return it != m_Values.end() ? it->second : "";
    }
    
    int getInt(const std::string& name) const {
        auto it = m_Values.find(name);
        return it != m_Values.end() ? std::stoi(it->second) : 0;
    }
    
    float getFloat(const std::string& name) const {
        auto it = m_Values.find(name);
        return it != m_Values.end() ? std::stof(it->second) : 0.0f;
    }
    
    double getDouble(const std::string& name) const {
        auto it = m_Values.find(name);
        return it != m_Values.end() ? std::stod(it->second) : 0.0;
    }
    
    const std::vector<std::string>& getPositional() const { return m_Positional; }
    
    /// @brief Print help message
    void printHelp() const {
        std::cout << m_ProgramName;
        if (!m_Description.empty()) {
            std::cout << " - " << m_Description;
        }
        std::cout << "\n\nUsage:\n  " << m_ProgramName << " [options]\n\nOptions:\n";
        
        for (const auto& arg : m_Arguments) {
            std::cout << "  --" << arg.name;
            if (arg.shortName != 0) {
                std::cout << ", -" << arg.shortName;
            }
            
            if (arg.type != ArgType::Flag) {
                std::cout << " <value>";
            }
            
            std::cout << "\n      " << arg.description;
            if (!arg.defaultValue.empty() && arg.type != ArgType::Flag) {
                std::cout << " (default: " << arg.defaultValue << ")";
            }
            if (arg.required) {
                std::cout << " [required]";
            }
            std::cout << "\n";
        }
        
        std::cout << "  --help, -h\n      Show this help message\n";
    }
    
private:
    const ArgDef* findArg(const std::string& name) const {
        for (const auto& arg : m_Arguments) {
            if (arg.name == name) return &arg;
        }
        return nullptr;
    }
    
    std::string m_ProgramName;
    std::string m_Description;
    std::vector<ArgDef> m_Arguments;
    std::map<char, std::string> m_ShortMap;
    std::map<std::string, std::string> m_Values;
    std::vector<std::string> m_Positional;
};

} // namespace Sirius
