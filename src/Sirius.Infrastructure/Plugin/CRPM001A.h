// CRPM001A.h - Plugin Manager
#ifndef CRPM001A_H
#define CRPM001A_H

#include <vector>
#include <string>
#include <memory>
#include "PHMT000A.h"

using namespace Sirius;

class PluginManager {
public:
    void loadPlugins();
    IMetric* getMetric(const std::string& name);
    IMetric* createMetric(const std::string& name) { return getMetric(name); } // Alias for compatibility
    const std::vector<std::string>& getMetricNames() const { return m_MetricNames; }

private:
    std::vector<std::unique_ptr<IMetric>> m_Metrics;
    std::vector<std::string> m_MetricNames;
};

#endif // CRPM001A_H