// CRPM001A.cpp - Plugin Manager Implementation
#include "CRPM001A.h"
#include "PHMT000A.h"

// Unified Metric Families (Dec 2025)
#include <PHMT100A.h>  // Kerr-Schild Family
#include <PHMT101A.h>  // Morris-Thorne Family
#include <PHMT102A.h>  // Warp Drive Family

using namespace Sirius;

void PluginManager::loadPlugins() {
    // =========================================================================
    // Kerr-Schild Family - covers 9 black hole spacetimes
    // =========================================================================
    // Default: Schwarzschild (M=1, a=Q=Λ=0)
    m_Metrics.push_back(std::make_unique<KerrSchildFamily>());
    m_MetricNames.push_back("Schwarzschild");
    
    // Kerr (M=1, a=0.9, Q=Λ=0)
    m_Metrics.push_back(std::make_unique<KerrSchildFamily>(KerrSchildParams::Kerr(1.0, 0.9)));
    m_MetricNames.push_back("Kerr");
    
    // Minkowski (M=a=Q=Λ=0)
    m_Metrics.push_back(std::make_unique<KerrSchildFamily>(KerrSchildParams::Minkowski()));
    m_MetricNames.push_back("Minkowski");
    
    // Reissner-Nordstrom (M=1, a=0, Q=0.5, Λ=0)  
    m_Metrics.push_back(std::make_unique<KerrSchildFamily>(KerrSchildParams::ReissnerNordstrom(1.0, 0.5)));
    m_MetricNames.push_back("Reissner-Nordstrom");
    
    // Kerr-Newman (M=1, a=0.6, Q=0.3, Λ=0)
    m_Metrics.push_back(std::make_unique<KerrSchildFamily>(KerrSchildParams::KerrNewman(1.0, 0.6, 0.3)));
    m_MetricNames.push_back("Kerr-Newman");
    
    // de Sitter (M=a=Q=0, Λ=0.01)
    m_Metrics.push_back(std::make_unique<KerrSchildFamily>(KerrSchildParams::DeSitter(0.01)));
    m_MetricNames.push_back("de-Sitter");
    
    // =========================================================================
    // Morris-Thorne Family - traversable wormholes
    // =========================================================================
    // Ellis Drainhole (zero-tidal wormhole)
    m_Metrics.push_back(std::make_unique<MorrisThorneFamily>(MorrisThorneParams::Ellis(1.0)));
    m_MetricNames.push_back("Ellis-Drainhole");
    
    // =========================================================================
    // Warp Drive Family - Alcubierre-class metrics
    // =========================================================================
    // Alcubierre (v=2, R=1)
    m_Metrics.push_back(std::make_unique<WarpDriveFamily>(WarpDriveParams::Superluminal(2.0, 1.0)));
    m_MetricNames.push_back("Alcubierre");
}

IMetric* PluginManager::getMetric(const std::string& name) {
    for (size_t i = 0; i < m_MetricNames.size(); ++i) {
        if (m_MetricNames[i] == name) {
            return m_Metrics[i].get();
        }
    }
    return nullptr;
}