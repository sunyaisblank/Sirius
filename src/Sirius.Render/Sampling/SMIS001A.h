// SMIS001A.h - Multiple Importance Sampling
// Component ID: SMIS001A (Sampling/ImportanceSampling)
//
// Implements Multiple Importance Sampling (MIS) for combining
// different sampling strategies with optimal variance reduction.
//
// MATHEMATICAL BASIS:
// The MIS estimator combines N sampling strategies using weights:
//   F_MIS = Σ_i (1/n_i) Σ_j [w_i(x_ij) × f(x_ij) / p_i(x_ij)]
//
// Balance heuristic weights (proven optimal):
//   w_i(x) = n_i × p_i(x) / Σ_k(n_k × p_k(x))
//
// Power heuristic (reduces variance in practice):
//   w_i(x) = (n_i × p_i(x))^β / Σ_k(n_k × p_k(x))^β
//
// REFERENCE: Veach & Guibas (1995) "Optimally Combining Sampling Techniques"

#pragma once

#include <vector>
#include <cmath>
#include <functional>

namespace Sirius::Sampling {

//==============================================================================
// MIS Weight Heuristics
//==============================================================================
enum class MISHeuristic : uint8_t {
    Balance,    ///< Balance heuristic: w = p_i / Σp_k
    Power,      ///< Power heuristic: w = p_i^β / Σp_k^β (default β=2)
    Maximum     ///< Maximum heuristic: w = 1 if p_i = max
};

//==============================================================================
// Sample Contribution
//==============================================================================
template<typename T>
struct Sample {
    T value;                ///< Sampled value (radiance, etc.)
    double pdf = 1.0;       ///< PDF of this sample
    int strategyIndex = 0;  ///< Which strategy generated this
    double weight = 1.0;    ///< MIS weight after computation
};

//==============================================================================
// Importance Sampler
//==============================================================================
template<typename T>
class ImportanceSampler {
public:
    using PDFFunction = std::function<double(const T& sample)>;
    
    struct Strategy {
        std::string name;
        PDFFunction pdf;
        int sampleCount = 1;
    };
    
    ImportanceSampler(MISHeuristic heuristic = MISHeuristic::Power, double beta = 2.0)
        : m_Heuristic(heuristic), m_Beta(beta) {}
    
    /// @brief Register a sampling strategy
    void addStrategy(const std::string& name, PDFFunction pdf, int sampleCount = 1) {
        m_Strategies.push_back({name, pdf, sampleCount});
    }
    
    /// @brief Compute MIS weight for a sample
    double computeWeight(const Sample<T>& sample) const {
        if (m_Strategies.empty()) return 1.0;
        
        double numerator = m_Strategies[sample.strategyIndex].sampleCount * sample.pdf;
        double denominator = 0.0;
        
        for (size_t i = 0; i < m_Strategies.size(); ++i) {
            double pdf_i = m_Strategies[i].pdf(sample.value);
            double n_i = m_Strategies[i].sampleCount;
            
            switch (m_Heuristic) {
                case MISHeuristic::Balance:
                    denominator += n_i * pdf_i;
                    break;
                case MISHeuristic::Power:
                    denominator += std::pow(n_i * pdf_i, m_Beta);
                    break;
                case MISHeuristic::Maximum:
                    denominator = std::max(denominator, n_i * pdf_i);
                    break;
            }
        }
        
        if (denominator <= 0) return 0.0;
        
        switch (m_Heuristic) {
            case MISHeuristic::Balance:
                return numerator / denominator;
            case MISHeuristic::Power:
                return std::pow(numerator, m_Beta) / denominator;
            case MISHeuristic::Maximum:
                return (numerator >= denominator) ? 1.0 : 0.0;
        }
        
        return 0.0;
    }
    
    /// @brief Compute weights for all samples
    void computeWeights(std::vector<Sample<T>>& samples) {
        for (auto& s : samples) {
            s.weight = computeWeight(s);
        }
    }
    
private:
    std::vector<Strategy> m_Strategies;
    MISHeuristic m_Heuristic;
    double m_Beta;
};

} // namespace Sirius::Sampling
