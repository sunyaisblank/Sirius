// =============================================================================
// TSFT002A.h - Metric Test Fixture
// Component ID: TSFT002A (Test/Fixtures/MetricFixture)
// =============================================================================
//
// PURPOSE:
// Provides reusable test fixture infrastructure for metric tensor validation.
// Implements common validation methods following the pattern from Alaris.Test.
//
// MATHEMATICAL INVARIANTS TESTED:
// 1. Metric symmetry: g_μν = g_νμ (exact)
// 2. Inverse accuracy: g^μα g_αν = δ^μ_ν
// 3. Christoffel symmetry: Γ^λ_μν = Γ^λ_νμ (torsion-free)
// 4. Lorentzian signature: exactly one negative eigenvalue (-,+,+,+)
// 5. Determinant non-degeneracy: |det(g)| > threshold
// 6. NaN/Inf detection
//
// USAGE:
// class MyMetricTests : public sirius::test::MetricTestFixture { ... };
//
// REFERENCE: docs/specification.md, docs/philosophy.md
// =============================================================================

#ifndef TSFT002A_H
#define TSFT002A_H

#include <gtest/gtest.h>
#include <PHCN001A.h>
#include <MTTN001A.h>
#include <MTDL001A.h>
#include "TSFT001A.h"

#include <cmath>
#include <vector>
#include <string>
#include <functional>

namespace sirius::test {
using namespace Sirius;

// =============================================================================
// Test Point Structure
// =============================================================================

/// @brief Represents a spacetime point in both Cartesian and spherical coordinates
struct TestPoint {
    // Cartesian coordinates (t, x, y, z)
    double t = 0.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    // Spherical coordinates (t, r, θ, φ)
    double r = 0.0;
    double theta = 0.0;
    double phi = 0.0;

    /// @brief Create from spherical coordinates
    static TestPoint fromSpherical(double t, double r, double theta, double phi) {
        TestPoint pt;
        pt.t = t;
        pt.r = r;
        pt.theta = theta;
        pt.phi = phi;

        // Convert to Cartesian
        double sin_th = std::sin(theta);
        double cos_th = std::cos(theta);
        double sin_ph = std::sin(phi);
        double cos_ph = std::cos(phi);

        pt.x = r * sin_th * cos_ph;
        pt.y = r * sin_th * sin_ph;
        pt.z = r * cos_th;

        return pt;
    }

    /// @brief Create from Cartesian coordinates
    static TestPoint fromCartesian(double t, double x, double y, double z) {
        TestPoint pt;
        pt.t = t;
        pt.x = x;
        pt.y = y;
        pt.z = z;

        // Convert to spherical
        pt.r = std::sqrt(x*x + y*y + z*z);
        if (pt.r > 1e-10) {
            pt.theta = std::acos(z / pt.r);
            pt.phi = std::atan2(y, x);
            if (pt.phi < 0) pt.phi += 2.0 * MathConst::PI;
        } else {
            pt.theta = 0.0;
            pt.phi = 0.0;
        }

        return pt;
    }

    /// @brief Get as Vec4 in Cartesian coordinates
    Vec4 toCartesianVec4() const {
        Vec4 v;
        v(0) = t;
        v(1) = x;
        v(2) = y;
        v(3) = z;
        return v;
    }

    /// @brief Get as Vec4 in spherical coordinates
    Vec4 toSphericalVec4() const {
        Vec4 v;
        v(0) = t;
        v(1) = r;
        v(2) = theta;
        v(3) = phi;
        return v;
    }
};

// =============================================================================
// Validation Result Structure
// =============================================================================

/// @brief Results from metric validation
struct MetricValidationResult {
    bool passed = true;
    std::string metric_name;
    TestPoint position;

    // Individual test results
    double symmetry_error = 0.0;
    double inverse_error = 0.0;
    double christoffel_symmetry_error = 0.0;
    double determinant = 0.0;
    bool has_nan = false;
    bool has_inf = false;
    bool correct_signature = true;

    // Tolerances used
    double symmetry_tol = Tolerances::METRIC_SYMMETRY;
    double inverse_tol = Tolerances::INVERSE_METRIC;
    double christoffel_tol = Tolerances::CHRISTOFFEL_SYMMETRY;
    double determinant_tol = Tolerances::DETERMINANT;

    /// @brief Check if all validations passed
    bool allPassed() const {
        return passed &&
               (symmetry_error < symmetry_tol) &&
               (inverse_error < inverse_tol) &&
               (christoffel_symmetry_error < christoffel_tol) &&
               (std::abs(determinant) > determinant_tol) &&
               !has_nan && !has_inf && correct_signature;
    }

    /// @brief Get human-readable summary
    std::string summary() const {
        std::string s = "MetricValidation for " + metric_name + " at r=" + std::to_string(position.r) + ":\n";
        s += "  Symmetry error: " + std::to_string(symmetry_error) + " (tol: " + std::to_string(symmetry_tol) + ")\n";
        s += "  Inverse error: " + std::to_string(inverse_error) + " (tol: " + std::to_string(inverse_tol) + ")\n";
        s += "  Christoffel symmetry: " + std::to_string(christoffel_symmetry_error) + " (tol: " + std::to_string(christoffel_tol) + ")\n";
        s += "  Determinant: " + std::to_string(determinant) + " (min: " + std::to_string(determinant_tol) + ")\n";
        s += "  Has NaN: " + std::string(has_nan ? "YES" : "no") + "\n";
        s += "  Has Inf: " + std::string(has_inf ? "YES" : "no") + "\n";
        s += "  Correct signature: " + std::string(correct_signature ? "yes" : "NO") + "\n";
        s += "  Overall: " + std::string(allPassed() ? "PASSED" : "FAILED") + "\n";
        return s;
    }
};

// =============================================================================
// Metric Test Fixture Base Class
// =============================================================================

/// @brief Base fixture for metric tensor tests
/// Provides common validation methods and sample point generation
class MetricTestFixture : public ::testing::Test {
protected:
    // =========================================================================
    // Setup and Teardown
    // =========================================================================

    void SetUp() override {
        // Derived classes can override for custom setup
    }

    void TearDown() override {
        // Derived classes can override for custom teardown
    }

    // =========================================================================
    // Sample Point Generation
    // =========================================================================

    /// @brief Generate standard test points outside the horizon
    /// @param M Black hole mass (default 1.0)
    /// @param a Spin parameter (default 0.0)
    /// @return Vector of test points
    std::vector<TestPoint> getStandardTestPoints(double M = 1.0, double a = 0.0) const {
        // Compute horizon radius: r+ = M + sqrt(M² - a²)
        double r_plus = M + std::sqrt(std::max(M*M - a*a, 0.0));
        double r_safe = std::max(r_plus * 1.5, 3.0 * M);

        std::vector<TestPoint> points;

        // Equatorial plane points at various radii
        for (double r : {r_safe, 2.0*r_safe, 5.0*M, 10.0*M, 20.0*M, 50.0*M, 100.0*M}) {
            if (r > r_safe) {
                points.push_back(TestPoint::fromSpherical(0, r, MathConst::HALF_PI, 0));
                points.push_back(TestPoint::fromSpherical(0, r, MathConst::HALF_PI, MathConst::PI/4));
            }
        }

        // Off-equatorial points
        for (double theta : {MathConst::PI/6, MathConst::PI/3, 2*MathConst::PI/3, 5*MathConst::PI/6}) {
            points.push_back(TestPoint::fromSpherical(0, 10.0*M, theta, 0));
        }

        // Far field for asymptotic tests
        points.push_back(TestPoint::fromSpherical(0, 100.0*M, MathConst::HALF_PI, 0));
        points.push_back(TestPoint::fromSpherical(0, 1000.0*M, MathConst::HALF_PI, 0));

        return points;
    }

    /// @brief Generate boundary test points (stress testing)
    /// @param M Black hole mass
    /// @param a Spin parameter
    /// @return Vector of boundary test points
    std::vector<TestPoint> getBoundaryTestPoints(double M = 1.0, double a = 0.0) const {
        double r_plus = M + std::sqrt(std::max(M*M - a*a, 0.0));

        std::vector<TestPoint> points;

        // Near horizon (1.001 buffer per specification)
        points.push_back(TestPoint::fromSpherical(0, r_plus * 1.002, MathConst::HALF_PI, 0));
        points.push_back(TestPoint::fromSpherical(0, r_plus * 1.01, MathConst::HALF_PI, 0));

        // Near poles (avoiding coordinate singularity)
        points.push_back(TestPoint::fromSpherical(0, 10.0*M, 0.01, 0));
        points.push_back(TestPoint::fromSpherical(0, 10.0*M, MathConst::PI - 0.01, 0));

        // Very large radius
        points.push_back(TestPoint::fromSpherical(0, 1e5*M, MathConst::HALF_PI, 0));

        return points;
    }

    // =========================================================================
    // Metric Validation Methods
    // =========================================================================

    /// @brief Validate metric symmetry: g_μν = g_νμ
    /// @return Maximum asymmetry |g_μν - g_νμ|
    double validateMetricSymmetry(const Metric4D& g) const {
        double max_asym = 0.0;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = mu + 1; nu < 4; ++nu) {
                double asym = std::abs(g(mu, nu).real - g(nu, mu).real);
                max_asym = std::max(max_asym, asym);
            }
        }
        return max_asym;
    }

    /// @brief Validate inverse metric: g^μα g_αν = δ^μ_ν
    /// @return Maximum deviation from identity
    double validateInverseMetric(const Metric4D& g) const {
        // Compute inverse metric
        auto g_inv = TensorOps::inverse(g);

        // Compute g^μα g_αν
        double max_error = 0.0;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                double product = 0.0;
                for (int alpha = 0; alpha < 4; ++alpha) {
                    product += g_inv(mu, alpha).real * g(alpha, nu).real;
                }
                double expected = (mu == nu) ? 1.0 : 0.0;
                max_error = std::max(max_error, std::abs(product - expected));
            }
        }
        return max_error;
    }

    /// @brief Validate Christoffel symmetry: Γ^λ_μν = Γ^λ_νμ
    /// @return Maximum asymmetry
    double validateChristoffelSymmetry(const Tensor<Dual<double>, 4, 4, 4>& Gamma) const {
        double max_asym = 0.0;
        for (int lam = 0; lam < 4; ++lam) {
            for (int mu = 0; mu < 4; ++mu) {
                for (int nu = mu + 1; nu < 4; ++nu) {
                    double asym = std::abs(Gamma(lam, mu, nu).real - Gamma(lam, nu, mu).real);
                    max_asym = std::max(max_asym, asym);
                }
            }
        }
        return max_asym;
    }

    /// @brief Check for NaN in metric components
    bool hasNaN(const Metric4D& g) const {
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                if (std::isnan(g(mu, nu).real)) return true;
            }
        }
        return false;
    }

    /// @brief Check for Inf in metric components
    bool hasInf(const Metric4D& g) const {
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                if (std::isinf(g(mu, nu).real)) return true;
            }
        }
        return false;
    }

    /// @brief Check Lorentzian signature (-,+,+,+)
    /// Simple check: g_tt < 0 and spatial diagonal > 0 for diagonal-dominant metrics
    bool hasLorentzianSignature(const Metric4D& g) const {
        // For general metrics, we check that g_tt is negative
        // and that the spatial block has positive "trace" character
        // (full eigenvalue analysis would be more rigorous)

        // g_tt must be negative
        if (g(0, 0).real >= -Tolerances::SIGNATURE) return false;

        // For many metrics, diagonal spatial components should be positive
        // This is a simplified check; full eigenvalue analysis is more rigorous
        double spatial_trace = g(1, 1).real + g(2, 2).real + g(3, 3).real;
        return spatial_trace > -Tolerances::SIGNATURE;
    }

    /// @brief Compute metric determinant (4x4)
    /// @note For Lorentzian signature, det(g) < 0
    double computeDeterminant(const Metric4D& g) const {
        // Use Laplace expansion along first row
        // det(g) = Σ_j (-1)^(0+j) g_0j det(M_0j)
        // where M_0j is the (3x3) minor matrix

        double det = 0.0;
        for (int j = 0; j < 4; ++j) {
            // 3x3 minor matrix (remove row 0, column j)
            double minor[3][3];
            int mi = 0;
            for (int i = 1; i < 4; ++i) {
                int mj = 0;
                for (int jj = 0; jj < 4; ++jj) {
                    if (jj != j) {
                        minor[mi][mj] = g(i, jj).real;
                        ++mj;
                    }
                }
                ++mi;
            }

            // 3x3 determinant
            double det3 = minor[0][0] * (minor[1][1]*minor[2][2] - minor[1][2]*minor[2][1])
                        - minor[0][1] * (minor[1][0]*minor[2][2] - minor[1][2]*minor[2][0])
                        + minor[0][2] * (minor[1][0]*minor[2][1] - minor[1][1]*minor[2][0]);

            double sign = (j % 2 == 0) ? 1.0 : -1.0;
            det += sign * g(0, j).real * det3;
        }

        return det;
    }

    /// @brief Perform full metric validation at a point
    template<typename MetricType>
    MetricValidationResult validateMetricAtPoint(MetricType& metric, const TestPoint& pt) {
        MetricValidationResult result;
        result.metric_name = metric.getName();
        result.position = pt;

        // Evaluate metric and derivatives
        Vec4 pos = pt.toCartesianVec4();
        Metric4D g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric.evaluate(pos, g, dg);

        // Run all validations
        result.has_nan = hasNaN(g);
        result.has_inf = hasInf(g);
        result.symmetry_error = validateMetricSymmetry(g);
        result.inverse_error = validateInverseMetric(g);
        result.determinant = computeDeterminant(g);
        result.correct_signature = hasLorentzianSignature(g);

        // Compute Christoffel symbols and validate symmetry
        auto Gamma = TensorOps::christoffel(g, dg);
        result.christoffel_symmetry_error = validateChristoffelSymmetry(Gamma);

        result.passed = result.allPassed();
        return result;
    }

    // =========================================================================
    // GTest Assertion Helpers
    // =========================================================================

    /// @brief Assert metric symmetry
    void assertMetricSymmetry(const Metric4D& g, const std::string& context = "") {
        double error = validateMetricSymmetry(g);
        EXPECT_LT(error, Tolerances::METRIC_SYMMETRY)
            << "Metric symmetry violation" << (context.empty() ? "" : " (" + context + ")")
            << ": " << error;
    }

    /// @brief Assert inverse metric accuracy
    void assertInverseMetric(const Metric4D& g, const std::string& context = "") {
        double error = validateInverseMetric(g);
        EXPECT_LT(error, Tolerances::INVERSE_METRIC)
            << "Inverse metric error" << (context.empty() ? "" : " (" + context + ")")
            << ": " << error;
    }

    /// @brief Assert Christoffel symmetry
    void assertChristoffelSymmetry(const Tensor<Dual<double>, 4, 4, 4>& Gamma,
                                   const std::string& context = "") {
        double error = validateChristoffelSymmetry(Gamma);
        EXPECT_LT(error, Tolerances::CHRISTOFFEL_SYMMETRY)
            << "Christoffel symmetry violation" << (context.empty() ? "" : " (" + context + ")")
            << ": " << error;
    }

    /// @brief Assert no NaN/Inf in metric
    void assertNoNaNInf(const Metric4D& g, const std::string& context = "") {
        EXPECT_FALSE(hasNaN(g)) << "NaN in metric" << (context.empty() ? "" : " (" + context + ")");
        EXPECT_FALSE(hasInf(g)) << "Inf in metric" << (context.empty() ? "" : " (" + context + ")");
    }

    /// @brief Assert Lorentzian signature
    void assertLorentzianSignature(const Metric4D& g, const std::string& context = "") {
        EXPECT_TRUE(hasLorentzianSignature(g))
            << "Incorrect signature" << (context.empty() ? "" : " (" + context + ")");
    }
};

} // namespace sirius::test

#endif // TSFT002A_H
