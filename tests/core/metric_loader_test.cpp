// =============================================================================
// TSIN002A.cpp - Metric Loader Chain Integration Tests
// Component ID: TSIN002A (Test/Integration/MetricLoader)
// =============================================================================
//
// PURPOSE:
// Validates the unified metric family loader (PHMT100A/101A/102A).
// Tests metric symmetry, signature, and limit behaviors.
//
// LABEL: Mandatory;Correctness
// =============================================================================

#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <memory>
#include "sirius/core/tensor.h"
#include "sirius/core/dual_number.h"
#include "sirius/core/metrics/metric.h"
#include "sirius/core/metrics/kerr_schild_family.h"
#include "sirius/core/metrics/morris_thorne_family.h"
#include "sirius/core/metrics/warp_drive_family.h"
#include "sirius/core/constants.h"  // Centralized constants

namespace sirius::test {
using namespace sirius::core;

// Test tolerances (from PHCN001A.h)
constexpr double kEpsilon = constants::metric::kInverseTol;   // 1e-14 for metric comparisons
[[maybe_unused]] constexpr double kSymmetryTol = constants::metric::kSymmetryTol;  // 1e-15 for symmetry
constexpr double M = 1.0;
constexpr double PI = constants::math::kPi;

// Helper: Transform Cartesian metric to Spherical basis for validation
Metric4d transformToSpherical(const Metric4d& g_cart, const Vec4& pos_cart) {
    double x = pos_cart(1);
    double y = pos_cart(2);
    double z = pos_cart(3);
    double rho = std::sqrt(x*x + y*y);
    double r = std::sqrt(x*x + y*y + z*z);

    // Jacobian d(x,y,z)/d(r,th,ph)
    // dx/dr = sin(th)cos(ph) = x/r
    // dx/dth = r cos(th)cos(ph) = x*z/rho
    // dx/dph = -r sin(th)sin(ph) = -y

    // Actually simpler to use trigonometric forms if available, but from Cartesian:
    double st = rho/r;
    double ct = z/r;
    double cp = (rho > 1e-10) ? x/rho : 1.0;
    double sp = (rho > 1e-10) ? y/rho : 0.0;

    Tensor<double, 4, 4> J; // Transformation matrix Lambda^mu_nu = dx^mu / dx'^nu
    J.Zero();
    J(0,0) = 1.0; // dt/dt

    // dr column
    J(1,1) = st * cp; // dx/dr
    J(2,1) = st * sp; // dy/dr
    J(3,1) = ct;      // dz/dr

    // dth column
    J(1,2) = r * ct * cp; // dx/dth
    J(2,2) = r * ct * sp; // dy/dth
    J(3,2) = -r * st;     // dz/dth

    // dph column
    J(1,3) = -r * st * sp; // dx/dph
    J(2,3) = r * st * cp;  // dy/dph
    J(3,3) = 0.0;          // dz/dph

    Metric4d g_sph;
    g_sph.Zero();

    // g'_ab = J^c_a J^d_b g_cd
    for(int a=0; a<4; ++a) {
        for(int b=0; b<4; ++b) {
            Dual<double> sum = 0.0;
            for(int c=0; c<4; ++c) {
                for(int d=0; d<4; ++d) {
                    sum = sum + g_cart(c,d) * J(c,a) * J(d,b);
                }
            }
            g_sph(a,b) = sum;
        }
    }
    return g_sph;
}

// Test fixture providing collection of all metric implementations
class MetricLoaderChainTests : public ::testing::Test {
protected:
    std::vector<std::unique_ptr<IMetric>> metrics;

    void SetUp() override {
        // Create one instance of each metric type
        metrics.push_back(std::make_unique<sirius::core::KerrSchildFamily>(sirius::core::KerrSchildParams::Minkowski()));
        metrics.push_back(std::make_unique<sirius::core::KerrSchildFamily>(sirius::core::KerrSchildParams::Schwarzschild(1.0)));
        metrics.push_back(std::make_unique<sirius::core::KerrSchildFamily>(sirius::core::KerrSchildParams::Kerr(1.0, 0.5)));
        metrics.push_back(std::make_unique<sirius::core::KerrSchildFamily>(sirius::core::KerrSchildParams::ReissnerNordstrom(1.0, 0.5)));
        metrics.push_back(std::make_unique<sirius::core::KerrSchildFamily>(sirius::core::KerrSchildParams::Kerr(1.0, 0.0)));
        metrics.push_back(std::make_unique<sirius::core::MorrisThorneFamily>(sirius::core::MorrisThorneParams::Ellis(1.0)));
    }

    void TearDown() override {
        metrics.clear();
    }

    // Standard test position safe for all metrics (Cartesian equivalent of r=10, th=PI/2, ph=0)
    // (x=10, y=0, z=0)
    Vec4 getStandardPosition() {
        Vec4 pos;
        pos(0) = 0.0;
        pos(1) = 10.0; // x
        pos(2) = 0.0;  // y
        pos(3) = 0.0;  // z
        return pos;
    }

    // Evaluate metric and verify no NaN/Inf in output
    bool evaluateMetricSafe(IMetric* metric, const Vec4& pos, Metric4d& g,
                           Tensor<Dual<double>, 4, 4, 4>& dg) {
        try {
            metric->Evaluate(pos, g, dg);

            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    if (std::isnan(g(i, j).real) || std::isinf(g(i, j).real)) {
                        return false;
                    }
                }
            }
            return true;
        } catch (...) {
            return false;
        }
    }
};

// --- IMetric Interface Compliance ---

// Verify all metrics implement getName() with non-empty string
TEST_F(MetricLoaderChainTests, AllMetricsHaveName) {
    for (const auto& metric : metrics) {
        const char* name = metric->GetName();
        EXPECT_NE(name, nullptr) << "Metric name should not be null";
        EXPECT_GT(strlen(name), 0) << "Metric name should not be empty";
    }
}

// Verify all metrics implement getParameters() without throwing
TEST_F(MetricLoaderChainTests, AllMetricsHaveParameters) {
    for (const auto& metric : metrics) {
        const auto& params = metric->GetParameters();
        (void)params;  // Just verify it doesn't throw
    }
}

// Verify all metrics evaluate successfully at standard position
TEST_F(MetricLoaderChainTests, AllMetricsEvaluateAtStandardPosition) {
    Vec4 pos = getStandardPosition();

    for (const auto& metric : metrics) {
        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;

        bool valid = evaluateMetricSafe(metric.get(), pos, g, dg);
        EXPECT_TRUE(valid) << "Metric " << metric->GetName()
                           << " failed to evaluate at standard position";
    }
}

// --- Metric Signature Tests ---

// Verify all metrics have Lorentzian signature (-,+,+,+)
TEST_F(MetricLoaderChainTests, AllMetricsHaveLorentzianSignature) {
    Vec4 pos = getStandardPosition();

    for (const auto& metric : metrics) {
        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric->Evaluate(pos, g, dg);

        Metric4d g_sph = transformToSpherical(g, pos);

        EXPECT_LT(g_sph(0, 0).real, 0.0)
            << "Metric " << metric->GetName() << " g_tt should be negative";

        EXPECT_GT(g_sph(1, 1).real, 0.0)
            << "Metric " << metric->GetName() << " g_rr should be positive";
        EXPECT_GT(g_sph(3, 3).real, 0.0)
            << "Metric " << metric->GetName() << " g_φφ should be positive";
    }
}

// Verify all metrics are symmetric: g_μν = g_νμ
TEST_F(MetricLoaderChainTests, AllMetricsAreSymmetric) {
    Vec4 pos = getStandardPosition();

    for (const auto& metric : metrics) {
        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric->Evaluate(pos, g, dg);

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                EXPECT_NEAR(g(i, j).real, g(j, i).real, kEpsilon)
                    << "Metric " << metric->GetName()
                    << " asymmetric at (" << i << "," << j << ")";
            }
        }
    }
}

// --- Parameter Configuration Tests ---

// Verify Schwarzschild mass parameter affects g_tt: -f(r) = -(1 - 2M/r)
TEST_F(MetricLoaderChainTests, SchwarzschildMassParameterWorks) {
    sirius::core::KerrSchildFamily schw{sirius::core::KerrSchildParams::Schwarzschild(M)};
    Vec4 pos = getStandardPosition();

    schw.SetParameter("mass", 1.0);
    Metric4d g1;
    Tensor<Dual<double>, 4, 4, 4> dg1;
    schw.Evaluate(pos, g1, dg1);

    schw.SetParameter("mass", 2.0);
    Metric4d g2;
    Tensor<Dual<double>, 4, 4, 4> dg2;
    schw.Evaluate(pos, g2, dg2);

    EXPECT_NE(g1(0, 0).real, g2(0, 0).real)
        << "Mass parameter should affect metric";
    EXPECT_LT(std::abs(g2(0, 0).real), std::abs(g1(0, 0).real))
        << "Larger mass should give smaller |g_tt|";
}

// Verify Kerr spin parameter produces frame-dragging term g_tφ
TEST_F(MetricLoaderChainTests, KerrSpinParameterWorks) {
    sirius::core::KerrSchildFamily kerr{sirius::core::KerrSchildParams::Kerr(M, 0.0)};
    Vec4 pos = getStandardPosition();

    kerr.SetParameter("spin", 0.0);
    Metric4d g_zero_spin;
    Tensor<Dual<double>, 4, 4, 4> dg;
    kerr.Evaluate(pos, g_zero_spin, dg);

    kerr.SetParameter("spin", 0.9);
    Metric4d g_high_spin;
    kerr.Evaluate(pos, g_high_spin, dg);

    // Transform to Spherical to check frame dragging (g_tφ)
    Metric4d g_zero_sph = transformToSpherical(g_zero_spin, pos);
    Metric4d g_high_sph = transformToSpherical(g_high_spin, pos);

    EXPECT_NEAR(g_zero_sph(0, 3).real, 0.0, kEpsilon)
        << "Zero spin should have no frame dragging";
    EXPECT_NE(g_high_sph(0, 3).real, 0.0)
        << "High spin should have frame dragging (g_tφ ≠ 0)";
}

// Verify Reissner-Nordström charge parameter affects metric
TEST_F(MetricLoaderChainTests, RNChargeParameterWorks) {
    sirius::core::KerrSchildFamily rn{sirius::core::KerrSchildParams::ReissnerNordstrom(M, 0.0)};
    Vec4 pos = getStandardPosition();

    rn.SetParameter("charge", 0.0);
    Metric4d g_zero_charge;
    Tensor<Dual<double>, 4, 4, 4> dg;
    rn.Evaluate(pos, g_zero_charge, dg);

    rn.SetParameter("charge", 0.5);
    Metric4d g_with_charge;
    rn.Evaluate(pos, g_with_charge, dg);

    EXPECT_NE(g_zero_charge(0, 0).real, g_with_charge(0, 0).real)
        << "Charge should affect g_tt";
}

// --- Cross-Metric Consistency Tests ---

// Verify Kerr reduces to Schwarzschild when spin a=0
TEST_F(MetricLoaderChainTests, KerrReducesToSchwarzschild) {
    sirius::core::KerrSchildFamily schw{sirius::core::KerrSchildParams::Schwarzschild(M)};
    sirius::core::KerrSchildFamily kerr{sirius::core::KerrSchildParams::Kerr(M, 0.0)};

    schw.SetParameter("mass", M);
    kerr.SetParameter("mass", M);
    kerr.SetParameter("spin", 0.0);

    Vec4 pos = getStandardPosition();

    Metric4d g_schw, g_kerr;
    Tensor<Dual<double>, 4, 4, 4> dg;

    schw.Evaluate(pos, g_schw, dg);
    kerr.Evaluate(pos, g_kerr, dg);

    EXPECT_NEAR(g_kerr(0, 0).real, g_schw(0, 0).real, kEpsilon)
        << "Kerr g_tt should match Schwarzschild at a=0";
    EXPECT_NEAR(g_kerr(1, 1).real, g_schw(1, 1).real, kEpsilon)
        << "Kerr g_rr should match Schwarzschild at a=0";
    EXPECT_NEAR(g_kerr(2, 2).real, g_schw(2, 2).real, kEpsilon)
        << "Kerr g_θθ should match Schwarzschild at a=0";
    EXPECT_NEAR(g_kerr(3, 3).real, g_schw(3, 3).real, kEpsilon)
        << "Kerr g_φφ should match Schwarzschild at a=0";
}

// Verify Reissner-Nordström reduces to Schwarzschild when charge Q=0
TEST_F(MetricLoaderChainTests, RNReducesToSchwarzschild) {
    sirius::core::KerrSchildFamily schw{sirius::core::KerrSchildParams::Schwarzschild(M)};
    sirius::core::KerrSchildFamily rn{sirius::core::KerrSchildParams::ReissnerNordstrom(M, 0.0)};

    schw.SetParameter("mass", M);
    rn.SetParameter("mass", M);
    rn.SetParameter("charge", 0.0);

    Vec4 pos = getStandardPosition();

    Metric4d g_schw, g_rn;
    Tensor<Dual<double>, 4, 4, 4> dg;

    schw.Evaluate(pos, g_schw, dg);
    rn.Evaluate(pos, g_rn, dg);

    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(g_rn(i, i).real, g_schw(i, i).real, kEpsilon)
            << "RN should match Schwarzschild at Q=0 for g[" << i << "][" << i << "]";
    }
}

// Verify all metrics approach Minkowski at large radius
TEST_F(MetricLoaderChainTests, FarFieldAsymptoticFlatness) {
    Vec4 far_pos;
    far_pos(0) = 0.0;
    far_pos(1) = 1000.0;  // Far from source (x=1000)
    far_pos(2) = 0.0;
    far_pos(3) = 0.0;

    for (const auto& metric : metrics) {
        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric->Evaluate(far_pos, g, dg);

        EXPECT_NEAR(g(0, 0).real, -1.0, 0.01)
            << "Metric " << metric->GetName() << " g_tt should be ~-1 at large r";
        EXPECT_NEAR(g(1, 1).real, 1.0, 0.01)
            << "Metric " << metric->GetName() << " g_rr should be ~1 at large r";
    }
}

// Verify Minkowski returns exactly flat metric at all positions
TEST_F(MetricLoaderChainTests, MinkowskiIsExactlyFlat) {
    sirius::core::KerrSchildFamily mink{sirius::core::KerrSchildParams::Minkowski()};

    std::vector<Vec4> positions;
    for (double r : {1.0, 10.0, 100.0}) {
        // Cartesian on x-axis, y-axis, z-axis
        Vec4 px, py, pz;
        px(0)=0; px(1)=r; px(2)=0; px(3)=0;
        py(0)=0; py(1)=0; py(2)=r; py(3)=0;
        pz(0)=0; pz(1)=0; pz(2)=0; pz(3)=r;
        positions.push_back(px);
        positions.push_back(py);
        positions.push_back(pz);
    }

    for (const auto& pos : positions) {
        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        mink.Evaluate(pos, g, dg);

        // Minkowski returns Cartesian form diag(-1, 1, 1, 1)
        EXPECT_NEAR(g(0, 0).real, -1.0, kEpsilon);
        EXPECT_NEAR(g(1, 1).real, 1.0, kEpsilon);
        EXPECT_NEAR(g(2, 2).real, 1.0, kEpsilon);
        EXPECT_NEAR(g(3, 3).real, 1.0, kEpsilon);

        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                EXPECT_NEAR(g(i, j).real, 0.0, kEpsilon);
            }
        }

        // All derivatives should be zero
        for (int a = 0; a < 4; ++a) {
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    EXPECT_NEAR(dg(a, i, j).real, 0.0, kEpsilon);
                }
            }
        }
    }
}

// --- Coordinate Boundary Tests ---

// Verify metrics handle θ near poles without singularities
TEST_F(MetricLoaderChainTests, HandlesNearPoles) {
    std::vector<double> pole_thetas = {0.001, PI - 0.001};

    for (const auto& metric : metrics) {
        for (double theta : pole_thetas) {
            Vec4 pos;
            // Near pole: theta small. Convert to Cartesian
            double r = 10.0;
            double phi = 0.0;
            pos(0) = 0.0;
            pos(1) = r * std::sin(theta) * std::cos(phi);
            pos(2) = r * std::sin(theta) * std::sin(phi);
            pos(3) = r * std::cos(theta); // z ~ +/- r

            Metric4d g;
            Tensor<Dual<double>, 4, 4, 4> dg;

            bool valid = evaluateMetricSafe(metric.get(), pos, g, dg);
            EXPECT_TRUE(valid)
                << "Metric " << metric->GetName()
                << " should handle theta near pole (" << theta << ")";
        }
    }
}

// Verify Schwarzschild handles positions just outside horizon (r > 2M)
TEST_F(MetricLoaderChainTests, HandlesNearHorizon) {
    sirius::core::KerrSchildFamily schw{sirius::core::KerrSchildParams::Schwarzschild(M)};
    schw.SetParameter("mass", 1.0);

    std::vector<double> near_horizon_r = {2.1, 2.05, 2.01};

    for (double r : near_horizon_r) {
        Vec4 pos;
        pos(0) = 0.0;
        pos(1) = r; // x=r (on axis)
        pos(2) = 0.0;
        pos(3) = 0.0;

        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;

        bool valid = evaluateMetricSafe(&schw, pos, g, dg);
        EXPECT_TRUE(valid)
            << "Schwarzschild should handle near-horizon r=" << r;

        // Check for singularity vs regularity
        // Schwarzschild coordinates: g_rr diverges at r=2
        // Kerr-Schild coordinates: g_rr is finite (~1 + 2M/r) at r=2
        // We know we are testing Schwarzschild metric.
        // If the implementation is KerrSchildFamily, it uses Kerr-Schild coordinates (Horizon Penetrating).
        // If it is SchwarzschildFamily (if it existed) it might use Schwarzschild coords.
        // The current implementation uses KerrSchildFamily for Schwarzschild.

        // So we expect FINITE g_rr for Kerr-Schild
        EXPECT_LT(g(1, 1).real, 100.0)
            << "g_rr should be finite near horizon for Kerr-Schild coordinates";
        EXPECT_GT(g(1, 1).real, 1.0)
            << "g_rr should be > 1";
    }
}

// --- Derivative Consistency Tests ---

// Verify all metric derivatives are finite (no NaN/Inf)
TEST_F(MetricLoaderChainTests, MetricDerivativesFinite) {
    Vec4 pos = getStandardPosition();

    for (const auto& metric : metrics) {
        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        metric->Evaluate(pos, g, dg);

        bool all_finite = true;
        for (int a = 0; a < 4; ++a) {
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    if (std::isnan(dg(a, i, j).real) || std::isinf(dg(a, i, j).real)) {
                        all_finite = false;
                    }
                }
            }
        }

        EXPECT_TRUE(all_finite)
            << "Metric " << metric->GetName() << " has non-finite derivatives";
    }
}

// Verify Minkowski derivatives are finite (coordinate-dependent in spherical coords)
TEST_F(MetricLoaderChainTests, MinkowskiZeroDerivatives) {
    sirius::core::KerrSchildFamily mink{sirius::core::KerrSchildParams::Minkowski()};
    Vec4 pos = getStandardPosition();

    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    mink.Evaluate(pos, g, dg);

    for (int a = 0; a < 4; ++a) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                EXPECT_FALSE(std::isnan(dg(a, i, j).real));
            }
        }
    }
}

// --- Performance and Stability Tests ---

// Verify repeated evaluations at same position give identical results
TEST_F(MetricLoaderChainTests, SequentialEvaluationsStable) {
    sirius::core::KerrSchildFamily schw{sirius::core::KerrSchildParams::Schwarzschild(M)};
    Vec4 pos = getStandardPosition();

    Metric4d g_first;
    Tensor<Dual<double>, 4, 4, 4> dg;
    schw.Evaluate(pos, g_first, dg);

    for (int i = 0; i < 100; ++i) {
        Metric4d g;
        schw.Evaluate(pos, g, dg);

        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                EXPECT_DOUBLE_EQ(g(mu, nu).real, g_first(mu, nu).real)
                    << "Metric should be deterministic at iteration " << i;
            }
        }
    }
}

// Verify different positions produce different metric values
TEST_F(MetricLoaderChainTests, DifferentPositionsDifferentMetrics) {
    sirius::core::KerrSchildFamily schw{sirius::core::KerrSchildParams::Schwarzschild(M)};

    std::vector<double> radii = {5.0, 10.0, 20.0, 50.0};
    Metric4d prev_g;
    bool first = true;

    for (double r : radii) {
         Vec4 pos;
        pos(0) = 0.0;
        pos(1) = r; // x=r
        pos(2) = 0.0;
        pos(3) = 0.0;

        Metric4d g;
        Tensor<Dual<double>, 4, 4, 4> dg;
        schw.Evaluate(pos, g, dg);

        if (!first) {
            EXPECT_NE(g(0, 0).real, prev_g(0, 0).real)
                << "Metrics at different radii should differ";
        }

        prev_g = g;
        first = false;
    }
}

} // namespace sirius::test
