// Double-precision beam (ray-bundle) integrator for the validation oracle:
// central geodesic and covariant Jacobi columns share one RK4 tableau. Drives
// the analytic beam-propagation benchmarks. Off the render path.
//   D^2 xi^mu / d lambda^2 + R^mu_nu_rho_sigma k^nu xi^rho k^sigma = 0
// Reference: James et al. (2015), Appendix B.
// Ported from INBI001A.h; numeric content bit-identical.

#pragma once

#include "sirius/base/contracts.h"
#include "sirius/oracle/metric_interface.h"
#include "sirius/oracle/transport_types.h"

#include <algorithm>
#include <cmath>
#include <numbers>

// CPU/GPU compatibility macro
#ifdef __CUDACC__
#define SIRIUS_HOST_DEVICE __host__ __device__
#else
#define SIRIUS_HOST_DEVICE
#endif

namespace sirius::oracle {

//==============================================================================
// BeamStateD: Complete state for beam propagation
//==============================================================================

struct BeamStateD {
    // Central geodesic state (64 bytes)
    Vec4d x;  // Position: (t, r, θ, φ)
    Vec4d k;  // Wave 4-vector (covariant): (k_t, k_r, k_θ, k_φ)

    // Jacobian matrix J^i_j = ∂x^i/∂x^j_0 (128 bytes)
    // Maps initial conditions to current position
    double J[4][4];

    // Covariant Jacobian velocity V^i_j = D J^i_j/dλ (128 bytes).
    // Storing V rather than the coordinate derivative is what makes the
    // Jacobi evolution tensorial in a non-Cartesian chart.
    double dJ[4][4];

    // Affine parameter
    double lambda;

    // Conserved quantities (computed at initialisation)
    double E;   // Energy: E = -k_t
    double Lz;  // Angular momentum: L_z = k_φ
    double Q;   // Carter constant (Kerr only)

    // Derived beam properties (cached)
    double solid_angle;    // Beam solid angle on sky
    double magnification;  // |det(J_angular)|^{-1}
    double major_axis;     // Ellipse semi-major axis [rad]
    double minor_axis;     // Ellipse semi-minor axis [rad]
    double orientation;    // PA of major axis [rad]

    // Beam status
    bool terminated;  // Hit horizon or escaped
    bool at_caustic;  // det(J) ≈ 0 (magnification → ∞)

    // Initial pixel solid angle (set at creation)
    double initial_pixel_solid_angle;

    //--------------------------------------------------------------------------
    // Initialisation
    //--------------------------------------------------------------------------

    SIRIUS_HOST_DEVICE void Initialise() {
        // Identity Jacobian (beam starts as a point)
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                J[i][j] = (i == j) ? 1.0 : 0.0;
                dJ[i][j] = 0.0;
            }
        }

        lambda = 0;
        solid_angle = 0;
        magnification = 1.0;
        major_axis = 0;
        minor_axis = 0;
        orientation = 0;
        terminated = false;
        at_caustic = false;
        initial_pixel_solid_angle = 0;
    }

    //--------------------------------------------------------------------------
    // Extract beam geometry from Jacobian
    //--------------------------------------------------------------------------

    SIRIUS_HOST_DEVICE void UpdateGeometry() {
        // Extract 2×2 angular submatrix (θ,φ) → (θ₀,φ₀) mapping
        // This gives the beam ellipse on the sky

        double a = J[2][2];  // ∂θ/∂θ₀
        double b = J[2][3];  // ∂θ/∂φ₀
        double c = J[3][2];  // ∂φ/∂θ₀
        double d = J[3][3];  // ∂φ/∂φ₀

        // Determinant of angular submatrix
        double det = a * d - b * c;

        // Check for caustic (det → 0)
        if (std::abs(det) < 1e-12) {
            at_caustic = true;
            magnification = 1e12;  // Cap at large value
        } else {
            at_caustic = false;
            magnification = 1.0 / std::abs(det);
        }

        // Singular value decomposition for ellipse axes
        // For matrix [[a,b],[c,d]], singular values are:
        // σ = √[(p ± √(p² - 4q²))/2] where p = a² + b² + c² + d², q = det

        double p = a * a + b * b + c * c + d * d;
        double q = det;
        double disc = p * p - 4 * q * q;

        if (disc < 0) disc = 0;  // Numerical protection
        double s = std::sqrt(disc);

        major_axis = std::sqrt(std::max(0.0, (p + s) / 2));
        minor_axis = std::sqrt(std::max(0.0, (p - s) / 2));

        // Orientation of the output ellipse: the major eigenvector of J J^T.
        // For J=[[a,b],[c,d]], tan(2φ)=2(ac+bd)/(a²+b²-c²-d²).
        // The former ab+cd expression is the right-singular-vector angle in the
        // input plane and does not orient an ellipse on the rendered sky.
        double num = 2 * (a * c + b * d);
        double den = a * a + b * b - c * c - d * d;
        orientation = 0.5 * std::atan2(num, den);

        // Solid angle = π × σ_max × σ_min × (initial pixel solid angle)
        solid_angle = std::numbers::pi * major_axis * minor_axis * initial_pixel_solid_angle;
    }

    //--------------------------------------------------------------------------
    // Convert from GeodesicStateD
    //--------------------------------------------------------------------------

    SIRIUS_HOST_DEVICE void FromGeodesic(const GeodesicStateD& geo) {
        x = geo.x;
        k = geo.k;
        lambda = geo.lambda;
        E = geo.E;
        Lz = geo.Lz;
        Q = geo.Q;
    }

    //--------------------------------------------------------------------------
    // Convert to GeodesicStateD
    //--------------------------------------------------------------------------

    SIRIUS_HOST_DEVICE GeodesicStateD ToGeodesic() const {
        GeodesicStateD geo;
        geo.x = x;
        geo.k = k;
        geo.lambda = lambda;
        geo.E = E;
        geo.Lz = Lz;
        geo.Q = Q;
        return geo;
    }
};

//==============================================================================
// BeamIntegratorD: Simultaneous geodesic and Jacobian integration
//==============================================================================

class BeamIntegratorD {
  public:
    struct Config {
        double escape_radius = 1e6;
    };

    [[nodiscard]] static bool IsRepresentedConfig(const Config& config) noexcept {
        return std::isfinite(config.escape_radius) && config.escape_radius > 0.0;
    }

    explicit BeamIntegratorD(const IMetricD* metric) : BeamIntegratorD(metric, Config()) {}

    BeamIntegratorD(const IMetricD* metric, const Config& config)
        : metric_(metric), config_(config) {
        SIRIUS_PRE(metric != nullptr);
        SIRIUS_PRE(IsRepresentedConfig(config));
    }

    //--------------------------------------------------------------------------
    // Single integration Step
    //--------------------------------------------------------------------------

    SIRIUS_HOST_DEVICE bool Step(BeamStateD& beam, double h) const {
        SIRIUS_PRE(std::isfinite(h) && h > 0.0);
        if (beam.terminated) return false;

        // Check if outside valid region
        if (!metric_->IsValid(beam.x)) {
            beam.terminated = true;
            return false;
        }

        // Check if escaped
        if (beam.x.r > config_.escape_radius) {
            beam.terminated = true;
            return false;
        }

        // Integrate the central Hamiltonian ray and every Jacobi column in the
        // same RK4 tableau. The former implementation leapfrogged the ray while
        // freezing curvature for an unrelated RK4 update, and treated dJ as an
        // ordinary derivative. In Boyer-Lindquist coordinates that drops the
        // connection terms and is not the covariant Jacobi equation.
        struct Derivative {
            Vec4d x;
            Vec4d k;
            double J[4][4];
            double V[4][4];
        };

        auto right_hand_side = [&](const Vec4d& x, const Vec4d& k, const double J[4][4],
                                   const double V[4][4], Derivative& out) {
            if (!metric_->IsValid(x)) {
                return false;
            }
            out.x = metric_->dHdp(x, k);
            out.k = metric_->dHdq(x, k);

            double Gamma[4][4][4];
            double R[4][4][4][4];
            metric_->Christoffel(x, Gamma);
            metric_->Riemann(x, R);

            for (int mu = 0; mu < 4; ++mu) {
                for (int column = 0; column < 4; ++column) {
                    double connection_j = 0.0;
                    double connection_v = 0.0;
                    for (int nu = 0; nu < 4; ++nu) {
                        for (int rho = 0; rho < 4; ++rho) {
                            connection_j += Gamma[mu][nu][rho] * out.x[nu] * J[rho][column];
                            connection_v += Gamma[mu][nu][rho] * out.x[nu] * V[rho][column];
                        }
                    }

                    double tidal = 0.0;
                    for (int nu = 0; nu < 4; ++nu) {
                        for (int rho = 0; rho < 4; ++rho) {
                            for (int sigma = 0; sigma < 4; ++sigma) {
                                tidal += R[mu][nu][rho][sigma] * out.x[nu] * J[rho][column] *
                                         out.x[sigma];
                            }
                        }
                    }

                    // V = D J / dλ, hence ordinary coordinate derivatives are:
                    // dJ/dλ = V - Γ(k,J), dV/dλ = -Γ(k,V) - R(k,J)k.
                    out.J[mu][column] = V[mu][column] - connection_j;
                    out.V[mu][column] = -connection_v - tidal;
                }
            }
            return true;
        };

        Derivative k1{}, k2{}, k3{}, k4{};
        double stage_J[4][4];
        double stage_V[4][4];
        SIRIUS_ASSERT(right_hand_side(beam.x, beam.k, beam.J, beam.dJ, k1));

        auto make_stage = [&](const Derivative& derivative, double scale) {
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    stage_J[i][j] = beam.J[i][j] + scale * derivative.J[i][j];
                    stage_V[i][j] = beam.dJ[i][j] + scale * derivative.V[i][j];
                }
            }
        };

        make_stage(k1, 0.5 * h);
        if (!right_hand_side(beam.x + k1.x * (0.5 * h), beam.k + k1.k * (0.5 * h), stage_J, stage_V,
                             k2)) {
            beam.terminated = true;
            return false;
        }

        make_stage(k2, 0.5 * h);
        if (!right_hand_side(beam.x + k2.x * (0.5 * h), beam.k + k2.k * (0.5 * h), stage_J, stage_V,
                             k3)) {
            beam.terminated = true;
            return false;
        }

        make_stage(k3, h);
        if (!right_hand_side(beam.x + k3.x * h, beam.k + k3.k * h, stage_J, stage_V, k4)) {
            beam.terminated = true;
            return false;
        }

        beam.x += (k1.x + k2.x * 2.0 + k3.x * 2.0 + k4.x) * (h / 6.0);
        beam.k += (k1.k + k2.k * 2.0 + k3.k * 2.0 + k4.k) * (h / 6.0);
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                beam.J[i][j] +=
                    (k1.J[i][j] + 2.0 * k2.J[i][j] + 2.0 * k3.J[i][j] + k4.J[i][j]) * (h / 6.0);
                beam.dJ[i][j] +=
                    (k1.V[i][j] + 2.0 * k2.V[i][j] + 2.0 * k3.V[i][j] + k4.V[i][j]) * (h / 6.0);
            }
        }
        beam.lambda += h;

        // Update beam geometry
        beam.UpdateGeometry();

        return true;
    }

  private:
    const IMetricD* metric_;
    Config config_;
};

}  // namespace sirius::oracle
