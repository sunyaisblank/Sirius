#pragma once

// Live-path polarisation transport (sirius::core): parallel transport of a
// polarisation 4-vector along the Cartesian Kerr-Schild null geodesics the CPU
// tracer integrates, and the map from a transported polarisation vector at the
// observer to the electric-vector position angle (EVPA) that rotates the
// Stokes Q/U pair. Specification E2 (docs/SPECIFICATION.md section 4).
//
// The oracle (oracle/polarisation_transport.h) validates this path against the
// Walker-Penrose complex constant in Boyer-Lindquist coordinates; here the same
// physics runs in the render-native Kerr-Schild Cartesian chart, so the
// transport equation is stated coordinate-freely and the Christoffel symbols
// come from the live metric family:
//   dx^mu/dl = k^mu,  df^mu/dl = -Gamma^mu_ab k^a f^b,   f.k = 0.
// Parallel transport preserves g_uv f^u f^v and g_uv f^u k^v exactly, so those
// inner products are the on-path conservation monitors.
//
// Precision follows core/geodesic_integrator: Core computes in double (Vec4 =
// Tensor<double,4>), per docs/STYLE.md section 7. Christoffel symbols are built
// from the family's g and dg through TensorOps, the same route the geodesic
// integrator uses. Reference for the transport-to-EVPA map: Connors, Piran &
// Stark (1980), ApJ 235, 224; Walker & Penrose (1970), Commun. Math. Phys. 18,
// 265.

#include "sirius/base/contracts.h"
#include "sirius/core/metrics/metric.h"
#include "sirius/core/polarisation/stokes.h"
#include "sirius/core/tensor.h"

#include <cmath>

namespace sirius::core {

// A null ray carrying a parallel-transported polarisation vector. position and
// velocity match core/geodesic_integrator's Lightray fields (Cartesian
// Kerr-Schild (t, x, y, z), contravariant k^mu); polarisation is the
// contravariant f^mu with f.k = 0.
struct PolarisedRay {
    Vec4 position;      // 4-position (t, x, y, z).
    Vec4 velocity;      // Contravariant tangent k^mu = dx^mu/dl.
    Vec4 polarisation;  // Contravariant polarisation f^mu, f.k = 0, spacelike.
    double affine = 0.0;
};

//==============================================================================
// Metric-projection helpers (live path uses metric inner products directly)
//==============================================================================

// -Gamma^mu_ab u^a v^b from a built Christoffel stack; shared by the geodesic
// (u = v = k) and transport (u = k, v = f) right-hand sides.
[[nodiscard]] inline Vec4 NegConnectionContraction(const ChristoffelSymbols& gamma, const Vec4& u,
                                                   const Vec4& v) {
    Vec4 out;
    for (int mu = 0; mu < 4; ++mu) {
        double s = 0.0;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) s += gamma.gamma(mu, a, b).real * u(a) * v(b);
        out(mu) = -s;
    }
    return out;
}

//==============================================================================
// ParallelTransportStep: one RK4 step of (x, k, f) on the live metric family
//==============================================================================

namespace detail {

struct PolarisedDeriv {
    Vec4 dx, dk, df;
};

// Right-hand side of the coupled system at a ray state. One metric evaluation
// and one Christoffel build serve both dk and df, so the tangent and the
// polarisation see an identical connection.
[[nodiscard]] inline PolarisedDeriv PolarisedRhs(IMetric& metric, const Vec4& x, const Vec4& k,
                                                 const Vec4& f) {
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric.Evaluate(x, g, dg);
    const ChristoffelSymbols gamma = TensorOps::Christoffel(g, dg);

    PolarisedDeriv d;
    d.dx = k;
    d.dk = NegConnectionContraction(gamma, k, k);
    d.df = NegConnectionContraction(gamma, k, f);
    return d;
}

}  // namespace detail

// Fourth-order Runge-Kutta step of size h advancing position, tangent, and
// polarisation together. RK4 (not the tracer's Dormand-Prince pair) because the
// polarisation carries no error-estimation state of its own; transporting k and
// f with one shared connection evaluation is what keeps f.k and f.f invariant to
// integrator order. Precondition: metric is non-null (enforced by the caller's
// contract). Postcondition: ray advanced by h in affine parameter.
inline void ParallelTransportStep(IMetric& metric, PolarisedRay& ray, double h) {
    const Vec4 x0 = ray.position;
    const Vec4 k0 = ray.velocity;
    const Vec4 f0 = ray.polarisation;

    const detail::PolarisedDeriv s1 = detail::PolarisedRhs(metric, x0, k0, f0);
    const detail::PolarisedDeriv s2 = detail::PolarisedRhs(
        metric, x0 + s1.dx * (0.5 * h), k0 + s1.dk * (0.5 * h), f0 + s1.df * (0.5 * h));
    const detail::PolarisedDeriv s3 = detail::PolarisedRhs(
        metric, x0 + s2.dx * (0.5 * h), k0 + s2.dk * (0.5 * h), f0 + s2.df * (0.5 * h));
    const detail::PolarisedDeriv s4 =
        detail::PolarisedRhs(metric, x0 + s3.dx * h, k0 + s3.dk * h, f0 + s3.df * h);

    ray.position = x0 + (s1.dx + s2.dx * 2.0 + s3.dx * 2.0 + s4.dx) * (h / 6.0);
    ray.velocity = k0 + (s1.dk + s2.dk * 2.0 + s3.dk * 2.0 + s4.dk) * (h / 6.0);
    ray.polarisation = f0 + (s1.df + s2.df * 2.0 + s3.df * 2.0 + s4.df) * (h / 6.0);
    ray.affine += h;
}

//==============================================================================
// Polarisation gauge fixing and inner products
//==============================================================================

// Orthonormalise a trial polarisation against a null tangent: f = trial -
// (trial.k / e_t.k) e_t with e_t = (1,0,0,0), giving f.k = 0 because k is null,
// then scale to g_uv f^u f^v = 1. The residual freedom f -> f + alpha k does not
// change the physical polarisation or the Walker-Penrose constant.
[[nodiscard]] inline Vec4 MakeOrthonormalPolarisation(IMetric& metric, const Vec4& x,
                                                      const Vec4& k, const Vec4& trial) {
    Metric4d g;
    Tensor<Dual<double>, 4, 4, 4> dg;
    metric.Evaluate(x, g, dg);

    Vec4 e_t;
    e_t(0) = 1.0;
    const double trial_dot_k = TensorOps::InnerProduct(trial, k, g);
    const double et_dot_k = TensorOps::InnerProduct(e_t, k, g);
    SIRIUS_PRE(std::abs(et_dot_k) > 0.0);

    Vec4 f = trial - e_t * (trial_dot_k / et_dot_k);
    const double norm2 = TensorOps::InnerProduct(f, f, g);
    SIRIUS_PRE(norm2 > 0.0);  // Spacelike polarisation.
    return f / std::sqrt(norm2);
}

//==============================================================================
// Transport -> EVPA -> Stokes Q/U rotation
//==============================================================================
//
// A linearly polarised ray's Stokes Q, U are set by the electric-vector
// position angle chi in a chosen sky basis (e_ref, e_perp) spanning the plane
// transverse to the line of sight in the observer frame: Q = p I cos 2chi,
// U = p I sin 2chi. Parallel transport rotates the physical polarisation
// vector, so chi at the observer differs from chi at emission; the difference
// rotates (Q, U). The projection uses metric inner products, and the rotation
// reuses the existing Stokes algebra (stokes.h MuellerMatrix::Rotation), so no
// second rotation implementation is introduced.

// EVPA of a polarisation vector f projected onto the observer sky basis:
// chi = atan2(f.e_perp, f.e_ref), in radians. e_ref and e_perp are the two
// spacelike screen vectors transverse to the ray in the observer frame.
[[nodiscard]] inline double ElectricVectorPositionAngle(const Vec4& f, const Vec4& e_ref,
                                                        const Vec4& e_perp, const Metric4d& g) {
    const double f_ref = TensorOps::InnerProduct(f, e_ref, g);
    const double f_perp = TensorOps::InnerProduct(f, e_perp, g);
    return std::atan2(f_perp, f_ref);
}

// Rotate a Stokes vector so its electric-vector position angle carries the
// transported polarisation direction. chi is the EVPA of f in the sky basis, and
// the result's EVPA is the emitted EVPA raised by chi (so reference-aligned
// emitted light acquires EVPA = chi, matching StokesVector::Evpa). I and V are
// unchanged. stokes.h MuellerMatrix::Rotation(theta) is a reference-frame
// rotation, which lowers the EVPA by theta, so Rotation(-chi) raises it by chi.
[[nodiscard]] inline StokesVector RotateStokesToObserverFrame(const StokesVector& emitted,
                                                              const Vec4& f_transported,
                                                              const Vec4& e_ref, const Vec4& e_perp,
                                                              const Metric4d& g) {
    const double chi = ElectricVectorPositionAngle(f_transported, e_ref, e_perp, g);
    return MuellerMatrix::Rotation(static_cast<float>(-chi)).Apply(emitted);
}

}  // namespace sirius::core
