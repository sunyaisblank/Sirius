# Mathematical Foundations

*From First Principles to Geodesic Ray Tracing*

> "The laws of physics must be of such a nature that they apply to systems of reference in any kind of motion."
> — Albert Einstein

## Overview

This document develops the mathematical theory underlying Sirius. We proceed from differential geometry through the geodesic equation to the specific spacetimes implemented in the renderer. Each section builds upon the previous; the reader who works through the material sequentially will acquire a complete understanding of the computational machinery.

Prerequisites: familiarity with multivariate calculus, linear algebra, and basic differential equations. Prior exposure to general relativity is helpful but not essential; the relevant concepts are developed from first principles.

---

## Part I: Differential Geometry

### 1.1 Manifolds and Coordinates

A **manifold** is a space that locally resembles Euclidean space. Spacetime is a four-dimensional manifold, meaning that near any point, we can assign four coordinates $(x^0, x^1, x^2, x^3)$ that smoothly parameterise the neighbourhood.

Coordinates are labels, not physical quantities. The choice of coordinates is arbitrary; physics must be independent of this choice. General relativity achieves this coordinate independence through tensorial formulations.

### 1.2 The Metric Tensor

The **metric tensor** $g_{\mu\nu}$ defines the geometry of spacetime. It is a symmetric $(0,2)$ tensor field that determines the spacetime interval between nearby events:

$$
ds^2 = g_{\mu\nu} \, dx^\mu dx^\nu
$$

where summation over repeated indices is implied (Einstein convention).

The Lorentzian signature $(-, +, +, +)$ distinguishes timelike intervals ($ds^2 < 0$), spacelike intervals ($ds^2 > 0$), and null intervals ($ds^2 = 0$). Light travels along null geodesics; massive particles move along timelike geodesics.

### 1.3 Covariant Derivative

The **covariant derivative** $\nabla_\mu$ generalises partial differentiation to curved spaces. For a vector field $V^\nu$:

$$
\nabla_\mu V^\nu = \partial_\mu V^\nu + \Gamma^\nu_{\mu\lambda} V^\lambda
$$

The Christoffel symbols $\Gamma^\nu_{\mu\lambda}$ encode how the coordinate basis vectors change across the manifold. They are not tensors; they transform in a specific way that ensures covariant derivatives of tensors are themselves tensors.

### 1.4 Christoffel Symbols

The **Christoffel symbols** (also called connection coefficients) are defined by:

$$
\Gamma^\lambda_{\mu\nu} = \frac{1}{2} g^{\lambda\sigma} \left( \partial_\mu g_{\sigma\nu} + \partial_\nu g_{\sigma\mu} - \partial_\sigma g_{\mu\nu} \right)
$$

This definition ensures metric compatibility ($\nabla_\lambda g_{\mu\nu} = 0$) and torsion-freeness ($\Gamma^\lambda_{\mu\nu} = \Gamma^\lambda_{\nu\mu}$).

For a 4-dimensional spacetime, there are 40 independent Christoffel symbols (64 components with the lower-index symmetry).

---

## Part II: The Geodesic Equation

### 2.1 Geodesics as Extremal Paths

A **geodesic** is a curve that extremises the proper length (for timelike curves) or satisfies the parallel transport condition (for null curves). Physically, geodesics represent the paths of freely falling particles.

The geodesic equation in affine parameterisation is:

$$
\frac{d^2 x^\mu}{d\lambda^2} + \Gamma^\mu_{\alpha\beta} \frac{dx^\alpha}{d\lambda} \frac{dx^\beta}{d\lambda} = 0
$$

where $\lambda$ is an affine parameter along the curve.

### 2.2 Null Geodesics

For light rays (massless particles), the tangent vector $k^\mu = dx^\mu/d\lambda$ satisfies the null condition:

$$
g_{\mu\nu} k^\mu k^\nu = 0
$$

This condition must be preserved by the geodesic equation. Analytically, it is preserved exactly; numerically, errors accumulate and explicit renormalisation is required.

### 2.3 First-Order Formulation

For numerical integration, we reformulate the second-order geodesic equation as a system of first-order equations. Define the state vector:

$$
\mathbf{y} = (x^0, x^1, x^2, x^3, k^0, k^1, k^2, k^3)
$$

The equations of motion are:

$$
\frac{dx^\mu}{d\lambda} = k^\mu
$$

$$
\frac{dk^\mu}{d\lambda} = -\Gamma^\mu_{\alpha\beta} k^\alpha k^\beta
$$

### 2.4 Conservation Laws

For spacetimes with symmetries, there exist conserved quantities along geodesics.

**Killing vectors**: If $\xi^\mu$ is a Killing vector (satisfying $\nabla_{(\mu} \xi_{\nu)} = 0$), then $g_{\mu\nu} \xi^\mu k^\nu$ is constant along geodesics.

For stationary spacetimes, $\partial_t$ is a Killing vector, giving conserved energy:

$$
E = -g_{t\mu} k^\mu
$$

For axisymmetric spacetimes, $\partial_\phi$ is a Killing vector, giving conserved angular momentum:

$$
L = g_{\phi\mu} k^\mu
$$

These conservation laws provide error monitors for numerical integration.

---

## Part III: Black Hole Spacetimes

### 3.1 Schwarzschild Metric

The Schwarzschild metric describes a static, spherically symmetric vacuum spacetime:

$$
ds^2 = -\left(1 - \frac{2M}{r}\right)dt^2 + \left(1 - \frac{2M}{r}\right)^{-1}dr^2 + r^2 d\Omega^2
$$

where $d\Omega^2 = d\theta^2 + \sin^2\theta \, d\phi^2$ and $M$ is the mass in geometric units ($G = c = 1$).

| Feature | Value | Significance |
|---------|-------|--------------|
| Event horizon | $r = 2M$ | Point of no return |
| Photon sphere | $r = 3M$ | Unstable circular photon orbits |
| ISCO | $r = 6M$ | Innermost stable circular orbit |

**Domain**: The exterior solution requires $r > 2M$.

### 3.2 Kerr Metric

The Kerr metric describes a stationary, axisymmetric vacuum spacetime (rotating black hole):

$$
ds^2 = -\left(1 - \frac{2Mr}{\Sigma}\right)dt^2 - \frac{4Mar\sin^2\theta}{\Sigma}dt\,d\phi + \frac{\Sigma}{\Delta}dr^2 + \Sigma \, d\theta^2 + \left(r^2 + a^2 + \frac{2Ma^2r\sin^2\theta}{\Sigma}\right)\sin^2\theta \, d\phi^2
$$

where:

- $\Sigma = r^2 + a^2\cos^2\theta$
- $\Delta = r^2 - 2Mr + a^2$
- $a = J/M$ is the spin parameter

| Feature | Value | Constraint |
|---------|-------|------------|
| Outer horizon | $r_+ = M + \sqrt{M^2 - a^2}$ | Requires $|a| \leq M$ |
| Inner horizon | $r_- = M - \sqrt{M^2 - a^2}$ | |
| Ergosphere | $r = M + \sqrt{M^2 - a^2\cos^2\theta}$ | |
| Ring singularity | $r = 0$, $\theta = \pi/2$ | Physical singularity |

**Domain**: Exterior solution requires $r > r_+$ and $|a| < M$ (sub-extremal).

### 3.3 Reissner-Nordström Metric

The Reissner-Nordström metric describes a static, spherically symmetric electrically charged black hole:

$$
ds^2 = -\left(1 - \frac{2M}{r} + \frac{Q^2}{r^2}\right)dt^2 + \left(1 - \frac{2M}{r} + \frac{Q^2}{r^2}\right)^{-1}dr^2 + r^2 d\Omega^2
$$

where $Q$ is the electric charge. The horizons are at:

$$
r_\pm = M \pm \sqrt{M^2 - Q^2}
$$

### 3.4 Kerr-Schild Form

The Kerr-Schild representation writes the metric as:

$$
g_{\mu\nu} = \eta_{\mu\nu} + H \, l_\mu l_\nu
$$

where $\eta_{\mu\nu}$ is the Minkowski metric, $H$ is a scalar function, and $l_\mu$ is a null vector field.

This form has computational advantages: the Christoffel symbols become polynomial expressions in Cartesian coordinates, eliminating coordinate singularities at the poles.

---

## Part IV: Exotic Geometries

### 4.1 Ellis Drainhole (Wormhole)

The Ellis drainhole is a traversable wormhole with metric:

$$
ds^2 = -dt^2 + d\ell^2 + (b_0^2 + \ell^2)(d\theta^2 + \sin^2\theta \, d\phi^2)
$$

where $\ell \in (-\infty, +\infty)$ is the proper radial coordinate and $b_0$ is the throat radius.

| Feature | Value |
|---------|-------|
| Throat | $\ell = 0$, radius $b_0$ |
| Topology | Two asymptotically flat regions connected |
| Horizon | None (traversable) |

### 4.2 Morris-Thorne Family

The general Morris-Thorne wormhole metric is:

$$
ds^2 = -e^{2\Phi(r)}dt^2 + \frac{dr^2}{1 - b(r)/r} + r^2 d\Omega^2
$$

where $\Phi(r)$ is the redshift function and $b(r)$ is the shape function. The throat is at the minimum of $r$ where $b(r_0) = r_0$.

### 4.3 Alcubierre Warp Drive

The Alcubierre metric describes a "warp bubble" moving faster than light:

$$
ds^2 = -dt^2 + (dx - v_s f(r_s) dt)^2 + dy^2 + dz^2
$$

where $v_s$ is the velocity of the bubble, $r_s$ is the distance from the bubble centre, and $f(r_s)$ is a shaping function.

---

## Part V: Numerical Integration

### 5.1 Fourth-Order Runge-Kutta (RK4)

The classical RK4 method for the system $\dot{\mathbf{y}} = \mathbf{f}(\mathbf{y})$:

$$
\mathbf{k}_1 = h \, \mathbf{f}(\mathbf{y}_n)
$$

$$
\mathbf{k}_2 = h \, \mathbf{f}(\mathbf{y}_n + \frac{1}{2}\mathbf{k}_1)
$$

$$
\mathbf{k}_3 = h \, \mathbf{f}(\mathbf{y}_n + \frac{1}{2}\mathbf{k}_2)
$$

$$
\mathbf{k}_4 = h \, \mathbf{f}(\mathbf{y}_n + \mathbf{k}_3)
$$

$$
\mathbf{y}_{n+1} = \mathbf{y}_n + \frac{1}{6}(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4)
$$

Local truncation error is $O(h^5)$; global error is $O(h^4)$.

### 5.2 Adaptive Stepping (DOPRI)

The Dormand-Prince method (DOPRI5) provides embedded error estimation by computing both fourth-order and fifth-order approximations. The step size is adjusted to maintain a target error tolerance:

$$
h_{new} = h_{old} \cdot \left(\frac{\epsilon_{tol}}{|\epsilon_{est}|}\right)^{1/5}
$$

with safety factors to prevent oscillation.

### 5.3 Null Constraint Preservation

Numerical errors cause the null condition $g_{\mu\nu} k^\mu k^\nu = 0$ to drift. Sirius renormalises the wavevector periodically:

Given spatial components $k^i$, solve for $k^t$ from:

$$
g_{tt}(k^t)^2 + 2g_{ti}k^t k^i + g_{ij}k^i k^j = 0
$$

Choose the future-directed root.

---

## Part VI: Radiative Transfer

### 6.1 Specific Intensity

The specific intensity $I_\nu$ measures energy per unit time, area, solid angle, and frequency. Along a geodesic with affine parameter $\lambda$:

$$
\frac{dI_\nu}{d\lambda} = j_\nu - \alpha_\nu I_\nu
$$

where $j_\nu$ is the emission coefficient and $\alpha_\nu$ is the absorption coefficient.

### 6.2 Relativistic Beaming

A source moving relative to the observer experiences relativistic beaming. The observed intensity scales as:

$$
I_{obs} = \delta^4 I_{emit}
$$

where $\delta$ is the Doppler factor.

### 6.3 Accretion Disk Models

For thin accretion disks, the Novikov-Thorne temperature profile is:

$$
T(r) \propto \left[\frac{M}{\dot{M}} \frac{1}{r^3} \mathcal{C}(r)\right]^{1/4}
$$

where $\mathcal{C}(r)$ encodes the relativistic corrections depending on the metric.

---

## Part VII: Automatic Differentiation

### 7.1 Dual Numbers

A **dual number** has the form $a + b\varepsilon$ where $\varepsilon^2 = 0$. Operations on dual numbers:

| Operation | Result |
|-----------|--------|
| $(a + b\varepsilon) + (c + d\varepsilon)$ | $(a+c) + (b+d)\varepsilon$ |
| $(a + b\varepsilon) \cdot (c + d\varepsilon)$ | $ac + (ad + bc)\varepsilon$ |
| $\sin(a + b\varepsilon)$ | $\sin a + b\cos a \cdot \varepsilon$ |
| $\sqrt{a + b\varepsilon}$ | $\sqrt{a} + \frac{b}{2\sqrt{a}}\varepsilon$ |

### 7.2 Forward-Mode AD

By evaluating a function with dual number inputs, we obtain both the function value and its derivative simultaneously. For $f: \mathbb{R}^n \to \mathbb{R}$ and input $\mathbf{x}$:

$$
f(\mathbf{x} + \mathbf{e}_i \varepsilon) = f(\mathbf{x}) + \frac{\partial f}{\partial x_i}\varepsilon
$$

This provides exact derivatives (to machine precision) without numerical differentiation error.

### 7.3 Application to Metrics

Sirius uses dual numbers to compute Christoffel symbols analytically. The metric tensor is evaluated with dual number coordinates; the dual parts give the derivatives $\partial_\mu g_{\nu\rho}$ needed for the Christoffel formula.

---

## Appendix A: Notation Reference

| Symbol | Meaning |
|--------|---------|
| $g_{\mu\nu}$ | Metric tensor (covariant) |
| $g^{\mu\nu}$ | Inverse metric tensor (contravariant) |
| $\Gamma^\lambda_{\mu\nu}$ | Christoffel symbols |
| $x^\mu$ | Spacetime coordinates |
| $k^\mu$ | Photon wavevector |
| $\lambda$ | Affine parameter |
| $M$ | Black hole mass (geometric units) |
| $a$ | Spin parameter ($a = J/M$) |
| $Q$ | Electric charge |

---

## Appendix B: Unit Conventions

Sirius uses **geometric units** where $G = c = 1$. In these units:

| Quantity | Dimension | SI Conversion |
|----------|-----------|---------------|
| Mass | Length | $M_{geom} = GM_{SI}/c^2$ |
| Time | Length | $t_{geom} = ct_{SI}$ |
| Angular momentum | Length² | $J_{geom} = GJ_{SI}/c^3$ |

---

## References

**General Relativity**

- Misner, C. W., Thorne, K. S., & Wheeler, J. A. (1973). *Gravitation*. W. H. Freeman.
- Wald, R. M. (1984). *General Relativity*. University of Chicago Press.

**Black Hole Geodesics**

- Bardeen, J. M. (1973). "Timelike and null geodesics in the Kerr metric." In *Black Holes*.
- Chandrasekhar, S. (1983). *The Mathematical Theory of Black Holes*. Oxford University Press.

**Numerical Methods**

- Hairer, E., Nørsett, S. P., & Wanner, G. (1993). *Solving Ordinary Differential Equations I*. Springer.
- Press, W. H. et al. (2007). *Numerical Recipes*. Cambridge University Press.

---

*End of Foundations*
