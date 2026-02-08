# Design Rationale

This document explains why Sirius is built the way it is: the physical assumptions it rests on, the computational strategies it employs, and the engineering constraints that shaped each choice. None of this is required reading for operating the renderer. Anyone modifying it will benefit from understanding the reasoning behind the decisions they inherit.

---

## Physical Assumptions

### The Test Particle Approximation

Sirius traces photons through a fixed spacetime geometry. The photons do not affect the geometry; they are test particles, probes that move through a gravitational field without perturbing it. This is standard in astrophysical ray tracing: the energy of a photon is negligible compared to the mass-energy of a black hole or other compact object.

The approximation simplifies the problem enormously. Rather than solving Einstein's field equations (a coupled system of nonlinear partial differential equations describing how matter and geometry evolve together), we solve the geodesic equation, an ordinary differential equation, on a pre-specified background. This is the difference between numerical relativity and geodesic ray tracing. Sirius does the latter.

### Null Geodesics

Light travels along **null geodesics**: paths through spacetime where the spacetime interval vanishes. The **metric tensor** $g_{\mu\nu}$ defines the geometry at each point. The **geodesic equation**

$$
\frac{d^2 x^\mu}{d\lambda^2} + \Gamma^\mu_{\alpha\beta} \frac{dx^\alpha}{d\lambda} \frac{dx^\beta}{d\lambda} = 0
$$

describes how a ray's trajectory bends in response to curvature, encoded in the **Christoffel symbols** $\Gamma^\mu_{\alpha\beta}$ (connection coefficients derived from the metric tensor). Here $x^\mu$ denotes the spacetime coordinates and $\lambda$ is an affine parameter along the ray.

The **null condition**

$$
g_{\mu\nu} k^\mu k^\nu = 0
$$

(where $k^\mu = dx^\mu/d\lambda$ is the tangent vector to the ray) must hold throughout the integration. If it drifts, the ray no longer represents a photon.

---

## Computational Strategy

### Why Analytic Christoffel Symbols

This is the single most consequential design decision in Sirius.

The Christoffel symbols appear in every evaluation of the geodesic equation. A typical ray requires hundreds of integration steps; the renderer evaluates millions of rays per frame. Any error in the Christoffel symbols propagates through every step and every ray.

**Numerical differentiation** of the metric tensor, the obvious low-effort approach, introduces truncation error proportional to the finite difference step size. This error is small on any single evaluation, but it accumulates across hundreds of steps. For a ray spiralling near the **photon sphere** (the surface of unstable circular photon orbits, at $r = 3M$ for a Schwarzschild black hole), where the geometry is steepest and the ray spends the most time, the accumulated error produces visible artefacts.

Sirius instead computes Christoffel symbols analytically, using **dual number automatic differentiation**. A dual number has the form $a + b\varepsilon$ where $\varepsilon^2 = 0$; evaluating a function with dual number inputs yields both the function value and its exact derivative simultaneously. The metric function is evaluated with dual number coordinates; the dual parts carry the derivatives $\partial_\mu g_{\nu\rho}$ needed for the Christoffel formula. This eliminates truncation error entirely for analytic spacetimes, confining numerical error to the integrator alone.

The cost is that every metric implementation must be written in terms of arithmetic operations on a templated scalar type (so that it works with both `double` and `Dual<double>`). This constrains the implementation style. The payoff is correct derivatives by construction.

For spacetimes specified on a numerical grid (from numerical relativity simulations, for instance), analytic differentiation is not available. In those cases, Sirius precomputes Christoffel symbols on the grid and uses high-order interpolation to reconstruct values at arbitrary points, trading memory for accuracy.

### Numerical Stability

A computation that produces a number is not necessarily a computation that produced the right number.

Sirius monitors mathematical invariants during integration:

- The **null constraint** $g_{\mu\nu} k^\mu k^\nu = 0$ must hold throughout. Tolerance: $< 10^{-6}$.
- In **stationary spacetimes** (those with a time-translation symmetry), the **energy** $E = -g_{t\mu} k^\mu$ is conserved along geodesics.
- In **axisymmetric spacetimes** (those with a rotational symmetry), the **angular momentum** $L = g_{\phi\mu} k^\mu$ is conserved along geodesics.
- Conservation law drift tolerance: $< 10^{-4}$ relative.

These conservation laws follow from Noether's theorem applied to the spacetime symmetries. They provide independent verification that the integrator is performing correctly; each check costs little to compute and catches errors that would otherwise be invisible.

When a ray's invariants drift beyond tolerance, Sirius terminates it and marks it as failed. Continuing the integration and producing a pixel that looks plausible but is wrong would be worse than an obvious error. A black pixel where a star should be tells the developer something went wrong. A subtly displaced star tells no one anything.

### Determinism

Given identical inputs, Sirius produces identical outputs. Across runs, across platforms, across builds.

This is an engineering requirement, not an aspiration. Without determinism:

- Regression tests cannot use exact comparison.
- Debugging by bisection becomes unreliable (the same commit may produce different outputs on successive runs).
- Rendered images cannot be cited in scientific publications, because the next run might differ.

Determinism in floating-point code requires sustained discipline. Reduction order must be fixed: parallel summation must follow a deterministic schedule. Random number generators must be seeded. Fused multiply-add (FMA) settings must be consistent across compilation units. These constraints are specified in the [Specification](specification.md) and enforced through testing.

---

## Design Principles

**Mathematical traceability.** Every computation traces to established physics. The geodesic equation, metric tensor expressions, and Christoffel symbols derive from standard general relativity (Misner, Thorne & Wheeler 1973; Wald 1984). Where an approximation or simplification is introduced, the documentation records what was assumed and where the approximation breaks down.

**Verifiability.** Every mathematical calculation has a corresponding test. Known analytic results (the Schwarzschild photon sphere at $r = 3M$, the ISCO at $r = 6M$, the weak-field deflection angle $\Delta\phi = 4M/b$) serve as reference values. Conservation laws provide continuous monitoring during integration. The target test-to-calculation ratio is 1.0: one test for every calculation.

**Bounded error.** Numerical operations have defined precision bounds. Error accumulation is modelled and monitored. When errors exceed tolerance, the system signals failure rather than producing unreliable results. The precision requirements are specified in the [Specification](specification.md).

**Simplicity over cleverness.** Complexity hides bugs. Straightforward implementations of well-understood algorithms are preferred over optimisations that obscure correctness. When optimisation is necessary, correctness testing precedes and accompanies it. Three readable lines are better than one clever one.

---

## Accuracy and Rendering

Sirius is designed for scientific visualisation. The default settings prioritise physical accuracy over rendering speed. When a trade-off exists between fidelity and performance, the default favours fidelity; lower-precision modes are available for interactive exploration.

General relativity produces inherently striking phenomena: the warped reflections of Einstein rings, the asymmetric Doppler brightening of a relativistic accretion disk, double images of stars behind a wormhole throat. Physical accuracy and visual impact reinforce each other; the geometry does the work.

---

## References

- Misner, C. W., Thorne, K. S., & Wheeler, J. A. (1973). *Gravitation*. W. H. Freeman.
- Wald, R. M. (1984). *General Relativity*. University of Chicago Press.
- Carroll, S. M. (2004). *Spacetime and Geometry*. Addison-Wesley.
- James, O., von Tunzelmann, E., Franklin, P., & Thorne, K. S. (2015). "Gravitational lensing by spinning black holes in astrophysics, and in the movie Interstellar." *Classical and Quantum Gravity*.
- Luminet, J.-P. (1979). "Image of a spherical black hole with thin accretion disk." *Astronomy and Astrophysics*.
- Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes*. Cambridge University Press.
