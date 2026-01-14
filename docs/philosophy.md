# Philosophy of Sirius

*On Spacetime, Light, and the Nature of Observation*

> "Space tells matter how to move; matter tells space how to curve."
> — John Archibald Wheeler

## Preface

This document explores the philosophical foundations of Sirius. Reading it is not required for operating the renderer; however, it may illuminate the reasoning behind design decisions and the assumptions that underlie the system.

Every scientific visualisation embeds a worldview. Assumptions about physics, computation, and observation shape the code as surely as the mathematics. Most rendering software conceals these assumptions beneath layers of abstraction. Sirius makes its philosophy explicit, rendering it available for scrutiny and refinement.

---

## Part I: The Nature of Spacetime

### What Is Spacetime?

Spacetime is the arena of physics. It is not a passive backdrop upon which events occur; it is an active participant that shapes the motion of matter and the propagation of light. General relativity, formulated by Einstein in 1915, reveals that spacetime is a four-dimensional pseudo-Riemannian manifold whose curvature encodes gravity.

The geometry of spacetime is described by the metric tensor $g_{\mu\nu}$, a symmetric tensor field that defines distances and angles at each point. The signature $(-, +, +, +)$ distinguishes time from space, ensuring that light travels along null geodesics and massive particles along timelike geodesics.

### Curvature and Gravity

In general relativity, gravity is not a force but a manifestation of curved geometry. A freely falling particle follows a geodesic, the straightest possible path through curved spacetime. The geodesic equation:

$$
\frac{d^2 x^\mu}{d\lambda^2} + \Gamma^\mu_{\alpha\beta} \frac{dx^\alpha}{d\lambda} \frac{dx^\beta}{d\lambda} = 0
$$

describes how the particle's worldline bends in response to the Christoffel symbols $\Gamma^\mu_{\alpha\beta}$, which encode the curvature of spacetime derived from the metric tensor.

This geometric perspective transforms our understanding of phenomena near compact objects. Light bending around a black hole is not deflected by a force; it follows the natural geometry of spacetime. The apparent bending is a consequence of Euclidean intuition failing in curved geometry.

### The Test Particle Approximation

Sirius operates in the test particle regime: the photons traced through spacetime are assumed to have negligible backreaction on the geometry. The spacetime is fixed; the light moves through it. This approximation is valid for astronomical observation, where photon energies are insignificant compared to the mass-energy of the sources.

---

## Part II: Light and Observation

### Null Geodesics

Light travels along null geodesics, paths for which the spacetime interval vanishes:

$$
g_{\mu\nu} k^\mu k^\nu = 0
$$

where $k^\mu$ is the tangent vector to the light ray. This condition must be preserved throughout the numerical integration of the geodesic equation. Sirius enforces this constraint through explicit renormalisation, maintaining the null condition to machine precision.

### Gravitational Lensing

The bending of light near massive objects produces gravitational lensing, one of the most striking predictions of general relativity. Strong lensing near black holes creates multiple images when a background source aligns with the gravitational lens. At the photon sphere ($r = 3M$ for Schwarzschild geometry), light can orbit unstably, producing the characteristic dark silhouette of a black hole shadow.

Sirius traces rays backward from the observer to their sources, naturally capturing all lensing effects. A ray that spirals close to the photon sphere may orbit multiple times before escaping, encoding information about the geometry in its final direction.

### Redshift and Blueshift

Light propagating through curved spacetime experiences gravitational redshift as it climbs out of a potential well and blueshift as it falls in. For an observer at infinity receiving light from a source near a gravitational mass, the observed frequency $\nu_{obs}$ relates to the emitted frequency $\nu_{emit}$ through:

$$
\frac{\nu_{obs}}{\nu_{emit}} = \sqrt{\frac{g_{tt}(r_{emit})}{g_{tt}(r_{obs})}}
$$

Sirius computes these frequency shifts to render accretion disks with physically accurate temperature distributions and spectral properties.

---

## Part III: Computational Philosophy

### Why Analytic Christoffel Symbols?

The Christoffel symbols appear in every evaluation of the geodesic equation. Numerical differentiation of the metric tensor introduces discretisation error that accumulates over thousands of integration steps per ray. For a renderer evaluating millions of rays per frame, small errors compound into visible artefacts.

Sirius computes Christoffel symbols analytically from closed-form expressions derived from the metric tensor. This eliminates numerical differentiation error entirely for analytic spacetimes, ensuring that the dominant source of error is the integrator itself rather than derivative approximation.

For numerical spacetimes (such as those from numerical relativity simulations), Sirius precomputes and caches Christoffel symbols on a grid, using high-order interpolation to reconstruct values at arbitrary points. This trades memory for accuracy, with the goal of achieving comparable precision to analytic metrics.

### Numerical Stability as Truth

A computation that produces a result is not necessarily correct. Sirius treats numerical stability as a measure of truthfulness: a result is trustworthy only if the mathematical invariants of the underlying physics are preserved.

The null constraint $g_{\mu\nu} k^\mu k^\nu = 0$ must be maintained throughout integration. Energy and angular momentum conservation (for stationary and axisymmetric spacetimes respectively) provide additional verification. Sirius monitors these invariants and terminates rays that exceed tolerance, preferring honesty about failure over plausible-looking garbage.

### Determinism and Reproducibility

Scientific visualisation must be reproducible. Given identical inputs, Sirius produces identical outputs across runs, platforms, and builds. This determinism enables:

- Debugging through bisection
- Regression testing with exact comparisons
- Scientific citation of rendered images

The commitment to determinism requires care with floating-point operations, thread ordering, and random number generation. Sirius documents these constraints and enforces them through testing.

---

## Part IV: Design Principles

### Mathematical Rigor

Every computation must be traceable to established physics. The geodesic equation, metric tensor expressions, and Christoffel symbols derive from textbook general relativity. Approximations and simplifications are documented with their assumptions and limits.

### Verifiability

All calculations must have corresponding validation tests. Mathematical invariants (null constraint, conservation laws, metric symmetry) are tested explicitly. Known analytic results (photon sphere radius, deflection angle, orbital periods) provide reference values for validation.

### Bounded Certainty

Numerical operations have defined precision bounds. Error accumulation is modelled and monitored. When errors exceed tolerance, Sirius signals failure rather than producing unreliable results.

### Simplicity

Complexity is a liability. Sirius prefers straightforward implementations of well-understood algorithms over clever optimisations that obscure correctness. When optimisation is necessary, correctness testing precedes and accompanies it.

---

## Part V: The Role of Visualisation

### Making the Invisible Visible

General relativity predicts phenomena that challenge human intuition. Black holes bend light so severely that they possess no visible surface, only a shadow against the background sky. Wormholes connect distant regions of spacetime through geometric tunnels. Warp drives contract space ahead and expand it behind.

These phenomena are not directly observable in everyday experience. Sirius serves to make them visible, translating the mathematics of differential geometry into images that can be comprehended, studied, and shared. The renderer is a telescope into theoretical spacetimes.

### Scientific Accuracy

Sirius is designed for scientific visualisation, not entertainment. The default settings prioritise physical accuracy over aesthetic appeal. When trade-offs exist between speed and accuracy, the default favours accuracy; performance options are available for interactive exploration when warranted.

This does not preclude beautiful images. General relativity produces inherently striking phenomena: the warped reflections of Einstein rings, the asymmetric brightness of relativistic accretion disks, the double imaging of stars behind a wormhole throat. Accuracy and beauty are not opposed; physical truth has its own aesthetic.

---

## Conclusion

Sirius rests on philosophical commitments:

Spacetime is geometry. Light follows the geometry. The renderer traces that geometry to produce images.

Computation must preserve mathematical truth. Invariants are tested, errors are bounded, and failures are acknowledged.

Determinism enables science. Reproducible results support verification, testing, and scientific communication.

Simplicity serves correctness. Complex code hides bugs; straightforward code reveals them.

Whether these principles produce accurate and useful visualisations is an empirical question that testing and comparison with observation can answer. The foundations are defensible, the mathematics sound, and the implementation rigorous.

*That is the philosophy of Sirius.*

---

## Further Reading

**General Relativity**

- Misner, C., Thorne, K., & Wheeler, J. (1973). *Gravitation*. W. H. Freeman.
- Wald, R. M. (1984). *General Relativity*. University of Chicago Press.
- Carroll, S. M. (2004). *Spacetime and Geometry*. Addison-Wesley.

**Black Hole Visualisation**

- James, O., von Tunzelmann, E., Franklin, P., & Thorne, K. S. (2015). "Gravitational lensing by spinning black holes in astrophysics, and in the movie Interstellar." *Classical and Quantum Gravity*.
- Luminet, J.-P. (1979). "Image of a spherical black hole with thin accretion disk." *Astronomy and Astrophysics*.

**Numerical Methods**

- Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes*. Cambridge University Press.
