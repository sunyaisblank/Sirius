# Sirius Type System

> "The beginning of wisdom is the definition of terms."
> — Socrates

## Preface

Types are the vocabulary of Sirius. They express meaning, enforce constraints, and make code auditable. This document builds from fundamental ideas about values and identity, through numeric and coordinate primitives, and into the domain types used in relativistic ray tracing. It also carries the philosophical position that every type is a claim about reality; a careless type is a careless claim.

The goal is clarity: a type should state what a value means, how it behaves, and what it forbids. The goal is integrity: a type should make invalid states hard to represent and easy to detect.

---

## Part I: First Principles

### 1.1 Value, Identity, and State

**Value** is a datum that stands on its own. Examples include a metric tensor component, a coordinate value, and a step size.

**Identity** refers to an entity that persists across time. A ray has an identity; its position and wavevector are values that change as the ray propagates.

**State** is the collection of values that describe an identity at a given time. State must be explicit and traceable; hidden state erodes auditability.

### 1.2 Invariants and Boundaries

An invariant is a rule that must remain true for a type. The system enforces invariants at boundaries. Inside a validated boundary, code can assume validity and remain efficient.

Example invariants:

- A radial coordinate $r > r_{horizon}$.
- A polar angle $\theta \in (0, \pi)$.
- The null condition $g_{\mu\nu} k^\mu k^\nu = 0$.

### 1.3 Meaning and Commitment

A type is a commitment to meaning. It is a promise about units, coordinate system, and allowed operations. This commitment has physical weight in a scientific system; untrue promises can become visualisation errors that mislead interpretation.

---

## Part II: Primitive Types

### 2.1 Numeric Types

| Type | Usage | Rationale |
|------|-------|-----------|
| `float` | GPU computation, real-time rendering | Performance; ~7 decimal digits |
| `double` | CPU computation, high-precision integration | Accuracy; ~15 decimal digits |
| `int` | Indices, loop counters, discrete quantities | Exact integer arithmetic |
| `size_t` | Array sizes, memory allocation | Platform-sized unsigned |

### 2.2 Precision Trade-offs

GPU shader code uses `float` for performance. CPU integration code uses `double` for precision. The boundary between GPU and CPU code must handle precision conversion explicitly.

| Context | Precision | Tolerance |
|---------|-----------|-----------|
| GPU ray tracing | Single | $10^{-5}$ relative |
| CPU geodesic integration | Double | $10^{-12}$ relative |
| Conservation law monitoring | Double | $10^{-4}$ relative |

### 2.3 Boolean

Boolean values represent predicates and termination conditions:

```cpp
bool isTerminated = ray.status != RayStatus::Propagating;
bool isNullPreserved = std::abs(nullNorm) < tolerance;
```

---

## Part III: Coordinate Types

### 3.1 Coordinate Systems

Sirius uses multiple coordinate systems with explicit type distinctions:

| Coordinate System | Notation | Domain |
|-------------------|----------|--------|
| Boyer-Lindquist | $(t, r, \theta, \phi)$ | $r > r_+$, $\theta \in (0, \pi)$ |
| Cartesian (Kerr-Schild) | $(t, x, y, z)$ | All $\mathbb{R}^3$ except singularity |
| Spherical (numerical) | $(t, r, \theta, \phi)$ | Grid bounds |

### 3.2 Vec4

The fundamental 4-vector type represents spacetime positions and wavevectors:

```cpp
struct Vec4 {
    double t, x, y, z;  // or t, r, theta, phi depending on context
};
```

**Invariant**: The coordinate system must be tracked externally or by a companion type.

### 3.3 Coordinate Transformations

Transformations between coordinate systems are explicit functions:

```cpp
Vec4 boyerLindquistToCartesian(const Vec4& bl, double a);
Vec4 cartesianToBoyerLindquist(const Vec4& cart, double a);
```

---

## Part IV: Tensor Types

### 4.1 Metric Tensor

The metric tensor $g_{\mu\nu}$ is a symmetric $4 \times 4$ matrix:

```cpp
using Metric4D = std::array<std::array<double, 4>, 4>;
```

**Invariants**:

- Symmetry: $g_{\mu\nu} = g_{\nu\mu}$
- Lorentzian signature: exactly one negative eigenvalue
- Invertibility: $\det(g) \neq 0$

### 4.2 Inverse Metric

The inverse metric $g^{\mu\nu}$ satisfies $g^{\mu\alpha}g_{\alpha\nu} = \delta^\mu_\nu$:

```cpp
Metric4D invertMetric(const Metric4D& g);
```

### 4.3 Christoffel Symbols

The Christoffel symbols form a 3-index array with symmetry in the lower indices:

```cpp
using Christoffel4D = std::array<std::array<std::array<double, 4>, 4>, 4>;
// christoffel[lambda][mu][nu] = Γ^λ_μν
```

**Invariant**: $\Gamma^\lambda_{\mu\nu} = \Gamma^\lambda_{\nu\mu}$

---

## Part V: Metric Interface

### 5.1 IMetric Interface

All spacetime metrics implement a common interface:

```cpp
class IMetric {
public:
    virtual Metric4D evaluate(const Vec4& position) const = 0;
    virtual Metric4D evaluateInverse(const Vec4& position) const = 0;
    virtual Christoffel4D christoffel(const Vec4& position) const = 0;
};
```

### 5.2 Metric Families

Metrics are organised into parameterised families:

| Family | Parameters | Members |
|--------|------------|---------|
| Kerr-Schild | $(M, a, Q, \Lambda)$ | Schwarzschild, Kerr, RN, Kerr-Newman, de Sitter |
| Morris-Thorne | $(b_0, \Phi)$ | Ellis, general wormholes |
| Warp Drive | $(v_s, R, \sigma)$ | Alcubierre, variants |

### 5.3 Parameter Bounds

Each metric family defines valid parameter ranges:

| Parameter | Range | Constraint |
|-----------|-------|------------|
| Mass $M$ | $[0.1, 100]$ | Practical rendering scale |
| Spin $a/M$ | $[0, 0.998]$ | Sub-extremal |
| Charge $Q/M$ | $[0, 0.999]$ | Sub-extremal |
| Throat $b_0$ | $> 0$ | Physical wormhole |

---

## Part VI: Ray Types

### 6.1 Ray State

A ray is characterised by its position, wavevector, and integration state:

```cpp
struct RayState {
    Vec4 position;      // Spacetime coordinates
    Vec4 wavevector;    // Null tangent vector k^μ
    double lambda;      // Affine parameter
    RayStatus status;   // Propagating, Terminated, Escaped, etc.
};
```

### 6.2 Ray Status

```cpp
enum class RayStatus {
    Propagating,    // Still integrating
    Escaped,        // Reached far boundary
    Captured,       // Fell into horizon
    HitDisk,        // Intersected accretion disk
    Terminated,     // Integration limit reached
    Error           // Numerical failure
};
```

### 6.3 Termination Conditions

| Condition | Detection | Action |
|-----------|-----------|--------|
| Escape | $r > r_{max}$ | Map to background texture |
| Capture | $r < r_{horizon} + \epsilon$ | Return black |
| Disk intersection | $\theta \approx \pi/2$ within disk radius | Compute emission |
| Step limit | steps $>$ maxSteps | Mark as terminated |
| NaN detection | `isnan(position)` or `isnan(wavevector)` | Mark as error |

---

## Part VII: Dual Number Types

### 7.1 Dual Number Template

```cpp
template<typename T>
struct Dual {
    T real;
    T dual;
    
    // Arithmetic operations preserving ε² = 0
};
```

### 7.2 Derivative Computation

Dual numbers enable automatic differentiation:

```cpp
Dual<double> x(r, 1.0);  // dx/dr = 1
auto result = f(x);       // result.dual = df/dr
```

### 7.3 Metric Derivatives

The Christoffel symbols require metric derivatives:

$$
\partial_\mu g_{\nu\rho}
$$

Using dual numbers with seeded directions provides these derivatives exactly.

---

## Part VIII: Observer Types

### 8.1 Observer State

```cpp
struct Observer {
    Vec4 position;      // Worldline position
    Vec4 velocity;      // 4-velocity (timelike, normalised)
    Tetrad frame;       // Local orthonormal frame
};
```

### 8.2 Tetrad

A tetrad is an orthonormal frame adapted to the observer:

```cpp
struct Tetrad {
    Vec4 e0;  // Timelike, parallel to 4-velocity
    Vec4 e1;  // Spacelike
    Vec4 e2;  // Spacelike
    Vec4 e3;  // Spacelike
};
```

**Invariants**:

- $g_{\mu\nu} e_a^\mu e_b^\nu = \eta_{ab}$ (Minkowski metric)
- $e_0$ is future-directed timelike

---

## Part IX: Result and Error Types

### 9.1 Failure Representation

Operations that can fail return explicit result types:

```cpp
std::optional<Metric4D> tryEvaluateMetric(const Vec4& position);
```

### 9.2 Error Conditions

| Error | Cause | Detection |
|-------|-------|-----------|
| OutOfDomain | Position outside valid region | Coordinate bounds check |
| SingularMetric | Determinant near zero | $|\det(g)| < \epsilon$ |
| NaNProduced | Invalid operation | `isnan()` check |
| ConvergenceFailed | Iteration limit exceeded | Iteration counter |

---

## Part X: Type Rules

These rules align with the coding standard:

| Rule | Level | Statement |
|------|-------|-----------|
| Precision selection | MUST | GPU uses `float`; CPU integration uses `double` |
| Boundary validation | MUST | Coordinates validated at metric entry |
| Invariant preservation | MUST | Null constraint enforced throughout integration |
| Coordinate explicitness | SHOULD | Coordinate system tracked with values |
| Failure expression | MUST | Failures return `optional` or set status; silence is forbidden |
| Immutable values | SHOULD | Metric evaluations return new objects |

---

## Appendix A: Type Naming Conventions

- `Vec4`: 4-component vector (position or velocity)
- `Metric4D`: $4 \times 4$ symmetric tensor
- `Christoffel4D`: $4 \times 4 \times 4$ connection array
- `Dual<T>`: Dual number with real type `T`
- `IMetric`: Abstract metric interface

---

## Appendix B: Checklist

Before adding or modifying a type, verify:

- [ ] The name expresses meaning and coordinates.
- [ ] Invariants are enforced at the boundary.
- [ ] Precision is appropriate for the context.
- [ ] Coordinate system is explicit.
- [ ] Failure modes are encoded in the return type.

---

*End of Type System*
