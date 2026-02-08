# Specification

This document defines the performance targets, mathematical invariants, precision requirements, and operational constraints for Sirius. Every requirement here is testable; the test suite verifies compliance.

---

---

## 2. Scope

Sirius covers four areas: metric evaluation (analytic and numerical spacetimes), geodesic integration (null geodesic propagation with adaptive stepping), radiative transfer (emission and absorption for accretion disk rendering), and visualisation (real-time and offline image generation).

Explicitly outside scope: numerical relativity (spacetime evolution), gravitational wave generation, and quantum effects (Hawking radiation).

---

## 3. Performance Requirements

### 3.1 Frame Rate Targets

| Metric Type | Target FPS | Minimum FPS | Resolution |
|-------------|------------|-------------|------------|
| Minkowski (flat) | 120 | 60 | 1920×1080 |
| Schwarzschild | 60 | 15 | 1920×1080 |
| Kerr | 30 | 10 | 1920×1080 |
| Numerical (ET) | 10 | 1 | 1920×1080 |

### 3.2 Integration Budget

| Metric | Steps per Ray | Reads per Step |
|--------|---------------|----------------|
| Analytic | 500 (target) | 0 (computed) |
| Numerical | 500 (target) | 1 (cached Christoffel) |

### 3.3 Latency Requirements

| Operation | Maximum Latency |
|-----------|-----------------|
| Metric evaluation | 1 μs |
| Christoffel computation (analytic) | 5 μs |
| Single integration step | 10 μs |
| Full ray trace | 10 ms |

---

## 4. Mathematical Invariants

### 4.1 Fundamental Invariants

These are the mathematical properties that must hold for a computation to be considered correct. The following MUST be preserved throughout integration:

| Invariant | Mathematical Form | Tolerance |
|-----------|-------------------|-----------|
| Null condition | $g_{\mu\nu} k^\mu k^\nu = 0$ | $< 10^{-6}$ |
| Energy conservation | $\frac{d}{d\lambda}(-g_{t\mu} k^\mu) = 0$ | $< 10^{-4}$ relative |
| Angular momentum | $\frac{d}{d\lambda}(g_{\phi\mu} k^\mu) = 0$ | $< 10^{-4}$ relative |
| Metric symmetry | $g_{\mu\nu} = g_{\nu\mu}$ | Exact |
| Christoffel symmetry | $\Gamma^\lambda_{\mu\nu} = \Gamma^\lambda_{\nu\mu}$ | Exact |

### 4.2 Physical Reference Values

All implementations must reproduce known results:

| Test Case | Expected Value | Tolerance | Reference |
|-----------|----------------|-----------|-----------|
| Schwarzschild photon sphere | $r = 3M$ | Exact | MTW §25.3 |
| Schwarzschild ISCO | $r = 6M$ | Exact | MTW §25.5 |
| Kerr ISCO (prograde, $a=M$) | $r = M$ | Exact | Bardeen 1972 |
| Weak field deflection | $\Delta\phi = 4M/b$ | $< 1\%$ | MTW §40.3 |

---

## 5. Precision Requirements

### 5.1 Floating-Point Characteristics

| Property | Double | Single (GPU) |
|----------|--------|--------------|
| Significand bits | 52 | 23 |
| Decimal digits | ~15 | ~7 |
| Machine epsilon | $2.2 \times 10^{-16}$ | $1.2 \times 10^{-7}$ |

### 5.2 Error Accumulation Model

For $n$ sequential operations with relative error $\epsilon$:

$$
\text{Accumulated error} \approx n \cdot \epsilon \cdot |x|
$$

| Path | Operations | GPU Error Bound | CPU Error Bound |
|------|------------|-----------------|-----------------|
| RK4 step | ~100 | $\approx 10^{-5}$ | $\approx 10^{-14}$ |
| Full ray (500 steps) | ~50,000 | $\approx 10^{-2}$ | $\approx 10^{-11}$ |

### 5.3 Tolerance Selection

| Context | Tolerance |
|---------|-----------|
| Null constraint (GPU) | $10^{-5}$ |
| Null constraint (CPU) | $10^{-10}$ |
| Conservation laws | $10^{-4}$ relative |
| Coordinate comparison | $10^{-12}$ |

---

## 6. Determinism Requirements

### 6.1 Requirement

Identical inputs MUST produce identical outputs across runs, platforms, and builds. This is necessary for regression testing, debugging by bisection, and scientific reproducibility.

### 6.2 Constraints

| Requirement | Specification |
|-------------|---------------|
| Floating-point consistency | Use consistent FMA settings |
| No uninitialised reads | All variables initialised before use |
| Fixed RNG seeds | Randomness must be seedable |
| Ordered reductions | Reduction order must be deterministic |
| Thread-safe access | No race conditions in parallel code |

---

## 7. Domain Constraints

### 7.1 Coordinate Bounds

| Coordinate | Constraint | Enforcement |
|------------|------------|-------------|
| $r$ (Schwarzschild) | $r > 2M \cdot 1.001$ | Clamp at boundary |
| $r$ (Kerr) | $r > r_+ + \epsilon$ | Clamp at boundary |
| $\theta$ | $(\epsilon, \pi - \epsilon)$ | Clamp at poles |
| $\phi$ | $[0, 2\pi)$ | Wrap with modulo |

### 7.2 Parameter Ranges

| Parameter | Range | Default |
|-----------|-------|---------|
| Black hole mass $M$ | $[0.1, 100]$ | 1.0 |
| Spin parameter $a/M$ | $[0, 0.998]$ | 0.0 |
| Charge $Q/M$ | $[0, 0.999]$ | 0.0 |
| Observer distance | $[5M, 1000M]$ | 20M |
| Integration steps | $[100, 10000]$ | 1500 |

---

## 8. Test Coverage Requirements

### 8.1 Test-to-Calculation Ratio

Every mathematical calculation MUST have a corresponding test. Target ratio: 1.0.

### 8.2 Test Categories

| Category | Code Prefix | Gating |
|----------|-------------|--------|
| Unit tests | `TSXX` | Mandatory |
| Integration tests | `TSIN` | Mandatory |
| Diagnostic tests | `TSDG` | Mandatory |
| Benchmark tests | `TSBM` | Soft |

### 8.3 Mandatory Tests

The following tests MUST pass for any valid build:

| Test | Validates |
|------|-----------|
| `TSPH001A` | Schwarzschild metric properties |
| `TSPH002A` | Christoffel symbol symmetry |
| `TSDG001A` | Numerical stability |
| `TSDG002A` | NaN/Inf detection |
| `TSDG003A` | Conservation laws |
| `TSDG004A` | Determinism |

---

## 9. Error Handling

### 9.1 Error Categories

| Category | Severity | Response |
|----------|----------|----------|
| Domain error | Warning | Clamp and continue |
| Numerical failure | Error | Terminate ray |
| Configuration error | Critical | Fail startup |

### 9.2 Failure Modes

| Operation | Failure Condition | Detection | Response |
|-----------|-------------------|-----------|----------|
| Division | Divisor = 0 | Check before | Return NaN, log |
| Square root | Argument < 0 | Check before | Return NaN, log |
| Metric evaluation | Outside domain | Coordinate check | Clamp |
| Integration step | NaN in state | `isnan()` check | Terminate ray |

### 9.3 NaN Policy

- NaN MUST NOT propagate silently
- All NaN occurrences MUST be logged with context
- Affected rays MUST be terminated with error status

---

## 10. Configuration

### 10.1 Required Parameters

| Parameter | Type | Default | Range |
|-----------|------|---------|-------|
| `resolution.width` | int | 1920 | [128, 8192] |
| `resolution.height` | int | 1080 | [128, 8192] |
| `integration.maxSteps` | int | 1500 | [100, 10000] |
| `integration.stepSize` | float | 0.1 | [0.001, 1.0] |
| `integration.tolerance` | float | $10^{-6}$ | [$10^{-12}$, $10^{-3}$] |
| `metric.type` | enum | Schwarzschild | - |
| `metric.mass` | float | 1.0 | [0.1, 100] |

### 10.2 Configuration Validation

All configuration MUST be validated at startup. Invalid configuration MUST prevent system start.

---

## 11. Logging and Audit

### 11.1 Required Log Events

| Event | Log Level | Required Fields |
|-------|-----------|-----------------|
| Frame rendered | INFO | Frame number, time, ray count |
| Metric loaded | INFO | Metric type, parameters |
| Ray terminated (error) | WARNING | Position, reason |
| Conservation violated | WARNING | Quantity, drift magnitude |
| Integration failed | ERROR | Position, NaN location |

---

## Appendix A: Error Codes

| Code | Name | Description |
|------|------|-------------|
| E001 | OutOfDomain | Position outside valid region |
| E002 | SingularMetric | Metric determinant near zero |
| E003 | NaNProduced | NaN detected in computation |
| E004 | StepLimitExceeded | Maximum integration steps reached |
| E005 | ConvergenceFailed | Iterative method did not converge |
| E006 | InvalidConfiguration | Configuration parameter invalid |

---

## Appendix B: Glossary

| Term | Definition |
|------|------------|
| Affine parameter | Path parameter along a geodesic, preserving the linear relationship between tangent vector magnitude and proper time (or distance) |
| Christoffel symbols | Connection coefficients $\Gamma^\lambda_{\mu\nu}$, derived from the metric tensor and encoding spacetime curvature |
| Geodesic | Path through spacetime that extremises proper length (timelike) or satisfies parallel transport (null) |
| ISCO | Innermost stable circular orbit; the smallest circular orbit stable against radial perturbations |
| Null | Zero spacetime interval ($ds^2 = 0$), the condition satisfied by light rays |
| Photon sphere | Surface of unstable circular photon orbits ($r = 3M$ for Schwarzschild) |

