# Sirius Refactor Specification

*Disciplined Change in High-Integrity Code*

> "We are what we repeatedly do. Excellence, then, is not an act, but a habit."
> — Aristotle

## Preface

Refactoring in Sirius is a disciplined act of stewardship. The system computes geodesics in curved spacetime, renders gravitational lensing, and produces scientific visualisations. This demands clarity of intent, precision of meaning, and respect for numerical stability.

This document provides a repeatable method for changing code whilst preserving integrity, performance, and auditability. It complements the [Coding Standard](standard.md) and [Type System](types.md).

---

## Part I: First Principles

### 1.1 Purpose and Consequence

A refactor is justified when it improves correctness, performance, maintainability, or auditability without altering intended behaviour. Each change must trace to a concrete purpose. Refactoring is not aesthetic decoration; it is a method of preserving long-term system integrity.

### 1.2 Meaning and Accountability

Every function is a statement about physics. Every type is a commitment to meaning. A refactor must preserve those commitments or replace them with stronger ones. The developer inherits accountability for any change that weakens meaning.

### 1.3 Boundaries and Invariants

Invariants are enforced at boundaries, not in every line of interior logic. A refactor must preserve this boundary discipline. Validation is explicit and local; internal logic assumes validated inputs and remains efficient.

### 1.4 Determinism and Auditability

Scientific software requires deterministic behaviour. Refactors must reduce ambiguity, remove hidden costs, and increase traceability. A refactor that obscures the path of data is a regression, even if it shortens a file.

---

## Part II: Refactor Contract

### 2.1 Scope Statement

Define the target components, files, and boundaries. Record the allowed change type.

| Change Type | Description |
|-------------|-------------|
| Mechanical | Code shape changes only; behaviour preserved |
| Structural | Reorganisation within the same behaviour and contract |
| Semantic | Clarifies meaning with stronger types or explicit invariants |
| Behavioural | Changes output or timing; requires explicit approval |

### 2.2 Behaviour Baseline

Document the current contract in practical terms:

- Inputs, outputs, and failure modes
- Expected latency in hot paths, where known
- External dependencies and data sources
- Mathematical invariants preserved

### 2.3 Performance Budget

Define the performance envelope for hot paths. A refactor SHOULD keep the envelope stable or improve it. Any measured regression MUST be surfaced as a risk.

### 2.4 Evidence Requirements

Each refactor SHOULD provide evidence proportional to risk:

| Risk Level | Evidence Required |
|------------|-------------------|
| Low | Unit tests passing |
| Medium | Unit tests, invariant checks |
| High | Full test suite, benchmarks, review |

---

## Part III: Workflow

### 3.1 Orientation

Identify the location of truth. Trace from data entry to output. Summarise the current flow in plain language. This step prevents refactoring the wrong surface.

**Questions to answer:**

- Where does the data enter the system?
- What transformations occur?
- What invariants are assumed?
- What outputs are produced?

### 3.2 Boundary Audit

List every boundary and its invariants:

| Boundary | Invariant | Enforcement |
|----------|-----------|-------------|
| Metric entry | $r > r_{horizon}$ | Coordinate clamp |
| Integration step | Null constraint | Renormalisation |
| Ray output | Valid pixel | NaN check |

Confirm each validation is explicit, local, and complete. Confirm internal logic assumes validated inputs.

### 3.3 Hot Path Mapping

Identify:

- Loops (geodesic integration, per-pixel computation)
- Pricing kernels (Christoffel evaluation, metric computation)
- Frequently called functions (vector operations, coordinate transforms)

Mark any allocations, virtual calls, and external I/O. Use this map to guide mechanical sympathy work.

### 3.4 Refactor Design

Plan the transformation in a minimal sequence:

1. Each step SHOULD be reversible
2. Each step SHOULD be independently testable
3. Avoid broad changes that mingle behaviour with structure

### 3.5 Implementation

Apply changes in small passes:

| Pass | Focus |
|------|-------|
| First | Clarify intent: explicit names, explicit types |
| Second | Localise logic: keep computation near its data |
| Third | Performance alignment: resource control, allocation reduction |

### 3.6 Verification

Run tests or provide equivalent evidence. If tooling cannot run in the current environment, record the constraint and provide a precise command for the operator.

```bash
cd Sirius.Build
ctest -C Release --output-on-failure
```

### 3.7 Review and Audit Trail

Record:

- Changes made
- Risks identified
- Follow-up actions required
- Audit note for future investigation

---

## Part IV: Refactor Principles

### 4.1 Mechanical Sympathy

Principles:

- Prefer contiguous data structures and span-based access
- Remove hidden allocations in hot paths
- Reuse buffers and keep GC pressure predictable

**Example:**

```cpp
// Before: repeated allocations inside loop
for (int i = 0; i < n; i++) {
    auto buffer = new double[m];
    Compute(buffer, i);
    delete[] buffer;
}

// After: reuse buffer with explicit lifetime
double* buffer = new double[m];
for (int i = 0; i < n; i++) {
    Compute(buffer, i);
}
delete[] buffer;
```

### 4.2 Logic Proximity

Principles:

- Keep the calculation close to its data
- Favour linear flow in complex mathematics and state transitions
- Defer abstraction until a second concrete use exists

**Example:**

```cpp
// Before: fragmented across multiple helpers
auto a = ComputeA(x);
auto b = ComputeB(a);
return ComputeC(b, y);

// After: localised computation
double a = x * x + k;
double b = a / (y + epsilon);
return b * factor;
```

### 4.3 Explicit Intent

Principles:

- Use named locals for predicates and thresholds
- Replace sentinel values with types or explicit results
- Reserve comments for external references and exceptions

**Example:**

```cpp
bool isConverged = error < tolerance;
if (isConverged) {
    return result;
}
```

### 4.4 Failure First

Principles:

- Validate inputs at the boundary
- Use deterministic fallbacks with explicit limits
- Prefer explicit error types or exceptions to silent nulls

**Example:**

```cpp
if (steps >= maxSteps) {
    return IntegrationResult::MaxStepsExceeded;
}
```

### 4.5 Deterministic Resources

Principles:

- Dispose created resources in the same scope
- Avoid ambiguous ownership
- Use clear lifetime boundaries for pooled or shared objects

---

## Part V: GPU-Specific Considerations

### 5.1 Shader Refactoring

GPU code requires additional care:

- No dynamic allocation; fixed-size arrays only
- Explicit precision (`float` vs `double`)
- Avoid divergent branches in tight loops
- Memory access patterns affect performance significantly

### 5.2 Host-Device Boundary

Changes to data structures crossing the host-device boundary require:

- Verification that layout matches
- Testing on multiple GPU architectures
- Performance measurement before and after

---

## Part VI: Evidence and Testing

### 6.1 Unit Tests

Every refactor that touches logic requires tests for:

- Invariants at boundaries
- Failure cases
- Stable outputs for known inputs

### 6.2 Regression Safety

When behaviour must remain stable, capture it with:

- Fixed input fixtures
- Deterministic random seeds
- Numeric tolerances that reflect model sensitivity

### 6.3 Performance Evidence

For hot paths, add or update benchmarks and compare against the previous envelope. Record any expected deviation with a rationale.

---

## Part VII: Documentation and Knowledge

### 7.1 Alignment with Standards

Confirm that the refactor respects:

- [Coding Standard](standard.md)
- [Type System](types.md)
- Component-specific requirements in this documentation

### 7.2 Philosophical Integrity

Before finalising, ask: does this change make the system more truthful?

Truth here means correspondence between code, meaning, and physical execution. A refactor that improves this correspondence is valuable; one that obscures it is harmful regardless of other benefits.

---

## Appendix: Refactor Checklist

- [ ] Purpose and scope defined
- [ ] Boundary invariants explicit
- [ ] Hot paths mapped and protected
- [ ] Logic close to data
- [ ] Explicit naming and types
- [ ] Resource lifetimes deterministic
- [ ] Evidence recorded
- [ ] Documentation aligned

---

*End of Refactor Specification*
