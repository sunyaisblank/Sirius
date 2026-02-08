# Refactoring

A repeatable method for changing code in Sirius while preserving correctness, performance, and auditability. The domain (numerical integration, tensor calculus, GPU kernels) makes careless refactoring expensive: a change that introduces a subtle numerical error may not produce a visible artefact until the ray passes close to a photon sphere, and by then the bug is far from the diff that introduced it.

This document complements the [Coding Standard](standard.md) and [Type System](type.md).

---

## Part I: Ground Rules

### 1.1 Justification

A refactor is justified when it improves correctness, performance, maintainability, or auditability without altering intended behaviour. Each change traces to a concrete purpose. "The code looks nicer" is not a purpose; "the boundary validation is now explicit rather than implicit" is.

### 1.2 Preserving Contracts

A refactor must preserve the contracts of the code it touches, or replace them with strictly stronger ones. If a function guaranteed $r > 2M$ at exit before the refactor, it must guarantee at least that afterward. If a type encoded a coordinate system assumption, the refactored type must encode the same assumption or make it more explicit.

### 1.3 Boundary Discipline

Sirius validates invariants at boundaries (module entry points, public function parameters) and assumes validity in internal logic. A refactor must preserve this pattern. Moving a validation check from a boundary into interior logic degrades performance; removing it from the boundary removes a guarantee that downstream code relies on.

### 1.4 Determinism

A refactor that changes the order of floating-point operations may change the output. In scientific software that requires bit-exact reproducibility, this is a behavioural change, not a mechanical one. Any refactor that reorders arithmetic in a hot path must be classified as behavioural and tested accordingly.

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

```mermaid
flowchart LR
    A["Orientation"] --> B["Boundary\nAudit"]
    B --> C["Hot Path\nMapping"]
    C --> D["Refactor\nDesign"]
    D --> E["Implementation"]
    E --> F["Verification"]
    F --> G["Review"]
    F -->|tests fail| E
```

### 3.1 Orientation

Trace the data flow from entry to output and summarise it in plain language before touching any code. This step prevents refactoring the wrong surface, changing the shape of code that was not the problem.

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

Align data layout and access patterns with the hardware:

- Prefer contiguous data structures and span-based access (cache lines reward locality)
- Remove hidden allocations in hot paths (heap allocation serialises threads)
- Reuse buffers rather than allocating and freeing per iteration

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

Keep computation close to its data:

- Favour linear control flow in mathematical code so the reader can follow the derivation top to bottom
- Defer abstraction until a second concrete use exists (premature extraction scatters logic across files without reducing complexity)

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

Make the code say what it means:

- Use named locals for predicates and thresholds (a reader should not have to decode `if (r < 2.002)`; call it `isNearHorizon`)
- Replace sentinel values with types or explicit results (`std::optional` over magic return values)
- Reserve comments for external references and non-obvious constraints; do not narrate the code

**Example:**

```cpp
bool isConverged = error < tolerance;
if (isConverged) {
    return result;
}
```

### 4.4 Failure First

Handle failure before the happy path:

- Validate inputs at the boundary and reject early
- Use deterministic fallbacks with explicit limits (e.g., max iteration counts)
- Prefer explicit error types over silent null returns

**Example:**

```cpp
if (steps >= maxSteps) {
    return IntegrationResult::MaxStepsExceeded;
}
```

### 4.5 Deterministic Resources

Manage lifetimes explicitly:

- Dispose created resources in the same scope that created them
- Avoid ambiguous ownership (if two objects hold a pointer, document which one frees it)
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
- [Type System](type.md)
- Component-specific requirements in this documentation

### 7.2 Consistency Check

Before finalising, verify that the refactored code matches its documentation: variable names correspond to documented quantities, types match documented invariants, and the data flow described in the architecture document still holds. A refactor that makes the code better but makes the documentation wrong has traded one problem for another.

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

