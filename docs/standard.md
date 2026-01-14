# Sirius Coding Standard

**Based on:** JPL Institutional Coding Standard, MISRA C++:2023, GPU Programming Best Practices  
**Applicability:** High-Performance Relativistic Ray Tracing

> "Programs must be written for people to read, and only incidentally for machines to execute."
> — Harold Abelson

---

## Preface

This standard establishes coding conventions for the Sirius codebase. The system involves numerical integration of geodesic equations, GPU compute shaders, and tensor calculus; these domains demand predictable execution, numerical stability, and hardware compatibility.

The standard uses three requirement levels:

| Level | Meaning |
|-------|---------|
| MUST | Mandatory; violations are defects |
| SHOULD | Recommended; deviations require justification |
| MAY | Optional; permitted when beneficial |

---

## Part I: Language Compliance

### 1.1 Language Version

All code MUST conform to ISO C++17. Use of compiler-specific extensions requires documented justification.

### 1.2 Compiler Warnings

All code MUST compile with maximum warnings enabled:

```cmake
target_compile_options(${TARGET} PRIVATE -Wall -Wextra -Wpedantic)
```

Third-party code in `Sirius.Dependencies/` is excluded from this requirement.

### 1.3 No Undefined Behaviour

The following are prohibited:

- Reading uninitialised variables
- Signed integer overflow
- Strict aliasing violations
- Use-after-free, double-free
- Null pointer dereference

---

## Part II: Predictable Execution

### 2.1 Bounded Loops

All loops MUST have a determinable upper bound or be bounded by a maximum iteration count.

**Allowed:**

```cpp
for (int i = 0; i < 4; ++i) { ... }                      // Fixed bound
for (int i = 0; i < maxSteps && !converged; ++i) { ... } // Bounded iteration
```

**Prohibited:**

```cpp
while (true) { ... }  // Unbounded
```

### 2.2 No Recursion

Direct and indirect recursion are prohibited. Stack overflow in numerical code corrupts results silently. Use explicit iteration with stack data structures when necessary.

### 2.3 Stack Allocation in Hot Paths

Numerical kernels MUST NOT allocate heap memory. Use:

- Stack-allocated arrays (`std::array<T, N>`)
- Pre-allocated buffers passed as parameters
- Object pools for complex types

```cpp
// Preferred: Stack allocation
std::array<double, 4> christoffel;

// Prohibited in hot paths: Heap allocation
std::vector<double> christoffel(4);
```

### 2.4 Const Correctness

Use `const` for:

- Parameters that are not modified
- Member functions that do not mutate state
- Local variables that are initialised once

---

## Part III: Defensive Coding

### 3.1 Null Safety

Use references, `std::optional`, or smart pointers instead of raw pointers.

Exception: C API interop and performance-critical GPU code may use raw pointers with documented lifetime.

### 3.2 Limited Scope

Variables MUST be declared at the smallest possible scope:

```cpp
// Preferred
for (int i = 0; i < n; ++i) {
    double value = compute(i);  // Declared where used
    process(value);
}

// Avoid
double value;  // Declared far from use
for (int i = 0; i < n; ++i) {
    value = compute(i);
    process(value);
}
```

### 3.3 Guard Clauses

Public functions MUST validate parameters at entry:

```cpp
void integrate(const Metric& g, double dt) {
    if (dt <= 0.0) {
        throw std::invalid_argument("dt must be positive");
    }
    // Proceed with valid parameters
}
```

### 3.4 Error Handling in Hot Paths

Use return codes or `std::optional` instead of exceptions in numerical code:

```cpp
// Preferred for hot paths
std::optional<double> computeChristoffel(const Metric& g, int mu, int alpha, int beta);

// Exceptions allowed for initialisation/setup
void loadShader(const std::string& path);  // May throw
```

---

## Part IV: Code Clarity

### 4.1 Limited Preprocessor

Preprocessor macros MUST be limited to:

- Include guards
- Platform detection (`#ifdef _WIN32`)
- Debug assertions

**Prohibited:**

```cpp
#define SQUARE(x) ((x) * (x))  // Use inline function
```

### 4.2 Function Size

Functions SHOULD NOT exceed 60 lines. Extract helper functions for complex logic.

### 4.3 Naming Conventions

| Element | Convention | Example |
|---------|------------|---------|
| Component ID | `XXYY###V` | `PHMT001A` |
| Class | `PascalCase` | `MetricTensor` |
| Function/method | `camelCase` | `evaluateChristoffel` |
| Constant | `SCREAMING_SNAKE_CASE` | `MAX_STEPS` |
| Variable | `camelCase` | `stepSize` |
| Template parameter | `PascalCase` | `Precision` |

### 4.4 Explicit Intent

Use named locals for predicates and thresholds:

```cpp
bool isConverged = error < tolerance;
if (isConverged) {
    return result;
}
```

---

## Part V: GPU Considerations

### 5.1 GLSL Compatibility

Shader code MUST:

- Use explicit `layout` qualifiers for all inputs/outputs
- Avoid vendor-specific extensions without fallbacks
- Be tested on AMD, NVIDIA, and Intel drivers

### 5.2 No Dynamic Allocation in Shaders

GLSL code MUST use fixed-size arrays:

```glsl
// Allowed
float christoffel[64];  // 4x4x4

// Dynamic allocation is not possible in GLSL
```

### 5.3 Precision Documentation

Document precision requirements for numerical operations:

```cpp
// Uses double precision for numerical stability
// GPU shader uses float with compensated summation
double integrateGeodesic(const Vec4& position, const Vec4& velocity);
```

---

## Part VI: Testing Requirements

### 6.1 Test Categories

| Category | Purpose | Gating |
|----------|---------|--------|
| Unit tests | Component logic and invariants | Mandatory |
| Boundary tests | Validation and failure paths | Mandatory |
| Diagnostic tests | Numerical stability | Mandatory |
| Benchmark tests | Performance regression | Soft gate |

### 6.2 Test Coverage

Every mathematical calculation MUST have a corresponding test. Target test-to-calculation ratio: 1.0.

### 6.3 Test Independence

Tests MUST be independent: no shared mutable state, deterministic execution, reproducible results.

---

## Part VII: Documentation Standards

### 7.1 Header Documentation

Every public function MUST have documentation:

```cpp
/// Compute Christoffel symbols at a given position.
/// @param position Spacetime coordinates (MUST be in valid domain)
/// @param[out] christoffel Output array (4x4x4)
/// @return true if successful, false if position is outside domain
bool computeChristoffel(const Vec4& position, Christoffel4D& christoffel);
```

### 7.2 Mathematical Notation

Document the mathematical basis for calculations:

```cpp
/// Compute the geodesic acceleration.
/// Uses the geodesic equation:
///   d²x^μ/dλ² = -Γ^μ_αβ dx^α/dλ dx^β/dλ
```

### 7.3 Reference Citations

Cite sources for non-trivial formulas:

```cpp
/// Novikov-Thorne temperature profile.
/// Reference: Page & Thorne (1974), ApJ 191, 499
```

---

## Part VIII: Performance

### 8.1 Allocation Awareness

Hot paths MUST avoid allocations. Profile before optimising; measure after.

### 8.2 Cache Efficiency

Prefer contiguous data structures. Access memory linearly where possible. Mark performance-critical code for review.

### 8.3 Benchmarking

Performance-critical code MUST have benchmarks. Regressions exceeding 10% require investigation.

---

## Appendix A: Rule Summary

| Rule | Level | Statement |
|------|-------|-----------|
| 1 | MUST | Conform to C++17 |
| 2 | MUST | Compile with maximum warnings |
| 3 | MUST | All loops have bounds |
| 4 | MUST | No recursion |
| 5 | MUST | No heap allocation in hot paths |
| 6 | SHOULD | Use const where applicable |
| 7 | SHOULD | Prefer references over pointers |
| 8 | MUST | Validate parameters at entry |
| 9 | SHOULD | Use std::optional for optional results |
| 10 | MUST | No undefined behaviour |
| 11 | SHOULD | Functions under 60 lines |
| 12 | MUST | Follow naming conventions |
| 13 | MUST | Test every calculation |

---

## References

1. **JPL D-60411**: JPL Institutional Coding Standard for C
2. **MISRA C++:2023**: Guidelines for the use of C++17 in critical systems
3. **OpenGL 4.6 Specification**: Khronos Group
4. **The Power of 10**: Rules for Developing Safety-Critical Code, Gerard J. Holzmann
5. **C++ Core Guidelines**: Bjarne Stroustrup and Herb Sutter

---

*End of Coding Standard*
