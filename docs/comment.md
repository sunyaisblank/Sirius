# Commentary Standard

---

## 1. Scope and Purpose

### 1.1 Scope

This specification defines requirements for inline documentation and commentary within source code files. It applies to all textual annotations embedded within source code that are intended for human readers and are syntactically ignored by compilers, interpreters, or other language processors.

This specification does not govern:
- External documentation (wikis, manuals, README files)
- Automatically generated documentation output
- Commit messages or version control annotations
- Issue tracker content
- Code formatting and organisation (e.g., include ordering, whitespace)

**Note:** Commentary standards operate alongside code formatting standards. Projects typically define formatting conventions (e.g., Google C++ Style Guide, Linux kernel style) and enforce them via automated formatters. This specification assumes such formatting standards exist and focuses exclusively on the content and purpose of comments.

### 1.2 Purpose

This specification exists to:
1. Establish measurable compliance criteria for code commentary
2. Ensure commentary serves its intended function of aiding comprehension
3. Prevent commentary that degrades maintainability
4. Provide a framework adaptable across programming languages and abstraction levels

### 1.3 Conformance Language

The key words "MUST", "MUST NOT", "REQUIRED", "SHALL", "SHALL NOT", "SHOULD", "SHOULD NOT", "RECOMMENDED", "MAY", and "OPTIONAL" in this document are to be interpreted as described in RFC 2119.

### 1.4 Relationship to External Style Guides

This specification is designed to complement, not replace, language-specific or organisation-specific style guides. Where this specification and an external style guide conflict, the external style guide takes precedence for language-specific syntax matters, while this specification takes precedence for semantic content requirements.

---

## 2. Definitions

### 2.1 Commentary Classifications

**Comment**: Any sequence of characters within a source file that is syntactically designated as non-executable annotation, regardless of the specific delimiter syntax used by the language.

**Documentation Comment**: A comment that uses a structured format recognised by documentation generation tools (e.g., Doxygen, JSDoc, docstrings, rustdoc).

**Inline Comment**: A comment appearing on the same line as executable code, typically following the code.

**Block Comment**: A comment occupying one or more complete lines, with no executable code on those lines.

**Section Comment**: A block comment that delineates a logical grouping of code elements within a file.

**Phase Label**: A brief inline marker identifying a logical stage within a function or code block.

**Tombstone Comment**: A comment marking code that has been removed or disabled, or explaining why something is intentionally absent.

**Instrumentation Comment**: A comment or annotation that serves both documentation and runtime tooling purposes (e.g., profiling scope markers).

### 2.2 Abstraction Level Classifications

For the purposes of this specification, programming contexts are classified into three abstraction levels:

**Level 0 (Hardware-Proximate)**: Code that directly manipulates hardware registers, memory addresses, interrupt vectors, or processor-specific instructions. Examples: kernel modules, device drivers, embedded firmware, assembly language, MMIO operations.

**Level 1 (Systems)**: Code that manages system resources, implements data structures, or provides foundational services without direct hardware manipulation. Examples: operating system services, runtime libraries, memory allocators, network protocol implementations, compilers, code generators.

**Level 2 (Application)**: Code that implements business logic, user interfaces, or high-level algorithms using abstractions provided by lower levels. Examples: web applications, data processing pipelines, configuration scripts, domain-specific tools.

### 2.3 Visibility Classifications

**Public Interface**: Any symbol, function, type, or module that is intended for use by code outside its defining compilation unit, package, or module boundary.

**Internal Interface**: Any symbol intended for use within a defined subsystem boundary but not by external consumers.

**Private Implementation**: Any symbol whose use is confined to its immediate defining context (e.g., file-static functions, private class members, closure-local bindings).

### 2.4 Domain Classifications

**General-Purpose Code**: Code that uses standard programming constructs without specialised domain knowledge.

**Domain-Specific Code**: Code that implements algorithms or concepts from a specialised field (e.g., quantum computing, signal processing, financial mathematics, computational physics, machine learning) where domain expertise is required to understand the implementation.

---

## 3. File-Level Requirements

### 3.1 License and Copyright Header

Every source file MUST begin with a license and copyright header unless the project's build system or packaging mechanism provides this information through an alternative verified mechanism.

**Requirement F-LIC-001**: The header MUST appear before any executable code or import statements.

**Requirement F-LIC-002**: The header SHOULD use SPDX identifiers where applicable.

**Requirement F-LIC-003**: The header MUST NOT exceed 20 lines excluding blank lines within the header block.

**Compliant Example**:
```
// SPDX-FileCopyrightText: Copyright (c) 2024 Organisation Name
// SPDX-License-Identifier: MIT
```

**Compliant Example (Extended)**:
```
/*
 * Copyright (c) 2024 Organisation Name
 * 
 * This file is part of ProjectName.
 * Licensed under the Apache License, Version 2.0.
 * See LICENSE file for details.
 */
```

### 3.2 File Purpose Statement

**Requirement F-PUR-001**: Source files exceeding 100 lines SHOULD contain a brief statement of the file's purpose, placed immediately after the license header.

**Requirement F-PUR-002**: The purpose statement MUST NOT duplicate information readily apparent from the filename and containing namespace or module path.

**Requirement F-PUR-003**: The purpose statement SHOULD identify the primary abstraction, subsystem, or responsibility implemented by the file.

**Non-Compliant Example**:
```
// user.cpp - Implementation of user functions
```
This merely restates the filename.

**Compliant Example**:
```
// Implements credential validation and session token generation for the
// authentication subsystem. Depends on the cryptographic primitives defined
// in crypto/hash.h.
```

**Compliant by Omission**: A file named `python_bindings.cpp` within namespace `nvfuser::python_frontend` needs no purpose statement if its role is evident from this context.

### 3.3 File Structure Comments

**Requirement F-STR-001**: Files exceeding 500 lines SHOULD use section comments to delineate logical groupings.

**Requirement F-STR-002**: Section comments MUST use a consistent visual format throughout the file.

**Requirement F-STR-003**: Section comments MUST NOT be used to compensate for poor file organisation; files requiring more than five section comments SHOULD be evaluated for decomposition.

**Compliant Example (Equals Style)**:
```
// ============================================================================
// Public Interface
// ============================================================================

// ============================================================================
// Internal Helpers
// ============================================================================
```

**Compliant Example (Box Style)**:
```
/******************************************************************************
 * Public Interface
 ******************************************************************************/
```

**Compliant Example (Minimal Style)**:
```
// --- Public Interface ---

// --- Internal Helpers ---
```

---

## 4. Declaration-Level Requirements

### 4.1 Public Interface Documentation

**Requirement D-PUB-001**: Every public interface element MUST have a documentation comment.

**Requirement D-PUB-002**: Documentation comments for public interfaces MUST include:
- A single-sentence summary (the "brief")
- Parameter descriptions for each parameter whose purpose is not self-evident from its name and type
- Return value description if the return type is not `void`, `unit`, or equivalent
- Exception or error condition descriptions if applicable

**Requirement D-PUB-003**: Documentation comments MUST be placed immediately preceding the declaration, with no intervening blank lines or code.

**Requirement D-PUB-004**: Documentation comments MUST use the documentation comment syntax appropriate to the language's ecosystem tooling (e.g., `///` or `/** */` for Doxygen, `"""` for Python docstrings, `///` for rustdoc).

**Compliant Example (C++)**:
```cpp
/// Computes the SHA-256 digest of the input buffer.
///
/// @param data     Pointer to the input byte sequence.
/// @param length   Number of bytes to process.
/// @param out      Output buffer; must have capacity for 32 bytes.
///
/// @return True if computation succeeded; false if any pointer is null.
///
/// @note This function is not thread-safe with respect to the global
///       entropy pool.
bool sha256_digest(const uint8_t* data, size_t length, uint8_t* out);
```

**Compliant Example (Python)**:
```python
def sha256_digest(data: bytes) -> bytes:
    """Compute the SHA-256 digest of the input bytes.

    Args:
        data: The byte sequence to hash.

    Returns:
        A 32-byte digest.

    Raises:
        ValueError: If data is empty.
    """
```

### 4.2 Internal Interface Documentation

**Requirement D-INT-001**: Internal interface elements SHOULD have documentation comments.

**Requirement D-INT-002**: Internal documentation MAY omit parameter descriptions where names and types are sufficient.

**Requirement D-INT-003**: Internal documentation MUST document non-obvious preconditions, postconditions, or invariants.

**Compliant Example**:
```cpp
// Set of local functions that are used to compose python FusionDefinition
// bindings. Ideally, these would be templated lambda functions but those
// are not available without C++20.
namespace {
Vector define_vector_base_fn(FusionDefinition& fd, std::vector<Scalar>& args) {
    ...
}
}
```

### 4.3 Private Implementation Documentation

**Requirement D-PRI-001**: Private implementation elements MAY omit documentation comments if their purpose is evident from context.

**Requirement D-PRI-002**: Private implementation elements MUST have documentation comments if:
- The implementation uses a non-obvious algorithm
- The implementation has subtle correctness requirements
- The implementation interacts with external state in non-apparent ways

### 4.4 Type and Structure Documentation

**Requirement D-TYP-001**: Type definitions (classes, structs, enums, interfaces, traits) with public visibility MUST include a documentation comment describing the type's role and responsibilities.

**Requirement D-TYP-002**: Enumeration values SHOULD have individual comments if the enumeration represents a non-trivial domain concept.

**Requirement D-TYP-003**: Struct or class fields that represent invariants or have validity constraints MUST document those constraints.

**Compliant Example**:
```cpp
/// Represents a half-open interval [start, end) in a one-dimensional space.
///
/// Invariant: start <= end
struct Interval {
    /// Inclusive lower bound.
    int start;
    
    /// Exclusive upper bound. Must be >= start.
    int end;
};
```

### 4.5 Constant and Configuration Documentation

**Requirement D-CON-001**: Named constants representing domain-specific values MUST include a comment explaining the value's origin or derivation.

**Requirement D-CON-002**: Magic numbers (numeric literals with non-obvious meaning) MUST NOT appear in code without either:
- Extraction to a named constant with appropriate documentation, OR
- An inline comment explaining the value

**Non-Compliant Example**:
```cpp
if (retry_count > 3) { ... }
```

**Compliant Example**:
```cpp
// Maximum retries before circuit breaker opens, per SLA agreement §4.2
constexpr int kMaxRetries = 3;

if (retry_count > kMaxRetries) { ... }
```

**Compliant Example (Inline Alternative)**:
```cpp
if (retry_count > 3) { ... }  // SLA mandates max 3 retry attempts
```

---

## 5. Implementation-Level Requirements

### 5.1 The "Why" Principle

**Requirement I-WHY-001**: Comments within implementation bodies MUST explain *why* code exists or *why* an approach was chosen, not *what* the code does.

**Requirement I-WHY-002**: Comments that restate what code does in natural language MUST NOT be present unless the code itself uses domain terminology that requires translation for the reader.

**Non-Compliant Example**:
```cpp
i++;  // increment i
```

**Non-Compliant Example**:
```cpp
// Loop through the array
for (int i = 0; i < n; i++) {
    // Check if element is greater than threshold
    if (arr[i] > threshold) {
        // Add to result
        result.push_back(arr[i]);
    }
}
```

**Compliant Example**:
```cpp
// Filter to elements exceeding threshold; order preservation is required
// by downstream correlation analysis.
for (int i = 0; i < n; i++) {
    if (arr[i] > threshold) {
        result.push_back(arr[i]);
    }
}
```

### 5.2 Domain-Specific Notation

**Requirement I-DOM-001**: Code operating in specialised domains SHOULD include comments that:
- Define domain-specific notation on first use within a file or function
- Map implementation variables to their mathematical or domain counterparts
- Reference authoritative sources for algorithms (textbooks, papers, specifications)

**Requirement I-DOM-002**: Mathematical expressions or domain formulas implemented in code SHOULD include a comment showing the formula in standard notation.

**Compliant Example (Quantum Computing)**:
```cpp
// Compute (U) Action (V) produces U V U†
// where † denotes the adjoint (conjugate transpose)
cudaq::compute_action(unitary_u, action_v, register);

// Perform ctrl-U^j for j = 2^i (phase estimation algorithm)
// Reference: Nielsen & Chuang, Chapter 5.2
for (int i = 0; i < n_counting_qubits; ++i) {
    for (int j = 0; j < (1 << i); ++j) {
        cudaq::control(oracle, counting_qubits[i], state_register);
    }
}
```

**Compliant Example (Signal Processing)**:
```cpp
// Hann window: w[n] = 0.5 * (1 - cos(2πn / N))
// Reduces spectral leakage compared to rectangular window
for (int n = 0; n < N; ++n) {
    window[n] = 0.5 * (1.0 - cos(2.0 * M_PI * n / N));
}
```

### 5.3 Algorithm and Complexity Comments

**Requirement I-ALG-001**: Non-trivial algorithms MUST include a comment identifying:
- The algorithm by name (if it has a standard name)
- The time complexity (using Big-O notation)
- The space complexity (if non-constant)
- A reference to source material (if applicable)

**Requirement I-ALG-002**: Optimisation techniques that sacrifice readability for performance MUST include a comment justifying the trade-off.

**Compliant Example**:
```cpp
// Quickselect (Hoare's algorithm) for O(n) expected-time k-th element.
// Worst case O(n²) but median-of-medians pivot selection makes this rare.
// Reference: CLRS 3rd ed., Chapter 9.
int quickselect(int* arr, int n, int k) {
    ...
}
```

**Compliant Example**:
```cpp
// Manual loop unrolling: benchmarks show 2.3x throughput improvement
// on the target architecture (ARM Cortex-A72) for this hot path.
// See perf/results/2024-01-15/unroll_comparison.txt
for (int i = 0; i < n; i += 4) {
    process(arr[i]);
    process(arr[i + 1]);
    process(arr[i + 2]);
    process(arr[i + 3]);
}
```

### 5.4 Performance Instrumentation

**Requirement I-PERF-001**: Performance-critical functions or code sections MAY include instrumentation markers that serve dual documentation and profiling purposes.

**Requirement I-PERF-002**: When instrumentation macros or annotations are used, their presence SHOULD indicate that the annotated code is performance-sensitive.

**Requirement I-PERF-003**: Projects using instrumentation markers SHOULD document the marker conventions in project-level documentation.

**Compliant Example**:
```cpp
Vector define_vector_base_fn(FusionDefinition& fd, std::vector<Scalar>& args) {
    FUSER_PERF_SCOPE("python_frontend::define_vector_base_fn");
    // The PERF_SCOPE macro indicates this function's performance is monitored
    // and that changes should be benchmarked.
    ...
}
```

### 5.5 Defensive Code Comments

**Requirement I-DEF-001**: Assertions, precondition checks, and invariant verifications SHOULD use assertion message strings that serve as documentation.

**Requirement I-DEF-002**: Defensive code guarding against "impossible" conditions MUST include a comment explaining why the condition is believed impossible and why the guard exists.

**Compliant Example**:
```cpp
NVF_CHECK(!fd.completed(), "Attempting to add to a completed definition!");
```

**Compliant Example**:
```cpp
// Should be unreachable: all enum cases handled above. Guard retained
// because external data deserialisation could produce invalid values.
default:
    return ErrorCode::InvalidState;
```

### 5.6 Control Flow Comments

**Requirement I-CTL-001**: Fall-through in switch statements (where the language permits and the fall-through is intentional) MUST be marked with an explicit comment.

**Requirement I-CTL-002**: Early returns, breaks, or continues that are not visually proximate to their controlling condition SHOULD have a comment explaining the exit condition.

**Requirement I-CTL-003**: Empty loop bodies or intentionally empty blocks MUST contain a comment confirming the emptiness is intentional.

**Compliant Example**:
```cpp
switch (state) {
    case State::Pending:
    case State::Queued:
        // fall through: both states require initialisation
    case State::Initialising:
        performInit();
        break;
    ...
}
```

**Compliant Example**:
```cpp
while (device.isBusy()) {
    // Spin-wait: polling is required; sleep would miss the 10µs window
}
```

### 5.7 Phase Labels for Complex Functions

**Requirement I-PHS-001**: Functions exceeding 100 lines or containing distinct processing phases SHOULD use inline phase labels to identify logical stages.

**Requirement I-PHS-002**: Phase labels MUST use a consistent format within a file.

**Compliant Example**:
```cpp
void processTransaction(Transaction& txn) {
    // [Validation Phase]
    if (!validateInputs(txn)) return;
    ...
    
    // [Transformation Phase]
    normaliseAmounts(txn);
    applyExchangeRates(txn);
    ...
    
    // [Persistence Phase]
    beginTransaction();
    ...
}
```

### 5.8 Concurrency and Synchronisation Comments

**Requirement I-CON-001**: Code that relies on specific synchronisation properties MUST document:
- Which locks protect which data
- The required lock ordering (if multiple locks are involved)
- Memory ordering requirements (for lock-free code)

**Requirement I-CON-002**: Functions that must be called from specific threading contexts (e.g., "main thread only", "must hold mutex X") MUST document this requirement.

**Compliant Example**:
```cpp
/// Updates the cached metrics.
///
/// Thread safety: Must be called with metrics_mutex_ held.
/// Must NOT be called from the metrics collection callback.
void updateMetricsCache() {
    ...
}
```

**Compliant Example**:
```cpp
// Acquire locks in canonical order: connection_lock_ < session_lock_ < buffer_lock_
// to prevent deadlock. See docs/threading.md for the complete ordering.
std::lock_guard<std::mutex> conn_guard(connection_lock_);
std::lock_guard<std::mutex> sess_guard(session_lock_);
```

### 5.9 Conditional Compilation Comments

**Requirement I-CND-001**: Conditional compilation blocks (preprocessor directives, feature flags, platform detection) MUST document:
- The condition being tested
- The behaviour difference between branches
- Any implications for testing or deployment

**Requirement I-CND-002**: The closing directive of a conditional block spanning more than 10 lines SHOULD include a comment identifying the corresponding opening condition.

**Compliant Example**:
```cpp
#ifdef TORCH_ENABLE_LLVM
// LLVM backend available: use JIT compilation for CPU kernels
bool cpu_fuser_enabled = true;
#else
// No LLVM: fall back to interpreter; performance will be reduced
bool cpu_fuser_enabled = false;
#endif

#if defined(__aarch64__)
    // Use ARM's native CRC32 instruction (ARMv8+).
    return __crc32d(crc, data);
#else
    // Software fallback for platforms without hardware CRC.
    return crc32_software(crc, data);
#endif // __aarch64__
```

---

## 6. Abstraction-Level-Specific Requirements

### 6.1 Level 0: Hardware-Proximate Code

The following requirements apply in addition to all general requirements when operating at Abstraction Level 0.

**Requirement L0-HW-001**: Register accesses MUST document:
- The register name and its function
- The significance of the value being written or the expected value being read
- Any timing constraints or sequencing requirements

**Compliant Example**:
```c
// Clear the interrupt pending flag by writing 1 to bit 7 of the status register.
// Must be done before re-enabling interrupts; see datasheet §5.2.3.
WRITE_REG(DEVICE_STATUS, 1 << 7);
```

**Requirement L0-HW-002**: Memory-mapped I/O operations MUST document volatility assumptions and compiler barrier requirements.

**Requirement L0-HW-003**: Bit manipulation operations MUST include comments identifying the bits by name or function, not merely by position.

**Non-Compliant Example**:
```c
config |= (1 << 4);
```

**Compliant Example**:
```c
config |= (1 << 4);  // Set ENABLE_DMA bit
```

**Compliant Example (Preferred)**:
```c
config |= ENABLE_DMA_BIT;  // Enable DMA transfers for this channel
```

**Requirement L0-HW-004**: Timing-sensitive code (spin loops, delays, memory barriers) MUST document the timing requirements and the basis for any delay values.

**Compliant Example**:
```c
// Delay required for PLL lock; minimum 100µs per datasheet §3.4.
// Using 150µs for margin against clock variance.
delay_microseconds(150);
```

### 6.2 Level 1: Systems Code

The following requirements apply in addition to all general requirements when operating at Abstraction Level 1.

**Requirement L1-SYS-001**: Functions that allocate resources MUST document ownership transfer semantics.

**Compliant Example**:
```cpp
/// Creates a new connection handle.
///
/// @return A handle that the caller owns. Caller must call
///         connection_destroy() when finished.
ConnectionHandle connection_create();
```

**Requirement L1-SYS-002**: Functions that may block MUST document the blocking behaviour and conditions for returning.

**Requirement L1-SYS-003**: Error handling paths that perform cleanup MUST document what resources are being released.

**Requirement L1-SYS-004**: Caching, pooling, or lazy initialisation strategies MUST be documented at the point of implementation.

**Compliant Example**:
```cpp
// Connection pool: maintains up to kMaxPoolSize idle connections.
// Connections idle for more than kIdleTimeoutMs are closed by the
// background reaper thread. Pool is thread-safe; access serialised
// by pool_mutex_.
class ConnectionPool { ... };
```

### 6.3 Level 2: Application Code

The following requirements apply in addition to all general requirements when operating at Abstraction Level 2.

**Requirement L2-APP-001**: Business logic that implements specific rules or policies MUST reference the source of those rules (requirements document, specification, user story, or similar).

**Compliant Example**:
```python
# Discount calculation per pricing policy v2.3:
# - 10% for orders over $100
# - Additional 5% for loyalty members
# See: docs/policies/pricing-v2.3.md
def calculate_discount(order, customer):
    ...
```

**Requirement L2-APP-002**: Configuration parameters with non-obvious values MUST document the rationale for the default value.

**Requirement L2-APP-003**: Integration points with external systems MUST document the expected interface contract and version.

---

## 7. Temporal and Maintenance Comments

### 7.1 TODO Comments

**Requirement T-TODO-001**: TODO comments MUST include at least one of:
- An identifier linking to an issue tracker (e.g., `TODO: bug 12345678`)
- An assignee or responsible party (e.g., `TODO(alice):`)
- A date or condition specifying when the TODO becomes actionable

**Requirement T-TODO-002**: TODO comments MUST NOT be used for tasks that should block the current change; such tasks must be completed before merge.

**Compliant Examples**:
```cpp
// TODO(#1234): Replace with constant-time comparison after security review.

// TODO: Remove after Q4 2024 migration completes.

// TODO(alice@): Refactor when upstream library exposes async API.

// TODO: bug 12345678 - Remove this after the 2047q4 compatibility window expires.
```

**Non-Compliant Example**:
```cpp
// TODO: Fix this
```

### 7.2 Deprecation Comments

**Requirement T-DEP-001**: Deprecated interfaces MUST include a comment specifying:
- The replacement interface or approach
- The timeline or version for removal (if known)

**Compliant Example**:
```cpp
/// @deprecated Use encryptAES256() instead. Will be removed in v3.0.
[[deprecated("Use encryptAES256() instead")]]
void encryptLegacy(const char* data);
```

### 7.3 Compatibility Comments

**Requirement T-CMP-001**: Workarounds for bugs in external dependencies MUST document:
- The external component and version affected
- A reference to the upstream bug report (if one exists)
- The conditions under which the workaround can be removed

**Compliant Example**:
```python
# Workaround for numpy #12345: einsum produces incorrect results for
# tensors with zero-length dimensions in versions < 1.22.
# Remove when minimum numpy version is bumped to 1.22+.
if has_zero_dimension(tensor):
    return fallback_implementation(tensor)
```

---

## 8. Comment Syntax and Formatting

### 8.1 General Formatting

**Requirement S-FMT-001**: Comments MUST be grammatically correct sentences with appropriate capitalisation and punctuation, except for:
- Single-word annotations (e.g., `// fallthrough`)
- Structured annotation tags (e.g., `@param`, `@return`)
- Phase labels (e.g., `// [Validation Phase]`)

**Requirement S-FMT-002**: Multi-line comments MUST have consistent indentation aligned with the code they describe.

**Requirement S-FMT-003**: Comments MUST NOT extend beyond the project's defined line length limit. If no limit is defined, comments SHOULD NOT exceed 100 characters per line.

### 8.2 Language-Specific Syntax

**Requirement S-SYN-001**: Within a single file, comment delimiter style MUST be consistent. Mixing styles (e.g., `//` and `/* */`) is permitted only when:
- One style is used for documentation comments and another for implementation comments, AND
- This usage is consistent throughout the project

**Requirement S-SYN-002**: The project MUST define and document which comment syntax is used for which purpose.

### 8.3 Namespace and Scope Closing Comments

**Requirement S-CLS-001**: Closing braces for namespaces, classes, or other scopes spanning more than 20 lines SHOULD include a comment identifying the scope being closed.

**Compliant Example**:
```cpp
namespace nvfuser {
namespace python_frontend {
...
} // namespace python_frontend
} // namespace nvfuser
```

### 8.4 Prohibited Content

**Requirement S-PRO-001**: Comments MUST NOT contain:
- Profanity or derogatory language
- Personal complaints or blame attribution
- Confidential information not appropriate for the repository
- Commented-out code exceeding 10 lines (use version control instead)

**Requirement S-PRO-002**: Comments MUST NOT include ASCII art exceeding 5 lines unless serving a genuine diagrammatic purpose.

---

## 9. Project-Level Configuration

### 9.1 Style Guide Layering

**Requirement P-STY-001**: Projects MAY adopt different commenting conventions for different subsystems or visibility layers.

**Requirement P-STY-002**: When multiple conventions are used within a project, the project MUST document:
- The boundary between convention domains (e.g., "public API vs internal implementation")
- A mapping from subsystem or layer to applicable convention
- The rationale for the differentiation

**Compliant Example (Project Documentation)**:
```markdown
## Commenting Conventions

This project uses layered commenting conventions:

| Layer | Convention | Rationale |
|-------|------------|-----------|
| Public API (`include/`) | Google C++ Style Guide | Consistency with ecosystem |
| Internal (`src/`) | LLVM/MLIR Style Guide | Alignment with MLIR codebase |
| Tests (`test/`) | Minimal | Test names should be self-documenting |
```

### 9.2 Toolchain Documentation

**Requirement P-TOOL-001**: Projects declaring conformance to this specification SHOULD document:
- The specific tools used to verify comment-related requirements
- The tool configuration (e.g., clang-tidy checks enabled, documentation coverage thresholds)
- The CI/CD integration point where verification occurs

**Compliant Example (Project Documentation)**:
```markdown
## Comment Verification

The following tools enforce commenting standards:

- **clang-format**: Enforces comment line length and alignment
- **clang-tidy**: `readability-*` checks for comment quality
- **Doxygen**: Generates documentation; warnings treated as errors for public APIs
- **lintrunner**: CI integration via `lintrunner --take CLANGFORMAT,CLANGTIDY`

Verification runs on every pull request via GitHub Actions.
```

---

## 10. Compliance Verification

### 10.1 Automated Verification

The following requirements are amenable to automated verification:

| Requirement | Verification Method |
|-------------|---------------------|
| F-LIC-001, F-LIC-003 | Static analysis of file header |
| D-PUB-001 | Documentation coverage tools |
| D-CON-002 | Regex patterns for magic numbers |
| T-TODO-001 | Regex pattern: `TODO` without identifier |
| S-FMT-003 | Line length linting |
| S-PRO-002 | Line count in comment blocks |
| S-CLS-001 | AST analysis for scope length |

### 10.2 Review-Based Verification

The following requirements require human review:

| Requirement | Review Criteria |
|-------------|-----------------|
| I-WHY-001 | Does comment explain *why*, not *what*? |
| I-ALG-001 | Is algorithm complexity documented? |
| I-DOM-001 | Is domain notation explained? |
| L0-HW-001 | Are hardware interactions explained? |
| L2-APP-001 | Are business rules traceable? |

### 10.3 Compliance Levels

Projects MAY declare conformance at one of three levels:

**Level A (Full Conformance)**: All MUST and SHOULD requirements are satisfied.

**Level B (Core Conformance)**: All MUST requirements are satisfied.

**Level C (Partial Conformance)**: A documented subset of MUST requirements is satisfied, with a roadmap to achieve Level B.

---

## 11. Exceptions and Waivers

### 11.1 Project-Level Exceptions

**Requirement E-PRJ-001**: Projects MAY define exceptions to specific requirements, provided:
- The exception is documented in a project-level configuration file
- The rationale for the exception is recorded
- The exception does not contradict a MUST NOT requirement

### 11.2 File-Level Exceptions

**Requirement E-FIL-001**: Individual files MAY be excluded from specific requirements by including a machine-readable exception marker.

**Example**:
```cpp
// comment-standard: disable=D-PUB-001 (generated file)
```

### 11.3 Inline Suppression

**Requirement E-INL-001**: Individual comment requirement violations MAY be suppressed inline, provided:
- The suppression marker immediately precedes the violation
- A brief justification accompanies the suppression

**Example**:
```cpp
// comment-standard: suppress=D-CON-002 (value from external protocol spec)
const int MAGIC_HANDSHAKE = 0x47525043;
```

---

## Annex A: Language-Specific Syntax Reference (Informative)

| Language | Single-line | Multi-line | Documentation |
|----------|-------------|------------|---------------|
| C/C++ | `//` | `/* */` | `///`, `/** */` |
| Python | `#` | `""" """` (strings) | `""" """` (docstrings) |
| Rust | `//` | `/* */` | `///`, `//!` |
| Java | `//` | `/* */` | `/** */` (Javadoc) |
| JavaScript/TypeScript | `//` | `/* */` | `/** */` (JSDoc) |
| Go | `//` | `/* */` | `//` (godoc convention) |
| Haskell | `--` | `{- -}` | `-- |`, `{- | -}` (Haddock) |
| Shell | `#` | `: ' '` (heredoc) | N/A |

---

## Annex B: Comment Density Guidelines (Informative)

The following table provides guidance on expected comment density by context:

| Context | Expected Density | Primary Focus |
|---------|------------------|---------------|
| Public API headers | High | Complete Doxygen/JSDoc coverage |
| Implementation files | Moderate | Explain "why", not "what" |
| Algorithmic sections | High | Step-by-step with complexity |
| Domain-specific code | High | Notation mapping, references |
| Obvious/trivial code | Minimal or none | Trust self-documenting code |
| Build configuration | Moderate | Options and their effects |
| Test files | Low | Test names should be descriptive |

---

## Annex C: Section Divider Formats (Informative)

Projects should select one format and use it consistently:

**Equals Style** (common in systems code):
```cpp
// ============================================================================
// Section Name
// ============================================================================
```

**Box Style** (common in C):
```c
/******************************************************************************
 * Section Name
 ******************************************************************************/
```

**Dash Style** (minimal):
```cpp
// --- Section Name ---
```

**Hash Style** (common in Python/Shell):
```python
# =============================================================================
# Section Name
# =============================================================================
```

---

## Annex D: Compliance Checklist (Informative)

### File-Level
- [ ] License header present and correctly formatted (F-LIC-*)
- [ ] File purpose statement present if file exceeds 100 lines (F-PUR-001)
- [ ] Section comments used if file exceeds 500 lines (F-STR-001)

### Declaration-Level
- [ ] All public interfaces have documentation comments (D-PUB-001)
- [ ] Documentation includes brief, parameters, return value, exceptions (D-PUB-002)
- [ ] Type invariants documented (D-TYP-003)
- [ ] Constants and magic numbers explained (D-CON-*)

### Implementation-Level
- [ ] Comments explain "why", not "what" (I-WHY-001)
- [ ] Domain notation defined and mapped (I-DOM-*)
- [ ] Algorithms identified by name and complexity (I-ALG-001)
- [ ] Concurrency requirements documented (I-CON-*)
- [ ] Conditional compilation explained (I-CND-001)

### Abstraction-Specific
- [ ] (L0) Hardware registers and timing documented
- [ ] (L1) Resource ownership documented
- [ ] (L2) Business rules traceable to requirements

### Maintenance
- [ ] TODOs have identifiers, assignees, or dates (T-TODO-001)
- [ ] Deprecations specify replacements (T-DEP-001)
- [ ] Workarounds have upstream references (T-CMP-001)

### Project Configuration
- [ ] Style guide layering documented if applicable (P-STY-*)
- [ ] Verification toolchain documented (P-TOOL-001)

---

## Annex E: Example Compliant File (Informative)

```cpp
// SPDX-FileCopyrightText: Copyright (c) 2024 Example Corporation
// SPDX-License-Identifier: Apache-2.0

// Implements rate limiting for the API gateway using a sliding window
// algorithm. Integrates with the metrics subsystem for limit breach
// reporting. See docs/architecture/rate-limiting.md for design rationale.

#ifndef GATEWAY_RATE_LIMITER_H
#define GATEWAY_RATE_LIMITER_H

#include <chrono>
#include <mutex>
#include <string>

namespace gateway {

/// Configuration for a rate limiting policy.
struct RateLimitConfig {
    /// Maximum requests permitted within the window.
    int max_requests;
    
    /// Duration of the sliding window. Must be positive.
    std::chrono::milliseconds window_duration;
};

/// Enforces rate limits on incoming requests using a sliding window algorithm.
///
/// Thread safety: All public methods are thread-safe. Internal state is
/// protected by a mutex.
///
/// Memory: O(n) where n is the maximum number of tracked clients.
class RateLimiter {
public:
    /// Constructs a rate limiter with the specified policy.
    ///
    /// @param config  The rate limiting policy to enforce.
    /// @throws std::invalid_argument if window_duration is not positive.
    explicit RateLimiter(RateLimitConfig config);
    
    /// Attempts to acquire permission for a request.
    ///
    /// @param client_id  Identifier for the requesting client.
    /// @return True if the request is permitted; false if rate limited.
    ///
    /// @note This method updates internal state; even "checking" the limit
    ///       has side effects on the sliding window.
    bool tryAcquire(const std::string& client_id);
    
    /// Returns the number of requests remaining for a client in the current window.
    ///
    /// @param client_id  Identifier for the client to query.
    /// @return Remaining request count; may be zero if limit reached.
    int remainingQuota(const std::string& client_id) const;

private:
    // Sliding window algorithm: track timestamps of requests within the window.
    // Older entries are lazily evicted on each tryAcquire call.
    // O(1) amortised insertion, O(k) query where k = requests in window.
    // Reference: https://konghq.com/blog/how-to-design-a-scalable-rate-limiting-algorithm
    
    RateLimitConfig config_;
    
    // Protects all mutable state below. Acquired for both read and write
    // operations to ensure consistency of the sliding window invariant.
    mutable std::mutex mutex_;
    
    // Maps client_id to the timestamps of their requests within the window.
    // Invariant: all timestamps are within [now - window_duration, now].
    std::unordered_map<std::string, std::vector<TimePoint>> request_history_;
    
    // Removes expired entries from a client's history.
    // Precondition: mutex_ must be held.
    void evictExpiredEntries(std::vector<TimePoint>& history);
};

}  // namespace gateway

#endif  // GATEWAY_RATE_LIMITER_H
```