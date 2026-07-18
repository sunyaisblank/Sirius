# Sirius C++26 Style Guide

This guide governs all C++ in the rebuilt Sirius. It layers three sources in priority order: the rules in this document, then the C++ Core Guidelines for semantic questions this document does not answer, then the Google C++ Style Guide for mechanical questions neither answers. Formatting is not debated in review at all; `.clang-format` at the repository root is the authority, and a file that round-trips through clang-format unchanged is correctly formatted. The measured toolchain baseline is GCC 14.2, Clang 21.1, and MSVC 19.4x, all in `-std=c++2c` mode or equivalent; every rule here compiles on that baseline today.

## 1. The safety posture

The purpose of adopting C++26 early is compile-time safety, so the default posture is: make invalid states unrepresentable at compile time where the type system can, assert them at runtime where it cannot, and never accept them silently.

- Every function that documents a precondition states it as a contract macro (section 2), not as a comment.
- Narrowing conversions are illegal outside a named cast; initialise with braces so the compiler enforces this.
- `enum class` always; a plain `enum` requires a recorded justification (interop with a C API is the only anticipated one).
- A `bool` parameter that selects between behaviours is forbidden; use a two-value `enum class` whose name reads at the call site. `render(scene, Tiling::kEnabled)` reviews itself; `render(scene, true)` does not.
- Numeric identifiers with distinct semantics take distinct types. `MetricId` stays a closed enum; tile indices, pixel indices, and byte offsets do not interchange.
- Public template parameters are constrained by concepts. An unconstrained `typename` in a public interface is a review error.
- `[[nodiscard]]` on every function whose return value is the point of calling it, which in this codebase is nearly every function in Core.

## 2. Contracts

C++26 contracts (P2900) are the intended end state, but the measured toolchain provides them only as GCC's experimental attribute syntax, so contract expressions are written once behind macros and the macros track the toolchain. `sirius/base/contracts.h` defines:

- `SIRIUS_PRE(expr)` and `SIRIUS_POST(expr)` for preconditions and postconditions, placed at the top of the function body (P2900's declaration syntax is adopted mechanically when both compilers ship it; the macro spelling is chosen so that migration is a script, not a review).
- `SIRIUS_ASSERT(expr)` for invariants mid-function, the analogue of `contract_assert`.
- `SIRIUS_AXIOM(expr)` for conditions too expensive to evaluate even in checked builds; documentation that a checker may exploit, never executed.

Semantics follow build mode: Debug and test builds enforce (violation terminates with location and expression text), Release builds of the CPU oracle observe (log and continue, because the oracle's job is to measure error, not to hide it), and Release hot paths ignore. The mode is a CMake-level decision, never a per-file one. When `__cpp_contracts` reports P2900 on both release compilers, the macros compile to native syntax and this section shrinks to the mapping table.

Contract expressions state the caller's observable obligations, not implementation detail: preconditions are the minimum the function needs, and a precondition stronger than necessary is itself a defect, because it restricts callers without benefit. Data crossing a module boundary is validated at that boundary once; a downstream re-check of an upstream guarantee duplicates authority and is removed.

## 3. C++26 feature rules

Tier one features are used directly wherever they apply.

- Pack indexing (`Ts...[0]`) replaces recursive template peeling. If you are writing a recursive template to reach the Nth element, stop.
- Deducing this replaces CRTP for static polymorphism and replaces const/non-const method duplication.
- `= delete("reason")` is the required form of every deleted overload; the reason names the alternative ("use MetricRegistry::lookup; string names drift").
- `static_assert` always carries a user message stating the violated requirement, not restating the condition.
- Placeholder variables (`auto _ =`) for RAII guards whose name would never be read.

Tier two, contracts, is section 2. Tier three, static reflection (P2996), is absent from the toolchain: no code in the repository may require it. The two subsystems that will eventually want it, the kernel constant tables and CLI/config binding, are structured as generated code (`*_gen.h`, emitted at build time by scripts under `scripts/` from a single declared source of truth) so that when reflection arrives it replaces a generator, not an architecture. Handwritten edits to generated files are build errors; the generator re-runs and CI diffs the output.

Library features: `std::expected` (section 4), `std::span` and `std::mdspan` at every boundary that today would take pointer plus length or a raw 2D index computation, `std::print` in tools and tests (the render CLI keeps FTXUI), `<numbers>` constants over literals. `std::simd` and senders/receivers are not yet in either standard library; the threading and SIMD abstractions live behind seams in `sirius/base/` so they can adopt the standard forms when they ship.

## 4. Errors

Errors classify by recoverability, and each class has exactly one channel.

- Fallible operations (file IO, device discovery, configuration parsing, kernel compilation) return `std::expected<T, Error>` where `Error` carries origin, operation, and offending input. Exceptions do not cross module boundaries and are not used on the physics or render paths at all.
- Invariant violations (a null direction vector inside the integrator, an unnormalised camera frame) are contract violations: they halt in checked builds and are never converted to error returns, because continuing past a broken invariant produces indeterminate output that looks like a render.
- Declined work is not an error channel misuse: on a live path with absent data, the unit of work declines loudly (the July failure philosophy), and fabricating a stand-in is permitted only in explicitly sandboxed modes.

## 5. Names

- Files: `snake_case.h` / `snake_case.cpp`; kernels `snake_case.slang`; a test file is `<subject>_test.cpp`. The codename scheme (`PHMT200A.h`) is retired; `docs/ARCHITECTURE.md` holds the mapping.
- Types and functions: `PascalCase`. Cheap accessors may be `snake_case` per Google practice.
- Variables `snake_case`; data members `snake_case_`; compile-time constants `kPascalCase`.
- Namespaces: everything lives under `sirius::`, one nested level per layer (`sirius::core`, `sirius::render`, `sirius::backend`, `sirius::app`). No `using namespace` in headers, no aliases that flatten layers.
- Physics symbols may use domain notation where the referenced equation makes single letters clearer than words (`g`, `Gamma`, `E`, `Lz`, `Q`), provided the declaration cites the equation or reference. `MTW` and paper citations belong at the declaration, not in every use.

## 6. Headers and dependencies

Headers are self-contained, include what they use, and prefer forward declarations where the full definition is not structurally required. `#pragma once` replaces include guards (a recorded deviation from Google: every target compiler supports it, and it removes a class of copy-paste guard collisions). Include order: matching header, project headers, third-party, standard library; clang-format enforces the grouping. A target links only what it uses; the dependency graph is acyclic with call direction from orchestration toward primitives, and a lateral include between peer modules is a layering violation to resolve, not to wave through. C++ modules are not adopted this engagement: MSVC/GCC/Clang module interop with CMake and with Slang-adjacent codegen is not yet stable enough, and headers with strict IWYU discipline are the fallback that loses least.

## 7. Floating point and numerics

- Core and the oracle compute in `double`. Kernels compute in `float` with the precision ladder the architecture document defines; the choice per path is a design decision recorded there, never an accident of literal suffixes.
- Every tolerance, epsilon, and physical constant comes from the constants authority (`sirius/core/constants.h`); a numeric literal standing in for a named constant is a review error, and a new unexplained constant is a defect.
- Accumulations whose length scales with image size or step count use compensated (Kahan or pairwise) summation; the naive form is permitted only with a recorded bound on accumulated error.
- Subtractions of nearly equal quantities are restructured analytically where possible; where not, the site carries a comment naming the cancellation and the retained precision.
- No fast-math flags anywhere in Core or the oracle. Kernel builds may enable finite-math-only per backend when the parity suites pass under it.

## 8. Concurrency

Prefer eliminating shared mutable state (per-tile ownership, message passing) over synchronising it. Where synchronisation is necessary use the narrowest primitive that suffices, and document lock acquisition order in the owning module's header. Render work is compute-bound and decomposes by tile; IO (progress reporting, image writing) stays on its owner thread. Data races are undefined behaviour to be eliminated by design, and ThreadSanitizer runs in CI on the threaded suites.

## 9. Comments and exposition

A comment states why, or states a constraint the code cannot express (units, frames, coordinate charts, citation of the governing equation); a comment that restates the code is noise and is removed in review. Every physics routine's declaration cites its source (paper, equation number, or textbook section) so a reviewer can check the mathematics against the reference rather than reverse-engineering intent. Decision records from the July remediation and from this rebuild (deliberate constraints that look like bugs) stay at the site as one-line comments naming the decision.

## 10. Tests

Tests are named by the postcondition they verify, not the method they call (`ShadowBoundaryMatchesBardeenAnalytic`, not `TestComputeShadow`). Every suite carries a generated ctest label; the generator script is the labelling authority and CI rejects stale label files. A test that asserts nothing is a defect. Physics suites distinguish software correctness (the code does what the design says) from model validity (the design matches the reference solution), and every numerical suite states its tolerance from the constants authority with a comment naming why that tolerance is achievable.
