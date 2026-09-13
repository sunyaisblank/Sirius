# Source closeout — 13 September 2026

The canonical workspace is `Sirius/`. The accumulated development is preserved
on `development/source-closeout-2026-09-13`, based on upstream main at
`8c7ba5a`. This is an unfinished development checkpoint, not a release.

The recovered history includes `2f7189b` (governed evidence and Vulkan session
output), `4a685bb` (subpixel point-source filtering), and `b5ef145` (the surviving
uncommitted coupled-transport and source-sky implementation). The latter
preserves all 90 changed/new files from the former completion workspace,
including CPU/source-sky and Kerr-infinity work, camera differentials,
four-column GPU continuation, dispatch changes, and their existing tests.

The old canonical checkout had remained at `072794b` with an earlier dirty
snapshot. That snapshot is retained locally as
`refs/archive/original-checkout-2026-09-13`; it was not reapplied over newer work.
Recovered branch tips also remain under `refs/archive/recovered-2026-09-13/`.
A verified recovery bundle, historical session notes/scripts, and the cleanup
manifest are retained locally in `.git/closeout/`. They are recovery material,
not build inputs or qualification evidence.

All 42 sibling Sirius/test-toolchain directories were removed after checking
that their Sirius commits and uncommitted source had been preserved. Obsolete
builds and render outputs were removed. New build output stays in
`bin/linux-gcc`; other projects were left untouched.

## Remaining implementation

The September 10 handoff reported incomplete retained-precision observer-frame
and camera transport, joint admission of the central ray and four physical
film/pupil derivatives, dense event handling, and physical detector integration.
It also reported unresolved CPU critical refusals and output workflows,
full-workload image quality/performance, and native platform qualification.
Those obligations remain open; this closeout does not claim to resolve them.

The handoff referenced staged retained-pair arithmetic and camera-direction
modules under the former `.sirius-release-work/evidence/` directory. That
directory and those modules were already absent when this closeout began.
Their reported numerical results cannot be reproduced from the surviving
checkout alone. Recover the original artifacts, or reimplement and independently
validate that stage, before relying on it. Historical success counts are not
current validation.

No external operating domain was admitted at configure time (0/8). Physical
Radeon, WSL2/Dozen, native Windows/macOS build and runtime, native viewer input,
and the exact IMAX workload still require independent qualification on a single
final revision. Release packaging remains disabled.

## Closeout validation

The recovered checkpoint compiled every application, shader (with SPIR-V validation),
and test target with the Linux GCC preset and warnings treated as errors.
Repository structure, operating-model validation and negative controls, build-policy
negative controls, and generated CTest labels passed. The baseline Mandatory run
passed its first 678 tests; it was stopped during the Vulkan parity tests before
resuming implementation. This is not a complete Mandatory pass. CI formatting
was then applied to 22 first-party files. Logs remain in `.git/closeout/`.

The owner expanded the task to complete the remaining renderer development.
The gaps above are the starting point for that continuing work, not a final
delivery declaration.
