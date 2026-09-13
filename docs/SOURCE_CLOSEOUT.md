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

## Implementation completed during closeout

The CPU launch now differentiates the actual smooth film and Cartesian pupil
coordinates, including the metric frame rebuilt at the displaced pupil event.
All four physical columns enter joint transport with their units preserved.
Kerr-family metric Hessians now come from the same geometric authorities as
the nominal metric, replacing repeated numerical derivative sampling where
the analytic Hessian is represented. Endpoint sampling preserves small
covariant variations directly, and the shared Hamiltonian contractions reduce
repeated work across the four columns.

Polarisation reconditioning now uses the regular live chart's Eulerian gauge.
This resolves the reproduced past-horizon capture failure without changing the
central trajectory. The live camera family is also used by the numerical
failure diagnostic reference.

The complete 487-test core executable passed. Focused transport/source tests
passed, as did the original CPU PPM, PNG/EXR and polarisation output workflows.
The PNG/EXR run tested `0b2520e`; the PPM run tested `4514adb`. These results do
not constitute a full Mandatory pass on the eventual final revision.

The CPU point-source detector is now connected to the actual continuous camera
and worker-local tracer for uncharged, zero-cosmological-constant Kerr-family
metrics. It preserves the original Gaussian packet and pupil through adaptive
subdivision, retraces candidate image roots, retains multiple images, and applies
their individual frequency and transmission. Failed or exhausted packets cannot
publish partial radiance. Other metric families retain their existing route.
Forward source-coordinate bounds avoid unnecessary refinement around empty fold
regions without increasing the work cap or relaxing the radiance tolerance.

Seven detector tests pass, including rotated folds, a close polynomial image
pair, disconnected visibility, shared ownership and cancellation/exhaustion.
A moving ThinLens Kerr frame with 100,000 stars produces finite, nonconstant
linear EXR output, identical with one and two rendering workers. Display grading
suppressed this deliberately faint fixture; its linear radiance is nonzero.
All 488 core tests pass after the camera tangent-basis correction that preserves
exact geometric zeros in continuous detector offsets. ThinLens and celestial
basis probes also pass on the pinned Radeon/Dozen route. These checks establish
the tested CPU connection, not full-scene convergence or throughput.

The saved CPU critical-ray failures now have a governed regression. Portable
twofold working arithmetic retains central momentum, the four Hamiltonian
variations, metric/gradient values and covariant conversions. Each derivative
uses its actual central stage tangent and shares that stage's private geometry.
The defining DP tableau's constant-field row sums are preserved explicitly.
The public state types and the original error budgets remain unchanged.

All fourteen original launch packets are preserved and replayed without changing
their directions, observer boosts, tolerances or 30,000-attempt limit. Twelve
reach physical outcomes; the other two explicitly exhaust that limit. None
terminates with the reproduced interpolation failure. Ten separately derived
minimum-step states agree with independent 75-digit RK4 refinement, and 960
retained metric/inverse/gradient fields agree with independent 75/110-digit
witnesses. The latter also detect scalar-only substitutions. These are finite
numerical witnesses, not certified interval enclosures.

All 491 core tests and all fourteen coupled-transport tests pass on this source,
including physical neighbouring-ray, event and refinement controls. The full
coupled suite takes about eight minutes on the current machine. These checks do
not constitute a complete Mandatory or full-workload qualification.

## Remaining implementation

The September 10 handoff reported incomplete retained-precision observer-frame
and camera transport, joint admission of the central ray and four physical
film/pupil derivatives, dense event handling, and physical detector integration.
The GPU production migration, GPU detector integration, full-workload image
quality/performance and native platform qualification remain open. The CPU
corrections above resolve the reproduced interpolation failures in the saved
launch regression; two original packets still require more work than their
unchanged attempt limit allows.

The handoff referenced staged retained-pair arithmetic and camera-direction
modules under the former `.sirius-release-work/evidence/` directory. That
directory and those modules were already absent when this closeout began.
Nine retained camera shaders, including the later complete metric/frame/launch
prototype and recorded arithmetic repairs, were subsequently recovered from
session tool output into `tests/support/retained_camera`. All three probes
compile in both narrow modes and pass SPIR-V and arithmetic-control checks.
The recovered independent eight-case reference also passes its high-precision
stability and null/frequency identity checks. The twelve original input packets
were subsequently recovered exactly and their helper-context coefficients were
observed on the pinned device. The monolithic complete camera exceeded 21 GiB
of host memory during pipeline preparation and was stopped before dispatch.

A bounded replacement camera stage now executes 5,292 retained arithmetic
instructions using 227 live registers. Both narrow-mode backend tests pass all
20 complete 104-value camera fixtures, 16 total invalid-request refusals,
4,160 component mutation controls and 40 whole-packet low-part deletion controls.
Its 100/180-digit independent fixture regeneration is byte-for-byte reproducible.
The observed stage needs 111,560 explicit buffer bytes and prepares in under
one second on the pinned Radeon/Dozen route. This resolves the camera prototype's
compiler expansion for these tests. Production transport, dense events,
continuation and detector integration still require their own implementation
and validation; historical success counts are not current evidence.

The camera and joint seven-stage Hamiltonian DP pair now have fixed-capacity
batched device interfaces and build-embedded programs. The transport program
uses 2,587 instructions and 177 registers; all four physical columns use their
actual central stage. Camera upload retains continuous film and pupil offsets,
including independently checked offsets below one binary32 spacing. Invalid
rows cannot publish previous or partial candidates, and repeated submissions
reuse their governed device buffers.

Transport uses three binary32 terms and a separate arithmetic radius. The
two-term prototype exceeded the original projected critical-ray budgets by
factors of 1.5 to 9.44. The three-term results agree with the independent
projected refinement to at most 1.3e-7 of those unchanged budgets. Optimized
three-term submission/readback takes 163 ms for 2,048 private candidates on the
pinned Radeon/Dozen route after pipeline preparation. This is a stage timing,
not complete renderer throughput. The camera uses two working terms; no input
correction is silently dropped on conversion.

All six retained backend tests pass, covering the recovered camera cases,
continuous inputs, twelve independent phase-space fixtures, low/tail deletion,
arithmetic enclosures and invalid-row reuse. Both new frozen reference sets
regenerate byte for byte. All application and test targets build with warnings
as errors. These device stages still require endpoint projection, accepted
continuation, dense events and detector connection before replacing the live
GPU trace path. Their arithmetic radii do not claim global ODE enclosures.

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
was then applied to 22 first-party files. A later Mandatory run at `f76db09`
passed its first 687 tests, including physical Radeon camera and metric probes,
before being stopped to resume development. Logs remain in `.git/closeout/`.

The owner expanded the task to complete the remaining renderer development.
The gaps above are the starting point for that continuing work, not a final
delivery declaration.
