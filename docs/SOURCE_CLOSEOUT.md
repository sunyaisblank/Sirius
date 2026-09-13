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

## Renderer integration

Kerr-family Vulkan frames now use bounded retained stages for camera launch,
physical phase initialization, joint seven-stage Hamiltonian transport, null
projection and dense admission. All four physical film/pupil derivatives share
the central stages. Full, embedded, midpoint and refined candidates remain
private until the common component budgets admit them. The shared host tracer
then localizes events and owns source and physical detector evaluation. Each
worker's retained phase survives accepted intervals and explicit rollback;
public binary64 fields remain views of that phase.

Both narrow precision selections use three binary32 transport terms and a
separate arithmetic radius. The camera uses two working terms. FP64 selection
uses the same retained records with exact binary64 products; the high part is
rounded by an explicit ordered binary32 multiplication. A native conversion
from the wide product failed the complete camera's residual checks on the
Radeon/Dozen route. The explicit split passes the independent references.
Arithmetic radii are not global ODE certificates.

The generated programs and validated SPIR-V are embedded in the product. All
six stages reuse fixed buffers. The renderer checks its allocation plan against
actual device buffer bytes and coalesces at most 64 worker requests. Timing
feedback adjusts subsequent batches within any safety ceiling; the existing safety duration remains in
force when the soft target is disabled. Source textures and catalogues stay
with the shared host source owner. The external frame publishes only after all
samples complete, preserving its previous image on cancellation or failure.
Tiny images reserve the minimum scratch tile without requiring eight actual
image pixels along each axis.

Exact Minkowski formulas avoid curved-metric evaluation while preserving the
DP tableau, differentiated null projection and general Hermite polynomial.
Independent nonlinear flat dense fixtures prevent treating an arbitrary flat
segment as a straight line. The detector applies its existing forward-range
bound to catalogue queries and refines coarse mixed-visibility cells before
spending catalogue visits. Neither change increases a work cap or relaxes a
radiance tolerance.

The device reference sets include 23 complete smooth-camera cases, twelve
Hamiltonian phase fixtures, fourteen endpoint/projection cases and 39 dense
samples. Camera and dense fixture regeneration is byte-for-byte reproducible.
The FP64 stage/controller test passes, as do the narrow joint-controller and
shared-tracer continuation/rollback checks. All seven detector regressions
pass. A complete moving ThinLens flat frame passes finite radiance, actual
allocation, single publication, cancellation and insufficient-budget checks.
The moving ThinLens Kerr/100,000-star CPU–device image comparison is undergoing
validation; this document does not yet declare that test or full Mandatory
complete.

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
