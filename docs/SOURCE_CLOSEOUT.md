# Source closeout — 13 September 2026

The canonical workspace is `Sirius/`. The accumulated development is preserved
on `development/source-closeout-2026-09-13`, based on upstream main at
`8c7ba5a`. This closeout consolidates source, renderer integration and recovery
history. Release packaging remains disabled pending operating-domain admission.

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
Arithmetic radii are not global ODE certificates or truncation estimators.
Embedded and independent-refinement comparisons of complete expansion centers
retain the existing local component budgets. After admission, the unchanged
three-term expansion defines the next local numerical initial value. Cached
candidate records keep their bounds and every private substage propagates them;
the next local problem starts a fresh arithmetic enclosure. This avoids treating
repeated coordinate-transform interval wrapping as global trajectory uncertainty
without rounding away retained limbs.

Host comparisons subtract all six input limbs before reducing the difference
to binary64. Lower-order increments use the same exact expansion accumulator
and include discarded terms in an outward arithmetic radius. This preserves
sparse tails even when their exponent span exceeds two doubles; a regression
checks a represented difference of 2^-120 beside common leading terms.

A strengthened continuation regression exposed two host handoff defects:
reconstruction of displacement columns from rounded increments and inclusion of
the rounded affine clock in the phase cache key. Accepted full endpoints now
remain authoritative, and the autonomous phase cache excludes that clock. Four
concurrent flat/Kerr traces pass with exactly one initialization per trace;
rollback and cancellation retain the correct phase.

The generated programs and validated SPIR-V are embedded in the product. All
six stages reuse fixed buffers. The renderer checks its allocation plan against
actual device buffer bytes and coalesces at most 64 worker requests. Timing
feedback adjusts subsequent batches within any safety ceiling; the existing safety duration remains in
force when the soft target is disabled. Source textures and catalogues stay
with the shared host source owner. The external frame publishes only after all
samples complete, preserving its previous image on cancellation or failure.
Tiny images reserve the minimum scratch tile without requiring eight actual
image pixels along each axis.

Each ray now owns a cooperative workgroup. The generated arithmetic DAG
schedules up to 64 independent expressions per layer without changing their
operands or summation order. Its allocator checks input and output ownership
across barriers; shared registers and status stay below the portable 16 KiB
limit. Independent RK state components run in parallel as well. Null projection
executes only its metric/tangent prefix before selecting a root, then evaluates
the complete physical output. Fixed buffer strides stay unchanged while only
requested rows are dispatched. All eleven retained backend/value tests pass
with cooperative arithmetic, RK component evaluation and corrected continuation,
including FP64 (160.743 seconds). The 491 core tests also pass. The moving
ThinLens Kerr CPU reference completes with exact serial/two-worker equality
(4,768.036 seconds).

The complete moving ThinLens Kerr frame agrees with the independently completed
CPU reference at a relative linear-RGB L1 error of 8.5437703874154113e-9, against
the unchanged .02 limit. All 32 RGBA values are finite and the reference signal
is nonzero. On the Radeon 780M through WSL2/Dozen, its eight pixels take
11,228.8 seconds across 1,338,417 retained dispatches. The maximum measured
submission is 469.421 ms, with seven soft subdivisions and zero safety fallbacks.
This focused measurement uses the renderer source committed in `bf30294` and
the preceding development build configuration; it is not a final-revision
Mandatory receipt. Throughput on this route is limited, and this small image
does not establish full-workload performance.
The original image comparison test also completed its own CPU reference and
passed in 14,574.406 seconds; the GPU measurement above is its device phase.

Exact Minkowski formulas avoid curved-metric evaluation while preserving the
DP tableau, differentiated null projection and general Hermite polynomial.
Independent nonlinear flat dense fixtures prevent treating an arbitrary flat
segment as a straight line. The detector applies its existing forward-range
bound to catalogue queries and refines coarse mixed-visibility cells before
spending catalogue visits. Neither change increases a work cap or relaxes a
radiance tolerance.

The device reference sets include 23 complete smooth-camera cases, fifteen
Hamiltonian phase fixtures, seventeen endpoint/projection cases and 48 dense
samples. Three independent nonradial weak-field cases add half-unit and unit
steps in both regular charts, with angular and pupil columns. The original
critical and flat witnesses remain unchanged. Fixture regeneration is
byte-for-byte reproducible.
All eleven retained backend/value tests pass with these expanded fixtures and
the sparse-tail correction (164.613 seconds), including FP64, independent
midpoint/refinement admission and the exact flat lower-increment regression.
The FP64 stage/controller test passes, as do the narrow joint-controller and
shared-tracer continuation/rollback checks. All seven detector regressions
pass. A complete moving ThinLens flat frame passes finite radiance, actual
allocation, single publication, cancellation and insufficient-budget checks.
The moving ThinLens Kerr/100,000-star CPU–device image comparison described
above measures the complete physical detector path. These focused results do
not substitute for the complete revision-bound Mandatory gate.

No external operating domain was admitted at configure time (0/8). Physical
Radeon, WSL2/Dozen, native Windows/macOS build and runtime, native viewer input,
and the exact IMAX workload still require independent qualification on a single
final revision. Release packaging remains disabled.

## Closeout validation

The complete-gate rehearsal exposed an exact finite-causal endpoint that was
left unowned by a strict outside-only predicate. The tracer now admits contact
at the accepted endpoint; the finite-pupil regression escapes in one interval.
The causal locator and its caller use the same radius evaluation. All thirteen
CPU boundary tests pass, including live disk-profile, bundle, and horizon checks.

The attestation exporter now takes its generated-input inventory from the build
gate authority, so retained-camera fixtures accompany every qualification
bundle. Its complete false-evidence controls pass. CPU arithmetic avoids an FMA
only for an exact zero factor, preserving signed zero and invalid-input behavior;
all 491 core tests pass, including the independent retained metric witnesses.
The installed 128x128 CPU smoke image keeps its original numerical settings and
600-second timeout. Smaller CPU tiles expose enough independent work to the
available threads, reducing this check from 598.19 to 230.12 seconds on the
current machine. Installation, relocation and missing-resource refusal all pass.

CPU session cancellation now reaches each pixel, camera sample and trace
interval. The outgoing horizon-chart worker receives the same cancellation
predicate. An interrupted trace discards its private physical result and
rejects its pending device interval before releasing trace ownership. Regressions
cancel an accepted private interval in flat and Schwarzschild charts, reuse the
tracer successfully, and interrupt an active render before any tile publication.
The active-session check passes in 0.401 seconds including scene initialization.
The headless CPU viewer completes and publishes its full 64x64 Schwarzschild
frame in 536.487 seconds. It uses smaller tiles to distribute the unchanged
physical samples among workers. Its completion deadline is separate from the
interval-level cancellation regressions.

The recovered checkpoint compiled every application, shader (with SPIR-V validation),
and test target with the Linux GCC preset and warnings treated as errors.
Repository structure, operating-model validation and negative controls, build-policy
negative controls, and generated CTest labels passed. The baseline Mandatory run
passed its first 678 tests; it was stopped during the Vulkan parity tests before
resuming implementation. This is not a complete Mandatory pass. CI formatting
was then applied to 22 first-party files. A later Mandatory run at `f76db09`
passed its first 687 tests, including physical Radeon camera and metric probes,
before being stopped to resume development. Logs remain in `.git/closeout/`.

The configured volume's `mandatory_gate.json` is the authority for the final
revision's complete zero-failure, zero-skip test estate. It binds the exact
source identity and product/test artifacts. Historical and focused test logs
must not be substituted for that receipt or for independent release evidence.
