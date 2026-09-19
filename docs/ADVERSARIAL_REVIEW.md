# Sirius adversarial operational review

Corrective review: 2026-09-07; consolidated 2026-09-19 (Australia/Sydney)

This is a review record, not a qualification receipt. The September corrective
work below supersedes conflicting historical dispositions. The retained earlier
ledger describes its own source snapshots; its counts, timings, “Fixed” entries
and physical outputs do not establish the state of a later revision.
`SPECIFICATION.md` and `tests/operating_model.json` define the obligations; the
clean source, executed artifacts and independently verified same-revision
receipts establish completion. The later local renderer result and remaining
performance limits are recorded in `SOURCE_CLOSEOUT.md`. Detailed historical
defects, session logs and the retired engagement report are preserved at
commit `9e20150`; use `git show 9e20150:docs/ADVERSARIAL_REVIEW.md` for that ledger.

## September corrective work

At this review checkpoint, focused diagnostics establish only the boundaries
they exercise. They are not complete Mandatory qualification, native-platform
admission, full-resolution physical evidence or installed-release readiness.
Later completion must be established by the exact revision's verified receipts.

| Finding | Current correction and evidence boundary |
|---|---|
| F03: a gate hashed candidate products only after testing | Both full and native build gates snapshot tested/product bytes before CTest, then reject differences in those bytes, source identity or full registration before issuing a receipt. Focused mutation controls pass, including restored obstruction fixtures. These boundary comparisons do not detect arbitrary transient tampering restored between observations. |
| F04: a generic pointer event could substitute for required input delivery | Viewer production and independent admission now require cursor and scroll observations alongside keyboard callbacks and a newly published Vulkan frame. Negative controls reject incomplete callback transcripts. Actual host delivery still requires a new native viewer campaign. |
| F06: viewer evidence was permanently reported as pending | Readiness distinguishes admitted exact-revision viewer evidence from current-host window/input availability. Admission does not assert that a window has been created on the querying host. Final installed-release reporting remains an execution obligation. |
| F08: a rejected CPU step could invent horizon capture by backward extrapolation and radial snapping | The CPU black-hole path now transforms the original observer, tangent and screens into a regular outgoing Kerr-Schild chart and localises capture on accepted finite segments. Exterior results return to the public chart; horizon results identify the outgoing chart. Focused outward-ray and shadow diagnostics support the correction; chart, coupled-transport and complete scientific acceptance remain subject to final verification. |
| F10: subtracting accumulated float affine times could erase or distort a short accepted step | Coupled CPU consumers use the accepted step's actual interval, with trace-local double accumulation and published affine length. Horizon and other terminal clips shorten that same interval. A live radial-ray witness checks the analytic horizon interval and Jacobi size after steps smaller than the accumulated float spacing; restoring the old subtraction fails that witness. |
| F09: strict fp64-specific skips rejected supported precision-limited routes | The Vulkan kernel boundary rejects unsupported Float64 modules and malformed instruction framing. Existing precision tests execute supported fp32/compensated science and require actual fp64 refusal when unsupported; capable devices retain their fp64 numerical comparisons. The render refusal checks its precise diagnostic. No-device boundary execution and compilation checks pass; native precision-limited execution and Radeon fp64 qualification remain outstanding. |
| F12: parallel CTest could overlap independent Vulkan consumers | Backend and render discovery now share the `sirius_vulkan_device` resource lock, including portability, render-capable CLI/viewer and end-to-end cases. Existing operational Vulkan and burn-in wrappers retain `RUN_SERIAL`; base, core, oracle and non-rendering app discovery remain parallel. This establishes device ownership within one CTest invocation, not exclusion against a separate CTest process or external GPU workload. Generated inventory checks establish registration coverage; physical execution remains a separate obligation. |
| F01/F02: caps could be bypassed and the cold controller could discard useful work without reducing latency | Independent width, row and area caps remain active with adaptation disabled. Heavy fp32 uses a measured 64×4 maximum and 750 ms soft target; other profiles retain their original target and conservative precision limits. Discrete growth avoids the cold one-row trap. The full governed-scene 128×128 diagnostic remains byte-identical, completing in 288 submissions with a 778 ms observed maximum and no safety fallback. A 512-pixel probe exceeded the stop threshold and was rejected. Full-resolution runtime and stability remain separate qualification obligations; a work cap cannot guarantee wall time. |

The final gate must run on the accepted clean revision. Eight exact-revision
domains, both governed physical image sizes, and explicit post-gate installation,
relocation and render/view initialisation remain separate requirements. The
strict install test's pre-gate rejection branch does not prove the later
installed workflow. Any further source or documentation correction precedes
freeze and requires new revision-bound evidence afterward.

## 1. Method

The review adapts the workflow in the appendix of
<https://arxiv.org/pdf/2607.13335>: define the exact model and success
conditions first; state what does not count; maintain distinct construction and
obstruction ledgers; demand concrete witnesses; audit every candidate for
fidelity, precision, consistency, extension, normalisation, bounds, domain,
quantifiers, imported assumptions, and circularity; then freeze the candidate
and repeat the audit.

For Sirius, the “model” is the complete operating path:

`source -> configure -> compile -> build gate -> install -> relocate volume ->
initialise -> select capability -> render -> write output -> diagnose`

The minimum operational profile is a relocatable CPU render. Vulkan, physical
GPU, and native-platform profiles are reported separately and cannot borrow
evidence from the CPU profile.

Success requires:

1. Every operator input is either represented exactly or rejected before work.
2. Every runtime resource is in the installed volume and resolves without the
   source/build tree.
3. Every required operating dimension has named, build-failing evidence.
4. Physics claims identify an independent oracle and a numerical tolerance.
5. A capability unavailable in the active environment is `UNEXECUTED`, never
   inferred green.

### Results that do not count

- A successful build that did not execute the Mandatory label.
- A registered, skipped, disabled, tautological, or assertion-free test.
- A source file, reference image, or historical log with no current consumer.
- “Compiled in” as evidence that a device, kernel, or resource is usable.
- A fallback that changes the requested spacetime, samples, features, or
  background.
- CPU evidence used to close a Vulkan claim, or Lavapipe evidence used to close
  a physical/native-driver claim.
- A test count used as a proxy for coverage.

## 2. Adversarial audit matrix

Each review pass asks:

| Audit | Sirius question |
|---|---|
| Fidelity | Is the requested metric, feature, sampling count, and asset the one rendered? |
| Precision | Are numeric parses complete and finite, and are physics tolerances explicit? |
| Adaptivity | Can backend auto-selection silently cross a capability boundary? |
| Fixed consistency | Do schema, CLI, session, kernel, docs, and tests name the same behaviour? |
| Global extension | Does a build-tree success remain valid after install and relocation? |
| Normalisation | Are camera frames, null rays, transfer encodes, and units applied exactly once? |
| Bounds | Are memory, dispatch duration, dimensions, samples, and output types bounded? |
| Output | Is success withheld when an output or mandatory resource is incomplete? |
| Domain | Are CPU, software Vulkan, physical Vulkan, Windows, WSL2, and macOS distinguished? |
| Quantifiers | Does “all backends/platforms” have evidence for every named member? |
| Imported results | Are papers, historical logs, and reference tapes independently connected to a live gate? |
| Non-circularity | Is the oracle independent of the implementation it judges? |

## 3. Operating-model evidence boundaries

`tests/operating_model.json` is the machine-readable claim ledger.
`scripts/verify-operating-model.py` proves that all ten P1–P6/E1–E4 acceptance
criteria, every required dimension, and every capability contract name existing
Mandatory evidence. The build also verifies that the generated label file and
repository structure are current. Seven P/E criteria are build-gated; P3, P5,
and E3 remain `attestation_required` because software evidence cannot establish
physical IMAX/780M operation. The 24 required dimensions cover attestation
admission/release alignment, compile contracts,
non-skipping and registration-complete evidence, input/config,
install/relocation, operator-script status, runtime resources, CPU rendering,
session lifecycle/cancellation, required Vulkan dispatch, exact sampling,
device/allocation identity, near-extremal/oracle physics and burn-in,
natural-scale trace-domain geometry, polarised
transport through film, Page-Thorne/volumetric transfer, ray bundles/filtered
stars, camera/lens/film, interactive-viewer projection, output encoding,
memory/dispatch, kernel portability, and metric/decline behaviour.

The same ledger now contains 30 explicit capability contracts. Each has one
state: `supported`, `bounded`, `fail_closed`, `substituted`, or
`attestation_required`. The installed volume includes the exact ledger,
`info capabilities` exposes it, and readiness fails when it is absent.

| Requested capability | Model disposition |
|---|---|
| Revision-bound release alignment | `supported`; complete clean-revision admission is enforced at configure/build/compile/runtime boundaries |
| Schwarzschild/Kerr thin-disk CPU polarisation | `supported` |
| Polarised volume, temporal blur, or Vulkan | `fail_closed` |
| Scalar CPU temporal disk blur | `fail_closed`; requires time-dependent emissivity and a covariant shutter integral |
| Scalar Vulkan temporal disk blur | `fail_closed` before device dispatch |
| Inverse-Compton corona and narrowband line transfer | `fail_closed`; frequency-dependent source/opacity models are absent |
| Grey volumetric disk | `bounded`; stationary Gaussian, vertically isothermal phenomenology with declared optical-depth law and optional procedural modulation |
| Cartesian Kerr-Schild vacuum/electrovac family | `bounded`; live vacuum, Einstein-Maxwell, and source-free Maxwell equations are independently gated, while rotating/charged nonzero-Lambda sectors decline |
| Spherical de Sitter/Kottler sector | `bounded`; exact positive-Lambda metric and both Kottler roots, with pure de Sitter horizonless for capture, Schwarzschild-de Sitter restricted below Nariai, observer at or below `0.99*r_c`, and directional boundary no later than `r_c` |
| Disk emission outside Schwarzschild/Kerr | `fail_closed`; requires `--no-disk` |
| Morris-Thorne one-sheet dark-throat scene | `supported` on CPU/Vulkan |
| Morris-Thorne two-sheet continuation | `supported` on CPU/Vulkan with an asymptotically normalised lapse, independently checked exact curvature/NEC violation, inversion-matched opposite boundary, and symmetric catalogue |
| Vulkan volumetric samples | `bounded` to 1..128; auto selects CPU above the bound |
| Viewer input-state logic | `supported`; native window delivery is `attestation_required` |
| P2900/P2996 | `substituted` by checked macros/explicit schemas and reported non-native |
| Physical Radeon, WSL2/Dozen, native Windows/macOS/MoltenVK, physical 5616x4096 | `attestation_required` |

The staged-volume test installs Sirius, moves the prefix, runs
`info readiness` from an unrelated directory, verifies the exact ten-criterion
P/E summary, mutates P1 to prove semantic rejection, renders a CPU Minkowski frame,
then removes the starfield while a valid hostile-working-directory decoy is
present and requires a non-zero diagnostic. Installation checks the exact
capability-specific file set and rejects empty artefacts. Separate negative
controls prove the cross-platform render workflow propagates a renderer failure
and rejects an exit-zero renderer that emits no file.

`linux-ci` is the non-skipping runtime profile. Configuration fails without the
Vulkan backend and compiled kernels; CTest fails without a ready device, a
governed multi-sample dispatch, or the repeated CPU/Vulkan P1 classifiers.
Portable source profiles may omit Vulkan, but cannot be cited as evidence for
this profile.

## 4. Evidence semantics by profile

Live attestation status is deliberately not frozen into this source file. An
attestation names the exact source revision that it verifies, so committing a
new status count would create a different revision and immediately invalidate
that count. The deterministic alignment receipt generated from independently
verified, same-revision bundles is the only authority for the live admitted and
pending partition. The table records durable implementation and evidence
boundaries instead.

| Profile | Status | Evidence boundary |
|---|---|---|
| Revision-bound release alignment | EXTERNAL LEDGER AUTHORITATIVE | A source-only qualification configure intentionally admits no external domains. Development artifacts remain inadmissible; release configure reports every domain absent from the supplied same-revision ledger, and packaging/initialisation remain fail-closed until all eight domains have verified evidence. |
| Pull-request/integration boundary | GOVERNED NON-RENDER PATH PRESENT | Linux, Windows, and macOS compile the complete strict topology and execute exactly nine authority controls without creating a Mandatory receipt. Pull requests cannot publish evidence; explicit dispatch may issue only the precisely scoped Windows/macOS compilation domains. |
| Configure/compile/build, GCC 14 | STRICT NON-RENDER PATH PRESENT | Qualification binds all six test executables, the candidate and nine live products under warnings as errors, then runs the exact nine authority controls without issuing a Mandatory receipt. Execution outcomes belong to the workflow record for the tested revision. |
| Configure/compile/build, Clang 21 | STRICT NON-RENDER PATH PRESENT | Qualification also emits and validates every Slang kernel, binds all six test executables, the candidate and nine live products under warnings as errors, and runs the exact nine authority controls without issuing a Mandatory receipt. Execution outcomes belong to the workflow record for the tested revision. |
| Relocatable CPU volume | LOCAL VALIDATION AT `cdcc254` | Complete installation, relocation and CPU image checks passed in the 934-test local selection; this is not a later-revision or release receipt. |
| Viewer-disabled build/install | PRE-ALIGNMENT FULL-PROFILE SNAPSHOT | Capability-specific install verification, readiness, and fail-closed `view` were previously exercised; no current-delta full profile is claimed |
| CPU physics/render path | LOCAL VALIDATION AT `cdcc254` | The complete 934-test local selection included scientific calculation, shadow, Doppler, polarisation, output writers and the CPU viewer. See `SOURCE_CLOSEOUT.md` for exact scope and exclusions. |
| Vulkan on WSL2 software and physical devices | PRE-ALIGNMENT FULL-PROFILE SNAPSHOT | Prior llvmpipe/Radeon-Dozen evidence covered dispatch and bounded scene semantics, but historical records are not admitted for a later source revision. |
| Interactive viewer refinement | PRE-ALIGNMENT PREFLIGHT | The prior source opened a GLFW/OpenGL XWayland window, published a Radeon Vulkan frame, and received host-delivered input; native-viewer admission is asserted only by a verified same-revision external record. |
| GCC ASan + UBSan + LSan | PRE-ALIGNMENT FULL-PROFILE SNAPSHOT | The historical 904-test sanitizer result is not promoted to evidence for a later source revision. |
| Physical Radeon 780M | PRE-ALIGNMENT WORKTREE SNAPSHOT | The historical 907-test and 697-test Radeon/Dozen runs plus both exact sparse-star frames are revision-specific and cannot be admitted for later source. |
| Native Windows build/runtime | BUILD PRODUCER PRESENT; LIVE RECORD EXTERNAL | Explicit Windows dispatch can emit a clean-Git compilation record binding the full registration, exact non-render authority estate, all executables/products, and receipt. Native Vulkan remains a separate full-estate physical-host domain and rejects Dozen. The alignment receipt, not this table, states whether either record exists for a revision. |
| macOS build/MoltenVK runtime | BUILD PRODUCER PRESENT; LIVE RECORD EXTERNAL | Explicit macOS dispatch can emit a clean-Git compilation record binding the full registration, exact non-render authority estate, all executables/products, and receipt. MoltenVK remains a separate full-estate physical-host domain. The alignment receipt, not this table, states whether either record exists for a revision. |

## 5. Previously recorded limitations

These boundaries remain distinct from release qualification. Historical
measurements and their detailed diagnoses remain in Git at `9e20150`; reassess
implementation and evidence against the selected release revision.

- P2 now has the specification's exact radial and circular Schwarzschild
  congruence pair at 1e-6, in addition to literal dual-vector ellipses,
  Riemann/orientation gates, live CPU/Vulkan behaviour, and the rotating-star
  flicker witness. This closes the earlier P2 obstruction; it does not supply
  evidence for unexecuted native operating systems.
- E2 is represented on the CPU Schwarzschild/Kerr thin-disk path. The live
  Kerr–Schild transport conserves the Walker–Penrose constant and agrees with
  the independent Boyer–Lindquist oracle; each disk crossing applies the
  flux-normalised Chandrasekhar-Sobolev semi-infinite pure electron-scattering
  atmosphere in the transported observer screen, and `--color-mode
  Polarisation` reaches the film buffer. Absorption, finite optical depth,
  returning radiation, magnetic/Faraday effects, polarised volumetric transfer,
  temporal disk blur, and Vulkan polarisation remain unrepresented and decline
  at configuration and typed-session boundaries.
- Full Page–Thorne disk emission is represented only for Schwarzschild and Kerr.
  Every other metric requires `--no-disk`; charged/cosmological disk models are
  absent, while Minkowski, de Sitter, Morris–Thorne, and Alcubierre have no
  accretion-disk semantics. Morris–Thorne `OneSheetCapture` renders one
  asymptotic sheet with a dark throat; `TwoSheet` traverses to the opposite
  asymptotic cutoff and uses the same catalogue through the inversion-related
  sky frame. Distinct content assigned to the second universe is absent.
- Alcubierre is an exact bounded kinematic metric family, not a matter-field,
  formation, stability, causality, or propulsion model. Its independently
  reconstructed Einstein tensor now proves the required negative Eulerian
  energy on the off-axis wall. That closes the prior geometry-versus-source
  ambiguity without claiming that a physical stress-energy construction is
  represented.
- P3/P5 physical qualification requires both 1920×1080 and 5616×4096 images
  of the governed moving ThinLens scene with beam-filtered 100,000-star
  sampling under the 2048 MiB cap. Admission requires sparse disk-free
  morphology, hashed typed request and actual-source records proving host
  catalogue/beam use, and the same transcript's measured retained stages,
  region coverage, device, budget, completion and output identity.
  Earlier full-resolution frames predate the physical pupil and sampling
  corrections and cannot satisfy these requirements. The complete small
  CPU/device comparison recorded in `SOURCE_CLOSEOUT.md` is also insufficient.
- Earlier retained-device timings exposed a substantial throughput limit:
  eight Kerr pixels at `bf30294` took more than three hours on Radeon/Dozen,
  and the larger moving ThinLens scene at `c0634c7` did not complete a pixel
  before interruption. Current discovery shares canonical 32×32 regions and
  batches independent probes; the 192×128, three-sample CPU workload completed
  in 386.29 seconds at `9661999`. These results and their exact scenes are in
  `SOURCE_CLOSEOUT.md`. Current GPU workload throughput remains unmeasured;
  correctness witnesses and component arithmetic timings do not establish it.
- Vulkan volumetric transfer deliberately caps `volumetric.samples` at 128 per
  geodesic segment to protect the first dispatch from an unbounded watchdog
  exposure. The CPU accepts the schema maximum of 4096; explicit Vulkan requests
  above 128 decline and `auto` logs the boundary before selecting CPU.
- `backend.enableDenoiser` and nonzero legacy `backend.cudaDevice` are retained
  as compatibility sentinels. They decline instead of implying unavailable
  denoiser/CUDA capabilities; neither is a DNGR parity criterion.
- The viewer's strict parsing, projection, refinement render, frame publication,
  cancellation, and press/repeat/release/pointer/scroll state transitions are
  gated. An earlier source snapshot opened the GLFW/OpenGL window and
  received host-delivered callbacks on WSLg/XWayland. That preflight is
  not promoted into a revision attestation: the current record contract also
  requires a published progressive Vulkan frame on the selected physical
  Radeon plus hashed readiness, inventory, and non-skipping JUnit evidence. It
  does not claim native-window coverage on unexecuted platforms.
- P2900 contracts and P2996 reflection remain absent from the measured
  compilers. Enforced contract macros and explicit serializers/parameter
  bindings are the active, gated substitutes. `info system` reports both
  native feature states as false.
- Native Windows and macOS build outcomes, native Windows Vulkan, and
  macOS/MoltenVK are revision-bound external facts and are not asserted by
  source prose. The verifier rejects llvmpipe for physical domains, Dozen for
  native Windows, and non-MoltenVK macOS runtime evidence. Native CI can emit
  build-only attestations after its test gate; the shared physical runtime
  producer is available but cannot substitute for execution on those hosts.
- Historical reference images have no current byte-identity test. Any new
  identity claim requires a checked-in manifest and an executable comparator.

## 6. Scope note

The system in scope is the complete Sirius repository and its build, install,
CLI, operator-script, CPU, Vulkan, viewer, physics, and output surfaces. External
drivers, physical devices, and native operating systems remain distinct
attestation domains and are never inferred from source or Lavapipe evidence.
