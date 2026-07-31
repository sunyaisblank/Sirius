# Sirius adversarial operational review

Date: 2026-07-31 (Australia/Sydney)

This document is the current ground truth for Sirius. `SPECIFICATION.md` remains
the target-state mandate; `ENGAGEMENT_REPORT.md` is historical evidence, not an
authority for current counts or closure.

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

## 3. Confirmed defects and dispositions

| Finding | Ground-truth witness | Disposition |
|---|---|---|
| Installed binary omitted SPIR-V and viewer shaders and embedded a build kernel path | Staged install contained only the binary and starfield | Fixed: one resource locator, staged build tree, complete install tree, install-time invariant |
| Source/build tree or an escaping symlink could mask an incomplete volume | Runtime searched working-directory parents and accepted canonical files outside the selected root | Fixed: only executable-volume roots are searched; absolute/traversing names and escaping symlinks decline; relocation removes an asset beside a valid working-directory decoy |
| Vulkan renderer could compile without Slang kernels | Backend target existed before `slangc` detection | Fixed: live Vulkan render is compiled only when backend and kernels both exist |
| CPU-only build hosts still required OpenGL development files | app configuration used unconditional `find_package(OpenGL REQUIRED)` | Fixed: a viewer-disabled profile compiles and installs the complete CPU/Vulkan CLI, reports the absent capability, and makes `view` fail closed |
| Normal builds did not run Mandatory tests | Gate default was off and had a skip option | Fixed: gate defaults on; test omission conflicts with an operational build |
| A reused ordinary-preset cache could retain the Mandatory gate as off, with the same exposure for warning/contract policy | shared preset relied on option defaults rather than setting policy | Fixed: every supported preset pins Mandatory tests on, warnings as errors, and enforce-mode contracts; explicit overrides remain cache-visible non-operational choices |
| Release contracts only observed violations and continued | default followed `NDEBUG` rather than the operational policy | Fixed: first-party targets compile with enforce mode in every build type by default; weaker modes are explicit cache-visible non-operational choices |
| Unknown test suites silently became non-gating Correctness suites | Label generator default branch | Fixed: unknown suites are hard errors |
| Parameterised/typed GoogleTests could escape the static name-to-label authority | `TEST_P` discovery creates names the generator never parsed | Fixed: unsupported dynamic registration macros are governance errors, with a build-gated negative control |
| Explicit missing/malformed config fell back to defaults | Loader caught errors and continued | Fixed: explicit files, JSON shape, validation, and set environment values fail closed |
| JSON typos and trailing/non-finite numbers were accepted or ignored | Unknown keys merged; `stod("1junk")` consumed a prefix | Fixed with strict shape and complete finite parsers |
| Unknown top-level command returned help with exit 0 | unrecognised command collapsed to an empty token | Fixed and Mandatory-gated |
| `optix`/`cuda` backend names silently selected Vulkan | GPU-name alias predicate | Fixed: retired names fail validation |
| Non-square sample counts rendered fewer samples | floor(sqrt(spp)) squared | Fixed by the exact-count `ForEachPixelSample` authority |
| Vulkan rendered one sample while reporting arbitrary requested SPP | no SPP kernel parameter | Fixed: every requested sample is dispatched and accumulated, including non-square counts |
| Hardware attestation, readiness, and rendering could name different Vulkan devices | renderer created device zero while the runbook matched only a device name | Fixed: strict zero-based selection is shared by inventory, readiness, dispatch, memory governance, and attestation |
| Vulkan memory governance used a device-local heap while buffers require host-visible coherent memory | budget and allocation domains differed | Fixed: enumeration reports the renderable heap, memory types are ranked by compatible-heap capacity, and the governor uses that exact domain |
| Vulkan ignored camera aberration, finite aperture, point stars, temperature-model choice, and Doppler suppression | fields absent from kernel parameter mapping | Fixed: all controls reach the Slang kernel and have live effect witnesses |
| Vulkan had no volumetric disk, turbulence, or corona | explicit capability rejection | Fixed for up to 128 midpoint samples per segment; larger requests decline before dispatch |
| Observer azimuth and CPU disk controls were accepted but ignored/hard-coded | session used `phi=0`; CPU colour used 30000 K | Fixed: azimuth, temperature scale, and model reach the live CPU path; Vulkan azimuth basis is wired |
| Typed camera API advertised unimplemented lens modes and silently made them pinhole; invalid worldlines were clamped | `Orthographic`/`Panoramic` enum values and the aberration clamp | Fixed: only represented lens types exist, malformed enum values contract-fail, and non-finite/superluminal internal worldlines contract-fail after operator validation |
| Session conversion inferred pinhole, digital film, or Novikov–Thorne from unknown values; Vulkan inferred pinhole/Novikov–Thorne from malformed typed values | conversion ternaries/default branch and Vulkan parameter packing | Fixed: disk temperature is a typed post-boundary enum; conversion explicitly parses every lens, film, and temperature value; CPU and Vulkan reject malformed typed values |
| Typed session callers could reach allocation or Vulkan capability selection with invalid dimensions, output types, metric identity, dependent features, or malformed enums | validation existed only in the JSON/CLI layer and output fell through to PPM | Fixed: one `SessionConfigIssue` authority validates the complete typed boundary before allocation, backend selection, or dispatch |
| Session `Start()` was synchronous and cancellation could publish partial tiles/output while progress and callbacks raced | work ran on the caller; shared counters/callbacks were unsynchronised; Vulkan committed bands directly | Fixed: asynchronous terminal lifecycle, cancellation checkpoints, transactional commits, mutex-safe snapshots/callbacks, and cancellation-without-output evidence |
| Reinitialising a tile scheduler retained its old completion count | `completed_count_` was not reset | Fixed with reset at initialisation and a repeated-use ledger witness |
| Typed point-star configuration accepted values that its generator silently clamped | validation stopped at the JSON boundary | Fixed: the session boundary rejects every non-finite/out-of-domain field; parallax is wired and unused model fields were removed |
| Direct metric construction admitted parameter combinations outside the implemented Kerr-Schild ansatz, while unknown setter keys were ignored | rotating/charged nonzero-Lambda instances and misspelled parameter names could be constructed below the config layer | Fixed: `MetricParameterIssue` is shared by config, session, and Vulkan; core constructors and setters contract-fail on unrepresented combinations and unknown keys |
| The metric-signature check inspected diagonal signs and accepted `g00` near zero | off-diagonal Lorentzian and degenerate synthetic tensors exposed the false proxy | Fixed: a deterministic symmetric Jacobi eigensolve verifies exact inertia (one negative, three positive, no zero eigenvalues), with mixed-chart fixtures corrected |
| Metrics without Page–Thorne emission either inherited the Kerr approximation or accepted `diskEnabled` and silently rendered no disk | metric lensing and emission capability were not separated; `NotApplicable` bypassed rejection | Fixed: exact live disk support is restricted to Schwarzschild and Kerr; every other metric requires `--no-disk` at config, typed-session, and Vulkan boundaries |
| Rendered polarisation controls could suggest Stokes transport while the live output only formed a local emission angle | render flags and colour mode were disconnected from transported Stokes state | Fixed on CPU: one operator `colorMode` enables a transported observer screen basis, Page–Thorne emitter frame, thermal Thomson Stokes emission, crossing accumulation, and film visualisation; core/live/oracle Walker–Penrose gates and an end-to-end render witness are Mandatory. Vulkan, volume, and temporal-blur combinations decline |
| Thin-disk shading applied relativistic beaming more than once | the caller formed `(g T)^4`, `ApplyColorMode` applied `g^4`, and the caller multiplied by `g^4` again | Fixed: emitted `T^4` and observed `g^4 T^4` are distinct values, TrueColor owns exactly one `g^4`, and a numeric unit witness pins the formula |
| Kerr disk calculations used Kerr–Schild cylindrical radius as Boyer–Lindquist orbital radius | at the equator `sqrt(x^2+y^2)=sqrt(r^2+a^2)`, so ISCO clipping, Page–Thorne temperature, volume density, redshift, and polarised emitter motion used the wrong radius | Fixed: the shared spheroidal-radius authority now drives CPU and Slang thin/volume disk paths; an exact coordinate witness distinguishes the two radii |
| Malformed vector indices, corona geometry, colour modes, tonemappers, trace outcomes, and film enums silently aliased or fell back | internal typed APIs retained permissive default branches after operator validation | Fixed: every malformed internal enum/index/state contracts or produces an explicit diagnostic rather than selecting another represented mode |
| Fisheye existed in the core camera factory but could not be selected through configuration | schema and session projection only named pinhole/thin lens | Fixed: Fisheye is a represented CPU mode; explicit Vulkan declines and `auto` selects CPU at the capability boundary |
| Viewer flags were parsed but the progressive session used hard-coded controls; frame/session state raced and shipped shaders were never consumed | CLI/config values did not reach the session; inline GLSL bypassed volume assets | Fixed: explicit projection, synchronised asynchronous refinement/cancellation, file-backed checked shaders, and a real headless frame witness |
| The viewer liveness test embedded local-host performance assumptions | hosted ASan/UBSan/LSan completed one refinement but teardown canceled a healthy second render at the fixed deadline | Fixed: the witness keeps a valid 64x64 Schwarzschild render, asynchronous completion, frame-size, state, and error assertions without duplicating multi-sample performance coverage; the exact sanitizer rerun passes |
| Missing starfield changed the image to grey/analytic background | CPU warning and Vulkan fallback | Fixed: the starfield is a mandatory scene input and absence declines |
| Kerr coordinate round-trip tolerated an intentionally wrong azimuth | forward transform omitted the Kerr phase while inverse subtracted it | Fixed with the exact oblate transform and 1e-12 round-trip gates |
| Horizon-capture test asserted an identity | `x + (total-x) >= total` | Fixed with a positive captured/disk-ray postcondition |
| One permanently disabled unstable-orbit test inflated registration counts | `DISABLED_PhotonSphereRadius` | Removed; the enabled instability test is the represented claim |
| `GTEST_SKIP()` could still make required-profile evidence look green | CMake maps GoogleTest skips to successful CTest skips | Fixed: the required profile force-includes a test policy that converts every `GTEST_SKIP()` to a fatal assertion; portable capability-domain skips remain explicit |
| P1 status had no Bardeen boundary or burn-in | no test referenced Bardeen screen coordinates | Fixed: CPU and Vulkan independently classify ten samples spanning the visible upper curve at 1920x1080, a/M=0.998, below one pixel; the required runtime profile repeats both three times |
| P1 CPU classifier counted any long-orbit ray as captured | cinematic `Spiraling` work bound was included in the shadow predicate | Fixed: the analytic classifier disables that policy, uses a strong-field-only fine step cap, and rejects numerical/max-step outcomes |
| Page-Thorne existed as an isolated model while the live paths used a legacy correction | rendered temperature did not equal the full model | Fixed: CPU and Slang render paths use the full Page-Thorne flux shape and cross-check against the independent Core model |
| Volumetric transfer restarted and overwrote at every RK45 segment | optical depth depended on the final traversed segment | Fixed: emission and optical depth form one recurrence across the ray; deterministic turbulence and inverse-Compton corona are live on CPU and Vulkan |
| “Beam ellipse” collapsed to a circular Gaussian; Vulkan propagated one vector | star filter accepted one sigma and kernel alpha carried one norm | Fixed: both paths propagate two vectors, extract both singular axes and output-plane orientation, and apply an anisotropic tangent-plane Gaussian |
| Beam orientation used the input right-singular vector | formula used `ab+cd` while documenting output position angle | Fixed in live and oracle geometry to the `MM^T` expression `ac+bd`, with a rotated-SVD witness |
| The oracle beam step treated a Boyer-Lindquist coordinate derivative as a covariant derivative and froze curvature outside the ray's integration tableau | the exact radial and circular Schwarzschild null-congruence solutions exposed the missing connection/stage coupling | Fixed: central ray and four covariant Jacobi columns share one RK4 tableau; radial and photon-sphere screen axes agree with independent closed forms to one part in 1e6 |
| A composed symplectic substep evaluated Boyer–Lindquist derivatives after crossing the chart's horizon | enforce-mode metric contracts terminated two conservation fixtures | Fixed: every substep checks the chart domain before its second derivative evaluation, termination never evaluates an invalid metric, and the energy/angular-momentum witnesses use long scattering rays with a minimum-step postcondition rather than vacuous plunges |
| The flat CPU/Vulkan point-catalogue parity fixture changed the metric to Minkowski but retained Kerr spin and disk state | the shared typed session boundary correctly declined both sessions | Fixed in evidence: the fixture now projects a coherent massless, spinless, diskless Minkowski scene instead of weakening product validation |
| Point-star sampling scanned 100,000 stars per escaped pixel | catalogue had no live acceleration structure | Fixed: deterministic latitude/longitude CSR index, exact exhaustive-oracle agreement, shared CPU/Vulkan candidate semantics, and bounded residency |
| Doppler toggle test attributed Kerr lensing/frame-dragging asymmetry to emitter motion | image-half ratio was the sole oracle | Fixed: the gate separately measures observed asymmetry and the isolated `(g/g_grav)^4` emitter factor; off is exactly zero in the isolated measure |
| Leak checking had no reproducible operational profile | ad-hoc ASan runs disabled leak detection around Vulkan | Fixed: GCC ASan/UBSan/LSan preset and CI job; one narrow, printed suppression names the repeatable 128-byte Vulkan-loader/driver process-lifetime allocation |
| Hardware script accepted software devices and called 4096x2864 “IMAX-class” | no physical-device precondition or attestation | Fixed: exact device matching, software rejection, readiness check, full log/hash/JSON attestation, and the specification's exact 5616x4096 frame |
| External-domain records had no enforced schema and device inventory omitted driver identity | a handwritten JSON result could call llvmpipe physical, call Dozen native Windows, or name an IMAX output without decoding it | Fixed: a Mandatory negative-control verifier rejects software devices, domain mismatch, wrong dimensions, false test state, and tampered hashes; Vulkan inventory includes driver/vendor/device/API identity; the physical runbook decodes 5616x4096 and verifies its own record |
| CPU disk motion blur existed only below the operator boundary | no JSON or CLI request could exercise it, so the polarisation/motion combination had only an internal decline | Fixed: explicit motion-blur schema and CLI controls reach the CPU session; polarised temporal integration and every Vulkan motion request decline at the operator and typed boundaries |
| Two-sheet wormholes were described only in comments | the operator could request a Morris-Thorne scene but could not state the required topology, allowing one-sheet output to be misunderstood | Fixed: `OneSheetCapture` and `TwoSheet` are explicit schema/CLI/typed identities; the former is represented, while the latter declines before render |
| GLFW repeat actions were interpreted as key releases and non-finite pointer input could poison camera state | the handler treated only action 1 as pressed and accepted NaN/Inf | Fixed: press/repeat/release semantics, finite pointer/scroll guards, camera update, and refinement restart are Mandatory-gated; native event delivery remains a separate attestation |
| P2900/P2996 and external-profile boundaries were prose-only at runtime | C++26 mode could be mistaken for native contracts/reflection and installed volumes carried no capability ledger | Fixed: compile-time language facts distinguish native features from checked-macro/explicit-schema substitutes; the 19-contract model is build-gated, installed, readiness-required, and exposed by `info capabilities` |
| Operator render scripts could report success from `tail` after a failed renderer and choose an arbitrary stale binary | pipeline status belonged to `tail`; `find ... | head -1` selected an unspecified build | Fixed: exact/explicit binary selection, status propagation, non-empty artefact checks, and failure/no-output negative controls |
| EXR metadata was constructed but discarded, finite fp32 radiance could overflow half storage, and direct writers accepted NaN/Inf | writer header and advertised color-space claims did not match the file | Fixed: fp32 channels, embedded Sirius metadata, conservative color naming, and finite checks at display, session, Vulkan, PNG, and EXR boundaries |
| Native CI passed its full test estate but could not issue the promised build attestation | CTest placed relative JUnit outputs under each `--test-dir`, while the verifier looked in the repository root | Fixed: Windows and macOS attestations consume the exact binary-tree JUnit paths; publication CI is the executable witness |
| The compiled capability authority exceeded MSVC's narrow string-literal limit | the 18,303-byte operating model was emitted as one raw literal and failed with C2026 | Fixed: CMake emits independently bounded 8,000-byte chunks; runtime validation reconstructs and semantically compares the complete authority |
| Sanitizer publication limits assumed the faster review host | both pushed streams reproduced a 120-second relocated-volume timeout and a viewer teardown cancellation even though the same evidence passed locally | Fixed without dropping evidence: the specification-minimum 128x128 relocated render and every obstruction remain, the viewer witness removes redundant sampling load, and bounded test/job liveness windows admit slow instrumented hosted CPUs |
| Historical byte-identity claims had no executable consumer | no test referenced the baseline directory | Claim withdrawn; tapes remain forensic evidence only |

## 4. Enforced operating model

`tests/operating_model.json` is the machine-readable claim ledger.
`scripts/verify-operating-model.py` proves that every required dimension and
capability contract names existing Mandatory evidence. The build also verifies
that the generated label file is current. The 22 required dimensions cover compile contracts,
non-skipping and registration-complete evidence, input/config,
install/relocation, operator-script status, runtime resources, CPU rendering,
session lifecycle/cancellation, required Vulkan dispatch, exact sampling,
device/allocation identity, near-extremal/oracle physics and burn-in, polarised
transport through film, Page-Thorne/volumetric transfer, ray bundles/filtered
stars, camera/lens/film, interactive-viewer projection, output encoding,
memory/dispatch, kernel portability, and metric/decline behaviour.

The same ledger now contains 19 explicit capability contracts. Each has one
state: `supported`, `bounded`, `fail_closed`, `substituted`, or
`attestation_required`. The installed volume includes the exact ledger,
`info capabilities` exposes it, and readiness fails when it is absent.

| Requested capability | Model disposition |
|---|---|
| Schwarzschild/Kerr thin-disk CPU polarisation | `supported` |
| Polarised volume, temporal blur, or Vulkan | `fail_closed` |
| Scalar CPU temporal disk blur | `supported` through JSON/CLI/session |
| Scalar Vulkan temporal disk blur | `fail_closed` before device dispatch |
| Disk emission outside Schwarzschild/Kerr | `fail_closed`; requires `--no-disk` |
| Morris-Thorne one-sheet dark-throat scene | `supported` on CPU/Vulkan |
| Morris-Thorne two-sheet continuation | `fail_closed` through `TwoSheet` request identity |
| Vulkan volumetric samples | `bounded` to 1..128; auto selects CPU above the bound |
| Viewer input-state logic | `supported`; native window delivery is `attestation_required` |
| P2900/P2996 | `substituted` by checked macros/explicit schemas and reported non-native |
| Physical Radeon, WSL2/Dozen, native Windows/macOS/MoltenVK, physical 5616x4096 | `attestation_required` |

The staged-volume test installs Sirius, moves the prefix, runs
`info readiness` from an unrelated directory, renders a CPU Minkowski frame,
then removes the starfield while a valid hostile-working-directory decoy is
present and requires a non-zero diagnostic. Installation checks the exact
capability-specific file set and rejects empty artefacts. Separate negative
controls prove both operator scripts propagate a renderer failure and reject an
exit-zero renderer that emits no file.

`linux-ci` is the non-skipping runtime profile. Configuration fails without the
Vulkan backend and compiled kernels; CTest fails without a ready device, a
governed multi-sample dispatch, or the repeated CPU/Vulkan P1 classifiers.
Portable source profiles may omit Vulkan, but cannot be cited as evidence for
this profile.

## 5. Current status by profile

| Profile | Status | Evidence boundary |
|---|---|---|
| Configure/compile/build, GCC 14 | READY | `linux-ci`, warnings as errors, governance, Mandatory gate, required Vulkan dispatch and P1 burn-in |
| Configure/compile/build, Clang 21 | READY | Warnings as errors, Slang portability emission, Mandatory gate |
| Relocatable CPU volume | READY | Install-time invariant plus relocation/render/failure-injection CTest |
| Viewer-disabled build/install | READY | GCC 14 warnings-as-errors compile, capability-specific install verification, readiness report, and fail-closed `view` command |
| CPU physics/render path | READY WITH CLAIM-SCOPED TOLERANCES | Full estate, exact Schwarzschild/Kerr Page-Thorne/volume path, near-extremal conservation, ten-point Bardeen gate, and live thin-disk polarisation transported through Stokes film output |
| Vulkan on current software device | READY WITH EXPLICIT BOUNDS | WSL2 llvmpipe; required dispatch, all precision binaries, multi-sampling, camera/lens/disk/volume/ellipse-filtered point stars; volume samples <=128; non-TrueColor, motion blur, and polarisation explicitly excluded |
| Interactive viewer refinement | HEADLESS READY; WINDOW RUNTIME UNEXECUTED | Strict parsing, projection, asynchronous render, synchronised frame publication, and cancellation are gated; no GLFW/OpenGL window or live input event was exercised |
| GCC ASan + UBSan + LSan | READY WITH ONE EXTERNAL SUPPRESSION | Broad Mandatory estate plus repaired/changed-surface reruns; one named 128-byte Vulkan loader/driver allocation suppressed and printed |
| Physical Radeon 780M | UNEXECUTED IN THIS REVIEW | Historical report names a run, but current environment exposes llvmpipe and no fresh physical-device evidence was produced |
| Native Windows build/runtime | BUILD ATTESTED; RUNTIME UNEXECUTED | Windows Server 2022 CI configures, builds, passes the full source-available estate, and emits a revision-bound build attestation; native Vulkan/device execution remains separate |
| macOS build/MoltenVK runtime | BUILD ATTESTED; RUNTIME UNEXECUTED | macOS 15 CI configures, builds, passes the full source-available estate, and emits a revision-bound build attestation; MoltenVK/device execution remains separate |

## 6. Remaining limitations (not hidden as “green”)

- P2 now has the specification's exact radial and circular Schwarzschild
  congruence pair at 1e-6, in addition to literal dual-vector ellipses,
  Riemann/orientation gates, live CPU/Vulkan behaviour, and the rotating-star
  flicker witness. This closes the earlier P2 obstruction; it does not supply
  evidence for unexecuted physical platforms.
- E2 is represented on the CPU Schwarzschild/Kerr thin-disk path. The live
  Kerr–Schild transport conserves the Walker–Penrose constant and agrees with
  the independent Boyer–Lindquist oracle; each disk crossing emits a physical
  thermal Thomson Stokes state in the transported observer basis, and
  `--color-mode Polarisation` reaches the film buffer. `motionBlur` is now an
  operator-visible scalar CPU capability. Polarised volumetric transfer,
  temporal disk blur, and the Vulkan kernel remain unrepresented and decline
  at configuration and typed-session boundaries rather than approximating the
  thin-disk CPU result.
- Full Page–Thorne disk emission is represented only for Schwarzschild and Kerr.
  Every other metric requires `--no-disk`; charged/cosmological disk models are
  absent, while Minkowski, de Sitter, Morris–Thorne, and Alcubierre have no
  accretion-disk semantics. Morris–Thorne `OneSheetCapture` renders one
  asymptotic sheet with a dark throat; an explicit `TwoSheet` request declines.
- The exact 5616x4096 memory plan, catalogue/index residency, and sublinear
  candidate query are Mandatory-gated, and the physical-hardware runbook now
  renders that exact frame. A fresh end-to-end 23-million-pixel filtered-star
  render was not executed on the current software device. P3/P5 IMAX runtime
  closure remains a physical-profile item.
- Vulkan volumetric transfer deliberately caps `volumetric.samples` at 128 per
  geodesic segment to protect the first dispatch from an unbounded watchdog
  exposure. The CPU accepts the schema maximum of 4096; explicit Vulkan requests
  above 128 decline and `auto` logs the boundary before selecting CPU.
- `backend.enableDenoiser` and nonzero legacy `backend.cudaDevice` are retained
  as compatibility sentinels. They decline instead of implying unavailable
  denoiser/CUDA capabilities; neither is a DNGR parity criterion.
- The viewer's strict parsing, projection, refinement render, frame publication,
  cancellation, and press/repeat/release/pointer/scroll state transitions are
  gated. The review did not open a GLFW/OpenGL window or receive native input;
  that exact domain requires a `viewer-native-window-input` attestation.
- P2900 contracts and P2996 reflection remain absent from the measured
  compilers. Enforced contract macros and explicit serializers/parameter
  bindings are the active, gated substitutes. `info system` reports both
  native feature states as false.
- Physical 780M, WSL2/Dozen, native Windows Vulkan, and macOS/MoltenVK require
  fresh attestations in their own domains. The verifier rejects llvmpipe for
  physical domains, Dozen for native Windows, and non-MoltenVK macOS runtime
  evidence. Native CI now emits build-only attestations after its test gate.
- Historical reference images have no current byte-identity test. Any new
  identity claim requires a checked-in manifest and an executable comparator.

## 7. Validation record

The counts below are refreshed only after the final source state completes the
corresponding profile. They are not projected from registration counts.

- GCC 14 strict required-Vulkan profile and repeated CPU/Vulkan P1 burn-in:
  889/889 passed; 678 Mandatory, 152 Operational, 47 Performance, and 27
  Stability tests passed. The isolated P1 burn-in ran for 14.42 seconds with
  required Vulkan and volume evidence.
- GCC 14 Release portable profile: 886/886 passed; 675 Mandatory, 149
  Operational, 47 Performance, and 26 Stability tests passed.
- Clang 21 Release with warnings as errors: 886/886 passed; 675 Mandatory, 149
  Operational, 47 Performance, and 26 Stability tests passed.
- GCC 14 ASan/UBSan/LSan: the complete 675-test Mandatory estate passed in
  1400.96 seconds. The run used the one documented, printed external
  Vulkan-loader/driver suppression.
- Slang `trace.slang`: SPIR-V validated for fp32, compensated fp32, and fp64;
  CUDA and Metal emission gates pass.
- Install relocation, hostile-working-directory initialisation, and resource
  failure injection passed in the estates above. A viewer-disabled GCC build
  compiled and installed, reported the capability absent, remained CPU/Vulkan
  ready, and rejected `view`. Reduced real operator runs produced 2/2 smoke and
  13/13 demonstration artefacts. The installed operating model reports 22
  review dimensions and 19 capability contracts. Attestation negative controls,
  shell syntax, JSON governance, label freshness, formatting, and diff
  whitespace were checked separately.

## 8. Scope note

The system in scope is the complete Sirius repository and its build, install,
CLI, operator-script, CPU, Vulkan, viewer, physics, and output surfaces. External
drivers, physical devices, and native operating systems remain distinct
attestation domains and are never inferred from source or Lavapipe evidence.
