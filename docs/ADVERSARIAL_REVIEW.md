# Sirius adversarial operational review

Date: 2026-08-31 (Australia/Sydney)

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
| Vulkan had no volumetric disk, turbulence, or corona | explicit capability rejection | The declared phenomenological grey volume and optional procedural density modulation are represented for up to 128 midpoint samples per segment; larger volume requests and every corona request decline before dispatch |
| Observer azimuth and CPU disk controls were accepted but ignored/hard-coded | session used `phi=0`; CPU colour used 30000 K | Fixed: azimuth, temperature scale, and model reach the live CPU path; Vulkan azimuth basis is wired |
| Typed camera API advertised unimplemented lens modes and silently made them pinhole; invalid worldlines were clamped | `Orthographic`/`Panoramic` enum values and the aberration clamp | Fixed: only represented lens types exist, malformed enum values contract-fail, and non-finite/superluminal internal worldlines contract-fail after operator validation |
| Session conversion inferred pinhole, digital film, or Novikov–Thorne from unknown values; Vulkan inferred pinhole/Novikov–Thorne from malformed typed values | conversion ternaries/default branch and Vulkan parameter packing | Fixed: disk temperature is a typed post-boundary enum; conversion explicitly parses every lens, film, and temperature value; CPU and Vulkan reject malformed typed values |
| Typed session callers could reach allocation or Vulkan capability selection with invalid dimensions, output types, metric identity, dependent features, or malformed enums | validation existed only in the JSON/CLI layer and output fell through to PPM | Fixed: one `SessionConfigIssue` authority validates the complete typed boundary before allocation, backend selection, or dispatch |
| Lens, disk, bloom, volumetric, motion, output, backend-specific CPU scheduling, and film parameters could be accepted while their consuming feature was different or disabled | parameter validation checked numeric ranges but not ownership, and several CLI parameter flags never selected their consumer | Fixed: shared camera/disk/feature defaults define the only neutral inactive state; parameter flags select their owner, absent film overrides inherit their selected preset, projection names are parsed once, and config plus typed-session boundaries reject every non-default latent control before work. Vulkan also rejects CPU tile/thread/parallel controls instead of accepting work it does not schedule |
| Session `Start()` was synchronous and cancellation could publish partial tiles/output while progress and callbacks raced | work ran on the caller; shared counters/callbacks were unsynchronised; Vulkan committed bands directly | Fixed: asynchronous terminal lifecycle, cancellation checkpoints, transactional commits, mutex-safe snapshots/callbacks, and cancellation-without-output evidence |
| Reinitialising a tile scheduler retained its old completion count | `completed_count_` was not reset | Fixed with reset at initialisation and a repeated-use ledger witness |
| Typed point-star configuration accepted a broader starfield model, including controls that the point renderer did not consume, while its generator silently clamped malformed values | validation stopped at the JSON boundary and the session type overstated the request | Fixed: `PointStarfieldConfig` contains only the five consumed catalogue controls; the session boundary rejects every non-finite/out-of-domain value and both backends expand the admitted request through one deterministic generator authority |
| Core metric, camera, disk, turbulence, integrator, oracle, post-process, and film setters rewrote malformed typed input or retained phantom controls | lower-level callers could bypass operator validation and obtain a different represented request through clamps, fallbacks, ignored fields, or unused toggles | Fixed: each retained typed boundary now has a total represented-domain predicate and contract-fails malformed input without mutation; unused overloads and fields were deleted instead of documented as capabilities |
| The point-star spatial index borrowed a caller vector while advertising reusable topology | caller destruction or mutation could invalidate entries or make topology disagree with the catalogue | Fixed: the index owns a validated immutable catalogue snapshot and all indexed CPU/Vulkan consumers address that same owned storage |
| Circular fisheye masking encoded an unused zero-weight ray whose direction violated the tracer's assumptions | `CameraRay::weight` was never consumed and masked rays could reach a unit-direction-only API | Fixed: camera rays carry explicit active/inactive state, both forms have complete representation invariants, and the session writes inactive samples as black without tracing them |
| The display buffer exposed unlocked mutable storage while render, snapshot, encoding, and viewer threads shared it | raw vector access escaped the buffer mutex; a dead gamma-cache API returned another pointer after releasing the lock; malformed dimensions and tiles were silently ignored or clipped | Fixed: mutation is callback-scoped under the buffer lock, readers take stable float snapshots, the unused byte cache is deleted, shape is contract-enforced, and owned output encoders diagnose non-finite radiance before transfer encoding |
| CPU and Slang trace policy drifted in retry ceilings, turbulence coordinates, point-star fallback colour, volumetric finite checks, capture finalisation, and endpoint event order | two live implementations had accumulated independent literals and terminal assumptions | Fixed: one shared host/device attempt ceiling and capture factor feed packing; turbulence uses the same reconstructed disk coordinates; malformed radiance inputs fail closed; captured rays avoid invalid terminal frames; and both paths accumulate disk events along an accepted segment before classifying its endpoint |
| Film preset selection was overwritten by concrete default-valued fields and operator text described a stock simulation | serialisation could not distinguish an absent override from an explicit default, so non-default presets collapsed | Fixed: film finishing is explicitly bounded and non-stock; optional overrides inherit the selected preset, JSON emits only present overrides, and internal/session names no longer claim physical stock simulation |
| Vulkan session initialisation constructed unused CPU tracers, thread pools, tiles, textures, and a duplicate 100,000-star catalogue before dispatch | backend selection happened after CPU-only scene construction and CPU scheduling controls were still admitted | Fixed: typed preflight enforces backend ownership and the Vulkan branch initialises only shared lifecycle/evidence state before entering its independently governed device path |
| Direct metric construction admitted parameter combinations outside the implemented Kerr-Schild ansatz, while unknown setter keys were ignored | rotating/charged nonzero-Lambda instances and misspelled parameter names could be constructed below the config layer | Fixed: `MetricParameterIssue` is shared by config, session, and Vulkan; core constructors and setters contract-fail on unrepresented combinations and unknown keys |
| Operator output labelled the raw observer coordinate radius as a number of `M`, EXR metadata labelled geometric mass as solar masses, and metrics accepted mass/throat/warp parameters their constructors ignored | validation passed the raw radius directly to the camera while scaling only its allowed interval; metric identity and length-unit semantics diverged across schema, CLI, viewer, session, trace, and metadata boundaries | Fixed: the registry owns mass and exotic-parameter applicability; identity changes select a compatible omitted-value mass default while explicit irrelevant mass, throat/topology, and warp/bubble changes decline at config and typed-session boundaries. Distance remains the established coordinate radius `r`, mass-bearing output derives `r/M`, focus and viewer motion are labelled as geometric coordinates, and EXR declares a geometric-coordinate length unit without inventing a physical scale |
| Schwarzschild-de Sitter capture reused the asymptotically flat Kerr-Newman horizon, pure de Sitter reported a zero-radius black-hole horizon, and the renderer could launch or sample its directional sky beyond the cosmological horizon | the public cosmological sector ignored both Kottler roots, admitted negative Lambda under de Sitter names, allowed Nariai/super-Nariai scenes without a distinct black-hole exterior, and treated positive-Lambda space as though its arbitrary far sphere were asymptotically flat | Fixed: one scale-safe spherical positive-Lambda authority solves the smaller black-hole and larger cosmological roots of `1-2M/r-Lambda*r^2/3=0` for the metric, validation, and tracer; capture uses only the black-hole root, pure de Sitter never captures, the ADM/Eulerian launch radius is bounded by `0.99*r_c`, and the finite directional sky boundary is clamped to `r_c`. Its fp32 radius rounds inward and an accepted crossing segment is localised back onto the boundary before sky sampling. Named Schwarzschild-de Sitter requires `9 Lambda M^2 < 1`; negative/rotating/charged cosmological requests, Nariai/super-Nariai scenes, and observers outside the governed causal patch decline. No static-observer or cosmological-transfer claim is inferred |
| The metric-signature check inspected diagonal signs and accepted `g00` near zero | off-diagonal Lorentzian and degenerate synthetic tensors exposed the false proxy | Fixed: a deterministic symmetric Jacobi eigensolve verifies exact inertia (one negative, three positive, no zero eigenvalues), with mixed-chart fixtures corrected |
| Metrics without Page–Thorne emission either inherited the Kerr approximation or accepted `diskEnabled` and silently rendered no disk | metric lensing and emission capability were not separated; `NotApplicable` bypassed rejection | Fixed: exact live disk support is restricted to Schwarzschild and Kerr; every other metric requires `--no-disk` at config, typed-session, and Vulkan boundaries |
| Rendered polarisation controls could suggest Stokes transport while the live output only formed a local emission angle | render flags and colour mode were disconnected from transported Stokes state | Fixed on CPU: one operator `colorMode` transports an observer-screen basis on the accepted central ray, projects into the Page–Thorne circular-emitter screen, and applies the flux-normalised Chandrasekhar-Sobolev semi-infinite pure electron-scattering atmosphere; absorption, finite optical depth, returning radiation, and magnetic/Faraday effects are excluded. Vulkan, volume, and temporal-blur combinations decline |
| Thin-disk shading applied relativistic beaming more than once | the caller formed `(g T)^4`, `ApplyColorMode` applied `g^4`, and the caller multiplied by `g^4` again | Fixed: emitted `T^4` and observed `g^4 T^4` are distinct values, TrueColor owns exactly one `g^4`, and a numeric unit witness pins the formula |
| Kerr disk calculations used Kerr–Schild cylindrical radius as Boyer–Lindquist orbital radius | at the equator `sqrt(x^2+y^2)=sqrt(r^2+a^2)`, so ISCO clipping, Page–Thorne temperature, volume density, redshift, and polarised emitter motion used the wrong radius | Fixed: the shared spheroidal-radius authority now drives CPU and Slang thin/volume disk paths; an exact coordinate witness distinguishes the two radii |
| Malformed vector indices, corona geometry, colour modes, tonemappers, trace outcomes, and film enums silently aliased or fell back | internal typed APIs retained permissive default branches after operator validation | Fixed: every malformed internal enum/index/state contracts or produces an explicit diagnostic rather than selecting another represented mode |
| Fisheye existed in the core camera factory but could not be selected through configuration | schema and session projection only named pinhole/thin lens | Fixed: Fisheye is a represented CPU mode; explicit Vulkan declines and `auto` selects CPU at the capability boundary |
| Viewer flags were parsed but the progressive session used hard-coded controls; its later projection still overwrote validated tile, parallelism, bloom, and exposure policy, silently disabled a volume requested without its disk owner, orbit input could leave the typed azimuth domain, frame/session state raced, and shipped shaders were never consumed | CLI/config values did not reach the session; preview-only literals and dependent-feature coercion replaced the session template; unbounded azimuth accumulation made a previously admitted viewer fail on a later refinement; inline GLSL bypassed volume assets | Fixed: one projection preserves every control owned by the chosen backend and every viewer-requested feature, neutralises CPU-only scheduling state under Vulkan, overrides only live viewer state, and passes the complete typed-session authority before initialization, so impossible feature ownership declines; orbit input is reduced to its periodic equivalent, asynchronous refinement/cancellation is synchronised, shaders are file-backed and checked, and pure projection plus real headless-frame witnesses gate the boundary |
| The viewer liveness test embedded local-host performance assumptions | hosted ASan/UBSan/LSan completed one refinement but teardown canceled a healthy second render at the fixed deadline | Fixed: the witness keeps a valid 64x64 Schwarzschild render, asynchronous completion, frame-size, state, and error assertions without duplicating multi-sample performance coverage; the exact sanitizer rerun passes |
| Missing starfield changed the image to grey/analytic background | CPU warning and Vulkan fallback | Fixed: the starfield is a mandatory scene input and absence declines |
| Kerr coordinate round-trip tolerated an intentionally wrong azimuth | forward transform omitted the Kerr phase while inverse subtracted it | Fixed with the exact oblate transform and 1e-12 round-trip gates |
| Horizon-capture test asserted an identity | `x + (total-x) >= total` | Fixed with a positive captured/disk-ray postcondition |
| One permanently disabled unstable-orbit test inflated registration counts | `DISABLED_PhotonSphereRadius` | Removed; the enabled instability test is the represented claim |
| `GTEST_SKIP()` could still make required-profile evidence look green | CMake maps GoogleTest skips to successful CTest skips | Fixed: the required profile force-includes a test policy that converts every `GTEST_SKIP()` to a fatal assertion; portable capability-domain skips remain explicit |
| P1 status had no Bardeen boundary or burn-in | no test referenced Bardeen screen coordinates | Fixed: CPU and Vulkan independently classify ten samples spanning the visible upper curve at 1920x1080, a/M=0.998, below one pixel; the required runtime profile repeats both three times |
| P1 conservation evidence omitted Carter Q and did not span a terminated full ray; the legacy split oracle missed the 1e-10 target | oracle/live witnesses measured only E and L_z, the fixed 20-unit oracle segment stopped near 40M, and the rendered Slang RK4 had no invariant witness | Fixed: an independent adaptive double Hamiltonian oracle crosses the 200M escape surface at a/M=0.998 with E/L_z/Q below 1e-10; the CPU RK45 and actual Slang Cartesian RK4 each traverse a terminated full ray and gate all three live drifts below 1e-4 |
| P1 CPU classifier counted any long-orbit ray as captured | a cinematic winding heuristic was included in the shadow predicate | Fixed: the heuristic and its termination outcome are deleted; the classifier uses a strong-field-only fine step cap, rejects numerical and max-step outcomes, and counts only physical horizon termination as capture |
| FP64 image agreement concealed a controller-induced false termination | the fp64 null threshold exactly equalled the integration target, so a converged minimum step with residual `1.00002e-9` was rejected and a valid ray stopped near 34.75M; the image test then compared black output with escaped fp32 rays while falsely claiming identical step schedules | Fixed: the shared fp64 null policy carries a documented 5e-9 guard band; fp32, compensated fp32, and fp64 compile as distinct non-image precision probes, SPIR-V capability inspection proves only fp64 declares `Float64`, and the same near-extremal state must reach physical escape while conserving E, L_z, Carter Q, and the null constraint |
| Page-Thorne existed as an isolated model while the live paths used a legacy correction | rendered temperature did not equal the full model | Fixed: CPU and Slang render paths use the Page-Thorne flux shape; independent midpoint quadrature with finite-difference derivatives covers retrograde, Schwarzschild, and prograde spins, and a separate large-radius witness requires the Newtonian r^-3 flux tail so a shared determinant-factor error cannot self-validate |
| The public Page-Thorne model accepted a configurable inner edge but always integrated from, and buffered, the ISCO | a truncated disk reported its declared edge while emitting nonzero flux there, and an edge inside the last stable orbit was accepted | Fixed: represented inner edges must be at or outside the ISCO; the zero-torque integral and numerical edge guard both use that exact edge, with independent truncated-disk quadrature as direct P4 evidence |
| Thin-disk coordinate helpers mapped an invalid radius or polar-axis singularity to plausible equatorial coordinates | invalid geometry could become `theta=pi/2` or `z=0` and impersonate a disk intersection | Fixed: both transforms now expose an optional owned domain, reject non-finite/axis inputs, and pass bidirectional round-trip witnesses across positive radii and signed heights |
| The Boyer-Lindquist Kerr oracle clamped spin after copying horizons derived from the unclamped request, and accepted charge although its equations contain no charge terms | exact-extremal and charged requests could combine one metric tensor with another spacetime's identity or horizons; the production ISCO authority separately clamped exact extremality to `0.9999` | Fixed: duplicated derived horizon fields are removed; the oracle accepts only finite uncharged Kerr with `M>0` and `|a|<=M` plus the exact `M=a=0` Minkowski limit, computes its own exact horizons without rewriting, and keeps black-hole-only radii unavailable on the flat limit; the production Kerr ISCO evaluates both exact extremal orientations while charged, massless, cosmological, and super-extremal ISCO requests decline |
| Direct Boyer-Lindquist tensor and Hamiltonian calls could enter the polar singularity and replace `sin(theta)` plus the metric determinant with finite floors | the oracle could return a plausible tensor for an event that its own `IsValid` rejected | Fixed: one finite-event predicate owns the numerical chart for metric, analytic derivatives, curvature, Hamiltonian phase space, and integration validity; polar-margin and non-finite events contract-fail, with no pole-clamped substitute |
| Polarised and canonical Boyer-Lindquist integrators still clamped axis-crossing states, classified every invalid symplectic substep as a horizon hit, silently retained a state when null projection had no real solution, and carried a mutable renormalization counter across independent calls | oracle evidence could mutate an unrepresented trajectory, report a coordinate failure as physical capture, claim an enforced null constraint that was not enforced, or make an identical step depend on call history | Fixed: null/polarisation initial-data constructors are fallible total functions, every RK stage and final state shares the metric's chart predicate, the last valid state is retained, mutually exclusive outcomes distinguish escape/capture/chart exit/constraint failure/exhaustion, implicit nonconvergence and failed null projection are explicit constraint failures, and renormalization cadence is an explicit step argument; the stationary-limit linear null equation is covered directly |
| The oracle interface advertised unused TTESI and frequency-transfer hooks; the latter described a ZAMO as static and returned an invented redshift of one on a zero emitter contraction | dead APIs could be mistaken for evidence of a live regularisation or transfer capability | Fixed: both zero-consumer hooks are removed; live frequency transfer remains owned and independently gated by the actual CPU/Slang Killing-field authority |
| The Slang Kerr ISCO mapped every represented `0 < |a/M| < 0.001` request to exactly `6M` | the device Page-Thorne path silently became Schwarzschild over an arbitrary admitted interval while the CPU evaluated Kerr | Fixed: zero spin is the only Schwarzschild branch, every other admitted spin evaluates the Bardeen formula, and a live SPIR-V parity witness at `a/M=0.0009` rejects interval aliasing |
| Obsolete Novikov-Thorne unit tests mirrored an unused Slang approximation | five green tests exercised only a test-local copy of dead kernel functions and never judged the production Page-Thorne path | Fixed: the unused approximation and its mirror are removed; the independent Page-Thorne comparison lives in the non-render oracle suite and is direct P4 evidence |
| `ShakuraSunyaev` was a public CPU/Vulkan mode but absent from the capability ledger, while its implementation was only `T proportional to r^-3/4` | the named model emitted nonzero flux at its declared zero-torque inner edge, and a Mandatory test merely repeated the same standalone power expression without touching production | Fixed: the CPU authority and Slang kernel use the standard Newtonian zero-torque shape `r^-3[1-sqrt(r_in/r)]`, normalised at the shared `1.5 r_in` temperature scale; core boundary, live CPU, and kernel-parity gates cover it, and the operating model records it as an explicit non-relativistic substitute with no alpha-disk vertical-structure claim |
| Mandatory analytic tests asserted assigned constants or unrelated bounds, and support suites contained unconditional success assertions | photon-sphere cases could not fail the Hamiltonian implementation; a Vulkan placeholder was counted as evidence; colour and memory cases ended in `true` | Fixed: analytic Kerr/Schwarzschild photon-orbit parameters now satisfy the production Hamiltonian null and radial-stationarity equations, ISCO pins use exact values, unconditional/placeholder registrations are removed, and surviving colour/memory/parameter cases assert observable postconditions |
| Test sources could remain green while simulating an unrelated implementation, restating textbook arithmetic, estimating hypothetical resource use, or printing observations without a postcondition | FPS suites fabricated metrics/RK4 or extrapolated single evaluations through assumed CUDA parallelism; geodesic/lensing “benchmarks” compared formulas to themselves; deviation, ADM, and memory suites assigned expected/modelled state directly; Schwarzschild, Kerr, and Reissner–Nordström suites never instantiated the Cartesian Kerr–Schild product; a Christoffel benchmark timed a class expressly documented as non-Kerr; live azimuth/adaptive-step/null probes were advisory | Fixed: disconnected/self-oracle registrations are removed and replaced only where a production-backed independent oracle exists; the black-hole suites now exercise exact Cartesian Kerr–Schild forms, finite-difference derivatives, oblate geometry, horizon/capture authorities, limits, and inverse identities; surviving live probes enforce branch-cut continuity, adaptive-step change, and the named CPU null tolerance; every governed source GoogleTest must contain a direct postcondition, operating-model evidence is checked again at its claim boundary, and negative controls reject assertion-free, obvious no-op, disabled, or `#if 0` tests |
| Spectral correctness largely retested fixture-local Planck/Wien/redshift equations, and a “GPU parity” suite compared two local copies of deleted CUDA-era arithmetic | production Planck and binned-radiance code could regress while independent-looking formulas stayed green | Fixed: the false parity suite is removed; production Planck, numerical Wien maximisation, independent spectral quadrature, Lorentz redshift, exact g-fourth deposition, wavelength boundaries, and integrated colour progression replace it. The independent gates exposed and fixed cancellation in `exp(x)-1` and non-finite/bin-boundary leakage |
| Starfield property loops could pass over an empty catalogue, while reusable phenomenological configurations retained NaN or infinity and fake jet/corona models exposed unsupported physics | loop bodies were the only postconditions and `std::clamp`/`std::max` do not sanitize NaN | Fixed: catalogue witnesses assert non-empty preconditions, retained core configurations reject values outside finite ordered domains without rewriting them, and the local jet plus grey/tinted corona implementations and their self-tests are deleted behind explicit fail-closed operator contracts |
| Registered correctness tests were omitted from the build-time `Mandatory` gate | CI and external full-estate runs included them, but an ordinary build could pass after running only a subset | Fixed: every enabled source GoogleTest is labelled `Mandatory`; the operational verifier compares all parsed source identities with CTest's live JSON inventory, rejects missing or undeclared registrations and non-Mandatory labels, and permits only the explicitly conditional Vulkan session source as an all-or-nothing suite; category labels still support focused runs |
| Deterministic physics witnesses could skip or report zero drift when their own fixture made no progress | angular-momentum and motion-blur tests treated a broken deterministic launch/camera scan as unavailable environment, while five Killing-vector gates initialised drift to zero and did not require an accepted step | Fixed: deterministic fixture preconditions and accepted-progress counts are fatal assertions; optional disk-hit, photon-ring, halation, and escape outcomes are now required by tests that claim them; skips remain only at actual compile/device capability boundaries and become failures in required profiles |
| Volumetric disk emission omitted gravitational/Doppler shift, while the near-extremal gravitational square root could become NaN inside r=2M | only thin crossings called the g-factor authority and both CPU/Slang evaluated an unguarded ergosphere expression | Fixed: thin and volumetric CPU/Slang sources call one exactly-once g^4 authority; finite near-extremal inner-disk, isolated toggle, and end-to-end Vulkan thin/volume witnesses pin the live branches |
| Volumetric transfer restarted and overwrote at every RK45 segment | optical depth depended on the final traversed segment | Fixed: emission and optical depth form one observer-to-source recurrence across accepted segments using invariant comoving path length. The bounded Gaussian grey atmosphere and procedural modulation are explicit; inverse-Compton corona transfer declines |
| “Beam ellipse” collapsed to a circular Gaussian; Vulkan propagated one vector | star filter accepted one sigma and kernel alpha carried one norm | Fixed: both paths propagate two vectors, extract both singular axes and output-plane orientation, and apply an anisotropic tangent-plane Gaussian |
| Beam orientation used the input right-singular vector | formula used `ab+cd` while documenting output position angle | Fixed in live and oracle geometry to the `MM^T` expression `ac+bd`, with a rotated-SVD witness |
| Vulkan beam integration exposed only scalar area, so the device output-orientation formula could regress while P2 remained green | the live kernel alpha carried geometric-mean expansion and only the CPU/oracle rotated-SVD witnesses observed orientation | Fixed: trace and parity probe consume one Slang ellipse-projection authority; a device rotated-SVD witness pins both singular axes, determinant, and output-plane orientation while the live trace gate proves propagated deviation reaches it |
| The oracle beam step treated a Boyer-Lindquist coordinate derivative as a covariant derivative and froze curvature outside the ray's integration tableau | the exact radial and circular Schwarzschild null-congruence solutions exposed the missing connection/stage coupling | Fixed: central ray and four covariant Jacobi columns share one RK4 tableau; radial and photon-sphere screen axes agree with independent closed forms to one part in 1e6 |
| The claimed symplectic-structure gate bounded nearby-ray separation instead of testing the canonical two-form | bounded separation is not implied by symplecticity, and the explicit kick-drift-kick map is not symplectic for the non-separable Kerr Hamiltonian | Fixed: the fixed-step oracle composes symmetric implicit-midpoint maps; a finite-difference gate enforces `D(Phi)^T J D(Phi) = J`, order accuracy is measured against the independent adaptive Hamiltonian solver, and variable-step/null-projection stabilisers are expressly outside the symplectic claim |
| The flat CPU/Vulkan point-catalogue parity fixture changed the metric to Minkowski but retained Kerr spin and disk state | the shared typed session boundary correctly declined both sessions | Fixed in evidence: the fixture now projects a coherent massless, spinless, diskless Minkowski scene instead of weakening product validation |
| Point-star sampling scanned 100,000 stars per escaped pixel | catalogue had no live acceleration structure | Fixed: deterministic latitude/longitude CSR index, exact exhaustive-oracle agreement, shared CPU/Vulkan candidate semantics, and bounded residency |
| Point-star mode was numerically nonzero but display-invisible at its operator default | scale 1 produced a disk-free Radeon PNG with channel ranges 1/1/0 and only three colours; the Vulkan test passed because near-black point mode differed from the texture | Fixed: the relative-flux zero point is display-calibrated at the anti-flicker-tested scale 100, disk-free live gates require a bounded sparse lit-pixel fraction, and the IMAX verifier enforces sparse dark-field morphology rather than accepting disk structure |
| P3 attestation admitted only the IMAX frame although the criterion quantifies both 1080p and 5616x4096 | no independent 1920x1080 artifact, transcript segment, wall time, dimensions, or morphology was required | Fixed: the physical runbook renders the identical governed sparse-star scene at both exact resolutions and the verifier independently binds, decodes, and morphology-checks each artifact and typed-session event |
| P5 aberration evidence covered only motion along the view axis | the analytic witness reduced the Lorentz transform to one dimension | Fixed: an independent photon-four-vector Lorentz oracle covers mixed-sign three-axis velocities, varied ray directions, and near-luminal finite worldlines in addition to both lens models and the governed physical scene |
| Doppler toggle test attributed Kerr lensing/frame-dragging asymmetry to emitter motion | image-half ratio was the sole oracle | Fixed: the gate separately measures observed asymmetry and the isolated `(g/g_grav)^4` emitter factor; off is exactly zero in the isolated measure |
| Leak checking had no reproducible operational profile | ad-hoc ASan runs disabled leak detection around Vulkan | Fixed: GCC ASan/UBSan/LSan preset and CI job; one narrow, printed suppression names the repeatable 128-byte Vulkan-loader/driver process-lifetime allocation |
| Hardware script accepted software devices and called 4096x2864 “IMAX-class” | no physical-device precondition or attestation | Fixed: exact device matching, software rejection, readiness check, full log/hash/JSON attestation, and the specification's exact 5616x4096 frame |
| External-domain records had no enforced schema and device inventory omitted driver identity | a handwritten JSON result could call llvmpipe physical, call Dozen native Windows, or name an IMAX output without decoding it | Fixed: a Mandatory negative-control verifier rejects software devices, domain mismatch, wrong dimensions, false test state, and tampered hashes; Vulkan inventory includes driver/vendor/device/API identity; the physical runbook decodes 5616x4096 and verifies its own record |
| Physical test, readiness, and revision claims were unbound result booleans | the transcript was hashed, but the verifier did not require its CTest result or inventory to support `test_estate_passed`, `runtime_ready`, and `source_revision` | Fixed: the runbook emits a hashed non-skipping JUnit report and explicit revision/readiness markers; the verifier binds its case count to CTest, requires `linux-ci`, and cross-checks the selected device and platform against the hashed inventory |
| Physical/runtime records could name clean `HEAD` while executing a dirty or stale binary tree | revision text and artifacts did not prove which configure/compiled authority the JUnit and device run consumed | Fixed: every runtime producer now bundles the clean configure-time receipt and its JUnit must contain the test that compared that receipt to the compiled authority; revision and model-domain bindings are revalidated during admission |
| “Full test estate” meant only a green boolean and case count | a partial JUnit file could omit registered obligations while retaining plausible metadata; a declared floor could drift below the real source estate unnoticed | Fixed: producers hash CTest's JSON registration inventory; the local build independently parses 732 source GoogleTests and refuses to fall below the policy floor, while external admission requires that same floor and exact equality between every enabled registered test name and the non-skipping JUnit case set |
| Full-frame scene claims were trusted from runbook metadata rather than bound to the executed session | structured unrelated PNGs and claimed P3/P5 scenes could satisfy the record without proving that bundles, stars, motion, and ThinLens reached dispatch | Fixed: the typed session emits canonical JSON; the verifier requires one exact event at 1920x1080 and one at 5616x4096, cross-checks every scene field, proves the actual 100,000-star catalogue reached Vulkan with beams enabled, and binds device, budget, dispatch completion, wall time, output name, and terminal state in each hashed transcript segment |
| E1 claimed all nine CPU-renderable metrics from registry metadata while live render evidence exercised only Morris-Thorne | registry completeness did not prove every advertised CPU factory/session path completed | Fixed: a registration-driven product witness constructs legal parameters and completes a frame for every CPU-supported registry row without a hand-maintained metric list |
| E2's machine ledger omitted its explicit decline witnesses | unsupported polarised volume/temporal/Vulkan and two-sheet requests were tested but not required evidence for the criterion | Fixed: configuration and typed-session decline witnesses are direct E2 evidence, alongside live Kerr, Schwarzschild-limit oracle, physical Stokes, and film-output gates |
| E2's cross-chart live/oracle witness exercised Kerr while the Schwarzschild rigid-transport witness remained oracle-only | the live thin-disk Stokes path did run Schwarzschild, but did not carry the same invariant comparison through a physically completed ray | Fixed: the existing cross-chart gate now completes escaping Schwarzschild and Kerr rays and compares their initial/final Walker-Penrose invariants between live Kerr-Schild Cartesian and independent Boyer-Lindquist paths |
| CPU disk motion blur rotated an axisymmetric stationary crossing through a Euclidean line-of-sight factor | the synthetic phase change was not a covariant time-dependent emissivity or shutter integral | Fixed: the operator identity remains parseable but all CPU/Vulkan motion-blur requests decline until a time-dependent emissivity and covariant temporal transfer are represented |
| CPU motion-blur evidence reimplemented a synthetic temporal formula disconnected from physical emissivity | source inspection and exact sample counts could not make the Euclidean phase rotation covariant | Fixed: the sampler and its mirror tests are deleted; schema, CLI, configuration, and typed-session witnesses require every temporal-blur request to decline before tracing |
| Two-sheet wormholes were described only in comments | the operator could request a Morris-Thorne scene but could not state the required topology, allowing one-sheet output to be misunderstood | Fixed: `OneSheetCapture` and `TwoSheet` are explicit schema/CLI/typed identities; the former is represented, while the latter declines before render |
| GLFW repeat actions were interpreted as key releases and non-finite pointer input could poison camera state | the handler treated only action 1 as pressed and accepted NaN/Inf | Fixed: press/repeat/release semantics, finite pointer/scroll guards, camera update, and refinement restart are Mandatory-gated; native event delivery remains a separate attestation |
| Viewer state exposed an unwritten restart flag, unused camera velocities, and an attached GLFW window pointer that the CLI—not the viewer—owned | downstream callers could infer observable refinement/window capabilities from dead fields | Fixed: the actual atomic restart lifecycle now publishes and consumes the synchronised restart state under the input gate; dead velocity and window attachment surfaces are removed |
| P2900/P2996, P1–P6/E1–E4, and external-profile boundaries were prose-only at runtime | C++26 mode could be mistaken for native contracts/reflection and installed volumes carried no complete acceptance ledger | Fixed: compile-time language facts distinguish native features from checked-macro/explicit-schema substitutes; all ten P/E criteria, 24 dimensions, and 29 capability contracts are build-gated, installed, readiness-required, and exposed by `info capabilities` |
| CPU/Vulkan trace termination and step lengths were fixed in absolute coordinates while public black-hole and exotic geometry scales varied by orders of magnitude | at large M the 200-unit escape sphere could coincide with the horizon, while Morris-Thorne `b0` and Alcubierre `R`/`1/sigma` were ignored entirely; observers could start inside the one-sheet capture throat, the kernel could classify distant observers as escaped before one step, and unresolved walls could be advertised as rendered | Fixed: one registry authority governs observer bounds and both trace paths with `M`, `b0`, or Alcubierre `max/min(R,1/sigma)`; the escape sphere encloses the exterior observer in the asymptotic region, device stepping is scale-covariant, and `0.1 <= sigma*R <= 100` keeps the radius and inverse-wall scales inside the fp32-resolved envelope |
| A valid Mandatory test could be substituted for an unrelated P/E or operating-model witness | evidence shape and existence checks did not preserve the reviewed semantic mapping | Fixed: the ten P/E evidence sets are explicit verifier policy, and canonical digests pin the complete acceptance, dimension, and capability sections; negative controls mutate both evidence identity and section semantics and require rejection |
| Release packaging could retain a complete attestation receipt while weakening the product built from that revision | `BUILD_TESTS`, the Mandatory gate, warnings-as-errors, contract enforcement, Vulkan/Slang, SPIR-V validation, and the native viewer remained independently optional | Fixed: a two-phase build-policy verifier makes all of them non-negotiable, permits only the Release configuration, checks both configured flags and realised target/tool topology, and supplies a Mandatory negative-control CTest that rejects every weakened variant |
| Qualification CI could execute different build inputs for the same source revision | remote Actions and FetchContent used movable tags, the Linux Slang archive had no checked digest, and Windows selected the latest third-party SwiftShader at runtime | Fixed: every Action and fetched source has an immutable commit identity, direct Slang/SwiftShader archives have versioned URLs and checked-in SHA-256 values, native Vulkan toolchain/device preflights fail early, and repository governance rejects movable actions/dependencies, latest rasterizer selection, and unchecked direct downloads |
| Pull-request integration was inseparable from resource-heavy rendering qualification | opening any pull request launched the full Linux, sanitizer, Windows, and macOS runtime estates even when the needed evidence was compilation and authority governance | Fixed: pull requests and explicit integration dispatches now compile the complete strict topology on Linux, Windows, and macOS and run exactly nine non-render authority controls while proving no promotion receipt exists; full Mandatory/runtime and native attestation jobs are push-only, and repository governance rejects either side crossing that boundary |
| Release install/package creation could bypass execution of the configured Mandatory target | `cmake --install` and CPack consumed an already-built tree without proving that its current tests or artifacts were the ones that passed | Fixed: the Mandatory target now owns CTest execution and writes a deterministic receipt only when its zero-skip JUnit name set exactly equals live registration; the receipt hashes seven test executables and the complete release product, and install/CPack rehash both the build tree and installed volume before succeeding |
| Requiring the final Mandatory receipt inside the tests that issue it made first release promotion self-dependent and allowed an old staged receipt to influence the next run | the gate deleted its canonical stamp but runtime tests and install checks consumed a receipt that could exist only after they passed | Fixed: each gate removes both canonical and staged receipts before CTest; release pre-gate operational branches require install/readiness to fail closed and prove the resource override is ignored, the C++ release validator is exercised against an explicit isolated volume without changing packaged lookup, and only a newly passed exact-estate receipt is staged afterward |
| An installed release could pass packaging and then lose or alter its Mandatory receipt, executable, or resources before startup | build/install verification ended before the product's operational trust boundary | Fixed: release readiness and render/view initialisation parse the installed receipt and rehash the running executable plus all 12 installed resources; exact seven-test/13-product inventories, clean revision, zero skips, model/alignment digests, sizes, and SHA-256 values are required, while missing receipt and same-size tamper controls fail closed |
| Install-time Python verification and runtime C++ verification could each be correct without ever agreeing on the same relocated volume | the installer rehashed files but did not execute the installed product authority | Fixed: release install/CPack now invokes the relocated binary's non-rendering `info capabilities` path after external rehashing and requires schema 3 plus aligned ultimate-ideal state, so installation fails unless the executable accepts the exact installed authority/product set |
| `SIRIUS_RESOURCE_DIR` could redirect a strict executable to a forged receipt and matching operator-controlled resources | development resource injection remained live in the release/qualification trust boundary | Fixed: the override is compiled out of qualification and release resource candidate selection, so all alignment/build-gate checks and subsequent renderer/viewer loads remain anchored beside the running executable; development retains explicit-root fail-closed behavior only for local diagnosis |
| Resource-consuming test executables lived outside the strict candidate volume and depended on the now-disabled development override | clean qualification correctly ignored the override, so application authority tests failed before the gate could issue evidence and render tests would have inspected a different tree | Fixed: application and render suites execute beside the candidate binary, consume its exact 12-resource volume, and share the gate's remove-test-stage receipt lifecycle; source governance rejects relocation or restoration of the override |
| Hardware admission required at least 891 JUnit cases after the governed strict estate had become 784 | a historical test count survived as policy, so an exact, complete, zero-skip current qualification run could never be admitted | Fixed: hardware admission uses the operating model's governed source-available floor, exact JUnit-to-live-CTest identity equality, the qualification gate's matching identity set, and named hardware semantic witnesses; no historical estate size can veto a complete current revision |
| Individually valid external attestations never became a configure/build/runtime authority | the verifier could approve files, but no complete-set admission existed; overall readiness could not distinguish a development preflight from an aligned release | Fixed: release configuration requires all eight domains exactly once at the current clean revision; a deterministic receipt is checked on every build, statically required by the compiled release authority, installed, reported, and compared byte-semantically before render/view initialisation; development remains explicit and non-packageable |
| External attestations could be generated by a weaker development artifact and later promote a different release product from the same source revision | the verifier accepted development receipts and carried no actual candidate executable or gate-generated test artifacts in the bundle | Fixed: only clean strict `qualification` receipts and gates are admissible; qualification enforces the full release-equivalent product but allows pending domains and cannot package. Each record carries the copied candidate bytes, exact seven-test/13-product gate, gate JUnit/log, rerun JUnit/live inventory, alignment receipt, model digest, and revision under one cross-checked hash/test-identity chain; release separately retests and gates its final product |
| Development readiness said `ready: true` while the same response said the strict ideal was unsatisfied | local resource readiness and aligned-system readiness shared one ambiguous boolean | Fixed: top-level `ready` and its exit status are reserved for a satisfied revision-bound ideal; development exposes a separate `evidence_generation_ready` preflight so external evidence can still be produced without calling the system aligned |
| Runtime alignment reporting collapsed the receipt to an unactionable `0/8` count | operators could not identify the exact missing authorities from the installed volume | Fixed: the validated runtime authority carries sorted admitted, pending, and required IDs; JSON and human readiness expose the partition, and relocated-runtime plus unit controls prove it is the exact eight-domain model set |
| Operating-model external profiles and release admission were separate policy lists | the viewer capability said `viewer-native-window` while the verifier admitted `viewer-native-window-input`; both authorities could validate independently | Fixed: exact capability-to-domain mappings are governance- and runtime-validated, admission derives the eight-profile set from the model, and every receipt carries the model SHA-256 checked by the compiled authority |
| Native Windows Vulkan and macOS/MoltenVK records had a weaker evidence chain than Radeon/Dozen | physical-device fields plus arbitrary PNG/log files and result booleans could pass without a hashed JUnit estate, selected-device inventory, revision/readiness transcript, or required native preset | Fixed: every runtime-device domain now requires the same non-skipping JUnit/inventory/revision/readiness/transcript chain; domain-specific presets are enforced and positive plus missing-artifact controls cover both native routes |
| Native Windows Vulkan and macOS/MoltenVK had verifier contracts but no executable evidence producer | operators had to hand-assemble a record and could not follow one source-bound route from native build through device selection to verification | Fixed: one clean-tree native producer chooses exactly one physical device, rejects Dozen/non-MoltenVK substitutions, binds configure/build/JUnit/inventory/readiness/probe output, rechecks revision identity, and verifies the record before publication |
| Non-IMAX physical runtime records accepted any hashed file with a `.png` suffix | native Vulkan success could be asserted with renamed text or spatially collapsed bytes | Fixed: every physical runtime PNG is fully decoded and must have represented dimensions, dynamic range, colour structure, and more than one row signature; IMAX retains its stricter sparse-point morphology gate |
| Application JSON types and conversion logic lived in the render layer and threw across the app/render seam | `render_config.h` mixed wire-schema strings with film and session-domain types | Fixed: app owns explicit JSON codecs, render owns film/session types, and one `std::expected` adapter projects validated strings into closed enums without an exception crossing the boundary |
| CPU backend source was compiled by the render target and the repository retained two unused vendor trees | CMake carried a deferred re-home comment; no first-party build/include referenced GLM or Dear ImGui | Fixed: `sirius_backend_cpu` is always available and singly owns the tracer; render links it as a target; GLM/ImGui are removed; executable structure governance prevents ownership, layer, boundary-name, and vendor-set regression |
| First-party physics and oracle code depended on non-standard `M_PI` macros and exposed mixed camel/snake data records | local compatibility defines and Java-style configuration/result members persisted after the C++26 port | Fixed: `<numbers>` is the sole standard pi source; app/render/backend/oracle boundary records use governed snake-case identifiers, public templates at the JSON seam are concept-constrained, and raw fixed-size tracer arrays are `std::array` |
| Render configuration returned an ambiguous boolean, accepted invalid state, and reported ordinary resource/output failures by throwing after a worker launched | validation and error ownership were deferred into asynchronous initialisation | Fixed: `Configure` returns `std::expected`, validates synchronously, and the render state machine consumes typed configuration/IO/physics errors without deliberate exceptions on the live path |
| The interactive viewer compiled its 450-line threaded implementation through every including translation unit | the legacy header-only port had never acquired target ownership | Fixed: the declaration-only header and singly owned `.cpp` separate the public surface from refinement implementation; copy/move semantics and discarded-return intent are explicit |
| C translation units and standard size/integer names depended on compiler defaults and global aliases | only the C++ dialect was pinned and first-party code mixed `size_t`/`uint32_t` spellings | Fixed: mixed targets pin C17, all first-party C++ uses `std::size_t`/fixed-width types, and both GCC 14 and Clang 21 compile under warnings-as-errors |
| A Dozen render completed and wrote output, then the render worker called a TLS destructor from an unloaded WSL D3D12 runtime | GDB placed the thread-exit program counter in the unmapped `libd3d12core.so` range; retaining both D3D12 mappings made the same dispatch exit cleanly | Fixed: the Vulkan backend detects `VK_DRIVER_ID_MESA_DOZEN`, retains the D3D12 wrapper and core for process lifetime, and gates dispatch plus device destruction on a worker thread |
| Operator render scripts could report success from `tail` after a failed renderer and choose an arbitrary stale binary | pipeline status belonged to `tail`; `find ... | head -1` selected an unspecified build | Fixed: exact/explicit binary selection, status propagation, non-empty artefact checks, and failure/no-output negative controls |
| EXR metadata was constructed but discarded, finite fp32 radiance could overflow half storage, and direct writers accepted NaN/Inf | writer header and advertised color-space claims did not match the file | Fixed: fp32 channels, embedded Sirius metadata, conservative color naming, and finite checks at display, session, Vulkan, PNG, and EXR boundaries |
| Product PPM output bypassed the reusable writer and applied a separate gamma-2.2 loop | the direct writer could be correct while the session route transfer-encoded differently | Fixed: the session delegates RGBA display-linear output to the owned PPM writer, the writer shares the sRGB boundary authority, and byte-exact direct plus completed product-route witnesses are P6 evidence |
| Native CI passed its full test estate but could not issue the promised build attestation | CTest placed relative JUnit outputs under each `--test-dir`, while the verifier looked in the repository root | Fixed: Windows and macOS attestations consume the exact binary-tree JUnit paths; publication CI is the executable witness |
| Native build records trusted a caller-supplied revision and could relabel a stale JUnit file | the record contained no configure/compiled-authority witness tying the test binary to the clean checkout | Fixed: the producer derives and cross-checks the clean Git revision, bundles the configure receipt, and requires the JUnit case that compared that receipt with the compiled authority |
| Native compilation claims were coupled to the rendering-inclusive promotion gate | a build-only domain could not be evidenced without executing unrelated runtime/image tests, while its JSON gate did not bundle the other compiled test binaries | Fixed: a distinct non-promoting native-build gate binds the full governed registration, exact nine non-render authority cases, all seven compiled executables, all thirteen products, clean revision, and alignment receipt; dispatch-only publication states runtime false, while every physical/runtime domain retains the complete Mandatory estate |
| External runbooks discovered a wrong host, software/ambiguous device, missing GUI/input route, or stale candidate only after starting expensive qualification | environmental mistakes could waste the very rendering/test resources needed for physical evidence, and an informal preflight could itself be mistaken for attestation | Fixed: one exact-candidate preflight permits only `info system`/`info readiness`, cross-checks the compiled receipt/resources/host/device/driver/memory/tools, rejects render/view/test commands, writes only outside the clean tree, and emits an explicitly non-promoting report with no admissible domains or artifacts |
| The WSL2 hardware and native-viewer producers rebuilt and reran the same complete qualification estate on the same candidate and Radeon/Dozen device | the second estate added resource cost but no independent source, product, registration, or physical-device authority; skipping it informally would permit stale or cross-device substitution | Fixed: viewer reuse accepts only an independently verified same-revision physical-Radeon record, rehashes its candidate/receipts/gate/JUnit/registration/products against the live volume, requires exact WSL2 device identity, and carries the upstream transcript; only the newly published viewer frame and native callbacks are executed separately |
| Native viewer input had a verifier contract but no evidence producer | headless state-machine tests could pass, but no runbook could prove GLFW callback delivery | Fixed: an opt-in callback transcript and clean-tree X11/XWayland runbook deliver keyboard, drag, scroll, and Escape through the host window system and hash the resulting record |
| E3 viewer attestation could bind callbacks to a CPU frame or unrelated device without a complete estate | the record lacked a published Vulkan-frame event, selected physical inventory/readiness, and non-skipping JUnit artifact | Fixed: the runbook selects the named Radeon, gates the full estate, waits for a native-window Vulkan refinement frame before injecting host input, and the verifier binds frame, callbacks, readiness, inventory, JUnit, log, and revision by hash |
| The compiled capability authority exceeded MSVC's narrow string-literal limit | the 18,303-byte operating model was emitted as one raw literal and failed with C2026 | Fixed: CMake emits independently bounded 8,000-byte chunks; runtime validation reconstructs and semantically compares the complete authority |
| Sanitizer publication limits assumed the faster review host | both pushed streams reproduced a 120-second relocated-volume timeout and a viewer teardown cancellation even though the same evidence passed locally | Fixed without dropping evidence: the specification-minimum 128x128 relocated render and every obstruction remain, the viewer witness removes redundant sampling load, and bounded test/job liveness windows admit slow instrumented hosted CPUs |
| Compiler estates interfered through fixed `/tmp` render filenames | concurrent GCC and Clang runs could remove the other process's completed EXR between write and assertion | Fixed: render-session probes own collision-resistant RAII temporary directories; simultaneous GCC/Clang reproduction now passes |
| Historical byte-identity claims had no executable consumer | no test referenced the baseline directory | Claim withdrawn; tapes remain forensic evidence only |
| The upstream specification still required an unimplemented tape-identity gate, called macOS CI-only, and tied correctness to a stale branch name | current runtime admission, statistical/oracle evidence, and repository state followed different authorities | Fixed: the target now requires executable semantic/oracle/output evidence, external MoltenVK admission, and one exact clean revision; branch integration remains an owner workflow rather than correctness evidence |

## 4. Enforced operating model

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

The same ledger now contains 28 explicit capability contracts. Each has one
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
| Spherical de Sitter/Kottler sector | `bounded`; exact positive-Lambda metric and both Kottler roots, with pure de Sitter horizonless for capture, Schwarzschild-de Sitter restricted below Nariai, observer at or below `0.99*r_c`, and directional boundary no later than `r_c` |
| Disk emission outside Schwarzschild/Kerr | `fail_closed`; requires `--no-disk` |
| Morris-Thorne one-sheet dark-throat scene | `supported` on CPU/Vulkan |
| Morris-Thorne two-sheet continuation | `fail_closed` through `TwoSheet` request identity |
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
controls prove both operator scripts propagate a renderer failure and reject an
exit-zero renderer that emits no file.

`linux-ci` is the non-skipping runtime profile. Configuration fails without the
Vulkan backend and compiled kernels; CTest fails without a ready device, a
governed multi-sample dispatch, or the repeated CPU/Vulkan P1 classifiers.
Portable source profiles may omit Vulkan, but cannot be cited as evidence for
this profile.

## 5. Evidence semantics by profile

Live attestation status is deliberately not frozen into this source file. An
attestation names the exact source revision that it verifies, so committing a
new status count would create a different revision and immediately invalidate
that count. The deterministic alignment receipt generated from independently
verified, same-revision bundles is the only authority for the live admitted and
pending partition. The table records durable implementation and evidence
boundaries instead.

| Profile | Status | Evidence boundary |
|---|---|---|
| Revision-bound release alignment | EXTERNAL LEDGER AUTHORITATIVE | A source-only qualification configure intentionally admits no external domains. Development artifacts remain inadmissible; release configure reports every domain absent from the supplied same-revision ledger, and packaging/initialisation remain fail-closed until all eight verified records exist. |
| Pull-request/integration boundary | GOVERNED NON-RENDER PATH PRESENT | Linux, Windows, and macOS compile the complete strict topology and execute exactly nine authority controls without creating a Mandatory receipt. Pull requests cannot publish evidence; explicit dispatch may issue only the precisely scoped Windows/macOS compilation domains. |
| Configure/compile/build, GCC 14 | STRICT NON-RENDER PATH PRESENT | Qualification binds all seven test executables and 13 products under warnings as errors, then runs the exact nine authority controls without issuing a Mandatory receipt. Execution outcomes belong to the workflow record for the tested revision. |
| Configure/compile/build, Clang 21 | STRICT NON-RENDER PATH PRESENT | Qualification also emits and validates every Slang kernel, binds all seven test executables and 13 products under warnings as errors, and runs the exact nine authority controls without issuing a Mandatory receipt. Execution outcomes belong to the workflow record for the tested revision. |
| Relocatable CPU volume | CURRENT INITIALISATION CHECKED; FULL PROFILE SNAPSHOT ONLY | Installed authority/resource readiness and tamper rejection pass after relocation; the relocation render snapshot predates this delta |
| Viewer-disabled build/install | PRE-ALIGNMENT FULL-PROFILE SNAPSHOT | Capability-specific install verification, readiness, and fail-closed `view` were previously exercised; no current-delta full profile is claimed |
| CPU physics/render path | ORACLE DELTA READY; RENDER PROFILE SNAPSHOT ONLY | Current Page-Thorne and canonical-symplectic oracle corrections pass on both compilers; near-extremal, shadow, volume, and film-output executions remain the immediately preceding snapshot under the no-render validation constraint |
| Vulkan on WSL2 software and physical devices | PRE-ALIGNMENT FULL-PROFILE SNAPSHOT | Prior llvmpipe/Radeon-Dozen evidence covered dispatch and bounded scene semantics, but historical records are not admitted for a later source revision. |
| Interactive viewer refinement | PRE-ALIGNMENT PREFLIGHT | The prior source opened a GLFW/OpenGL XWayland window, published a Radeon Vulkan frame, and received host-delivered input; native-viewer admission is asserted only by a verified same-revision external record. |
| GCC ASan + UBSan + LSan | PRE-ALIGNMENT FULL-PROFILE SNAPSHOT | The historical 904-test sanitizer result is not promoted to evidence for a later source revision. |
| Physical Radeon 780M | PRE-ALIGNMENT WORKTREE SNAPSHOT | The historical 907-test and 697-test Radeon/Dozen runs plus both exact sparse-star frames are revision-specific and cannot be admitted for later source. |
| Native Windows build/runtime | BUILD PRODUCER PRESENT; LIVE RECORD EXTERNAL | Explicit Windows dispatch can emit a clean-Git compilation record binding the full registration, exact non-render authority estate, all executables/products, and receipt. Native Vulkan remains a separate full-estate physical-host domain and rejects Dozen. The alignment receipt, not this table, states whether either record exists for a revision. |
| macOS build/MoltenVK runtime | BUILD PRODUCER PRESENT; LIVE RECORD EXTERNAL | Explicit macOS dispatch can emit a clean-Git compilation record binding the full registration, exact non-render authority estate, all executables/products, and receipt. MoltenVK remains a separate full-estate physical-host domain. The alignment receipt, not this table, states whether either record exists for a revision. |

## 6. Remaining limitations (not hidden as “green”)

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
  asymptotic sheet with a dark throat; an explicit `TwoSheet` request declines.
- The exact 5616x4096 memory plan, catalogue/index residency, and sublinear
  candidate query are Mandatory-gated. The physical runbook now requires the
  same scene at both 1920x1080 and 5616x4096 to combine beam-filtered
  100,000-star sampling, nonzero camera velocity, and a finite-aperture lens
  under the 2048 MiB cap. A preflight of
  that exact scene completed on the Radeon but exposed a finite-aperture
  projection defect: all channels spanned at most 14 code values and only 62
  colours survived at full resolution. That frame was rejected. The shared
  CPU/Slang projection now preserves the requested field of view, a physical
  Radeon regression retains a 56.0 per cent shadow with a 0--4.54 linear
  luminance span, and the attestation verifier decodes the PNG to reject
  collapsed output. The corrected current-worktree render then completed on the
  Radeon in 16,841.9 seconds and 32,736 governed dispatches. Its 5,013,725-byte
  5616x4096 RGBA PNG has SHA-256
  `0a83eaf51a1d9c6293cc0ae1d11fe6ffcc941859a618ec23dedcaf3852343bc8`;
  every RGB channel spans 216 code values and 66 distinct colours survive.
  Visual review resolves the central shadow, asymmetric lensing, photon-ring
  structure, and warped emission surface rather than the rejected bands, so it
  closes the ThinLens projection diagnosis. It does not close P3: disk emission
  supplies the image structure, with 45.5 per cent bright channels and only
  54.0 per cent dark channels. A disk-free Radeon probe then exposed that the
  point catalogue's scale-1 default quantised to channel ranges 1/1/0 and three
  colours. The display-calibrated scale-100 fix produces ranges 151/140/129,
  at least 64 colours, 99.96 per cent dark channels, a 0.020 per cent bright
  channel fraction, and bright sources across 8.2 per cent of rows. That exact
  moving ThinLens scene completed in 20.0 seconds and 88 governed dispatches;
  both strengthened Vulkan gates pass on the Radeon. The verifier now requires
  this sparse disk-free morphology plus one hashed typed-session event, proves
  the actual catalogue reached Vulkan with beams enabled, and binds device,
  budget, completion, and output filename in the same transcript segment. The
  existing 5616x4096 diagnostic predates those requirements and is rejected as
  P3 evidence. Admission remains pending until the clean-revision runbook
  produces both exact sparse-star records.
- Vulkan volumetric transfer deliberately caps `volumetric.samples` at 128 per
  geodesic segment to protect the first dispatch from an unbounded watchdog
  exposure. The CPU accepts the schema maximum of 4096; explicit Vulkan requests
  above 128 decline and `auto` logs the boundary before selecting CPU.
- `backend.enableDenoiser` and nonzero legacy `backend.cudaDevice` are retained
  as compatibility sentinels. They decline instead of implying unavailable
  denoiser/CUDA capabilities; neither is a DNGR parity criterion.
- The viewer's strict parsing, projection, refinement render, frame publication,
  cancellation, and press/repeat/release/pointer/scroll state transitions are
  gated. The current worktree has also opened the GLFW/OpenGL window and
  received host-delivered callbacks on WSLg/XWayland. That older preflight is
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

## 7. Validation record

The counts below are refreshed only after the final source state completes the
corresponding profile. They are not projected from registration counts.

- 2026-08-27 alignment delta: GCC 14 and Clang 21 compiled the new
  receipt/admission/runtime-authority path under warnings-as-errors. All three
  non-rendering alignment tests passed on both compilers (incomplete/ambiguous-
  set rejection, installed-receipt tamper rejection, compiled/staged authority
  equality). Operating-model evidence/profile governance passed on both; the
  expanded per-domain attestation controls and native-producer host/device
  negative controls passed locally because they are compiler-independent.
  A real three-test non-rendering selection proved CTest JSON inventory/JUnit
  name-set equality. Label/structure governance, receipt installation, and a
  relocated initialization probe passed with strict readiness false,
  evidence-generation readiness true, the declared 700-case external floor,
  and 0/8 domains admitted; dirty release admission failed as required. Per
  owner direction, no rendering test was run to validate this delta. The full-
  estate counts below therefore remain the immediately preceding implementation
  snapshot and are not promoted as evidence for the new alignment code.

- 2026-08-27 non-render oracle delta: the Page-Thorne production profile passed
  an independent midpoint-quadrature/finite-difference comparison for three
  spins on GCC 14 and Clang 21. A formerly non-probative symplectic test first
  exposed an `8.2e-6` canonical-form defect in the explicit non-separable map;
  after replacement with implicit midpoint/Yoshida composition, all nine
  symplectic-oracle cases passed on both compilers, including the direct
  canonical-form, `1e-10` conservation, and independent order-accuracy gates.

- 2026-08-29 physical-model closure delta: back-propagating the Page-Thorne
  determinant exposed a shared production/oracle defect: the Newtonian
  `r^-3` prefactor had been multiplied directly by a radial factor that itself
  tends as `r^-2`, producing a false `r^-5` tail. Production now forms the
  dimensionless `(2/3) r^2 q(r)` correction, the independent oracle uses the
  physical `q(r)/r` flux shape, and a separate large-radius slope witness
  prevents recurrence. The same delta removes the unused kernel tidal
  approximation, makes the accepted central-ray Hermite segment authoritative
  for disk events, volume samples, polarisation, and Jacobi stages (including
  tangent and equal-endpoint-sign cubic plane roots), introduces
  scale-aware/null-gated adaptive GPU step doubling, and declares every
  retained emission closure or fail-closed boundary in the 26-capability
  operating model. On GCC 14 and Clang 21, all 365 core and 121 oracle tests
  passed together with the named non-render Page-Thorne, integrated-blackbody,
  and adaptive-integrator storage-buffer probes; all SPIR-V precision rungs plus
  CUDA and Metal emissions compiled and validated. Rendering test binaries were
  compiled; no rendering test completed or contributed evidence. Development
  inventory is 741 Mandatory
  registrations: 732 governed GoogleTests and nine operational CTests.

- 2026-08-27 spectral/configuration/build-gate delta: 26 production spectral
  validation/utility cases, nine binned-radiance cases, and 74
  starfield/corona/turbulence/disk/jet cases passed on GCC 14 and Clang 21. Both
  warnings-as-errors profiles compiled the application plus core, backend, and
  rendering test targets; rendering tests were compiled but not executed. The
  source/live-inventory, release-policy, and artifact-bound build-gate
  operational controls passed on both profiles. The complete Mandatory runner
  was intentionally not invoked because it includes owner-only rendering
  validation, so no build-gate receipt or full-estate pass is claimed.
  The current 741 development registrations (732 GoogleTests plus nine operational CTests)
  remain an inventory count, not a full-estate pass claim. Disconnected
  simulated/self-oracle registrations were removed; independent product-backed
  replacements, finite-domain controls, and build/install receipt controls were
  added where they establish a real postcondition. No rendering test was
  executed for this delta.

- 2026-08-27 release-runtime receipt delta: GCC 14 and Clang 21 compiled the
  executable-path and SHA-256 authorities plus release initialisation checks
  under warnings-as-errors. Five NIST/streaming SHA-256 witnesses plus eight
  authority/path witnesses passed on both compilers; the release-volume
  fixture accepted an exact 13-product receipt and rejected a same-size resource
  alteration, a skipped test, and a missing receipt. The source/live-inventory,
  release-policy, and artifact-bound build-gate operational controls also passed
  on both profiles, as did the development installed-receipt tamper obstruction.
  The Mandatory command was inspected to confirm remove, run, then stage order.
  These were non-render selections; no full-estate pass or release receipt is
  claimed, and the strict-mode pre-gate branches remain pending the external
  clean-revision configuration.

- 2026-08-27 qualification-promotion delta: an end-to-end trace found that the
  external runbooks still built ordinary development artifacts, while admission
  accepted their receipts and a later release build produced a different
  product. The sole evidence producer is now clean `qualification` mode: it
  enforces the release-equivalent configure/product policy, packaged-resource
  lookup, and Mandatory gate while allowing pending domains and disabling CPack.
  Every record directly carries and cross-checks the candidate executable, all
  12 runtime resources, exact seven-test/13-product gate, gate-generated
  JUnit/log, independently rerun JUnit/live inventory, alignment receipt,
  operating-model digest, and source revision. Development and release-mode
  evidence producers, substituted candidates/products, missing gates, divergent
  test identities, and dirty qualification are negative controls. GCC 14 and
  Clang 21 compiled the executable and application authority under
  warnings-as-errors; the same nine selected non-render governance, alignment,
  build-policy, build-gate, attestation, producer, and C++ authority tests passed
  on both. A disposable clean snapshot then configured strict qualification
  with both compilers, emitted and validated all Slang kernels, and compiled all
  seven gate-bound executables plus the 13-product volume. That execution
  exposed a real trust-boundary defect: the application and render test binaries
  were not beside the candidate, so strict lookup correctly ignored their
  development resource override. Both suites now run from the exact candidate
  volume; governance rejects either relocation or a restored override. The nine
  non-render controls passed again on both qualification builds. Static
  governance reports 732 governed GoogleTests and 23 dimensions plus 26
  capability contracts. Development registers 741 Mandatory tests; strict
  qualification registers 744 because it additionally forbids GoogleTest skips
  and requires the Vulkan dispatch and P1 burn-in controls. Those three strict
  runtime controls were compiled but not executed in this no-render pass. Native CI now
  pins every Action and fetched dependency by commit, verifies the directly
  downloaded Slang and SwiftShader archives by SHA-256, and preflights Slang,
  SPIR-V validation, and a dispatchable SwiftShader/MoltenVK device before
  qualification. Repository governance rejects a regression to movable or
  unchecked inputs. No rendering test or full Mandatory runner was executed.
  This local mirror has a synthetic audit revision and provides implementation
  proof only; the external qualification branches and all eight admissible
  clean-revision domain attestations remain unexecuted.

- 2026-08-30 spherical Kottler correction: the advertised CPU cosmological
  sector now distinguishes the exact smaller black-hole and larger
  cosmological roots, treats pure de Sitter as horizonless for capture, and
  rejects nonpositive Lambda plus Nariai/super-Nariai public scenes. GCC 14 and
  Clang 21 each compiled the application and all six GoogleTest executables
  under warnings-as-errors, including the render executable without executing
  it; both passed all 368 core tests and the selected config boundaries. Live
  CTest inventories contain exactly 746 Mandatory registrations (737 governed
  GoogleTests plus nine operational CTests), and both inventories satisfy the
  10-criterion, 23-dimension, 28-capability model. Runtime `info capabilities`
  independently reports all 28 contracts. No rendering test was executed.
- 2026-08-30 mass-scale trace-domain correction: the CPU and Vulkan paths now
  share M-scaled step bounds and an escape surface that encloses the observer.
  The M=1 production defaults remain byte-compatible by construction, while
  the public 0.1 and 100 mass endpoints plus a 1000M observer are pinned by a
  non-rendering authority test. The live model now contains 24 dimensions and
  the source estate contains 738 governed GoogleTests (747 development CTests
  including the nine operational controls). No rendering test was executed.
- 2026-08-30 coordinate-unit and metric-parameter correction: public observer
  distance remains the established coordinate radius rather than being
  reinterpreted as a multiplier of mass; mass-bearing displays now derive
  `r/M`, while EXR metadata no longer invents a solar-mass scale. The registry
  also makes parameter applicability explicit: Minkowski, de Sitter,
  Morris-Thorne, and Alcubierre require zero mass; non-default throat/topology
  and warp/bubble fields are accepted only by their owning metric.
  Identity-aware omitted-value defaults and explicit irrelevant input fail
  closed at both configuration and typed-session boundaries. A non-render
  metadata witness brings the current estate to 739 governed source GoogleTests
  and 748 development CTests including nine operational controls. No rendering
  test was executed.
- 2026-08-30 feature-ownership correction: lens, disk-temperature, Doppler,
  diagnostic-colour, bloom, volumetric, temporal, film, and typed point-star
  parameters now require their consuming identity or feature at both public
  configuration and typed-session boundaries. CLI parameter flags select the
  owner they configure, and `MakeSessionConfig` can no longer be used to bypass
  schema validation. Shared defaults remove the former app/session bloom and
  volume drift. GCC 14 and Clang 21 each compile the complete explicit product
  and six-test-executable topology; the focused pure config/projection controls
  pass on both. The estate now contains 743 governed source GoogleTests and 752
  development CTests. During initial witness selection the then-combined
  `ConfigurationConversionPreservesObserverAndDiskControls` test unexpectedly
  executed its embedded 8x8 in-memory CPU preview once per compiler; that
  preview was then split under the explicit rendering name
  `InMemoryPreviewCompletesWithoutWritingOutput` and was not executed again; the
  later ownership audit retired that rendering witness in favour of a pure
  configure-time rejection for an inactive output path.
- 2026-08-30 exotic natural-scale correction: Morris-Thorne observer and trace
  bounds now scale with `b0`; Alcubierre uses `max(R,1/sigma)` for scene extent
  and `min(R,1/sigma)` for the smallest integration feature. The public and
  typed boundaries require an exterior observer and reject `sigma*R` outside
  `[0.1,100]`, while the device step schedule no longer hard-codes one
  coordinate unit. CLI, registry, schema comments, and session diagnostics now
  identify `sigma` as an inverse wall length. GCC 14 and Clang 21 each compiled
  the complete explicit product and six-test-executable topology plus every
  Slang precision/portability output; the three pure scale/config/CLI probes and
  all 24 registry/Alcubierre mathematical probes pass on both. The estate now
  contains 744 governed source GoogleTests, 753 development CTests, 24 operating
  dimensions, and 29 capability contracts. No rendering test was executed.
- 2026-08-30 cosmological causal-patch correction: de Sitter and spherical
  Kottler observers now intersect the ordinary natural-scale interval with
  `r <= 0.99*r_c`. One scale-safe root authority is shared by metric capture,
  public and typed validation, and CPU trace construction; the fp32 outer
  boundary rounds inward from `r_c`, and the accepted crossing segment is
  localised back onto it before the directional sky is sampled. GCC 14 and
  Clang 21 each compiled the executable and all six test binaries under
  warnings-as-errors, including the render test executable without executing
  it. The 11 pure registry/Kottler probes and two pure config/trace-domain
  probes pass on both. The estate now contains 745 governed source GoogleTests,
  754 development CTests, 24 operating dimensions, and 29 capability contracts.
  No rendering test was executed.
- 2026-08-30 interactive-projection correction: progressive sessions now
  preserve template-owned tile, thread, parallelism, bloom, and exposure
  policy instead of replacing it with viewer literals. The same projection
  preserves requested volumetric state and rejects it without its disk owner,
  passes `SessionConfigIssue` before initialization, and keeps orbit azimuth in
  its periodic typed domain. GCC 14 and Clang 21 each compiled the complete
  product and six-test-executable topology, including the render test
  executable without executing it; both passed the two pure viewer
  projection/input controls and the exact nine non-render authority CTests.
  The estate remains 745 governed source GoogleTests and 754 development
  CTests. No rendering test was executed.
- 2026-08-31 typed-boundary and live-backend alignment delta: lower metric,
  camera, disk, turbulence, RK45, oracle, post-process, film-finish, point-star,
  pixel-sampling, tracer, and display-buffer APIs now reject malformed requests
  without clamping or fallback. The indexed point catalogue owns the exact
  validated snapshot used by CPU and Vulkan; fisheye masking is explicit;
  absent film overrides inherit their preset; CPU-only scheduling and inactive
  output controls have typed owners; and Vulkan initialization no longer builds
  unused CPU scene state or requires the texture asset in point-catalogue mode.
  CPU and Slang share the 20,000-attempt ceiling, disk-segment/terminal-event
  precedence, reconstructed turbulence coordinates, and safe captured-endpoint
  handling. The explicit device ABI is dense at 65 consumed float slots with no
  future-corona padding. GCC 14 and Clang 21 each compiled the application, all
  six test executables, and fp32/fp32-comp/fp64 SPIR-V plus CUDA and Metal
  emissions under contracts and warnings-as-errors; the render executable was
  compiled but not run. On each compiler, 20 focused core boundary cases, three
  oracle boundary cases, and eight app/configuration/projection cases passed.
  Source, label, live-inventory, semantic-policy, alignment, and strict-build
  governance accept 753 governed source GoogleTests and 762 development CTests,
  with 24 operating dimensions and 29 capability contracts. No rendering test
  was executed.
- 2026-08-31 Boyer-Lindquist integration-domain correction: polarised null and
  screen initial-data construction now returns absence for invalid, impossible,
  or ambiguous requests, including an exact linear solve on the Kerr
  stationary-limit surface. Every polarised RK stage uses the metric's finite
  off-axis chart predicate and preserves the last valid state on exit. The
  canonical integrator no longer clamps `theta`, no longer calls every chart
  exit a horizon hit, and reports an impossible null projection separately
  instead of accepting the unchanged state. GCC 14 and Clang 21 each pass both
  new non-render witnesses and all 129 oracle tests. Source, label, semantic
  policy, repository, and operating-model governance accept 765 governed source
  GoogleTests, 24 operating dimensions, and 29 capability contracts. No
  rendering test was executed.
- GCC 14 strict required-Vulkan profile: the pre-alignment snapshot passed the complete
  907/907 estate on the Radeon 780M through Dozen in 210.36 seconds. A separate
  required-only run passed 697/697 on the Radeon in 281.34 seconds and the same
  697/697 on isolated Lavapipe in 119.48 seconds. The two physical JUnit files
  contain exactly 907 and 697 cases respectively, each with zero failures, zero
  errors, and zero skips; both physical runs include the required dispatch and
  repeated CPU/Vulkan P1 burn-in. Their timings are validation records, not
  performance comparisons; the later full run reused warmed driver state.
- GCC 14 Release portable profile (pre-alignment snapshot): 904/904 passed through an explicitly isolated
  Lavapipe ICD; 694 Mandatory, 160 Operational, 47 Performance, and 26 Stability
  tests passed. The concurrent run completed in approximately 39 seconds; that
  contended validation timing is not a performance result.
- Clang 21 Release with warnings as errors (pre-alignment snapshot): 904/904 passed through the same
  isolated Lavapipe route; 694 Mandatory, 160 Operational, 47 Performance, and
  26 Stability tests passed. Its concurrent run likewise completed in
  approximately 39 seconds and is not a performance result.
- GCC 14 ASan/UBSan/LSan (pre-alignment snapshot): the complete 904-test estate passed in 1578.53
  seconds with no Sirius defect and the one
  documented, printed external Vulkan-loader/driver suppression. Vulkan tests
  were pinned to Lavapipe so the physical Radeon render remained isolated.
- Physical current-worktree P3/P5 preflight: the exact 1920x1080 scene completed
  in 1575.7 seconds (SHA-256 `72d079f7dfa1f637e15e4773114fa91ed1256d093d4d96f4fa117f82824eaa8c`)
  and the exact 5616x4096 scene completed in 17105.1 seconds (SHA-256
  `094ae5b1379d39d9c6dfdddc7199f9d14d09f50c674505decd0bb364b0d768e6`).
  Both decoded at the exact dimensions and passed the sparse dynamic-range,
  resolved-source, row-distribution, and two-dimensional-structure gate. They
  remain non-admissible until reproduced by the clean-revision runbook.
- Native viewer current-worktree preflight: a real XWayland/GLFW window used the
  Radeon Vulkan backend, published a 256x192 progressive frame, and recorded
  host-delivered keyboard press/release, pointer drag, scroll, and Escape
  callbacks. The transcript and viewer log were hashed, but remain
  non-admissible until reproduced by the clean-revision viewer runbook.
- Slang `trace.slang`: SPIR-V validated for fp32, compensated fp32, and fp64;
  CUDA and Metal emission gates pass.
- Install relocation, hostile-working-directory initialisation, and resource
  failure injection passed in the estates above. A viewer-disabled GCC build
  compiled and installed, reported the capability absent, remained CPU/Vulkan
  ready, and rejected `view`. Reduced real operator runs produced 2/2 smoke and
  13/13 demonstration artefacts. The installed operating model reports all ten
  P/E criteria, 24 review dimensions, and 29 capability contracts. The final
  clean-revision physical record must bind the selected 780M/Dozen identity,
  three precision rungs, and exact decoded 1920x1080 and 5616x4096 P3/P5
  outputs under the 2048 MiB cap; the next native viewer record must bind a
  physical-Radeon progressive Vulkan frame and hashed WSLg/XWayland callback
  evidence. Attestation negative controls, shell syntax, JSON
  governance, label freshness, formatting, and diff whitespace were checked
  separately.

## 8. Scope note

The system in scope is the complete Sirius repository and its build, install,
CLI, operator-script, CPU, Vulkan, viewer, physics, and output surfaces. External
drivers, physical devices, and native operating systems remain distinct
attestation domains and are never inferred from source or Lavapipe evidence.
