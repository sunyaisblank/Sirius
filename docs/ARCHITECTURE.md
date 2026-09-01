# Sirius Architecture

This document is the design authority for the rebuilt system. It records the layer structure, the single-authority seams, the kernel and backend strategy, the memory and precision design that reaches 2 GB integrated GPUs, the dependency decisions, and the complete mapping from the retired codename files to their successors. The specification (`docs/SPECIFICATION.md`) says what must be true; this document says how the system makes it true.

## 1. Layers

The system keeps the July layering, orchestration above primitives, dependencies pointing downward only, and adds one layer the old system lacked: a device abstraction that makes every GPU vendor an adapter rather than an architecture.

```
sirius::app        CLI, external configuration, typed session projection, viewer
sirius::render     session lifecycle, tiles, memory governor, output writers
sirius::backend    device abstraction: CPU tracer, Vulkan compute, future adapters
sirius::kernels    Slang single-source device physics (+ generated tables)
sirius::core       metrics, geodesics, tensors, disk, spectral, camera, constants
sirius::oracle     double-precision Boyer-Lindquist validation stack (test-only)
sirius::base       contracts, error types, threading and SIMD seams
```

`core` and `oracle` compile and test with no GPU, no window, and no Slang toolchain present. `kernels` is source consumed by `backend` at build time; nothing above `backend` knows which device executes a render.

Each compiled implementation has one CMake owner. `sirius_backend_cpu` always
owns the CPU tracer; `sirius_backend` owns only the optional Vulkan device; and
`sirius_render` consumes those targets without compiling their source files.
Application JSON structs stay in `app/config`, retain the public camel-case JSON
wire keys for compatibility, and project once into closed render-domain enums
and snake-case values through `session_config_adapter`. Exceptions do not cross
that boundary: projection returns `std::expected<SessionConfig, Error>`.

The interactive viewer follows the same ownership rule. Its public header holds
declarations and state; `viewer/interactive_viewer.cpp` owns the threaded
refinement loop. CLI and test translation units therefore do not compile the
viewer implementation through an oversized header.

The source repository contains only these durable top-level concerns:

| Path | Ownership |
|---|---|
| `src/sirius/` | First-party C++ and Slang, split by the layers above |
| `tests/` | Unit, oracle, backend, render, application, and operational gates |
| `scripts/` | Test-label, operating-model, structure, attestation governance, and native evidence producers |
| `cmake/` | Toolchain, dependency, embedded-model, and install-volume logic |
| `assets/` | Required runtime assets installed with the executable |
| `lib/` | Only the four live header/source vendors enumerated in section 8 |
| `docs/` | Specification, architecture, style, current review, and history |

`scripts/verify-repository-structure.py` makes this more than a diagram: normal
builds fail if a translation unit gains two owners, a lower layer reaches up, a
retired vendor/configuration path returns, or boundary identifiers regress.

External evidence has the same one-authority shape. Individual bundles are
validated by `verify-attestation.py`; `verify-alignment.py` admits exactly one
witness for each model-derived required domain at one clean revision and emits
the sole alignment receipt, including the exact operating-model digest. An
individual bundle's gated operating-model product must first byte-match that
independently hashed source model; agreement among a substituted bundle's own
receipt, gate, and product is insufficient. An external witness is admissible
only from a strict, non-packageable qualification
build: its copied candidate executable, qualification gate, gate-generated
JUnit/log, independently rerun JUnit/live inventory, and configure-time
alignment receipt form one cross-checked byte and test-identity chain.
Development artifacts cannot enter that chain, while qualification may retain
pending domains so evidence production does not circularly require the release
it is intended to enable. Before that expensive authority executes,
`attestation_preflight.py` queries only the exact candidate's system and
readiness interfaces, cross-checks its compiled receipt/resources against the
clean worktree and physical host/device route, and emits a deliberately
non-promoting plan outside the source tree. It cannot replace any runtime,
image, or input witness. CMake
rechecks and embeds that receipt, while
`app/alignment_authority` compares the installed copy with the compiled bytes
before operational initialisation and exposes the exact admitted/pending/required
domain partition through readiness. Thus upstream evidence cannot drift from the
downstream runtime claim or collapse into an unactionable count.

External test-estate claims follow the same rule: CTest's JSON registration is
hashed beside JUnit, the enabled-name sets must be identical, and a local
source parser plus the external verifier independently enforce the model's
700-case floor so a trivially self-consistent subset cannot become an authority.
Runtime records additionally bind the clean
configure receipt to the JUnit case that compared it with the compiled binary.

## 2. Single authorities

The July remediation's central lesson was that every duplicated authority eventually drifts. The rebuild carries each existing authority forward and closes the one hole the inventory found.

| Concern | Authority | Change from July |
|---|---|---|
| Metric identity and natural scale | `core/metrics/registry.h`, closed `MetricId` enum | Carried and strengthened. The GPU kernel's former physics catalogue (including dead Gödel, Taub-NUT, and `Custom` branches) is deleted. A compact device-family ABI remains, but `MapMetric` exhaustively translates every registry identity to a represented family or an explicit decline, and parity probes pin those ABI values. Observer and trace-domain scaling share the same `M`, `b0`, or Alcubierre `max/min(R,1/sigma)` authority; the same registry also owns the scale-safe positive-$\Lambda$ Kottler roots, the `0.99*r_c` observer bound, and the cosmological trace clamp. |
| Kerr--Schild field-equation boundary | live `core/metrics/kerr_schild_family.h`; independent `tests/core/metric_curvature_oracle.h` and `kerr_schild_field_equation_test.cpp` | Strengthened. A fourth-order connection-difference Ricci oracle uses a generic inversion of the sampled live covariant metric, while an independently differentiated Kerr--Newman potential supplies $F_{\mu\nu}$ and $T^{\rm EM}_{\mu\nu}$. Vacuum, Einstein--Maxwell, source-free Maxwell, and spherical $\Lambda$ equations are Mandatory across signs, limits, and scales. Structural negative controls reject production-$H$ reuse, specialised-inverse coupling, oracle-only metric substitution, or removal of the Maxwell/cosmological branches. |
| Optional-feature parameter ownership | `core/camera.h`, `core/disk/disk_defaults.h`, `core/feature_defaults.h`, `app/config/config_schema.h`, and `render::SessionConfigIssue` | Lens, disk-emission, bloom, volumetric, temporal, point-catalogue, output, backend-specific CPU scheduling, and film-finish parameters have one explicit owner and neutral inactive state. CLI parameter flags select their owning feature; absent film overrides inherit the selected preset, while JSON and typed-session callers that retain non-default parameters under a different or disabled owner decline before work. Interactive projection preserves controls owned by the selected backend, changes only live viewer state, neutralises CPU-only scheduling controls for Vulkan, and preflights the resulting typed session before initialization. |
| Host tolerances and constants | `core/constants.h` | Carried. Device-only fp32/fp64 controller tolerances have different accumulation semantics and therefore live once in `gr_integrator.slang`; production and the live P1 parity probe call the same policy functions, while the probe binds that accepted envelope to CPU/oracle invariants. |
| Host blackbody laws | `core/spectral/blackbody_laws.h`, consuming `core/constants.h` | Strengthened. Host colour integration and binned spectral radiance delegate to one Planck spectral-radiance-per-metre function; only the bin facade converts density to per nanometre. Wien peak wavelength and Stefan-Boltzmann hemispheric radiant exitance are unit-explicit fallible quantities. Repository governance rejects a second host Planck implementation, copied Wien constant, detached consumer, or the former exitance-as-radiance name; independent CODATA, quadrature, numerical-maximum, and Slang compute witnesses remain comparators. |
| Celestial tangent orientation | `core/celestial_tangent_basis.h`, mirrored by `gr_deviation.slang` and live parity-gated | The observer Sachs screen, terminal beam ellipse, and anisotropic point-catalogue filter select the same least-aligned coordinate axis, so an ellipse angle never crosses an implicit basis rotation. |
| Kerr circular orbits and frame transfer | `core/kerr_orbits.h`, `core/relativistic_transfer.h` | Strengthened. Metric identity, Page-Thorne disk setup, invariant frequency transfer, and CPU tracer emitter motion share one exact signed-spin ISCO/angular-velocity authority. Thin and volumetric frame frequencies use $E$, $L_z$, and the actual Boyer-Lindquist latitude; the diagnostic branch is a locally non-rotating ZAMO, never the radially moving Kerr-Schild slicing normal. The Boyer-Lindquist oracle and Slang implementation remain independent comparators, with cross-consumer, arbitrary-latitude, structural, and device-parity gates preventing reuse or drift. |
| Capture, topology, and causal-patch surfaces | `IMetric::InsideCaptureSurface`, `IMetric::IsotropicEllisThroatRadius`, analytic Kerr-Schild inverses, registry-owned exact Kottler black-hole/cosmological roots, `core/trace_boundary.h`, `kernels/gr_trace_event.slang`, `kernels/gr_ellis_topology.slang` | Strengthened. Kottler's smaller root owns capture and its larger root owns observer admission and the outer directional boundary. The horizonless Ellis metric exposes its regular isotropic throat without calling it capture. Its live host/device surface admits only the asymptotically normalised zero-tidal member; an independent connection-difference oracle recovers the exact Ricci scalar and radial null-energy-condition violation on both ends. `OneSheetCapture` uses complete accepted-segment root isolation, including hidden and tangent contact, and publishes a distinct dark throat event. `TwoSheet` crosses into $0<\rho<b_0/2$, reaches the inversion-matched second infinity, records the terminal sheet, and transforms its sky direction by the exact inversion Jacobian. Production rendering and the compute parity probe call the same device two-sheet trajectory authority. |
| Alcubierre geometry and stress-energy boundary | `core/metrics/registry.h`, `core/metrics/warp_drive_family.h`, mirrored by `gr_metrics.slang` and fallible `gr_integrator.slang` stages | Strengthened. One complete host domain governs finite velocity/radius/inverse-wall parameters and the jointly resolved $0.1\le\sigma R\le100$ wall; direct CPU setters cannot bypass it, and the Slang mirror rejects the same domain before metric, connection, or any RK stage can turn invalid data into flat motion. The exact unit-determinant $t$--$x$ block owns one analytic inverse on each backend. A shared test-only fourth-order connection-difference framework independently reconstructs the Einstein tensor and pins the negative Eulerian wall density, while repository negative controls reject a future domain bypass, determinant repair, or infallible stage. |
| Test labels | `scripts/generate-ctest-labels.py` | Carried. |
| CIE 1931 observer fit | `core/cie1931_observer.h`, mirrored by `gr_disk.slang` and live parity-gated | Strengthened. One explicitly approximate Wyman-Sloan-Shirley 2-degree observer fit serves both host spectral facades at their actual wavelengths; it is bounded across 380--780 nm against checksum-identified official CIE 018:2019 samples, and out-of-band input fails closed. The device disk path owns one fp32 mirror, direct compute parity pins it, and repository governance rejects the former wavelength-mislabeled table, detached facades, copied Gaussian lobes, or a missing parity route. |
| XYZ-D65 to linear sRGB | `core/xyz_srgb.h`, mirrored by `gr_disk.slang` and live parity-gated | Strengthened. The explicitly normalised host blackbody-colour path consumes one exact rational transform derived from the declared sRGB primaries and D65 white point. The binned physical-radiance facade stops at unnormalised XYZ. The device disk path owns one fp32 mirror, direct compute parity pins it independently, extended-gamut values remain linear instead of being hidden by this stage, and repository governance rejects rounded legacy coefficients, a second exact coefficient copy, or a detached consumer. |
| Binned spectral-to-display boundary | deliberately absent | Strengthened by removal. `SpectralRadiance::ToXyz()` retains physical radiometric scale, so a primary transform, clipping, and transfer curve alone cannot define an sRGB preview. The unused `ToSrgb()` facade and its arbitrary-scale range test are retired; repository governance rejects both that API and direct display encoding inside the physical-radiance facade until exposure or normalisation and display intent are explicitly represented. |
| Spectral-to-ACES boundary | deliberately absent | Strengthened by removal. `SpectralRadiance` carries physical radiance density and its integrated tristimulus values inherit that radiometric scale; an AP0 primary matrix alone cannot establish ACES2065-1 relative exposure. The unused `ToAces()` facade and detached matrix are retired, while repository governance rejects either returning without an explicitly represented exposure/encoding authority. |
| Fitted display tone curve | `core/postprocess.h` | Strengthened. `ACESFit` names the five-coefficient Narkowicz scalar fit exactly and owns its config spelling, parser, enum, coefficients, and dispatch. It makes no Academy ACES Output Transform or colour-management claim; bare `ACES` declines, while repository governance rejects a copied coefficient set, detached parser/session consumers, the former falsely broad internal name, or misleading CLI help. Source-derived samples and the positive analytic derivative pin the represented rational function. |
| Display transfer encode | `core/srgb_transfer.h` for host writers; `RDSD003A.frag` as the sole live-view mirror | Strengthened. One host IEC 61966-2-1 curve and clipped round-to-nearest 8-bit quantiser serve explicitly normalised blackbody previews, explicit byte views, PNG, and PPM. The live viewer consumes the already tone-mapped, graded, bounded display-linear frame and mirrors only that IEC curve in its installed fragment shader; hardware framebuffer sRGB is disabled to prevent a second encode. PNG/PPM still own their file encoding, EXR remains linear, non-finite byte input declines, and repository governance rejects a second production breakpoint, viewer tone map, approximate gamma, detached consumer, or direct encoding of binned physical radiance. Independent numerical, shader-source, and output decode/round-trip tests remain the comparators. |
| Device physics | one Slang kernel set under `src/sirius/kernels/` | New. Replaces the 6,561-line CUDA monolith that inlined its own drifted copies of metric and disk physics. |
| Dispatch duration | `render/dispatch_governor.h` (`BandController`) | New (2026-07-28, physical-GPU validation). Bounds how long any single compute dispatch holds the device, the constraint OS GPU watchdogs enforce; the memory governor remains the sole authority for residency. |

P2996 reflection is not available and is not simulated. The metric registry is
the C++ identity authority; configuration has explicit serializers and strict
shape validation; Slang constants and dispatch parameters have explicit
packing plus CPU/kernel parity witnesses. The compile-time language capability
authority reports this as `explicit-schema`, never as native reflection.

## 3. Kernel strategy

All device physics is written once, in Slang, as modules with interfaces per concern (metric family, integrator, disk emission, beam propagation). Slang was chosen over three alternatives: raw GLSL lacks the module and generics system a physics library needs; SYCL single-source C++ has no Metal path and weak Windows support outside oneAPI; WGSL has no 64-bit floating point at all. Slang compiles the same source to SPIR-V for Vulkan today, and to CUDA, Metal, and HLSL when a native adapter earns its keep by profiling; the compiler is Khronos-hosted, Apache-licensed, and ships prebuilt for all three desktop platforms.

Two live integrators exist by design, and the difference is methodological rather than accidental. The kernel integrates Cartesian geodesics with adaptive RK4 step doubling and compensated or fp64 state accumulation on the corresponding precision rungs; the CPU reference tracer uses adaptive Dormand-Prince RK45. Both methods reject a candidate step when either its scale-aware truncation estimate or relative null residual exceeds the owned tolerance, preserving the prior state for retry, and both share the same 20,000-attempt render ceiling. The kernel additionally applies an accepted-state constraint method: after the unprojected RK4 candidate passes those local bounds, a stable quadratic/linear solve changes only its temporal tangent component, selecting the root nearest the incoming component so the light-cone branch is continuous; an unrepresented projection rejects the attempt. This prevents independently accumulated roundoff from becoming an irreducible global null defect while leaving the adaptive controller responsible for local integration error. One host authority scales both paths' step bounds and far boundary with the scene's natural geometry and expands that boundary to enclose the observer: `M` for black-hole sectors, `b0` for Morris-Thorne, `max(R,1/sigma)` for the Alcubierre extent, and `min(R,1/sigma)` for its narrowest integration feature. Minkowski and de Sitter retain the conservative unit scale. This prevents the kernel's pre-step escape test from accepting an observer outside a fixed sphere and keeps exotic metrics scale-covariant instead of silently reverting to unit geometry. A fixed-step Yoshida composition of implicit-midpoint maps is deliberately confined to the double-precision oracle stack; its canonical two-form is tested directly, while state-dependent step selection and optional oracle-only null projection are not misreported as symplectic operations. Two independent live methods agreeing within stated tolerance is stronger evidence than one method executed twice, so backend parity gates are statistical (per-pixel relative radiance bounds, conserved-quantity drift bounds) rather than bitwise. The reference tapes in `renders/` remain CPU-produced.

The explicit host/device trace-parameter ABI has 66 occupied slots, numbered 0
through 65. Packing, bounds, and portability builds pin that literal boundary;
unused padding is not presented as a capability.

Precision inside kernels follows the ladder in section 6. Full-ray Mandatory
diagnostics independently measure energy, axial angular momentum, Carter Q,
and the null residual on the CPU RK45 and Vulkan Cartesian-RK4 paths. Constraint
drift participates in live step acceptance on both paths; neither live method
re-normalises, projects, or otherwise changes the null branch after a step.

## 4. Backends and platforms

`sirius::backend` exposes one interface: enumerate devices, report capabilities (memory budget, fp64 support, queue families), compile or load kernels, dispatch tiles, read back results. Selection order is explicit and user-overridable, never silent.

| Backend | Reaches | Status in this engagement |
|---|---|---|
| CPU tracer | every machine, no GPU or loader required | Ported live; the reference implementation |
| Vulkan compute | AMD, Intel, NVIDIA on Windows and Linux; Apple via MoltenVK; WSL2 via Dozen or Lavapipe; CPU-only machines via Lavapipe/SwiftShader | Required and non-skipping on Lavapipe in the current Linux operational profile. A 2026-07-28 Dozen/780M run is historical evidence; the current review does not present it as a fresh physical-driver attestation. |
| Slang→CUDA native | NVIDIA, if profiling shows Vulkan leaves performance unclaimed | Emission is a standing build gate: trace.slang → trace.cu every build, entry point pinned by the KernelPortability suite; adapter deferred |
| Slang→Metal native | Apple silicon | Emission is a standing build gate: trace.slang → trace.metal every build, entry point pinned by the KernelPortability suite; adapter deferred |

WSL2 is a named platform tier, not an afterthought: no native GPU kernel driver exists there, so the backend probes in order Dozen (Vulkan over the D3D12 paravirtual device), then Lavapipe, then declines to the CPU tracer with a clear message. The OptiX and CUDA host code, the PTX loader, and the vendored OptiX headers are retired with the kernel monolith; the Radeon 780M-class target is served by the Vulkan backend through the standard AMD drivers on Windows and RADV on Linux.

On WSL, Dozen's dynamically loaded `libd3d12.so` wrapper and
`libd3d12core.so` register a pthread TLS destructor. Sirius retains both
mappings for process lifetime after and only after the Vulkan driver identifies
itself as `VK_DRIVER_ID_MESA_DOZEN`; without the lease, destroying the instance
on a render worker can unmap the destructor before the worker exits. The backend
worker-thread dispatch gate exercises this ordering, while other Vulkan drivers
never enter the WSL-specific path.

## 5. Render orchestration

The session orchestrator, progress tracking, and viewer port forward with their roles intact. On CPU, the tile scheduler threads over spiral-ordered operator tiles. Vulkan does not initialise that unused CPU scheduler or its duplicate scene objects: the backend derives device-budget tiles, then subdivides them into watchdog-governed row bands, while the viewer consumes complete progressive frames. Both live paths propagate two deviation vectors, extract the singular axes and output-plane orientation, and feed that literal ellipse to an index that owns the exact validated point-catalogue snapshot used to build its topology. The terminal Sachs screen and catalogue tangent plane call one least-aligned-axis basis authority on each host/device path, with a live cross-backend probe pinning the mirrored implementations; anisotropic filtering therefore interprets orientation in the basis that produced it. Thin and volumetric disk contributions are accumulated along each accepted observer-to-scene segment before terminal horizon classification; captured endpoints do not construct an invalid terminal Eulerian frame. For each accepted central-ray segment, the Jacobi RK4 tableau evaluates connection and full Riemann curvature at the accepted start, cubic-Hermite midpoint, and accepted end, sharing each stage across both deviation columns. The double-precision beam integrator remains off the render path as an oracle and is gated against the exact radial and circular Schwarzschild null-congruence solutions to one part in \(10^6\).

The CPU polarisation path carries two observer-screen vectors through the same
accepted central-ray segments, reconditions them within the local observer
screen only to remove integrator drift, and projects the circular-emitter disk
screen into that transported basis. Emission uses the flux-normalised
Chandrasekhar-Sobolev semi-infinite pure electron-scattering atmosphere; it
does not claim absorption, finite optical depth, returning radiation, or
magnetic/Faraday effects.
The live Kerr-Schild implementation conserves the Walker-Penrose constant and
agrees across charts with the independent Boyer-Lindquist oracle. That oracle
owns one finite off-axis chart domain end to end: null and screen initial-data
construction return absence for impossible or ambiguous requests, every
polarised RK stage is checked without angular clamping, and integration retains
the last represented state while distinguishing capture, escape, chart exit,
constraint failure, and exhaustion. The canonical oracle likewise validates
its chart before admitting a radial capture event, is history-free under an
explicit renormalization cadence, and reports midpoint nonconvergence or an
impossible null projection instead of accepting an unchanged non-null state. Stokes
crossings are accumulated before the bounded film-inspired display finish. Vulkan,
volumetric transfer, and temporal disk blur do not cross this capability
boundary.

The relativistic thin-disk default and physical parity authority remains the
zero-torque Page-Thorne profile integrated from the declared stable-orbit
inner edge, which must be at or outside the ISCO. CPU and Slang evaluate every
represented nonzero spin rather than substituting the Schwarzschild ISCO over
a low-spin interval. The CPU tracer constructs and invalidates its cached disk
profile from that declared edge, and both thin and volumetric Slang consumers
pass the same edge into their independent quadrature rather than recomputing an
ISCO at the shading boundary. A non-render live CPU trace witness covers a
larger truncated edge after `SetConfig`; source governance rejects either host
cache retention or a host/device ISCO substitution. `ShakuraSunyaev` is an
explicit, narrower substitution:
CPU and Slang evaluate the Newtonian zero-torque shape
$F(r)\propto r^{-3}[1-\sqrt{r_\mathrm{in}/r}]$ and define the operator
temperature at $1.5r_\mathrm{in}$. The shared normalisation and CPU/Slang
parity prevent the former bare power law from reappearing; this mode does not
represent relativistic flux or alpha-disk density, opacity, or vertical
structure.

## 6. Memory governor and precision ladder

The governor exists because the 780M-class target has a 2 GB budget that a naive full-frame HDR pipeline exhausts (a 5616 by 4096 IMAX frame at RGBA32F is 368 MB per buffer before ray state, which at 96 bytes per ray for position, momentum, deviation vectors, and accumulators is another 2.2 GB full-frame). The design bounds device residency by construction rather than by hope.

- At Vulkan startup the backend reports usable budget (`VK_EXT_memory_budget` where present, else a conservative fraction of device-local memory, else a user cap) and independently derives the largest device tile whose complete working set (ray state, accumulation, readback staging) fits a fixed fraction of that budget. CPU tile/thread controls neither govern nor initialise this path. Full-frame buffers live host-side only.
- The precision ladder has three rungs, selected per device and recorded in the render metadata: fp64 kernels where the device offers usable `shaderFloat64`; fp32 with Kahan compensation on Cartesian position and velocity state accumulation elsewhere; and plain fp32 for preview quality. Full-ray diagnostics derive Hamiltonian/null residual, E, L_z, and Carter-Q drift from those states. Separate storage-buffer probe artefacts compile the same trajectory and controller authorities for all three rungs, prove the fp64 module actually declares SPIR-V `Float64`, and require physical capture or escape plus the stated invariant envelope without shading an image. Image comparisons remain render-behaviour evidence, not the precision authority: capture boundaries are discontinuous and each precision policy owns a distinct adaptive schedule.
- Progressive refinement is the interactivity strategy on weak devices: samples per pixel grow monotonically across passes, every pass is a complete image, and the viewer stays responsive at one pass per tile budget regardless of device speed.
- The dispatch governor (`render/dispatch_governor.h`) is the memory governor's time-domain companion, added when physical-GPU validation showed the two constraints are independent: a tile can fit the memory budget and still exceed the operating system's GPU watchdog in a single dispatch (Windows TDR removes the device after ~2 s; software Vulkan has no watchdog, so only silicon exposed this). Tiles are submitted as adaptive row bands targeting a wall-time budget per dispatch (default 250 ms, `SIRIUS_DISPATCH_TARGET_MS` overrides, 0 disables). The trace kernel addresses pixels absolutely, so banding is value-transparent by construction; the controller's growth is damped and its shrink unbounded because overshooting risks device removal while undershooting only costs submission overhead.

## 7. The port: codename to descriptive mapping

The port renames every file per the style guide. The table is the traceability record from July's evidence ledger and project memory to the new tree; the A/B variant split (production versus double-precision oracle twin) survives as the `core`/`oracle` layer split.

| Old (src/…) | New (src/sirius/…) |
|---|---|
| Sirius.Core/Constants/PHCN001A.h | core/constants.h |
| Sirius.Core/Tensor/MTTN001A.{h,cpp} | core/tensor.{h,cpp} |
| Sirius.Core/Autodiff/MTDL001A.h | core/dual_number.h |
| Sirius.Core/Coordinate/PHCT002A.h | core/coordinates.h |
| Sirius.Core/Metric/PHMT000A.h | core/metrics/metric.h |
| Sirius.Core/Metric/PHMT100A.h | core/metrics/kerr_schild_family.h |
| Sirius.Core/Metric/PHMT101A.h | core/metrics/morris_thorne_family.h |
| Sirius.Core/Metric/PHMT102A.h | core/metrics/warp_drive_family.h |
| Sirius.Core/Metric/PHMT200A.h | core/metrics/registry.h |
| Sirius.Core/Metric/PHSM001A.h | retired: unused preset-only SI scaling; constants remain central in core/constants.h |
| Sirius.Core/Geodesic/PHGD001A.{h,cpp} | core/geodesic_integrator.{h,cpp} |
| Sirius.Core/Camera/CMBS001A.h | core/camera.h |
| Sirius.Core/Disk/PHAD000A.h | core/disk/disk_model.h |
| Sirius.Core/Disk/PHAD001A.h | core/disk/novikov_thorne_disk.h |
| Sirius.Core/Disk/PHTR001A.h | core/disk/turbulence.h |
| Sirius.Core/Environment/PHSF001A.h | core/starfield.h |
| Sirius.Core/Spectral/PHSP001A.h | core/spectral/blackbody.h |
| Sirius.Core/Spectral/MTSB001A.h | core/spectral/spectral_radiance.h |
| Sirius.Core/Spectral/PHSC001A.h | core/spectral/colour_modes.h |
| Sirius.Core/Polarisation/PHPL001A.h | core/polarisation/stokes.h |
| Sirius.Core/PostProcess/PPOP001A.h | core/postprocess.h |
| Sirius.Core/Metric/PHMT000B.h | oracle/metric_interface.h |
| Sirius.Core/Metric/PHMT100B.h | oracle/kerr_boyer_lindquist.h |
| Sirius.Core/Symplectic/PHSI001A.h | oracle/symplectic_integrator.h |
| Sirius.Core/Transport/MTTP001A.h | oracle/transport_types.h |
| Sirius.Render/Integration/INBI001A.h | oracle/beam_integrator.h |
| Sirius.Render/Integration/GTRC001A.{h,cpp} | backend/cpu/geodesic_tracer.{h,cpp} |
| Sirius.Render/Acceleration/Backend/ACIB001A.h, ACBM001A.h, ACBF001A.cpp | backend/device.h, backend/device_selector.{h,cpp} |
| Sirius.Render/Acceleration/OptiX/* (RDOP/RDOX/RDPTX) | retired; physics reauthored in kernels/*.slang, host in backend/vulkan/ |
| Sirius.Render/Output/OUIB001A.h | render/image_buffer.h |
| Sirius.Render/Output/OUEW001A.{h,cpp} | render/exr_writer.{h,cpp} |
| Sirius.Render/Output/OUPN001A.{h,cpp} | render/png_writer.{h,cpp} |
| Sirius.Render/Output/RDFL001A.h | render/film_pipeline.h |
| Sirius.Render/Pipeline/RDRT001A.{h,cpp} | render/renderer.{h,cpp} |
| Sirius.Render/Shader/RDSD003A.{vert,frag} | app/viewer/shaders/RDSD003A.{vert,frag}; later unreachable RDSD004A/RDSD005A effects retired |
| Sirius.Infrastructure/Application/CREP001A.cpp | app/main.cpp |
| Sirius.Infrastructure/Cli/CRCL001A–006A | app/cli/{command_router,render_command,info_command,config_command,cli_output,view_command}.{h,cpp} |
| Sirius.Infrastructure/Configuration/CRCF002A, CRCF003A | app/config/{config_loader,config_schema,session_config_adapter}.{h,cpp} |
| Sirius.Infrastructure/Configuration/CRFM001A | render/film_config.h |
| Sirius.Infrastructure/Platform/CRPF001A.{h,cpp} | app/platform_paths.{h,cpp} |
| Sirius.Infrastructure/Session/SMSM001A, SNST001A, SNEV001A, SNPR001A, SNDP001A, SNTL001A, SNRS001A | render/session/{state_machine,session_states,session_events,progress_tracker,display_buffer,tile_scheduler,render_session}.h/.cpp |
| Sirius.Infrastructure/Viewer/UIVW001A.h | app/viewer/interactive_viewer.{h,cpp} |
| Sirius.Test/TS*.cpp | tests/, one file per subject, names per style guide |

New files with no predecessor: `base/contracts.h`, `base/error.h`, `base/threading.h`, `kernels/*.slang`, `backend/vulkan/*`, `render/memory_governor.{h,cpp}`, `core/polarisation/walker_penrose.h` (specification E2), and the beam propagation additions to the CPU tracer (P2).

## 8. Dependencies

The vendored source set is deliberately small and closed: GLFW and glad serve
the optional OpenGL viewer, while stb and tinyexr/miniz serve image IO. FTXUI,
nlohmann_json, and GoogleTest are fetched at configure time with immutable
commit pins. Slang is a build-time compiler found at `/opt/slang` or
`SLANG_ROOT`; Vulkan's
loader and headers are system dependencies. GLM and Dear ImGui were removed
after a source/build audit proved they had no consumer; the viewer uses its own
small camera state and direct OpenGL blit. CUDA, OptiX, PTX host machinery, and
glad's former compute duties remain retired. GLFW's disabled upstream examples,
self-tests, and demo-only dependency copies are also omitted; its runtime,
platform, Wayland, and MinGW compatibility sources remain intact. No scene graph, ECS, renderer
framework, or unused prebuilt library is carried in the repository.

## 9. Validation architecture

Application and rendering tests run beside the strict candidate and inspect the same staged resources and receipt lifecycle, never an environment-selected substitute. Their execution boundary is structural: `sirius_app_tests` owns only parsing, validation, configuration/session projection, and input-state checks; every case that can enter a render session, dispatch CPU/Vulkan work, write a rendered output, or publish a completed frame is compiled only into `sirius_render_tests` and labelled `Rendering`. Shared source governance rejects the six governed rendering identities outside `tests/render`, any app-side viewer/session start, and every new unclassified `RenderCommand` or `ViewCommand` execution. This makes an unfiltered application-suite invocation non-rendering while preserving the complete rendering-inclusive Mandatory estate for external qualification and release.

Pull-request and explicitly dispatched integration are a deliberately non-render compile boundary. Linux, Windows, and macOS each configure the strict qualification topology, compile the product and all six test executables, and run only nine named governance/authority controls. Pull requests cannot upload evidence. An explicit dispatch may issue only the Windows/macOS compilation domains through a separate non-promoting gate that binds the complete registration and compiled artifact topology while declaring runtime false. No integration path creates a Mandatory receipt; the complete runtime estate remains mandatory for physical/runtime domains and final release.

Six gate families. Unit and property gates cover textbook values, exact identities, registry round-trips, live-path conservation, and determinism. Oracle gates compare the live path with the independent double-precision Boyer-Lindquist stack. Backend gates compare CPU and Vulkan statistics and enforce the Vulkan capability boundary. Operational gates install to a staged volume, relocate it, initialise from a hostile working directory, render, remove a mandatory resource to prove fail-closed behaviour, require a real software-Vulkan dispatch, and repeat the P1 near-extremal classifier as burn-in. Attestation-admission gates keep physical and native evidence revision- and qualification-candidate-bound and make complete admission a configure/build/compile/runtime invariant for releases. Qualification and release build gates remove their prior runtime receipt, then emit a deterministic zero-skip JUnit receipt bound to live CTest registration, every test executable, and every strict product; pre-gate strict-mode install/readiness are negative controls, and install rejects an absent, stale, or product-mismatched final receipt. Only after CTest succeeds is the new receipt staged beside the build-tree executable. Qualification cannot package or report top-level readiness, but it locks resource discovery to its executable volume and is the sole mode accepted by external-evidence admission. Release separately regenerates its own gate for the final product; release install/package creation, readiness, and render/view initialisation parse the installed receipt and independently rehash the running executable plus the complete installed product set. Qualification and release resource discovery ignore the development `SIRIUS_RESOURCE_DIR` override, preventing an environment-selected tree from substituting for that volume. Development initialisation remains receipt-optional and override-capable for local diagnosis, but its artifacts cannot be promoted to evidence. The sanitizer profile runs the Mandatory estate under ASan, UBSan, and LSan with one named software-Vulkan process-lifetime suppression. Historical image tapes remain useful forensic material, but no current executable byte-identity gate consumes them; `docs/ADVERSARIAL_REVIEW.md` records that claim as unverified rather than treating files alone as evidence.
