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
external witness is admissible only from a strict, non-packageable qualification
build: its copied candidate executable, qualification gate, gate-generated
JUnit/log, independently rerun JUnit/live inventory, and configure-time
alignment receipt form one cross-checked byte and test-identity chain.
Development artifacts cannot enter that chain, while qualification may retain
pending domains so evidence production does not circularly require the release
it is intended to enable. CMake
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
| Metric identity | `core/metrics/registry.h`, closed `MetricId` enum | Carried. The GPU kernel's private `MetricType` enum (which had grown Gödel and Taub-NUT dead branches and a `Custom` escape hatch) is deleted; the kernel-side enum and dispatch tables are generated from the registry at build time, so a second catalogue can no longer exist. |
| Tolerances and constants | `core/constants.h` | Carried; the kernel-side constant block is generated from it. |
| ISCO | `AccretionDisk::ComputeIsco` | Carried. |
| Capture surfaces | `IMetric::InsideCaptureSurface`, analytic Kerr-Schild inverses | Carried. |
| Test labels | `scripts/generate-ctest-labels.py` | Carried. |
| Transfer encodes | each writer applies its own, once; EXR stays linear | Carried. |
| Device physics | one Slang kernel set under `src/sirius/kernels/` | New. Replaces the 6,561-line CUDA monolith that inlined its own drifted copies of metric and disk physics. |
| Dispatch duration | `render/dispatch_governor.h` (`BandController`) | New (2026-07-28, physical-GPU validation). Bounds how long any single compute dispatch holds the device, the constraint OS GPU watchdogs enforce; the memory governor remains the sole authority for residency. |

P2996 reflection is not available and is not simulated. The metric registry is
the C++ identity authority; configuration has explicit serializers and strict
shape validation; Slang constants and dispatch parameters have explicit
packing plus CPU/kernel parity witnesses. The compile-time language capability
authority reports this as `explicit-schema`, never as native reflection.

## 3. Kernel strategy

All device physics is written once, in Slang, as modules with interfaces per concern (metric family, integrator, disk emission, beam propagation). Slang was chosen over three alternatives: raw GLSL lacks the module and generics system a physics library needs; SYCL single-source C++ has no Metal path and weak Windows support outside oneAPI; WGSL has no 64-bit floating point at all. Slang compiles the same source to SPIR-V for Vulkan today, and to CUDA, Metal, and HLSL when a native adapter earns its keep by profiling; the compiler is Khronos-hosted, Apache-licensed, and ships prebuilt for all three desktop platforms.

Two live integrators exist by design, and the difference is methodological rather than accidental. The kernel integrates Cartesian geodesics with fixed-form RK4 and compensated or fp64 state accumulation on the corresponding precision rungs; the CPU reference tracer uses adaptive Dormand-Prince RK45. A fixed-step Yoshida composition of implicit-midpoint maps is deliberately confined to the double-precision oracle stack; its canonical two-form is tested directly, while state-dependent step selection and optional null projection are not misreported as symplectic operations. Two independent live methods agreeing within stated tolerance is stronger evidence than one method executed twice, so backend parity gates are statistical (per-pixel relative radiance bounds, conserved-quantity drift bounds) rather than bitwise. The reference tapes in `renders/` remain CPU-produced.

Precision inside kernels follows the ladder in section 6. Full-ray Mandatory
diagnostics independently measure energy, axial angular momentum, Carter Q,
and the null residual on the CPU RK45 and Vulkan Cartesian-RK4 paths. The CPU
integrator additionally re-normalises a ray when its null residual crosses the
owned threshold; the Vulkan render loop remains the fixed RK4 method being
measured and is not silently altered by its diagnostic.

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

The session orchestrator, tile scheduler, progress tracking, and viewer port forward with their roles intact. Tiles remain the unit of work everywhere: the CPU tracer threads over tiles, the Vulkan backend dispatches a compute grid per tile batch, and the viewer's progressive refinement consumes tiles in the existing spiral order. Both live paths propagate two deviation vectors, extract the singular axes and output-plane orientation, and feed that literal ellipse to the indexed point-star filter. The kernel constructs the Riemann tensor once per step and contracts it with both vectors. The double-precision beam integrator remains off the render path as an oracle; its covariant Jacobi state shares one RK4 tableau with the central ray and is gated against the exact radial and circular Schwarzschild null-congruence solutions to one part in \(10^6\).

The CPU polarisation path carries two observer screen vectors through the same
accepted geodesic steps, reconditions them only to remove integrator drift, and
projects the thin-disk thermal scattering state into that transported basis.
The live Kerr-Schild implementation conserves the Walker-Penrose constant and
agrees across charts with the independent Boyer-Lindquist oracle. Stokes
crossings are accumulated before false-colour film visualisation. Vulkan,
volumetric transfer, and temporal disk blur do not cross this capability
boundary.

## 6. Memory governor and precision ladder

The governor exists because the 780M-class target has a 2 GB budget that a naive full-frame HDR pipeline exhausts (a 5616 by 4096 IMAX frame at RGBA32F is 368 MB per buffer before ray state, which at 96 bytes per ray for position, momentum, deviation vectors, and accumulators is another 2.2 GB full-frame). The design bounds device residency by construction rather than by hope.

- At startup the backend reports usable budget (`VK_EXT_memory_budget` where present, else a conservative fraction of device-local memory, else a user cap). The governor derives the largest tile whose complete working set (ray state, accumulation, readback staging) fits a fixed fraction of that budget, and the tile scheduler simply receives a smaller tile size on smaller devices. Full-frame buffers live host-side only.
- The precision ladder has three rungs, selected per device and recorded in the render metadata: fp64 kernels where the device offers usable `shaderFloat64`; fp32 with Kahan compensation on Cartesian position and velocity state accumulation elsewhere; and plain fp32 for preview quality. Full-ray diagnostics derive Hamiltonian/null residual, E, L_z, and Carter-Q drift from those states. The oracle's job is to quantify what each rung costs: per-rung error against the double-precision reference is a recorded test artefact, not a guess.
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
| Sirius.Core/Metric/PHSM001A.h | core/metrics/astrophysical_scaling.h |
| Sirius.Core/Geodesic/PHGD001A.{h,cpp} | core/geodesic_integrator.{h,cpp} |
| Sirius.Core/Camera/CMBS001A.h | core/camera.h |
| Sirius.Core/Disk/PHAD000A.h | core/disk/disk_model.h |
| Sirius.Core/Disk/PHAD001A.h | core/disk/novikov_thorne_disk.h |
| Sirius.Core/Disk/PHCR001A.h | core/disk/corona.h |
| Sirius.Core/Disk/PHTR001A.h | core/disk/turbulence.h |
| Sirius.Core/Environment/PHSF001A.h | core/starfield.h |
| Sirius.Core/Jet/PHJT001A.h | core/jet.h |
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
| Sirius.Render/Shader/RDSD00*.{vert,frag} | app/viewer/shaders/*.{vert,frag} |
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

Resource-consuming application and render test executables run beside the strict candidate. They inspect the same staged resources and receipt lifecycle, never an environment-selected substitute; repository governance pins that topology.

Pull-request integration is a deliberately non-promoting compile boundary. It configures the strict qualification topology, compiles the product and all six test executables, and runs only nine named non-render governance/authority controls. It proves that no Mandatory-gate receipt is created and cannot upload an attestation. The complete runtime estate remains mandatory for qualification and release on protected-branch or external-domain producers; repository governance checks both sides of this separation.

Six gate families. Unit and property gates cover textbook values, exact identities, registry round-trips, live-path conservation, and determinism. Oracle gates compare the live path with the independent double-precision Boyer-Lindquist stack. Backend gates compare CPU and Vulkan statistics and enforce the Vulkan capability boundary. Operational gates install to a staged volume, relocate it, initialise from a hostile working directory, render, remove a mandatory resource to prove fail-closed behaviour, require a real software-Vulkan dispatch, and repeat the P1 near-extremal classifier as burn-in. Attestation-admission gates keep physical and native evidence revision- and qualification-candidate-bound and make complete admission a configure/build/compile/runtime invariant for releases. Qualification and release build gates remove their prior runtime receipt, then emit a deterministic zero-skip JUnit receipt bound to live CTest registration, every test executable, and every strict product; pre-gate strict-mode install/readiness are negative controls, and install rejects an absent, stale, or product-mismatched final receipt. Only after CTest succeeds is the new receipt staged beside the build-tree executable. Qualification cannot package or report top-level readiness, but it locks resource discovery to its executable volume and is the sole mode accepted by external-evidence admission. Release separately regenerates its own gate for the final product; release install/package creation, readiness, and render/view initialisation parse the installed receipt and independently rehash the running executable plus the complete installed product set. Qualification and release resource discovery ignore the development `SIRIUS_RESOURCE_DIR` override, preventing an environment-selected tree from substituting for that volume. Development initialisation remains receipt-optional and override-capable for local diagnosis, but its artifacts cannot be promoted to evidence. The sanitizer profile runs the Mandatory estate under ASan, UBSan, and LSan with one named software-Vulkan process-lifetime suppression. Historical image tapes remain useful forensic material, but no current executable byte-identity gate consumes them; `docs/ADVERSARIAL_REVIEW.md` records that claim as unverified rather than treating files alone as evidence.
