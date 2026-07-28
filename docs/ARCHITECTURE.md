# Sirius Architecture

This document is the design authority for the rebuilt system. It records the layer structure, the single-authority seams, the kernel and backend strategy, the memory and precision design that reaches 2 GB integrated GPUs, the dependency decisions, and the complete mapping from the retired codename files to their successors. The specification (`docs/SPECIFICATION.md`) says what must be true; this document says how the system makes it true.

## 1. Layers

The system keeps the July layering, orchestration above primitives, dependencies pointing downward only, and adds one layer the old system lacked: a device abstraction that makes every GPU vendor an adapter rather than an architecture.

```
sirius::app        CLI, configuration, session orchestration, viewer
sirius::render     tracing orchestration, tiles, memory governor, output writers
sirius::backend    device abstraction: CPU tracer, Vulkan compute, future adapters
sirius::kernels    Slang single-source device physics (+ generated tables)
sirius::core       metrics, geodesics, tensors, disk, spectral, camera, constants
sirius::oracle     double-precision Boyer-Lindquist validation stack (test-only)
sirius::base       contracts, error types, threading and SIMD seams
```

`core` and `oracle` compile and test with no GPU, no window, and no Slang toolchain present. `kernels` is source consumed by `backend` at build time; nothing above `backend` knows which device executes a render.

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

The generated-table mechanism is the reflection stand-in the style guide names: a build-time script reads the registry and constants headers (their declarative sections carry structured comments) and emits both the Slang-side tables and a C++ mirror test that fails if generation is stale. When P2996 reflection lands in release compilers, the generator is replaced and nothing downstream moves.

## 3. Kernel strategy

All device physics is written once, in Slang, as modules with interfaces per concern (metric family, integrator, disk emission, beam propagation). Slang was chosen over three alternatives: raw GLSL lacks the module and generics system a physics library needs; SYCL single-source C++ has no Metal path and weak Windows support outside oneAPI; WGSL has no 64-bit floating point at all. Slang compiles the same source to SPIR-V for Vulkan today, and to CUDA, Metal, and HLSL when a native adapter earns its keep by profiling; the compiler is Khronos-hosted, Apache-licensed, and ships prebuilt for all three desktop platforms.

Two live integrators exist by design, and the difference is methodological rather than accidental. The kernel integrates with the time-transformed Yoshida symplectic family the GPU path already uses, horizon-regularised and conservation-preserving in single precision; the CPU reference tracer keeps the Dormand-Prince RK45 Hamiltonian formulation. Two independent methods agreeing within stated tolerance is stronger evidence than one method executed twice, so backend parity gates are statistical (per-pixel relative radiance bounds, conserved-quantity drift bounds) rather than bitwise. The reference tapes in `renders/` remain CPU-produced.

Precision inside kernels follows the ladder in section 6; conserved quantities (E, L_z, Carter Q, the null condition) are monitored per ray and renormalised on the schedule the constants authority sets, exactly as the CPU path does.

## 4. Backends and platforms

`sirius::backend` exposes one interface: enumerate devices, report capabilities (memory budget, fp64 support, queue families), compile or load kernels, dispatch tiles, read back results. Selection order is explicit and user-overridable, never silent.

| Backend | Reaches | Status in this engagement |
|---|---|---|
| CPU tracer | every machine, no GPU or loader required | Ported live; the reference implementation |
| Vulkan compute | AMD, Intel, NVIDIA on Windows and Linux; Apple via MoltenVK; WSL2 via Dozen or Lavapipe; CPU-only machines via Lavapipe/SwiftShader | Built and gated on Lavapipe; validated on physical 780M silicon via Dozen (2026-07-28) |
| Slang→CUDA native | NVIDIA, if profiling shows Vulkan leaves performance unclaimed | Emission is a standing build gate: trace.slang → trace.cu every build, entry point pinned by the KernelPortability suite; adapter deferred |
| Slang→Metal native | Apple silicon | Emission is a standing build gate: trace.slang → trace.metal every build, entry point pinned by the KernelPortability suite; adapter deferred |

WSL2 is a named platform tier, not an afterthought: no native GPU kernel driver exists there, so the backend probes in order Dozen (Vulkan over the D3D12 paravirtual device), then Lavapipe, then declines to the CPU tracer with a clear message. The OptiX and CUDA host code, the PTX loader, and the vendored OptiX headers are retired with the kernel monolith; the Radeon 780M-class target is served by the Vulkan backend through the standard AMD drivers on Windows and RADV on Linux.

## 5. Render orchestration

The session orchestrator, tile scheduler, progress tracking, and viewer port forward with their roles intact. Tiles remain the unit of work everywhere: the CPU tracer threads over tiles, the Vulkan backend dispatches a compute grid per tile batch, and the viewer's progressive refinement consumes tiles in the existing spiral order. Beam propagation (specification P2) joins the live path in both backends: the kernel already carries deviation-vector state per ray from the GPU implementation, and the CPU tracer gains the same propagation, both validated against the double-precision oracle beam integrator, which remains off the render path.

## 6. Memory governor and precision ladder

The governor exists because the 780M-class target has a 2 GB budget that a naive full-frame HDR pipeline exhausts (a 5616 by 4096 IMAX frame at RGBA32F is 368 MB per buffer before ray state, which at 96 bytes per ray for position, momentum, deviation vectors, and accumulators is another 2.2 GB full-frame). The design bounds device residency by construction rather than by hope.

- At startup the backend reports usable budget (`VK_EXT_memory_budget` where present, else a conservative fraction of device-local memory, else a user cap). The governor derives the largest tile whose complete working set (ray state, accumulation, readback staging) fits a fixed fraction of that budget, and the tile scheduler simply receives a smaller tile size on smaller devices. Full-frame buffers live host-side only.
- The precision ladder has three rungs, selected per device and recorded in the render metadata: fp64 kernels where the device offers usable `shaderFloat64`; fp32 with double-single compensated arithmetic for the four integration-critical accumulations (Hamiltonian constraint, E, L_z, Q drift) elsewhere; plain fp32 for preview quality. The oracle's job is to quantify what each rung costs: per-rung error against the double-precision reference is a recorded test artefact, not a guess.
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
| Sirius.Infrastructure/Configuration/CRCF002A, CRCF003A, CRFM001A | app/config/{config_loader,config_schema,film_config}.h/.cpp |
| Sirius.Infrastructure/Platform/CRPF001A.{h,cpp} | app/platform_paths.{h,cpp} |
| Sirius.Infrastructure/Session/SMSM001A, SNST001A, SNEV001A, SNPR001A, SNDP001A, SNTL001A, SNRS001A | render/session/{state_machine,session_states,session_events,progress_tracker,display_buffer,tile_scheduler,render_session}.h/.cpp |
| Sirius.Infrastructure/Viewer/UIVW001A.h | app/viewer/interactive_viewer.h |
| Sirius.Test/TS*.cpp | tests/, one file per subject, names per style guide |

New files with no predecessor: `base/contracts.h`, `base/error.h`, `base/threading.h`, `kernels/*.slang`, `backend/vulkan/*`, `render/memory_governor.{h,cpp}`, `core/polarisation/walker_penrose.h` (specification E2), and the beam propagation additions to the CPU tracer (P2).

## 8. Dependencies

Kept, vendored in `lib/` as now: GLFW, glad, Dear ImGui, stb, tinyexr with miniz, glm (viewer mathematics only; Core keeps its own tensor types because glm is float-first and 3/4-component, the wrong shape for 4x4 double tensors with derivative stacks). Kept, fetched at configure time with pinned tags: FTXUI, nlohmann_json, GoogleTest. Added: Slang (build-time compiler, found at `/opt/slang` or `SLANG_ROOT`, plus optional runtime library for on-device kernel specialisation later), Vulkan loader and headers (system packages), and Vulkan Memory Allocator (vendored header) for the governor's budget-aware allocation. Removed: CUDA toolkit, OptiX SDK headers, the PTX multi-arch machinery, and glad's OpenGL-for-compute duties (GL remains only in the viewer blit path). Nothing else enters; in particular no scene-graph, ECS, or renderer-framework dependency, because the system's shape (analytic surfaces, one scene, physics-driven) gains nothing from them.

## 9. Validation architecture

Four gate families, layered exactly as the governing specification requires. Unit and property gates: the ported Mandatory suites (textbook metric values, Christoffel symmetry, exact Kerr-Schild identities, registry round-trips, conservation on the live path, determinism). Oracle gates: live-path agreement with the double-precision Boyer-Lindquist stack on conserved quantities, beam cross-sections against the oracle beam integrator, per-rung precision-ladder error quantification. Backend parity gates: CPU tracer versus Vulkan-on-Lavapipe on radiance statistics and conservation, kernel blackbody and Christoffel agreement suites carried from July, all executable in this environment because Lavapipe needs no hardware. Image gates: behaviour-preserving port stages reproduce the CPU reference tapes pixel-identically; behaviour-changing stages (beams on the CPU path, the kernel integrator on Vulkan) land with named re-baselined references and the difference explained in the commit.
