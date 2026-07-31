# Sirius

Sirius is a general relativistic ray tracing engine. It renders black holes, wormholes, and warp-drive spacetimes by integrating photon geodesics through analytically specified metric tensor fields, producing gravitational lensing, black hole shadows, and accretion disk emission from the geometry itself rather than from approximations layered on flat space. The engine is written in C++26 with a single Slang kernel source for GPU execution, and it targets the fidelity standard of DNGR, the renderer DNEG and Kip Thorne built for Interstellar; `docs/SPECIFICATION.md` states that standard as measurable criteria and `docs/ADVERSARIAL_REVIEW.md` records the current scorecard and limits.

## The physics

Light in curved spacetime follows null geodesics, and the renderer exists to solve their equation of motion. For coordinates $x^\mu$ and an affine parameter $\lambda$ along the ray,

$$
\frac{d^2 x^\mu}{d\lambda^2} + \Gamma^\mu_{\nu\rho} \frac{dx^\nu}{d\lambda} \frac{dx^\rho}{d\lambda} = 0,
$$

where the Christoffel symbols derive from the metric and its first derivatives. Sirius supplies those derivatives analytically. The CPU reference path integrates the Cartesian geodesic equation with adaptive Dormand-Prince RK45 and null normalisation; the live Vulkan kernel uses Cartesian RK4. Backend parity is statistical rather than bitwise, and conservation on the live path is a build-failing test.

The black hole family uses the Kerr-Schild form $g_{\mu\nu} = \eta_{\mu\nu} + H\,l_\mu l_\nu$ with null $l$, which has no coordinate singularity at the poles and yields two exact identities the tests exploit: the closed-form inverse $g^{\mu\nu} = \eta^{\mu\nu} - H\,l^\mu l^\nu$ and $\det g = -1$. One four-parameter family $(M, a, Q, \Lambda)$ covers Minkowski, Schwarzschild, Kerr, Reissner-Nordstrom, Kerr-Newman, de Sitter, and Schwarzschild-de Sitter. Schwarzschild and Kerr accretion-disk emission uses the Novikov-Thorne temperature profile with the full Page-Thorne (1974) flux function, blackbody spectral colour, and exactly one $g^4$ relativistic-beaming factor. A double-precision Boyer-Lindquist oracle stack (metric with analytic Christoffels and Riemann, symplectic integrator, geodesic-deviation beam integrator, Walker-Penrose polarisation transport) lives deliberately off the render path as the validation reference. On CPU, `Polarisation` mode parallel-transports an observer screen basis along the live Kerr-Schild ray, forms thermal Thomson Stokes emission in the Page-Thorne disk frame, and carries the accumulated state into the film buffer. Live and oracle transport are independently pinned by the Walker-Penrose constant; Vulkan and volumetric polarisation requests decline explicitly.

## What it renders

| Metric | Parameters | Backends | Accretion disk |
|--------|------------|----------|-----------------|
| Minkowski | none | CPU, Vulkan | Not applicable; use `--no-disk` |
| Schwarzschild | mass | CPU, Vulkan | Page-Thorne |
| Kerr | mass, spin | CPU, Vulkan | Page-Thorne |
| Reissner-Nordstrom | mass, charge | CPU only | Unavailable; use `--no-disk` |
| Kerr-Newman | mass, spin, charge | CPU only | Unavailable; use `--no-disk` |
| de-Sitter | lambda | CPU only | Not applicable; use `--no-disk` |
| Schwarzschild-de-Sitter | mass, lambda | CPU only | Unavailable; use `--no-disk` |
| Morris-Thorne | throat radius | CPU, Vulkan | Not applicable; use `--no-disk`; one sheet, dark throat |
| Alcubierre | warp velocity, bubble radius, wall thickness | CPU, Vulkan | Not applicable; use `--no-disk` |

`sirius info metrics` prints this catalogue and disk support from the same registry the validator and routers use. Requesting a metric, disk, or backend combination that cannot be represented is a clear error; the engine never substitutes a different spacetime or disk model. Output format follows the file extension: `.png` and `.ppm` receive the display pipeline (tonemap, grade, exactly one transfer encode applied by the owning writer); `.exr` receives finite fp32 linear HDR radiance untouched plus Sirius render metadata, for compositing.

## Backends

The CPU path is the reference implementation and runs without Vulkan or OpenGL. It supports pinhole, thin-lens, and fisheye projection plus the physical thin-disk polarisation output. The GPU path is one Slang kernel source (`src/sirius/kernels/`) compiled to SPIR-V and dispatched through Vulkan; the same source is also compiled to CUDA and Metal as a portability gate, not as a CUDA or Metal runtime backend. A memory governor sizes tiles against the host-visible coherent heap the adapter actually allocates from, and a dispatch governor bounds device occupancy. `SIRIUS_VULKAN_DEVICE=<zero-based-index>` selects one enumerated device strictly; the readiness report, renderer, tests, and hardware attestation use the same identity. Backend `auto` selects Vulkan only when a device, metric, and requested scene semantics are all represented; otherwise it logs the reason and selects CPU. Vulkan represents exact arbitrary sample counts, camera aberration, pinhole and finite-aperture thin lenses, Page-Thorne and Shakura-Sunyaev temperature modes, Doppler suppression, dual-vector ray bundles, ellipse-filtered 100,000-star catalogues, and volumetric disk/turbulence/corona transfer. Vulkan volumetric marching is deliberately bounded at 128 midpoint samples per geodesic segment; larger explicit requests decline and `auto` uses CPU. Diagnostic colour and polarisation modes are CPU-only. The precision ladder offers plain fp32, `SIRIUS_PRECISION=fp32-comp`, and `SIRIUS_PRECISION=fp64` on devices exposing `shaderFloat64`.

## Building

Requirements: CMake 3.28+, and a compiler with a C++26 mode: GCC 14+, Clang 17+, AppleClang 16+, or MSVC 19.40+ (Visual Studio 2022 17.10). The Vulkan backend additionally wants the Vulkan headers/loader and the Slang compiler (`slangc`, found at `/opt/slang` or `SLANG_ROOT`). The interactive window uses OpenGL and GLFW; when OpenGL development files are absent, or `SIRIUS_BUILD_VIEWER=OFF`, Sirius builds the complete CPU/Vulkan CLI and the `view` command declines explicitly. Portable builds likewise degrade clearly to a complete CPU system when Vulkan is absent; the `linux-ci` operational profile sets `SIRIUS_REQUIRE_VULKAN_RUNTIME=ON` and rejects configuration, missing kernels, runtime skips, or failed dispatch.

```bash
# Linux / WSL2
cmake --preset linux-gcc && cmake --build --preset linux-gcc
./bin/linux-gcc/src/sirius/app/sirius

# Windows
cmake --preset windows-msvc && cmake --build --preset windows-msvc
```

Presets exist for GCC and Clang (Release and Debug), MSVC, macOS, required Linux Vulkan CI, and GCC ASan/UBSan/LSan. Every supported preset pins warnings as errors, enforce-mode contracts, Mandatory tests, and source governance; `linux-ci` additionally runs the non-skipping Vulkan/P1 operational profile. `cmake --install` creates a relocatable volume containing the binary, starfield, installed operating model, every compiled backend resource, and viewer shaders when the viewer is compiled, then verifies that exact capability-specific structure before succeeding.

## Running

```bash
# A Kerr black hole with its disk, 256 samples per pixel
sirius render -m Kerr -a 0.9 -s 256 -o kerr.png

# The same scene as linear HDR for compositing
sirius render -m Kerr -a 0.9 -s 256 -o kerr.exr

# Transported thermal-disk polarisation (CPU selected automatically)
sirius render -m Kerr -a 0.9 --color-mode Polarisation -o kerr-polarisation.png

# Scalar CPU temporal disk integration
sirius render -m Kerr -a 0.9 --cpu --motion-blur --shutter-time 0.25 \
  --motion-samples 7 -o kerr-motion.png

# Explicitly require Vulkan with a volumetric turbulent disk and corona
sirius render -m Kerr -a 0.9 -s 4 --backend vulkan --volumetric \
  --turbulence --corona -o kerr_vk.png

# Interactive progressive viewer, machine-readable system report
sirius view --spin 0.9 --backend auto
sirius info system --json
sirius info readiness --json
sirius info capabilities --json
```

Configuration layers in a fixed order: struct defaults, then a JSON config file, then `SIRIUS_*` environment variables, then command-line flags; later layers win, every parameter is validated at startup against the ranges the physics supports, and invalid configuration stops the run rather than being clamped. `sirius render --help` lists the full flag set.

`render_test.sh` and `render_demo.sh` select only the exact preset binary (or an explicit `SIRIUS_BINARY`), propagate renderer failures, and require every declared output to exist and be non-empty. They never search for an arbitrary stale build.

The DNGR-technique features are explicit flags, each default-off so the pinned reference render is stable: `--beams` propagates two ray-bundle deviation vectors and derives the oriented beam ellipse; `--starfield point` replaces the background texture with a spatially indexed 100,000-star point catalogue filtered through both ellipse axes; `--volumetric`, `--turbulence`, and `--corona` select live radiative transfer; `--motion-blur` selects scalar CPU temporal disk integration; `--doppler-beaming off` suppresses the disk's orbital Doppler asymmetry while retaining gravitational redshift; `--color-mode Polarisation` selects transported CPU Stokes output; and a camera velocity applies special-relativistic aberration across all lens models. Polarisation combined with volume, temporal blur, or Vulkan declines. `--wormhole-topology TwoSheet` is a named but fail-closed request; `OneSheetCapture` is the represented dark-throat topology.

## Testing

```bash
ctest --test-dir bin/linux-gcc -j"$(nproc)"          # everything
ctest --test-dir bin/linux-gcc -L Mandatory          # the build gate

# Instrumented mandatory gate
cmake --preset linux-gcc-sanitize
cmake --build --preset linux-gcc-sanitize -j"$(nproc)"
```

Labels carry gating semantics: Mandatory and Operational suites fail every normal build by default; Correctness suites cover additional feature behaviour; Performance suites measure without gating. New suites without an explicit policy are a hard generator error. Parameterised/typed GoogleTest macros are also rejected until the generator can prove their exact discovered names, preventing tests from escaping the Mandatory label. `tests/operating_model.json` maps 22 required operating dimensions and 19 explicit capability contracts to Mandatory evidence. The build rejects a missing state, renamed witness, or de-gated claim; the same model is installed and available through `info capabilities`.

## Documentation

`docs/SPECIFICATION.md` is the target, while `docs/ADVERSARIAL_REVIEW.md` is the current evidence-based status and limitation ledger. `docs/ENGAGEMENT_REPORT.md` is retained as historical closure evidence and is not current ground truth.

## References

- Misner, Thorne and Wheeler (1973), *Gravitation*. W. H. Freeman.
- Visser (2007), "The Kerr spacetime: a brief introduction", arXiv:0706.0622.
- Page and Thorne (1974), "Disk-accretion onto a black hole", ApJ 191, 499.
- Morris and Thorne (1988), "Wormholes in spacetime and their use for interstellar travel", Am. J. Phys. 56, 395.
- James, von Tunzelmann, Franklin and Thorne (2015), "Gravitational lensing by spinning black holes in astrophysics, and in the movie Interstellar", Class. Quantum Grav. 32, 065001.
- Walker and Penrose (1970), "On quadratic first integrals of the geodesic equations for type {22} spacetimes", Commun. Math. Phys. 18, 265.
- Bardeen, Press and Teukolsky (1972), "Rotating black holes", ApJ 178, 347.
