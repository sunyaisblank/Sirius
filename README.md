# Sirius

Sirius is a general relativistic ray tracing engine. It renders black holes, wormholes, and warp-drive spacetimes by integrating photon geodesics through analytically specified metric tensor fields, producing gravitational lensing, black hole shadows, and accretion disk emission from the geometry itself rather than from approximations layered on flat space. The engine is written in C++26 with a single Slang kernel source for GPU execution, and it targets the fidelity standard of DNGR, the renderer DNEG and Kip Thorne built for Interstellar; `docs/SPECIFICATION.md` states that standard as measurable criteria and `docs/ENGAGEMENT_REPORT.md` tracks the scorecard.

## The physics

Light in curved spacetime follows null geodesics, and the renderer exists to solve their equation of motion. For coordinates $x^\mu$ and an affine parameter $\lambda$ along the ray,

$$
\frac{d^2 x^\mu}{d\lambda^2} + \Gamma^\mu_{\nu\rho} \frac{dx^\nu}{d\lambda} \frac{dx^\rho}{d\lambda} = 0,
$$

where the Christoffel symbols derive from the metric and its first derivatives. Sirius supplies those derivatives analytically: each metric family implements closed-form $g_{\mu\nu}$ and $\partial_\sigma g_{\mu\nu}$, which removes the truncation error finite differencing would accumulate over the hundreds of steps a ray takes near a photon sphere. The CPU reference path integrates a Dormand-Prince RK45 scheme in Hamiltonian form with covariant momenta and per-step renormalisation of the null condition; the GPU kernels integrate with the time-transformed Yoshida symplectic family. Two independent methods agreeing within stated tolerance is stronger evidence than one method twice, so backend parity is gated statistically, and conservation of energy and angular momentum on the live path is enforced by build-failing test rather than assumed.

The black hole family uses the Kerr-Schild form $g_{\mu\nu} = \eta_{\mu\nu} + H\,l_\mu l_\nu$ with null $l$, which has no coordinate singularity at the poles and yields two exact identities the tests exploit: the closed-form inverse $g^{\mu\nu} = \eta^{\mu\nu} - H\,l^\mu l^\nu$ and $\det g = -1$. One four-parameter family $(M, a, Q, \Lambda)$ covers Minkowski, Schwarzschild, Kerr, Reissner-Nordstrom, Kerr-Newman, de Sitter, and Schwarzschild-de Sitter. Accretion disk emission uses the Novikov-Thorne temperature profile with the Page-Thorne (1974) flux function, blackbody spectral colour, and $g^4$ relativistic beaming; polarisation is transported physically along geodesics and certified by conservation of the Walker-Penrose complex constant. A double-precision Boyer-Lindquist oracle stack (metric with analytic Christoffels and Riemann, symplectic integrator, geodesic-deviation beam integrator, polarisation transport) lives deliberately off the render path as the validation reference the live paths are tested against.

## What it renders

| Metric | Parameters | Backends |
|--------|------------|----------|
| Minkowski | none | CPU, Vulkan |
| Schwarzschild | mass | CPU, Vulkan |
| Kerr | mass, spin | CPU, Vulkan |
| Reissner-Nordstrom | mass, charge | CPU only |
| Kerr-Newman | mass, spin, charge | CPU only |
| de-Sitter | lambda | CPU only |
| Schwarzschild-de-Sitter | mass, lambda | CPU only |
| Morris-Thorne | throat radius | deferred (declines on both; kernel physics parity-tested) |
| Alcubierre | warp velocity, bubble radius, wall thickness | CPU, Vulkan |

`sirius info metrics` prints this catalogue from the same registry the validator and routers use. Requesting a metric on a backend that cannot represent it is a clear error; the engine never substitutes a different spacetime. Output format follows the file extension: `.png` and `.ppm` receive the display pipeline (tonemap, grade, exactly one transfer encode applied by the owning writer); `.exr` receives linear HDR radiance untouched, for compositing.

## Backends

The CPU path is the reference implementation and runs everywhere, threads over spiral-ordered tiles, and produced the pinned reference images the test estate compares against byte-for-byte. The GPU path is one Slang kernel source (`src/sirius/kernels/`) compiled to SPIR-V and dispatched through a Vulkan compute seam, which reaches AMD, Intel, and NVIDIA hardware on Windows and Linux, Apple silicon through MoltenVK, and WSL2 or GPU-less machines through Lavapipe software Vulkan; the same source emits CUDA and Metal for future native adapters. A memory governor sizes tiles to the device budget (a 2 GB integrated GPU is a supported render target, not merely a supported install), and kernel-versus-oracle parity suites run on Lavapipe in CI, so the GPU physics is tested on machines with no GPU at all. GPU rendering is explicit opt-in (`--backend vulkan` or `SIRIUS_RENDER_BACKEND=vulkan`); `auto` resolves to the CPU until the go-live decision flips it.

## Building

Requirements: CMake 3.28+, and a compiler with a C++26 mode: GCC 14+, Clang 17+, AppleClang 16+, or MSVC 19.40+ (Visual Studio 2022 17.10). The Vulkan backend additionally wants the Vulkan headers/loader and the Slang compiler (`slangc`, found at `/opt/slang` or `SLANG_ROOT`); both degrade with a clear message when absent, leaving the CPU system complete.

```bash
# Linux / WSL2
cmake --preset linux-gcc && cmake --build --preset linux-gcc
./bin/linux-gcc/src/sirius/app/sirius

# Windows
cmake --preset windows-msvc && cmake --build --preset windows-msvc
```

Presets exist for GCC and Clang (Release and Debug), MSVC, and macOS; `linux-ci` enables the Mandatory build gate. On WSL2 no extra setup is needed beyond `libvulkan-dev` and `mesa-vulkan-drivers`: Lavapipe provides the Vulkan device the backend suites use.

## Running

```bash
# A Kerr black hole with its disk, 256 samples per pixel
sirius render -m Kerr -a 0.9 -s 256 -o kerr.png

# The same scene as linear HDR for compositing
sirius render -m Kerr -a 0.9 -s 256 -o kerr.exr

# Through the Vulkan backend (explicit opt-in)
sirius render -m Kerr -a 0.9 --backend vulkan -o kerr_vk.png

# Interactive progressive viewer, machine-readable system report
sirius view -m Kerr -a 0.9
sirius info system --json
```

Configuration layers in a fixed order: struct defaults, then a JSON config file, then `SIRIUS_*` environment variables, then command-line flags; later layers win, every parameter is validated at startup against the ranges the physics supports, and invalid configuration stops the run rather than being clamped. `sirius render --help` lists the full flag set.

The DNGR-technique features are explicit flags, each default-off so the pinned reference render is bit-stable: `--beams` propagates ray bundles by the geodesic deviation equation and reports the beam ellipse per pixel; `--starfield point` replaces the background texture with a 100,000-star point catalogue filtered through the beam footprint (the anti-flicker approach of the Interstellar renderer); `--doppler-beaming off` suppresses the disk's orbital Doppler asymmetry while retaining gravitational redshift, mirroring the film's artistic choice (full physics is the default); a camera velocity applies special-relativistic aberration across all lens models.

## Testing

```bash
ctest --test-dir bin/linux-gcc -j"$(nproc)"          # everything
ctest --test-dir bin/linux-gcc -L Mandatory          # the build gate
```

Labels carry gating semantics: Mandatory suites (textbook metric values, exact Kerr-Schild identities, conservation on the live path, kernel parity, Walker-Penrose conservation, the registry) fail the build under `-DSIRIUS_MANDATORY_TESTS=ON`; Correctness suites cover feature behaviour; Performance suites warn only. The label file is generated from the test sources by `scripts/generate-ctest-labels.py` and CI rejects a stale copy, so a renamed suite cannot silently fall out of the gate.

## Documentation

`docs/SPECIFICATION.md` (the DNGR-parity mandate and criteria), `docs/ARCHITECTURE.md` (layers, single authorities, kernel and backend design, the legacy-name mapping), `docs/STYLE.md` (the C++26 style guide), `docs/ENGAGEMENT_REPORT.md` (the rebuild's closing report and scorecard).

## References

- Misner, Thorne and Wheeler (1973), *Gravitation*. W. H. Freeman.
- Visser (2007), "The Kerr spacetime: a brief introduction", arXiv:0706.0622.
- Page and Thorne (1974), "Disk-accretion onto a black hole", ApJ 191, 499.
- Morris and Thorne (1988), "Wormholes in spacetime and their use for interstellar travel", Am. J. Phys. 56, 395.
- James, von Tunzelmann, Franklin and Thorne (2015), "Gravitational lensing by spinning black holes in astrophysics, and in the movie Interstellar", Class. Quantum Grav. 32, 065001.
- Walker and Penrose (1970), "On quadratic first integrals of the geodesic equations for type {22} spacetimes", Commun. Math. Phys. 18, 265.
- Bardeen, Press and Teukolsky (1972), "Rotating black holes", ApJ 178, 347.
