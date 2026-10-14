# Sirius

Sirius is a general relativistic ray tracing engine. It renders black holes, wormholes, and warp-drive spacetimes by integrating photon geodesics through analytically specified metric tensor fields, producing gravitational lensing, black hole shadows, and accretion disk emission from the geometry itself rather than from approximations layered on flat space.

## The physics

Light in curved spacetime follows null geodesics, and the entire renderer exists to solve their equation of motion. For coordinates $x^\mu$ and an affine parameter $\lambda$ along the ray,

$$
\frac{d^2 x^\mu}{d\lambda^2} + \Gamma^\mu_{\nu\rho} \frac{dx^\nu}{d\lambda} \frac{dx^\rho}{d\lambda} = 0,
$$

where the Christoffel symbols $\Gamma^\mu_{\nu\rho}$ encode the curvature. They derive from the metric tensor and its first derivatives,

$$
\Gamma^\mu_{\nu\rho} = \frac{1}{2} g^{\mu\sigma} \left( \partial_\nu g_{\sigma\rho} + \partial_\rho g_{\sigma\nu} - \partial_\sigma g_{\nu\rho} \right),
$$

so the accuracy of every rendered pixel reduces to the accuracy of the metric derivatives and of the integrator. Sirius supplies the derivatives analytically: each metric family implements closed-form expressions for $g_{\mu\nu}$ and $\partial_\sigma g_{\mu\nu}$, which removes the truncation error that finite differencing would accumulate over the hundreds of integration steps a ray takes near a photon sphere. Integration uses a Dormand-Prince RK45 scheme in a Hamiltonian formulation with covariant momenta, with the null condition $g_{\mu\nu} k^\mu k^\nu = 0$ monitored and re-normalised each step. Conservation of energy and angular momentum in stationary axisymmetric spacetimes is enforced by test rather than assumed: the mandatory suite integrates rays on the live path and fails the build if drift exceeds $10^{-4}$ relative.

The black hole family uses the Kerr-Schild form,

$$
g_{\mu\nu} = \eta_{\mu\nu} + H \, l_\mu l_\nu,
$$

with $\eta_{\mu\nu}$ the Minkowski metric, $H$ a scalar field, and $l_\mu$ a null vector. This Cartesian representation was chosen over Boyer-Lindquist coordinates because it has no coordinate singularity at the rotation poles and because nullity of $l$ gives two exact identities the code exploits: the inverse metric is $g^{\mu\nu} = \eta^{\mu\nu} - H l^\mu l^\nu$ in closed form, and $\det g = -1$ everywhere, both of which serve as test oracles. One four-parameter family $(M, a, Q, \Lambda)$ covers Minkowski, Schwarzschild, Kerr, Reissner-Nordstrom, Kerr-Newman, de Sitter, and Schwarzschild-de Sitter; the cosmological constant folds into $H$ exactly for $a = 0$, and the rotating de Sitter forms are rejected at the configuration boundary rather than approximated.

Accretion disk emission uses the Novikov-Thorne relativistic thin-disk temperature profile with the Page and Thorne (1974) flux function, or a Shakura-Sunyaev power law when selected. Brightness scales as the physical $T^4$ against a single reference temperature. A volumetric disk mode ray-marches a Gaussian vertical structure with optional turbulence and corona; these, together with film-stock simulation, are aesthetic features and are labelled as such rather than presented as physics.

## What it renders

| Metric | Parameters | Backends |
|--------|------------|----------|
| Minkowski | none | CPU, GPU |
| Schwarzschild | mass | CPU, GPU |
| Kerr | mass, spin | CPU, GPU |
| Reissner-Nordstrom | mass, charge | CPU only |
| Kerr-Newman | mass, spin, charge | CPU only |
| de-Sitter | lambda | CPU only |
| Schwarzschild-de-Sitter | mass, lambda | CPU only |
| Morris-Thorne | throat radius | GPU only |
| Alcubierre | warp velocity, bubble radius, wall thickness | CPU, GPU |

`sirius info metrics` prints this table from the same registry the validator and routers use, so it cannot drift from the code. Names are case-insensitive and common aliases are accepted (`Wormhole`, `WarpDrive`, `KerrNewman`). Requesting a metric on a backend that cannot represent it is an error with a clear message; the engine never substitutes a different spacetime. The two asymmetries are deliberate: charge and the cosmological constant are not yet plumbed into the GPU kernel, and the Morris-Thorne implementation evaluates in spherical coordinates that the Cartesian CPU tracer cannot drive (its Cartesian embedding is recorded follow-up work). The wormhole is the zero-tidal-force subclass with constant redshift function, a physically legitimate choice from Morris and Thorne (1988), not an approximation.

Output format follows the file extension. `.ppm` and `.png` receive the display pipeline (tonemapping, grading, then exactly one transfer encode, applied by the writer that owns the format); `.exr` receives the linear HDR radiance untouched, for downstream compositing. Tonemappers: ACES (default), Reinhard, Filmic, None.

## Building

Requirements: CMake 3.16 or newer and a C++17 compiler. CUDA 12 and the OptiX 7/8/9 SDK enable the GPU backend; without them the build degrades to the CPU path with a warning, and CI builds both configurations. OpenGL and GLFW (vendored) support the interactive viewer; FTXUI and nlohmann-json are fetched at configure time.

```bash
# Linux / WSL2
cmake -B bin/Sirius.Build && cmake --build bin/Sirius.Build
./bin/Sirius.Build/bin/Sirius

# Windows
cmake -B bin/Sirius.Build -G "Visual Studio 17 2022" -A x64
cmake --build bin/Sirius.Build --config Release
```

OptiX headers are searched in `lib/optix/`, `/opt/nvidia/optix/`, `~/NVIDIA-OptiX-SDK/`, the `OptiX_INSTALL_DIR` environment variable, and the standard Windows install paths. For WSL2, `scripts/setup-wsl2-deps.sh`, `setup-wsl2-cuda.sh`, and `setup-wsl2-optix.sh` install the toolchain end to end.

Build options: `BUILD_TESTS` (default on), `SIRIUS_PTX_ALL_ARCHS` (PTX for SM 75/80/86/89, default on), and `SIRIUS_MANDATORY_TESTS`, which makes the build itself fail unless every test labelled Mandatory passes; `SIRIUS_SKIP_MANDATORY` overrides it for development.

## Running

The binary has four commands. `render` produces an image in batch, `view` opens an interactive OpenGL window with progressive refinement (WASD orbit, mouse look, scroll zoom), `info` reports the system, metric catalogue, backends, and effective configuration, and `config` manages configuration files (`show`, `validate`, `init`, `paths`).

```bash
# A Kerr black hole with its disk, 256 samples per pixel
sirius render -m Kerr -a 0.9 -s 256 -o kerr.png

# The same scene as linear HDR for compositing
sirius render -m Kerr -a 0.9 -s 256 -o kerr.exr

# A wormhole (GPU backend required), and a warp bubble on the CPU
sirius render -m Morris-Thorne --throat-radius 1.5 -o wormhole.ppm
sirius render -m Alcubierre --warp-velocity 0.8 --cpu -o warp.png

# Machine-readable system report
sirius info system --json
```

Configuration layers in a fixed order: struct defaults, then a JSON config file (`--config` or the standard search paths), then `SIRIUS_*` environment variables (`SIRIUS_METRIC`, `SIRIUS_SPIN`, `SIRIUS_WIDTH`, `SIRIUS_BACKEND`, and the rest listed in `sirius --help`), then command-line flags. Later layers win. Global flags (`--json`, `--verbose`, `--no-color`, `--config`) are recognised anywhere on the command line. Every parameter is validated at startup against the ranges the physics supports (spin up to the 0.998 Thorne limit, sub-extremal spin-charge combinations, lambda only with zero spin); invalid configuration stops the run rather than being silently clamped.

The main render options: `-w/-h` resolution, `-s` samples per pixel, `-m` metric, `-a` spin, `-d` observer distance, `-i` inclination, `--fov`, `--temperature-model` (NovikovThorne or ShakuraSunyaev), `--tonemapper`, `--exposure`, `--bloom`, `--volumetric` with `--turbulence` and `--corona`, `--film` with `--film-preset`, and `--cpu`/`--no-gpu` to force the CPU backend. `sirius render --help` gives the full list.

`render_test.sh` renders a two-image smoke test and `render_demo.sh` a thirteen-image suite across metrics, inclinations, and feature combinations; `scripts/benchmark-gpu-cpu.sh` compares backend performance across resolutions.

## Architecture

```
Sirius.Infrastructure  (CLI, configuration, render session, viewer)
        ↓
Sirius.Render          (OptiX kernel and host, CPU tracer, output writers)
        ↓
Sirius.Core            (metrics, geodesics, tensors, disk physics, constants)
```

Dependencies point downward only; the physics layer compiles and tests without a GPU or a window. Two decisions shape the layout beyond the layering. First, metric identity is a closed enum with a single registry (`PHMT200A`) that every consumer parses through, because the previous free-form strings drifted into six divergent catalogues and a wormhole request once silently rendered a black hole. Second, alongside the live single-precision Cartesian path, `Sirius.Core` keeps a double-precision Boyer-Lindquist stack (metric, symplectic integrator, beam integrator) that is deliberately not on the render path: it is the validation oracle the test suite integrates against, and the live path is required by test to agree with it on conserved quantities.

One caveat is worth stating honestly: the GPU kernel (`RDOP002A.cu`) carries its own inlined copies of the metric and disk physics, written for single-precision CUDA. Parity between kernel and Core is enforced where tests can reach it (blackbody colour, Christoffel agreement) and is otherwise a known maintenance surface; a harness that executes the kernel under test needs a GPU in the loop and is recorded follow-up work.

Source files follow a `[Domain][Category][Sequence][Variant]` naming convention (`PHMT100A.h` reads as Physics, Metric Tensor, family 100, production variant). The scheme is terse deliberately: many small, tightly scoped files would otherwise produce long near-identical descriptive names. Each file's header comment states its component identity and role.

## Testing

```bash
ctest --test-dir bin/Sirius.Build --output-on-failure   # everything
ctest --test-dir bin/Sirius.Build -L Mandatory          # the build gate
```

Tests carry labels with distinct gating semantics: Mandatory suites (metric properties against textbook values, Christoffel symmetry, conservation on the live path, determinism, the registry, the exact Kerr-Schild inverse identities) fail the build when `SIRIUS_MANDATORY_TESTS` is on; Correctness suites cover feature behaviour; Performance suites warn without failing. The label file is generated from the test sources by `scripts/generate-ctest-labels.py` and CI rejects a stale copy, so a renamed suite cannot silently fall out of the gate.

## References

- Misner, Thorne and Wheeler (1973), *Gravitation*. W. H. Freeman.
- Visser (2007), "The Kerr spacetime: a brief introduction", arXiv:0706.0622.
- Page and Thorne (1974), "Disk-accretion onto a black hole", ApJ 191, 499.
- Morris and Thorne (1988), "Wormholes in spacetime and their use for interstellar travel", Am. J. Phys. 56, 395.
- James, von Tunzelmann, Franklin and Thorne (2015), "Gravitational lensing by spinning black holes in astrophysics, and in the movie Interstellar", Class. Quantum Grav. 32, 065001.
- Bardeen, Press and Teukolsky (1972), "Rotating black holes", ApJ 178, 347.
