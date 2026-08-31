# Sirius

Sirius is a general relativistic ray tracing engine. It renders black holes, wormholes, and warp-drive spacetimes by integrating photon geodesics through analytically specified metric tensor fields, producing gravitational lensing, black hole shadows, and accretion disk emission from the geometry itself rather than from approximations layered on flat space. The engine is written in C++26 with a single Slang kernel source for GPU execution, and it targets the fidelity standard of DNGR, the renderer DNEG and Kip Thorne built for Interstellar; `docs/SPECIFICATION.md` states that standard as measurable criteria and `docs/ADVERSARIAL_REVIEW.md` records the current scorecard and limits.

## The physics

Light in curved spacetime follows null geodesics, and the renderer exists to solve their equation of motion. For coordinates $x^\mu$ and an affine parameter $\lambda$ along the ray,

$$
\frac{d^2 x^\mu}{d\lambda^2} + \Gamma^\mu_{\nu\rho} \frac{dx^\nu}{d\lambda} \frac{dx^\rho}{d\lambda} = 0,
$$

where the Christoffel symbols derive from the metric and its first derivatives. Sirius supplies those derivatives analytically. The CPU reference path integrates the Cartesian Hamilton equations with adaptive Dormand-Prince RK45. The live Vulkan kernel uses adaptive Cartesian RK4 step doubling, with compensated or fp64 state accumulation on the corresponding precision rungs. Both reject steps whose scale-aware truncation error or relative null residual exceeds tolerance. The device additionally projects only the accepted temporal tangent onto the nearest represented null root, preserving the spatial tangent and light-cone branch; the CPU path does not project. Their shared trace domain scales step bounds and the far boundary with the selected metric's natural geometry—black-hole mass `M`, Morris-Thorne throat `b0`, or the Alcubierre radius/inverse-wall scales—and always encloses the observer. In the CPU positive-$\Lambda$ sector the ordinary finite far sphere is additionally clamped to the exact cosmological horizon and the accepted crossing segment is localised back onto that boundary, so a directional background is never sampled beyond the represented causal patch. Backend parity is statistical rather than bitwise, and full-ray conservation on each live path is a build-failing test.

The black hole family uses the Kerr-Schild form $g_{\mu\nu} = \eta_{\mu\nu} + H\,l_\mu l_\nu$ with null $l$, which has no coordinate singularity at the poles and yields two exact identities the tests exploit: the closed-form inverse $g^{\mu\nu} = \eta^{\mu\nu} - H\,l^\mu l^\nu$ and $\det g = -1$. One four-parameter family $(M, a, Q, \Lambda)$ covers Minkowski, Schwarzschild, Kerr, Reissner-Nordstrom, Kerr-Newman, de Sitter, and the spherical positive-$\Lambda$ Schwarzschild-de Sitter (Kottler) sector. There the black-hole and cosmological horizons are the two positive roots of $1-2M/r-\Lambda r^2/3=0$; capture uses the smaller root, pure de Sitter has no capture surface, and named Schwarzschild-de Sitter scenes require $9\Lambda M^2<1$. Negative $\Lambda$ and rotating or charged cosmological sectors decline. The camera worldline is the metric-aware ADM/Eulerian observer in the horizon-penetrating chart, not an implied static observer, but its launch radius must remain at or below $0.99r_c$ inside the cosmological horizon. The star texture is a finite directional boundary condition applied no later than $r_c$, not a cosmological radiative-transfer model. Schwarzschild and Kerr thin-disk emission defaults to the zero-torque Page-Thorne (1974) flux; the explicitly selected `ShakuraSunyaev` mode is a bounded Newtonian substitute with $F(r)\propto r^{-3}[1-\sqrt{r_\mathrm{in}/r}]$, not a relativistic or complete alpha-disk model. Both modes use Planck spectral radiance, invariant $I_\nu/\nu^3$ transfer, and the corresponding single bolometric $g^4$ factor. A double-precision Boyer-Lindquist oracle stack (metric with analytic Christoffels and Riemann, a fixed-step implicit-midpoint/Yoshida symplectic map, geodesic-deviation beam integrator, Walker-Penrose polarisation transport) lives deliberately off the render path as the validation reference. Its canonical two-form is tested directly; optional null projection is an oracle-only stabiliser, not a live method. On CPU, `Polarisation` mode parallel-transports an observer-screen basis along the accepted Kerr-Schild ray, projects it into the circular emitter screen, and applies the flux-normalised Chandrasekhar-Sobolev semi-infinite pure electron-scattering atmosphere. That closure excludes absorption, finite optical depth, returning radiation, and magnetic/Faraday effects. Live and oracle transport are independently pinned by the Walker-Penrose constant; Vulkan and volumetric polarisation requests decline explicitly.

## What it renders

| Metric | Parameters | Backends | Accretion disk |
|--------|------------|----------|-----------------|
| Minkowski | none | CPU, Vulkan | Not applicable; use `--no-disk` |
| Schwarzschild | mass | CPU, Vulkan | Page-Thorne |
| Kerr | mass, spin | CPU, Vulkan | Page-Thorne |
| Reissner-Nordstrom | mass, charge | CPU only | Unavailable; use `--no-disk` |
| Kerr-Newman | mass, spin, charge | CPU only | Unavailable; use `--no-disk` |
| de-Sitter | positive lambda | CPU only | Not applicable; use `--no-disk` |
| Schwarzschild-de-Sitter | mass, positive lambda; $9\Lambda M^2<1$ | CPU only | Unavailable; use `--no-disk` |
| Morris-Thorne | Ellis throat radius | CPU, Vulkan | Exact isotropic chart; no disk; one output sheet, dark throat |
| Alcubierre | warp velocity, bubble radius, inverse wall scale `sigma` | CPU, Vulkan | Not applicable; use `--no-disk` |

`sirius info metrics` prints this catalogue and disk support from the same registry the validator and routers use. Requesting a metric, disk, or backend combination that cannot be represented is a clear error; the engine never substitutes a different spacetime or disk model. Output format follows the file extension: `.png` and `.ppm` receive the display pipeline (tonemap, grade, exactly one transfer encode applied by the owning writer); `.exr` receives finite fp32 linear HDR radiance untouched plus Sirius render metadata, for compositing. The interactive viewer consumes the same display-linear result and its sole live shader applies exactly one IEC sRGB encode—no second tone map or approximate gamma curve. The default `ACESFit` tonemapper is the scalar [Narkowicz fit](https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/) to sampled ACES 1 output, independently applied to working-RGB channels. It is not an [Academy ACES Output Transform](https://docs.acescentral.com/system-components/output-transforms/): Sirius does not claim the corresponding ACES working-space, chroma/gamut, display-primary, white-point, viewing-condition, or display-encoding semantics, and a bare `ACES` request declines rather than silently selecting the fit.

The binned spectral facade integrates physical spectral radiance into tristimulus values whose scale remains radiometric and deliberately stops at `ToXyz()`. It does not apply a D65-to-sRGB matrix and clipped transfer encoding to an arbitrary absolute scale: [sRGB is display-referred, D65, and has a defined in-gamut reference range](https://www.w3.org/TR/css-color-4/#predefined-sRGB), so a preview first needs an explicitly represented exposure or normalisation and display intent. The former unused `ToSrgb()` method relied on callers inventing that scale and was removed; the distinct live blackbody-colour path retains its explicit brightest-channel normalisation. The binned facade likewise does not expose its values as [ACES2065-1](https://docs.acescentral.com/encodings/aces2065-1/) relative exposure. Its former `ToAces()` method merely multiplied by an AP0 coordinate matrix and was removed: a primary matrix can change coordinates, but cannot supply the scene-relative exposure normalisation required to establish an ACES encoding.

Metric mass `M`, observer distance, thin-lens focus distance, throat radius, and bubble radius use one scale-free geometric coordinate-length system; Sirius does not infer metres or solar masses. Alcubierre `sigma` has inverse-coordinate-length units, so increasing it makes the wall thinner, not thicker. `observer.distance` is a coordinate radius, not a dimensionless ratio: mass-bearing scenes require `5M <= r <= 1000M`; the exact isotropic one-sheet Ellis profile requires `0.1 <= b0 <= 1000` and `5b0 <= rho <= 1000b0`; and Alcubierre requires `5L <= r <= 1000L` for `L=max(R,1/sigma)`, which also places the far boundary in the asymptotic region. Positive-$\Lambda$ scenes additionally require `r <= 0.99*r_c`, intersecting the ordinary interval with the exact cosmological causal patch. The governed Alcubierre profile requires `0.1 <= sigma*R <= 100`, keeping the radius and inverse-wall scales within the fp32-resolved envelope. Minkowski, de Sitter, Morris-Thorne, and Alcubierre have no mass parameter and require `mass=0`; Minkowski and de Sitter otherwise retain the 5-to-1000 coordinate-radius interval. Changing the metric identity without supplying mass preserves a compatible lower-layer value or supplies the compatible default. Non-default throat/topology values apply only to Morris-Thorne, non-default warp/bubble values apply only to Alcubierre, and any explicit irrelevant parameter declines rather than being ignored.

## Backends

The CPU path is the reference implementation and runs without Vulkan or OpenGL. It supports pinhole, thin-lens, and fisheye projection plus the physical thin-disk polarisation output. Thin-lens rays preserve the requested field of view at the aperture centre; because the scale-free spacetime has no physical M-to-millimetre conversion, 50 mm-equivalent defines one virtual lens unit for the f-number calculation. The GPU path is one Slang kernel source (`src/sirius/kernels/`) compiled to SPIR-V and dispatched through Vulkan; the same source is also compiled to CUDA and Metal as a portability gate, not as a CUDA or Metal runtime backend. A memory governor sizes tiles against the host-visible coherent heap the adapter actually allocates from, and a dispatch governor bounds device occupancy. `SIRIUS_VULKAN_DEVICE=<zero-based-index>` selects one enumerated device strictly; the readiness report, renderer, tests, and hardware attestation use the same identity. Backend `auto` selects Vulkan only when a device, metric, and requested scene semantics are all represented; otherwise it logs the reason and selects CPU. Vulkan represents exact arbitrary sample counts, camera aberration, pinhole and finite-aperture thin lenses, Page-Thorne and the explicitly bounded zero-torque Shakura-Sunyaev temperature substitute, the explicit ZAMO Doppler-suppression diagnostic, dual-vector full-Riemann ray bundles, ellipse-filtered 100,000-star catalogues, and the bounded phenomenological grey volume with optional procedural density modulation. Volume marching is limited to 128 midpoint samples per accepted geodesic segment; larger explicit requests decline and `auto` uses CPU. Corona, narrowband line transfer, temporal motion blur, diagnostic colour, and polarisation are not represented on Vulkan; requests decline rather than borrowing another model. The precision ladder offers plain fp32, `SIRIUS_PRECISION=fp32-comp`, and `SIRIUS_PRECISION=fp64` on devices exposing `shaderFloat64`.

## Building

Requirements: CMake 3.28+, Python 3.10+, and a compiler with a C++26 mode: GCC 14+, Clang 17+, AppleClang 16+, or MSVC 19.40+ (Visual Studio 2022 17.10). First-party targets use C++26; their isolated vendored C translation units are pinned to portable C17. The Vulkan backend additionally wants the Vulkan headers/loader and the Slang compiler (`slangc`, found at `/opt/slang` or `SLANG_ROOT`). The interactive window uses OpenGL and GLFW; when OpenGL development files are absent, or `SIRIUS_BUILD_VIEWER=OFF`, Sirius builds the complete CPU/Vulkan CLI and the `view` command declines explicitly. Portable builds likewise degrade clearly to a complete CPU system when Vulkan is absent; the `linux-ci` operational profile sets `SIRIUS_REQUIRE_VULKAN_RUNTIME=ON` and rejects configuration, missing kernels, runtime skips, or failed dispatch.

```bash
# Linux / WSL2
cmake --preset linux-gcc && cmake --build --preset linux-gcc
./bin/linux-gcc/src/sirius/app/sirius

# Windows
cmake --preset windows-msvc && cmake --build --preset windows-msvc
```

Presets exist for GCC and Clang (Release and Debug), MSVC, macOS, required Linux Vulkan CI, and GCC ASan/UBSan/LSan. Every supported preset pins warnings as errors, enforce-mode contracts, Mandatory tests, and source governance; `linux-ci` additionally runs the non-skipping Vulkan/P1 operational profile. Source governance also enforces one CMake owner per first-party translation unit, downward-only layer dependencies, boundary naming, and the closed live-vendor set. `cmake --install` creates a relocatable volume containing the binary, starfield, installed operating model, every compiled backend resource, and viewer shaders when the viewer is compiled, then verifies that exact capability-specific structure before succeeding. A release install additionally requires the deterministic Mandatory-gate receipt emitted only after live CTest registration exactly equals a zero-skip JUnit estate; the receipt binds the revision, model, alignment authority, all test executables, and every installed product hash.

The application test executable is structurally non-rendering: it owns parsing, validation, projection, and input-state checks only. Render-session entry, CLI dispatch, CPU/Vulkan output, and progressive-frame cases live in the rendering test executable and carry an explicit `Rendering` label. Both executables remain beside the candidate so strict modes inspect that same product volume without a development resource override; source governance rejects relocation, an unclassified app-side command execution, or restoration of the override.

Remote CI Actions and fetched source dependencies are pinned to immutable commits. Direct qualification-tool downloads use versioned URLs and checked-in SHA-256 values; source governance rejects movable tags, latest-release selection, or unchecked direct downloads.

Normal presets use `SIRIUS_ALIGNMENT_MODE=development`: they compile and test
the system while reporting every absent external domain as pending, but their
artifacts are deliberately inadmissible as external evidence. Hardware, native
runtime, viewer, and native CI producers reconfigure the exact clean revision in
`qualification` mode. Qualification enforces the complete release-equivalent
build/product policy, the zero-skip Mandatory gate, executable-volume resource
locking, and candidate hashing while allowing external domains to remain
pending; it cannot package or report top-level readiness. Each record includes
the copied qualification executable, alignment receipt, product/test gate,
gate-generated JUnit and log, and an independently rerun JUnit/inventory pair.
The verifier cross-checks their byte hashes and exact test identities before the
record can become release input. A distributable package is available only from
`release` mode. That mode requires
a clean Git revision plus verified evidence for the physical Radeon, WSL2/Dozen,
native Windows build and Vulkan, native macOS build and MoltenVK, native viewer
input, and exact IMAX domain. Configure generates one deterministic receipt;
every build revalidates it, the compiler embeds it, installation carries it,
and render/view initialisation rejects a missing or altered copy:

```bash
cmake --preset linux-gcc \
  -DSIRIUS_ALIGNMENT_MODE=release \
  -DSIRIUS_REQUIRE_VULKAN_RUNTIME=ON \
  -DSIRIUS_ATTESTATION_ROOT=/absolute/path/to/attestation-bundles
cmake --build --preset linux-gcc
```

Qualification and release policy require tests and the Mandatory build gate, enforce-mode
contracts, warnings-as-errors, the sole `Release` configuration, compiled
Vulkan/Slang kernels validated by `spirv-val`, and the native viewer. Configure
fails rather than packaging a reduced product. Installation and CPack fail if
the complete gate was not run after the last artifact or registration change.
Release readiness, render, and viewer initialisation then parse that installed
gate receipt and rehash the running executable plus every installed product;
missing, skipped, stale, or same-size altered evidence fails closed. Development
builds do not require a gate receipt and remain local diagnostic builds; they
cannot generate admissible external evidence. `SIRIUS_RESOURCE_DIR` is
likewise development-only; qualification and release binaries ignore it and resolve the receipt,
executable, and resources exclusively from their packaged volume. To avoid a
self-proving gate, the Mandatory target removes any prior staged receipt before
CTest, requires pre-gate strict-mode install/readiness to decline, and stages only
the receipt issued after the exact zero-skip estate succeeds.
Release install/CPack then invokes the relocated binary's non-rendering
`info capabilities` path, so external byte verification and the product's own
runtime authority must agree before installation succeeds.

## Running

```bash
# A Kerr black hole with its disk, 256 samples per pixel
sirius render -m Kerr -a 0.9 -s 256 -o kerr.png

# The same scene as linear HDR for compositing
sirius render -m Kerr -a 0.9 -s 256 -o kerr.exr

# Transported thermal-disk polarisation (CPU selected automatically)
sirius render -m Kerr -a 0.9 --color-mode Polarisation -o kerr-polarisation.png

# Explicitly require Vulkan with the bounded grey volume and procedural density modulation
sirius render -m Kerr -a 0.9 -s 4 --backend vulkan --volumetric \
  --turbulence -o kerr_vk.png

# Interactive progressive viewer, machine-readable system report
sirius view --spin 0.9 --backend auto
sirius info system --json
sirius info readiness --json
sirius info capabilities --json
```

`info readiness` reserves top-level `ready`/`release_ready` for the strict
aligned-system claim. Development and qualification builds with valid resources
but pending external domains return non-zero with `ready: false` and
`evidence_generation_ready: true`; only the qualification artifacts emitted by
hardware/native runbooks are admissible. The `ultimate_ideal` object also reports the exact
sorted admitted, pending, and required domain IDs, so an operator can act on a
blocked count without consulting the source tree. A release build can report
ready only after all eight revision-bound domains have been admitted.

Configuration layers in a fixed order: struct defaults, then a JSON config file, then `SIRIUS_*` environment variables, then command-line flags; later layers win, every parameter is validated at startup against the ranges the physics supports, and invalid configuration stops the run rather than being clamped. `sirius render --help` lists the full flag set.

`render_test.sh` and `render_demo.sh` select only the exact preset binary (or an explicit `SIRIUS_BINARY`), propagate renderer failures, and require every declared output to exist and be non-empty. They never search for an arbitrary stale build.

The DNGR-technique features are explicit flags, each default-off so the pinned reference render is stable: `--beams` propagates two full-Riemann Jacobi vectors and derives the oriented beam ellipse; `--starfield point` replaces the background texture with a spatially indexed 100,000-star point catalogue filtered through both ellipse axes, using a display-calibrated relative-flux zero point so the catalogue survives PNG quantisation; `--volumetric` selects the declared stationary, vertically isothermal Gaussian grey atmosphere, midpoint-samples every accepted integration segment without using endpoint membership as an intersection proxy, and preserves every finite positive optically thin layer; `--turbulence` adds explicitly procedural density modulation; `--doppler-beaming off` substitutes the locally non-rotating ZAMO frequency at the actual disk latitude as a nonphysical diagnostic while retaining gravitational and frame-dragging transfer; `--color-mode Polarisation` selects transported CPU Stokes output; and a camera velocity applies special-relativistic aberration across all lens models. Corona, narrowband lines, temporal motion blur, two-sheet wormhole traversal, and every incomplete polarisation combination have named fail-closed boundaries.

## Testing

```bash
ctest --test-dir bin/linux-gcc -j"$(nproc)"          # everything
ctest --test-dir bin/linux-gcc -L Mandatory          # the build gate

# Instrumented mandatory gate
cmake --preset linux-gcc-sanitize
cmake --build --preset linux-gcc-sanitize -j"$(nproc)"
```

Labels carry gating semantics: Mandatory and Operational suites fail every normal build by default; Correctness suites cover additional feature behaviour; Performance suites measure without gating. New suites without an explicit policy are a hard generator error. Parameterised/typed GoogleTests are also rejected until the generator can prove their exact discovered names, preventing tests from escaping the Mandatory label. `tests/operating_model.json` maps all ten required P1–P6/E1–E4 acceptance criteria, 24 operating dimensions, and 29 explicit capability contracts to Mandatory evidence. Seven criteria are source/build gated; P3, P5, and E3 remain explicitly attestation-required because present software evidence cannot prove physical IMAX/780M operation. The build rejects a missing state, renamed witness, de-gated claim, or incomplete revision-bound release receipt; the same authorities are installed and available through `info capabilities` and `info readiness`.

External native/runtime evidence additionally carries CTest's JSON inventory:
its non-skipping JUnit names must equal every enabled registration and meet the
operating model's conservative 700-case source-available floor. A green subset
cannot be relabelled as the complete estate.

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
