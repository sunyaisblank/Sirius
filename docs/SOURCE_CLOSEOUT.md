# Renderer validation and integration status

This records scoped renderer validation and remaining qualification limits.
Development and retention instructions are in [DEVELOPMENT.md](DEVELOPMENT.md).
Earlier experiment details remain in Git history, including `9e20150` and
`9661999`; historical results do not qualify later revisions.

## Completed local validation

Renderer revision `cdcc254` passed all 934 selected tests, with zero failures
or skips, in 5,923.72 seconds. The selected identities, clean revision and seven
executable hashes were checked before and after the run. All explicit
application and test targets built successfully. The 97 Vulkan-related cases
were excluded at the user's request; the inventory contained 1,031 tests.

The complete selection covered Kerr PPM/PNG/EXR frames, Doppler and polarisation
comparisons, the CPU viewer, both analytic shadow references, coupled transport,
application behavior, installation/relocation and source/build authorities.
The moving ThinLens Kerr point detector produced identical linear radiance with
one and two workers in 1,210.11 seconds. The original 256-ray throughput test
passed in 46.52 seconds against its unchanged 60-second limit.

The exclusion covered `VulkanBackend`, `KernelParity`, `RetainedCameraProgram`,
`RetainedComputeTest`, `KernelInfinityDevice`, `KernelTrace`, `KernelBeam`, and
`VulkanRenderSession`, plus
`RenderCommandParse.ExplicitGpuRequestRunsVulkanWhenDevicePresent` and
`ViewCommandOperational.VulkanRefinementPublishesProgressiveFrames`.
CPU-only retained arithmetic, kernel portability emission, Vulkan configuration
parsing and mocked cancellation remained included. This result cannot issue the
complete Mandatory receipt or qualify a later source revision.

## Implemented numerical path

The CPU camera differentiates physical film and pupil coordinates through the
metric frame at the displaced observer. Four columns enter joint Hamiltonian
transport. Analytic Kerr-family Hessians share the metric authorities, and
retained twofold arithmetic preserves critical rays and covariant variations.
Original launch packets and independently generated high-precision metric,
transport and projection fixtures remain in `tests/support/`.

The physical point-source detector uses the live continuous camera and tracer,
preserves packets through subdivision, retains multiple images, and withholds
partial radiance on failure or exhaustion. CPU cancellation reaches private
trace intervals before image publication. ThinLens, endpoint, projection,
shadow-bracketing and gauge regressions preserve their physical samples,
numerical tolerances and work limits.

Kerr-family Vulkan frames use bounded retained camera, initialization, joint
transport, projection and dense stages. The host tracer owns event localization
and source/detector evaluation. Accepted phase records survive continuation and
rollback. Cooperative workgroups, fixed scratch buffers, batched requests and
bounded dispatch adaptation reduce device overhead. These numerical working
representations and their arithmetic radii are not global ODE certificates.

## Performance and qualification limits

A complete moving ThinLens Kerr/100,000-star CPU–device comparison at renderer
revision `bf30294` measured relative linear-RGB L1 error
`8.5437703874154113e-9` against the unchanged `0.02` limit. Its eight device
pixels took 11,228.8 seconds and 1,338,417 retained dispatches on Radeon 780M
through WSL2/Dozen. The maximum measured submission was 469.421 ms, with seven
soft subdivisions and no safety fallback. This proves that focused image
agreement, not practical full-workload throughput or final-source qualification.

The larger 192×128, three-sample moving ThinLens catalogue scene was interrupted
at `c0634c7` after about 113 minutes with no completed pixels. It had executed
2,808,742 accepted device intervals, 7,443 phase initializations and no rejected
intervals or safety reductions. A regular detector packet alone evaluates 601
coordinates: across 73,728 camera samples that implies about 44.3 million probes
before nonlinear refinement and image-root searches. That observation motivated shared detector sampling and bounded CPU/device
work ownership; it remains the historical device result for that workload.

Subsequent development hands outward vacuum detector rays to the checked Kerr
infinity map earlier (`df95428`). Current point-source discovery shares canonical
32×32 regions at each original pupil, bisects declined regions, and falls back to
original scalar footprints within a bounded extra-work allowance. Every pixel
retains its Gaussian, source map, frequency and image transmission. Workers cache
only completed region RGB. Sky-only point scenes need no additional centre ray;
scenes with disk or volume emission retain their foreground traces. Device jobs
own whole regions so separate pixel workers do not repeat their discovery.
CPU tiles smaller than a region now form one worker-owned scheduling group,
retaining their requested bounds and separate publication. Each region is
evaluated once, including partial edges. Tile completion uses its direct ID
index instead of scanning the full tile list.

Independent detector cell probes now execute through a bounded ray queue. Each
worker owns its tracer; root refinements and foreground rays use that same pool,
so coordinator threads do not add competing trace work. Probe coordinates,
cache identity and result order remain unchanged. Retained-device capacity now
accounts for up to 13 independent probes within each region. This provides work
to batch even when an image contains only one region; the existing device memory
and dispatch limits still apply.

A 32×32 region of the moving ThinLens Kerr scene used 607 shared probes for 1,024
footprints, collected in 87 sampler calls of at most 13 probes. Its inner and
infinity-tail work remained 137,580 and 44,454 attempts. Its two nonzero reference
pixels agreed with separately evaluated original footprints to maximum relative
RGB difference `1.5695371152294026e-10`, unchanged from `889e11f`. That physical
check took 141.24 seconds including the independent scalar evaluations. It
measures available probe concurrency, not GPU dispatch throughput. The analytic
1,024-footprint check uses 792 probes and matches its independent Gaussian flux
oracle within 2e-6 relative.

After grouping small tiles, the same 32×16, one-sample, 100,000-star
serial/two-worker comparison passed in 83.92 seconds with identical linear pixels
and unchanged requested tile counts (106.75 seconds at `78886be`). The 33×1
edge/three-sample comparison likewise retained identical pixels in 30.03 seconds
(previously 35.00). These are local comparison-test wall times, not portable
throughput guarantees. Eighteen focused scheduler, queue and session checks
passed. The explicit Linux GCC render-test target built with warnings as errors;
format, source ownership and the 1,059-case live CTest inventory checks passed.

The full 192×128, three-sample CPU workload completed at `9661999` in
386.29 seconds (63.62 pixels/second), using 15 workers and all 24 requested
32×32 tiles. It used the scene in
[`moving_kerr_detector_scene.h`](../tests/support/moving_kerr_detector_scene.h):
Kerr mass 1 and spin 0.7, observer radius 50 and inclination 1.5708 radians,
2-degree field of view, camera beta (0.1, 0.8, 0), and a 50 mm-equivalent thin
lens at f/2.8 focused at 50 coordinate units. Disk, bloom, film and ray-bundle
flags were off. The catalogue contained 100,000 stars, seed 42, distances
1–10,000 parsecs and brightness 1e-5. The app's JSON configuration preserved
these controls; its typed scene transcript confirmed them. The completed linear
EXR decoded to finite, nonnegative RGB with 88 nonzero pixels and maximum
channel value `3.2453112908115145e-7`. This measures CPU workload completion,
not full-image convergence or current GPU throughput.

The linear output is retained locally under `renders/moving-kerr-192x128/`.
Its SHA-256 is `a65b191d83684a0958bb39d7f681b4d310c0712ed1897a171ca2ed9b41b80bb4`.
The PNG beside it is an inspection preview with linear values multiplied by
1e7 before the normal display pipeline; it is not a second rendered exposure.
The producing executable hash was
`2029e0b98b124bdf1ca6f8a1c29287ee4c8b302a9030416243cb275d6d340f85`.

That workload exposed an ETA error: the first 15 completed tiles arrived as a
burst after minutes of tracing, causing the old short-window rate to predict
two seconds remaining. Progress now uses monotonic total elapsed throughput,
with deterministic burst, quiet-period and restart checks. Scheduling and
radiance are unchanged. Configuration, catalogue-domain, scene transcript,
Gaussian normalization and disk-intersection checks passed (64 GoogleTests),
as did the attestation verifier's missing/mismatched catalogue controls.
The two deterministic progress checks, strict Linux GCC render-test build,
format/source checks and current 1,065-case live CTest inventory also passed.

PR integration found three native portability issues: MSVC could not prove
disk-event locals initialized, AppleClang rejected `constexpr std::abs`, and
Apple's standard library lacked the Doppler test's `std::jthread`. The fixes
preserve the disk calculations, Gaussian formula, original raster reduction and
independent per-worker tracers. Native integration runs compile all product and
test targets and execute non-render authority controls; they do not issue the
complete runtime or scientific receipt.

No external domain was admitted in the closeout build (0/8). Physical Radeon,
WSL2/Dozen, native Windows/macOS build and runtime, native viewer input and the
exact IMAX workload retain their independent qualification requirements. Release
packaging remains disabled. The exact revision's `mandatory_gate.json` and
verified external bundles, rather than historical logs, establish qualification.

## Recovery and repository state

The canonical checkout preserved the recovered source on
`development/source-closeout-2026-09-13`, based on `main` at `8c7ba5a`.
Forty-two sibling Sirius/test-toolchain directories were removed after their
commits and uncommitted source were preserved. Redundant local topic branches
and the already-archived stash were removed. Recovery refs remain under
`refs/archive/`; a verified bundle and recovery records remain in
`.git/closeout/`. They are not build inputs or qualification evidence.
