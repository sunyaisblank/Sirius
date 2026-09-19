# Renderer source closeout — 14 September 2026

This records the completed local renderer validation from the previous session,
not completion of the subsequent repository consolidation or release
qualification. Detailed recovery and experiment history remains at commit
`9e20150` (`git show 9e20150:docs/SOURCE_CLOSEOUT.md`). Current development and
retention instructions are in [DEVELOPMENT.md](DEVELOPMENT.md).

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
before nonlinear refinement and image-root searches. Component arithmetic
speedups do not resolve that algorithmic workload. Detector sampling and
CPU/device work ownership require further performance work while preserving
physical acceptance criteria.

Subsequent development hands outward vacuum detector rays to the checked Kerr
infinity map earlier (`df95428`). Current point-source discovery shares canonical
32×32 regions at each original pupil, bisects declined regions, and falls back to
original scalar footprints within a bounded extra-work allowance. Every pixel
retains its Gaussian, source map, frequency and image transmission. Workers cache
only completed region RGB. Sky-only point scenes need no additional centre ray;
scenes with disk or volume emission retain their foreground traces. Device jobs
own whole regions so separate pixel workers do not repeat their discovery.

A 32×32 region of the moving ThinLens Kerr scene used 607 shared probes for 1,024
footprints. Its two nonzero reference pixels agreed with separately evaluated
original footprints to maximum relative RGB difference `1.5695371152294026e-10`.
That check took 134.57 seconds including the independent scalar evaluations.
The complete 32×16, one-sample, 100,000-star serial/two-worker frame comparison
passed in 107.65 seconds with identical linear pixels; the former 4×2 comparison
at `e1755eb` took 116.46 seconds. This is increased image resolution, not a timing
comparison of identical workloads. The analytic 1,024-footprint check used 792
probes and matched its independent Gaussian flux oracle within 2e-6 relative.
Twenty-two focused detector/session tests passed across the implementation
steps, including bounded fallback, disconnected visibility, folds, private
failure, cancellation, image edges, three-sample accumulation and untouched
linear EXR output. These are local source/CPU checks. The complete 192×128
three-sample workload and current device throughput remain unmeasured; GPU
probe scheduling and exact-domain qualification remain separate work.

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
