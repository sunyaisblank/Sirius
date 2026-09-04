# Sirius external attestation domains

Local CPU and software-Vulkan evidence cannot close a native platform,
physical-device, native-window, or full-frame hardware claim. These domains
are `attestation_required` in `tests/operating_model.json`; an absent record is
`UNEXECUTED`, not a failure and not a pass.

Every admitted record uses schema version 1 and passes:

```bash
python3 scripts/verify-attestation.py \
  --source-root /path/to/clean/exact-revision-checkout \
  path/to/attestation.json
```

The verifier requires a passing domain-appropriate test estate, source revision,
UTC completion time, exact domain identity, and byte count plus SHA-256 for every
artifact. Standalone verification derives both the clean Git revision and model
identity from `--source-root`; omitting the option selects the verifier's own
checkout. It rejects an attested revision that differs from that selected
worktree. It independently hashes the selected checkout's canonical
`tests/operating_model.json` and rejects a bundle whose gated operating-model
product differs, even when the bundle's receipt and gate are internally
self-consistent. Aggregate admission passes that independently derived digest
into every record verification before admitting its domains; qualification and
release admission also require `--source-root` and refuse a selected model that
differs from that exact Git worktree. Native-build records bind the full
governed CTest registration, the
exact nine-case non-render authority estate, every compiled test executable, all
nine products, the alignment receipt, and the copied candidate. Runtime
records must additionally identify the physical device, driver, Vulkan API,
memory budget, and fp64 capability, and must hash the rendered PNG, execution
transcript, device inventory, configure-time alignment receipt, qualification
build-gate receipt, copied qualification executable, gate-generated JUnit/log,
and non-skipping independently rerun JUnit report plus CTest's machine-readable
registration inventory. Development and release artifacts are rejected as
evidence producers: the receipt and gate must both name `qualification`, the
source must be clean, and qualification must have enforced the same complete
build/product topology as release without claiming that pending external domains
were already satisfied. The local
operating-model gate independently parses the source estate and rejects fewer
than 700 enabled GoogleTests. This governed floor is not a historical fixed
estate count. Build-domain evidence must register the full governed estate and
execute the fixed non-render authority selection; runtime-domain JUnit case
names must equal every enabled registered test and meet that same floor. Both
contain the compiled-receipt authority witness;
the receipt must name the same clean revision and exact model-derived domains.
The gate-generated and rerun JUnit case-name sets and summaries must agree; the
gate must bind the submitted alignment receipt, operating-model digest, all
seven tested executables, all nine candidate products, and the actual copied
candidate bytes admitted as `qualification-sirius.bin`.
The report, inventory, revision, selected-device readiness, and required profile
are cross-checked rather than trusted from result booleans. Artifact paths must be
unique and remain inside the record bundle. Its Mandatory self-test proves that
software devices, missing hardware identity or upstream evidence, false test
state, wrong 1080p or IMAX dimensions, ambiguous or unbound revisions, aliased
artifacts, tampered hashes, and substituted models or revisions are rejected.
Repository source governance separately rejects a detached source-root CLI.

## Zero-render preflight

Before any expensive test estate, renderer, or native viewer is started,
`scripts/attestation_preflight.py` can prove that the exact qualification
candidate, compiled alignment receipt, governed resources, executing host,
selected physical Vulkan device, driver route, free output space, and required
operator tools agree. The script is non-promoting by construction: it executes
only `info system` and `info readiness`, emits `planned_domains` rather than
`domains`, carries no admissible artifact set, states `promotable=false` and
`external_execution_completed=false`, and refuses to write inside the clean
source worktree. A successful preflight is necessary preparation, never runtime
or image evidence.

For the consolidated WSL2/Radeon/Dozen/IMAX/viewer route, compile the candidate
without invoking the Mandatory target and keep the report outside the source:

```bash
cmake --preset linux-ci \
  -DSIRIUS_ALIGNMENT_MODE=qualification \
  -DSIRIUS_REQUIRE_VULKAN_RUNTIME=ON
cmake --build --preset linux-ci --target SiriusAlignmentGate sirius
VK_DRIVER_FILES=/absolute/path/to/dzn_icd.x86_64.json \
  python3 scripts/attestation_preflight.py \
    --profile wsl2-radeon-viewer \
    --candidate bin/linux-ci/src/sirius/app/sirius \
    --expected-revision "$(git rev-parse HEAD)" \
    --output /tmp/sirius-wsl2-preflight.json
```

The profile requires the real executing host to be WSL2 with `/dev/dxg`, an
unambiguous integrated/discrete Radeon 780M reporting Dozen/dzn, at least the
governed 2048 MiB render budget, a compiled viewer, `DISPLAY`, `timeout`, and
`xdotool`. Native hosts use `--profile windows-native-vulkan` or
`--profile macos-moltenvk`; they respectively reject Dozen and require
MoltenVK. Omitting `--output` prints the same report to standard output.

## Physical Radeon and WSL2/Dozen

`scripts/validate-hardware.sh` selects one exact non-software Vulkan device and
uses that selection for inventory, readiness, the complete test estate, the
1920x1080 and 5616x4096 governed renders, and every precision rung. Both scenes
disable disk emission so their structure must come from the ray-bundle-filtered,
display-calibrated 100,000-star point catalogue; it combines that catalogue with
a nonzero timelike camera velocity and the finite-aperture thin lens. Film
jitter and pupil position come from separate deterministic camera-sample
dimensions; reusing the film pair for the pupil cannot satisfy the build gate.
Mandatory host/device projection controls first require that this lens moves the metric
launch event across its bounded rest-frame pupil and reconverges at the focus
plane; they do not replace the physical frame. The verifier
requires those P3/P5 semantics in the record, fully decodes both governed PNGs,
and admits only structured RGB images at 1920 by 1080 and 5616 by 4096 under a
2048 MiB budget. It also requires one unambiguous hashed typed-session scene
event at each resolution to agree with every claim,
the actual 100,000-star catalogue to reach Vulkan dispatch with beams enabled,
the logged device and memory budget to match the record, and that same render to
write the admitted PNG before reaching `Complete`. A hashed JUnit report must
contain the complete green estate with no skips, agree with the CTest transcript,
and come from the `linux-ci` profile; the hashed inventory must select the same
device named by readiness and rendering. This rejects dimensionally valid but
visually collapsed output as well as metadata that describes a scene or test
state the runbook did not execute. The runbook refuses a dirty source tree so
the attested full commit revision identifies the code that was built.

Native Linux/RADV:

```bash
SIRIUS_ATTESTATION_DOMAINS=physical-radeon-780m,physical-imax-5616x4096 \
  scripts/validate-hardware.sh
```

WSL2 through Dozen:

```bash
SIRIUS_ATTESTATION_DOMAINS=physical-radeon-780m,wsl2-dozen,physical-imax-5616x4096 \
  VK_DRIVER_FILES=/absolute/path/to/dzn_icd.json \
  scripts/validate-hardware.sh
```

The Dozen domain additionally requires WSL2 to be reported and the selected
device's Vulkan driver name or information to contain `Dozen` or `dzn`.
Lavapipe is rejected because its device kind is software.

## Windows and macOS

On explicit workflow dispatches and protected-branch pushes, the native CI
matrix can emit downloadable `windows-native-build` and `macos-native-build`
records after a clean strict qualification compilation and the exact non-render
build-authority estate succeed. The producer derives the recorded commit
from the tested Git checkout, requires that checkout to remain clean after the
tests, and cross-checks it against the workflow revision; caller-provided text
cannot name a different source. Each record also carries the configure-time
alignment receipt, non-promoting native-build gate, gate JUnit/log, all seven
compiled executables, all nine products, copied candidate, and current CTest
registration. The verifier requires at least 700 governed registered tests and
the exact nine zero-skip authority cases; it rejects any claim that this proves
runtime or the unexecuted estate. These records do not close native Vulkan.

Pull requests and explicit integration dispatches use the same non-render Linux,
Windows, and macOS jobs. Pull requests compile and exercise the authority
controls but cannot publish evidence. An explicit dispatch may additionally
issue only the Windows/macOS compilation domains through the separate
native-build gate. Neither path creates a Mandatory promotion receipt;
physical/runtime attestations remain push-only or external-runbook evidence.

`windows-native-vulkan` separately requires a native Windows host, a physical
Vulkan device, readiness, and runtime artifacts; Dozen is expressly rejected
for that domain. `macos-moltenvk` separately requires a native macOS host,
physical device evidence, and a driver identity containing `MoltenVK`. Both
runtime domains must also bind a non-skipping JUnit estate, the selected-device
inventory, exact source revision, readiness marker, and their required native
preset in hashed evidence. A PNG and a success boolean alone are rejected.

Run the shared producer from a clean checkout on the native host:

```text
# Windows PowerShell
python scripts/validate-native-runtime.py --expected-revision $env:GITHUB_SHA

# macOS shell
python3 scripts/validate-native-runtime.py --expected-revision "$GITHUB_SHA"
```

When more than one physical Vulkan device is present, pass
`--device-pattern` with an unambiguous substring. The producer chooses only one
physical device, rejects Dozen on native Windows and non-MoltenVK drivers on
macOS, binds that selection through inventory/readiness/tests/a small native
Vulkan probe, rechecks the clean revision, and verifies its own record before
reporting success.

## Native viewer input

Handler semantics for press, repeat, release, mouse drag, scroll, finite input,
and progressive restart are Mandatory-gated without a window. The
`viewer-native-window-input` domain remains separate: a record must show that a
GLFW/OpenGL window was created on the named physical Radeon target, that a
progressive Vulkan frame was published, and that native keyboard and pointer
callbacks were both observed. The verifier independently binds the selected
device, readiness marker, complete non-skipping JUnit report, inventory,
execution log, and callback/frame transcript by hash. Direct calls to handlers,
a CPU-only frame, software Vulkan, or unsupported booleans do not satisfy this
domain.

On an X11 or XWayland desktop, install `xdotool` and run:

```bash
scripts/validate-viewer.sh
```

When the same clean WSL2 checkout has just completed the consolidated Radeon,
Dozen, and IMAX hardware run, avoid rebuilding and rerunning the full estate:

```bash
SIRIUS_REUSE_QUALIFICATION_ATTESTATION=/absolute/path/to/hardware-attestation.json \
  scripts/validate-viewer.sh
```

`reuse_qualification_evidence.py` first runs the normal attestation verifier on
that hardware record, then requires its exact source revision, WSL2 platform,
complete Radeon 780M/Dozen device object, qualification candidate bytes,
alignment/Mandatory receipts, CTest registration/JUnit, and all governed product
hashes to equal the live viewer volume. Only those authority artifacts and the
verified upstream transcript are copied; stale, circular viewer, native-Linux,
wrong-driver, cross-device, or altered-resource reuse fails before a window is
created. The viewer still must publish a new Vulkan frame and receive new native
keyboard and pointer callbacks. Thus one physical qualification estate can
support both records without reusing either domain's physical observation.

The runbook refuses a dirty source tree, selects and verifies the named physical
Vulkan device, builds and runs the complete test estate or reuses the exact
verified WSL2 hardware authority described above, opens the real viewer,
waits for a published Vulkan refinement frame, delivers keyboard, pointer-drag,
and scroll events through the host window system, then exits through an Escape
event. The viewer writes the transcript only when
`SIRIUS_VIEWER_INPUT_TRANSCRIPT` is explicitly set by the runbook.

## Record disposition

Attestations are evidence for their exact source revision and named domain
only. A workflow definition, historical report, record from another commit,
software Vulkan run, or successful memory-plan unit test cannot be promoted
into a physical/native pass.

The current admitted-domain count is external state and is deliberately not
hardcoded in source documentation. Editing source to record that count creates
a new source revision for which the records do not apply. Operators obtain the
live admitted and pending partition only by running `verify-alignment.py` over
the independently verified, same-revision bundles; its deterministic receipt
is the authority consumed by release configuration and the installed runtime.

An individual passing record also does not make a release aligned. Sirius
admits a set only when all eight domains are covered exactly once and every
record passes this verifier for the same clean source revision. The release
configuration command is:

```bash
cmake --preset linux-gcc \
  -DSIRIUS_ALIGNMENT_MODE=release \
  -DSIRIUS_REQUIRE_VULKAN_RUNTIME=ON \
  -DSIRIUS_ATTESTATION_ROOT=/absolute/path/to/attestation-bundles
```

`scripts/verify-alignment.py` derives those eight domains from the exact
attestation-required capability profiles and binds the operating-model SHA-256
into the deterministic admission receipt. The governed model is checked out as
LF-terminated UTF-8 on every host, qualification rejects any platform-local
CRLF model identity, and the receipt is written as explicit UTF-8/LF bytes;
therefore the same revision, model, and admitted record set has one byte-level
authority across Linux, Windows, and macOS.
Configuration fails on missing, duplicate, dirty, wrong-revision, or
non-qualification evidence;
it also refuses disabled tests, a weakened Mandatory/Werror/contract policy,
non-Release configurations, missing Vulkan/Slang or `spirv-val`, and a product
without the native viewer. The build rechecks it, the executable embeds it,
and the Mandatory target writes a deterministic build-gate receipt only after
live CTest registration exactly equals a zero-failure, zero-error, zero-skip
JUnit estate. It removes both prior canonical and staged receipts before CTest;
the release pre-gate controls require install/readiness to decline, and the new
receipt is staged only after the estate passes. That receipt hashes every test executable and release product;
install and CPack reject a missing or stale receipt and rehash the installed
volume, then execute the relocated binary's non-rendering `info capabilities`
initialisation and require it to report the aligned schema. Runtime readiness plus render/view initialisation compare the installed
alignment receipt with the compiled authority, parse the installed Mandatory
receipt, and rehash the running executable plus every installed product. A
missing receipt, any skipped case, or even a same-size product alteration fails
closed in release mode. Release lookup also ignores `SIRIUS_RESOURCE_DIR`, so an
environment-selected receipt/resource tree cannot replace the packaged volume.
Development mode retains the override and does not require a gate receipt,
records pending domains explicitly, and cannot produce a release package or an
admissible attestation. Qualification mode is the sole evidence-production
authority: it requires a clean source, strict release-equivalent product, and
fresh Mandatory gate while allowing pending domains and disabling packaging.
It exposes `evidence_generation_ready` while top-level `ready` remains false and
returns non-zero until the strict ideal is aligned; evidence generation is
therefore possible without either a circular release build or a development
artifact being promoted as release evidence.
