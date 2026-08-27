# Sirius external attestation domains

Local CPU and software-Vulkan evidence cannot close a native platform,
physical-device, native-window, or full-frame hardware claim. These domains
are `attestation_required` in `tests/operating_model.json`; an absent record is
`UNEXECUTED`, not a failure and not a pass.

Every admitted record uses schema version 1 and passes:

```bash
python3 scripts/verify-attestation.py path/to/attestation.json
```

The verifier requires a passing test estate, source revision, UTC completion
time, exact domain identity, and byte count plus SHA-256 for every artifact.
Runtime records must also identify the physical device, driver, Vulkan API,
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
estate count: the external JUnit case-name set must equal every enabled
registered test, meet that same conservative floor, and contain the
compiled-receipt authority witness;
the receipt must name the same clean revision and exact model-derived domains.
The gate-generated and rerun JUnit case-name sets and summaries must agree; the
gate must bind the submitted alignment receipt, operating-model digest, all
seven tested executables, all thirteen candidate products, and the actual copied
candidate bytes admitted as `qualification-sirius.bin`.
The report, inventory, revision, selected-device readiness, and required profile
are cross-checked rather than trusted from result booleans. Artifact paths must be
unique and remain inside the record bundle. Its Mandatory self-test proves that
software devices, missing hardware identity or upstream evidence, false test
state, wrong 1080p or IMAX dimensions, ambiguous or unbound revisions, aliased
artifacts, and tampered hashes are rejected.

## Physical Radeon and WSL2/Dozen

`scripts/validate-hardware.sh` selects one exact non-software Vulkan device and
uses that selection for inventory, readiness, the complete test estate, the
1920x1080 and 5616x4096 governed renders, and every precision rung. Both scenes
disables disk emission so its structure must come from the ray-bundle-filtered,
display-calibrated 100,000-star point catalogue; it combines that catalogue with
a nonzero timelike camera velocity and the finite-aperture thin lens. The verifier
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

On protected-branch pushes, the CI matrix emits downloadable
`windows-native-build` and `macos-native-build` records only after a clean
strict qualification compilation and the complete
source-available CTest estate succeed. The producer derives the recorded commit
from the tested Git checkout, requires that checkout to remain clean after the
tests, and cross-checks it against the workflow revision; caller-provided text
cannot name a different source. Each record also carries the configure-time
alignment receipt, candidate/product build gate, gate JUnit/log, copied tested
executable, and current CTest registration. It requires exact inventory/JUnit
equality and requires the JUnit witness that compared that receipt to the
compiled authority, preventing a partial, development-mode, or stale
report/binary tree from being relabelled as the current revision. These build
records do not close native Vulkan.

Pull requests and explicit integration dispatches use separate non-render
Linux, Windows, and macOS jobs. Each compiles the same product and test topology
and exercises governance/authority negative controls, but neither runs the
Mandatory estate nor creates, uploads, or admits an attestation. Repository
governance requires the full native jobs to remain push-only and requires every
integration job to prove that no promotion receipt appears.

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

The runbook refuses a dirty source tree, selects and verifies the named physical
Vulkan device, builds and runs the complete test estate, opens the real viewer,
waits for a published Vulkan refinement frame, delivers keyboard, pointer-drag,
and scroll events through the host window system, then exits through an Escape
event. The viewer writes the transcript only when
`SIRIUS_VIEWER_INPUT_TRANSCRIPT` is explicitly set by the runbook.

## Record disposition

Attestations are evidence for their exact source revision and named domain
only. A workflow definition, historical report, record from another commit,
software Vulkan run, or successful memory-plan unit test cannot be promoted
into a physical/native pass.

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
into the deterministic admission receipt.
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
