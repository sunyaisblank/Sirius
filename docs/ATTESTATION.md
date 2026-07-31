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
memory budget, and fp64 capability, and must hash both a rendered PNG and the
execution transcript. Artifact paths must be unique and remain inside the
record bundle. Its Mandatory self-test proves that software devices, missing
hardware identity, false test state, wrong IMAX dimensions, ambiguous
revisions, aliased artifacts, and tampered hashes are rejected.

## Physical Radeon and WSL2/Dozen

`scripts/validate-hardware.sh` selects one exact non-software Vulkan device and
uses that selection for inventory, readiness, the complete test estate, the
5616x4096 governed render, and every precision rung. It decodes the governed
PNG header and admits only 5616 by 4096 under a 2048 MiB budget. It refuses a
dirty source tree so the attested full commit revision identifies the code that
was built.

Native Linux/RADV:

```bash
SIRIUS_ATTESTATION_DOMAINS=physical-radeon-780m,physical-imax-5616x4096 \
  scripts/validate-hardware.sh
```

WSL2 through Dozen:

```bash
SIRIUS_ATTESTATION_DOMAINS=physical-radeon-780m,wsl2-dozen,physical-imax-5616x4096 \
  VK_ICD_FILENAMES=/absolute/path/to/dzn_icd.json \
  scripts/validate-hardware.sh
```

The Dozen domain additionally requires WSL2 to be reported and the selected
device's Vulkan driver name or information to contain `Dozen` or `dzn`.
Lavapipe is rejected because its device kind is software.

## Windows and macOS

The CI matrix emits downloadable `windows-native-build` and
`macos-native-build` records only after native compilation and the complete
source-available CTest estate succeed. These build records do not close native
Vulkan.

`windows-native-vulkan` separately requires a native Windows host, a physical
Vulkan device, readiness, and runtime artifacts; Dozen is expressly rejected
for that domain. `macos-moltenvk` separately requires a native macOS host,
physical device evidence, and a driver identity containing `MoltenVK`.

## Native viewer input

Handler semantics for press, repeat, release, mouse drag, scroll, finite input,
and progressive restart are Mandatory-gated without a window. The
`viewer-native-window-input` domain remains separate: a record must show that a
GLFW/OpenGL window was created and that native keyboard and pointer callbacks
were both observed, with the input transcript included as a hashed artifact.
Direct calls to the handlers do not satisfy this domain.

## Record disposition

Attestations are evidence for their exact source revision and named domain
only. A workflow definition, historical report, record from another commit,
software Vulkan run, or successful memory-plan unit test cannot be promoted
into a physical/native pass.
