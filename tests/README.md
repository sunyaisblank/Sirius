# Test harness

Each layer owns one GoogleTest executable. `base` covers contracts and resource
primitives; `core` covers geometry, integration and physical models; `oracle`
contains independent mathematical references; `backend` covers live CPU/device
transport; `render` covers sessions, images and writers; `app` covers parsing,
configuration and input-state logic without rendering. `operational` covers
installation, evidence admission, build policy and operator workflows.

The backend executable always includes the CPU tracer, coupled transport and
source-map tests. Vulkan development files add device cases; Slang enables
their kernel inputs. CTest inventory governance checks those configured
capabilities explicitly: absent device tools cannot remove CPU coverage, and
an enabled device suite cannot silently lose tests.

The source authorities are `operating_model.json`, the policy in
`scripts/generate-ctest-labels.py`, and the operational registrations in
`tests/CMakeLists.txt`. `labels/CTestLabels.cmake` is generated: do not edit it
by hand. Reference generators and high-precision fixtures in `support/` are
test inputs, not disposable logs. Keep their independent derivations intact.
Repository source-contract checks and their negative controls live by subject
in `scripts/governance/`; `scripts/verify-repository-structure.py` remains the
single command used by CMake and CI.

## Choose validation by the change

| Change | First useful check | Broader check when warranted |
|---|---|---|
| Documentation or workflow instructions | Review commands and references | No renderer run required |
| Python governance or test registration | Relevant verifier/self-test | Live CTest inventory after configuration |
| CMake ownership or presets | Configure, inspect target graph and test registration | Compile affected targets when flags, dependencies or generated inputs change |
| CLI/configuration logic | `sirius_app_tests` or named cases | A render only when execution semantics change |
| Metric, arithmetic or integration | Affected core/backend cases and independent references | Complete local scientific selection for changes across numerical paths |
| Session, detector or image writer | Affected render cases | Complete relevant images after the implementation stabilizes |
| Driver, device, native input or release | Exact-domain runbook | Complete required qualification on the final clean revision |

For source-only changes these checks need no compiler or renderer:

```sh
python3 -B scripts/generate-ctest-labels.py --check
python3 -B scripts/verify-operating-model.py
python3 -B scripts/verify-repository-structure.py
python3 -B tests/operational/render_workflow_test.py
```

The operator workflow check uses real stub child processes on Linux, Windows,
and macOS. It checks exit status, missing/empty images, stale output refusal,
separate logs, explicit selection, and interrupted partial output cleanup.
It does not invoke Sirius. Sanitizer suppression lives in `tests/sanitizers/`; its
single Vulkan-loader allocation pattern does not suppress Sirius allocations.

Compile only the affected test target, inspect the selection, then execute it:

```sh
cmake --build --preset linux-gcc --target sirius_core_tests -j4
ctest --test-dir bin/linux-gcc -N -R '^CameraFilmDifferential\.'
ctest --test-dir bin/linux-gcc --output-on-failure --no-tests=error \
  -R '^CameraFilmDifferential\.'
```

Add `-C Release` for multi-configuration generators. Names are exact GoogleTest
identities. An unfiltered CTest or render/backend executable can perform costly
image work. `Rendering` is a useful selection label, not a runtime estimate:
individual transport and oracle cases can also be expensive. CPU-only work is
not necessarily short. Inspect test bodies before selecting unfamiliar suites.
Device consumers share a CTest resource lock, but independent CTest processes
do not share that lock; avoid overlapping device runs.

Do not lower image sizes, physical samples, numerical tolerances or work caps
just to shorten validation. Improve test scheduling when it preserves the
asserted behavior, and keep performance measurements separate from correctness.

## Complete gates and evidence

`cmake --build --preset <preset> --target RunMandatoryTests` runs the governed
estate and verifies exact inventory, zero failures/skips, source identity and
product/test hashes before issuing `generated/sirius/mandatory_gate.json`.
Normal builds include this target. Direct `ctest -L Mandatory` runs the selected
tests but does not issue that receipt. Focused tests cannot satisfy it.

`RunNativeBuildEvidence` binds a strict build and the fixed non-render authority
selection; it cannot stand in for runtime or image evidence. Qualification
and release modes retain their complete requirements. Development profiles can
skip unavailable Vulkan cases, so a green local selection is not a Vulkan pass.

The attestation authority control invokes only
`RenderEvidence.RetainedWireRecordsFeedAttestationControls` from the render test
executable. That case serializes synthetic scene/completion records without
rendering or opening a GPU. The verifier consumes those native records before
mutating source ownership, work coverage, stage counts, budgets and completion
ordering; a handwritten valid fixture alone cannot establish producer compatibility.

The last completed local renderer run is recorded in
[SOURCE_CLOSEOUT.md](../docs/SOURCE_CLOSEOUT.md), including its Vulkan exclusions
and exact tested revision. Reuse that historical evidence only for its stated
scope. Qualification receipts are revision-bound and cannot be reused after
changes by renaming them.
