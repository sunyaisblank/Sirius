# Development workflow

Work in the canonical checkout. Use the existing CMake presets, with generated
build files under `bin/<preset>`. Do not create sibling development workspaces.

## Files and ownership

| Path | Contents | Retention |
|---|---|---|
| `src/sirius/` | Production layers and the independent mathematical oracle | Tracked source |
| `tests/` | Behavioral tests, reference generators/fixtures, operating model, sanitizer suppression | Tracked source |
| `cmake/` | Toolchain, dependencies, policy, alignment, test gates, installation | Tracked build definitions |
| `scripts/` | Source governance, generated-kernel tooling, render examples, qualification producers | Tracked tools |
| `assets/`, `lib/` | Runtime assets and the closed live vendor set | Tracked inputs |
| `bin/<preset>/` | Build cache, compiled products, generated fixtures and build-gate records | Ignored; keep reusable builds |
| `out/<task>/` | Disposable diagnostics, logs, inventories and temporary scripts | Ignored; remove after recording useful findings |
| `renders/<run>/` | Requested rendered images | Ignored; retain only useful outputs |
| `attestations/` | Verified qualification bundles | Ignored; keep separately from disposable diagnostics |

The two top-level shell render scripts are replaced by `scripts/render.py`.
The former `sanitizers/` directory is now `tests/sanitizers/`, consumed by both
the sanitizer test preset and the receipt-producing build gate. Historical
engagement reports and per-session chronology belong in Git history, not beside
current operating instructions. Recovery bundles already stored in
`.git/closeout/` are local recovery material, not build or qualification inputs.

## Build only what changed

Configure when changing CMake, dependencies, options, or presets, and when the
alignment check requests a new source-revision receipt. Compile the affected
explicit target, then select tests as described in [the harness guide](../tests/README.md).

```sh
cmake --preset linux-gcc
cmake --build --preset linux-gcc --target sirius_core_tests -j4
```

Targets include `sirius`, `sirius_base_tests`, `sirius_core_tests`,
`sirius_oracle_tests`, `sirius_backend_tests`, `sirius_render_tests`, and
`sirius_app_tests`. Building an explicit target compiles its dependencies and
discovers its GoogleTests, without running the complete gate. A normal build
without `--target` retains the complete Mandatory gate; use that at an
appropriate integration boundary, not after every edit. Explicit targets do
not issue qualification receipts.

The root CMake file sequences the layers. `sirius_build_policy.cmake` checks
configuration and product policy, `sirius_alignment.cmake` owns revision-bound
alignment, `sirius_testing.cmake` owns source/build gates, and
`sirius_install.cmake` owns installation and packaging. Each source layer owns
its translation units; `tests/CMakeLists.txt` owns test registration.

## Render examples deliberately

```sh
python3 scripts/render.py
python3 scripts/render.py --scene kerr --dry-run
python3 scripts/render.py --scene kerr --backend cpu
python3 scripts/render.py --preset windows-msvc --scene wormhole --dry-run
python3 scripts/render.py --scene kerr --scene schwarzschild --format exr \
  --width 512 --height 512 --samples 32
```

Use `python` instead of `python3` where that is the Python 3 executable.
No selection lists scenes and exits. `--all` explicitly selects the whole
catalogue. Defaults are 128×128, one sample per pixel, PNG, and CPU; these are
small operator examples, not scientific convergence or performance evidence.
Larger resolution and sample counts must be requested. Scene arguments still
participate in the application's normal configuration layering; parameters
they do not set retain the selected binary's configuration and environment.

The runner never configures, builds, discovers alternate binaries, or chooses
a different backend after an explicit backend failure. `--preset` selects one
build directory; `--config` selects its multi-configuration build type;
`--binary` selects an explicit executable. `--dry-run` prints argument arrays
without invoking the executable or creating files.

Each execution has a unique run directory. A case succeeds only if the child
exits zero and creates a nonempty image. Failure stops the batch, removes the
failed partial image, and retains completed images and diagnostic logs.
`out/render/<run>/run.json` records commands and completion state, but is not a
qualification manifest. A hard process kill may leave partial files; a run
still marked `running` is incomplete. Review or remove those outputs before
keeping them.

## Finish a task

Preserve useful source and concise findings in Git. Remove task diagnostics
from `out/`; remove obsolete renders by their run directory after preserving
anything worth keeping. Retain working build caches to avoid recompilation.
Do not use a blanket ignored-file deletion: it would also remove useful builds,
images, toolchains, and qualification bundles.

Check `git status --short`, `git branch -vv`, and `git worktree list`. Identify
unfinished work on a development branch, and state whether commits are pushed.
Final integration belongs on `main` once the requested source work is complete.
Local source and focused tests do not establish full scientific, native-device,
or release qualification; those receipts follow [ATTESTATION.md](ATTESTATION.md).
