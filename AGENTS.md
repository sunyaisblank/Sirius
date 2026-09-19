# Sirius workspace

Keep Sirius development inside this checkout. Do not create sibling clones,
worktrees, toolchains, evidence folders, or test output in `.project`.

Use the existing CMake presets and their ignored `bin/<preset>` directories.
Keep disposable diagnostics under `out/`, rendered output under `renders/`, and
qualification bundles under `attestations/`. Remove temporary artifacts when
their task finishes. Preserve useful source and concise findings in Git before
removing a development workspace.

Before concluding a task, check the working tree, branch tracking, and any
temporary worktrees. State whether changes are committed and pushed. Keep
unfinished development on an identified development branch and distinguish
source/build checks from full scientific and release qualification.
