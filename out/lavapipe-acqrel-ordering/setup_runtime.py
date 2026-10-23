"""Reuse original controls/readers and exact matched ELF for a provider-only trial."""
from pathlib import Path
import hashlib, json, shutil

R = Path.cwd(); W = Path(__file__).resolve().parent
OLD = R / 'out/native-original-input-ownership'
assert json.loads((W / 'ordinary/ordinary-code-review.json').read_text())['structural_premise_observed']
assert not (W / 'runtime').exists() and not (W / 'readers').exists()
def pin(p): return {'bytes': p.stat().st_size,
    'sha256': hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()}
providers = json.loads((OLD / 'providers.json').read_text())
identity = json.loads((W / 'provider-identity.json').read_text())
candidate = providers['providers']['candidate']
candidate.update(library=str(R / identity['library']), ICD=str(R / identity['ICD']),
    library_identity=identity['identity'], ICD_identity=identity['ICD_identity'])
candidate['historical_identity_only'] = 'historical_library binds the unchanged expected hardware identity from the original provider; new actual library/features/properties/code/boundaries are independently bound below, and observed C8/ABBA identities must match.'
(W / 'providers.json').write_text(json.dumps(providers, indent=2) + '\n')
READ = W / 'readers'; RUN = W / 'runtime'; READ.mkdir(); RUN.mkdir()
for name in ['native_evidence.py', 'numeric_gates.py']:
    (READ / name).symlink_to('../../native-original-input-ownership/readers/' + name)
shutil.copyfile(OLD / 'readers/expected-readback-layout.json', READ / 'expected-readback-layout.json')
wrappers = {
    'read_results.py': [
        ("RUN=R/'out/native-original-input-ownership/runtime'", "RUN=R/'out/lavapipe-acqrel-ordering/runtime'"),
        ("doc(W.parent/'compiler-result.json')", "doc(R/'out/native-original-input-ownership/compiler-result.json')"),
        ("doc(W.parent/'previous-control-8-acceptance.json')", "doc(R/'out/native-original-input-ownership/readers/control-8-acceptance.json')")],
    'read_controls.py': [
        ("doc(W.parent / 'compiler-result.json')", "doc(R / 'out/native-original-input-ownership/compiler-result.json')"),
        ("doc(W.parent / 'provider-history/candidate-default/device-features.json')", "doc(W.parent / 'ordinary/transport/device-features.json')"),
        ("doc(W.parent / 'provider-history/candidate-default/float-controls-properties.json')", "doc(W.parent / 'ordinary/transport/float-controls-properties.json')")],
    'run_abba.py': [
        ("RUN=R/'out/native-original-input-ownership/runtime'", "RUN=R/'out/lavapipe-acqrel-ordering/runtime'"),
        ("RUN/'run_case.py'", "RUN.parent/'run_case.py'")],
    'decimal_check.py': [("RUN=R/'out/native-original-input-ownership/runtime'", "RUN=R/'out/lavapipe-acqrel-ordering/runtime'")]
}
for name, replacements in wrappers.items():
    basis = OLD / 'readers' / name
    source = basis.read_text()
    for before, after in replacements:
        assert source.count(before) == 1, (name, before)
        source = source.replace(before, after)
    if name == 'read_results.py':
        prefix = "from pathlib import Path\nimport importlib.util\n_binding = Path(__file__).resolve().parents[1] / 'provider_binding.py'\n_spec = importlib.util.spec_from_file_location('actual_provider_binding', _binding)\n_module = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_module)\n_module.admit()\n"
    else: prefix = ''
    # Preserve all original predicates; alter only artifact roots and new provider binding.
    wrapper = prefix + 'from pathlib import Path\nimport hashlib\n'
    wrapper += '_basis = Path.cwd() / ' + repr(str(basis.relative_to(R))) + '\n'
    wrapper += 'assert hashlib.sha256(_basis.read_bytes()).hexdigest() == ' + repr(pin(basis)['sha256']) + '\n'
    wrapper += '_source = _basis.read_text()\n'
    for before, after in replacements:
        wrapper += 'assert _source.count(' + repr(before) + ') == 1\n'
        wrapper += '_source = _source.replace(' + repr(before) + ', ' + repr(after) + ')\n'
    wrapper += 'exec(compile(_source, str(_basis), "exec"))\n'
    (READ / name).write_text(wrapper)
for name in ['run_owned.py', 'owned_guard.py']:
    if name == 'owned_guard.py':
        (RUN / name).symlink_to('../../native-cpu-row-ownership/cold-frame/owned_guard.py')
    else:
        (RUN / name).write_text('"""Reuse the preserved owner with this action directory."""\nfrom pathlib import Path\n_preserved_runner = Path(__file__).resolve().parents[2] / "native-cpu-row-ownership/cold-frame/run_owned.py"\nexec(compile(_preserved_runner.read_text(), str(_preserved_runner), "exec"))\n')
print(json.dumps({'original_predicates_reused': True, 'same_existing_ELF_both_arms': True,
    'new_shader_or_CPP_compilations': 0, 'new_workload_executions': 0}))
