"""Freeze exact original build evidence, consumers, controllers and provider."""
from pathlib import Path
import hashlib, importlib.util, json, subprocess
from native_contracts import ROOT, WORK, sha

trial = json.loads((WORK / 'candidate-source.json').read_text())
head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert head == trial['source_revision'] and not subprocess.check_output(['git', 'status', '--porcelain'])
output = WORK / 'execution-bindings.json'
assert not output.exists()
spec = importlib.util.spec_from_file_location('attestation', ROOT / 'scripts/verify-attestation.py')
v = importlib.util.module_from_spec(spec); spec.loader.exec_module(v); g = v.load_build_gate_verifier()
files = {}
def bind(path, expected=None):
    path = path.resolve()
    record = {'bytes': path.stat().st_size, 'sha256': sha(path)}
    if expected is not None:
        assert record == {k: expected[k] for k in ('bytes', 'sha256')}, str(path)
    key = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
    assert key not in files or files[key] == record
    files[key] = record

producers = {}
for label, revision in [('baseline', trial['baseline_revision']), ('candidate', head)]:
    stage = ROOT / 'bin/windows-msvc' / ('native-' + revision[:7])
    bundle = ROOT / 'attestations/native-build' / revision[:7] / 'windows-build'
    gate_path = stage / 'generated/sirius/native_build_gate.json'
    assert gate_path.read_bytes() == (bundle / 'native_build_gate.json').read_bytes()
    gate = g.validate_native_build_document(json.loads(gate_path.read_text()))
    assert gate['source'] == {'revision': revision, 'clean': True}
    evidence = bundle / 'windows-build.json'; data = json.loads(evidence.read_text())
    if revision == head:
        v.verify_path(evidence, ROOT)
    else:
        model = subprocess.check_output(['git', 'show', revision + ':tests/operating_model.json'])
        assert model == (ROOT / 'tests/operating_model.json').read_bytes()
        v.verify_document_against_authority(data, evidence, revision, hashlib.sha256(model).hexdigest())
    g.verify_recorded_files(gate, ROOT, stage)
    assert len(v.copy_qualification_test_inputs(gate_path, stage)) == 20 and len(gate['tested_artifacts']) == 7
    pe = json.loads((WORK / (label + '-native-array-bindings.json')).read_text())
    assert pe['pass'] and pe['source_revision'] == revision and len(pe['arrays']) == 40
    for record in pe['products'].values():
        bind(ROOT / record['path'], record)
    for path in stage.rglob('*'):
        assert not path.is_symlink()
        if path.is_file():
            bind(path)
    for group in ('tested_artifacts', 'product_artifacts', 'test_input_artifacts'):
        for record in gate[group].values():
            bind((stage if record['root'] == 'build' else ROOT) / record['path'], record)
    bind(evidence)
    producers[label] = {'source_revision': revision, 'stage': str(stage.relative_to(ROOT)), 'gate_sha256': sha(gate_path), 'canonical_and_consumed_inputs_verified': 20, 'complete_recorded_files_verified': True}
paths = subprocess.check_output(['git', 'ls-files', 'src/sirius/kernels'], text=True).splitlines()
paths += ['src/sirius/backend/vulkan/vulkan_device.cpp', 'src/sirius/backend/retained_compute.cpp',
          'src/sirius/backend/retained_integrator.cpp', 'tests/backend/retained_compute_test.cpp', 'tests/backend/retained_dopri_test.cpp',
          'scripts/verify-attestation.py', 'scripts/verify-build-gate.py', 'scripts/build-retained-kernels.py', '.github/workflows/ci.yml']
for path in paths:
    bind(ROOT / path)
for name in ['run_native_control.ps1','run_native.py','run_native_sequence.py','native_contracts.py',
             'frozen-retention-comparator.py','prepare_native_bindings.py','bind_candidate_native.py',
             'baseline-terminal-audit.ps1','emergency_native_cleanup.ps1','run_terminal_audit.py',
             'candidate-source.json','candidate-build-and-arrays.json','candidate-source-review.json',
             'candidate-software-owner-review.json','independent-software-controls-review.json',
             'candidate-controls-complete.json','candidate-linux-readonly-payloads.json',
             'candidate-native-apparatus-review.json','candidate-powershell-parser.json',
             'baseline-retained_kernels.h','candidate-retained_kernels.h',
             'baseline-native-array-bindings.json','candidate-native-array-bindings.json',
             'reconstruction-windows-eef7baf.json','reconstruction-windows-bd82e52.json',
             'candidate-native-preregistration.txt','candidate-native-preregistration.json',
             'candidate-plan-binding.json','retention-gate-frozen.json',
             'terminal-baseline-baseline-mixed-accepted.json','baseline-native-apparatus-review.json',
             'controls/timestamps/tests.xml','baseline-timestamp-controls/timestamps/result.json']:
    bind(WORK / name)
assert json.loads((WORK/'candidate-controls-complete.json').read_text())['existing_tests']==9
assert len(json.loads((WORK/'candidate-controls-complete.json').read_text())['results'])==10
assert json.loads((WORK/'independent-software-controls-review.json').read_text())['pass']
assert json.loads((WORK/'candidate-native-apparatus-review.json').read_text())['pass']
assert json.loads((WORK/'terminal-baseline-baseline-mixed-accepted.json').read_text())['pass']
assert producers['baseline']['source_revision']=='eef7baf23f2e8c9906a467a6eef9ba7b6317d5ee'
assert producers['candidate']['source_revision']==head
baseline_pe=json.loads((WORK/'baseline-native-array-bindings.json').read_text())
candidate_pe=json.loads((WORK/'candidate-native-array-bindings.json').read_text())
assert not baseline_pe['changed_arrays'] and len(candidate_pe['changed_arrays'])==6 and candidate_pe['unchanged_arrays']==34
assert candidate_pe['changed_arrays']==sorted(json.loads((WORK/'candidate-plan-binding.json').read_text())['changed_expected_arrays'])
provider = {
    '/mnt/c/Windows/System32/vulkan-1.dll': {'bytes': 1719408, 'sha256': '45ee2efc5c6d3986da8181775610f308122c752e55e703340d9733bbbc4b9f1a'},
    '/mnt/c/Windows/System32/DriverStore/FileRepository/amdvlk.inf_amd64_914ba89eaaafdf60/amdvlk64.dll': {'bytes': 103213064, 'sha256': '5c0561968155f4b0aa80fe9a5bbd17072d5942ecb906c687e749dc92058a8650'},
}
for path, record in provider.items():
    bind(Path(path), record)
ps = Path('/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe')
assert sha(ps) == '8bb6fa8c283b4d92120b1ef249a9b311b0f804d4cabbe9981159976c8be76a5e'; bind(ps)
result = {'candidate_revision': head, 'baseline_revision': trial['baseline_revision'], 'producers': producers,
          'files': files, 'file_count': len(files), 'source_tree_clean': True, 'qualification_claimed': False,
          'scope': 'All recorded stages/products/consumed inputs, actual 40-array PE bindings, relevant committed source, apparatus/controller/preregistration and provider bytes. Historical baseline retains its original authority.'}
output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'files': len(files), 'producers': producers, 'dispatches': 0}))
