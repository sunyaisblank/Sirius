"""Validate the existing CI export and reconstruct its exact ignored native stage."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, importlib.util, json, shutil, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
HEAD = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
SOURCE = json.loads((WORK / 'source.json').read_text())
assert HEAD == SOURCE['source_revision']
assert HEAD == subprocess.check_output(['git', 'rev-parse', '@{u}'], text=True).strip()
assert not subprocess.check_output(['git', 'status', '--porcelain'], text=True)
ARTIFACT = WORK / 'artifact'
stage_relative = Path(SOURCE['stage_root'])
expected_stage = Path('bin/windows-msvc') / ('native-' + HEAD[:7] + '-stage-seven-word-witness')
assert stage_relative == expected_stage and not stage_relative.is_absolute()
assert '..' not in stage_relative.parts and ROOT == ROOT.resolve()
STAGE = ROOT / stage_relative
for ancestor in (ROOT / 'bin', ROOT / 'bin/windows-msvc', STAGE):
    assert not ancestor.is_symlink(), 'Native stage ancestry must not redirect writes'
assert not STAGE.exists()

def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

def safe(root, relative):
    path = root / relative
    assert not Path(relative).is_absolute() and '..' not in Path(relative).parts
    assert path.resolve().is_relative_to(root.resolve())
    return path

document = json.loads((ARTIFACT / 'windows-build.json').read_text())
assert document['status'] == 'pass' and document['source_revision'] == HEAD
assert document['domains'] == ['windows-native-build']
assert document['platform'] == {'os': 'windows', 'wsl2': False}
assert document['claims']['compiled_topology_complete']
assert document['claims']['runtime_ready'] is False
payloads = {}
for name, record in document['artifacts'].items():
    path = safe(ARTIFACT, record['path'])
    assert seal(path) == {k: record[k] for k in ('bytes', 'sha256')}
    payloads.setdefault((record['bytes'], record['sha256']), path)

sys.path.insert(0, str(ROOT / 'scripts'))
spec = importlib.util.spec_from_file_location('native_gate_verifier', ROOT / 'scripts/verify-build-gate.py')
verifier = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verifier)
gate = json.loads((ARTIFACT / 'native_build_gate.json').read_text())
verifier.validate_native_build_document(gate)
assert gate['source'] == {'revision': HEAD, 'clean': True}
STAGE.mkdir(parents=True)

def copy_record(record, destination):
    source = payloads[(record['bytes'], record['sha256'])]
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        assert destination.read_bytes() == source.read_bytes()
    else:
        shutil.copyfile(source, destination)
    assert seal(destination) == {k: record[k] for k in ('bytes', 'sha256')}

tested = {}
for name, record in gate['tested_artifacts'].items():
    assert record['root'] == 'build'
    tested[name] = safe(STAGE, record['path'])
    copy_record(record, tested[name])
for record in list(gate['product_artifacts'].values()) + list(gate['test_input_artifacts'].values()):
    if record['root'] == 'build':
        copy_record(record, safe(STAGE, record['path']))
    else:
        assert record['root'] == 'source'
        assert seal(safe(ROOT, record['path'])) == {k: record[k] for k in ('bytes', 'sha256')}

for name, record in gate['product_artifacts'].items():
    if name == 'sirius': continue
    relative = verifier.INSTALLED_PRODUCTS[name].removeprefix('share/sirius/')
    copy_record(record, safe(tested['sirius'].parent / 'resources', relative))
for consumer in ('sirius_backend_tests', 'sirius_render_tests'):
    resource = tested[consumer].parent / 'resources'
    for name, record in gate['test_input_artifacts'].items():
        copy_record(record, safe(resource, verifier.TEST_INPUT_PATHS[name]))
    for name in ('trace_spv', 'trace_fp32comp_spv', 'trace_fp64_spv'):
        copy_record(gate['product_artifacts'][name], safe(resource, verifier.INSTALLED_PRODUCTS[name].removeprefix('share/sirius/')))
verifier.verify_test_input_copies(tested, gate['product_artifacts'], gate['test_input_artifacts'])
gate_path = STAGE / 'generated/sirius/native_build_gate.json'
gate_path.parent.mkdir(parents=True, exist_ok=True)
shutil.copyfile(ARTIFACT / 'native_build_gate.json', gate_path)
exe = tested['sirius_backend_tests'].read_bytes()
embedded = {}
for name in ('retained_transport.spv', 'retained_endpoint.spv'):
    path = WORK / 'kernels' / name
    stage = 'transport' if name == 'retained_transport.spv' else 'endpoint'
    reviewed = json.loads((WORK / 'compiled-review.json').read_text())
    assert reviewed['accepted'] and seal(path) == reviewed['modules'][stage]
    data = path.read_bytes()
    assert data in exe, 'Exact observer module missing from exported executable'
    embedded[name] = {**seal(path), 'exact_bytes_present_in_executable': True}

providers = json.loads((ROOT / 'attestations/native-vulkan/2053ef7/native-shader-clock-placement/files/transport/plan.json').read_text())['providers']
for provider in providers:
    assert seal(Path(provider['linux_path'])) == {k: provider[k] for k in ('bytes', 'sha256')}
inputs = {str(p.relative_to(ROOT)): seal(p) for p in STAGE.rglob('*') if p.is_file()}
tracked = subprocess.check_output(['git', 'ls-files', '-z']).split(b'\0')
for raw in tracked:
    if raw:
        path = raw.decode(); inputs[path] = seal(ROOT / path)
for path in list(WORK.iterdir()) + list((WORK / 'kernels').iterdir()):
    if path.is_file() and path.name not in ('execution-plan.json', 'controller-review.json'):
        inputs[str(path.relative_to(ROOT))] = seal(path)
host = json.loads((WORK / 'host-source-review.json').read_text())
assert host['accepted'] and host['source_revision'] == SOURCE['baseline_revision']
for name in ('word_value_observer.h', 'retained_compute_test.cpp'):
    assert seal(ROOT / 'tests/backend' / name) == host['outputs'][name]
for name, digest in host['inputs'].items():
    assert seal(WORK / name)['sha256'] == digest
compiled = json.loads((WORK / 'compiled-review.json').read_text())
assert compiled['accepted'] and compiled['source_revision'] == SOURCE['baseline_revision']
assert compiled['schema_seal'] == seal(WORK / 'value-schema.json')
run = json.loads((WORK / 'ci-run.json').read_text())
ci = json.loads((WORK / 'ci.json').read_text())
assert run['source_revision'] == ci['source_revision'] == HEAD
assert ci['accepted'] and ci['run_id'] == run['run_id']
assert ci['metadata']['id'] == run['run_id'] and ci['metadata']['head_sha'] == HEAD
assert ci['metadata']['head_branch'] == SOURCE['branch']
assert ci['metadata']['event'] == 'workflow_dispatch'
assert ci['metadata']['status'] == 'completed' and ci['metadata']['conclusion'] == 'success'
success = {'integration-no-render', 'integration-windows-no-render', 'integration-macos-no-render'}
skipped = {'linux-gate', 'linux-sanitizers', 'windows-build', 'macos-build'}
jobs = ci['jobs']['jobs']
assert len(jobs) == 7 and {j['name'] for j in jobs} == success | skipped
assert all(j['status'] == 'completed' and j['conclusion'] == ('success' if j['name'] in success else 'skipped') for j in jobs)
extraction = json.loads((WORK / 'artifact-extraction.json').read_text())
assert extraction['accepted'] and extraction['source_revision'] == HEAD and extraction['run_id'] == run['run_id']
artifact = extraction['artifact']
assert artifact['name'] == 'sirius-windows-native-attestation'
assert artifact['workflow_run']['id'] == run['run_id'] and artifact['workflow_run']['head_sha'] == HEAD
assert not artifact['expired']
assert seal(WORK / 'artifact.zip') == {'bytes': artifact['size_in_bytes'], 'sha256': artifact['digest'].removeprefix('sha256:')}
api_artifacts = json.loads((WORK / 'ci-artifacts.json').read_text())['artifacts']
assert [a for a in api_artifacts if a['id'] == artifact['id']] == [artifact]
assert len(extraction['extracted_files']) == 40
for record in extraction['extracted_files']:
    assert seal(safe(ARTIFACT, record['path'])) == {k: record[k] for k in ('bytes', 'sha256')}
plan = {
    'source_revision': HEAD, 'build_root': str(STAGE.relative_to(ROOT)),
    'executable_directory': str(tested['sirius_backend_tests'].parent.relative_to(ROOT)),
    'gate_sha256': seal(gate_path)['sha256'], 'input_seals': inputs,
    'providers': providers, 'exact_observer_modules': embedded,
    'cases': {'device-timestamps': 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults',
              'mixed-transport': 'RetainedComputeTest.MixedIndependentTransportLayersPreserveWordsAndRefusal',
              'factored-endpoints': 'RetainedComputeTest.FactoredEndpointsPreserveIndependentRootsAndBoundaryAdmission'},
    'purpose': 'Original numerical/ordered-stream/refusal/recovery controls with native-default word-value sidecars and complete actual input/output snapshots; original inputs/modes/assertions/8MiB/governors intact. Existing 320s/4GiB native tool guard, 400s bridge. No full science/release or ordinary-benefit claim.'}
(WORK / 'execution-plan.json').write_text(json.dumps(plan, indent=2) + '\n')
print(json.dumps({'source_revision': HEAD, 'stage_files': sum(p.is_file() for p in STAGE.rglob('*')), 'native_backend': seal(tested['sirius_backend_tests']), 'observer_modules': embedded}))
