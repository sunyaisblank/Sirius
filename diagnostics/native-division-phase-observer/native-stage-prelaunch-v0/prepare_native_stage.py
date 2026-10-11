"""Validate the existing CI export and reconstruct its exact ignored native stage."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, importlib.util, json, shutil, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
HEAD = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert HEAD == '451e2cbf833b87f51b1819421d68163a4b5fd1e5'
assert HEAD == subprocess.check_output(['git', 'rev-parse', '@{u}'], text=True).strip()
assert not subprocess.check_output(['git', 'status', '--porcelain'], text=True)
ARTIFACT = WORK / 'artifact'
STAGE = ROOT / 'bin/windows-msvc/native-451e2cb-division-phase-observer'
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
    path = (WORK / 'phase/retained_transport.spv' if name == 'retained_transport.spv' else ROOT / 'attestations/native-vulkan/2053ef7/native-shader-clock-admission/files/phase/retained_endpoint.spv')
    data = path.read_bytes()
    assert data in exe, 'Exact observer module missing from exported executable'
    embedded[name] = {**seal(path), 'exact_bytes_present_in_executable': True}

providers = json.loads((ROOT / 'attestations/native-vulkan/2053ef7/native-shader-clock-placement/files/transport/plan.json').read_text())['providers']
for provider in providers:
    assert seal(Path(provider['linux_path'])) == {k: provider[k] for k in ('bytes', 'sha256')}
inputs = {str(p.relative_to(ROOT)): seal(p) for p in STAGE.rglob('*') if p.is_file()}
for path in ('src/sirius/backend/vulkan/vulkan_device.cpp', 'tests/backend/retained_compute_test.cpp', 'tests/backend/shader_clock_observer.h', 'out/native-division-phase-observer/host-source-review.json', 'out/native-division-phase-observer/prepare_native_stage.py', 'out/native-division-phase-observer/phase/retained_transport.spv', 'out/native-division-phase-observer/compiled-review.json', 'out/native-division-phase-observer/compiler.json', 'out/native-division-phase-observer/preparation.json'):
    inputs[path] = seal(ROOT / path)
plan = {
    'source_revision': HEAD, 'build_root': str(STAGE.relative_to(ROOT)),
    'executable_directory': str(tested['sirius_backend_tests'].parent.relative_to(ROOT)),
    'gate_sha256': seal(gate_path)['sha256'], 'input_seals': inputs,
    'providers': providers, 'exact_observer_modules': embedded,
    'cases': {'device-timestamps': 'RetainedComputeTest.DeviceTimestampsPreserveOriginalIntervalResults',
              'mixed-transport': 'RetainedComputeTest.MixedIndependentTransportLayersPreserveWordsAndRefusal',
              'factored-endpoints': 'RetainedComputeTest.FactoredEndpointsPreserveIndependentRootsAndBoundaryAdmission'},
    'purpose': 'Original numerical/ordered-stream/refusal/recovery controls with uncalibrated native-default shader-clock phase/division sidecars and exact input snapshots; original inputs/modes/assertions/8MiB/governors intact. Existing 320s/4GiB native tool guard, 400s bridge. No full science/release or ordinary-benefit claim.'}
(WORK / 'execution-plan.json').write_text(json.dumps(plan, indent=2) + '\n')
print(json.dumps({'source_revision': HEAD, 'stage_files': len(inputs) - 9, 'native_backend': seal(tested['sirius_backend_tests']), 'observer_modules': embedded}))
