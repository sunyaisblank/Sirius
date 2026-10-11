"""Reconstruct the preserved native producer whose whole source tree is current."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, importlib.util, json, shutil, subprocess
root = Path.cwd(); work = Path(__file__).resolve().parent
source = json.loads((work / 'source.json').read_text())
producer = source['baseline_revision']; live = source['source_revision']
def git(*args): return subprocess.check_output(['git', *args])
assert git('rev-parse', 'HEAD').decode().strip() == live == git('rev-parse', '@{u}').decode().strip()
assert not git('status', '--porcelain')
assert git('rev-parse', live + '^{tree}') == git('rev-parse', producer + '^{tree}')
bundle = root / 'attestations/native-build' / producer[:7] / 'windows-build'
stage = root / 'bin/windows-msvc' / ('native-' + producer[:7] + '-executed-cost')
receipt = work / 'baseline-reconstruction.json'
assert not stage.exists() and not receipt.exists()
spec = importlib.util.spec_from_file_location('sirius_attestation', root / 'scripts/verify-attestation.py')
v = importlib.util.module_from_spec(spec); spec.loader.exec_module(v)
g = v.load_build_gate_verifier()
data = json.loads((bundle / 'windows-build.json').read_text())
gate = g.validate_native_build_document(json.loads((bundle / 'native_build_gate.json').read_text()))
assert gate['source'] == {'revision': producer, 'clean': True} and data['source_revision'] == producer
model = git('show', producer + ':tests/operating_model.json')
assert model == (root / 'tests/operating_model.json').read_bytes()
v.verify_document_against_authority(data, bundle / 'windows-build.json', producer, hashlib.sha256(model).hexdigest())
stage.mkdir(parents=True)
source_view = stage / 'recorded-source'; records = []
def put(src, relative):
    rel = Path(relative); assert not rel.is_absolute() and '..' not in rel.parts
    dst = stage / rel; dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists(): assert dst.read_bytes() == src.read_bytes()
    else: shutil.copy2(src, dst)
    return dst
def bound(src, record):
    b = src.read_bytes()
    assert len(b) == record['bytes'] and hashlib.sha256(b).hexdigest() == record['sha256']
    if record['root'] == 'build': return put(src, record['path'])
    assert record['root'] == 'source'
    assert b == git('show', producer + ':' + record['path']) == (root / record['path']).read_bytes()
    records.append(record)
    return put(src, Path('recorded-source') / record['path'])
executables = {}
for name, record in gate['tested_artifacts'].items():
    src = bundle / ('qualification-sirius.bin' if name == 'sirius' else v.QUALIFICATION_TEST_EVIDENCE[name])
    executables[name] = bound(src, record)
for name, record in gate['product_artifacts'].items():
    src = bundle / ('qualification-sirius.bin' if name == 'sirius' else ('alignment_receipt.json' if name == 'alignment_receipt' else v.QUALIFICATION_PRODUCT_EVIDENCE[name]))
    bound(src, record)
for name, record in gate['test_input_artifacts'].items(): bound(bundle / v.QUALIFICATION_TEST_INPUT_EVIDENCE[name], record)
for key, name in [('junit', 'native_build_gate_junit.xml'), ('log', 'native-build-gate-log')]: bound(bundle / name, gate['ctest'][key])
put(bundle / 'native_build_gate.json', 'generated/sirius/native_build_gate.json')
for consumer in ['sirius', 'sirius_backend_tests', 'sirius_render_tests']:
    parent = executables[consumer].parent.relative_to(stage)
    for name, rel in v.QUALIFICATION_RUNTIME_RESOURCE_PATHS.items(): put(bundle / v.QUALIFICATION_PRODUCT_EVIDENCE[name], parent / 'resources' / rel)
    put(bundle / 'alignment_receipt.json', parent / 'resources/model/alignment_receipt.json')
    if consumer != 'sirius':
        for name, rel in g.TEST_INPUT_PATHS.items(): put(bundle / v.QUALIFICATION_TEST_INPUT_EVIDENCE[name], parent / 'resources' / rel)
g.verify_recorded_files(gate, source_view, stage)
checked = v.copy_qualification_test_inputs(stage / 'generated/sirius/native_build_gate.json', stage)
assert len(checked) == 20
assert git('rev-parse', 'HEAD').decode().strip() == live and not git('status', '--porcelain')
result = {'native_producer_revision': producer, 'live_source_revision': live, 'whole_source_trees_exact': True, 'stage': str(stage.relative_to(root)), 'official_export': str(bundle.relative_to(root)), 'whole_gate_records_verified': True, 'consumed_inputs_checked': len(checked), 'executables': {n: {'path': str(p.relative_to(root)), 'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for n, p in executables.items()}, 'gate_sha256': hashlib.sha256((bundle / 'native_build_gate.json').read_bytes()).hexdigest(), 'scope': 'Original historical producer reconstructed on byte-identical current source; no new native build, runtime dispatch, performance or scientific qualification.'}
receipt.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: result[k] for k in ('native_producer_revision', 'live_source_revision', 'whole_source_trees_exact', 'consumed_inputs_checked')}))
