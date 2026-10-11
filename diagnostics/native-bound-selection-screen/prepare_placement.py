"""Prepare the reviewed existing zero-dispatch driver query for one new module."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
PRIOR = ROOT / 'attestations/native-vulkan/24a8e14/native-rhs-layer-clock-execution/files/placement'
HEAD = '2053ef79c90752d37b0075b8137cfb9b571b6820'

def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}

def read(path):
    return json.loads(path.read_text())

assert WORK == ROOT / 'out/native-bound-selection'
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == HEAD
assert not subprocess.check_output(['git', 'status', '--porcelain'], text=True)
admission = read(WORK / 'compile-admission.json')
assert admission['accepted'] and admission['source_revision'] == HEAD
assert seal(WORK / 'phase/retained_transport.spv') == admission['spv']
assert seal(WORK / 'phase/retained_transport.spvasm') == admission['assembly']
assert seal(WORK / 'compiler.json') == admission['compiler_receipt']
compiled_review = read(WORK / 'compiled-review.json')
assert compiled_review['accepted']
assert seal(WORK / 'compile-admission.json') == compiled_review['admission_seal']
review = read(PRIOR / 'source-review.json')
assert review['accepted']
target = WORK / 'placement'
assert not target.exists()
target.mkdir()
records = {}
names = ('native_clock_placement.cs', 'placement-owner.ps1', 'placement-child.ps1',
         'run_placement_query.py', 'emergency_native_cleanup.ps1',
         'run_native_retirement.py', 'verify_native_retirement.ps1')
for name in names:
    source = PRIOR / name
    if name in review['controller_seals']:
        assert seal(source) == review['controller_seals'][name]
    text = source.read_text()
    if name == 'native_clock_placement.cs':
        old = 'e3a1589d852a4254312bb0e9e4b0343ab6aea32855c02c4f0f93b9fe1e18e5a3'
        assert text.count(old) == 1
        text = text.replace(old, seal(WORK / 'phase/retained_transport.spv')['sha256'])
    if name == 'run_placement_query.py':
        old = "WORK.parent/'prepare_layer_placement.py'"
        assert text.count(old) == 1
        text = text.replace(old, "WORK.parent/'prepare_placement.py'")
    (target / name).write_text(text)
    records[name] = {'basis': str(source.relative_to(ROOT)), 'basis_seal': seal(source),
                     'current_seal': seal(target / name),
                     'byte_identical': source.read_bytes() == (target / name).read_bytes()}
for name in ('retained_transport.spv', 'retained_endpoint.spv'):
    source = WORK / 'phase' / name if name == 'retained_transport.spv' else PRIOR / name
    (target / name).write_bytes(source.read_bytes())
providers = read(PRIOR / 'plan.json')['providers']
paths = [p for p in target.iterdir() if p.is_file()]
paths += [Path(__file__), WORK / 'compiler.json', WORK / 'compile-admission.json', WORK / 'compiled-review.json',
          WORK / 'source-review.json', WORK / 'prototype.patch',
          WORK / 'kernels/retained_pair.slang', WORK / 'phase/retained_transport.slang']
inputs = {str(p.relative_to(ROOT)): seal(p) for p in paths}
for provider in providers:
    assert seal(Path(provider['linux_path'])) == {k: provider[k] for k in ('bytes', 'sha256')}
plan = {'source_revision': HEAD, 'selected_stage': 'Transport',
        'input_seals': inputs, 'providers': providers,
        'purpose': 'One isolated native-default Transport integer-mask bound candidate '
        'with unchanged diagnostic observer. One ordinary pipeline/static shader-info '
        'query, zero dispatch. Original 50s/512MiB child guard,45s intrinsic watchdog, '
        '120s bridge,8MiB payload limit unchanged. No numerical, benefit, rendering '
        'or qualification verdict.'}
(target / 'plan.json').write_text(json.dumps(plan, indent=2) + '\n')
(target / 'preparation.json').write_text(json.dumps(
    {'sources': records, 'plan_seal': seal(target / 'plan.json')}, indent=2) + '\n')
print(json.dumps({'plan': seal(target / 'plan.json'), 'sources': records}))
