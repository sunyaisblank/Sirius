"""Prepare exact-source restoration helpers; never build or change production source."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, subprocess

ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
BASELINE = '2053ef79c90752d37b0075b8137cfb9b571b6820'
assert WORK == ROOT/'out/native-division-first-correction-screen'
OLD = ROOT/'attestations/native-vulkan/451e2cb/native-division-phase-observer/files'

def read(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))

def seal(p):
    return {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}

def put(p, value):
    assert not p.exists(), p
    p.write_text(value)

tracked = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', BASELINE], text=True).splitlines()
preparation = read(WORK/'runtime/preparation.json')
assert preparation['accepted'] and preparation['baseline_revision'] == BASELINE
originals = {name: preparation['input_seals'][name] for name in tracked}
assert len(originals) == 532
changed = {'scripts/build-retained-kernels.py', 'src/sirius/kernels/retained_triple.slang'}
for name, expected in originals.items():
    p = WORK/'runtime/baseline'/name if name in changed else ROOT/name
    assert seal(p) == expected, name
    if name in changed:
        blob = subprocess.check_output(['git', 'show', BASELINE+':'+name])
        assert expected == {'bytes': len(blob), 'sha256': hashlib.sha256(blob).hexdigest()}
header = WORK/'runtime/baseline/bin/linux-gcc/src/sirius/backend/retained/retained_kernels.h'
arrays = read(WORK/'runtime/baseline-arrays.json')
assert seal(header) == {'bytes': arrays['header_bytes'], 'sha256': arrays['header_sha256']}
put(WORK/'baseline-source.json', json.dumps({'revision': BASELINE,
    'tracked_source_seals': originals, 'production_header': seal(header),
    'scope': 'Original accepted source epoch: two source versions retained before ea6 application; other 530 files still exact. Current candidate ROOT is intentionally different.'}, indent=2)+'\n')

outputs = {}
for name in ('restore_build.py', 'verify_restored_products.py'):
    original = (OLD/name).read_text()
    target = original.replace("ROOT/'out/native-division-phase-observer'", "ROOT/'out/native-division-first-correction-screen'")
    assert original != target
    if name == 'verify_restored_products.py':
        old = "transport = WORK/'phase/retained_transport.spv'\nendpoint = ROOT/'attestations/native-vulkan/2053ef7/native-shader-clock-admission/files/phase/retained_endpoint.spv'"
        new = "transport = WORK/'runtime/candidate/retained_transport.spv'\nendpoint = WORK/'runtime/candidate/retained_endpoint.spv'"
        assert target.count(old) == 1
        target = target.replace(old, new).replace('diagnostic_modules_absent', 'rejected_candidate_modules_absent').replace('diagnostic modules absent', 'rejected candidate modules absent')
    put(WORK/name, target)
    outputs[name] = seal(WORK/name)
inputs = {str(p.relative_to(ROOT)): seal(p) for p in (Path(__file__), OLD/'restore_build.py', OLD/'verify_restored_products.py', WORK/'runtime/preparation.json', WORK/'runtime/baseline-arrays.json', header)}
put(WORK/'restoration-preparation.json', json.dumps({'accepted': True, 'input_seals': inputs,
    'output_seals': outputs, 'baseline_source_seal': seal(WORK/'baseline-source.json'),
    'production_changes_or_builds': False, 'controller_governors_unchanged': True,
    'reader_checks': 'All532 accepted source seals, exact original header and all40 readonly ELF payloads in app/backend; whole candidate T/E payloads absent.'}, indent=2)+'\n')
print(json.dumps({'accepted': True, 'outputs': outputs, 'baseline_source_seal': seal(WORK/'baseline-source.json')}))
