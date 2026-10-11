"""Bind a single zero-dispatch Transport query to the current exact layer module."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, subprocess
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
OLD = ROOT / 'attestations/native-vulkan/2053ef7/native-shader-clock-placement/files/transport'
HEAD = '24a8e145d1da5dca75a4fa26ec8e639020eecf37'
def seal(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip() == HEAD
assert subprocess.check_output(['git', 'rev-parse', '@{u}'], text=True).strip() == HEAD
assert not subprocess.check_output(['git', 'status', '--porcelain'], text=True)
prior_review = json.loads((OLD / 'source-review.json').read_text())
assert prior_review['accepted']
target = WORK / 'placement'
assert not target.exists()
target.mkdir()
records = {}
for name in ('native_clock_placement.cs', 'placement-owner.ps1', 'placement-child.ps1',
             'run_placement_query.py', 'emergency_native_cleanup.ps1'):
    source = OLD / name
    assert seal(source) == prior_review['controller_seals'][name]
    text = source.read_text()
    if name == 'native_clock_placement.cs':
        old_hash = '0eff0f81bf94582c5401e544b28bca9056ee77beb0264d4ece541aa8443dd076'
        assert text.count(old_hash) == 1
        text = text.replace(old_hash, seal(WORK / 'phase/retained_transport.spv')['sha256'])
    if name == 'run_placement_query.py':
        assert text.count("WORK/'prepare_single_stage.py'") == 1
        text = text.replace("WORK/'prepare_single_stage.py'", "WORK.parent/'prepare_layer_placement.py'")
        old = "record['accepted']=record['accepted'] and record['linux_child_group_absent']\nprefix.with_suffix"
        assert text.count(old) == 1
        text = text.replace(old, "record['accepted']=record['accepted'] and record['linux_child_group_absent']\nif signals:record['outer_stop']=True;record['accepted']=False\nrecord['controller_signals']=list(signals)\nprefix.with_suffix")
    (target / name).write_text(text)
    records[name] = {'basis': str(source.relative_to(ROOT)), 'basis_seal': seal(source),
                     'current_seal': seal(target / name), 'byte_identical': source.read_bytes() == (target / name).read_bytes()}
for name in ('run_native_retirement.py', 'verify_native_retirement.ps1'):
    source = WORK / name
    assert seal(source) == json.loads((WORK / 'native-controller-source-review.json').read_text())['controller_seals'][name]
    text = source.read_text()
    if name == 'run_native_retirement.py':
        old = "assert len(sys.argv)==2 and sys.argv[1] in ('device-timestamps','mixed-transport','factored-endpoints')\nWORK=HERE/sys.argv[1];action='retirement';signals=[]"
        assert text.count(old) == 1
        text = text.replace(old, "assert len(sys.argv)==1\nWORK=HERE;action='retirement';signals=[]")
    (target / name).write_text(text)
    records[name] = {'basis': str(source.relative_to(ROOT)), 'basis_seal': seal(source), 'current_seal': seal(target / name)}
for name in ('retained_transport.spv', 'retained_endpoint.spv'):
    source = WORK / 'phase' / name if name == 'retained_transport.spv' else ROOT / 'attestations/native-vulkan/2053ef7/native-shader-clock-admission/files/phase' / name
    (target / name).write_bytes(source.read_bytes())
providers = json.loads((OLD / 'plan.json').read_text())['providers']
inputs = {str(p.relative_to(ROOT)): seal(p) for p in target.iterdir() if p.is_file()}
for p in (Path(__file__), WORK / 'phase/retained_transport.slang', WORK / 'compiled-review.json', WORK / 'compiler.json', ROOT / 'tests/backend/shader_clock_observer.h'):
    inputs[str(p.relative_to(ROOT))] = seal(p)
for provider in providers:
    assert seal(Path(provider['linux_path'])) == {k:provider[k] for k in ('bytes', 'sha256')}
plan = {'source_revision': HEAD, 'selected_stage': 'Transport', 'input_seals': inputs, 'providers': providers,
        'purpose': 'One exact new Transport layer-observer module, ordinary pipeline/static shader-info query, zero dispatch. Existing50s/512MiB child guard,45s intrinsic watchdog,120s bridge,8MiB payload limit unchanged. No clock rate, exclusive layer cost, numerical or rendering benefit claim.'}
(target / 'plan.json').write_text(json.dumps(plan, indent=2) + '\n')
(target / 'preparation.json').write_text(json.dumps({'sources': records, 'plan_seal': seal(target / 'plan.json'), 'scope': 'Source preparation only; no helper/compiler/native invocation.'}, indent=2) + '\n')
print(json.dumps({'sources': {n:r['current_seal'] for n,r in records.items()}, 'plan': seal(target / 'plan.json')}))
