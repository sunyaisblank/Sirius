"""Bind one ordinary zero-dispatch driver query to the exact division observer."""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import hashlib, json, subprocess
ROOT = Path.cwd()
WORK = Path(__file__).resolve().parent
OLD = ROOT / 'attestations/native-vulkan/24a8e14/native-rhs-layer-clock-execution/files/placement'
BASELINE = '2053ef79c90752d37b0075b8137cfb9b571b6820'
def read(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))
def seal(p):
    return {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
def git(*a):
    return subprocess.check_output(['git', *a], text=True, timeout=30).strip()
def one(text, old, new):
    assert text.count(old) == 1, (old, text.count(old))
    return text.replace(old, new)
assert ROOT.resolve() == Path('/home/astra/.project/Sirius')
assert WORK == ROOT / 'out/native-division-phase-observer'
assert git('rev-parse', 'HEAD') == BASELINE == git('rev-parse', '@{u}') and not git('status', '--porcelain')
compiler = read(WORK / 'compiler.json')
assert compiler['accepted'] and compiler['baseline_revision'] == BASELINE
compiled = read(WORK / 'compiled-review.json')
assert compiled['accepted'] and compiled['baseline_revision'] == BASELINE
required = {WORK / 'compiler.json', WORK / 'phase/retained_transport.spv',
            WORK / 'phase/retained_transport.spvasm', WORK / 'phase/retained_transport.compiler.spv',
            WORK / 'phase/retained_transport.slang', WORK / 'preparation.json',
            WORK / 'source-review.json', WORK / 'compiler-controller-review.json'}
required |= {Path(name) for name in compiler['input_seals']}
assert all(str(p) in compiled['reviewed_seals'] for p in required)
assert seal(WORK / 'phase/retained_transport.spv') == compiler['product']['spv']
assert seal(WORK / 'phase/retained_transport.spvasm') == compiler['product']['assembly']
for name, expected in compiled['reviewed_seals'].items():
    assert seal(Path(name)) == expected
prior_review = read(OLD / 'source-review.json')
assert prior_review['accepted']
target = WORK / 'placement'
assert not target.exists()
target.mkdir()
records = {}
for name in ('native_clock_placement.cs', 'placement-owner.ps1', 'placement-child.ps1',
             'run_placement_query.py', 'emergency_native_cleanup.ps1',
             'run_native_retirement.py', 'verify_native_retirement.ps1'):
    source = OLD / name
    assert seal(source) == prior_review['controller_seals'][name]
    text = source.read_text()
    if name == 'native_clock_placement.cs':
        text = one(text, 'e3a1589d852a4254312bb0e9e4b0343ab6aea32855c02c4f0f93b9fe1e18e5a3',
                   seal(WORK / 'phase/retained_transport.spv')['sha256'])
        text = one(text, 'One fixed Transport unadopted native-default phase-observer module',
                   'One exact unadopted native-default division primitive observer module')
    if name == 'run_placement_query.py':
        text = one(text, "WORK.parent/'prepare_layer_placement.py'", "WORK.parent/'prepare_division_placement.py'")
        text = one(text, "review=read(WORK/'source-review.json');assert review['accepted']",
                   "review=read(WORK/'source-review.json');assert review['accepted'] and review['baseline_revision']==head\nassert ROOT.resolve()==Path('/home/astra/.project/Sirius') and WORK==ROOT/'out/native-division-phase-observer/placement'\nfor name,expected in review['reviewed_seals'].items():assert seal(Path(name))==expected")
        text = one(text, "seals={seal_key(p):seal(p) for p in paths}",
                   "paths=sorted(set(paths)|{Path(name) for name in review['reviewed_seals']})\nseals={seal_key(p):seal(p) for p in paths}")
    (target / name).write_text(text)
    records[name] = {'basis': str(source.relative_to(ROOT)), 'basis_seal': seal(source),
                     'current_seal': seal(target / name), 'byte_identical': source.read_bytes() == (target / name).read_bytes()}
for name in ('retained_transport.spv', 'retained_endpoint.spv'):
    source = (WORK / 'phase' / name if name == 'retained_transport.spv' else
              ROOT / 'attestations/native-vulkan/2053ef7/native-shader-clock-admission/files/phase' / name)
    (target / name).write_bytes(source.read_bytes())
providers = read(OLD / 'plan.json')['providers']
inputs = {str(p.relative_to(ROOT)): seal(p) for p in target.iterdir() if p.is_file()}
for p in (Path(__file__), WORK / 'phase/retained_transport.slang', WORK / 'compiled-review.json',
          WORK / 'compiler.json', WORK / 'preparation.json', WORK / 'source-review.json'):
    inputs[str(p.relative_to(ROOT))] = seal(p)
for provider in providers:
    assert seal(Path(provider['linux_path'])) == {k: provider[k] for k in ('bytes', 'sha256')}
plan = {'source_revision': BASELINE, 'selected_stage': 'Transport', 'input_seals': inputs,
        'providers': providers,
        'purpose': 'One exact new native-default Transport division primitive observer module, ordinary pipeline/static shader-info query, zero dispatch. Existing50s/512MiB child guard,45s intrinsic watchdog,120s bridge,8MiB payload limit unchanged. Actual native code placement prerequisite; no exclusive primitive cost, clock rate, numerical equivalence, rendering benefit or qualification claim.'}
(target / 'plan.json').write_text(json.dumps(plan, indent=2) + '\n')
(target / 'preparation.json').write_text(json.dumps({'sources': records, 'plan_seal': seal(target / 'plan.json'),
    'scope': 'Source preparation only; no helper/compiler/native invocation.'}, indent=2) + '\n')
print(json.dumps({'sources': {n: r['current_seal'] for n, r in records.items()}, 'plan': seal(target / 'plan.json')}))
