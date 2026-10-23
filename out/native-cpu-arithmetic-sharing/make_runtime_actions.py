"""Bind original cases, prerequisites and limits to this exact private product."""
from pathlib import Path
import json
import re
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
def doc(p): return json.loads(p.read_text())
assert doc(W / 'ordinary/ordinary-code-review.json')['structural_premise_observed']
matched = doc(W / 'matched-result.json')
assert matched['passed']
old = doc(R / 'out/native-original-input-ownership/runtime/actions.json')
providers = doc(W / 'providers.json')['providers']
inputs = {str(p.relative_to(R)) for p in W.rglob('*') if p.is_file()}
inputs |= {matched['executable']}
inputs |= {x['module'] for x in doc(W / 'compiler-result.json')['modules'].values()}
for folder in [R / 'bin/linux-gcc/src/sirius/backend/retained',
               R / 'out/native-original-input-ownership/provider-history',
               R / 'out/native-original-input-ownership/reused/prior-proof']:
    inputs |= {str(p.relative_to(R)) for p in folder.rglob('*') if p.is_file()}
inputs |= {'out/native-original-input-ownership/readers/control-8-acceptance.json',
           'out/native-cpu-row-ownership/cold-frame/run_owned.py',
           'out/native-cpu-row-ownership/cold-frame/owned_guard.py',
           '/etc/ld.so.cache', '/usr/bin/python3'}
for name in doc(W / 'provider/review-admission.json')['source_files']:
    inputs.add('out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8/' + name)
for p in [R / matched['executable'], Path('/usr/bin/python3'),
          *[Path(v['library']) for v in providers.values()]]:
    inputs |= {str(p), str(p.resolve(strict=True))}
    text = subprocess.check_output(['/usr/bin/ldd', str(p)], text=True)
    inputs |= {str(Path(x).resolve()) for x in re.findall(r'(/[^\s()]+)', text)
               if Path(x).is_file()}
for v in providers.values():
    inputs |= {v['library'], v['ICD']}
assert all((Path(p) if p.startswith('/') else R / p).is_file() for p in inputs)
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
actions = {}
for label, basis in old.items():
    action = {k: basis[k] for k in ['arm', 'case', 'prerequisites', 'kind',
                                   'expected_exit', 'limits']}
    provider = providers[action['arm']]
    action.update(
        revision=revision,
        argv=['/usr/bin/python3', '-B', str(W / 'run_case.py'), label, '{folder}'],
        frozen_inputs=sorted(inputs), software_provider=True,
        runtime_icd=provider['ICD'], runtime_provider=provider['library'],
        scope='Original ' + label + ': exact matched host, unchanged cooperative T, private native CPU E/DP opcode control sharing and exact-hash IEEE32/intergroup/bulk-zero Default provider. Original public35 modules, all40 fields,7RHS/3RK/4quartics, modes, refusals/rollback, LocalSize64, cases/counts, cache/governors/assertions and limits unchanged. C8/all32 first, then ONE A1/B1/B2/A2 with3 internal repeats and original22+Decimal; remaining7 only afterall22; original cold frame only afterall8. No repeated rejected AMDTransport screen, pooling/subtraction/resampling/adoption or scientific/release inference.')
    actions[label] = action
p = W / 'runtime/actions.json'
assert not p.exists()
p.write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'source': revision, 'frozen_inputs': len(inputs),
                  'actions': len(actions), 'order': ['a1', 'b1', 'b2', 'a2']}))
