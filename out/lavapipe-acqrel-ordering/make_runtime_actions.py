"""Keep the exact original control cases, order, limits and finite admission."""
from pathlib import Path
import json, re, subprocess

R = Path.cwd(); W = Path(__file__).resolve().parent
assert json.loads((W / 'ordinary/ordinary-code-review.json').read_text())['structural_premise_observed']
old = json.loads((R / 'out/native-original-input-ownership/runtime/actions.json').read_text())
providers = json.loads((W / 'providers.json').read_text())['providers']
inputs = set(old['control-8']['frozen_inputs'])
inputs |= {'out/native-cpu-row-ownership/cold-frame/run_owned.py',
           'out/native-cpu-row-ownership/cold-frame/owned_guard.py'}
for folder in [W, W / 'runtime', W / 'readers']:
    inputs |= {str(p.relative_to(R)) for p in folder.iterdir() if p.is_file()}
inputs |= {str(p.relative_to(R)) for p in (W / 'ordinary').rglob('*') if p.is_file()}
inputs |= {str(p.relative_to(R)) for p in (W / 'build').iterdir() if p.is_file()}
for value in providers.values():
    inputs |= {value['library'], value['ICD']}
    output = subprocess.check_output(['ldd', value['library']], text=True)
    inputs |= {str(Path(x).resolve()) for x in re.findall(r'(/[^\s()]+)', output) if Path(x).is_file()}
for name in inputs: assert (R / name).is_file(), name
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
actions = {}
for label, basis in old.items():
    action = {k: basis[k] for k in ['arm', 'case', 'prerequisites', 'kind', 'expected_exit', 'limits']}
    provider = providers[action['arm']]
    action.update(revision=revision, argv=['/usr/bin/python3', '-B', str(W / 'run_case.py'), label, '{folder}'],
        frozen_inputs=sorted(inputs), software_provider=True,
        runtime_icd=provider['ICD'], runtime_provider=provider['library'],
        scope='Original ' + label + ', original same admitted complete-input matched ELF and all exact shader/program bytes; only candidate provider exact acquire-release lowering changes. Same actual stock baseline, cooperative64-lane workgroups, coroutine rendezvous, arithmetic/fields/refusal/rollback, counts, default product mode, cache policy, governors, assertions and original limits. C8 then sole A1/B1/B2/A2 with3 repeats/slot and original22 thresholds; remaining7 only afterall22 PASS and original changed cold frame only afterall8 PASS. Historical hardware/precision expectations are composed with actual new provider features/properties/code/boundaries and fresh runtime identities. No unchanged retry, pooling, subtraction or scientific/release/adoption inference.')
    actions[label] = action
p = W / 'runtime/actions.json'; assert not p.exists()
p.write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'source': revision, 'frozen_inputs': len(inputs), 'same_existing_ELF': True,
    'new_CPP_or_shader_compilations': 0, 'fixed_benefit_order': ['a1', 'b1', 'b2', 'a2']}))
