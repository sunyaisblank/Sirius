"""Bind the changed provider to the original finite controls and sole ABBA."""
from pathlib import Path
import hashlib, json, re, subprocess
R = Path.cwd(); W = Path(__file__).resolve().parent
OLD = R / 'out/native-cooperative-transport-intergroup'
def doc(p): return json.loads(p.read_text())
assert doc(W / 'ordinary/ordinary-code-review.json')['structural_premise_observed']
old = doc(R / 'out/native-original-input-ownership/runtime/actions.json')
providers = doc(W / 'providers.json')['providers']
inputs = set(old['control-8']['frozen_inputs'])
matched = doc(W / 'matched-result.json')
inputs |= {matched['executable'], str((OLD / 'headers/retained_kernels.h').relative_to(R))}
inputs |= {d['module'] for d in doc(W / 'compiler-result.json')['modules'].values()}
inputs |= {d['candidate']['path'] for d in doc(W / 'source-contract.json')['modified_files'].values()}
inputs |= {'out/native-cpu-row-ownership/cold-frame/run_owned.py', 'out/native-cpu-row-ownership/cold-frame/owned_guard.py'}
inputs |= {str(p.relative_to(R)) for p in (OLD / 'matched').iterdir() if p.is_file()}
inputs |= {str((OLD / name).relative_to(R)) for name in ['source-scope.json', 'retained_compute.cpp', 'matched-result.json']}
for folder in [W, W / 'runtime', W / 'readers']:
    inputs |= {str(p.relative_to(R)) for p in folder.iterdir() if p.is_file()}
for folder in [W / 'ordinary', W / 'build']:
    inputs |= {str(p.relative_to(R)) for p in folder.rglob('*') if p.is_file()}
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
        frozen_inputs=sorted(inputs), software_provider=True, runtime_icd=provider['ICD'], runtime_provider=provider['library'],
        scope='Original ' + label + ' on the unchanged exact mixed f958 host and three immutable modules, with the newly scoped post-IR machine CodeGenNone provider for exact Endpoint/DopriPhase row-owner variants. Cooperative Transport keeps ordinary compilation. Same actual stock baseline, original logical groups/LocalSize64, arithmetic/fields/refusal/rollback/counts/default product mode/cache policy/governors/assertions and limits. C8 then exactly one A1/B1/B2/A2 with3 repeats/slot and original22 thresholds; remaining7 only afterall22 PASS, original cold frame only afterall8 PASS. Actual complete code, optimized-IR equality, allocator isolation, features/properties/MXCSR and fresh residency/readbacks required. No unchanged retry, pooling, subtraction, resampling, adoption or scientific/release inference.')
    actions[label] = action
p = W / 'runtime/actions.json'; assert not p.exists()
p.write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'source': revision, 'frozen_inputs': len(inputs), 'same_existing_ELF': True,
    'new_CPP_or_shader_compilations': 0, 'fixed_benefit_order': ['a1', 'b1', 'b2', 'a2']}))
