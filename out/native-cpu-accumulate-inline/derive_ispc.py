"""Copy sealed L sources; add exactly one inline qualifier per stage."""
from pathlib import Path
import hashlib
import json
import sys

root = Path.cwd()
work = Path(__file__).resolve().parent
prior = root / 'out/native-cpu-simd-endpoint-regions'
destination = root / 'bin/linux-gcc/native-cpu-accumulate-inline/ispc-source'
assert sys.argv[1:] in ([], ['--write'])


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def publish(path, data):
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        assert sys.argv[1:] == ['--write'], path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)


proof = json.loads((work / 'source-proof.json').read_text())
old_receipt = json.loads((prior / 'ispc-derivation.json').read_text())
assert proof['instructions'] == old_receipt['instructions'] == 9380 and not proof['patch_applied']
inputs = dict(proof['source_inputs'])
for name in ('source-proof.json', 'accumulate-inline.patch', 'coverage-plan.json'):
    path = work / name
    inputs[str(path.relative_to(root))] = pin(path.read_bytes())
stages, count = {}, 0
for slug in ('transport', 'endpoint', 'dopri_phase'):
    old = old_receipt['stages'][slug]
    basis = root / old['derived_path']
    text = basis.read_text()
    assert pin(text.encode()) == old['derived'] == proof['stages'][slug]['basis']
    header = 'static RetainedTriple_0 RTAccumulateTerm_0('
    replacement = 'static inline RetainedTriple_0 RTAccumulateTerm_0('
    assert text.count(header) == 1 and replacement not in text
    revised = text.replace(header, replacement, 1)
    assert revised.replace(replacement, header, 1) == text
    assert pin(revised.encode()) == proof['stages'][slug]['candidate']
    calls = [line for line in text.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
    assert calls == [line for line in revised.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
    count += len(calls)
    path = destination / (slug + '.ispc')
    stages[slug] = {**old, 'derived_path': str(path.relative_to(root)), 'derived': pin(revised.encode()),
        'inline_basis_path': old['derived_path'], 'inline_basis': old['derived'],
        'one_inline_qualifier_only': True, 'full_inverse_restores_L_source': True,
        'literal_instruction_calls_preserved': len(calls)}
    publish(path, revised.encode())
assert count == 9380
receipt = {'instructions': count, 'stages': stages, 'source_inputs': inputs,
    'derivation_basis': {'generator_path': str((prior / 'derive_ispc.py').relative_to(root)),
        'generator': pin((prior / 'derive_ispc.py').read_bytes()),
        'receipt_path': str((prior / 'ispc-derivation.json').relative_to(root)),
        'receipt': pin((prior / 'ispc-derivation.json').read_bytes())},
    'scope': 'All3 complete sealed L sources, one inline qualifier each. All other bytes, all9380 literal instruction calls, masks/math/control/counts/FP/eligibility retained. Source derivation only; no compilation or numerical/benefit/frame acceptance.'}
publish(work / 'ispc-derivation.json', (json.dumps(receipt, indent=2) + '\n').encode())
for name, expected in inputs.items():
    assert pin((root / name).read_bytes()) == expected, name
print(json.dumps({'source_copies': 3, 'inline_qualifiers': 3, 'instructions': count,
                  'compiler_or_runtime_executed': False}))
