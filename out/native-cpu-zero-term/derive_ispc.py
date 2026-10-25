"""Derive new complete source copies by inserting only the reviewed zero arm."""
from pathlib import Path
import hashlib
import json
import re
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
assert sys.argv[1:] in ([], ['--write'])


def pin(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def publish(path, data):
    if path.exists():
        assert path.read_bytes() == data, path
    else:
        assert sys.argv[1:] == ['--write']
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            stream.write(data)


proof = json.loads((W / 'source-proof.json').read_text())
previous = json.loads((R / 'out/native-cpu-accumulate-inline/ispc-derivation.json').read_text())
assert proof['instructions'] == previous['instructions'] == 9380
stages, calls = {}, 0
inputs = dict(proof['source_inputs'])
for name in ('prepare_proposal.py', 'source-proof.json', 'zero-term.patch', 'derive_ispc.py'):
    inputs[str((W / name).relative_to(R))] = pin((W / name).read_bytes())
for slug, entry in proof['stages'].items():
    basis = R / entry['basis_path']
    text = basis.read_text()
    assert pin(text.encode()) == entry['basis']
    header = re.search(r'static inline RetainedTriple_0 RTAccumulateTerm_0\(', text)
    assert header
    end = text.index('        return RTInvalid_0();\n    }', header.end()) + len('        return RTInvalid_0();\n    }')
    arm = entry['zero_arm']
    candidate = text[:end] + arm + text[end:]
    assert candidate.replace(arm, '', 1) == text
    assert pin(candidate.encode()) == entry['candidate']
    literal = lambda value: [line for line in value.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
    assert literal(text) == literal(candidate)
    calls += len(literal(text))
    destination = R / entry['future_candidate_path']
    publish(destination, candidate.encode())
    stages[slug] = {**previous['stages'][slug],
        'derived_path': str(destination.relative_to(R)), 'derived': pin(candidate.encode()),
        'zero_term_basis_path': entry['basis_path'], 'zero_term_basis': entry['basis'],
        'inverse_removes_one_zero_arm_only': True, 'zero_term_arithmetic_specialized': True,
        'literal_instruction_calls_preserved': len(literal(text))}
assert calls == 9380
result = {'instructions': calls, 'stages': stages, 'source_inputs': inputs,
          'scope': 'New full source copies for all3 stages. Add only the reviewed exact signed-zero arm after original validity/x-bound refusal; nonzero path unchanged, all9380calls and every field/column/RHS/trial/quartic/count retained. Source derivation only, no compilation or runtime qualification.'}
publish(W / 'ispc-derivation.json', (json.dumps(result, indent=2) + '\n').encode())
print(json.dumps({'new_source_copies': 3, 'instructions': calls, 'workloads_executed': 0}))
