"""Inline the two product helpers; inverse must restore exact S source."""
from pathlib import Path
import hashlib
import json
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

old_path = R / 'out/native-cpu-zero-term/ispc-derivation.json'
old = json.loads(old_path.read_text())
inputs = {str(old_path.relative_to(R)): pin(old_path.read_bytes()),
          str(Path(__file__).resolve().relative_to(R)): pin(Path(__file__).read_bytes())}
stages = {}
calls = 0
for slug, entry in old['stages'].items():
    basis = R / entry['derived_path']
    data = basis.read_bytes()
    assert pin(data) == entry['derived']
    source = data.decode()
    candidate = source
    for name in ('RTAccumulateProduct_0', 'RTAccumulateRoundedProduct_0'):
        before = 'static RetainedTriple_0 ' + name + '('
        assert candidate.count(before) == 1
        candidate = candidate.replace(before, 'static inline RetainedTriple_0 ' + name + '(', 1)
    inverse = candidate
    for name in ('RTAccumulateProduct_0', 'RTAccumulateRoundedProduct_0'):
        inverse = inverse.replace('static inline RetainedTriple_0 ' + name + '(', 'static RetainedTriple_0 ' + name + '(', 1)
    assert inverse.encode() == data
    literal = lambda text: [line for line in text.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
    assert literal(candidate) == literal(source)
    count = len(literal(candidate)); calls += count
    destination = R / 'bin/linux-gcc/native-cpu-product-inline/ispc-source' / (slug + '.ispc')
    publish(destination, candidate.encode())
    inputs[str(basis.relative_to(R))] = pin(data)
    stages[slug] = {**entry, 'product_basis_path': str(basis.relative_to(R)),
        'product_basis': pin(data), 'derived_path': str(destination.relative_to(R)),
        'derived': pin(candidate.encode()), 'two_inline_qualifiers_only': True,
        'inverse_restores_exact_S_bytes': True, 'literal_instruction_calls_preserved': count}
assert calls == 9380
result = {'instructions': calls, 'stages': stages, 'source_inputs': inputs,
    'scope': 'Only two static-to-static-inline product-helper qualifiers change in each complete S source. All arithmetic bodies, ordered literal instructions, guards, output fields, film/pupil columns, RHS/trials/quartics, rounding and ABI exports retained. Actual S objects contain six calls per helper per stage; removal is a source-supported premise, not a speed claim. No runtime qualification.'}
publish(W / 'ispc-derivation.json', (json.dumps(result, indent=2) + '\n').encode())
print(json.dumps({'new_source_copies': 3, 'instructions': calls, 'qualifiers_per_stage': 2}))
