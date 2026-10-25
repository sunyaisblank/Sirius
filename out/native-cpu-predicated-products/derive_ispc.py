"""Prepare predicated product accumulators; no compilation or workload."""
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


def function(source, name):
    match = re.search(r'^static (?:inline )?[^\n]+\b' + name + r'\(', source, re.M)
    assert match, name
    opening = source.index('{', match.end())
    end, depth = opening + 1, 1
    while depth:
        depth += (source[end] == '{') - (source[end] == '}')
        end += 1
    return source[match.start():end]


basis_path = R / 'out/native-cpu-static-registers/ispc-derivation.json'
basis = json.loads(basis_path.read_text())
inputs = {str(basis_path.relative_to(R)): pin(basis_path.read_bytes()),
          str(Path(__file__).resolve().relative_to(R)): pin(Path(__file__).read_bytes())}
stages = {}
for slug, entry in basis['stages'].items():
    path = R / entry['derived_path']
    data = path.read_bytes()
    assert pin(data) == entry['derived']
    original = data.decode()
    candidate = original
    replacements, proofs = [], {}
    for name in ('RTAccumulateProduct_0', 'RTAccumulateRoundedProduct_0'):
        old = function(original, name)
        signature, body = old.split('{', 1)
        sum_name, a, b = re.search(
            r'\(varying RetainedTriple_0 \* uniform (\w+), float (\w+), float (\w+)\)', signature).groups()
        assert re.findall(r'!RPBounded_0\((\w+), 1568669696U\)', body) == [a, b]
        assert re.findall(r'RPScalarEqual_0\((\w+), \w+\)', body) == [a, b]
        assert body.count('return RTInvalid_0();') == 1 and body.count('return *' + sum_name + ';') == 1
        zero_return = body.index('return *' + sum_name + ';')
        after_zero = body.index('}', zero_return) + 1
        first = re.search(r'^\s*float\s+\w+(?:;|\s*=)', body[after_zero:], re.M)
        assert first
        tail = body[after_zero + first.start():body.rfind('}')]
        returned, = re.findall(r'\breturn (\w+);', tail)
        last = tail.rfind('return ' + returned + ';')
        core = tail[:last]
        assert not re.search(r'\breturn\b', core)
        prefix = f'''{{
    bool SiriusInputGood = ({sum_name}->valid_0 != 0U) &
        RPBounded_0({a}, 1568669696U) & RPBounded_0({b}, 1568669696U);
    bool SiriusProductZero = RPScalarEqual_0({a}, RPZero_0()) |
        RPScalarEqual_0({b}, RPZero_0());
    // Preserve the cheap original zero-product return for all executing lanes.
    if (all(SiriusInputGood & SiriusProductZero))
        return *{sum_name};
'''
        suffix = f'''RetainedTriple_0 SiriusGeneral = {returned};
    RetainedTriple_0 SiriusResult;
'''
        for field in ('hi_0', 'lo_0', 'tail_0', 'error_0'):
            suffix += (f'    SiriusResult.{field} = select(SiriusInputGood, '
                       f'select(SiriusProductZero, {sum_name}->{field}, SiriusGeneral.{field}), RPZero_0());\n')
        suffix += (f'    SiriusResult.valid_0 = select(SiriusInputGood, '
                   f'select(SiriusProductZero, {sum_name}->valid_0, SiriusGeneral.valid_0), 0U);\n'
                   '    return SiriusResult;\n}')
        new = signature + prefix + core + suffix
        assert candidate.count(old) == 1
        candidate = candidate.replace(old, new, 1)
        replacements.append((old, new))
        proofs[name] = {'basis': pin(old.encode()), 'derived': pin(new.encode()),
                        'ordered_nonzero_core_copied_exactly': pin(core.encode()),
                        'invalid_precedes_zero': True, 'all_five_zero_product_sum_words_preserved': True}
    inverse = candidate
    for old, new in reversed(replacements):
        assert inverse.count(new) == 1
        inverse = inverse.replace(new, old, 1)
    assert inverse == original
    calls = lambda text: [line for line in text.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
    assert calls(candidate) == calls(original)
    destination = R / 'bin/linux-gcc/native-cpu-predicated-products/ispc-source' / (slug + '.ispc')
    publish(destination, candidate.encode())
    inputs[str(path.relative_to(R))] = pin(data)
    stages[slug] = {'basis_path': str(path.relative_to(R)), 'basis': pin(data),
                    'derived_path': str(destination.relative_to(R)), 'derived': pin(candidate.encode()),
                    'helpers': proofs, 'instructions': len(calls(original)), 'inverse_restores_exact_AA_source': True}
assert sum(v['instructions'] for v in stages.values()) == 9380
proof = {'status': 'contingent source only; no compilation or workload', 'instructions': 9380,
         'stages': stages, 'source_inputs': inputs,
         'scope': 'Replace only the two product-accumulator input/zero early-return masks with explicit result selection. Original complete nonzero calculation core, including conditional error update, is copied byte-for-byte. Invalid input selects original canonical invalid before zero selection; zero products preserve all five original sum words, including signed zeros and the original nonzero valid word. Full AA static metadata, ordered tuples, arithmetic helpers, row/refusal/rollback/publication remain.',
         'reasoning': [
             'For each executing lane, bad input returns the original canonical invalid. Good zero product returns the exact original sum. Good nonzero product executes and returns the original core. These three cases are disjoint and exhaustive; the uniform all-good-zero return preserves the original cheap common case.',
             'The original nonzero core reads sum and writes only local temporaries; RetainedProduct32, ordered FP/bit helpers and RTAccumulateTerm have no external or indexing effects. Extra calculations on unselected lanes can raise FP status, so this is suitable only inside the existing saved/masked IEEE32 worker scope with full restoration and no intervening status observation. Actual FP ownership/code inspection and runtime comparisons remain required; this source argument is not machine-checked proof.',
             'Do not compile or execute unless AA fails the decisive original frame and its terminal result has been preserved. Discard unused contingent output if AA succeeds.'
         ]}
publish(W / 'ispc-derivation.json', (json.dumps(proof, indent=2) + '\n').encode())
print(json.dumps({'source_stages': 3, 'changed_pure_helpers_each': 2, 'instructions': 9380,
                  'compiled': False, 'executed_workload': False}))
