"""Replace two private arithmetic-store bodies with the checked dense C++ ABI."""
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
    match = re.search(r'^static (?:inline |noinline )?[^\n]*\b' + name + r'\(', source, re.M)
    assert match, name
    opening = source.index('{', match.end())
    end, depth = opening + 1, 1
    while depth:
        depth += (source[end] == '{') - (source[end] == '}')
        end += 1
    return source[match.start():end]


basis_path = R / 'out/native-cpu-predicated-products/ispc-derivation.json'
basis = json.loads(basis_path.read_text())
inputs = {str(path.relative_to(R)): pin(path.read_bytes()) for path in
          (basis_path, Path(__file__).resolve(), W / 'retained8.cpp', W / 'cpp-code-inspection.json',
           W / 'cpp-build-result.json', W / 'cpp/build/owner.json', W / 'cpp/run/owner.json', W / 'cpp/run/stdout.log')}
probe = json.loads((W / 'cpp/run/stdout.log').read_text())
assert probe['passed'] and probe['calls'] == 11880 and probe['compared_words'] == 1425600
assert probe['input_alias_and_inactive_words_checked'] and probe['FP_controls_unchanged']
assert json.loads((W / 'cpp-code-inspection.json').read_text())['passed']
for phase in ('build', 'run'):
    owner = json.loads((W / 'cpp' / phase / 'owner.json').read_text())
    assert owner['status'] == 'terminal' and owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
declarations = '''
// Scalar C ABI into shared dense C++26 AVX2 arithmetic. The caller's live
// ISPC mask is explicit; inactive scratch words are merged unchanged.
extern "C" void sirius_native_dense_add8(uniform uint32_t * uniform scratch,
    uniform uint32_t destination, uniform uint32_t a, uniform uint32_t b, uniform uint32_t active);
extern "C" void sirius_native_dense_multiply8(uniform uint32_t * uniform scratch,
    uniform uint32_t destination, uniform uint32_t a, uniform uint32_t b, uniform uint32_t active);
'''
stages = {}
for slug, entry in basis['stages'].items():
    path = R / entry['derived_path']
    data = path.read_bytes()
    assert pin(data) == entry['derived']
    original = data.decode()
    candidate = original
    marker = 'typedef int32 int32_t;\n'
    assert candidate.count(marker) == 1
    candidate = candidate.replace(marker, marker + declarations, 1)
    changes = []
    for helper, target in (('RTAddStore_0', 'sirius_native_dense_add8'),
                           ('RTMultiplyStore_0', 'sirius_native_dense_multiply8')):
        old = function(candidate, helper)
        signature = old[:old.index('{')]
        assert all(name in signature for name in ('SiriusDestination', 'SiriusA', 'SiriusB', 'SiriusRowBase', 'SiriusContext'))
        new = signature + '''{
    ''' + target + '''((uniform uint32_t * uniform)&SiriusContext->retainedScratch_0[0],
        SiriusDestination, SiriusA, SiriusB, (uniform uint32_t)lanemask());
}'''
        assert candidate.count(old) == 1
        candidate = candidate.replace(old, new, 1)
        changes.append((old, new))
    address = function(original, 'ScratchAddressUniform_0')
    parameters = re.findall(r'uniform uint32_t (\w+)', address.split('(', 1)[1].split(')', 1)[0])
    assert len(parameters) == 4
    assert re.search(r'return\s+' + parameters[0] + r'\s*\*\s*5U\s*\+\s*' + parameters[1] + r'\s*;', address)
    assert original.count('memset((void * uniform)&context, 0, sizeof(context));') == 1
    scratch, = re.findall(r'uint32_t retainedScratch_0\[(\d+)\];', original)
    assert int(scratch) % 5 == 0
    inverse = candidate
    for old, new in reversed(changes):
        assert inverse.count(new) == 1
        inverse = inverse.replace(new, old, 1)
    inverse = inverse.replace(marker + declarations, marker, 1)
    assert inverse == original
    calls = lambda text: [line for line in text.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
    assert calls(candidate) == calls(original)
    destination = R / 'bin/linux-gcc/native-cpu-dense-arithmetic/ispc-source' / (slug + '.ispc')
    publish(destination, candidate.encode())
    inputs[str(path.relative_to(R))] = pin(data)
    stages[slug] = {'basis_path': str(path.relative_to(R)), 'basis': pin(data),
                    'derived_path': str(destination.relative_to(R)), 'derived': pin(candidate.encode()),
                    'scratch_columns': int(scratch), 'instructions': len(calls(original)),
                    'inverse_restores_exact_AB_source': True}
assert sum(v['instructions'] for v in stages.values()) == 9380
record = {'status': 'source bridge prepared; not compiled or renderer-qualified',
          'instructions': 9380, 'stages': stages, 'source_inputs': inputs,
          'scope': 'Only two shared private arithmetic-store bodies call the checked C++26 dense AVX2 implementation via scalar C ABI. Every ordered canonical tuple, all original ISPC arithmetic helpers/generic fallback, row/control/refusal/rollback/publication remain. C++ operands are loaded in full before writes; each destination field merges the explicit active lane mask. Actual shared scratch layout is five varying32 columns per register, eight physical words per column, with complete context initialization before use. Existing original index guards dominate calls. No eligibility/count/precision/physical operation changes.',
          'proof_limits': 'The bridge inverse restores exact AB source. Dense arithmetic equivalence is coordinator source reasoning plus one focused finite scalar-reference comparison, not machine-checked proof or independent scientific validation. ISPC uniform pointer reinterpretation/lanemask ABI, linked call targets and complete row integration still require actual code association and fresh full38/C8 before the sole original frame.'}
publish(W / 'ispc-derivation.json', (json.dumps(record, indent=2) + '\n').encode())
print(json.dumps({'source_stages': 3, 'dense_C_ABI_calls_each': 2, 'instructions': 9380,
                  'renderer_compiled_or_executed': False}))
