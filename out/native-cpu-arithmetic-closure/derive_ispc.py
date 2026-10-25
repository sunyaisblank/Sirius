"""Pack compound operands directly into one initialized scratch bank."""
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


basis_path = R / 'out/native-cpu-masked-arithmetic/ispc-derivation.json'
basis = json.loads(basis_path.read_text())
probe=json.loads((W/'cpp/run/stdout.log').read_text())
assert probe['passed'] and probe['compared_words']==2851200 and probe['calls']==23760
assert json.loads((W/'cpp/run/owner.json').read_text())['passed']
assert json.loads((W/'cpp-code-inspection.json').read_text())['passed']
inputs = {str(path.relative_to(R)): pin(path.read_bytes()) for path in
          (basis_path, Path(__file__).resolve(), W / 'retained8.cpp', W / 'cpp-code-inspection.json', W / 'cpp-build-result.json', W / 'cpp/run/owner.json', W / 'cpp/run/stdout.log',
           R/'src/sirius/kernels/retained_triple.slang',R/'src/sirius/kernels/retained_pair.slang',R/'src/sirius/kernels/retained_scalar.slang')}
stages = {}
for slug, entry in basis['stages'].items():
    path = R / entry['derived_path']
    data = path.read_bytes()
    assert pin(data) == entry['derived']
    original = data.decode()
    marker = 'typedef int32 int32_t;\n'
    assert original.count(marker) == 1
    declarations = '''
// One boundary per complete original division/sqrt refinement sequence.
extern "C" void sirius_native_dense_divide8(uniform uint32_t * uniform scratch,
    uniform uint32_t destination, uniform uint32_t a, uniform uint32_t b, uniform uint32_t active);
extern "C" void sirius_native_dense_sqrt8(uniform uint32_t * uniform scratch,
    uniform uint32_t destination, uniform uint32_t a, uniform uint32_t b, uniform uint32_t active);
'''
    candidate = original.replace(marker,marker+declarations,1)
    changes = []
    for helper, target in (('RTDivide_0', 'sirius_native_dense_divide8'),
                           ('RTSqrt_0', 'sirius_native_dense_sqrt8')):
        old = function(candidate, helper)
        signature = old[:old.index('{')]
        arguments=re.findall(r'varying RetainedTriple_0 \* uniform (\w+)',signature)
        assert len(arguments)==(2 if helper=='RTDivide_0' else 1)
        body = '''{
    uniform uint32_t SiriusBank[15 * programCount];
    memset(SiriusBank, 0, sizeof(SiriusBank));
'''
        for i, field in enumerate(('hi_0','lo_0','tail_0','error_0','valid_0')):
            for operand,source in enumerate(arguments):
                offset=5*operand
                value = source + '->' + field
                if i != 4:
                    value = 'RPBits_0(' + value + ')'
                body += f'    SiriusBank[{offset+i} * programCount + programIndex] = {value};\n'
        body += f'    {target}(SiriusBank, 2, 0, 1, (uniform uint32_t)lanemask());\n'
        body += '    RetainedTriple_0 SiriusResult;\n'
        for i, field in enumerate(('hi_0','lo_0','tail_0','error_0','valid_0')):
            value = f'SiriusBank[{10+i} * programCount + programIndex]'
            if i != 4:
                value = 'RPFromBits_0(' + value + ')'
            body += f'    SiriusResult.{field} = {value};\n'
        body += '    return SiriusResult;\n}'
        new = signature + body
        assert candidate.count(old) == 1
        candidate = candidate.replace(old,new,1)
        changes.append((old,new))
    inverse = candidate
    for old,new in reversed(changes):
        assert inverse.count(new) == 1
        inverse = inverse.replace(new,old,1)
    inverse = inverse.replace(marker+declarations,marker,1)
    assert inverse == original
    calls = lambda s: [line for line in s.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
    assert calls(candidate) == calls(original)
    destination = R / 'bin/linux-gcc/native-cpu-arithmetic-closure/ispc-source' / (slug+'.ispc')
    publish(destination,candidate.encode())
    inputs[str(path.relative_to(R))] = pin(data)
    stages[slug] = {'basis_path':str(path.relative_to(R)),'basis':pin(data),
                   'derived_path':str(destination.relative_to(R)),'derived':pin(candidate.encode()),
                   'instructions':len(calls(original)),'inverse_restores_exact_AE_source':True}
assert sum(v['instructions'] for v in stages.values()) == 9380
record = {'status':'source derived; not renderer-qualified', 'instructions':9380,
          'stages':stages,'source_inputs':inputs,
          'scope':'Only original pure RTDivide/RTSqrt bodies route one initialized scratch bank to the checked complete C++ refinement closure. The same three corrections, final residual, original seed/exponent/scaling/rounding/refusal and five-word results remain; eleven nested compound Add/Multiply pack/ABI/unpack boundaries become in-unit compositions on a full refinement path. Existing AE Store and compound Add/Multiply calls remain. All9380tuples and other row/control/refusal/rollback/publication/FP/scheduling contracts are unchanged; inverse restores exact AE source. Fresh23760call/2851200word comparison is supporting finite compatibility, not independent numerical proof; actual linked joins/full38/C8/original frame remain required.'}
publish(W/'ispc-derivation.json',(json.dumps(record,indent=2)+'\n').encode())
print(json.dumps({'source_stages':3,'instructions':9380,'compiled_or_executed':False}))
