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


basis_path = R / 'out/native-cpu-dense-arithmetic/ispc-derivation.json'
basis = json.loads(basis_path.read_text())
probe=json.loads((W/'cpp/run/stdout.log').read_text())
assert probe['passed'] and probe['compared_words']==1425600
assert json.loads((W/'cpp/run/owner.json').read_text())['passed']
assert json.loads((W/'cpp-code-inspection.json').read_text())['passed']
inputs = {str(path.relative_to(R)): pin(path.read_bytes()) for path in
          (basis_path, Path(__file__).resolve(), W / 'retained8.cpp', W / 'cpp-code-inspection.json', W / 'cpp-build-result.json', W / 'cpp/run/owner.json', W / 'cpp/run/stdout.log')}
stages = {}
for slug, entry in basis['stages'].items():
    path = R / entry['derived_path']
    data = path.read_bytes()
    assert pin(data) == entry['derived']
    original = data.decode()
    marker = 'typedef int32 int32_t;\n'
    assert original.count(marker) == 1
    candidate = original
    changes = []
    for helper, target in (('RTAdd_0', 'sirius_native_dense_add8'),
                           ('RTMultiply_0', 'sirius_native_dense_multiply8')):
        old = function(candidate, helper)
        signature = old[:old.index('{')]
        a,b = re.search(r'\(varying RetainedTriple_0 \* uniform (\w+), varying RetainedTriple_0 \* uniform (\w+)\)', signature).groups()
        body = '''{
    uniform uint32_t SiriusBank[15 * programCount];
    memset(SiriusBank, 0, sizeof(SiriusBank));
'''
        for i, field in enumerate(('hi_0','lo_0','tail_0','error_0','valid_0')):
            for source,offset in ((a,0),(b,5)):
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
    assert inverse == original
    calls = lambda s: [line for line in s.splitlines() if line.strip().startswith('if(!EvaluateOneStatic_0<')]
    assert calls(candidate) == calls(original)
    destination = R / 'bin/linux-gcc/native-cpu-masked-arithmetic/ispc-source' / (slug+'.ispc')
    publish(destination,candidate.encode())
    inputs[str(path.relative_to(R))] = pin(data)
    stages[slug] = {'basis_path':str(path.relative_to(R)),'basis':pin(data),
                   'derived_path':str(destination.relative_to(R)),'derived':pin(candidate.encode()),
                   'instructions':len(calls(original)),'inverse_restores_exact_AC_source':True}
assert sum(v['instructions'] for v in stages.values()) == 9380
record = {'status':'source derived; not renderer-qualified', 'instructions':9380,
          'stages':stages,'source_inputs':inputs,
          'scope':'AE replaces only the two pure compound Add/Multiply bodies of AC source with one initialized 15-column scratch bank and direct calls to the same explicit mask/index C ABI as Store calls. Both operands are packed before external writes; output occupies distinct initialized register2. New checked AE math eliminates discarded-lane FP work and private aggregate helper calls under unchanged strict FP flags. Full inverse restores exact AC source; all9380tuples and other row/control/refusal/publication contracts remain. No packing-only C++ wrapper or copies are needed. Source reasoning and finite scalar agreement are not independent mathematical/scientific acceptance; linked targets/full38/C8/original frame remain required.'}
publish(W/'ispc-derivation.json',(json.dumps(record,indent=2)+'\n').encode())
print(json.dumps({'source_stages':3,'instructions':9380,'compiled_or_executed':False}))
