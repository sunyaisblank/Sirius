"""Contingent value-buffer ABI for compound Add/Multiply calls; no workloads."""
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
inputs = {str(path.relative_to(R)): pin(path.read_bytes()) for path in
          (basis_path, Path(__file__).resolve(), W / 'values_bridge.cpp')}
declarations = '''
// Explicit initialized word buffers avoid assumptions about varying-struct layout.
extern "C" void sirius_native_dense_values_add8(uniform uint32_t * uniform a,
    uniform uint32_t * uniform b, uniform uint32_t * uniform output, uniform uint32_t active);
extern "C" void sirius_native_dense_values_multiply8(uniform uint32_t * uniform a,
    uniform uint32_t * uniform b, uniform uint32_t * uniform output, uniform uint32_t active);
'''
stages = {}
for slug, entry in basis['stages'].items():
    path = R / entry['derived_path']
    data = path.read_bytes()
    assert pin(data) == entry['derived']
    original = data.decode()
    marker = 'typedef int32 int32_t;\n'
    assert original.count(marker) == 1
    candidate = original.replace(marker, marker + declarations, 1)
    changes = []
    for helper, target in (('RTAdd_0', 'sirius_native_dense_values_add8'),
                           ('RTMultiply_0', 'sirius_native_dense_values_multiply8')):
        old = function(candidate, helper)
        signature = old[:old.index('{')]
        a,b = re.search(r'\(varying RetainedTriple_0 \* uniform (\w+), varying RetainedTriple_0 \* uniform (\w+)\)', signature).groups()
        body = '''{
    uniform uint32_t SiriusA[5 * programCount];
    uniform uint32_t SiriusB[5 * programCount];
    uniform uint32_t SiriusOutput[5 * programCount];
    memset(SiriusA, 0, sizeof(SiriusA));
    memset(SiriusB, 0, sizeof(SiriusB));
    memset(SiriusOutput, 0, sizeof(SiriusOutput));
'''
        for i, field in enumerate(('hi_0','lo_0','tail_0','error_0','valid_0')):
            for source,destination in ((a,'SiriusA'),(b,'SiriusB')):
                value = source + '->' + field
                if i != 4:
                    value = 'RPBits_0(' + value + ')'
                body += f'    {destination}[{i} * programCount + programIndex] = {value};\n'
        body += f'    {target}(SiriusA, SiriusB, SiriusOutput, (uniform uint32_t)lanemask());\n'
        body += '    RetainedTriple_0 SiriusResult;\n'
        for i, field in enumerate(('hi_0','lo_0','tail_0','error_0','valid_0')):
            value = f'SiriusOutput[{i} * programCount + programIndex]'
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
    destination = R / 'bin/linux-gcc/native-cpu-dense-values/ispc-source' / (slug+'.ispc')
    publish(destination,candidate.encode())
    inputs[str(path.relative_to(R))] = pin(data)
    stages[slug] = {'basis_path':str(path.relative_to(R)),'basis':pin(data),
                   'derived_path':str(destination.relative_to(R)),'derived':pin(candidate.encode()),
                   'instructions':len(calls(original)),'inverse_restores_exact_AC_source':True}
assert sum(v['instructions'] for v in stages.values()) == 9380
record = {'status':'contingent source only; not compiled or executed','instructions':9380,
          'stages':stages,'source_inputs':inputs,
          'scope':'Route remaining compound Add/Multiply calls to the exact checked AC dense arithmetic through initialized explicit five-column word buffers and the live lane mask. AC direct scratch Store calls remain. New C++ adapters only snapshot a/b/old-output, invoke the existing checked scratch ABI and copy its destination back. No arithmetic source or original dense object changes. All ordered tuples and other original row/control/refusal/rollback/publication remain; inverse restores exact AC source.',
          'reasoning':['Packing writes the original five operand words for each executing lane into fully initialized 8-lane buffers; no varying-struct layout assumption. All snapshots precede external writes.',
                       'The C++ adapter snapshots both inputs and old output into distinct local bank registers before invoking the unchanged AC helper. The existing helper contract changes only active destination words; copying that bank destination back preserves inactive output, including aliases.',
                       'Exact original dense object/opcodes and pinned FP scope must be reused. Fresh actual linked targets/full38/C8 and the original frame remain necessary. Source reasoning is not machine-checked proof or independent math/scientific qualification.',
                       'Compile/run only after a preserved AC original-frame rejection; discard unused contingent output if AC succeeds.']}
publish(W/'ispc-derivation.json',(json.dumps(record,indent=2)+'\n').encode())
print(json.dumps({'source_stages':3,'instructions':9380,'compiled_or_executed':False}))
