"""Inspect the exact validated optional owner modules, without executing them."""
import hashlib, json, re, sys
from pathlib import Path
sys.dont_write_bytecode = True
work = Path(__file__).resolve().parent
phase = work / 'compiler-screen'
results = {}
for kind, expected in (('transport', 9708), ('endpoint', 12616)):
    folder = phase / kind
    text = (folder / 'final.spvasm').read_text()
    assert re.findall(r'OpCapability (\S+)', text) == ['Shader', 'RoundingModeRTE']
    assert set(re.findall(r'OpTypeInt (\d+) [01]', text)) == {'32'}
    assert re.findall(r'OpTypeFloat (\d+)', text) == ['32']
    assert not re.search(r'FPFastMathMode|FPFastMathDefault|RelaxedPrecision|OpFma|ReadClock|DenormPreserve', text)
    entry = re.search(r'OpEntryPoint GLCompute (%\S+)', text)[1]
    assert f'OpExecutionMode {entry} LocalSize 16 1 1' in text
    assert f'OpExecutionMode {entry} RoundingModeRTE 32' in text
    operations = re.findall(r'(%\S+) = OpF(?:Add|Sub) ', text)
    decorated = set(re.findall(r'OpDecorate (%\S+) NoContraction', text))
    assert operations and all(value in decorated for value in operations)
    bindings = dict(re.findall(r'OpDecorate (%\S+) Binding (\d+)', text))
    sets = dict(re.findall(r'OpDecorate (%\S+) DescriptorSet (\d+)', text))
    assert sorted(bindings.values()) == ['0', '1', '2']
    assert all(sets[value] == '0' for value in bindings)
    names = dict(re.findall(r'OpName (%\S+) "([^"\n]+)"', text))
    assert names[next(value for value, binding in bindings.items() if binding == '2')] == 'stageSevenOwner'
    defs = dict(re.findall(r'^\s*(%\S+) = (.+)$', text, re.M))
    def size(value):
        parts = defs[value].split()
        op = parts[0]
        if op in ('OpTypeInt', 'OpTypeFloat'):
            assert parts[1] == '32'
            return 4
        if op == 'OpTypePointer': return size(parts[2])
        if op == 'OpTypeArray':
            count = defs[parts[2]].split()
            assert count[0] == 'OpConstant'
            return size(parts[1]) * int(count[2], 0)
        if op == 'OpTypeVector': return size(parts[1]) * int(parts[2])
        if op == 'OpTypeStruct': return sum(map(size, parts[1:]))
        raise AssertionError((value, parts))
    shared = {names.get(value, value): size(definition.split()[1])
              for value, definition in defs.items()
              if re.match(r'OpVariable %\S+ Workgroup(?: |$)', definition)}
    assert sum(shared.values()) == expected, shared
    assert {name: shared[name] for name in ('s7Enabled','s7RowBase','s7Family','s7Invalid','s7Captured')} == dict.fromkeys(('s7Enabled','s7RowBase','s7Family','s7Invalid','s7Captured'), 4)
    binary = folder / 'owner.spv'
    results[kind] = {'bytes': binary.stat().st_size, 'sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
                     'shared_variables_bytes': shared, 'shared_bytes': sum(shared.values()),
                     'bindings': {names.get(value,value): int(binding) for value,binding in bindings.items()},
                     'static_control_barriers': text.count('OpControlBarrier'),
                     'guarded_float_operations': len(operations), 'local_size': [16,1,1]}
result = {'accepted': True, 'modules': results,
          'scope': 'Validated compiler modules and structural IR only. No device allocation, runtime replay, gain, scientific or release qualification.'}
(work / 'owner-ir.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
