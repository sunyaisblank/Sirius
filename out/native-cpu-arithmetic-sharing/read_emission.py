"""Inspect emitted ownership/storage and unchanged scalar functions; no dispatch."""
from pathlib import Path
import hashlib, json, re, subprocess
R = Path.cwd(); W = Path(__file__).resolve().parent
C = R / 'bin/linux-gcc/src/sirius/backend/retained'
d = json.loads((W / 'compiler-result.json').read_text())
def function(s, name):
    body = re.search(r'^\s*%' + name + r' = OpFunction .*?OpFunctionEnd', s, re.M | re.S)[0]
    local = set(re.findall(r'^\s*(%\S+) =', body, re.M)); ids = {}
    def rename(m):
        v = m[0]
        if v not in ids: ids[v] = '%N' + str(len(ids))
        return ids[v]
    body = re.sub(r'%[A-Za-z0-9_]+', lambda m: rename(m) if m[0] in local or m[0][1:].isdigit() else m[0], body)
    return '\n'.join(x.strip() for x in body.splitlines())
out = {}
for kind, stem in [('Endpoint', 'endpoint'), ('DopriPhase', 'dopri_phase')]:
    row = d['modules'][kind]; p = R / row['module']; s = p.with_suffix('.spvasm').read_text()
    assert hashlib.sha256(p.read_bytes()).hexdigest() == row['sha256']
    old = subprocess.check_output([d['tools']['disassembler'], str(C / ('retained_' + stem + '.spv'))], text=True)
    assert re.findall(r'OpCapability (.*)', s) == re.findall(r'OpCapability (.*)', old)
    assert re.findall(r'OpExecutionMode (.*)', s) == re.findall(r'OpExecutionMode (.*)', old)
    assert 'LocalSize 64 1 1' in s and 'DenormPreserve 32' in s and 'RoundingModeRTE 32' in s
    assert 'OpControlBarrier' not in s and 'OpMemoryBarrier' not in s and 'OpAtomic' not in s
    assert not re.findall(r'OpVariable %\S+ Workgroup', s)
    constants = {k: int(v) for k, v in re.findall(r'(%\S+) = OpConstant %(?:u?int) (\d+)', s)}
    sizes = {k: 4 for k in re.findall(r'(%\S+) = OpTypeInt 32 [01]', s)}
    for k, element, length in re.findall(r'(%\S+) = OpTypeArray (%\S+) (%\S+)', s):
        if element in sizes and length in constants: sizes[k] = sizes[element] * constants[length]
    pointer = dict(re.findall(r'(%\S+) = OpTypePointer Private (%\S+)', s))
    private = {k: sizes[pointer[t]] for k, t in re.findall(r'(%\S+) = OpVariable (%\S+) Private', s)}
    required = {'%retainedScratch', '%retainedStatus'}
    required |= {'%retainedCoefficients'} if kind == 'Transport' else {'%nativeProjection'} if kind == 'Endpoint' else set()
    assert set(private) == required and sum(private.values()) <= 16384
    # No implicit zero initializer on these invocation-private arrays/status.
    assert not re.findall(r'OpVariable %\S+ Private %\S+', s)
    entry = re.search(r'^\s*%ComputeMain = OpFunction .*?OpFunctionEnd', s, re.M | re.S)[0]
    assert re.search(r'(%\S+) = OpLoad %uint %gl_LocalInvocationIndex\s+(%\S+) = OpIEqual %bool \1 %uint_0\s+OpSelectionMerge (%\S+) None\s+OpBranchConditional \2 (%\S+) \3\s+\4 = OpLabel\s+(%\S+) = OpLoad %v3uint %gl_WorkGroupID\s+%\S+ = OpFunctionCall %void %ComputeOwnedRow \5 %uint_0\s+OpBranch \3', entry)
    unchanged = ['Project'] if kind == 'Endpoint' else ['ReadOriginal']
    for name in unchanged: assert function(s, name) == function(old, name), (kind, name)
    helper = function(s, 'RTEvaluateArithmetic')
    assert helper.count('OpFunctionCall %RetainedTriple %RTMultiply ') == 1
    assert helper.count('OpFunctionCall %RetainedTriple %RTAdd ') == 1
    assert function(s,'EvaluateOne').count('OpFunctionCall %RetainedTriple %RTEvaluateArithmetic ') == 1
    out[kind] = {'module': row, 'actual_Private_variables_bytes': private,
        'actual_Private_total_bytes': sum(private.values()), 'Workgroup_variables': 0,
        'barriers': 0, 'atomics': 0, 'entry_only_lane0_calls_complete_row_with_literal_lane0': True,
        'LocalSize': [64, 1, 1], 'original_capabilities_and_modes_equal': True,
        'unchanged_scalar_functions': unchanged, 'implicit_Private_initializers': 0}
p = W / 'emission-result.json'; assert not p.exists()
p.write_text(json.dumps({'passed': True, 'modules': out,
    'scope': 'Emitted SPIRV/storage and scalar function identity only. Actual LLVM private lifetime, code, numerical/benefit/frame qualification still required.'}, indent=2) + '\n')
print(json.dumps({'passed': True, 'Private_bytes': {k: v['actual_Private_total_bytes'] for k, v in out.items()},
    'both_actual_lane0_guards': True, 'both_Workgroup_variables_and_barriers_and_atomics': 0}))
