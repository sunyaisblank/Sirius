"""Read the two already captured complete ordinary modules; never dispatch."""
from pathlib import Path
import collections, hashlib, json, re, subprocess
R = Path.cwd(); W = Path(__file__).resolve().parent
compiler = json.loads((W.parent / 'compiler-result.json').read_text())
def normalized_function(text, name):
    body = re.search(r'^[ \t]*%' + name + r' = OpFunction .*?OpFunctionEnd', text, re.M | re.S)[0]
    local_ids = set(re.findall(r'^[ \t]*(%\S+) =', body, re.M))
    ids = {}
    def rename(match):
        token = match[0]
        if token not in ids: ids[token] = '%N' + str(len(ids))
        return ids[token]
    renamed = re.sub(r'%[A-Za-z0-9_]+', lambda m: rename(m) if m[0] in local_ids or m[0][1:].isdigit() else m[0], body)
    return '\n'.join(line.strip() for line in renamed.splitlines())

def expanded_fences(text):
    functions = {}
    for match in re.finditer(r'^[ \t]*(%\S+) = OpFunction .*?OpFunctionEnd', text, re.M | re.S):
        body = match[0]
        functions[match[1]] = (body.count('OpControlBarrier'),
            re.findall(r'OpFunctionCall %\S+ (%\S+)', body))
    memo = {}
    def count(name, active):
        if name in memo: return memo[name]
        assert name not in active
        own, called = functions[name]
        memo[name] = own + sum(count(c, active | {name}) for c in called)
        return memo[name]
    return count('%ComputeMain', set())

rows = {}
for label, kind, stem in [('transport', 'Transport', 'transport'), ('endpoint', 'Endpoint', 'endpoint')]:
    F = W / label
    verification = json.loads((F / 'capture-verification.json').read_text())
    bindings = json.loads((F / 'bindings.json').read_text())
    event, = [x for x in bindings['jit_events'] if x.get('file') == verification['actual_entry']['file']]
    section, = [x for x in event['live_sections_after_pointer_lookup'] if x['name'] == '.ltext']
    raw = (F / section['path']).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == section['sha256']
    regions = [(x['address'], x['address'] + x['bytes']) for x in event['functions']]
    dump = subprocess.check_output(['/usr/bin/objdump', '-D', '-b', 'binary', '-m', 'i386:x86-64',
        '-M', 'intel', '--adjust-vma=' + str(section['address']), str(F / section['path'])], text=True)
    lines = []; ops = []
    for line in dump.splitlines():
        m = re.match(r'\s*([0-9a-f]+):\s+(?:[0-9a-f]{2}\s+)+\s*([a-z][a-z0-9]*)\s*(.*)', line)
        if m and any(a <= int(m[1], 16) < b for a, b in regions):
            lines.append(line); ops.append(m[2])
    (W / (label + '-called-functions.asm')).write_text('\n'.join(lines) + '\n')
    counts = collections.Counter(ops)
    ir = (F / verification['ordinary_ir_export']['path']).read_text()
    assert hashlib.sha256(ir.encode()).hexdigest() == verification['ordinary_ir_export']['sha256']
    shader = (R / 'bin/linux-gcc/native-original-input-ownership/modules' / ('owned-' + stem + '.spvasm')).read_text()
    record = compiler['modules'][kind]
    assert verification['actual_complete_module']['sha256'] == record['sha256']
    assert verification['actual_entry']['actual_mxcsr'] & 0xe040 == 0
    locked = sum(bool(re.search(r'\block or\s+DWORD PTR \[rsp[^]]*\],0x0\b', x)) for x in lines)
    ir_fences = len(re.findall(r'^\s*fence seq_cst\b', ir, re.M))
    expected_fences = expanded_fences(shader)
    assert counts['mfence'] + locked == ir_fences == expected_fences
    assert counts['vaddps'] > 0 and counts['vmulps'] > 0
    rows[kind] = {'complete_module': verification['actual_complete_module'],
        'functions': event['functions'], 'total_called_function_bytes': sum(x['bytes'] for x in event['functions']),
        'source_module_barriers': record['candidate_control_barriers'],
        'expanded_shader_barriers': expected_fences, 'ordinary_IR_fences': ir_fences,
        'machine_full_fences': counts['mfence'] + locked,
        'packed_add': counts['vaddps'], 'packed_multiply': counts['vmulps'],
        'packed_divide': counts['vdivps'], 'packed_sqrt': counts['vsqrtps'],
        'callback_state_restoration': verification['callback_restoration_observed']}

endpoint = (R / compiler['modules']['Endpoint']['module']).with_suffix('.spvasm').read_text()
tools = compiler['tools']
old = subprocess.check_output([tools['disassembler'], str(R / 'bin/linux-gcc/src/sirius/backend/retained/retained_endpoint.spv')], text=True)
assert normalized_function(endpoint, 'Project') == normalized_function(old, 'Project')
# The complete curved call is in the only true block of the built-in lane-zero
# test; the false branch joins after all seven published words and their barrier.
guards = list(re.finditer(r'(%\S+) = OpLoad %uint %gl_LocalInvocationIndex\s+(%\S+) = OpIEqual %bool \1 %uint_0\s+OpSelectionMerge (%\S+) None\s+OpBranchConditional \2 (%\S+) \3\s+\4 = OpLabel\s+(.*?)\s+\3 = OpLabel\s+OpControlBarrier', endpoint, re.S))
guards = [g for g in guards if re.search(r'OpFunctionCall %bool %Project %rowBase %component %root', g[5])]
assert len(guards) == 1
guard, = guards
assert re.search(r'OpFunctionCall %bool %Project %rowBase %component %root', guard[5])
for i in range(7): assert re.search(r'OpAccessChain %\S+ %nativeProjection %int_' + str(i) + r'\b', guard[5])
assert all(compiler['modules'][k]['LocalSize'] == [64, 1, 1] for k in rows)
assert all(compiler['modules'][k]['actual_shared_bytes'] <= 16384 for k in ['Transport', 'Endpoint'])
transport = (R / compiler['modules']['Transport']['module']).with_suffix('.spvasm').read_text()
read_original = re.search(r'^\s*%ReadOriginal = OpFunction .*?OpFunctionEnd', transport, re.M | re.S)[0]
assert 'originalTerms' in read_original and 'RTPack' not in read_original
read_input = re.search(r'^\s*%ReadInput = OpFunction .*?OpFunctionEnd', transport, re.M | re.S)[0]
assert '%ReadStored' in read_input and '%ReadOriginal' in read_input
assert 'OpAtomicUMin' in endpoint and '%firstInvalidOriginal' in endpoint
old_transport = subprocess.check_output([tools['disassembler'], str(R / 'bin/linux-gcc/src/sirius/backend/retained/retained_transport.spv')], text=True)
assert normalized_function(transport, 'ReadOriginalInput') == normalized_function(old_transport, 'ReadOriginal')
out = {'structural_premise_observed': True, 'rows': rows,
    'original_native_Project_SPIRV_function_equal_after_function_local_ID_renaming': True,
    'curved_projection_gated_by_actual_LocalInvocationIndex_zero': True,
    'all7_projection_words_join_before_consumption': True,
    'same_workgroups_and_programs_and_scalar_algorithms': True,
    'Transport_normalized_original_is_complete_immutable_shared_value': True,
    'Transport_mutable_current_state_ReadInput_preserved': True,
    'Endpoint_first_invalid_uses_workgroup_atomic_min': True,
    'actual_shared_bytes': {k: compiler['modules'][k]['actual_shared_bytes'] for k in ['Transport', 'Endpoint']},
    'scope': 'Complete actual changed modules, ordinary optimized LLVM and called code. Four-byte guard dispatch only; no retained numerical, benefit, cold-frame or release acceptance.'}
assert not (W / 'ordinary-code-review.json').exists()
(W / 'ordinary-code-review.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
