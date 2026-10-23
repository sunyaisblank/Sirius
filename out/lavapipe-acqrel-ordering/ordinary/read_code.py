"""Inspect captured ordinary code; never compile, dispatch or resample."""
from pathlib import Path
import collections, hashlib, json, re, subprocess

R = Path.cwd(); W = Path(__file__).resolve().parent
compiler = json.loads((R / 'out/native-original-input-ownership/compiler-result.json').read_text())
prior = json.loads((R / 'out/native-original-input-ownership/ordinary/ordinary-code-review.json').read_text())
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

rows = {}; observed = True
for label, kind in [('transport', 'Transport'), ('endpoint', 'Endpoint'), ('dopri_phase', 'DopriPhase')]:
    F = W / label
    verification = json.loads((F / 'capture-verification.json').read_text())
    bindings = json.loads((F / 'bindings.json').read_text())
    event, = [x for x in bindings['jit_events'] if x.get('file') == verification['actual_entry']['file']]
    section, = [x for x in event['live_sections_after_pointer_lookup'] if x['name'] == '.ltext']
    raw = (F / section['path']).read_bytes(); assert hashlib.sha256(raw).hexdigest() == section['sha256']
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
    record = compiler['modules'][kind]
    shader = (R / record['module']).with_suffix('.spvasm').read_text()
    assert verification['actual_complete_module']['sha256'] == record['sha256']
    assert verification['actual_entry']['actual_mxcsr'] & 0xe040 == 0
    expected = expanded_fences(shader)
    fences = collections.Counter(re.findall(r'^\s*fence (\w+)\b', ir, re.M))
    locked = sum(bool(re.search(r'\block or\s+DWORD PTR \[rsp[^]]*\],0x0\b', x)) for x in lines)
    full_fences = counts['mfence'] + locked
    coroutine_resumes = [x['name'] for x in event['functions'] if 'resume' in x['name']]
    row_observed = (dict(fences) == {'acq_rel': expected} and full_fences == 0 and
        bool(coroutine_resumes) and counts['vaddps'] > 0 and counts['vmulps'] > 0)
    observed = observed and row_observed
    rows[kind] = {'complete_module': verification['actual_complete_module'],
        'functions': event['functions'], 'total_called_function_bytes': sum(x['bytes'] for x in event['functions']),
        'expanded_original_shader_barriers': expected, 'ordinary_IR_fences_by_order': dict(fences),
        'machine_full_fences': full_fences, 'coroutine_resume_functions': coroutine_resumes,
        'packed_add': counts['vaddps'], 'packed_multiply': counts['vmulps'],
        'callback_state_restoration': verification['callback_restoration_observed'],
        'structural_premise_observed': row_observed}
    if kind in prior['rows']:
        rows[kind]['previous_same_module_seq_cst_fences'] = prior['rows'][kind]['ordinary_IR_fences']
        rows[kind]['previous_same_module_full_fences'] = prior['rows'][kind]['machine_full_fences']
        assert prior['rows'][kind]['complete_module']['sha256'] == record['sha256']
out = {'structural_premise_observed': observed, 'rows': rows,
    'exact_prior_modules_and_cooperative_workgroups_unchanged': True,
    'scope': 'Actual complete ordinary driver code from4B guard dispatch only. Acquire-release fences and coroutine resumes are retained; zero full hardware fences is a lowering observation, not a speedup, exclusive-cost measurement or numerical/scientific/runtime/release acceptance.'}
p = W / 'ordinary-code-review.json'; assert not p.exists()
p.write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({'structural_premise_observed': observed,
    'rows': {k: {x: v[x] for x in ['ordinary_IR_fences_by_order', 'machine_full_fences',
        'coroutine_resume_functions', 'total_called_function_bytes']} for k, v in rows.items()}}))
assert observed, 'No C8/ABBA: acquire-release/intact-rendezvous code premise failed.'
