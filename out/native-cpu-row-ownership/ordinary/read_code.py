"""Read complete captured ordinary IR and actual called code; never dispatch."""
from pathlib import Path
import collections, hashlib, json, re, subprocess
R = Path.cwd(); W = Path(__file__).resolve().parent
compiler = json.loads((W.parent / 'compiler-result.json').read_text())
emission = json.loads((W.parent / 'emission-result.json').read_text()); assert emission['passed']
rows = {}
for label, kind in [('transport', 'Transport'), ('endpoint', 'Endpoint'), ('dopri_phase', 'DopriPhase')]:
    f = W / label; v = json.loads((f / 'capture-verification.json').read_text())
    b = json.loads((f / 'bindings.json').read_text())
    e, = [x for x in b['jit_events'] if x.get('file') == v['actual_entry']['file']]
    assert v['actual_complete_module']['sha256'] == compiler['modules'][kind]['sha256']
    assert v['actual_entry']['actual_mxcsr'] & 0xe040 == 0 and v['callback_restoration_observed']
    section, = [x for x in e['live_sections_after_pointer_lookup'] if x['name'] == '.ltext']
    raw = (f / section['path']).read_bytes(); assert hashlib.sha256(raw).hexdigest() == section['sha256']
    regions = [(x['address'], x['address'] + x['bytes']) for x in e['functions']]
    dump = subprocess.check_output(['/usr/bin/objdump', '-D', '-b', 'binary', '-m', 'i386:x86-64',
        '-M', 'intel', '--adjust-vma=' + str(section['address']), str(f / section['path'])], text=True)
    lines = []; ops = []
    for line in dump.splitlines():
        m = re.match(r'\s*([0-9a-f]+):\s+(?:[0-9a-f]{2}\s+)+\s*([a-z][a-z0-9]*)\s*(.*)', line)
        if m and any(a <= int(m[1], 16) < z for a, z in regions): lines.append(line); ops.append(m[2])
    (W / (label + '-called-functions.asm')).write_text('\n'.join(lines) + '\n')
    ir = (f / v['ordinary_ir_export']['path']).read_text()
    assert hashlib.sha256(ir.encode()).hexdigest() == v['ordinary_ir_export']['sha256']
    mutable_globals = [x for x in ir.splitlines() if x.startswith('@') and re.search(r'\bglobal\b', x)]
    assert not mutable_globals
    registers = compiler['modules'][kind]['registers']
    scratch = re.findall(r'^\s*(%\S+) = alloca \[' + str(registers * 5) + r' x <8 x i32>\], align 32$', ir, re.M)
    assert len(scratch) == 1
    assert not re.findall(r'^\s*fence ', ir, re.M)
    assert len(e['functions']) == 1 and e['functions'][0]['name'] == 'cs_variant'
    counts = collections.Counter(ops)
    locked_full = sum(bool(re.search(r'\block or\s+DWORD PTR \[rsp[^]]*\],0x0\b', x)) for x in lines)
    assert counts['mfence'] + locked_full == 0
    assert counts['vaddps'] + counts['vaddss'] > 0 and counts['vmulps'] + counts['vmulss'] > 0
    rows[kind] = {'module': v['actual_complete_module'], 'functions': e['functions'],
        'total_called_function_bytes': sum(x['bytes'] for x in e['functions']),
        'ordinary_mutable_process_globals': 0, 'scratch_is_function_stack_alloca': scratch[0],
        'scratch_stack_vector_bytes': registers * 5 * 8 * 4,
        'actual_SIMD_lanes': 8, 'ordinary_IR_fences': 0, 'machine_full_fences': 0,
        'coroutine_resume_functions': 0, 'packed_add': counts['vaddps'], 'scalar_add': counts['vaddss'],
        'packed_multiply': counts['vmulps'], 'scalar_multiply': counts['vmulss'],
        'callback_state_restoration_observed': True}
out = {'structural_premise_observed': True, 'rows': rows,
    'actual_per_invocation_stack_storage_no_mutable_process_globals': True,
    'same_LocalSize64_and_group_row_correspondence_and_tables': True,
    'cost_limit': 'Ordinary code remains8-wide SIMD with one active logical invocation; stack scratch is8 copies. Code is not smaller. No speed, memory-win, numerical, benefit, frame or release conclusion follows.',
    'scope': 'Complete actual changed modules, ordinary optimized LLVM and called machine code. Four-byte guard dispatch only.'}
p = W / 'ordinary-code-review.json'; assert not p.exists(); p.write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({'passed': True, 'called_function_bytes': {k: v['total_called_function_bytes'] for k, v in rows.items()},
    'actual_mutable_process_globals_fences_coroutines': 0, 'SIMD_lanes': 8,
    'scratch_stack_bytes': {k: v['scratch_stack_vector_bytes'] for k, v in rows.items()}}))
