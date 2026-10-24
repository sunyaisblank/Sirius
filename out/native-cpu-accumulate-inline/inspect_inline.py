"""Verify actual inlining/code shape after the owned compile; no shader execution."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()}


record = doc(work / 'ispc-inline-compile-result.json')
owner = doc(work / 'compile/ispc-inline/owner.json')
assert record['all_objects_compiled'] and record['all_prepared_source_bytes_unchanged']
assert record['new_compilations'] == 3 and record['reused_compiled_objects'] == 0
assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
isa_path = work / 'object-code-inspection-inline.json'
isa = doc(isa_path)
assert set(isa['objects']) == {'transport', 'endpoint', 'dopri_phase'}
old = doc(root / 'out/native-cpu-simd-endpoint-regions/ispc-uniform-compile-result.json')
baseline = {row['stage']: row['object'] for row in old['objects']}
selected = ('RTAccumulateProduct_0', 'RTAccumulateRoundedProduct_0', 'RTAdd_0')


def decode(identity):
    path = root / identity['path']
    assert pin(path) == {key: identity[key] for key in ('bytes', 'sha256')}
    argv = ['/usr/bin/objdump', '-dr', '--no-show-raw-insn', str(path)]
    asm = subprocess.check_output(argv, text=True, timeout=60)
    sizes = {}
    for line in subprocess.check_output(['/usr/bin/nm', '-S', '--defined-only', str(path)], text=True).splitlines():
        words = line.split()
        if len(words) == 4 and words[2] in ('t', 'T'):
            sizes[words[3]] = int(words[1], 16)
    bodies = {}
    for part in re.split(r'\n(?=[0-9a-f]+ <)', asm):
        match = re.match(r'[0-9a-f]+ <([^>]+)>:', part)
        if match:
            bodies[match[1]] = part
    calls = [line.strip() for line in asm.splitlines() if
             re.search(r'\b(?:callq?|jmpq?)\b.*RTAccumulateTerm_0|R_X86_64_\w+\s+RTAccumulateTerm_0', line)]
    metrics = {}
    for prefix in selected:
        rows = []
        for name, body in bodies.items():
            if not name.startswith(prefix):
                continue
            ins = [line.strip() for line in body.splitlines() if re.match(r'^\s*[0-9a-f]+:\s+\S', line)]
            rows.append({'symbol': name, 'code_bytes': sizes[name], 'instructions': len(ins),
                'stack_reference_instructions': sum('%rsp' in line or '%rbp' in line for line in ins),
                'stack_vector_stores': sum(bool(re.search(r'\bv\w+\s+%[xy]mm\d+,.*\(%r(?:sp|bp)\)', line)) for line in ins),
                'stack_vector_loads': sum(bool(re.search(r'\bv\w+\s+.*\(%r(?:sp|bp)\).*%[xy]mm\d+', line)) for line in ins),
                'stack_pointer_operations': [line for line in ins if
                    re.search(r'\b(?:sub|add|and|mov)\b.*,%rsp|\b(?:push|pop)\b', line)],
                'referenced_vector_registers': sorted(set(re.findall(r'%[xy]mm\d+', body))),
                'term_call_sites': [line for line in ins if re.search(r'\b(?:callq?|jmpq?)\b.*RTAccumulateTerm_0', line)]})
        metrics[prefix] = rows  # An omitted outer helper is recorded, not invented.
    return {'object': identity, 'objdump_argv': argv, 'object_code_bytes': sum(sizes.values()),
            'term_symbols': [name for name in sizes if name.startswith('RTAccumulateTerm_0')],
            'all_term_call_or_jump_or_relocation_sites': calls, 'callers': metrics}


assert not (work / 'inline-code-inspection.json').exists()
result = {'passed': False, 'revision': record['revision'], 'stages': {},
    'compile_result': pin(work / 'ispc-inline-compile-result.json'), 'owner': pin(work / 'compile/ispc-inline/owner.json'),
    'ISA_FP_inspection': pin(isa_path), 'scope': 'Actual object source/code inspection only. Stack traffic and referenced register sets are proxies, not liveness/spill-cost or speedup measurements. All original runtime/scientific/release gates and new Transport compatibility remain outstanding.'}
for row in record['objects']:
    slug = row['stage']
    assert row['argv_executed'] and row['exit'] == 0 and isa['objects'][slug]['object'] == row['object']
    previous, candidate = decode(baseline[slug]), decode(row['object'])
    assert previous['all_term_call_or_jump_or_relocation_sites']
    assert not candidate['all_term_call_or_jump_or_relocation_sites'], (slug, candidate)
    for caller in selected:
        assert previous['callers'][caller]
        assert all(not item['term_call_sites'] for item in candidate['callers'][caller])
    result['stages'][slug] = {'baseline': previous, 'candidate': candidate,
        'object_code_byte_delta': candidate['object_code_bytes'] - previous['object_code_bytes'],
        'all_RTAccumulateTerm_ABI_calls_removed': True}
result['passed'] = True
(work / 'inline-code-inspection.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'passed': True, 'stages': 3, 'remaining_term_calls': 0, 'runtime_executed': False}))
