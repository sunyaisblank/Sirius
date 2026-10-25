"""Compare actual new code to R; retain code size/stack proxies, no speed claim."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()}


def code(identity):
    path = R / identity['path']
    assert pin(path) == {key: identity[key] for key in ('bytes', 'sha256')}
    sizes = {}
    for line in subprocess.check_output(['/usr/bin/nm', '-S', '--defined-only', str(path)], text=True).splitlines():
        words = line.split()
        if len(words) == 4 and words[2] in ('t', 'T'):
            sizes[words[3]] = int(words[1], 16)
    asm = subprocess.check_output(['/usr/bin/objdump', '-dr', '--no-show-raw-insn', str(path)], text=True)
    terms = [line.strip() for line in asm.splitlines() if re.search(r'\b(?:callq?|jmpq?)\b.*RTAccumulateTerm_0|R_X86_64_\w+\s+RTAccumulateTerm_0', line)]
    selected = {}
    for body in re.split(r'\n(?=[0-9a-f]+ <)', asm):
        match = re.match(r'[0-9a-f]+ <([^>]+)>:', body)
        if match and match[1].startswith(('RTAdd_0', 'RTAccumulateProduct_0', 'RTAccumulateRoundedProduct_0')):
            instructions = [line for line in body.splitlines() if re.match(r'^\s*[0-9a-f]+:\s+\S', line)]
            selected[match[1]] = {'code_bytes': sizes[match[1]], 'instructions': len(instructions),
                'stack_reference_instructions': sum('%rsp' in line or '%rbp' in line for line in instructions)}
    return {'object': identity, 'object_code_bytes': sum(sizes.values()), 'term_calls': terms,
            'selected_callers': selected}


record = doc(W / 'ispc-zero-compile-result.json')
owner = doc(W / 'compile/ispc-zero/owner.json')
isa_path = W / 'object-code-inspection-zero.json'
isa = doc(isa_path)
assert record['all_objects_compiled'] and record['all_prepared_source_bytes_unchanged'] and record['new_compilations'] == 3
assert owner['passed'] and owner['status'] == 'terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
old = {row['stage']: row['object'] for row in doc(R / 'out/native-cpu-accumulate-inline/ispc-inline-compile-result.json')['objects']}
result = {'passed': False, 'revision': record['revision'], 'compile_result': pin(W / 'ispc-zero-compile-result.json'),
          'ISA_FP_inspection': pin(isa_path), 'stages': {},
          'scope': 'Actual three-object source/code/ISA inspection only. Signed-zero semantics follow the reviewed argument; full native compatibility/original runtime gates still required. Code size and stack-reference proxies are not benefit/liveness measurements.'}
for row in record['objects']:
    slug = row['stage']
    assert row['argv_executed'] and row['exit'] == 0 and isa['objects'][slug]['object'] == row['object']
    baseline, candidate = code(old[slug]), code(row['object'])
    assert not baseline['term_calls'] and not candidate['term_calls']
    assert baseline['selected_callers'] != candidate['selected_callers'], 'No material new caller code shape'
    result['stages'][slug] = {'baseline': baseline, 'candidate': candidate,
        'all_RTAccumulateTerm_ABI_calls_removed': True,
        'object_code_byte_delta': candidate['object_code_bytes'] - baseline['object_code_bytes']}
result['passed'] = True
with (W / 'zero-code-inspection.json').open('x') as stream:
    json.dump(result, stream, indent=2); stream.write('\n')
print(json.dumps({'passed': True, 'stages': 3, 'runtime_executed': False}))
