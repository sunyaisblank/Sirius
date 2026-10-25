"""Inspect actual helper code shape; no runtime speed or equivalence claim."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

R = Path.cwd(); W = Path(__file__).resolve().parent
doc = lambda p: json.loads(p.read_text())
old = {row['stage']: row['object'] for row in doc(R / 'out/native-cpu-product-inline/ispc-product-compile-result.json')['objects']}
new = doc(W / 'ispc-predicated-compile-result.json')
owner = doc(W / 'compile/ispc-predicated/owner.json')
assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert new['all_objects_compiled'] and new['new_compilations'] == 3
assert doc(W / 'object-code-inspection-predicated.json')['objects']

def inspect(entry):
    p = R / entry['path']
    assert p.stat().st_size == entry['bytes'] and hashlib.file_digest(p.open('rb'), 'sha256').hexdigest() == entry['sha256']
    sizes = {}
    for line in subprocess.check_output(['nm', '-S', '--defined-only', str(p)], text=True).splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[2] in ('t', 'T'):
            sizes[parts[3]] = int(parts[1], 16)
    asm = subprocess.check_output(['objdump', '-dr', '--no-show-raw-insn', str(p)], text=True)
    calls = [line.strip() for line in asm.splitlines() if re.search(r'\b(?:callq?|jmpq?)\b.*RTAccumulateTerm_0|R_X86_64_\w+\s+RTAccumulateTerm_0', line)]
    helpers = {name: size for name, size in sizes.items() if name.startswith(('RTAdd_0___', 'RTMultiply_0___', 'RPUpAdd_0___', 'RetainedSum32_0___'))}
    assert not calls
    return {'object': entry, 'Term_ABI_calls': 0, 'code_bytes': sum(sizes.values()), 'helper_code_bytes': helpers}

result = {'passed': False, 'stages': {}, 'scope': 'Actual code shape and independent ISA/FP inspection only. All selected original arithmetic expressions and lane result mapping follow source reasoning; full38/C8 and original later gates still required. Static sizes are not runtime frequencies, causal shares or speed evidence.'}
for row in new['objects']:
    slug = row['stage']; before, after = inspect(old[slug]), inspect(row['object'])
    assert before['helper_code_bytes'] != after['helper_code_bytes']
    result['stages'][slug] = {'basis': before, 'candidate': after,
        'object_code_byte_delta': after['code_bytes'] - before['code_bytes']}
result['passed'] = True
with (W / 'predicated-code-inspection.json').open('x') as stream:
    json.dump(result, stream, indent=2); stream.write('\n')
print(json.dumps({'passed': True, 'stages': 3, 'runtime_executed': False}))
