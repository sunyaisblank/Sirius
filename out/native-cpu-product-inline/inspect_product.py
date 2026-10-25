"""Confirm both real product-helper call boundaries disappeared; no speed claim."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
doc = lambda p: json.loads(p.read_text())
old = {row['stage']: row['object'] for row in doc(R / 'out/native-cpu-zero-term/ispc-zero-compile-result.json')['objects']}
new = doc(W / 'ispc-product-compile-result.json')
owner = doc(W / 'compile/ispc-product/owner.json')
assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert new['all_objects_compiled'] and new['new_compilations'] == 3
assert doc(W / 'object-code-inspection-product.json')['objects']

def inspect(entry):
    p = R / entry['path']
    assert p.stat().st_size == entry['bytes'] and hashlib.file_digest(p.open('rb'), 'sha256').hexdigest() == entry['sha256']
    asm = subprocess.check_output(['objdump', '-dr', '--no-show-raw-insn', str(p)], text=True)
    calls = [line.strip() for line in asm.splitlines() if re.search(r'\b(?:callq?|jmpq?)\b.*RTAccumulate(?:Rounded)?Product_0|R_X86_64_\w+\s+RTAccumulate(?:Rounded)?Product_0', line)]
    sizes = [int(parts[1], 16) for line in subprocess.check_output(['nm', '-S', '--defined-only', str(p)], text=True).splitlines()
             if len(parts := line.split()) == 4 and parts[2] in ('t', 'T')]
    return {'object': entry, 'helper_calls': calls, 'code_bytes': sum(sizes)}

result = {'passed': False, 'stages': {}, 'scope': 'Actual object shape/ISA/FP inspection only; code size is not a speed or runtime acceptance measurement.'}
for row in new['objects']:
    slug = row['stage']; a, b = inspect(old[slug]), inspect(row['object'])
    assert a['helper_calls'] and not b['helper_calls']
    result['stages'][slug] = {'baseline': a, 'candidate': b, 'both_product_helper_ABI_calls_removed': True}
result['passed'] = True
with (W / 'product-code-inspection.json').open('x') as stream:
    json.dump(result, stream, indent=2); stream.write('\n')
print(json.dumps({'passed': True, 'stages': 3, 'runtime_executed': False}))
