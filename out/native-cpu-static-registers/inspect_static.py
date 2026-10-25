"""Require the two shared private-scratch helpers in each actual new object."""
from pathlib import Path
import hashlib
import json
import subprocess

R = Path.cwd(); W = Path(__file__).resolve().parent
doc = lambda p: json.loads(p.read_text())
compiled = doc(W / 'ispc-registers-compile-result.json')
owner = doc(W / 'compile/ispc-registers/owner.json')
assert compiled['all_objects_compiled'] and compiled['new_compilations'] == 3
assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert doc(W / 'object-code-inspection-registers.json')['objects']
rows = {}
for item in compiled['objects']:
    entry = item['object']; path = R / entry['path']
    assert path.stat().st_size == entry['bytes'] and hashlib.file_digest(path.open('rb'), 'sha256').hexdigest() == entry['sha256']
    functions = {}
    for line in subprocess.check_output(['nm', '-S', '--defined-only', str(path)], text=True).splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[2] in ('t', 'T'):
            functions[parts[3]] = int(parts[1], 16)
    helpers = {name: size for name, size in functions.items() if name.startswith(('RTAddStore_0___', 'RTMultiplyStore_0___'))}
    assert any(name.startswith('RTAddStore_0___') for name in helpers)
    assert any(name.startswith('RTMultiplyStore_0___') for name in helpers)
    rows[item['stage']] = {'object': entry, 'shared_store_helpers': helpers, 'total_code_bytes': sum(functions.values())}
result = {'passed': True, 'stages': rows,
    'scope': 'Object shape plus separate ISA/FP inspection only; code size is not a performance measure. Requires whole-program opcode/module joins, fresh full38/C8 and the original frame; Z broader results remain explicitly predecessor evidence. No runtime, equivalence, frame or scientific/release claim.'}
with (W / 'static-register-code-inspection.json').open('x') as stream:
    json.dump(result, stream, indent=2); stream.write('\n')
print(json.dumps({'passed': True, 'stages': 3, 'runtime_executed': False}))
