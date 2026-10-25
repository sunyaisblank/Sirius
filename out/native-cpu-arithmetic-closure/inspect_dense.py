"""Require both new C ABI references and exact merged-adapter function joins."""
from pathlib import Path
import hashlib
import json
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
doc = lambda path: json.loads(path.read_text())
compiled = doc(W / 'ispc-dense-compile-result.json')
owner = doc(W / 'compile/ispc-dense/owner.json')
assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert compiled['all_objects_compiled'] and compiled['new_compilations'] == 3
assert compiled['new_CPP_compilations'] == 0 and compiled['new_relocatable_links'] == 1
assert doc(W / 'object-code-inspection-dense.json')['objects']
assert compiled['adapter_input_code_association']['passed']
assert all(row['all_retained_function_opcodes_match'] for row in compiled['adapter_input_code_association']['objects'])
for row in [item['object'] for item in compiled['objects']] + [compiled['merged_adapter']]:
    path = R / row['path']
    assert path.stat().st_size == row['bytes'] and hashlib.file_digest(path.open('rb'), 'sha256').hexdigest() == row['sha256']
for item in compiled['objects']:
    path = R / item['object']['path']
    undefined = subprocess.check_output(['/usr/bin/nm', '-u', str(path)], text=True)
    relocations = subprocess.check_output(['/usr/bin/objdump', '-r', str(path)], text=True)
    for name in ('sirius_native_dense_add8', 'sirius_native_dense_multiply8', 'sirius_native_dense_divide8', 'sirius_native_dense_sqrt8'):
        assert name in undefined and name in relocations
result = {'passed': True, 'raw_objects': [row['object'] for row in compiled['objects']],
          'merged_adapter': compiled['merged_adapter'],
          'adapter_original_and_CPP_code_association': compiled['adapter_input_code_association'],
          'both_dense_ABI_references_in_all3': True,
          'scope': 'Source/ISA and relocatable object joins only. Every original adapter and checked C++ function opcode matches its combined object. Both new call relocations are present in each raw stage. Actual linked helper targets and full38/C8/frame still required; no current renderer performance/scientific acceptance.'}
with (W / 'dense-code-inspection.json').open('x') as stream:
    json.dump(result, stream, indent=2)
    stream.write('\n')
print(json.dumps({'passed': True, 'stages': 3, 'original_adapter_and_CPP_code_joined': True}))
