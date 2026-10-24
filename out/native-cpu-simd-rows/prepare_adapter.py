"""Seal the eight-row C ABI adapter using the exact original strict C++ argv."""
from pathlib import Path
import hashlib
import json
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent


def pin(path):
    return {'bytes': path.stat().st_size,
            'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()}


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not (work / 'prepared-adapter.json').exists()
assert not (work / 'adapter-compile-result.json').exists()
record = json.loads((work / 'ispc-uniform-compile-result.json').read_text())
owner = json.loads((work / 'compile/ispc-uniform/owner.json').read_text())
assert record['all_objects_compiled'] and owner['passed'] and owner['whole_inputs_unchanged']
inputs = dict(record['source_inputs'])
for item in record['objects']:
    value = item['object']
    inputs[value['path']] = {k: value[k] for k in ('bytes', 'sha256')}
for name, expected in inputs.items():
    assert pin(root / name) == expected, name
scalar = json.loads((root / 'out/native-cpu-static-program/batch-compile-result.json').read_text())
argv = list(scalar['objects'][0]['argv'])
assert all(flag in argv for flag in ('-std=c++2c', '-O2', '-march=x86-64',
    '-mfpmath=sse', '-fno-fast-math', '-frounding-math', '-ffp-contract=off'))
output = root / 'bin/linux-gcc/native-cpu-simd-rows/objects/rows_adapter.o'
assert not output.exists()
argv[argv.index('-c') + 1] = str(work / 'rows_adapter.cpp')
argv[argv.index('-o') + 1] = str(output)
for name in ('rows_adapter.cpp', 'rows_bridge.h', 'prepare_adapter.py', 'build_adapter.py',
             'ispc-uniform-compile-result.json', 'compile/ispc-uniform/owner.json', 'object-code-inspection-uniform.json'):
    path = work / name
    inputs[str(path.relative_to(root))] = pin(path)
prepared = {'inputs': inputs, 'argv': argv, 'output': str(output.relative_to(root))}
(work / 'prepared-adapter.json').write_text(json.dumps(prepared, indent=2) + '\n')
actions = json.loads((work / 'compile/actions.json').read_text())
assert 'adapter' not in actions
actions['adapter'] = {**actions['ispc-uniform'],
    'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'argv': ['/usr/bin/python3', '-B', str(work / 'build_adapter.py')],
    'frozen_inputs': sorted(set(inputs) | {str((work / 'prepared-adapter.json').relative_to(root))}),
    'scope': 'One C ABI adapter compile only, exact original strict C++ flags. All3 SIMD and all3 original scalar objects plus F pure admission remain unchanged. No provider or shader execution.'}
(work / 'compile/actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'prepared_adapter_objects': 1, 'compiled': False}))
