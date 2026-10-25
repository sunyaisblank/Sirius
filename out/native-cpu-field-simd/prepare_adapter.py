"""Seal one new strict adapter compile after the actual full38 comparison."""
from pathlib import Path
import hashlib
import json
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not (W / 'prepared-adapter.json').exists()
assert not (W / 'adapter-compile-result.json').exists()
compatibility = R / 'out/native-cpu-field-compatibility'
assert doc(compatibility / 'result.json')['passed'] and doc(compatibility / 'run/owner.json')['passed']
assert doc(compatibility / 'result.json')['comparison']['compared_output_words'] == 187960
inputs = dict(doc(W / 'prepared-compile.json')['inputs'])
for kind in ('cpp', 'ispc'):
    result = doc(W / (kind + '-compile-result.json'))
    assert result['all_objects_compiled'] and doc(W / 'compile' / kind / 'owner.json')['passed']
    for row in result['objects']:
        inputs[row['object']['path']] = {k: row['object'][k] for k in ('bytes', 'sha256')}
for path in (W / 'rows_adapter.cpp', W / 'prepare_adapter.py', W / 'build_adapter.py',
             W / 'object-code-inspection.json', compatibility / 'result.json', compatibility / 'run/owner.json',
             R / 'out/native-cpu-simd-packet-policy/rows_bridge.h'):
    inputs[str(path.relative_to(R))] = pin(path)
scalar = doc(R / 'out/native-cpu-static-program/batch-compile-result.json')
argv = list(scalar['objects'][0]['argv'])
output = R / 'bin/linux-gcc/native-cpu-field-simd/objects/rows_adapter.o'
assert not output.exists()
argv[argv.index('-c')+1], argv[argv.index('-o')+1] = str(W / 'rows_adapter.cpp'), str(output)
for name, expected in inputs.items():
    assert pin(R / name) == expected, name
(W / 'prepared-adapter.json').write_text(json.dumps({'inputs': inputs, 'argv': argv,
    'output': str(output.relative_to(R))}, indent=2) + '\n')
actions = doc(W / 'compile/actions.json'); assert 'adapter' not in actions
actions['adapter'] = {**actions['cpp'],
    'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'argv': ['/usr/bin/python3', '-B', str(W / 'build_adapter.py')],
    'frozen_inputs': sorted(set(inputs) | {str((W / 'prepared-adapter.json').relative_to(R))}),
    'scope': 'One field-row adapter compile with original strict scalar recipe; actual OS/ISA guard, pureF admission and D fallback retained. No runtime/benefit/scientific acceptance.'}
(W / 'compile/actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'adapter_compilations_planned': 1, 'executed': False}))
