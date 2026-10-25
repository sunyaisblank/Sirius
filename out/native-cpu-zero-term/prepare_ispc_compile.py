"""Seal the actual R compiler recipe for three new full-stage compilations."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'out/native-cpu-accumulate-inline'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


assert not os.environ.get('ISPC_ARGS')
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not any((W / name).exists() for name in ('prepared-ispc-zero-inputs.json', 'ispc-zero-compile-result.json', 'compile/actions.json'))
assert not (R / 'bin/linux-gcc/native-cpu-zero-term/objects').exists()
subprocess.run(['/usr/bin/python3', '-B', str(W / 'derive_ispc.py')], check=True)
derived, old = doc(W / 'ispc-derivation.json'), doc(B / 'ispc-inline-compile-result.json')
assert derived['instructions'] == 9380 and old['all_objects_compiled'] and old['new_compilations'] == 3
inputs = {**old['source_inputs'], **derived['source_inputs']}
objects = []
flags = ['-O2', '--target=avx2-i32x8', '--cpu=core-avx2', '--arch=x86-64', '--addressing=32', '--PIC', '--opt=disable-fma', '--math-lib=system']
for row in old['objects']:
    slug = row['stage']; entry = derived['stages'][slug]
    assert row['argv_executed'] and row['exit'] == 0 and row['argv'][2:-3] == flags
    assert row['argv'][-3:-1] == ['--emit-obj', '-o']
    assert row['argv'][1] == str(R / entry['zero_term_basis_path'])
    argv = list(row['argv'])
    destination = 'bin/linux-gcc/native-cpu-zero-term/objects/' + slug + '.o'
    argv[1], argv[-1] = str(R / entry['derived_path']), str(R / destination)
    inverse = list(argv); inverse[1], inverse[-1] = row['argv'][1], row['argv'][-1]
    assert inverse == row['argv']
    inputs[entry['derived_path']] = entry['derived']
    objects.append({'stage': slug, 'output': destination, 'argv': argv,
                    'original_R_argv': row['argv'], 'argv_inverse_exact': True})
compiler = Path(objects[0]['argv'][0])
assert pin(compiler)['sha256'] == '1813a6b83442dfcc23a1a82a641f22452fedbff7b429877e2f31957957e0168e'
resolution = subprocess.check_output(['/usr/bin/ldd', '-r', str(compiler)], text=True, stderr=subprocess.STDOUT)
assert 'undefined symbol' not in resolution and 'not found' not in resolution
for name in re.findall(r'(/[^\s()]+)', resolution):
    path = Path(name)
    if path.is_file():
        path = path.resolve(strict=True)
        current = pin(path)
        assert str(path) not in inputs or inputs[str(path)] == current
        inputs[str(path)] = current
for name in ('derive_ispc.py', 'ispc-derivation.json', 'prepare_ispc_compile.py', 'build_ispc.py',
             'inspect_objects.py', 'inspect_zero.py', 'compile/run_owned.py', 'compile/owned_guard.py'):
    path = W / name; inputs[str(path.relative_to(R))] = pin(path)
for path in (B / 'ispc-inline-compile-result.json', B / 'compile/ispc-inline/owner.json',
             Path('/usr/bin/nm'), Path('/usr/bin/objdump')):
    inputs[str(path.relative_to(R)) if path.is_relative_to(R) else str(path)] = pin(path)
for name, expected in inputs.items():
    assert pin(R / name) == expected, name
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
scope = 'Compile all3 complete S ISPC stages once with exact R compiler/flags, only reviewed signed-zero arithmetic arm added. No host/adapter/scalar/provider recompilation or runtime/scientific qualification; all9380calls/counts/precision/guards/FP/governors retained.'
prepared = {'revision': revision, 'inputs': inputs, 'objects': objects, 'fixed_flags': flags,
            'compiler_runtime_resolution': resolution, 'scope': scope}
(W / 'prepared-ispc-zero-inputs.json').write_text(json.dumps(prepared, indent=2) + '\n')
base = doc(B / 'compile/actions.json')['ispc-inline']
action = {**base, 'revision': revision, 'argv': ['/usr/bin/python3', '-B', str(W / 'build_ispc.py')],
          'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
          'frozen_inputs': sorted(set(inputs) | {str((W / 'prepared-ispc-zero-inputs.json').relative_to(R))}),
          'scope': scope}
(W / 'compile/actions.json').write_text(json.dumps({'ispc-zero': action}, indent=2) + '\n')
print(json.dumps({'new_compilations_planned': 3, 'executed': False}))
