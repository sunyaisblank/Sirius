"""Seal one exact-recipe compile with only the actual required input closure."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'out/native-cpu-product-inline'

def doc(path):
    return json.loads(path.read_text())

def pin(path):
    return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()}

assert not os.environ.get('ISPC_ARGS')
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not (W / 'prepared-ispc-predicated-inputs.json').exists()
assert not (R / 'bin/linux-gcc/native-cpu-predicated-arithmetic/objects').exists()
subprocess.run(['/usr/bin/python3', '-B', str(W / 'derive_ispc.py')], check=True)
derived = doc(W / 'ispc-derivation.json')
old = doc(B / 'ispc-product-compile-result.json')
assert derived['instructions'] == 9380 and old['all_objects_compiled'] and old['new_compilations'] == 3
inputs = dict(derived['source_inputs'])
objects = []
flags = ['-O2', '--target=avx2-i32x8', '--cpu=core-avx2', '--arch=x86-64', '--addressing=32', '--PIC', '--opt=disable-fma', '--math-lib=system']
for row in old['objects']:
    slug = row['stage']; entry = derived['stages'][slug]
    assert row['argv_executed'] and row['exit'] == 0 and row['argv'][2:-3] == flags
    assert row['argv'][-3:-1] == ['--emit-obj', '-o']
    assert row['argv'][1] == str(R / entry['predicated_basis_path'])
    argv = list(row['argv'])
    destination = 'bin/linux-gcc/native-cpu-predicated-arithmetic/objects/' + slug + '.o'
    argv[1], argv[-1] = str(R / entry['derived_path']), str(R / destination)
    inverse = list(argv); inverse[1], inverse[-1] = row['argv'][1], row['argv'][-1]
    assert inverse == row['argv']
    inputs[entry['derived_path']] = entry['derived']
    objects.append({'stage': slug, 'output': destination, 'argv': argv, 'original_W_argv': row['argv'], 'argv_inverse_exact': True})
compiler = Path(objects[0]['argv'][0])
assert pin(compiler)['sha256'] == '1813a6b83442dfcc23a1a82a641f22452fedbff7b429877e2f31957957e0168e'
inputs[str(compiler.relative_to(R))] = pin(compiler)
resolution = subprocess.check_output(['/usr/bin/ldd', '-r', str(compiler)], text=True, stderr=subprocess.STDOUT)
assert 'undefined symbol' not in resolution and 'not found' not in resolution
for name in re.findall(r'(/[^\s()]+)', resolution):
    path = Path(name)
    if path.is_file():
        path = path.resolve(strict=True); inputs[str(path)] = pin(path)
for name in ('derive_ispc.py', 'ispc-derivation.json', 'prepare_ispc_compile.py', 'build_ispc.py',
             'inspect_objects.py', 'inspect_predication.py', 'compile/run_owned.py', 'compile/owned_guard.py'):
    path = W / name; inputs[str(path.relative_to(R))] = pin(path)
for path in (B / 'ispc-product-compile-result.json', B / 'compile/ispc-product/owner.json',
             B / 'inspect_objects.py', R / 'out/native-cpu-simd-rows/inspect_objects.py',
             R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
             R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py',
             Path('/usr/bin/nm'), Path('/usr/bin/objdump')):
    inputs[str(path.relative_to(R)) if path.is_relative_to(R) else str(path)] = pin(path)
for name, expected in inputs.items():
    assert pin(R / name) == expected, name
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
scope = 'Compile three complete X ISPC stages once using exact W compiler/flags; only pure inner arithmetic lane selection/control changes. Ordered live arithmetic expressions and all literal instructions retained. No host/scalar/adapter rebuild, count/precision/guard/FP/governor changes or runtime/scientific claim. Historical receipts are provenance; retired R objects are not current compiler inputs.'
prepared = {'revision': revision, 'inputs': inputs, 'objects': objects, 'fixed_flags': flags,
            'compiler_runtime_resolution': resolution, 'scope': scope}
(W / 'prepared-ispc-predicated-inputs.json').write_text(json.dumps(prepared, indent=2) + '\n')
base = doc(B / 'compile/actions.json')['ispc-product']
action = {**base, 'revision': revision, 'argv': ['/usr/bin/python3', '-B', str(W / 'build_ispc.py')],
          'frozen_inputs': sorted(set(inputs) | {str((W / 'prepared-ispc-predicated-inputs.json').relative_to(R))}), 'scope': scope}
(W / 'compile/actions.json').write_text(json.dumps({'ispc-predicated': action}, indent=2) + '\n')
print(json.dumps({'new_compilations_planned': 3, 'executed': False}))
