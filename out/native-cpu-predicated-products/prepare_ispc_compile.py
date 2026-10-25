"""Seal one exact-recipe compile only after AA's preserved frame rejection."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'out/native-cpu-static-registers'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


failure = R / 'out/native-cpu-static-registers-cohort/cold-frame'
owner = doc(failure / 'frame/owner.json')
assert owner['status'] == 'terminal' and not owner['passed']
assert owner['stop_reason'] == 'diagnostic_timeout'
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
assert doc(failure / 'result.json')['attempted'] and not doc(failure / 'result.json')['accepted']
manifest = R / 'attestations/qualification' / owner['source']['revision'][:8] / 'native-cpu-static-registers-frame-timeout/MANIFEST.json'
assert not doc(manifest)['accepted']
assert not os.environ.get('ISPC_ARGS')
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not (W / 'prepared-ispc-products-inputs.json').exists()
assert not (R / 'bin/linux-gcc/native-cpu-predicated-products/objects').exists()
subprocess.run(['/usr/bin/python3', '-B', str(W / 'derive_ispc.py')], check=True)
derived = doc(W / 'ispc-derivation.json')
old = doc(B / 'ispc-registers-compile-result.json')
assert derived['instructions'] == 9380
assert old['all_objects_compiled'] and old['new_compilations'] == 3
inputs = dict(derived['source_inputs'])
for path in (failure / 'frame/owner.json', failure / 'result.json', manifest):
    inputs[str(path.relative_to(R))] = pin(path)
flags = ['-O2', '--target=avx2-i32x8', '--cpu=core-avx2', '--arch=x86-64', '--addressing=32', '--PIC', '--opt=disable-fma', '--math-lib=system']
objects = []
for row in old['objects']:
    slug = row['stage']
    entry = derived['stages'][slug]
    assert row['argv_executed'] and row['exit'] == 0 and row['argv'][2:-3] == flags
    assert row['argv'][-3:-1] == ['--emit-obj', '-o']
    assert row['argv'][1] == str(R / entry['basis_path'])
    argv = list(row['argv'])
    destination = 'bin/linux-gcc/native-cpu-predicated-products/objects/' + slug + '.o'
    argv[1], argv[-1] = str(R / entry['derived_path']), str(R / destination)
    inverse = list(argv)
    inverse[1], inverse[-1] = row['argv'][1], row['argv'][-1]
    assert inverse == row['argv']
    inputs[entry['derived_path']] = entry['derived']
    objects.append({'stage': slug, 'output': destination, 'argv': argv,
                    'basis_AA_argv': row['argv'], 'argv_inverse_exact': True})
compiler = Path(objects[0]['argv'][0])
assert pin(compiler)['sha256'] == '1813a6b83442dfcc23a1a82a641f22452fedbff7b429877e2f31957957e0168e'
inputs[str(compiler.relative_to(R))] = pin(compiler)
resolution = subprocess.check_output(['/usr/bin/ldd', '-r', str(compiler)], text=True, stderr=subprocess.STDOUT)
assert 'undefined symbol' not in resolution and 'not found' not in resolution
for name in re.findall(r'(/[^\s()]+)', resolution):
    path = Path(name)
    if path.is_file():
        path = path.resolve(strict=True)
        inputs[str(path)] = pin(path)
for name in ('derive_ispc.py', 'ispc-derivation.json', 'prepare_ispc_compile.py', 'build_ispc.py',
             'inspect_objects.py', 'inspect_products.py', 'compile/run_owned.py', 'compile/owned_guard.py'):
    path = W / name
    inputs[str(path.relative_to(R))] = pin(path)
for path in (B / 'ispc-registers-compile-result.json', B / 'compile/ispc-registers/owner.json',
             R / 'out/native-cpu-simd-rows/inspect_objects.py',
             R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
             R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py',
             Path('/usr/bin/nm'), Path('/usr/bin/objdump')):
    inputs[str(path.relative_to(R)) if path.is_relative_to(R) else str(path)] = pin(path)
fp_source = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8/src/gallium/drivers/llvmpipe/lp_state_cs.c'
assert pin(fp_source) == {'bytes': 109640, 'sha256': '34b1f093c09543c7f1c48873081d014849088b3839316dc6b3671d2b6add9cd1'}
inputs[str(fp_source.relative_to(R))] = pin(fp_source)
for name, expected in inputs.items():
    assert pin(R / name) == expected, name
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
scope = 'Compile complete AB raw3 once after preserved AA rejection. Only two product accumulators replace varying early input/zero returns with exact per-lane result selection; original complete nonzero cores and all9380ordered tuples retained. Invalid precedes zero; zero products preserve allfive original sum words. Extra unselected-lane pure FP calculations require the exact pinned masked/restored worker scope and subsequent ISA/linked-FP inspection. Same compiler/flags, full AA metadata/scratch ABI, eligibility/numerical/refusal/publication unchanged. Fresh component/C8 and original frame required; Z broader evidence stays predecessor evidence. No runtime/science/release claim.'
prepared = {'revision': revision, 'inputs': inputs, 'objects': objects,
            'fixed_flags': flags, 'compiler_runtime_resolution': resolution, 'scope': scope}
(W / 'prepared-ispc-products-inputs.json').write_text(json.dumps(prepared, indent=2) + '\n')
base = doc(B / 'compile/actions.json')['ispc-registers']
action = {**base, 'revision': revision, 'argv': ['/usr/bin/python3', '-B', str(W / 'build_ispc.py')],
          'frozen_inputs': sorted(set(inputs) | {str((W / 'prepared-ispc-products-inputs.json').relative_to(R))}),
          'scope': scope}
(W / 'compile/actions.json').write_text(json.dumps({'ispc-products': action}, indent=2) + '\n')
print(json.dumps({'new_compilations_planned': 3, 'executed': False}))
