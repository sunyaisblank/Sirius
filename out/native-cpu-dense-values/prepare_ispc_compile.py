"""Seal one exact-recipe compile only after AA's preserved frame rejection."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'out/native-cpu-dense-arithmetic'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


failure = R / 'out/native-cpu-dense-arithmetic-cohort/cold-frame'
owner = doc(failure / 'frame/owner.json')
assert owner['status'] == 'terminal' and not owner['passed']
assert owner['stop_reason'] == 'diagnostic_timeout'
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
assert doc(failure / 'result.json')['attempted'] and not doc(failure / 'result.json')['accepted']
manifest = R / 'attestations/qualification' / owner['source']['revision'][:8] / 'native-cpu-dense-arithmetic-frame-timeout/MANIFEST.json'
assert not doc(manifest)['accepted']
assert not os.environ.get('ISPC_ARGS')
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not (W / 'prepared-ispc-values-inputs.json').exists()
object_directory = R / 'bin/linux-gcc/native-cpu-dense-values/objects'
assert {p.name for p in object_directory.iterdir()} == {'values_bridge.o'}
bridge_identity = doc(W / 'bridge-build-result.json')['object']
assert pin(R / bridge_identity['path']) == {k: bridge_identity[k] for k in ('bytes', 'sha256')}
subprocess.run(['/usr/bin/python3', '-B', str(W / 'derive_ispc.py')], check=True)
derived = doc(W / 'ispc-derivation.json')
old = doc(B / 'ispc-dense-compile-result.json')
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
    destination = 'bin/linux-gcc/native-cpu-dense-values/objects/' + slug + '.o'
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
             'inspect_objects.py', 'inspect_values.py', 'compile/run_owned.py', 'compile/owned_guard.py'):
    path = W / name
    inputs[str(path.relative_to(R))] = pin(path)
for path in (B / 'ispc-dense-compile-result.json', B / 'compile/ispc-dense/owner.json',
             R / 'out/native-cpu-simd-rows/inspect_objects.py',
             R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
             R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py',
             Path('/usr/bin/nm'), Path('/usr/bin/objdump')):
    inputs[str(path.relative_to(R)) if path.is_relative_to(R) else str(path)] = pin(path)
fp_source = R / 'out/lavapipe-acqrel-ordering/upstream/mesa-25.2.8/src/gallium/drivers/llvmpipe/lp_state_cs.c'
assert pin(fp_source) == {'bytes': 109640, 'sha256': '34b1f093c09543c7f1c48873081d014849088b3839316dc6b3671d2b6add9cd1'}
inputs[str(fp_source.relative_to(R))] = pin(fp_source)
bridge = doc(W / 'bridge-build-result.json')
bridge_owner = doc(W / 'bridge/build/owner.json')
assert bridge['passed'] and bridge['new_CPP_compilations'] == 1 and bridge['arithmetic_compilations'] == 0
assert bridge_owner['passed'] and bridge_owner['source_unchanged'] and bridge_owner['whole_inputs_unchanged']
assert not bridge_owner['remaining_owned_processes'] and not bridge_owner['cleanup_errors']
helper = bridge['object']
previous_adapter = doc(B / 'provider/source-contract.json')['native_objects'][-1]
for row in (helper, previous_adapter):
    assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}
    inputs[row['path']] = {k: row[k] for k in ('bytes', 'sha256')}
for path in (W / 'bridge-build-result.json', W / 'bridge/build/owner.json',
             W / 'bridge-code-inspection.json', B / 'cpp-code-inspection.json',
             B / 'cpp/run/owner.json', B / 'cpp/run/stdout.log', B / 'provider/source-contract.json',
             R / 'out/native-cpu-field-compatibility/elf_proof.py',
             R / 'out/native-cpu-accumulate-compatibility/elf_proof.py'):
    inputs[str(path.relative_to(R))] = pin(path)
assert doc(B / 'cpp/run/stdout.log')['compared_words'] == 1425600
linker = Path('/usr/bin/ld').resolve(strict=True)
inputs[str(linker)] = pin(linker)
for name, expected in doc(W / 'prepared-bridge.json')['inputs'].items():
    assert pin(R / name) == expected
    inputs[name] = expected
merged_adapter = 'bin/linux-gcc/native-cpu-dense-values/objects/rows_adapter.o'
merged_argv = [str(linker), '-r', '-o', str(R / merged_adapter),
               str(R / previous_adapter['path']), str(R / helper['path'])]
for name, expected in inputs.items():
    assert pin(R / name) == expected, name
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
scope = 'Compile AD raw3 once after preserved AC frame rejection. Remaining compound Add/Multiply functions use explicit initialized word buffers and live lane masks through the small checked packing adapter. Existing AC direct-scratch calls remain. No arithmetic recompile or source change: exact prior combined I/checked C++ math object is reused by relocatable link with one new packing object. All9380tuples and other original row/control/refusal/rollback/publication remain; existing strict ISPC compiler/flags and Y FP/scheduling retained. Fresh linked targets/full38/C8/original frame required; Z broader results remain predecessor evidence. No independent math/performance/science/release claim.'
prepared = {'revision': revision, 'inputs': inputs, 'objects': objects,
            'fixed_flags': flags, 'compiler_runtime_resolution': resolution, 'scope': scope,
            'merged_adapter': merged_adapter, 'merged_adapter_argv': merged_argv,
            'adapter_input_objects': [previous_adapter, helper]}
(W / 'prepared-ispc-values-inputs.json').write_text(json.dumps(prepared, indent=2) + '\n')
base = doc(B / 'compile/actions.json')['ispc-dense']
action = {**base, 'revision': revision, 'argv': ['/usr/bin/python3', '-B', str(W / 'build_ispc.py')],
          'frozen_inputs': sorted(set(inputs) | {str((W / 'prepared-ispc-values-inputs.json').relative_to(R))}),
          'scope': scope}
(W / 'compile/actions.json').write_text(json.dumps({'ispc-values': action}, indent=2) + '\n')
print(json.dumps({'new_compilations_planned': 3, 'executed': False}))
