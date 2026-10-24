"""Seal one serial ISPC compile of the complete three-stage row translation."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import sys

root = Path.cwd()
work = Path(__file__).resolve().parent
destination = root / 'bin/linux-gcc/native-cpu-simd-rows/objects-uniform'
basis = root / 'out/native-cpu-static-program'


def pin(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def read(path):
    return json.loads(path.read_text())


assert not os.environ.get('ISPC_ARGS', '')
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert not any(p.exists() for p in (destination, work / 'prepared-ispc-uniform-inputs.json',
    work / 'ispc-uniform-compile-result.json'))
assert pin(basis / 'batch-compile-result.json')['sha256'] == 'f8434acad86353ab1a345cd35a9f4e6ef338805a2d87974ebacb49d2d9082be2'
scalar = read(basis / 'batch-compile-result.json')
assert scalar['all_prepared_source_bytes_unchanged'] and len(scalar['objects']) == 3
inputs = dict(scalar['source_inputs'])
for item in scalar['objects']:
    assert item['exit'] == 0
    value = item['object']
    inputs[value['path']] = {k: value[k] for k in ('bytes', 'sha256')}
for name, expected in inputs.items():
    assert pin(root / name) == expected, name

acquisition = read(work / 'compiler-acquisition.json')
environment = read(work / 'compiler-environment.json')
assert acquisition['release_tag'] == 'v1.32.0'
assert acquisition['asset']['digest'] == 'sha256:faf266bda2b702d68976ffd26bbd36c1ddba43424cc3c44f25b092cd1aa1f64a'
compiler_name = environment['compiler']
assert compiler_name == 'bin/linux-gcc/native-cpu-simd-rows/compiler/bin/ispc'
assert acquisition['extracted_members'][compiler_name]['sha256'] == '1813a6b83442dfcc23a1a82a641f22452fedbff7b429877e2f31957957e0168e'
assert '1.32.0' in environment['version_output']
for name, value in acquisition['extracted_members'].items():
    expected = {k: value[k] for k in ('bytes', 'sha256')}
    assert pin(root / name) == expected, name
    inputs[name] = expected
# Reuse the installer-recorded runtime resolution; never select another tool.
runtime_paths = set(re.findall(r'/[^\s()]+', environment['runtime_resolution']))
runtime_paths |= {'/etc/ld.so.cache', '/usr/bin/python3', str(Path(sys.executable).resolve())}
for name in runtime_paths:
    assert Path(name).is_file(), name
    inputs[name] = pin(Path(name))

# Default derivation mode only compares existing bytes; it cannot overwrite.
subprocess.run([sys.executable, '-B', str(work / 'derive_ispc.py')], check=True)
derived = read(work / 'ispc-uniform-derivation.json')
assert derived['instructions'] == 9380 and set(derived['stages']) == {'transport', 'endpoint', 'dopri_phase'}
flags = ['-O2', '--target=sse2-i32x8', '--cpu=x86-64', '--arch=x86-64',
         '--addressing=64', '--PIC', '--opt=disable-fma', '--math-lib=system']
objects = []
for slug in ('transport', 'endpoint', 'dopri_phase'):
    row = derived['stages'][slug]
    assert pin(root / row['basis_path']) == row['basis']
    assert pin(root / row['derived_path']) == row['derived']
    assert row['logical_invocations_per_active_row'] == 64 and not row['physical_math_rederived']
    assert scalar['source_inputs'][row['basis_path']] == row['basis']
    text = (root / row['derived_path']).read_text()
    assert not re.search(r'\wLOAD_WORD\b|\.Load\(|\.GetDimensions\(', text)
    assert 'uniform uint32_t * uniform input' in text
    assert 'uniform uint32_t * uniform output' in text
    assert 'varying KernelContext_0 context;' in text
    assert 'if (programIndex < groups) ComputeOwnedRow_0(id, invocation, &context);' in text
    inputs[row['derived_path']] = row['derived']
    output = destination / (slug + '.o')
    objects.append({'stage': slug, 'output': str(output.relative_to(root)),
                    'argv': [str(root / compiler_name), str(root / row['derived_path']),
                             *flags, '--emit-obj', '-o', str(output)]})
for name in ('compiler-acquisition.json', 'compiler-environment.json', 'derive_ispc.py',
             'ispc-uniform-derivation.json', 'prepare_ispc_compile.py', 'build_ispc.py',
             'compile/run_owned.py', 'compile/owned_guard.py'):
    path = work / name
    inputs[str(path.relative_to(root))] = pin(path)
for path in (basis / 'batch-compile-result.json', basis / 'header-preparation.json',
             root / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
             root / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py'):
    inputs[str(path.relative_to(root))] = pin(path)
admission = root / 'bin/linux-gcc/native-cpu-dispatch-admission/admission.o'
assert pin(admission) == {'bytes': 194704, 'sha256': 'b8587c2f6c834b72b8378c4634a92837e3644e4fd494ce283897711559c4c616'}
inputs[str(admission.relative_to(root))] = pin(admission)
prepared = {'revision': revision, 'inputs': inputs, 'objects': objects, 'fixed_flags': flags,
            'scope': 'Three complete per-row ISPC objects with uniform canonical operands/private scratch indices; preserved 9380 literal calls/all5canonicalprograms, original generic helpers and all arithmetic. Original three D scalar objects and F admission object remain unchanged. No provider, host, shader execution or qualification.'}
(work / 'prepared-ispc-uniform-inputs.json').write_text(json.dumps(prepared, indent=2) + '\n')
old = read(work / 'generation/actions.json')['c-source']
action = {**old, 'revision': revision, 'kind': 'compile',
          'argv': ['/usr/bin/python3', '-B', str(work / 'build_ispc.py')],
          'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
          'frozen_inputs': sorted(set(inputs) | {str((work / 'prepared-ispc-uniform-inputs.json').relative_to(root))}),
          'scope': prepared['scope']}
actions = read(work / 'compile/actions.json')
assert 'ispc-uniform' not in actions
actions['ispc-uniform'] = action
(work / 'compile/actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'prepared_objects': 3, 'preserved_scalar_objects': 3, 'preserved_admission_objects': 1,
                  'compiled': False, 'revision': revision}))
