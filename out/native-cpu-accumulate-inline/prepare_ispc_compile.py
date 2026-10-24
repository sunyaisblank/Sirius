"""Seal the exact L compiler recipe for three real qualifier-only compilations."""
from pathlib import Path
import hashlib
import json
import os
import subprocess

root = Path.cwd()
work = Path(__file__).resolve().parent
prior = root / 'out/native-cpu-simd-endpoint-regions'
destination = root / 'bin/linux-gcc/native-cpu-accumulate-inline/objects'


def read(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


assert not os.environ.get('ISPC_ARGS', '')
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert not any(path.exists() for path in (destination, work / 'prepared-ispc-inline-inputs.json',
    work / 'ispc-inline-compile-result.json', work / 'compile/actions.json'))
subprocess.run(['/usr/bin/python3', '-B', str(work / 'derive_ispc.py')], check=True)
derived = read(work / 'ispc-derivation.json')
old = read(prior / 'ispc-uniform-compile-result.json')
assert old['all_objects_compiled'] and old['all_prepared_source_bytes_unchanged']
old_owner = read(prior / 'compile/ispc-uniform/owner.json')
assert old_owner['passed'] and old_owner['source_unchanged'] and old_owner['whole_inputs_unchanged']
assert not old_owner['cleanup_errors'] and not old_owner['remaining_owned_processes']
assert derived['instructions'] == 9380 and set(derived['stages']) == {'transport', 'endpoint', 'dopri_phase'}
inputs = {**old['source_inputs'], **derived['source_inputs']}
acquisition, environment = read(work / 'compiler-acquisition.json'), read(work / 'compiler-environment.json')
assert acquisition == read(prior / 'compiler-acquisition.json')
assert environment == read(prior / 'compiler-environment.json')
assert acquisition['release_tag'] == 'v1.32.0'
compiler = root / environment['compiler']
assert pin(compiler)['sha256'] == '1813a6b83442dfcc23a1a82a641f22452fedbff7b429877e2f31957957e0168e'
q = read(root / 'out/native-cpu-transport-worker-occupancy/provider/source-contract.json')
assert len(q['native_objects']) == 8
for obj in q['native_objects']:
    inputs[obj['path']] = {key: obj[key] for key in ('bytes', 'sha256')}
objects, flags = [], ['-O2', '--target=avx2-i32x8', '--cpu=core-avx2', '--arch=x86-64',
    '--addressing=32', '--PIC', '--opt=disable-fma', '--math-lib=system']
for row in old['objects']:
    slug = row['stage']
    entry = derived['stages'][slug]
    assert row['exit'] == 0 and row['argv'][0] == str(compiler) and row['argv'][2:-3] == flags
    assert row['argv'][-3:-1] == ['--emit-obj', '-o']
    assert row['argv'][1] == str(root / entry['inline_basis_path'])
    actual = row
    if not row['argv_executed']:
        # L truthfully reused K Transport/Dopri; bind the actual earlier compile too.
        previous = root / 'out/native-cpu-simd-structural-uniform/ispc-uniform-compile-result.json'
        actual, = [item for item in read(previous)['objects'] if item['stage'] == slug]
        assert actual['argv_executed'] and actual['exit'] == 0
        assert actual['argv'][0] == row['argv'][0] and actual['argv'][2:-3] == flags
        assert {key: actual['object'][key] for key in ('bytes', 'sha256')} == {
            key: row['object'][key] for key in ('bytes', 'sha256')}
        inputs[str(previous.relative_to(root))] = pin(previous)
    assert entry['logical_invocations_per_active_row'] == 64 and not entry['physical_math_rederived']
    assert entry['one_inline_qualifier_only'] and entry['full_inverse_restores_L_source']
    output = destination / (slug + '.o')
    argv = list(row['argv'])
    argv[1], argv[-1] = str(root / entry['derived_path']), str(output)
    inverse = list(argv); inverse[1], inverse[-1] = row['argv'][1], row['argv'][-1]
    assert inverse == row['argv']
    inputs[entry['derived_path']] = entry['derived']
    inputs[entry['basis_path']] = entry['basis']
    objects.append({'stage': slug, 'output': str(output.relative_to(root)), 'argv': argv,
        'original_L_argv': row['argv'], 'original_L_argv_executed': row['argv_executed'],
        'actual_basis_compile_argv': actual['argv'], 'argv_inverse_exact': True})
assert [item['stage'] for item in objects] == ['transport', 'endpoint', 'dopri_phase']
for name in ('compiler-acquisition.json', 'compiler-environment.json', 'derive_ispc.py',
    'ispc-derivation.json', 'prepare_ispc_compile.py', 'build_ispc.py', 'inspect_objects.py',
    'inspect_inline.py', 'compile/run_owned.py', 'compile/owned_guard.py', 'prepare_proposal.py',
    'source-proof.json', 'accumulate-inline.patch', 'excerpts.asm', 'coverage-plan.json'):
    path = work / name
    inputs[str(path.relative_to(root))] = pin(path)
for path in (prior / 'ispc-uniform-compile-result.json', prior / 'compile/ispc-uniform/owner.json',
    prior / 'compile/actions.json', root / 'out/native-cpu-simd-rows/inspect_objects.py',
    root / 'out/native-cpu-transport-worker-occupancy/provider/source-contract.json',
    root / 'out/native-cpu-runtime-observe/decode_all_threads.py', Path('/usr/bin/objdump'), Path('/usr/bin/nm')):
    inputs[str(path.relative_to(root)) if path.is_relative_to(root) else str(path)] = pin(path)
for name, expected in inputs.items():
    assert pin(root / name) == expected, name
prepared = {'revision': revision, 'inputs': inputs, 'objects': objects, 'fixed_flags': flags,
    'scope': 'Compile all3 complete R ISPC sources once, only RTAccumulateTerm inline qualifiers changed. Exact L/K compiler/tool/flags/9380calls preserved; Q other5 objects/P host/ISA/scheduler/FP/guards unchanged. No provider/host/other-native compilation, shader runtime, benefit or science qualification.'}
(work / 'prepared-ispc-inline-inputs.json').write_text(json.dumps(prepared, indent=2) + '\n')
base = read(prior / 'compile/actions.json')['ispc-uniform']
action = {**base, 'revision': revision, 'argv': ['/usr/bin/python3', '-B', str(work / 'build_ispc.py')],
    'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'frozen_inputs': sorted(set(inputs) | {str((work / 'prepared-ispc-inline-inputs.json').relative_to(root))}),
    'scope': prepared['scope']}
(work / 'compile/actions.json').write_text(json.dumps({'ispc-inline': action}, indent=2) + '\n')
print(json.dumps({'prepared_objects': 3, 'new_native_compilations_planned': 3, 'executed': False}))
