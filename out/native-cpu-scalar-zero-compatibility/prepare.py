"""Reuse the original full38 driver/baseline; compile only the new ABI shim."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
OLD = R / 'out/native-cpu-accumulate-compatibility'
B = R / 'bin/linux-gcc/native-cpu-scalar-zero-compatibility'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


def add(path):
    inputs[str(path.relative_to(R)) if path.is_relative_to(R) else str(path)] = pin(path)


mode, = sys.argv[1:]
assert mode in ('build', 'run')
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
old = doc(OLD / 'compile-completion.json')
old_prep = doc(OLD / 'prepared.json')
assert old['passed'] and doc(OLD / 'result.json')['passed']
inputs = dict(old['actual_compiled_dependencies'])
objects = []
for kind in ('cpp',):
    receipt = R / ('out/native-cpu-scalar-zero/' + kind + '-compile-result.json')
    compiled = doc(receipt); owner = doc(receipt.parent / 'compile' / kind / 'owner.json')
    assert compiled['all_objects_compiled'] and compiled['new_compilations'] == 3
    assert owner['passed'] and owner['whole_inputs_unchanged'] and owner['source_unchanged']
    assert owner['remaining_owned_processes'] == [] and not owner['cleanup_errors']
    inputs.update(compiled['source_inputs'])
    for item in compiled['objects']:
        row = item['object']; add(R / row['path'])
        assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}
        if item['stage'] == 'transport':
            objects.append(row)
    add(receipt); add(receipt.parent / 'compile' / kind / 'owner.json')
admission = old_prep['objects'][3]; objects.append(admission)
for row in objects + [old['driver_object'], old['renamed_baseline_object'], old_prep['original_baseline_object']]:
    assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}
    add(R / row['path'])
for path in (OLD / 'compile-completion.json', OLD / 'prepared.json', OLD / 'result.json',
             OLD / 'task.py', OLD / 'elf_proof.py', OLD / 'elf_proof_complete.py', OLD / 'rename-proof.json', R / 'out/native-cpu-field-compatibility/elf_proof.py',
             R / 'out/native-cpu-accumulate-inline/coverage-plan.json',
             R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
             R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py', *W.glob('*.py'), W / 'shim.cpp'):
    add(path)
tools = doc(R / 'out/native-cpu-host-representation/probe/preparation.json')['tool_dependencies']
for name, expected in tools.items():
    assert pin(Path(name)) == expected, name
    inputs[name] = expected
if mode == 'build':
    assert not B.exists() and not (W / 'prepared.json').exists() and not (W / 'actions.json').exists()
    argv = list(old_prep['compile_argv'])
    argv[argv.index('-c')+1] = str(W / 'shim.cpp')
    argv[argv.index('-o')+1] = str(B / 'shim.o')
    argv[argv.index('-MF')+1] = str(W / 'build/shim.d')
    link = old_prep['link_argv'][:old_prep['link_argv'].index(str(R / old['driver_object']['path']))]
    link += [str(R / old['driver_object']['path']), str(R / old['renamed_baseline_object']['path'])]
    link += [str(R / row['path']) for row in objects]
    link += [str(B / 'shim.o'), '-o', str(B / 'compatibility'), '-Wl,--no-relax', '-lm']
    prepared = {'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        'source_inputs': inputs, 'objects': objects, 'candidate_transport_object': objects[0],
        'driver': old['driver_object']['path'], 'baseline': old['renamed_baseline_object']['path'],
        'executable': str((B / 'compatibility').relative_to(R)), 'compile_argv': argv, 'link_argv': link,
        'original_baseline_object': old_prep['original_baseline_object'],
        'scope': 'New scalar-zero Transport CPP plus original pure admission, unchanged full38 driver and renamed baseline. One new diagnostic ABI-shim compile and one link; zero driver/baseline/old math recompilations or objcopy. All38 original cases,40 capacity,187960words, input/FP/ISA/inactive checks retained. Finite compatibility only.'}
    for name, expected in inputs.items():
        assert pin(R / name) == expected, name
    B.mkdir(parents=True)
    with (W / 'prepared.json').open('x') as stream:
        json.dump(prepared, stream, indent=2); stream.write('\n')
    actions = {}
else:
    owner = doc(W / 'build/owner.json'); build = doc(W / 'build-result.json')
    assert owner['passed'] and owner['whole_inputs_unchanged'] and not owner['remaining_owned_processes']
    assert not owner['cleanup_errors'] and build['passed']
    prepared = doc(W / 'prepared.json')
    inputs.update(prepared['source_inputs']); inputs.update(build['runtime_libraries'])
    for key in ('driver_object', 'renamed_baseline_object', 'executable', 'link_proof', 'rename_proof'):
        row = build[key]; add(R / row['path'])
        assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}
    add(W / 'build-result.json'); add(W / 'build/owner.json')
    actions = doc(W / 'actions.json'); assert 'run' not in actions
add(W / 'prepared.json')
actions[mode] = {'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'kind': 'diagnostic', 'expected_exit': 0, 'software_provider': False,
    'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'frozen_inputs': sorted(inputs),
    'argv': ['/usr/bin/python3', '-B', str(W / ('build.py' if mode == 'build' else 'run_task.py')), mode, '{folder}'],
    'scope': prepared['scope']}
(W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'prepared': mode, 'new_shim_compilations': 1 if mode == 'build' else 0, 'executed': False}))
