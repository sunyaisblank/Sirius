"""Complete the successful driver compilation at link, without compiling again."""
from pathlib import Path
import json
import re
import shlex
import subprocess
import sys
from elf_proof import pin, rename_proof, associate

R = Path.cwd()
W = Path(__file__).resolve().parent


def doc(path):
    return json.loads(path.read_text())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def identity(path):
    return {'path': str(path.relative_to(R)), **pin(path)}


mode, *args = sys.argv[1:]
assert mode in ('prepare', 'link')
prepared = doc(W / 'prepared.json')
failed = doc(W / 'compile-result.json')
owner = doc(W / 'compile/owner.json')
assert owner['status'] == 'terminal' and not owner['passed'] and owner['child_exit'] == 1
assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
assert not owner['cleanup_errors'] and not owner['remaining_owned_processes']
assert all(row['absent'] for row in owner['observed_births'])
assert not failed['passed'] and len(failed['commands']) == 2
assert [row['argv'] for row in failed['commands']] == [prepared['rename_argv'], prepared['compile_argv']]
assert all(row['exit'] == 0 for row in failed['commands'])
for name, expected in prepared['source_inputs'].items():
    assert pin(R / name) == {key: expected[key] for key in ('bytes', 'sha256')}, name
baseline, driver, executable = [R / prepared[key] for key in ('baseline', 'driver', 'executable')]
assert not executable.exists() if mode == 'prepare' else True
if mode == 'prepare':
    assert not args and not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    dependencies, added = {}, {}
    text = (W / 'compile/driver.d').read_text().replace('\\\n', ' ')
    for name in shlex.split(text.split(':', 1)[1]):
        path = Path(name)
        path = (path if path.is_absolute() else R / path).resolve(strict=True)
        key = str(path.relative_to(R)) if path.is_relative_to(R) else str(path)
        dependencies[key] = pin(path)
        expected = prepared['source_inputs'].get(key, prepared['source_inputs'].get(str(path)))
        if expected is None:
            added[key] = {**pin(path), 'mtime_ns': path.stat().st_mtime_ns}
        else:
            assert dependencies[key] == {k: expected[k] for k in ('bytes', 'sha256')}
    assert set(added) == {'/usr/lib/gcc/x86_64-linux-gnu/14/include/stdbool.h'}
    assert added[next(iter(added))]['sha256'] == '7b9c438a5b2f082a560a72da8d752da51062cc6f84d1c4acec45ff0a68e8d037'
    assert added[next(iter(added))]['mtime_ns'] < (W / 'compile/inputs-before.json').stat().st_mtime_ns
    inputs = {**prepared['source_inputs'], **dependencies}
    for path in (W / 'prepared.json', W / 'compile-result.json', W / 'compile/owner.json',
                 W / 'compile/driver.d', W / 'compile/stderr.log', W / 'rename-proof.json',
                 W / 'complete_build.py', W / 'prepare_run.py', W / 'run_task.py', baseline, driver):
        inputs[str(path.relative_to(R))] = pin(path)
    plan = {'source_inputs': inputs, 'actual_compiled_dependencies': dependencies,
            'added_after_CPP_compilation_before_link_and_execution': added,
            'driver_object': identity(driver), 'original_compile_result': identity(W / 'compile-result.json'),
            'original_compile_owner': identity(W / 'compile/owner.json'),
            'link_argv': prepared['link_argv'], 'scope': 'Original rename and single CPP compile succeeded; dependency accounting then stopped on one previously unlisted GCC stdbool.h. Preserve that failed owner/receipt. Seal every actual dependency and the successful CPP object now; execute the original unrun link once. No recompile or raw native execution, and no claim that the added header was frozen before the original CPP compile.'}
    write(W / 'link-inputs.json', plan)
    actions = doc(W / 'actions.json')
    assert 'link' not in actions and 'run' not in actions
    actions['link'] = {**actions['compile'],
        'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        'argv': ['/usr/bin/python3', '-B', str(W / 'complete_build.py'), 'link', '{folder}'],
        'frozen_inputs': sorted(set(inputs) | {str((W / 'link-inputs.json').relative_to(R))}),
        'scope': plan['scope']}
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps({'prepared': 'link', 'new_CPP_compilations': 0, 'new_links_planned': 1}))
else:
    folder, = args
    assert Path(folder) == W / 'link'
    plan = doc(W / 'link-inputs.json')
    for name, expected in plan['source_inputs'].items():
        assert pin(R / name) == expected, name
    assert not executable.exists() and identity(driver) == plan['driver_object']
    subprocess.run(plan['link_argv'], cwd=R, check=True)
    proof = associate(executable, [{'path': str(baseline), **pin(baseline)},
                                  {**prepared['candidate_transport_object'], 'path': str(R / prepared['candidate_transport_object']['path'])}])
    proof['executable'] = identity(executable)
    proof['objects'][0]['object'] = identity(baseline)
    proof['objects'][1]['object'] = prepared['candidate_transport_object']
    write(W / 'link-proof.json', proof)
    assert rename_proof(R / prepared['original_baseline_object']['path'], baseline,
                        prepared['rename_mapping']) == doc(W / 'rename-proof.json')
    resolution = subprocess.check_output(['/usr/bin/ldd', '-r', str(executable)], text=True, stderr=subprocess.STDOUT)
    (W / 'link/runtime-resolution.txt').write_text(resolution)
    assert 'undefined symbol' not in resolution and 'not found' not in resolution
    libraries = {str(Path(p).resolve(strict=True)) for p in re.findall(r'(/[^\s()]+)', resolution) if Path(p).is_file()}
    for name, expected in plan['source_inputs'].items():
        assert pin(R / name) == expected, name
    result = {**failed, 'passed': True, 'prepared': pin(W / 'prepared.json'),
              'commands': [*failed['commands'], {'argv': plan['link_argv'], 'cwd': str(R), 'exit': 0}],
              'diagnostic_CPP_compilations': 1, 'diagnostic_links': 1,
              'additional_CPP_compilations_in_link_continuation': 0,
              'driver_object': plan['driver_object'], 'executable': identity(executable),
              'link_proof': identity(W / 'link-proof.json'), 'rename_proof': identity(W / 'rename-proof.json'),
              'runtime_libraries': {p: pin(Path(p)) for p in sorted(libraries)},
              'actual_compiled_dependencies': plan['actual_compiled_dependencies'],
              'source_inputs_unchanged': True, 'compile_continuation': plan}
    write(W / 'compile-completion.json', result)
    print(json.dumps({'link_passed': True, 'new_CPP_compilations': 0, 'total_CPP_compilations': 1}))
