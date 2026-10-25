"""Finish accounting for the existing single-compile/single-link diagnostic."""
from pathlib import Path
import json
import re
import subprocess
import sys
from elf_proof_complete import pin, rename_proof, associate

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
assert mode in ('prepare', 'read')
prepared = doc(W / 'prepared.json')
plan = doc(W / 'link-inputs.json')
for phase in ('compile', 'link'):
    owner = doc(W / phase / 'owner.json')
    assert owner['status'] == 'terminal' and not owner['passed'] and owner['child_exit'] == 1
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    assert all(row['absent'] for row in owner['observed_births'])
assert 'AssertionError: /usr/lib/gcc/x86_64-linux-gnu/14/include/stdbool.h' in (W / 'compile/stderr.log').read_text()
assert (W / 'link/stderr.log').read_text().rstrip().endswith('KeyError: 6')
for name, expected in plan['source_inputs'].items():
    assert pin(R / name) == expected, name
baseline, driver, executable = [R / prepared[key] for key in ('baseline', 'driver', 'executable')]
assert executable.is_file() and identity(driver) == plan['driver_object']
if mode == 'prepare':
    assert not args and not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    inputs = dict(plan['source_inputs'])
    for path in (W / 'link-inputs.json', W / 'link/owner.json', W / 'link/stderr.log',
                 W / 'elf_proof_complete.py', W / 'read_build.py', W / 'prepare_run_complete.py', executable):
        inputs[str(path.relative_to(R))] = pin(path)
    proof_inputs = {'source_inputs': inputs, 'executable': identity(executable),
                    'scope': 'Read-only completion after one successful CPP compilation and one successful link. Prior owned phases failed only in dependency accounting and executable dynamic-symbol parsing, before raw workload execution. Preserve both failed owners; no compile/link retries and no claim of their success as whole phases.'}
    write(W / 'proof-inputs.json', proof_inputs)
    actions = doc(W / 'actions.json')
    assert 'proof' not in actions and 'run' not in actions
    actions['proof'] = {**actions['link'],
        'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        'argv': ['/usr/bin/python3', '-B', str(W / 'read_build.py'), 'read', '{folder}'],
        'frozen_inputs': sorted(set(inputs) | {str((W / 'proof-inputs.json').relative_to(R))}),
        'scope': proof_inputs['scope']}
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
else:
    folder, = args
    assert Path(folder) == W / 'proof'
    proof_inputs = doc(W / 'proof-inputs.json')
    for name, expected in proof_inputs['source_inputs'].items():
        assert pin(R / name) == expected, name
    proof = associate(executable, [{'path': str(baseline), **pin(baseline)},
                                  {**prepared['candidate_transport_object'], 'path': str(R / prepared['candidate_transport_object']['path'])}])
    proof['executable'] = identity(executable)
    proof['objects'][0]['object'] = identity(baseline)
    proof['objects'][1]['object'] = prepared['candidate_transport_object']
    write(W / 'link-proof.json', proof)
    assert rename_proof(R / prepared['original_baseline_object']['path'], baseline,
                        prepared['rename_mapping']) == doc(W / 'rename-proof.json')
    resolution = subprocess.check_output(['/usr/bin/ldd', '-r', str(executable)], text=True, stderr=subprocess.STDOUT)
    (W / 'proof/runtime-resolution.txt').write_text(resolution)
    assert 'undefined symbol' not in resolution and 'not found' not in resolution
    libraries = {str(Path(p).resolve(strict=True)) for p in re.findall(r'(/[^\s()]+)', resolution) if Path(p).is_file()}
    for name, expected in proof_inputs['source_inputs'].items():
        assert pin(R / name) == expected, name
    original = doc(W / 'compile-result.json')
    result = {**original, 'passed': True, 'prepared': pin(W / 'prepared.json'),
              'commands': [*original['commands'], {'argv': plan['link_argv'], 'cwd': str(R), 'exit': 0,
                    'scope': 'Executed once before preserved link-owner reader failure; existing ELF now verified.'}],
              'diagnostic_CPP_compilations': 1, 'diagnostic_links': 1,
              'additional_CPP_compilations_in_link_continuation': 0,
              'additional_links_in_proof_completion': 0,
              'driver_object': identity(driver), 'executable': identity(executable),
              'link_proof': identity(W / 'link-proof.json'), 'rename_proof': identity(W / 'rename-proof.json'),
              'runtime_libraries': {p: pin(Path(p)) for p in sorted(libraries)},
              'actual_compiled_dependencies': plan['actual_compiled_dependencies'],
              'source_inputs_unchanged': True, 'compile_continuation': plan,
              'proof_completion': proof_inputs,
              'prior_failed_owners': [identity(W / phase / 'owner.json') for phase in ('compile', 'link')]}
    write(W / 'compile-completion.json', result)
    print(json.dumps({'read_only_proof_passed': True, 'CPP_compilations': 1, 'links': 1, 'raw_executions': 0}))
