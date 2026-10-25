"""Seal one new link and one comparison using the unchanged genuine R fixture driver."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
OLD = R / 'out/native-cpu-accumulate-compatibility'
B = R / 'bin/linux-gcc/native-cpu-local-identities-compatibility'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


def named(row):
    path = R / row['path']
    assert pin(path) == {key: row[key] for key in ('bytes', 'sha256')}
    return path


mode, = sys.argv[1:]
assert mode in ('link', 'run')
assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
compiled_path = R / 'out/native-cpu-local-identities/ispc-dense-compile-result.json'
compiled = doc(compiled_path)
assert compiled['all_objects_compiled'] and compiled['all_prepared_source_bytes_unchanged']
assert compiled['new_compilations'] == 3 and compiled['new_raw_compilations'] == 3 and all(row['argv_executed'] and row['exit'] == 0 for row in compiled['objects'])
contract_path = R / 'out/native-cpu-local-identities/provider/source-contract.json'
objects = doc(contract_path)['native_objects']
assert len(objects) == 8 and objects[4:7] == [row['object'] for row in compiled['objects']]
old_prepared = doc(OLD / 'prepared.json')
old_build = doc(OLD / 'compile-completion.json')
old_result = doc(OLD / 'result.json')
assert old_build['passed'] and old_build['diagnostic_CPP_compilations'] == 1 and old_build['diagnostic_links'] == 1
assert old_build['prepared'] == pin(OLD / 'prepared.json')
assert named(old_result['compile_result']) == OLD / 'compile-completion.json'
assert old_result['passed'] and old_result['comparison']['compared_output_words'] == 187960
for phase in ('proof', 'run'):
    owner = doc(OLD / phase / 'owner.json')
    assert owner['passed'] and owner['status'] == 'terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
for key in ('driver_object', 'renamed_baseline_object', 'rename_proof'):
    named(old_build[key])
inputs = {**compiled['source_inputs'], **doc(contract_path)['source_inputs'], **old_build['actual_compiled_dependencies']}
for row in objects + [old_build['driver_object'], old_build['renamed_baseline_object'], old_prepared['original_baseline_object']]:
    named(row); inputs[row['path']] = {key: row[key] for key in ('bytes', 'sha256')}
for path in (compiled_path, contract_path, OLD / 'prepared.json', OLD / 'compile-completion.json',
             OLD / 'result.json', OLD / 'proof/owner.json', OLD / 'run/owner.json',
             R / 'out/native-cpu-field-compatibility/elf_proof.py', OLD / 'elf_proof_complete.py', OLD / 'task.py', OLD / 'rename-proof.json',
             R / 'out/native-cpu-accumulate-inline/coverage-plan.json', *W.glob('*.py')):
    inputs[str(path.relative_to(R))] = pin(path)
tool_receipt = R / 'out/native-cpu-host-representation/probe/preparation.json'
tools = doc(tool_receipt)['tool_dependencies']
assert any(Path(name).name == 'collect2' for name in tools)
for name, expected in tools.items():
    assert pin(Path(name)) == expected, name
    inputs[name] = expected
inputs[str(tool_receipt.relative_to(R))] = pin(tool_receipt)
for name in ('/usr/bin/x86_64-linux-gnu-g++-14', '/usr/bin/ld', '/usr/bin/ldd'):
    path = Path(name).resolve(strict=True); inputs[str(path)] = pin(path)
if mode == 'link':
    assert not B.exists() and not (W / 'prepared.json').exists() and not (W / 'actions.json').exists()
    argv = list(old_prepared['link_argv'])
    for old, new in zip(old_prepared['objects'], objects):
        original, replacement = str(R / old['path']), str(R / new['path'])
        assert argv.count(original) == 1
        argv[argv.index(original)] = replacement
    assert argv.count(str(R / old_build['driver_object']['path'])) == 1
    assert argv.count(str(R / old_build['renamed_baseline_object']['path'])) == 1
    argv[argv.index('-o') + 1] = str(B / 'compatibility')
    assert '-Wl,--no-relax' in argv
    assert '-c' not in argv and not any(word.endswith(('.cpp', '.cc', '.cxx')) for word in argv)
    for name, expected in inputs.items():
        assert pin(R / name) == expected, name
    prepared = {'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        'source_inputs': inputs, 'objects': objects, 'candidate_transport_object': objects[4],
        'original_baseline_object': old_prepared['original_baseline_object'],
        'baseline': old_build['renamed_baseline_object']['path'], 'driver': old_build['driver_object']['path'],
        'executable': str((B / 'compatibility').relative_to(R)), 'link_argv': argv,
        'original_R_link_argv': old_prepared['link_argv'], 'reused_R_build': {'path': str((OLD / 'compile-completion.json').relative_to(R)), **pin(OLD / 'compile-completion.json')},
        'scope': 'Reuse exact genuine R fixture driver and renamed Q baseline object; both already compiled/verified. New AN raw3 with original integer identity stores before unchanged AM CPP calls plus unchanged U scalar3/pureF and combined original Iadapter/checked dense C++, one new link, zero CPP/old-math recompilations or objcopy operations. Original38 inputs/capacity40/187960words and scoped FP/ISA/admission checks unchanged. Historical R accounting failures stay historical.'}
    B.mkdir(parents=True)
    with (W / 'prepared.json').open('x') as stream:
        json.dump(prepared, stream, indent=2); stream.write('\n')
    actions = {}
else:
    owner = doc(W / 'link/owner.json')
    assert owner['passed'] and owner['status'] == 'terminal' and owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    build = doc(W / 'link-result.json'); assert build['passed'] and build['new_CPP_compilations'] == 0
    prepared = doc(W / 'prepared.json')
    inputs.update(prepared['source_inputs']); inputs.update(build['runtime_libraries'])
    for key in ('driver_object', 'renamed_baseline_object', 'rename_proof', 'link_proof', 'executable'):
        row = build[key]; named(row); inputs[row['path']] = {k: row[k] for k in ('bytes', 'sha256')}
    for path in (W / 'link-result.json', W / 'link/owner.json'):
        inputs[str(path.relative_to(R))] = pin(path)
    actions = doc(W / 'actions.json'); assert 'run' not in actions and not (W / 'run').exists()
inputs[str((W / 'prepared.json').relative_to(R))] = pin(W / 'prepared.json')
action = {'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    'kind': 'diagnostic', 'expected_exit': 0, 'software_provider': False,
    'limits': {'seconds': 300, 'rss_kib': 4194304, 'reserve_kib': 2097152},
    'frozen_inputs': sorted(inputs),
    'argv': ['/usr/bin/python3', '-B', str(W / ('link.py' if mode == 'link' else 'run_task.py')), mode, '{folder}'],
    'scope': prepared['scope'] if mode == 'link' else 'One fresh AN Transport compatibility workload, all38 unchanged cases/all187960words. No repeated study, tolerance, timing or original-gate substitution.'}
actions[mode] = action
(W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
print(json.dumps({'prepared': mode, 'new_CPP_compilations': 0, 'executed': False}))
