"""One bounded optional compile and unchanged original arithmetic compatibility."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-exact-product-residual'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def write(path, value):
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2) + '\n')


def owner(phase):
    value = doc(W / phase / 'owner.json')
    assert value['status'] == 'terminal' and value['passed']
    assert value['source_unchanged'] and value['whole_inputs_unchanged']
    assert not value['remaining_owned_processes'] and not value['cleanup_errors']
    return value


phase, = sys.argv[1:]
proof = doc(W / 'derivation.json')
for name, expected in proof['source_inputs'].items():
    assert pin(R / name) == expected, name
for name, expected in proof['derived'].items():
    assert pin(W / name) == expected, name
for row in proof['reused_objects']:
    assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}
old = R / 'out/native-cpu-value-fact-stores'
prepared, built = doc(old / 'prepared-cpp.json'), doc(old / 'cpp-build-result.json')
assert built['passed']
for name, expected in prepared['inputs'].items():
    assert pin(R / name) == expected, name
driver = next(row['object'] for row in built['compilations'] if row['name'] == 'driver')
assert pin(R / driver['path']) == {k: driver[k] for k in ('bytes', 'sha256')}

if phase in ('build', 'link', 'run'):
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    inputs = set(prepared['inputs']) | set(proof['source_inputs'])
    inputs |= {row['path'] for row in proof['reused_objects']} | {driver['path']}
    inputs |= {str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
    if phase == 'build':
        assert not B.exists()
        argv = ['/usr/bin/python3', '-B', str(W / 'verify.py'), 'compile']
        scope = 'One separately guarded optional-v4 compile; strict original flags and five-word/eight-lane ABI. Only bounded product32 residual changes; original-v3 and baseline dispatcher reused exactly. No workload or performance/adoption claim.'
        actions = {}
    else:
        actions = doc(W / 'actions.json')
        owner('build')
        inputs |= {str((B / 'optional.o').relative_to(R)), str((W / 'compile-result.json').relative_to(R))}
        assert pin(B / 'optional.o') == doc(W / 'compile-result.json')['optional_object']
        if phase == 'link':
            argv = [prepared['link'][0], str(R / driver['path']),
                    str(R / 'bin/linux-gcc/native-cpu-v4-arithmetic/dispatch.o'),
                    str(R / 'bin/linux-gcc/native-cpu-v4-arithmetic/v3.o'), str(B / 'optional.o'),
                    '-o', str(B / 'abi-check')]
            scope = 'Link exact already-compiled original555660-call tagged fixture to unchanged dispatch/v3 and new optional object; no fixture compilation, provider or acceptance transfer.'
        else:
            owner('link')
            inputs.add(str((B / 'abi-check').relative_to(R)))
            argv = [str(B / 'abi-check')]
            scope = 'One unchanged original555660-call arithmetic compatibility fixture; bank/tag/mask/alias/refusal/FP-control checks unchanged. Runtime fast-branch coverage and benefit remain unqualified. No provider/full38/C8/frame/scientific/release claim.'
    assert phase not in actions and not (W / phase).exists()
    actions[phase] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='diagnostic', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=argv, frozen_inputs=sorted(inputs), scope=scope)
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=phase, executed=False)))
elif phase == 'compile':
    B.mkdir(exist_ok=False)
    args = proof['compile_argv']
    subprocess.run(args, check=True)
    write(W / 'compile-result.json', dict(argv=args, exit=0, optional_object=pin(B / 'optional.o'),
        driver_recompiled=False, original_v3_recompiled=False, baseline_dispatcher_recompiled=False))
elif phase == 'read':
    value = owner('run')
    lines = (W / 'run/stdout.log').read_text().splitlines()
    result, = [json.loads(line) for line in lines if line.startswith('{"passed"')]
    assert result['passed'] and result['calls'] == 555660
    write(W / 'compatibility.json', dict(passed=True, comparison=result,
        optional_object=pin(B / 'optional.o'), executable=pin(B / 'abi-check'),
        scope='Original finite scalar compatibility only; no positive fast-branch coverage, benefit, provider, frame or scientific/release qualification.'))
    print(json.dumps(result))
else:
    raise ValueError(phase)
