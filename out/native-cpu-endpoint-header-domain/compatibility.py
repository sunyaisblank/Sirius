"""Build and execute the one focused real Endpoint routing comparison."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-endpoint-header-domain'


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def write(path, value):
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2) + '\n')


def owner(name):
    value = doc(W / name / 'owner.json')
    assert value['passed'] and value['status'] == 'terminal'
    assert value['source_unchanged'] and value['whole_inputs_unchanged']
    assert not value['remaining_owned_processes'] and not value['cleanup_errors']
    return value


phase, = sys.argv[1:]
build = doc(W / 'build-result.json')
owner('build')
assert build['new_compilations'] == 2
assert doc(W / 'code-proof.json')['passed']
for row in build['objects']:
    assert row['exit'] == 0
    assert pin(R / row['object']['path']) == {k: row['object'][k] for k in ('bytes', 'sha256')}
retained = doc(R / 'out/native-cpu-v4-provider/provider/source-contract.json')
for row in [retained['native_objects'][5]] + retained['adapter_input_objects'][1:]:
    assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}, row['path']

if phase in ('prepare-driver', 'prepare-run'):
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    inputs = set(doc(W / 'actions.json')['build']['frozen_inputs'])
    inputs |= {str(p.relative_to(R)) for p in W.rglob('*') if p.is_file()
               and p.name != 'actions.json'}
    inputs |= {str(p.relative_to(R)) for p in (R / 'tests/support/retained_transport').glob('*.h')}
    objects = [row['object']['path'] for row in build['objects']]
    objects += ['bin/linux-gcc/native-cpu-local-identities/objects/endpoint.o']
    objects += ['bin/linux-gcc/native-cpu-v4-arithmetic/' + name + '.o'
                for name in ('dispatch', 'v3', 'optional')]
    inputs |= set(objects)
    if phase == 'prepare-driver':
        assert (W / 'driver.cpp').is_file()
        args = list(build['objects'][1]['argv'])
        args = args[:args.index('-c')]
        # The reference consumes the original generic runtime program header,
        # before later static-program/scalar-zero source specializations.
        generic = str(R / 'bin/linux-gcc/native-cpu-slang-bridge/typed-batch')
        at = args.index(generic)
        assert args[at - 1] == '-isystem'
        del args[at - 1:at + 1]
        args[1:1] = ['-isystem', generic]
        args += ['-I', str(R), '-c', str(W / 'driver.cpp'), '-o', str(B / 'driver.o')]
        link = [args[0], str(B / 'driver.o')] + [str(R / n) for n in objects]
        link += ['-Wl,--no-relax', '-o', str(B / 'compatibility')]
        write(W / 'driver-plan.json', dict(compile=args, link=link, objects=objects))
        inputs.add(str((W / 'driver-plan.json').relative_to(R)))
        name, argv = 'driver-build', ['/usr/bin/python3', '-B', str(W / 'compatibility.py'), 'driver-build']
        scope = 'One focused driver compile and link of unchanged generated scalar reference body, real raw Endpoint and exact verified AO math; original packet ABI and FP flags. No provider or timing qualification.'
    else:
        owner('driver-build')
        assert doc(W / 'driver-build-result.json')['passed']
        inputs |= {str((B / n).relative_to(R)) for n in ('driver.o', 'compatibility')}
        name, argv = 'compatibility', [str(B / 'compatibility')]
        scope = 'One focused mixed23-row Endpoint complete-output comparison, canonical and malformed specialized register headers, real scalar and raw SIMD plus scalar tails, unchanged inputs/inactive words/canaries/FP controls and all-other-program-word refusal. Supporting scalar compatibility, not independent mathematical or provider/scientific/release qualification.'
    assert all((R / path).is_file() for path in inputs)
    actions = doc(W / 'actions.json')
    assert name not in actions and not (W / name).exists()
    actions[name] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='diagnostic', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=argv, frozen_inputs=sorted(inputs), scope=scope)
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=name, executed=False, frozen_inputs=len(inputs))))
elif phase == 'driver-build':
    plan = doc(W / 'driver-plan.json')
    for key in ('compile', 'link'):
        subprocess.run(plan[key], check=True)
    write(W / 'driver-build-result.json', dict(passed=True, plan=pin(W / 'driver-plan.json'),
        driver=pin(B / 'driver.o'), executable=pin(B / 'compatibility')))
elif phase == 'read':
    value = owner('compatibility')
    results = [json.loads(line) for line in (W / 'compatibility/stdout.log').read_text().splitlines()
               if line.startswith('{')]
    result, = results
    assert result['passed']
    for key, expected in dict(selection_rows=23, capacity=32, header_values=6,
            layouts=2, reference_rows=138, scalar_packets=360, raw_packets=24,
            other_program_word_refusals=21142, boundary_refusals=21).items():
        assert result[key] == expected, (key, result[key])
    for key in ('whole3829_words_and_inactive_capacity_and_canaries', 'input_immutable',
                'MXCSR_controls_unchanged', 'caller_MXCSR_state_restored'):
        assert result[key]
    write(W / 'compatibility.json', dict(passed=True, comparison=result,
        objects=[row['object'] for row in build['objects']],
        inputs_unchanged=value['whole_inputs_unchanged'], source_unchanged=value['source_unchanged'],
        executable=dict(path=str((B / 'compatibility').relative_to(R)), **pin(B / 'compatibility')),
        reference_header=dict(path='bin/linux-gcc/native-cpu-slang-bridge/typed-batch/endpoint.native.hpp',
            **pin(R / 'bin/linux-gcc/native-cpu-slang-bridge/typed-batch/endpoint.native.hpp')),
        scope='Focused routing and whole-row scalar compatibility only. The original independent Control1 and all provider/cohort/scientific/release gates remain required.'))
    print(json.dumps(result))
else:
    raise ValueError(phase)
