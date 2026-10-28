"""Build and apply the two reviewed, fixed strict-conversion passes once."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import re
import shlex
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-strict-ir-linking'
D = B / 'toolchain/bin'


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,
                    sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


phase, = sys.argv[1:]
outputs = {'remaining-plugin': 'constrain_remaining.so',
           'constrain-transport': 'transport.strict.bc', 'constrain-endpoint': 'endpoint.strict.bc'}
if phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    assert json.loads((W / 'machine-dp-features/owner.json').read_text())['passed']
    for stage in ('transport', 'endpoint'):
        owner = json.loads((W / ('emit-' + stage) / 'owner.json').read_text())
        assert owner['passed'] and owner['status'] == 'terminal'
    config = shlex.split(subprocess.check_output([str(D / 'llvm-config'), '--cxxflags'], text=True))
    commands = {'remaining-plugin': [str(D / 'clang++'), *config, '-O2', '-fPIC', '-shared',
                 str(W / 'constrain_remaining.cpp'), str(B / 'sha256-only.a'),
                 '-Wl,--exclude-libs,ALL', '-o', str(B / outputs['remaining-plugin'])]}
    for stage in ('transport', 'endpoint'):
        commands['constrain-' + stage] = [str(D / 'opt'),
            '-load-pass-plugin=' + str(B / outputs['remaining-plugin']),
            '-passes=sirius-constrain-' + stage + ',verify', '-verify-each',
            str(B / (stage + '.ll')), '-o', str(B / outputs['constrain-' + stage])]
    paths = [*W.glob('*.py'), W / 'constrain_remaining.cpp', W / 'remaining_scope_review.json',
             W / 'scope_contract.json', W / 'toolchain.json', B / 'sha256-only.a',
             R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
             R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py']
    for stage in ('transport', 'endpoint'):
        row = json.loads((W / (stage + '-inventory.json')).read_text())['output']
        assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}
        paths += [R / row['path'], W / (stage + '-inventory.json'),
                  R / ('bin/linux-gcc/native-cpu-inline-private-bank/ispc-source/' + stage + '.ispc')]
    for name in ('clang++', 'opt', 'llvm-dis', 'llvm-config'):
        row = json.loads((W / 'toolchain.json').read_text())['tools'][name]
        assert pin(R / row['resolved']) == {k: row[k] for k in ('bytes', 'sha256')}
        paths.append(R / row['resolved'])
    inputs = {str(p.relative_to(R)): pin(p) for p in paths}
    write(W / 'remaining-plan.json', dict(commands=commands, outputs=outputs, inputs=inputs))
    actions = json.loads((W / 'actions.json').read_text())
    for name in commands:
        assert name not in actions and not (B / outputs[name]).exists()
        actions[name] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
            kind='diagnostic', expected_exit=0, software_provider=False,
            limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
            argv=['/usr/bin/python3', '-B', str(W / 'remaining.py'), name],
            frozen_inputs=sorted(inputs) + ['out/native-cpu-strict-ir-linking/remaining-plan.json'] +
                          ([] if name == 'remaining-plugin' else [str((B / outputs['remaining-plugin']).relative_to(R))]),
            scope='First fixed Transport/Endpoint strict conversion, with full input fingerprint and reviewed scalar coefficient cast mapping. Existing strict FP/masks/layout/counts/work limits unchanged. No numerical run, native-provider/frame or release qualification.')
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=list(commands), compiler_run=False)))
else:
    plan = json.loads((W / 'remaining-plan.json').read_text())
    assert phase in plan['commands']
    for name, expected in plan['inputs'].items():
        assert pin(R / name) == expected, name
    if phase != 'remaining-plugin':
        owner = json.loads((W / 'remaining-plugin/owner.json').read_text())
        assert owner['passed'] and owner['status'] == 'terminal'
    target = B / plan['outputs'][phase]
    assert not target.exists()
    subprocess.run(plan['commands'][phase], check=True)
    result = dict(output=dict(path=str(target.relative_to(R)), **pin(target)), argv=plan['commands'][phase])
    if phase == 'remaining-plugin':
        undefined = {line.split()[-1] for line in subprocess.check_output(
            ['nm', '-D', str(target)], text=True).splitlines() if ' U ' in line and '_ZN4llvm' in line}
        exported = {line.split()[-1] for line in subprocess.check_output(
            ['nm', '-D', str(D / 'opt')], text=True).splitlines() if len(line.split()) == 3}
        assert not undefined - exported, sorted(undefined - exported)
    else:
        ll = target.with_suffix('.ll')
        assert not ll.exists()
        subprocess.run([str(D / 'llvm-dis'), str(target), '-o', str(ll)], check=True)
        text = ll.read_text()
        ordinary = re.findall(r'^\s+%[^=]+ = (?:fadd|fsub|fmul|fdiv|frem|fcmp|fpext|fptrunc|sitofp|uitofp|fptosi|fptoui)\b', text, re.M)
        assert not ordinary, ordinary[:3]
        calls = Counter(re.findall(r'\bcall\s+.*?@(llvm\.experimental\.constrained\.[^ (]+)\(', text))
        expected = ({'fadd.v8f32': 82, 'fsub.v8f32': 13, 'fcmp.v8f32': 349,
                     'fcmps.v8f32': 39, 'sitofp.f32.i32': 4} if phase == 'constrain-transport' else
                    {'fadd.v8f32': 99, 'fsub.v8f32': 13, 'fmul.v8f32': 1,
                     'fcmp.v8f32': 453, 'fcmps.v8f32': 31})
        assert dict(calls) == {'llvm.experimental.constrained.' + k: v for k, v in expected.items()}, calls
        result.update(ir=dict(path=str(ll.relative_to(R)), **pin(ll)),
            constrained_calls=dict(calls), ordinary_fp_operations=0,
            scalar_coefficient_casts=4 if phase == 'constrain-transport' else 0,
            scope='Exact-stage strict IR/verifier gate only; value-ABI integration, native ISA/retention, full finite compatibility, provider and renderer qualification remain pending.')
    write(W / (phase + '-result.json'), result)
    print(json.dumps({k: v for k, v in result.items() if k != 'argv'}))
