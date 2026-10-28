"""One source derivation, LLVM verification, and native adapter emission; no workload."""
from pathlib import Path
import json
import subprocess
import sys
from math_source import R, W, B, MATH, RECIPES, pin, sources

D = R / 'bin/linux-gcc/native-cpu-strict-ir-linking/toolchain/bin'
TOOLS = R / 'out/native-cpu-strict-ir-linking/toolchain.json'
RUN = W / 'math'


def doc(path):
    return json.loads(path.read_text())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def terminal(name):
    p = RUN / name / 'owner.json'
    d = doc(p)
    assert d['status'] == 'terminal' and d['passed'] and d['child_exit'] == 0
    assert d['source_unchanged'] and d['whole_inputs_unchanged']
    assert not d['remaining_owned_processes'] and not d['cleanup_errors']
    assert all(row['absent'] for row in d['observed_births'])


phase, = sys.argv[1:]
if phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    assert not (RUN / 'actions.json').exists() and not (B / 'wide-math.ll').exists()
    sources()
    commands = {
        'math-source': ['/usr/bin/python3', '-B', str(W / 'math_source.py'), 'derive'],
        'math-verify': [str(D / 'opt'), '-passes=verify', '-verify-each', str(B / 'wide-math.ll'), '-o', str(B / 'wide-math.bc')],
        'math-native': [str(D / 'llc'), '-O2', '-filetype=obj', '-relocation-model=pic', '-mcpu=x86-64-v4',
                        str(B / 'wide-math.bc'), '-o', str(B / 'wide-math.o')],
    }
    inputs = {str(p.relative_to(R)): pin(p) for p in [
        *(W / name for name in ('math_source.py', 'build_math.py', 'native_math_check.py')),
        RUN / 'run_owned.py', RUN / 'owned_guard.py',
        TOOLS, MATH / 'math.bc', MATH / 'math.ll', MATH / 'math.o',
        R / 'out/native-cpu-canonical-closed-remaining/native_check.py',
        R / 'out/native-cpu-field-compatibility/elf_proof.py',
        R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
        R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py']}
    for name in ('opt', 'llc'):
        value = doc(TOOLS)['tools'][name]
        path = R / value['resolved']
        assert pin(path) == {k:value[k] for k in ('bytes', 'sha256')}
        inputs[str(path.relative_to(R))] = pin(path)
    for tool in ('objdump', 'readelf'):
        path = Path('/usr/bin') / tool
        inputs[str(path)] = pin(path)
    write(W / 'math-plan.json', dict(inputs=inputs, commands=commands))
    actions = {}
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    extra = []
    for name in commands:
        if name == 'math-verify':
            extra = [str((B / 'wide-math.ll').relative_to(R)), str((W / 'math-source-contract.json').relative_to(R)),
                     str((RUN / 'math-source/owner.json').relative_to(R))]
        elif name == 'math-native':
            extra += [str((B / 'wide-math.bc').relative_to(R)), str((W / 'math-verify-result.json').relative_to(R)),
                      str((RUN / 'math-verify/owner.json').relative_to(R))]
        actions[name] = dict(revision=revision, kind='diagnostic', expected_exit=0, software_provider=False,
            limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
            argv=['/usr/bin/python3', '-B', str(W / 'build_math.py'), name],
            frozen_inputs=sorted(inputs) + [str((W / 'math-plan.json').relative_to(R))] + list(extra),
            scope='Source/LLVM/native adapter check only. Original strict8 math object stays exact. '
                  'Async FDE coverage required from first native emission; full GP/opcode/link/admission and raw16 '
                  'subgroup/initialization/context/runtime proofs remain required. No numerical, speed, provider or frame acceptance.')
    # Keep this compiler owner's action file fixed when raw stages advance.
    write(RUN / 'actions.json', actions)
    print(json.dumps(dict(prepared=list(commands), workload=False)))
else:
    plan = doc(W / 'math-plan.json')
    assert phase in plan['commands']
    for name, expected in plan['inputs'].items():
        assert pin(R / name) == expected, name
    order = list(plan['commands'])
    i = order.index(phase)
    if i:
        terminal(order[i - 1])
        previous = doc(W / (order[i - 1] + '-result.json'))['output']
        assert pin(R / previous['path']) == {k:previous[k] for k in ('bytes', 'sha256')}
    output = B / {'math-source':'wide-math.ll', 'math-verify':'wide-math.bc', 'math-native':'wide-math.o'}[phase]
    assert not output.exists()
    subprocess.run(plan['commands'][phase], check=True)
    result = dict(output=dict(path=str(output.relative_to(R)), **pin(output)), passed=True, workload=False)
    if phase == 'math-native':
        from native_math_check import all_adapters
        names = ['sirius_wide_value_' + name + '16' for name in RECIPES]
        checks = all_adapters(output, names)
        for recipe, name in zip(RECIPES, names):
            assert checks[name]['core_relocations'] == {'sirius_local_value_' + recipe + '8': 2}
        result['native'] = checks
    write(W / (phase + '-result.json'), result)
    print(json.dumps({k:result[k] for k in ('passed', 'output', 'workload')}))
