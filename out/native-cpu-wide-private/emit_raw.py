"""Derive private subgroup-preserving raw16 source and emit each stage once."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys
from raw_source import R, W, B, planned, pin

RUN = W / 'raw'
ORDER = ('derive', 'emit-transport', 'emit-endpoint', 'emit-dopri_phase')


def doc(path):
    return json.loads(path.read_text())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def terminal(folder):
    p = folder / 'owner.json'
    d = doc(p)
    assert d['passed'] and d['status'] == 'terminal' and d['child_exit'] == 0
    assert d['source_unchanged'] and d['whole_inputs_unchanged']
    assert d['observed_births'] and all(row['absent'] for row in d['observed_births'])
    assert not d['remaining_owned_processes'] and not d['cleanup_errors']
    assert p.with_name('inputs-before.json').read_bytes() == p.with_name('inputs-after.json').read_bytes()


def math_interface():
    terminal(W / 'math/math-native')
    result = doc(W / 'math-native-result.json')
    product = result['output']
    assert pin(R / product['path']) == {k:product[k] for k in ('bytes', 'sha256')}
    assert len(result['native']) == 9
    for name, row in result['native'].items():
        assert row['unwind_table']['full_function_FDE'] and not row['emitted_fused_instructions']
        assert sum(row['core_relocations'].values()) == 2
    return product


phase, = sys.argv[1:]
if phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    assert not (RUN / 'actions.json').exists() and not (W / 'raw-source-contract.json').exists()
    outputs, contract = planned()
    assert not any(path.exists() for path in outputs)
    math = math_interface()
    inputs = dict(contract['source_inputs'])
    for path in (W / 'emit_raw.py', RUN / 'run_owned.py', RUN / 'owned_guard.py',
                 W / 'math-native-result.json', W / 'math/math-native/owner.json', R / math['path'],
                 R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
                 R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py'):
        inputs[str(path.relative_to(R))] = pin(path)
    commands = {'derive': ['/usr/bin/python3', '-B', str(W / 'raw_source.py'), 'derive']}
    commands.update({'emit-' + stage: argv for stage, argv in contract['prospective_emit_commands'].items()})
    assert list(commands) == list(ORDER)
    write(W / 'raw-plan.json', dict(commands=commands, inputs=inputs,
        expected_source_outputs=contract['outputs'], math_interface=math,
        scope='Source/ordinary IR only. Strict16/value/context/native/private-admission/runtime gates remain open.'))
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    actions = {}
    for i, name in enumerate(ORDER):
        extra = []
        if i:
            extra += list(contract['outputs']) + [str((W / 'raw-source-contract.json').relative_to(R)),
                str((W / 'derive-result.json').relative_to(R)), str((RUN / 'derive/owner.json').relative_to(R))]
        if i > 1:
            previous = ORDER[i - 1]
            extra += [str((B / (previous.removeprefix('emit-') + '.ll')).relative_to(R)),
                      str((W / (previous + '-result.json')).relative_to(R)),
                      str((RUN / previous / 'owner.json').relative_to(R))]
        actions[name] = dict(revision=revision, kind='diagnostic', expected_exit=0, software_provider=False,
            limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
            argv=['/usr/bin/python3', '-B', str(W / 'emit_raw.py'), name],
            frozen_inputs=sorted(inputs) + [str((W / 'raw-plan.json').relative_to(R))] + extra,
            scope='One ordered raw16 source or ordinary ISPC emission, no retries. Exact9380 tuples, '
                  'eight-row subgroup decisions, original physical fields/RHS/trials/quartics, full initialization, '
                  'I/O layout, refusal and publication order retained. Only compiler-private width/packing changes. '
                  'No strict/native/provider/numerical/benefit/frame or scientific/release acceptance.')
    write(RUN / 'actions.json', actions)
    print(json.dumps(dict(prepared=list(commands), workload=False)))
else:
    plan = doc(W / 'raw-plan.json')
    assert phase in ORDER
    for name, expected in plan['inputs'].items():
        assert pin(R / name) == expected, name
    math_interface()
    i = ORDER.index(phase)
    if i:
        terminal(RUN / 'derive')
        source = doc(W / 'raw-source-contract.json')
        assert source['outputs'] == plan['expected_source_outputs']
        for name, expected in source['outputs'].items():
            assert pin(R / name) == expected, name
        if i > 1:
            previous = ORDER[i - 1]
            terminal(RUN / previous)
            product = doc(W / (previous + '-result.json'))['output']
            assert pin(R / product['path']) == {k:product[k] for k in ('bytes', 'sha256')}
        output = B / (phase.removeprefix('emit-') + '.ll')
        assert not output.exists()
    subprocess.run(plan['commands'][phase], check=True)
    if not i:
        contract = doc(W / 'raw-source-contract.json')
        assert contract['outputs'] == plan['expected_source_outputs']
        result = dict(source_only=True, outputs=contract['outputs'], compiler_run=False)
    else:
        text = output.read_text()
        assert not re.findall(r'^\s+.*\b(?:fadd|fsub|fmul|fdiv|fcmp|call)\s+(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b', text, re.M)
        assert not re.search(r'@llvm\.(?:fma|fmuladd)\.', text)
        calls = Counter(re.findall(r'\bcall\s+.*?@([^ (]+)\(', text))
        math = {key:value for key,value in calls.items() if key.startswith('sirius_wide_pointer_')}
        assert math and all(key.endswith('16') for key in math)
        result = dict(output=dict(path=str(output.relative_to(R)), **pin(output)),
            definitions=len(re.findall(r'^define ', text, re.M)), pointer_math_calls=math,
            ordinary_FP=dict(Counter(re.findall(r'^\s+%[^=]+ = (fadd|fsub|fmul|fdiv|frem|fneg|fcmp|sitofp|uitofp|fptosi|fptoui|fpext|fptrunc)\b', text, re.M))),
            comparison_predicates=dict(Counter(re.findall(r'= fcmp (\w+) ', text))),
            context_types=re.findall(r'^%v16_varying_KernelContext_0 = type .+$', text, re.M),
            scope='Actual ordinary raw16 IR only. Fresh strict/model/interface/context/native proofs are required; no numerical transfer.')
        assert len(result['context_types']) == 1
    write(W / (phase + '-result.json'), result)
    print(json.dumps(result))
