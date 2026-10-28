"""Run fixed strict/value stage passes once, with immutable source inventories."""
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
V = R / 'out/native-cpu-strict-ir-linking'
SDK = R / 'bin/linux-gcc/native-cpu-strict-ir-linking'
D = SDK / 'toolchain/bin'


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def doc(path):
    return json.loads(path.read_text())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def terminal(name):
    path = W / name / 'owner.json'
    value = doc(path)
    assert value['status'] == 'terminal' and value['passed'] and value['child_exit'] == 0
    assert value['source_unchanged'] and value['whole_inputs_unchanged']
    assert not value['remaining_owned_processes'] and not value['cleanup_errors']
    assert value['observed_births'] and all(row['absent'] for row in value['observed_births'])
    assert path.with_name('inputs-before.json').read_bytes() == path.with_name('inputs-after.json').read_bytes()


stage, phase = sys.argv[1:]
assert stage in ('transport', 'endpoint')
mode = phase.removeprefix('prepare-').removesuffix('-plugin')
assert mode in ('scope', 'value')
B = R / ('bin/linux-gcc/native-cpu-canonical-closed-' + stage)
contract_path = V / ('closed_' + stage + '_' + mode + '_contract.json')
contract = doc(contract_path)
assert contract['stage'] == stage and contract['phase'] == mode and contract['source_inverse_exact']
action = mode + '-' + stage
plugin_action = mode + '-plugin-' + stage
plan_path = W / (stage + '-' + mode + '-plan.json')
if phase.startswith('prepare-'):
    assert phase == 'prepare-' + mode
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    terminal('emit-' + stage)
    if mode == 'value':
        terminal('scope-' + stage)
    source = R / contract['output']['path']
    assert pin(source) == {key:contract['output'][key] for key in ('bytes', 'sha256')}
    assert pin(R / contract['input']['path']) == {key:contract['input'][key] for key in ('bytes', 'sha256')}
    config = shlex.split(subprocess.check_output([str(D / 'llvm-config'), '--cxxflags'], text=True))
    outputs = {plugin_action: mode + '.so', action: stage + ('.strict.bc' if mode == 'scope' else '.value.bc')}
    assert all(not (B / name).exists() for name in outputs.values())
    commands = {
        plugin_action: [str(D / 'clang++'), *config, '-O2', '-fPIC', '-shared', str(source),
            str(SDK / 'sha256-only.a'), '-Wl,--exclude-libs,ALL', '-o', str(B / outputs[plugin_action])],
        action: [str(D / 'opt'), '-load-pass-plugin=' + str(B / outputs[plugin_action]),
            '-passes=' + contract['pass_name'] + ',verify', '-verify-each',
            str(R / contract['input']['path']), '-o', str(B / outputs[action])]}
    paths = [*W.glob('*.py'), V / 'derive_closed_remaining_passes.py', source, contract_path,
        SDK / 'sha256-only.a', V / 'toolchain.json']
    inputs = dict(contract['source_inputs'])
    inputs.update({str(path.relative_to(R)):pin(path) for path in paths})
    for name in ('clang++', 'opt', 'llvm-config', 'llvm-dis'):
        tool = doc(V / 'toolchain.json')['tools'][name]
        path = R / tool['resolved']
        assert pin(path) == {key:tool[key] for key in ('bytes', 'sha256')}
        inputs[str(path.relative_to(R))] = pin(path)
    for name, expected in inputs.items():
        assert pin(R / name) == expected, name
    write(plan_path, dict(commands=commands, outputs=outputs, inputs=inputs, contract=pin(contract_path)))
    actions_path = W / 'actions.json'
    actions = doc(actions_path)
    prior = W / (stage + '-' + mode + '-prior-actions.json')
    with prior.open('xb') as stream:
        stream.write(actions_path.read_bytes())
    for name in commands:
        assert name not in actions
        extra = [] if name == plugin_action else [str((B / outputs[plugin_action]).relative_to(R)),
            str((W / (stage + '-' + mode + '-plugin-result.json')).relative_to(R)),
            str((W / plugin_action / 'owner.json').relative_to(R))]
        actions[name] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
            kind='diagnostic', expected_exit=0, software_provider=False,
            limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
            argv=['/usr/bin/python3', '-B', str(Path(__file__).resolve()), stage, mode + ('-plugin' if name == plugin_action else '')],
            frozen_inputs=sorted(inputs) + [str(plan_path.relative_to(R)), str(prior.relative_to(R))] + extra,
            scope='One exact closed mixed-stage strict/value pass build or execution. Source/IR/layout/cast/table/mask/rounding/exception/FP-effect inventories remain pinned. Value wrappers are counted before inlining, preserving raw call multiplicity. No context/native/FDE/provider/workload/numerical/benefit/frame/scientific/release qualification.')
    actions_path.write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=list(commands), runtime_run=False)))
else:
    assert phase in (mode, mode + '-plugin')
    plan = doc(plan_path)
    assert plan['contract'] == pin(contract_path)
    for name, expected in plan['inputs'].items():
        assert pin(R / name) == expected, name
    selected = plugin_action if phase.endswith('-plugin') else action
    if selected == action:
        terminal(plugin_action)
        product = doc(W / (stage + '-' + mode + '-plugin-result.json'))['output']
        assert pin(R / product['path']) == {key:product[key] for key in ('bytes', 'sha256')}
    output = B / plan['outputs'][selected]
    assert not output.exists()
    subprocess.run(plan['commands'][selected], check=True)
    result = dict(phase=selected, output=dict(path=str(output.relative_to(R)), **pin(output)),
        input=contract['input'], scope_contract=pin(V / ('closed_' + stage + '_scope_contract.json')),
        argv=plan['commands'][selected])
    if selected == action:
        ll = output.with_suffix('.ll')
        assert not ll.exists()
        subprocess.run([str(D / 'llvm-dis'), str(output), '-o', str(ll)], check=True)
        text = ll.read_text()
        headers = re.findall(r'^define .*$', text, re.M)
        attrs = dict(re.findall(r'^attributes #(\d+) = (.*)$', text, re.M))
        assert headers and all('strictfp' in attrs[re.search(r'#(\d+)', header)[1]] for header in headers)
        assert not re.findall(r'^\s+%[^=]+ = (?:fadd|fsub|fmul|fdiv|frem|fcmp|fpext|fptrunc|sitofp|uitofp|fptosi|fptoui)\b', text, re.M)
        assert not re.findall(r'^\s+.*\b(?:fneg|fadd|fsub|fmul|fdiv|fcmp|call)\s+(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b', text, re.M)
        constrained = Counter(re.findall(r'\bcall\s+.*?@(llvm.experimental.constrained.[^ (]+)\(', text))
        if mode == 'scope':
            shape = contract['actual_inventory']
            expected = Counter({'llvm.experimental.constrained.' + op + '.v8f32': shape['ops'].get(op, 0)
                for op in ('fadd', 'fsub', 'fmul')})
            expected['llvm.experimental.constrained.fcmp.v8f32'] = sum(shape['predicates'].get(p, 0) for p in ('oeq', 'une', 'uno'))
            expected['llvm.experimental.constrained.fcmps.v8f32'] = sum(shape['predicates'].get(p, 0) for p in ('oge', 'ole', 'ogt', 'olt'))
            if shape['ops'].get('sitofp', 0):
                expected['llvm.experimental.constrained.sitofp.v8f32.v8i32'] = shape['ops']['sitofp']
            assert constrained == +expected, (constrained, +expected)
            assert len(headers) == shape['definitions']
        else:
            assert constrained == Counter(doc(W / (stage + '-scope-result.json'))['constrained_calls'])
        calls = Counter(re.findall(r'\bcall\s+.*?@([^ (]+)\(', text))
        original_inventory = doc(W / (stage + '-inventory.json'))
        negates = len(re.findall(r'= fneg ', text))
        assert negates == original_inventory['ordinary_fp_operations'].get('fneg', 0)
        assert calls.get('llvm.fabs.v8f32', 0) == original_inventory['call_inventory'].get('llvm.fabs.v8f32', 0)
        for line in text.splitlines():
            if 'call ' in line and '@llvm.experimental.constrained.' in line:
                assert 'metadata !"fpexcept.strict"' in line
                assert '.fcmp' in line or 'metadata !"round.dynamic"' in line
        pointer = Counter({name:count for name,count in calls.items() if name.startswith(contract['pointer_math_prefix'])})
        typed = Counter({name:count for name,count in calls.items() if name.startswith('sirius_local_value_')})
        recipes = contract['exact_calls_by_recipe'] if mode == 'value' else {
            name.removeprefix('sirius_native_v4_'):count for name,count in contract['actual_inventory']['math_calls'].items()}
        assert pointer == Counter({contract['pointer_math_prefix'] + name:count for name,count in recipes.items()})
        assert not any(name.startswith('sirius_native_v4_') for name in calls)
        if mode == 'value':
            assert typed == Counter({'sirius_local_value_' + name:1 for name,count in recipes.items() if count})
            adapters = [header for header in headers if header.startswith('define internal ') and '@' + contract['pointer_math_prefix'] in header]
            assert len(adapters) == len(typed)
            assert all('alwaysinline' in attrs[re.search(r'#(\d+)', header)[1]] for header in adapters)
            result.update(source_value_contract=pin(contract_path),
                accepted_artifact_owner=str((W / action / 'owner.json').relative_to(R)))
        else:
            assert not typed
        result.update(ir=dict(path=str(ll.relative_to(R)), **pin(ll)), definitions=len(headers),
            constrained_calls=dict(constrained), pointer_calls=dict(pointer), value_calls=dict(typed),
            ordinary_fneg=negates)
    write(W / (stage + '-' + mode + ('-plugin' if selected == plugin_action else '') + '-result.json'), result)
    print(json.dumps({key:result[key] for key in ('phase', 'output', 'definitions') if key in result}))
