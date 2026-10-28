"""Fixed context promotion/cleanup/native gates for the two remaining stages."""
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


def terminal(action):
    path = W / action / 'owner.json'
    value = doc(path)
    assert value['status'] == 'terminal' and value['passed'] and value['child_exit'] == 0
    assert value['source_unchanged'] and value['whole_inputs_unchanged']
    assert not value['remaining_owned_processes'] and not value['cleanup_errors']
    assert value['observed_births'] and all(row['absent'] for row in value['observed_births'])
    assert path.with_name('inputs-before.json').read_bytes() == path.with_name('inputs-after.json').read_bytes()


stage, phase = sys.argv[1:]
assert stage in ('transport', 'endpoint')
B = R / ('bin/linux-gcc/native-cpu-canonical-closed-' + stage)
contract_path = V / ('closed_' + stage + '_context_contract.json')
contract = doc(contract_path)
assert contract['stage'] == stage and contract['phase'] == 'context'
plan_path = W / (stage + '-context-plan.json')
steps = ('context-plugin', 'context', 'cleanup', 'native')
if phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    terminal('value-' + stage)
    source = R / contract['output']['path']
    assert pin(source) == {key:contract['output'][key] for key in ('bytes', 'sha256')}
    value = doc(W / (stage + '-value-result.json'))
    assert contract['input'] == value['output']
    assert pin(R / value['output']['path']) == {key:value['output'][key] for key in ('bytes', 'sha256')}
    config = shlex.split(subprocess.check_output([str(D / 'llvm-config'), '--cxxflags'], text=True))
    pipeline = 'function(sroa,mem2reg,instcombine,simplifycfg,early-cse<memssa>,dse),globaldce,verify'
    outputs = dict(zip(steps, ('context.so', stage + '.context.bc', stage + '.closed.bc', stage + '.closed.o')))
    commands = {
        'context-plugin': [str(D / 'clang++'), *config, '-O2', '-fPIC', '-shared', str(source),
            str(SDK / 'sha256-only.a'), '-Wl,--exclude-libs,ALL', '-o', str(B / outputs['context-plugin'])],
        'context': [str(D / 'opt'), '-load-pass-plugin=' + str(B / outputs['context-plugin']),
            '-passes=' + contract['pass_name'] + ',verify', '-verify-each', str(R / contract['input']['path']), '-o', str(B / outputs['context'])],
        'cleanup': [str(D / 'opt'), '-passes=' + pipeline, '-verify-each', '-mcpu=skylake-avx512',
            str(B / outputs['context']), '-o', str(B / outputs['cleanup'])],
        'native': [str(D / 'llc'), '-O2', '-filetype=obj', '-relocation-model=pic', '-mcpu=skylake-avx512',
            str(B / outputs['cleanup']), '-o', str(B / outputs['native'])]}
    assert all(not (B / name).exists() for name in outputs.values())
    paths = [*W.glob('*.py'), source, contract_path, V / 'derive_closed_remaining_context.py',
        V / 'toolchain.json', SDK / 'sha256-only.a',
        R / 'out/native-cpu-field-compatibility/elf_proof.py',
        R / 'out/native-cpu-accumulate-compatibility/elf_proof.py',
        Path('/usr/bin/readelf').resolve(), Path('/usr/bin/objdump').resolve()]
    inputs = dict(contract['source_inputs'])
    inputs.update({str(path.relative_to(R)) if path.is_relative_to(R) else str(path):pin(path) for path in paths})
    for name in ('clang++', 'opt', 'llvm-config', 'llvm-dis', 'llc'):
        tool = doc(V / 'toolchain.json')['tools'][name]
        path = R / tool['resolved']
        assert pin(path) == {key:tool[key] for key in ('bytes', 'sha256')}
        inputs[str(path.relative_to(R))] = pin(path)
    for name, expected in inputs.items():
        assert pin(R / name) == expected, name
    value_contract = doc(V / ('closed_' + stage + '_value_contract.json'))
    staged_words = value_contract['context']['output_words'] * 8
    before = (B / (stage + '.value.ll')).read_text()
    arrays = re.findall(r'^\s+%staged_output = alloca \[(\d+) x i32\], align 64$', before, re.M)
    assert arrays == [str(staged_words), str(staged_words)]
    write(plan_path, dict(commands=commands, outputs=outputs, inputs=inputs, cleanup_pipeline=pipeline,
        contract=pin(contract_path), transactional_staging_words=staged_words))
    actions_path = W / 'actions.json'
    actions = doc(actions_path)
    prior_path = W / (stage + '-context-prior-actions.json')
    with prior_path.open('xb') as stream:
        stream.write(actions_path.read_bytes())
    for index, step in enumerate(steps):
        action = step + '-' + stage
        assert action not in actions
        extra = []
        if index:
            previous = steps[index - 1]
            extra += [str((B / outputs[previous]).relative_to(R)),
                str((W / (stage + '-' + previous + '-result.json')).relative_to(R)),
                str((W / (previous + '-' + stage) / 'owner.json').relative_to(R))]
        if step == 'native':
            extra += [str((B / (stage + '.closed.ll')).relative_to(R))]
        actions[action] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
            kind='diagnostic', expected_exit=0, software_provider=False,
            limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
            argv=['/usr/bin/python3', '-B', str(Path(__file__).resolve()), stage, step],
            frozen_inputs=sorted(inputs) + [str(plan_path.relative_to(R)), str(prior_path.relative_to(R))] + extra,
            scope='One exact remaining-stage private context compiler step. Preserve original transaction staging, coupled fields, mixed masks, initialization/refusal/order/FP/math and governors. Only proved reset bookkeeping expansion and bounded coefficient isolation are permitted. Final entry requires asynchronous unwind metadata from first native emission. No provider/runtime/numerical/benefit/frame/scientific/release qualification.')
    actions_path.write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=[step + '-' + stage for step in steps], transactional_staging_bytes=staged_words * 4)))
else:
    assert phase in steps
    plan = doc(plan_path)
    assert plan['contract'] == pin(contract_path)
    for name, expected in plan['inputs'].items():
        assert pin(R / name) == expected, name
    index = steps.index(phase)
    if index:
        previous = steps[index - 1]
        terminal(previous + '-' + stage)
        product = doc(W / (stage + '-' + previous + '-result.json'))['output']
        assert pin(R / product['path']) == {key:product[key] for key in ('bytes', 'sha256')}
    output = B / plan['outputs'][phase]
    assert not output.exists()
    subprocess.run(plan['commands'][phase], check=True)
    result = dict(phase=phase, stage=stage, output=dict(path=str(output.relative_to(R)), **pin(output)),
        argv=plan['commands'][phase], numerical_qualification_transferred=False)
    if output.suffix == '.bc':
        ll = output.with_suffix('.ll')
        assert not ll.exists()
        subprocess.run([str(D / 'llvm-dis'), str(output), '-o', str(ll)], check=True)
        text = ll.read_text()
        headers = re.findall(r'^define .*$', text, re.M)
        attrs = dict(re.findall(r'^attributes #(\d+) = (.*)$', text, re.M))
        assert len(headers) == 1 and '@sirius_closed_' + stage + '_packet(' in headers[0]
        attributes = attrs[re.search(r'#(\d+)', headers[0])[1]]
        assert 'strictfp' in attributes and 'uwtable' in attributes and 'uwtable(sync)' not in attributes
        assert not re.findall(r'^\s+%[^=]+ = (?:fadd|fsub|fmul|fdiv|frem|fcmp|fpext|fptrunc|sitofp|uitofp|fptosi|fptoui)\b', text, re.M)
        assert not re.findall(r'^\s+.*\b(?:fadd|fsub|fmul|fdiv|fcmp|call)\s+(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b', text, re.M)
        constrained = Counter(re.findall(r'\bcall\s+.*?@(llvm.experimental.constrained.[^ (]+)\(', text))
        for line in text.splitlines():
            if 'call ' in line and '@llvm.experimental.constrained.' in line:
                assert 'metadata !"fpexcept.strict"' in line
                assert '.fcmp' in line or 'metadata !"round.dynamic"' in line
        calls = Counter(re.findall(r'\bcall\s+.*?@(sirius_local_value_[^ (]+)\(', text))
        assert calls and not re.findall(r'\bcall\s+.*?@sirius_closed_' + stage + r'_scope_v4_', text)
        allocas = re.findall(r'^\s+%[^=]+ = alloca .*$', text, re.M)
        staged = [row for row in allocas if row.strip().startswith('%staged_output =')]
        assert len(staged) == 1 and staged[0].strip() == '%staged_output = alloca [' + str(plan['transactional_staging_words']) + ' x i32], align 64'
        context_allocas = [row for row in allocas if 'alloca %v8_varying_KernelContext_0' in row]
        other_allocas = [row.strip() for row in allocas if row not in staged]
        result.update(ir=dict(path=str(ll.relative_to(R)), **pin(ll)), definitions=len(headers),
            constrained_calls=dict(constrained), value_calls=dict(calls), owned_context_allocas=len(context_allocas),
            memory_instruction_counts=dict(Counter(re.findall(r'^\s+(?:%[^=]+ = )?(load|store|alloca)\b', text, re.M))),
            masked_store_calls=text.count('call void @llvm.masked.store'), allocas=[row.strip() for row in allocas],
            original_transactional_staging_bytes=plan['transactional_staging_words'] * 4)
        if phase == 'cleanup':
            assert not context_allocas, 'Required main-bank Context promotion failed; stop before native.'
            if stage == 'endpoint':
                assert not other_allocas, ('Unexpected Endpoint residual private storage', other_allocas)
            else:
                assert len(other_allocas) <= 1 and all(re.search(
                    r'alloca (?:\[125 x <8 x i32>\]|\[4000 x i8\]), align 64$', row)
                    for row in other_allocas), other_allocas
    elif phase == 'native':
        from native_check import check
        result.update(check(output, 'sirius_closed_' + stage + '_packet'))
        result['cleaned_IR_value_calls'] = doc(W / (stage + '-cleanup-result.json'))['value_calls']
    write(W / (stage + '-' + phase + '-result.json'), result)
    print(json.dumps({key:result[key] for key in ('phase', 'stage', 'output', 'owned_context_allocas', 'memory_instruction_counts', 'instructions') if key in result}))
