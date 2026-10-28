"""Build one shared strict value-math module and the three bounded raw modules."""
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
CORE_NAMES = ['sirius_shared_value_' + name for name in (
    'dense_add8', 'dense_multiply8', 'dense_divide8', 'dense_sqrt8',
    'tagged_add8', 'tagged_subtract8', 'tagged_multiply8', 'tagged_divide8', 'tagged_sqrt8')]
INPUTS = dict(dp='dopri_phase.strict.sha256.bc', transport='transport.strict.bc',
              endpoint='endpoint.strict.bc', math='private-math.bc')


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,
                    sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


phase, = sys.argv[1:]
if phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    assert json.loads((W / 'constrain-endpoint/owner.json').read_text())['passed']
    config = shlex.split(subprocess.check_output([str(D / 'llvm-config'), '--cxxflags'], text=True))
    commands = {'shared-plugin': [str(D / 'clang++'), *config, '-O2', '-fPIC', '-shared',
                 str(W / 'shared_value_abi.cpp'), str(B / 'sha256-only.a'),
                 '-Wl,--exclude-libs,ALL', '-o', str(B / 'shared_value_abi.so')]}
    outputs = {'shared-plugin': 'shared_value_abi.so'}
    exports = {}
    for stage, filename in INPUTS.items():
        if stage == 'math':
            exports[stage] = CORE_NAMES
        else:
            ll = B / filename.replace('.bc', '.ll')
            if stage == 'dp':
                # DP conversion was emitted as bitcode only. Its original two
                # export names were already sealed by the prior link plan.
                exports[stage] = json.loads((W / 'link-plan.json').read_text())['exports']
            else:
                exports[stage] = re.findall(r'^define (?!internal).*?@(sirius_scope_[^ (]+)\(', ll.read_text(), re.M)
                assert len(exports[stage]) == 2
        lower, opt, native = ('shared-' + stage + '-' + step for step in ('lower', 'opt', 'native'))
        outputs.update({lower: 'shared-' + stage + '.bc', opt: 'shared-' + stage + '.O2.bc',
                        native: 'shared-' + stage + '.o'})
        commands[lower] = [str(D / 'opt'), '-load-pass-plugin=' + str(B / outputs['shared-plugin']),
                           '-passes=sirius-shared-' + stage + ',verify', '-verify-each',
                           str(B / filename), '-o', str(B / outputs[lower])]
        commands[opt] = [str(D / 'opt'), '-passes=internalize,default<O2>,function(sroa),globaldce,verify',
                         '-verify-each', '-internalize-public-api-list=' + ','.join(exports[stage]),
                         '-mcpu=' + ('x86-64-v4' if stage == 'math' else 'skylake-avx512'),
                         str(B / outputs[lower]), '-o', str(B / outputs[opt])]
        commands[native] = [str(D / 'llc'), '-O2', '-filetype=obj', '-relocation-model=pic',
                            '-mcpu=' + ('x86-64-v4' if stage == 'math' else 'skylake-avx512'),
                            str(B / outputs[opt]), '-o', str(B / outputs[native])]
    paths = [*W.glob('*.py'), W / 'shared_value_abi.cpp', W / 'shared_value_contract.json',
             W / 'value_abi.cpp', W / 'scope_contract.json', W / 'remaining_scope_review.json',
             W / 'toolchain.json', W / 'link-plan.json', B / 'sha256-only.a',
             R / 'out/native-cpu-exact-product-zero-residual/retained8.cpp',
             R / 'out/native-cpu-inline-private-bank/packinghelper.isph',
             R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
             R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py']
    paths += [B / name for name in INPUTS.values()]
    paths += [R / ('bin/linux-gcc/native-cpu-inline-private-bank/ispc-source/' + s + '.ispc')
              for s in ('dopri_phase', 'transport', 'endpoint')]
    for name in ('clang++', 'opt', 'llvm-dis', 'llvm-config', 'llc'):
        row = json.loads((W / 'toolchain.json').read_text())['tools'][name]
        assert pin(R / row['resolved']) == {k: row[k] for k in ('bytes', 'sha256')}
        paths.append(R / row['resolved'])
    inputs = {str(p.relative_to(R)): pin(p) for p in paths}
    write(W / 'shared-plan.json', dict(commands=commands, outputs=outputs, inputs=inputs, exports=exports,
        core_names=CORE_NAMES, scope='Separate raw3 and one shared9-core module; fixed external hidden fastcc signatures, strict FP, no contraction, full original target feature closure.'))
    actions = json.loads((W / 'actions.json').read_text())
    for name in commands:
        assert name not in actions and not (B / outputs[name]).exists()
        dependencies = []
        if name != 'shared-plugin':
            dependencies.append('shared_value_abi.so')
            if not name.endswith('-lower'):
                prior = name.rsplit('-', 1)[0] + ('-lower' if name.endswith('-opt') else '-opt')
                dependencies.append(outputs[prior])
        actions[name] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
            kind='diagnostic', expected_exit=0, software_provider=False,
            limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
            argv=['/usr/bin/python3', '-B', str(W / 'shared.py'), name],
            frozen_inputs=sorted(inputs) + ['out/native-cpu-strict-ir-linking/shared-plan.json'] +
                          [str((B / n).relative_to(R)) for n in dependencies],
            scope='One shared strict value-math and exact raw-stage implementation gate; static compiler/IR/native evidence only. No numerical, provider, original-frame, scientific or release qualification.')
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=list(commands), compiler_run=False)))
else:
    plan = json.loads((W / 'shared-plan.json').read_text())
    assert phase in plan['commands']
    for name, expected in plan['inputs'].items():
        assert pin(R / name) == expected, name
    target = B / plan['outputs'][phase]
    assert not target.exists()
    if phase != 'shared-plugin':
        dependency = 'shared-plugin' if phase.endswith('-lower') else phase.rsplit('-', 1)[0] + ('-lower' if phase.endswith('-opt') else '-opt')
        owner = json.loads((W / dependency / 'owner.json').read_text())
        assert owner['passed'] and owner['status'] == 'terminal'
    subprocess.run(plan['commands'][phase], check=True)
    result = dict(output=dict(path=str(target.relative_to(R)), **pin(target)), argv=plan['commands'][phase])
    if phase == 'shared-plugin':
        undefined = {line.split()[-1] for line in subprocess.check_output(
            ['nm', '-D', str(target)], text=True).splitlines() if ' U ' in line and '_ZN4llvm' in line}
        exported = {line.split()[-1] for line in subprocess.check_output(
            ['nm', '-D', str(D / 'opt')], text=True).splitlines() if len(line.split()) == 3}
        assert not undefined - exported, sorted(undefined - exported)
    elif phase.endswith('-native'):
        listing = subprocess.check_output(['objdump', '-d', '-M', 'intel', '--no-show-raw-insn', str(target)], text=True)
        ops = Counter(re.findall(r'^\s*[0-9a-f]+:\s+(\S+)', listing, re.M))
        assert not any(op.startswith(('vfm', 'vfnm')) or op in ('(bad)', '.byte') for op in ops)
        result.update(no_FMA=True, decoded_instructions=sum(ops.values()), instruction_counts=dict(ops),
            defined_symbols=subprocess.check_output(['nm', '-S', '--defined-only', str(target)], text=True),
            undefined_symbols=subprocess.check_output(['nm', '-u', str(target)], text=True))
    else:
        ll = target.with_suffix('.ll')
        assert not ll.exists()
        subprocess.run([str(D / 'llvm-dis'), str(target), '-o', str(ll)], check=True)
        text = ll.read_text()
        headers = re.findall(r'^define .*$', text, re.M)
        attrs = dict(re.findall(r'^attributes #(\d+) = (.*)$', text, re.M))
        assert headers and all('strictfp' in attrs[re.search(r'#(\d+)', h)[1]] for h in headers)
        assert not re.search(r'^\s+%[^=]+ = (?:fadd|fsub|fmul|fdiv|frem|fcmp|fpext|fptrunc|sitofp|uitofp|fptosi|fptoui)\b', text, re.M)
        assert not re.search(r'^\s+.*\b(?:fneg|fadd|fsub|fmul|fdiv|fcmp|call)\s+(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b', text, re.M)
        assert not re.search(r'@llvm\.(?:experimental\.constrained\.)?(?:fma|fmuladd)\.', text)
        result.update(ir=dict(path=str(ll.relative_to(R)), **pin(ll)), definition_headers=headers,
            definition_attributes=attrs,
            adapter_calls=len(re.findall(r'\bcall\b[^\n]*@sirius_scope_v4_[^ (]+\(', text)),
            core_calls=Counter(re.findall(r'\bcall\b[^\n]*@(sirius_shared_value_[^ (]+)\(', text)),
            core_definitions=[h for h in headers if '@sirius_shared_value_' in h],
            core_declarations=re.findall(r'^declare .*@sirius_shared_value_.*$', text, re.M),
            payload_allocas=len(re.findall(r'alloca \[120 x i32\]', text)),
            tag_allocas=len(re.findall(r'alloca \[24 x i32\]', text)),
            constrained_calls=Counter(re.findall(r'\bcall\s+.*?@(llvm\.experimental\.constrained\.[^ (]+)\(', text)))
        signatures = {}
        for header in result['core_definitions'] + result['core_declarations']:
            name = re.search(r'@(sirius_shared_value_[^ (]+)\(', header)[1]
            tagged = '_tagged_' in name
            assert ' hidden ' in header and ' fastcc ' in header, header
            before, arguments = header.split('@' + name + '(', 1)
            arguments = arguments.split(')', 1)[0]
            assert before.count('<8 x i32>') == (6 if tagged else 5), header
            assert arguments.count('<8 x i32>') == (12 if tagged else 10), header
            parts = arguments.split(',')
            assert len(parts) == (13 if tagged else 11) and re.match(r'\s*i32(?:\s|$)', parts[-1]), header
            assert not re.search(r'\b(?:sret|byval|inreg|nest|inalloca|signext|zeroext|nonnull)\b', arguments), header
            signatures[name] = dict(vector_inputs=12 if tagged else 10, mask='i32',
                                    vector_returns=6 if tagged else 5, calling_convention='fastcc',
                                    linkage='hidden external', header=header)
        result['fixed_value_signatures'] = signatures
    write(W / (phase + '-result.json'), result)
    print(json.dumps({k: v for k, v in result.items() if k in ('output', 'ir', 'adapter_calls', 'payload_allocas', 'tag_allocas', 'no_FMA', 'decoded_instructions')}))
    if phase.endswith('-opt'):
        assert result['payload_allocas'] == result['tag_allocas'] == result['adapter_calls'] == 0
        if phase == 'shared-math-opt':
            assert len(result['core_definitions']) == 9 and not result['core_declarations']
            assert set(result['fixed_value_signatures']) == set(CORE_NAMES)
        else:
            assert not result['core_definitions'] and result['core_declarations']
