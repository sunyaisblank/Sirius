"""One private value-interface implementation, before any new runtime gate."""
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
outputs = {'value-plugin': 'value_abi.so', 'value-lower': 'dp-value.bc',
           'value-optimize': 'dp-value.O2.bc'}
if phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    owner = json.loads((W / 'optimize/owner.json').read_text())
    assert owner['passed'] and owner['status'] == 'terminal'
    previous = json.loads((W / 'link-plan.json').read_text())
    config = shlex.split(subprocess.check_output([str(D / 'llvm-config'), '--cxxflags'], text=True))
    commands = {
        'value-plugin': [str(D / 'clang++'), *config, '-O2', '-fPIC', '-shared',
                         str(W / 'value_abi.cpp'), str(B / 'sha256-only.a'),
                         '-Wl,--exclude-libs,ALL', '-o', str(B / outputs['value-plugin'])],
        'value-lower': [str(D / 'opt'), '-load-pass-plugin=' + str(B / outputs['value-plugin']),
                        '-passes=sirius-value-dp,verify', '-verify-each',
                        str(B / 'dp-linked.bc'), '-o', str(B / outputs['value-lower'])],
        'value-optimize': [str(D / 'opt'), '-passes=internalize,default<O2>,function(sroa),globaldce,verify',
                           '-verify-each', '-internalize-public-api-list=' + ','.join(previous['exports']),
                           '-mcpu=skylake-avx512', '-mattr=-fma', str(B / outputs['value-lower']),
                           '-o', str(B / outputs['value-optimize'])],
    }
    paths = [*W.glob('*.py'), W / 'value_abi.cpp', W / 'value_abi_contract.json',
             W / 'scope_contract.json', W / 'private_math.cpp', W / 'link-result.json',
             W / 'link-plan.json', W / 'optimize-result.json', W / 'toolchain.json',
             B / 'dp-linked.bc', B / 'sha256-only.a',
             R / 'out/native-cpu-exact-product-zero-residual/retained8.cpp',
             R / 'out/native-cpu-inline-private-bank/packinghelper.isph',
             R / 'bin/linux-gcc/native-cpu-inline-private-bank/ispc-source/dopri_phase.ispc',
             R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
             R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py']
    for name in ('clang++', 'opt', 'llvm-dis', 'llvm-config'):
        row = json.loads((W / 'toolchain.json').read_text())['tools'][name]
        path = R / row['resolved']
        assert pin(path) == {k: row[k] for k in ('bytes', 'sha256')}
        paths.append(path)
    inputs = {str(p.relative_to(R)): pin(p) for p in paths}
    write(W / 'value-plan.json', dict(commands=commands, outputs=outputs, inputs=inputs,
        exports=previous['exports'],
        hypothesis='A shared value core per arithmetic body plus tiny inlined adapters and initialized-local masked-memory normalization removes the observed packed-bank escape. No blanket arithmetic-body duplication or new workload.'))
    actions = json.loads((W / 'actions.json').read_text())
    dependencies = {'value-plugin': [], 'value-lower': ['value_abi.so'], 'value-optimize': ['dp-value.bc']}
    for name in commands:
        assert name not in actions and not (W / name).exists()
        actions[name] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
            kind='diagnostic', expected_exit=0, software_provider=False,
            limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
            argv=['/usr/bin/python3', '-B', str(W / 'value.py'), name],
            frozen_inputs=sorted(inputs) + ['out/native-cpu-strict-ir-linking/value-plan.json'] +
                          [str((B / n).relative_to(R)) for n in dependencies[name]],
            scope='Exact-module private value-interface compiler/structural gate only. One shared strict arithmetic body per operation; no numerical fixture, provider rebuild, original renderer or speed/adoption claim. Original public/scalar/fallback implementations unchanged.')
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=list(commands), compiler_run=False)))
else:
    plan = json.loads((W / 'value-plan.json').read_text())
    assert phase in plan['commands']
    for name, expected in plan['inputs'].items():
        assert pin(R / name) == expected, name
    prior = {'value-lower': 'value-plugin', 'value-optimize': 'value-lower'}.get(phase)
    if prior:
        owner = json.loads((W / prior / 'owner.json').read_text())
        assert owner['passed'] and owner['status'] == 'terminal'
        output = json.loads((W / (prior + '-result.json')).read_text())['output']
        assert pin(R / output['path']) == {k: output[k] for k in ('bytes', 'sha256')}
    output = B / plan['outputs'][phase]
    assert not output.exists()
    subprocess.run(plan['commands'][phase], check=True)
    result = dict(phase=phase, argv=plan['commands'][phase],
                  output=dict(path=str(output.relative_to(R)), **pin(output)))
    if phase == 'value-plugin':
        undefined = {line.split()[-1] for line in subprocess.check_output(
            ['nm', '-D', str(output)], text=True).splitlines() if ' U ' in line and '_ZN4llvm' in line}
        exported = {line.split()[-1] for line in subprocess.check_output(
            ['nm', '-D', str(D / 'opt')], text=True).splitlines() if len(line.split()) == 3}
        assert not undefined - exported, sorted(undefined - exported)
    else:
        text_path = output.with_suffix('.ll')
        assert not text_path.exists()
        subprocess.run([str(D / 'llvm-dis'), str(output), '-o', str(text_path)], check=True)
        text = text_path.read_text()
        headers = re.findall(r'^define .*$', text, re.M)
        attrs = dict(re.findall(r'^attributes #(\d+) = (.*)$', text, re.M))
        assert headers and all('strictfp' in attrs[re.search(r'#(\d+)', row)[1]] for row in headers)
        ordinary = Counter(re.findall(r'^\s+%[^=]+ = (fadd|fsub|fmul|fdiv|frem|fcmp|fpext|fptrunc|sitofp|uitofp|fptosi|fptoui)\b', text, re.M))
        assert not ordinary, ordinary
        constrained = Counter(re.findall(r'\bcall\s+.*?@(llvm\.experimental\.constrained\.[^ (]+)\(', text))
        if phase == 'value-lower':
            assert dict(constrained) == json.loads((W / 'link-result.json').read_text())['constrained_calls']
        result.update(ir=dict(path=str(text_path.relative_to(R)), **pin(text_path)),
            definitions=len(headers), definition_headers=headers,
            ordinary_fp_operations=dict(ordinary), constrained_calls=dict(constrained),
            adapter_calls=len(re.findall(r'\bcall\b[^\n]*@sirius_scope_v4_[^ (]+\(', text)),
            core_calls=len(re.findall(r'\bcall\b[^\n]*@sirius_value_core_[^ (]+\(', text)),
            packed_payload_allocas=len(re.findall(r'alloca \[120 x i32\]', text)),
            packed_tag_allocas=len(re.findall(r'alloca \[24 x i32\]', text)),
            core_definitions=sum('@sirius_value_core_' in h for h in headers),
            structural_gate='Adoption requires elimination of original480/96-byte packed banks and pointer adapters. IR aggregate values may still spill/use stack at machine level. No speed or numerical claim.')
    write(W / (phase + '-result.json'), result)
    print(json.dumps({k: v for k, v in result.items() if k not in ('argv', 'definition_headers', 'constrained_calls')}))
