"""Compile/link strict private math, then inspect one ordinary O2 pipeline."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import re
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
outputs = dict(math='private-math.bc', link='dp-linked.bc', optimize='dp-linked.O2.bc')
if phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    owner = json.loads((W / 'convert-sha256/owner.json').read_text())
    assert owner['passed'] and owner['status'] == 'terminal'
    original = (B / 'dopri_phase.ll').read_text()
    exports = re.findall(r'^define (?!internal).*?@(sirius_simd_v4_dopri_phase_packet[^ (]*)\(', original, re.M)
    assert len(exports) == 2
    exports = [n.replace('sirius_simd_v4_dopri_phase_packet', 'sirius_scope_dp_packet', 1) for n in exports]
    commands = {
        'math': [str(D / 'clang++'), '--target=x86_64-unknown-linux-gnu',
                 '--gcc-install-dir=/usr/lib/gcc/x86_64-linux-gnu/14', '-stdlib=libstdc++',
                 '-std=c++2c', '-O2', '-fPIC', '-DNDEBUG', '-march=x86-64-v4',
                 '-mprefer-vector-width=256', '-fno-fast-math', '-ffp-model=strict',
                 '-frounding-math', '-ffp-exception-behavior=strict', '-ffp-contract=off',
                 '-fdenormal-fp-math=ieee', '-ffp-eval-method=source',
                 '-Xclang', '-fdenormal-fp-math-f32=ieee', '-emit-llvm',
                 '-c', str(W / 'private_math.cpp'), '-o', str(B / outputs['math'])],
        'link': [str(D / 'llvm-link'), str(B / 'dopri_phase.strict.sha256.bc'),
                 str(B / outputs['math']), '-o', str(B / outputs['link'])],
        'optimize': [str(D / 'opt'), '-passes=internalize,default<O2>,verify', '-verify-each',
                     '-internalize-public-api-list=' + ','.join(exports),
                     '-mcpu=skylake-avx512', str(B / outputs['link']),
                     '-o', str(B / outputs['optimize'])],
    }
    paths = [*W.glob('*.py'), W / 'private_math.cpp', W / 'scope_contract.json',
             W / 'constrain_scope.cpp', W / 'convert-sha256-result.json',
             W / 'toolchain.json', B / 'dopri_phase.strict.sha256.bc',
             R / 'out/native-cpu-exact-product-zero-residual/retained8.cpp',
             R / 'bin/linux-gcc/native-cpu-exact-product-zero-residual/optional.o',
             R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
             R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py']
    for name in ('clang++', 'llvm-link', 'opt', 'llvm-dis'):
        row = json.loads((W / 'toolchain.json').read_text())['tools'][name]
        path = R / row['resolved']
        assert pin(path) == {k: row[k] for k in ('bytes', 'sha256')}
        paths.append(path)
    inputs = {str(p.relative_to(R)): pin(p) for p in paths}
    write(W / 'link-plan.json', dict(commands=commands, inputs=inputs, outputs=outputs,
        exports=exports, forced_inlining=False,
        scope='One strict math compilation and linked ordinary O2 implementation. Internalization permits constant-index propagation and ordinary inlining. Preserve slow_* noinline and all arithmetic; no runtime or speed study.'))
    actions = json.loads((W / 'actions.json').read_text())
    dependencies = dict(math=[], link=['private-math.bc'], optimize=['dp-linked.bc'])
    for name in commands:
        assert name not in actions and not (W / name).exists()
        actions[name] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
            kind='diagnostic', expected_exit=0, software_provider=False,
            limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
            argv=['/usr/bin/python3', '-B', str(W / 'link.py'), name],
            frozen_inputs=sorted(inputs) + ['out/native-cpu-strict-ir-linking/link-plan.json'] +
                          [str((B / n).relative_to(R)) for n in dependencies[name]],
            scope='Strict private C++/ISPC linking and structural optimization gate only. No original workload, numerical compatibility, speed or provider adoption established. Old public/fallback objects remain exact.')
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=list(commands), compiler_run=False)))
else:
    plan = json.loads((W / 'link-plan.json').read_text())
    assert phase in plan['commands']
    for name, expected in plan['inputs'].items():
        assert pin(R / name) == expected, name
    if phase != 'math':
        prior = 'math' if phase == 'link' else 'link'
        owner = json.loads((W / prior / 'owner.json').read_text())
        assert owner['passed'] and owner['status'] == 'terminal'
        output = json.loads((W / (prior + '-result.json')).read_text())['output']
        assert pin(R / output['path']) == {k: output[k] for k in ('bytes', 'sha256')}
    output = B / plan['outputs'][phase]
    assert not output.exists()
    subprocess.run(plan['commands'][phase], check=True)
    text_path = output.with_suffix('.ll')
    assert not text_path.exists()
    subprocess.run([str(D / 'llvm-dis'), str(output), '-o', str(text_path)], check=True)
    text = text_path.read_text()
    definitions = re.findall(r'^define .*$', text, re.M)
    attrs = dict(re.findall(r'^attributes #(\d+) = (.*)$', text, re.M))
    assert definitions and all('strictfp' in attrs[re.search(r'#(\d+)', row)[1]] for row in definitions)
    flags = re.findall(r'^\s+.*\b(?:fadd|fsub|fmul|fdiv|fcmp|call)\s+(?:fast|reassoc|nnan|ninf|nsz|arcp|contract|afn)\b', text, re.M)
    assert not flags, flags[:3]
    ordinary = Counter(re.findall(r'^\s+%[^=]+ = (fadd|fsub|fmul|fdiv|frem|fcmp|fpext|fptrunc|sitofp|uitofp|fptosi|fptoui)\b', text, re.M))
    assert not ordinary, ordinary
    constrained = Counter(re.findall(r'\bcall\s+.*?@(llvm\.experimental\.constrained\.[^ (]+)\(', text))
    math_calls = Counter(re.findall(r'\bcall\s+.*?@(sirius_scope_v4_[^ (]+)\(', text))
    result = dict(phase=phase, argv=plan['commands'][phase],
        output=dict(path=str(output.relative_to(R)), **pin(output)),
        ir=dict(path=str(text_path.relative_to(R)), **pin(text_path)), definitions=len(definitions),
        ordinary_fp_operations=dict(ordinary), constrained_calls=dict(constrained),
        private_math_calls=dict(math_calls), definition_headers=definitions,
        scope='Static compiler gate only. Strict definitions/no ordinary FP/no fast flags checked; numerical/primitive retention/predicate and emitted ISA validation still required before runtime adoption.')
    write(W / (phase + '-result.json'), result)
    print(json.dumps({k: result[k] for k in ('phase', 'output', 'definitions', 'private_math_calls')}))
