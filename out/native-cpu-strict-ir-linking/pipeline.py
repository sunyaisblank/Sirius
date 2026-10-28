"""Bounded exact-module implementation steps; no renderer or workload retries."""
from pathlib import Path
import hashlib
import json
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
fixed = phase.endswith('-sha256')
plan_path = W / ('pipeline-sha256-plan.json' if fixed else 'pipeline-plan.json')
if phase in ('prepare', 'prepare-sha256'):
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    verified = json.loads((W / 'verify/owner.json').read_text())
    assert verified['passed'] and verified['status'] == 'terminal'
    original = json.loads((W / 'inventory.json').read_text())['output']
    assert pin(R / original['path']) == {k: original[k] for k in ('bytes', 'sha256')}
    # opt exports its own LLVM implementation. Do not statically link another
    # copy into the plugin (duplicate registries/contexts would be incorrect).
    config = subprocess.check_output([str(D / 'llvm-config'), '--cxxflags'], text=True)
    assert '-D_GLIBCXX_USE_CXX11_ABI=1' in config and '-fno-rtti' in config
    commands = {
        'plugin': [str(D / 'clang++'), *shlex.split(config), '-O2', '-fPIC', '-shared',
                   str(W / 'constrain_scope.cpp'), '-o', str(B / 'constrain_scope.so')],
        'convert': [str(D / 'opt'), '-load-pass-plugin=' + str(B / 'constrain_scope.so'),
                    '-passes=sirius-constrain-dp,verify', '-verify-each',
                    str(B / 'dopri_phase.ll'), '-o', str(B / 'dopri_phase.strict.bc')],
    }
    if fixed:
        failed = json.loads((W / 'convert/owner.json').read_text())
        assert failed['status'] == 'terminal' and not failed['passed']
        assert not failed['remaining_owned_processes']
        # The first load exposed exactly three unexported SHA methods. Only
        # SHA256.cpp.o is imported, never the whole Support/Core registries.
        commands['plugin'][-2:-2] = [str(B / 'sha256-only.a'), '-Wl,--exclude-libs,ALL']
        commands['plugin'][-1] = str(B / 'constrain_scope.sha256.so')
        commands['convert'][1] = '-load-pass-plugin=' + str(B / 'constrain_scope.sha256.so')
        commands['convert'][-1] = str(B / 'dopri_phase.strict.sha256.bc')
        commands = {k + '-sha256': v for k, v in commands.items()}
    sources = [*W.glob('*.py'), W / 'constrain_scope.cpp', W / 'scope_contract.json',
               W / 'toolchain.json', W / 'inventory.json', W / 'plan.json',
               R / 'out/native-cpu-exact-product-zero-residual/retained8.cpp',
               R / 'out/native-cpu-exact-product-zero-residual/optional.cpp',
               R / original['path'],
               R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
               R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py']
    sources += [R / n for n in json.loads((W / 'plan.json').read_text())['source_inputs']]
    if fixed:
        sources += [B / 'toolchain/lib/libLLVMSupport.a', Path('/usr/bin/ar'),
                    W / 'convert/owner.json', W / 'convert/stderr.log']
    for name in ('clang++', 'opt', 'llvm-config'):
        tool = json.loads((W / 'toolchain.json').read_text())['tools'][name]
        assert pin(R / tool['resolved']) == {k: tool[k] for k in ('bytes', 'sha256')}
        sources.append((D / name).resolve())
    inputs = {str(p.relative_to(R)) if p.is_relative_to(R) else str(p): pin(p) for p in sources}
    write(plan_path, dict(commands=commands, inputs=inputs,
        plugin_linkage='Resolve LLVM symbols from exact SDK opt executable; matching libstdc++ C++11 ABI, RTTI/EH off. No duplicate static LLVM.',
        outputs=({'plugin-sha256': 'constrain_scope.sha256.so',
                  'convert-sha256': 'dopri_phase.strict.sha256.bc'} if fixed else
                 dict(plugin='constrain_scope.so', convert='dopri_phase.strict.bc'))))
    actions = json.loads((W / 'actions.json').read_text())
    for name in commands:
        assert name not in actions and not (W / name).exists()
        plugin_name = 'constrain_scope.sha256.so' if fixed else 'constrain_scope.so'
        extra = [str((B / plugin_name).relative_to(R))] if name.startswith('convert') else []
        actions[name] = dict(
            revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
            kind='diagnostic', expected_exit=0, software_provider=False,
            limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
            argv=['/usr/bin/python3', '-B', str(W / 'pipeline.py'), name],
            frozen_inputs=sorted(inputs) + [str(plan_path.relative_to(R))] + extra,
            scope='Build and apply the reviewed exact-DP strict conversion to a separate hidden copy. Full input pin/fingerprint and LLVM verifier; all cost included. No link/inlining, runtime, numerical, speed or adoption claim. Old public and fallback objects remain untouched.')
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=list(commands), compiler_run=False)))
else:
    plan = json.loads(plan_path.read_text())
    assert phase in plan['commands']
    for name, expected in plan['inputs'].items():
        assert pin(R / name) == expected, name
    if phase.startswith('convert'):
        predecessor = 'plugin-sha256' if fixed else 'plugin'
        owner = json.loads((W / predecessor / 'owner.json').read_text())
        assert owner['passed'] and owner['status'] == 'terminal'
        made = json.loads((W / (predecessor + '-result.json')).read_text())['output']
        assert pin(R / made['path']) == {k: made[k] for k in ('bytes', 'sha256')}
    output = B / plan['outputs'][phase]
    assert not output.exists()
    support = None
    if phase == 'plugin-sha256':
        member = subprocess.check_output(['/usr/bin/ar', 'p', str(B / 'toolchain/lib/libLLVMSupport.a'), 'SHA256.cpp.o'])
        single = B / 'SHA256.cpp.o'
        with single.open('xb') as stream:
            stream.write(member)
        subprocess.run(['/usr/bin/ar', 'cr', str(B / 'sha256-only.a'), str(single)], check=True)
        support = dict(member='SHA256.cpp.o', **pin(single))
    subprocess.run(plan['commands'][phase], check=True)
    if phase.startswith('plugin'):
        undefined = {line.split()[-1] for line in subprocess.check_output(
            ['nm', '-D', str(output)], text=True).splitlines() if ' U ' in line and '_ZN4llvm' in line}
        exported = {line.split()[-1] for line in subprocess.check_output(
            ['nm', '-D', str(D / 'opt')], text=True).splitlines() if len(line.split()) == 3}
        if fixed:
            assert not undefined - exported, sorted(undefined - exported)
    result = dict(phase=phase, argv=plan['commands'][phase],
                  output=dict(path=str(output.relative_to(R)), **pin(output)),
                  isolated_support_member=support,
                  scope='Compiler implementation gate only; no runtime/numerical/speed/adoption qualification.')
    write(W / (phase + '-result.json'), result)
    print(json.dumps(result))
