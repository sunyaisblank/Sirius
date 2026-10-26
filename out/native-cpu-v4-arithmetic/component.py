"""Build one separately guarded optional ISA variant of exact original math."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-v4-arithmetic'
OLD = R / 'out/native-cpu-value-fact-stores'


def doc(p):
    return json.loads(p.read_text())


def pin(p):
    with p.open('rb') as f:
        return dict(bytes=p.stat().st_size, sha256=hashlib.file_digest(f, 'sha256').hexdigest())


def write(p, x):
    with p.open('x') as f:
        f.write(json.dumps(x, indent=2) + '\n')


def closure():
    x = doc(OLD / 'prepared-cpp.json')
    for n, e in x['inputs'].items():
        assert pin(R / n) == e, n
    built = doc(OLD / 'cpp-build-result.json')
    assert built['passed']
    for row in built['compilations']:
        assert pin(R / row['object']['path']) == {k: row['object'][k] for k in ('bytes', 'sha256')}
    return x, built


phase, = sys.argv[1:]
if phase == 'derive':
    prepared, built = closure()
    source = (OLD / 'retained8.cpp').read_text()
    signatures = re.findall(r'void (sirius_native_(?:dense|tagged)_\w+8)\((.*?)\) noexcept', source, re.S)
    assert len(signatures) == 10 and len({n for n, _ in signatures}) == 10
    names = [n for n, _ in signatures]
    optional = '// Exact original math; only external names and guarded target differ.\n'
    optional += ''.join('#define ' + n + ' ' + n.replace('sirius_native_', 'sirius_native_v4_') + '\n' for n in names)
    optional += '#include "../native-cpu-value-fact-stores/retained8.cpp"\n'
    dispatch = '''// Compile without optional ISA flags. The provider retains its original v3 guard.
// No vector or Triple crosses this original pointer/u32 ABI.
#include <cstdint>
extern "C" __attribute__((visibility("hidden")))
bool sirius_native_v4_supported() noexcept {
    return __builtin_cpu_supports("x86-64-v4");
}
'''
    for name, parameters in signatures:
        assert len(re.findall('std::uint32_t', parameters)) in (5, 6)
        v3 = name.replace('sirius_native_', 'sirius_native_v3_')
        v4 = name.replace('sirius_native_', 'sirius_native_v4_')
        arguments = ','.join(part.strip().split()[-1].lstrip('*') for part in parameters.split(','))
        for n in (v3, v4):
            dispatch += 'extern "C" __attribute__((visibility("hidden"))) void ' + n + '(' + parameters + ') noexcept;\n'
        dispatch += 'extern "C" __attribute__((visibility("hidden"))) void ' + name + '(' + parameters + ') noexcept {\n'
        dispatch += '    if(sirius_native_v4_supported()) ' + v4 + '(' + arguments + ');\n'
        dispatch += '    else ' + v3 + '(' + arguments + ');\n}\n'
    for name, data in [('optional.cpp', optional), ('dispatch.cpp', dispatch)]:
        with (W / name).open('x') as f:
            f.write(data)
    write(W / 'derivation.json', dict(original_source=pin(OLD / 'retained8.cpp'),
        exports=names, optional=pin(W / 'optional.cpp'), dispatch=pin(W / 'dispatch.cpp'),
        scope='Same original five-word/eight-lane arithmetic source and FP operation order. '
              'Optional v4 object behind baseline-compiled OS-aware builtin guard; original v3 '
              'object, provider guard, raw stages, scalar tails, host, SPIR-V, public ABI, '
              'allocations/layouts/counts/precision/FP policy remain. '
              'No IFUNC, LTO, global initializer, buffer, setting or acceptance claim.'))
elif phase in ('build', 'link', 'run'):
    prepared, built = closure()
    revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    proof = doc(W / 'derivation.json')
    assert pin(W / 'optional.cpp') == proof['optional'] and pin(W / 'dispatch.cpp') == proof['dispatch']
    inputs = set(prepared['inputs']) | {str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
    inputs |= {row['object']['path'] for row in built['compilations']}
    if phase == 'build':
        assert not B.exists()
        argv = ['/usr/bin/python3', '-B', str(W / 'component.py'), 'compile']
        scope = 'One optional v4 exact-math compile plus baseline-only CPU dispatcher; no workload or ISA sweep.'
    elif phase == 'link':
        owner = doc(W / 'build/owner.json')
        assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
        inputs |= {str((B / n).relative_to(R)) for n in ('optional.o', 'dispatch.o', 'v3.o')}
        driver = next(row['object']['path'] for row in built['compilations'] if row['name'] == 'driver')
        argv = [prepared['link'][0], str(R / driver), str(B / 'dispatch.o'), str(B / 'v3.o'), str(B / 'optional.o'), '-o', str(B / 'abi-check')]
        scope = 'Link unchanged original tagged arithmetic driver with original v3, optional v4 and baseline dispatcher.'
    else:
        owner = doc(W / 'link/owner.json')
        assert owner['passed'] and owner['source_unchanged'] and owner['whole_inputs_unchanged']
        inputs.add(str((B / 'abi-check').relative_to(R)))
        argv = [str(B / 'abi-check')]
        scope = 'One original555660-call tagged arithmetic compatibility fixture on actual guarded target; exact bank/tag/mask/alias/refusal/FP-controls checks. No math-oracle, provider, benefit, frame or scientific/release admission.'
    actions = doc(W / 'actions.json') if (W / 'actions.json').exists() else {}
    assert phase not in actions
    actions[phase] = dict(revision=revision, kind='diagnostic', expected_exit=0,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=argv, frozen_inputs=sorted(inputs), software_provider=False, scope=scope)
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=phase, executed=False)))
elif phase == 'compile':
    prepared, built = closure()
    B.mkdir(exist_ok=False)
    flags = prepared['compilations'][0]['argv'][:prepared['compilations'][0]['argv'].index('-c')]
    assert flags.count('-march=x86-64-v3') == 1 and '-ffp-contract=off' in flags
    commands = []
    for name in ('optional', 'dispatch'):
        args = list(flags)
        args[args.index('-march=x86-64-v3')] = '-march=x86-64-v4' if name == 'optional' else '-march=x86-64'
        if name == 'optional':
            args += ['-mprefer-vector-width=256']
        args += ['-c', str(W / (name + '.cpp')), '-o', str(B / (name + '.o'))]
        code = subprocess.run(args).returncode
        commands.append(dict(argv=args, exit=code))
        assert code == 0, commands
    old = next(row['object']['path'] for row in built['compilations'] if row['name'] == 'retained8')
    args = ['objcopy']
    for name in doc(W / 'derivation.json')['exports']:
        args += ['--redefine-sym', name + '=' + name.replace('sirius_native_', 'sirius_native_v3_')]
    args += [str(R / old), str(B / 'v3.o')]
    subprocess.run(args, check=True)
    write(W / 'compile-result.json', dict(commands=commands, rename_argv=args,
        objects={n:pin(B / (n + '.o')) for n in ('optional', 'dispatch', 'v3')},
        original_v3=next(row['object'] for row in built['compilations'] if row['name'] == 'retained8'),
        driver_recompiled=False, original_v3_recompiled=False))
    print(json.dumps(dict(optional_compilations=1, baseline_dispatcher_compilations=1, original_v3_recompiled=False)))
else:
    raise ValueError(phase)
