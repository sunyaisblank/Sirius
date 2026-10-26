"""Derive and build one shared, bounded Endpoint program admission extension."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-endpoint-header-domain'


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def doc(path):
    return json.loads(path.read_text())


def write(path, value):
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2) + '\n')


phase, = sys.argv[1:]
if phase == 'derive':
    sources = {
        'native_packet.hpp': R / 'out/native-cpu-active-prefix/native_packet.hpp',
        'admission.cpp': R / 'out/native-cpu-dispatch-admission/admission.cpp',
        'endpoint.cpp': R / 'bin/linux-gcc/native-cpu-scalar-zero/typed-batch/endpoint.cpp',
    }
    original = {name: path.read_text() for name, path in sources.items()}
    header = original['native_packet.hpp']
    assert header.count('inline bool admit_packet(') == 1
    header = header.replace('#include <limits>', '#include <limits>\n#include <endpoint.program.hpp>', 1)
    header = header.replace('inline bool admit_packet(',
        'template<bool endpoint_header_domain = false>\ninline bool admit_packet(', 1)
    old = '    return std::memcmp(suffix, program, program_words * sizeof(std::uint32_t)) == 0;'
    assert header.count(old) == 1
    new = '''    if constexpr (endpoint_header_domain) {
        // The unchanged shader compares this word only for the specialized
        // zero-spin/charge/lambda family, refusing that family on a mismatch.
        // It must still execute: internal and public output words are observable.
        constexpr std::size_t specialized_base = 11240;
        constexpr std::size_t register_word = specialized_base + 1;
        static_assert(sizeof(endpoint_program) / sizeof(std::uint32_t) == 21143);
        static_assert(endpoint_program[1] == 612 &&
                      endpoint_program[register_word] == endpoint_program[1]);
        static_assert(endpoint_program[specialized_base] == 1943 &&
                      endpoint_program[specialized_base + 2] == 100);
        if (input_row_words != 225 || output_row_words != 3829 ||
            program != endpoint_program || program_words != 21143)
            return false;
        return std::memcmp(suffix, program, register_word * sizeof(std::uint32_t)) == 0 &&
               std::memcmp(suffix + (register_word + 1) * sizeof(std::uint32_t),
                           program + register_word + 1,
                           (program_words - register_word - 1) * sizeof(std::uint32_t)) == 0;
    }
''' + old
    header = header.replace(old, new, 1)
    admission = original['admission.cpp']
    assert admission.count('#include "../native-cpu-active-prefix/native_packet.hpp"') == 1
    admission = admission.replace('#include "../native-cpu-active-prefix/native_packet.hpp"',
                                  '#include "native_packet.hpp"', 1)
    old = 'sirius_cpu_bridge::admit_packet(ADMISSION_VALUES, 225U, 3829U,'
    assert admission.count(old) == 1
    admission = admission.replace(old, old.replace('admit_packet(', 'admit_packet<true>('), 1)
    scalar = original['endpoint.cpp']
    assert scalar.count('sirius_cpu_bridge::admit_packet(') == 1
    scalar = scalar.replace('sirius_cpu_bridge::admit_packet(', 'sirius_cpu_bridge::admit_packet<true>(', 1)
    for name, text in [('native_packet.hpp', header), ('admission.cpp', admission), ('endpoint.cpp', scalar)]:
        with (W / name).open('x') as stream:
            stream.write(text)
    basis = {str(path.relative_to(R)): pin(path) for path in sources.values()}
    for name in ('src/sirius/kernels/retained_endpoint.slang',
                 'bin/linux-gcc/native-cpu-scalar-zero/typed-batch/endpoint.native.hpp',
                 'bin/linux-gcc/native-cpu-slang-bridge/typed-batch/endpoint.program.hpp',
                 'out/native-cpu-dispatch-admission/batch-compile-result.json',
                 'out/native-cpu-scalar-zero/cpp-compile-result.json'):
        basis[name] = pin(R / name)
    write(W / 'derivation.json', dict(basis=basis,
        derived={name: pin(W / name) for name in sources},
        admitted_program_difference=dict(stage='Endpoint', relative_word=11241, bytes=4,
            canonical=612, allowed='Any uint32; the unchanged shader owns acceptance/refusal'),
        scope='One shared preflight/scalar-worker predicate. All original bounds/alignment/alias/grid/capacity checks remain. All other program bytes and Transport/DopriPhase stay strict; no shader, row/field/count/layout, arithmetic, FP, scheduler, input or output rewrite. Source premise does not identify the prior timeout cost.'))
elif phase == 'prepare-build':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    proof = doc(W / 'derivation.json')
    for name, expected in proof['basis'].items():
        assert pin(R / name) == expected, name
    for name, expected in proof['derived'].items():
        assert pin(W / name) == expected, name
    assert not B.exists()
    records = [doc(R / 'out/native-cpu-dispatch-admission/batch-compile-result.json')['objects'][0],
               next(row for row in doc(R / 'out/native-cpu-scalar-zero/cpp-compile-result.json')['objects']
                    if row['stage'] == 'endpoint')]
    commands = []
    for name, record in zip(('admission', 'endpoint'), records):
        args = list(record['argv'])
        args[args.index('-c') + 1] = str(W / (name + '.cpp'))
        args[args.index('-o') + 1] = str(B / (name + '.o'))
        commands.append(args)
    inputs = set(proof['basis']) | {str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
    for directory in ('bin/linux-gcc/native-cpu-scalar-zero/typed-batch',
                      'bin/linux-gcc/native-cpu-slang-bridge/typed-batch',
                      'bin/linux-gcc/native-cpu-static-program/typed-batch',
                      'out/native-cpu-slang-bridge', '/opt/slang/include'):
        inputs |= {str(p.relative_to(R)) if p.is_relative_to(R) else str(p)
                   for p in (R / directory).rglob('*') if p.is_file()}
    inputs |= {'/usr/bin/x86_64-linux-gnu-g++-14', '/usr/bin/ld'}
    write(W / 'build-plan.json', dict(commands=commands))
    inputs.add(str((W / 'build-plan.json').relative_to(R)))
    write(W / 'actions.json', {'build': dict(
        revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='build', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=['/usr/bin/python3', '-B', str(W / 'component.py'), 'build'],
        frozen_inputs=sorted(inputs), scope='Two baseline strict-FP C++ compilations for the same Endpoint admission extension: pure admission and unchanged scalar shader wrapper. No provider, workload or qualification claim.')} )
elif phase == 'build':
    B.mkdir(exist_ok=False)
    result = []
    for args in doc(W / 'build-plan.json')['commands']:
        subprocess.run(args, check=True)
        path = Path(args[args.index('-o') + 1])
        result.append(dict(argv=args, exit=0, object=dict(path=str(path.relative_to(R)), **pin(path))))
    write(W / 'build-result.json', dict(objects=result, new_compilations=2))
    print(json.dumps(dict(new_compilations=2, passed=True)))
else:
    raise ValueError(phase)
