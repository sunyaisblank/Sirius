"""Bounded raw-stage compilation and reuse of two sealed complete-bank fixtures.

This is a component check. No provider is changed or renderer launched.
"""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-v4-raw-packets'
O = B / 'objects'
sys.dont_write_bytecode = True
sys.path.insert(0, str(R / 'out/native-cpu-field-compatibility'))
from elf_proof import Elf, associate


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size,
                    sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def named(path):
    return dict(path=str(path.relative_to(R)), **pin(path))


def write(path, value):
    with path.open('x') as stream:
        stream.write(json.dumps(value, indent=2) + '\n')


def owner(phase):
    value = doc(W / phase / 'owner.json')
    assert value['status'] == 'terminal' and value['passed']
    assert value['source_unchanged'] and value['whole_inputs_unchanged']
    assert not value['remaining_owned_processes'] and not value['cleanup_errors']
    return value


def verify_inputs():
    proof = doc(W / 'derivation.json')
    for key in ('source_inputs', 'derived'):
        for name, expected in proof[key].items():
            assert pin(R / name) == expected, name
    for row in proof['reused_objects']:
        assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}
    return proof


def instructions(path):
    text = subprocess.check_output(['objdump', '-d', '-M', 'intel',
        '--no-show-raw-insn', str(path)], text=True)
    result = []
    for line in text.splitlines():
        match = re.match(r'^\s*[0-9a-f]+:\s+(\S+)(?:\s+(.*))?$', line)
        if match:
            result.append((match[1], match[2] or ''))
    assert result and not any(op in ('(bad)', '.byte') for op, _ in result)
    return result


def rename_reference(original, output, old, new, recovered=False):
    if recovered:
        assert output.exists()
    else:
        assert not output.exists()
        subprocess.run(['objcopy', '--redefine-sym', old + '=' + new,
                        str(original), str(output)], check=True)
    a, b = Elf(original), Elf(output)
    inverse = {new: old}
    assert any(not s['section'] and s['name'] == old for s in a.symbols)
    assert not any(s['name'] == new for s in a.symbols)
    assert sorted(a.descriptor(s, {}) for s in a.symbols) == sorted(
        b.descriptor(s, inverse) for s in b.symbols)
    assert a.names == b.names
    for i, section in enumerate(a.sections):
        if section[2] & 2:
            # File positions can move; allocated payloads and their runtime
            # properties, symbols and relocation semantics must remain exact.
            assert tuple(section[k] for k in (1, 2, 3, 5, 8, 9)) == tuple(
                b.sections[i][k] for k in (1, 2, 3, 5, 8, 9)), a.names[i]
            assert a.bytes(i) == b.bytes(i), a.names[i]
    def relocs(elf, names):
        return sorted((elf.section_name(section), offset, kind, addend,
                       elf.descriptor(target, names))
                      for section, offset, kind, addend, target in elf.relocations)
    assert relocs(a, {}) == relocs(b, inverse)
    return dict(original=named(original), renamed=named(output), old=old, new=new,
        all_allocated_sections_identical=True,
        all_symbols_and_relocations_inverse_normalized_equal=True)


phase, = sys.argv[1:]
proof = verify_inputs()
if phase.startswith('prepare-'):
    step = phase.removeprefix('prepare-')
    assert step in ('build', 'finish-build', 'link', 'run')
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    inputs = set(proof['source_inputs']) | set(proof['derived'])
    inputs |= {row['path'] for row in proof['reused_objects']}
    inputs |= {str(p.relative_to(R)) for p in W.iterdir() if p.is_file()}
    for tool in ('/usr/bin/objcopy', '/usr/bin/objdump', '/usr/bin/ld',
                 '/usr/bin/x86_64-linux-gnu-g++-14', '/usr/bin/ldd'):
        inputs.add(str(Path(tool).resolve(strict=True)))
    # These remain the original fixtures; only their undefined raw entry symbol
    # is renamed. The checked v4 constructor is linked into both executables.
    fixture_paths = [
        'out/native-cpu-accumulate-compatibility/prepared.json',
        'out/native-cpu-accumulate-compatibility/compile-completion.json',
        'bin/linux-gcc/native-cpu-accumulate-compatibility/driver.o',
        'bin/linux-gcc/native-cpu-accumulate-compatibility/baseline.o',
        'bin/linux-gcc/native-cpu-endpoint-header-domain/driver.o',
        'out/native-cpu-endpoint-header-domain/driver-build-result.json',
        'out/native-cpu-endpoint-header-provider/provider/source-contract.json',
        'out/native-cpu-field-compatibility/elf_proof.py',
        'out/native-cpu-accumulate-compatibility/elf_proof.py',
        'out/native-cpu-row-ownership/cold-frame/run_owned.py',
        'out/native-cpu-row-ownership/cold-frame/owned_guard.py',
    ]
    inputs |= set(fixture_paths)
    if step != 'build':
        if step == 'finish-build':
            failed = doc(W / 'build/owner.json')
            assert failed['status'] == 'terminal' and not failed['passed']
            assert failed['source_unchanged'] and failed['whole_inputs_unchanged']
            assert not failed['remaining_owned_processes'] and not failed['cleanup_errors']
            inputs.add('attestations/qualification/ed9bf462/native-cpu-v4-raw-fixture-copy-check-failure/MANIFEST.json')
        else:
            owner('finish-build')
            assert doc(W / 'code-proof.json')['passed']
        inputs |= {str(p.relative_to(R)) for p in O.glob('*.o')}
    if step == 'run':
        owner('link')
        inputs |= {str((B / name).relative_to(R)) for name in ('transport-check', 'endpoint-check')}
        inputs |= set(doc(W / 'link-result.json')['runtime_libraries'])
    actions = {} if step == 'build' else doc(W / 'actions.json')
    assert step not in actions and not (W / step).exists()
    actions[step] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='diagnostic', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=['/usr/bin/python3', '-B', str(W / 'qualify.py'), step],
        frozen_inputs=sorted(inputs),
        scope='One raw-v4 component phase; original 8-lane ABI, original AO arithmetic, AP admission and scalar fallback. Two original finite complete-bank fixtures only. No provider, performance, frame, full scientific or release acceptance.')
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=step, executed=False, frozen_inputs=len(inputs))))
elif phase == 'build':
    O.mkdir(exist_ok=False)
    commands = list(proof['raw_compilations'])
    for row in commands:
        assert not (R / row['output']).exists()
        subprocess.run(row['argv'], check=True)
    adapter = list(proof['adapter_compile_argv'])
    subprocess.run(adapter, check=True)
    gate = adapter[:adapter.index('-c')] + ['-c', str(W / 'cpu_gate.cpp'), '-o', str(O / 'cpu_gate.o')]
    subprocess.run(gate, check=True)
    renames = []
    for stage, original in (
        ('transport', R / 'bin/linux-gcc/native-cpu-accumulate-compatibility/driver.o'),
        ('endpoint', R / 'bin/linux-gcc/native-cpu-endpoint-header-domain/driver.o')):
        renames.append(rename_reference(original, O / (stage + '-driver.o'),
            'sirius_simd_' + stage + '_packet', 'sirius_simd_v4_' + stage + '_packet'))
    verify_inputs()
    write(W / 'build-result.json', dict(passed=True, raw_compilations=commands,
        adapter_compile_argv=adapter, gate_compile_argv=gate,
        fixture_reference_renames=renames,
        products=[named(p) for p in sorted(O.glob('*.o'))],
        original_math_and_raw3_recompiled=False))
elif phase == 'finish-build':
    sealed = doc(R / 'attestations/qualification/ed9bf462/native-cpu-v4-raw-fixture-copy-check-failure/MANIFEST.json')
    assert sealed['compiler_and_raw_arithmetic_not_failed'] and sealed['no_numerical_workload_run']
    assert len(sealed['products']) == 6
    for name, expected in sealed['products'].items():
        assert pin(R / name) == expected, name
    renames = []
    for stage, original in (
        ('transport', R / 'bin/linux-gcc/native-cpu-accumulate-compatibility/driver.o'),
        ('endpoint', R / 'bin/linux-gcc/native-cpu-endpoint-header-domain/driver.o')):
        renames.append(rename_reference(original, O / (stage + '-driver.o'),
            'sirius_simd_' + stage + '_packet', 'sirius_simd_v4_' + stage + '_packet',
            recovered=stage == 'transport'))
    write(W / 'build-result.json', dict(passed=True,
        original_compile_phase_passed=False, original_failure=sealed['failure'],
        fixture_copy_verification_completed=True, fixture_reference_renames=renames,
        recovered_products=sealed['products'], new_compilations=0,
        products=[named(p) for p in sorted(O.glob('*.o'))],
        original_math_and_raw3_recompiled=False))
elif phase == 'link':
    assert doc(W / 'code-proof.json')['passed']
    old = doc(R / 'out/native-cpu-accumulate-compatibility/compile-completion.json')
    ap = doc(R / 'out/native-cpu-endpoint-header-provider/provider/source-contract.json')
    objects = [R / row['path'] for row in ap['native_objects'][:7]]
    objects += [R / 'bin/linux-gcc/native-cpu-v4-arithmetic' / (name + '.o')
                for name in ('dispatch', 'v3', 'optional')]
    objects += [O / (stage + '.o') for stage in ('transport', 'endpoint', 'dopri_phase')]
    objects += [O / 'cpu_gate.o']
    rows, libraries = [], {}
    for stage in ('transport', 'endpoint'):
        selected = [O / (stage + '-driver.o'), *objects]
        if stage == 'transport':
            selected += [R / old['renamed_baseline_object']['path'], O / 'rows_adapter.o']
        executable = B / (stage + '-check')
        assert not executable.exists()
        argv = [proof['adapter_compile_argv'][0], *map(str, selected),
                '-Wl,--no-relax', '-lm', '-o', str(executable)]
        subprocess.run(argv, check=True)
        association = associate(executable, [dict(path=str(p), **pin(p)) for p in selected])
        write(W / (stage + '-link-proof.json'), association)
        resolution = subprocess.check_output(['ldd', '-r', str(executable)],
            text=True, stderr=subprocess.STDOUT)
        assert 'not found' not in resolution and 'undefined symbol' not in resolution
        for name in re.findall(r'(/[^\s()]+)', resolution):
            path = Path(name)
            if path.is_file():
                path = path.resolve(strict=True)
                libraries[str(path)] = pin(path)
        rows.append(dict(stage=stage, argv=argv, executable=named(executable)))
    write(W / 'link-result.json', dict(passed=True, executables=rows, runtime_libraries=libraries))
elif phase == 'run':
    results = []
    for stage in ('transport', 'endpoint'):
        result = subprocess.run([str(B / (stage + '-check'))], check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        write(W / (stage + '-output.json'), dict(stdout=result.stdout, stderr=result.stderr))
        print(result.stdout, end='')
        results.append(dict(stage=stage, exit=result.returncode))
    write(W / 'compatibility.json', dict(passed=True, fixtures=results,
        limitations='Original finite Transport38 and focused Endpoint fixture only. DopriPhase object linked but not executed; no performance, hit-rate, provider, C8/all8/ABBA/frame, scientific or release admission.'))
else:
    raise ValueError(phase)
