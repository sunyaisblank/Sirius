"""Owned preparation/execution boundary; this file does nothing without a phase.

Four production baseline compilations plus one diagnostic wrapper. Original
fixture objects/assertions, raw3, math and pure admission are never rebuilt.
"""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import sys
import zipfile

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-dispatch-admission-token'
sys.dont_write_bytecode = True
sys.path.insert(0, str(R / 'out/native-cpu-field-compatibility'))
from elf_proof import Elf


def doc(path):
    return json.loads(path.read_text())


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def named(path):
    return dict(path=str(path.relative_to(R)), **pin(path))


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def terminal(phase):
    owner = doc(W / phase / 'owner.json')
    assert owner['status'] == 'terminal' and owner['passed'] and owner['child_exit'] == 0
    assert owner['source_unchanged'] and owner['whole_inputs_unchanged']
    assert not owner['remaining_owned_processes'] and not owner['cleanup_errors']
    return owner


def verify_inputs():
    proof, finite = doc(W / 'derivation.json'), doc(W / 'finite-inputs.json')
    for name, expected in proof['source_inputs'].items():
        assert pin(R / name) == expected, name
    for name, expected in proof['derived'].items():
        assert pin(W / name) == expected, name
    for name, expected in finite['source_inputs'].items():
        assert pin(R / name) == expected, name
    for row in proof['reused_objects']:
        assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}, row['path']
    assert pin(R / finite['fixture_archive']['path']) == {
        k: finite['fixture_archive'][k] for k in ('bytes', 'sha256')}
    row = finite['baseline_object']
    assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}
    return proof, finite


def rename_proof(original, output, mapping):
    assert not output.exists()
    argv = ['/usr/bin/objcopy']
    for before, after in mapping.items():
        argv += ['--redefine-sym', before + '=' + after]
    subprocess.run(argv + [str(original), str(output)], check=True)
    a, b = Elf(original), Elf(output)
    inverse = {v: k for k, v in mapping.items()}
    assert len(inverse) == len(mapping)
    assert all(any(s['name'] == old for s in a.symbols) for old in mapping)
    assert a.names == b.names and a.kind == b.kind == 1
    assert sorted(a.descriptor(s, {}) for s in a.symbols) == sorted(b.descriptor(s, inverse) for s in b.symbols)
    for i, section in enumerate(a.sections):
        if section[2] & 2:
            assert tuple(section[k] for k in (1, 2, 3, 5, 8, 9)) == tuple(b.sections[i][k] for k in (1, 2, 3, 5, 8, 9))
            assert a.bytes(i) == b.bytes(i)
    def relocs(elf, names):
        return sorted((elf.section_name(sec), offset, kind, addend, elf.descriptor(target, names))
                      for sec, offset, kind, addend, target in elf.relocations)
    assert relocs(a, {}) == relocs(b, inverse)
    return dict(original=named(original), renamed=named(output), mapping=mapping,
        allocated_sections_byte_exact=True, symbols_and_relocations_inverse_exact=True)


def objects_for(stage, proof, finite):
    # New scalar3 + new adapter. Reuse pure admission, fallback raw3, AS raw3,
    # baseline dispatcher/original-v3, and the qualified optional math.
    paths = [B / (name + '.o') for name in ('transport', 'endpoint', 'dopri_phase', 'rows_adapter')]
    paths += [R / row['path'] for row in proof['reused_objects']]
    paths += [B / 'fixture_routes.o', B / 'cpu_gate.o']
    if stage == 'transport':
        paths.append(R / finite['baseline_object']['path'])
    assert len(paths) == len(set(paths))
    return [B / (stage + '-driver.o'), *paths]


def exact_route_joins(executable, selected, word):
    linked = Elf(executable)
    wanted = {
        'sirius_native_transport_admit_packet', 'sirius_native_endpoint_admit_packet',
        'sirius_native_transport_rows_try_packet', 'sirius_native_endpoint_rows_try_packet',
        'sirius_native_transport_rows_admitted_packet', 'sirius_native_endpoint_rows_admitted_packet',
        'sirius_native_endpoint_try_packet', 'sirius_native_endpoint_admitted_packet',
        'sirius_native_transport_admitted_packet', 'sirius_native_dopri_phase_admitted_packet',
    }
    joins, reached = [], set()
    for path in (B / 'fixture_routes.o', B / 'rows_adapter.o'):
        assert path in selected
        obj = Elf(path)
        for section, offset, kind, addend, target in obj.relocations:
            if target['name'] not in wanted:
                continue
            original, = [s for s in obj.symbols if s['section'] == section and
                          (s['info'] & 15) == 2 and s['value'] <= offset < s['value'] + s['size']]
            final, = [s for s in linked.symbols if s['section'] and s['name'] == original['name']]
            destination, = [s for s in linked.symbols if s['section'] and s['name'] == target['name']]
            field = final['value'] + offset - original['value']
            assert kind in (2, 4, 42)
            address = field + word(linked, field, '<i') - addend
            if kind == 42:
                slot = address
                address = word(linked, slot, '<Q')
                dynamic, = [(k, a, t) for _, location, k, a, t in linked.relocations if location == slot]
                assert dynamic[0] == 8 and dynamic[1] == destination['value'] and not dynamic[2]['name']
            assert address == destination['value'], (path, original['name'], target['name'])
            reached.add(target['name'])
            joins.append(dict(object=str(path.relative_to(R)), caller=original['name'],
                relocation_type=kind, target=target['name'], field_address=field,
                actual_target_address=address, expected_target_address=destination['value']))
    assert reached == wanted, (wanted - reached, reached - wanted)
    return joins


phase, = sys.argv[1:]
assert phase in ('prepare-build', 'prepare-link', 'prepare-run', 'build', 'link', 'run')
proof, finite = verify_inputs()
if phase.startswith('prepare-'):
    step = phase.removeprefix('prepare-')
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    inputs = set(proof['source_inputs']) | set(finite['source_inputs'])
    inputs |= {str((W / name).relative_to(R)) for name in proof['derived']}
    inputs |= {row['path'] for row in proof['reused_objects']}
    inputs |= {str(path.relative_to(R)) for path in W.iterdir() if path.is_file()}
    # Same include closure as the original scalar recipes; no toolchain copy.
    for directory in ('bin/linux-gcc/native-cpu-scalar-zero/typed-batch',
                      'bin/linux-gcc/native-cpu-static-program/typed-batch',
                      'bin/linux-gcc/native-cpu-slang-bridge/typed-batch',
                      'out/native-cpu-slang-bridge', '/opt/slang/include'):
        for path in (R / directory).rglob('*'):
            if path.is_file():
                inputs.add(str(path.relative_to(R)) if path.is_relative_to(R) else str(path))
    for tool in ('/usr/bin/x86_64-linux-gnu-g++-14', '/usr/bin/objcopy', '/usr/bin/objdump',
                 '/usr/bin/ld', '/usr/bin/ldd', '/usr/bin/readelf'):
        inputs.add(str(Path(tool).resolve(strict=True)))
    actions = {} if step == 'build' else doc(W / 'actions.json')
    if step == 'build':
        assert not B.exists() and not (W / 'actions.json').exists()
    else:
        terminal('build')
        assert doc(W / 'code-proof.json')['passed']
        inputs |= {str(path.relative_to(R)) for path in B.glob('*.o')}
    if step == 'run':
        terminal('link')
        inputs |= {str((B / (stage + '-check')).relative_to(R)) for stage in ('transport', 'endpoint')}
        inputs |= set(doc(W / 'link-result.json')['runtime_libraries'])
    assert step not in actions and not (W / step).exists()
    assert all((R / name).is_file() for name in inputs)
    actions[step] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='diagnostic', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=['/usr/bin/python3', '-B', str(W / 'qualify.py'), step], frozen_inputs=sorted(inputs),
        scope='One bounded dispatch-admission-token component phase. Four production baseline compiles plus one diagnostic-wrapper compile; sealed original fixture objects/raw3/math/pure-admission reused. Original Transport38 and Endpoint fixtures each run once in public mode and once in actual-pure-preflight/admitted mode, with all original assertions/counts unchanged per invocation. No additional cases/retries, benchmark, Mesa capability-branch observation, provider/frame/scientific/release admission.')
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared=step, executed=False, frozen_inputs=len(inputs))))
elif phase == 'build':
    B.mkdir(exist_ok=False)
    production = []
    for command in proof['prospective_compile_commands']:
        subprocess.run(command['argv'], check=True)
        production.append(dict(role=command['role'], argv=command['argv'],
                               object=named(Path(command['argv'][command['argv'].index('-o') + 1]))))
    diagnostic_argv = list(proof['prospective_compile_commands'][-1]['argv'])
    diagnostic_argv[diagnostic_argv.index('-c') + 1] = str(W / 'fixture_routes.cpp')
    diagnostic_argv[diagnostic_argv.index('-o') + 1] = str(B / 'fixture_routes.o')
    subprocess.run(diagnostic_argv, check=True)
    restored, renamed = [], []
    with zipfile.ZipFile(R / finite['fixture_archive']['path']) as archive:
        for name, item in finite['archive_objects'].items():
            data = archive.read(item['member'])
            assert dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest()) == item['pin']
            destination = B / name
            with destination.open('xb') as stream:
                stream.write(data)
            restored.append(dict(**named(destination), archive_member=item['member']))
    for stage, mapping in finite['fixture_symbol_rewrites'].items():
        renamed.append(rename_proof(B / (stage + '-original-driver.o'), B / (stage + '-driver.o'), mapping))
    verify_inputs()
    subprocess.run(['/usr/bin/python3', '-B', str(W / 'check_code.py')], check=True)
    assert doc(W / 'code-proof.json')['passed']
    write(W / 'build-result.json', dict(passed=True, production_objects=production,
        new_production_compilations=4, diagnostic_object=dict(argv=diagnostic_argv, object=named(B / 'fixture_routes.o')),
        new_diagnostic_compilations=1, original_fixture_test_raw_math_and_admission_rebuilds=0,
        restored_fixture_objects=restored, symbol_rename_proofs=renamed))
elif phase == 'link':
    terminal('build')
    assert doc(W / 'code-proof.json')['passed']
    # Exact original ELF association module; only the declared fixture GOT
    # target name changes from direct raw to this diagnostic route symbol.
    module = R / 'out/native-cpu-v4-raw-packets/fixture_elf.py'
    text = module.read_text()
    assert text.count("'sirius_simd_v4_transport_packet'") == 1
    text = text.replace("'sirius_simd_v4_transport_packet'", "'sirius_token_transport_raw'", 1)
    scope = {}
    exec(compile(text, str(module), 'exec'), scope)
    results, libraries = [], {}
    for stage in ('transport', 'endpoint'):
        selected = objects_for(stage, proof, finite)
        executable = B / (stage + '-check')
        assert not executable.exists()
        argv = [proof['prospective_compile_commands'][0]['argv'][0], *map(str, selected),
                '-Wl,--no-relax', '-lm', '-o', str(executable)]
        subprocess.run(argv, check=True)
        joined = scope['associate'](executable, [dict(path=str(path), **pin(path)) for path in selected])
        joined['exact_preflight_public_and_private_target_joins'] = exact_route_joins(executable, selected, scope['word'])
        joined['scope'] += ' The candidate fixture GOT target is the actual diagnostic routing function; original fixture instruction/data bytes and inverse-normalized relocations remain exact. Every selected object function is associated; this is not Mesa capability-branch execution.'
        write(W / (stage + '-link-proof.json'), joined)
        resolution = subprocess.check_output(['ldd', '-r', str(executable)], text=True, stderr=subprocess.STDOUT)
        assert 'not found' not in resolution and 'undefined symbol' not in resolution
        notes = subprocess.check_output(['readelf', '-n', str(executable)], text=True)
        assert not re.search(r'ISA needed:.*x86-64-v[234]', notes)
        for name in re.findall(r'(/[^\s()]+)', resolution):
            path = Path(name)
            if path.is_file():
                path = path.resolve(strict=True)
                libraries[str(path)] = pin(path)
        results.append(dict(stage=stage, argv=argv, executable=named(executable)))
    write(W / 'link-result.json', dict(passed=True, executables=results, runtime_libraries=libraries,
        fixture_GOT_targets_inverse_checked=True, all_linked_functions_associated=True,
        runtime_minimum_not_v3_or_v4=True))
elif phase == 'run':
    terminal('build')
    terminal('link')
    assert doc(W / 'code-proof.json')['passed']
    results = {}
    for mode in ('public', 'admitted'):
        results[mode] = {}
        for stage in ('transport', 'endpoint'):
            folder = W / 'run' / mode / stage
            folder.mkdir(parents=True, exist_ok=False)
            route_file = folder / 'routes.json'
            env = dict(os.environ, SIRIUS_TOKEN_FIXTURE_MODE=mode, SIRIUS_TOKEN_FIXTURE_ROUTES=str(route_file))
            argv = [str(B / (stage + '-check'))] + ([str(folder)] if stage == 'transport' else [])
            ran = subprocess.run(argv, check=True, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            assert not ran.stderr
            comparison, routes = json.loads(ran.stdout), doc(route_file)
            assert comparison['passed'] and routes['mode'] == mode
            for name, expected in finite['fixture_expected'][stage].items():
                assert comparison[name] == expected, (stage, mode, name, comparison[name], expected)
            for name, expected in finite['route_expected'][stage][mode].items():
                assert routes[name] == expected, (stage, mode, name, routes[name], expected)
            assert routes['trusted_call_after_failed_preflight'] is False
            fixture_environment = {key: env[key] for key in
                                   ('SIRIUS_TOKEN_FIXTURE_MODE', 'SIRIUS_TOKEN_FIXTURE_ROUTES')}
            results[mode][stage] = dict(result=comparison, routes=routes, argv=argv,
                                       child_only_fixture_environment=fixture_environment)
            write(folder / 'output.json', dict(stdout=ran.stdout, stderr=ran.stderr,
                argv=argv, child_only_fixture_environment=fixture_environment))
            print(json.dumps(dict(stage=stage, mode=mode, passed=True)))
    write(W / 'compatibility.json', dict(passed=True, fixtures=results,
        fixture_invocations=4, original_assertions_counts_and_sources_unchanged_per_invocation=True,
        actual_unchanged_pure_preflight_before_private_calls=True,
        trusted_calls_after_failed_preflight=0, public_and_private_routes_validated=True,
        limitations='Finite original Transport38 and focused Endpoint fixtures only. Public guards and diagnostic pure-preflight/private routing, full banks/inactive words/canaries, all21142 Endpoint malformed words/21 bounds, header-family refusals and original FP assertions retained. DopriPhase compiles/links but does not execute here. No actual Mesa job-token branch, speed, provider/C8/all8/ABBA/frame/scientific/release acceptance; fresh C8 must observe that branch and hidden entries.'))
