"""One baseline adapter compile and twelve-input relocation link; no workload."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

R = Path.cwd()
W = Path(__file__).resolve().parent
B = R / 'bin/linux-gcc/native-cpu-closed-dp-component'
AU = R / 'out/native-cpu-shared-value-component'
DP = R / 'out/native-cpu-canonical-closed-dp'
M = R / 'out/native-cpu-local-effects'
SUFFIXES = [kind + '_' + op + '8' for kind, ops in (
    ('dense', ('add', 'multiply', 'divide', 'sqrt')),
    ('tagged', ('add', 'subtract', 'multiply', 'divide', 'sqrt'))) for op in ops]


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def doc(path):
    return json.loads(path.read_text())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def terminal(path):
    row = doc(path)
    assert row['status'] == 'terminal' and row['passed'] and row['child_exit'] == 0
    assert row['source_unchanged'] and row['whole_inputs_unchanged']
    assert not row['remaining_owned_processes'] and not row['cleanup_errors']
    assert path.with_name('inputs-before.json').read_bytes() == path.with_name('inputs-after.json').read_bytes()


phase, = sys.argv[1:]
if phase == 'prepare':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    derivation = doc(W / 'derivation.json')
    assert pin(W / 'rows_adapter.cpp') == {k: derivation['derived'][k] for k in ('bytes', 'sha256')}
    terminal(AU / 'build/owner.json')
    old = doc(AU / 'build-plan.json')
    reused = old['reused']
    assert len(reused) == 6
    au_values = old['new_value_objects'][1:]
    assert [Path(row['path']).name for row in au_values] == ['shared-transport.o', 'shared-endpoint.cleanup.o', 'shared-math.o']
    for root, name in ((DP, 'native'), (M, 'math-native')):
        terminal(root / name / 'owner.json')
    added = [doc(DP / 'native-result.json')['output'], doc(M / 'math-native-result.json')['output']]
    products = reused + au_values + added
    for row in products:
        assert pin(R / row['path']) == {k: row[k] for k in ('bytes', 'sha256')}
    compile_argv = list(derivation['compile_argv'])
    merge_argv = ['/usr/bin/ld', '-r', '-o', str(B / 'rows_adapter.o'), str(B / 'adapter-only.o')]
    merge_argv += [str(R / row['path']) for row in products]
    sources = [*W.glob('*.py'), W / 'rows_adapter.cpp', W / 'derivation.json',
        AU / 'build-plan.json', AU / 'build-result.json', AU / 'build/owner.json',
        DP / 'native-result.json', DP / 'native/owner.json', M / 'math-native-result.json', M / 'math-native/owner.json',
        R / 'out/native-cpu-dispatch-admission-token/admitted_bridge.h',
        R / 'out/native-cpu-dispatch-admission-token/native_packet.hpp',
        R / 'out/native-cpu-row-ownership/cold-frame/run_owned.py',
        R / 'out/native-cpu-row-ownership/cold-frame/owned_guard.py']
    sources += [R / row['path'] for row in products]
    sources += [Path(compile_argv[0]).resolve(), Path('/usr/bin/ld').resolve()]
    inputs = {str(p.relative_to(R)) if p.is_relative_to(R) else str(p): pin(p) for p in sources}
    inputs.update(derivation['source_inputs'])
    # Protect the complete source/compiler proof of each genuinely new object.
    for plan in (DP / 'context-plan.json', M / 'plan.json'):
        inputs[str(plan.relative_to(R))] = pin(plan)
        inputs.update(doc(plan)['inputs'])
    write(W / 'build-plan.json', dict(inputs=inputs, compile_argv=compile_argv, merge_argv=merge_argv,
        reused=reused, reused_au_value_objects=au_values, new_closed_objects=added,
        scope='One baseline adapter compile and12-part relocation link. Six public AS/GCC/dispatcher constituents and AU scoped T/E/shared9 math remain exact; only admitted DP v4 selects the closed DP and separate existing local-effects regcall math9. Compiled interface only, no numerical transfer.'))
    assert not B.exists()
    write(W / 'actions.json', {'build': dict(
        revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='diagnostic', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=['/usr/bin/python3', '-B', str(W / 'build.py'), 'build'],
        frozen_inputs=sorted(inputs) + ['out/native-cpu-closed-dp-component/build-plan.json'],
        scope='One baseline adapter compile plus exact12-part merged object; public branches and AU T/E unchanged. Static interface only. Fresh original C8/all32 is required; no ABBA/frame gate is prepared or authorized here.')})
elif phase == 'build':
    plan = doc(W / 'build-plan.json')
    for name, expected in plan['inputs'].items():
        assert pin(R / name) == expected, name
    B.mkdir()
    subprocess.run(plan['compile_argv'], check=True)
    subprocess.run(plan['merge_argv'], check=True)
    target = B / 'rows_adapter.o'
    undefined = subprocess.check_output(['nm', '-u', str(target)], text=True)
    assert not any(prefix in undefined for prefix in ('sirius_shared_value_', 'sirius_local_value_', 'sirius_scope_', 'sirius_closed_')), undefined
    listing = subprocess.check_output(['readelf', '-Ws', str(target)], text=True)
    expected = {prefix + suffix for prefix in ('sirius_shared_value_', 'sirius_local_value_') for suffix in SUFFIXES}
    cores, closed = [], []
    for line in listing.splitlines():
        fields = line.split()
        if len(fields) != 8:
            continue
        if fields[-1] in expected or fields[-1].startswith('sirius_closed_dp_packet'):
            assert fields[3:6] == ['FUNC', 'GLOBAL', 'HIDDEN'] and fields[6] != 'UND' and int(fields[2]) > 0
            (cores if fields[-1] in expected else closed).append(fields[-1])
    assert set(cores) == expected and len(cores) == 18
    assert closed == ['sirius_closed_dp_packet']
    write(W / 'build-result.json', dict(passed=True, output=dict(path=str(target.relative_to(R)), **pin(target)),
        adapter=pin(B / 'adapter-only.o'), core_definitions=sorted(cores), closed_raw_exports=closed,
        all_raw_scope_and_core_references_bound=True, remaining_external_symbols=undefined,
        reused=plan['reused'], reused_au_value_objects=plan['reused_au_value_objects'], new_closed_objects=plan['new_closed_objects'],
        compiled_interface_passed=True, numerical_qualification_transferred=False,
        scope='Exact12-part merged object/static symbol binding only. Original genuine-job C8/all32 is the next numerical gate. No ABBA/frame/scientific/release acceptance.'))
    print(json.dumps(dict(passed=True, core_definitions=18, merged_bytes=target.stat().st_size)))
else:
    raise ValueError(phase)
