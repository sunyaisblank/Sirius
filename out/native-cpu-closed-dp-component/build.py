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
READER = 'out/native-cpu-closed-dp-component/build.py'
AUDIT_TOOLS = {name: Path('/usr/bin') / name for name in ('nm', 'readelf')}
FAILED_OWNER = dict(bytes=1840, sha256='33450b323f08bdea73d0ef3348f28798c19e690f674488728713adba26d1463b')
ORIGINAL_PLAN = dict(bytes=21530, sha256='dc19ffb3ab77ea20cea39a8dfabb806e7b69355b5a025f07db823a62714cbaa8')
RETAINED_PRODUCTS = {
    'bin/linux-gcc/native-cpu-closed-dp-component/adapter-only.o': dict(bytes=7680,
        sha256='884c89b1bc837009efa855cf4eb198c0f17caf622e3021ee4097cfeb9ed3bdf6'),
    'bin/linux-gcc/native-cpu-closed-dp-component/rows_adapter.o': dict(bytes=22724920,
        sha256='9ffdb4519cb8065b49b0b25b8ae539bef24d300fe14becfd86b55368c27bdbb7')}


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def doc(path):
    return json.loads(path.read_text())


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def terminal(path, success=True):
    row = doc(path)
    assert row['status'] == 'terminal' and row['passed'] is success and row['child_exit'] == (0 if success else 1)
    assert row['source_unchanged'] and row['whole_inputs_unchanged']
    assert not row['remaining_owned_processes'] and not row['cleanup_errors']
    assert row['observed_births'] and all(birth['absent'] for birth in row['observed_births'])
    assert path.with_name('inputs-before.json').read_bytes() == path.with_name('inputs-after.json').read_bytes()


def retained_failure():
    assert pin(W / 'build/owner.json') == FAILED_OWNER
    assert pin(W / 'build-plan.json') == ORIGINAL_PLAN
    terminal(W / 'build/owner.json', success=False)
    assert pin(W / 'build/stderr.log') == dict(bytes=400,
        sha256='ebbd36822d58b1a055b2c5ddd1189f30222f0265d9d06267aa66f07163f03c08')
    assert "ValueError: invalid literal for int() with base 10: '0x204cb1'" in (W / 'build/stderr.log').read_text()
    for name, expected in RETAINED_PRODUCTS.items():
        assert pin(R / name) == expected, name


def symbol_size(text):
    # GNU readelf uses decimal for short sizes and0x-prefixed hex for large ones.
    assert text and (text.isdecimal() or text.startswith('0x') and all(c in '0123456789abcdefABCDEF' for c in text[2:]))
    return int(text, 16 if text.startswith('0x') else 10)


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
elif phase == 'prepare-audit':
    assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
    retained_failure()
    assert not (W / 'build-result.json').exists() and not (W / 'build-audit').exists()
    plan = doc(W / 'build-plan.json')
    inputs = dict(plan['inputs'])
    override = dict(path=READER, before=inputs[READER], after=pin(R / READER))
    assert override['before'] != override['after']
    inputs[READER] = override['after']
    for name, expected in inputs.items():
        assert pin(R / name) == expected, name
    inputs.update(RETAINED_PRODUCTS)
    inputs.update({str(path.resolve()): pin(path.resolve()) for path in AUDIT_TOOLS.values()})
    history = [W / 'build-plan.json', W / 'build/owner.json', W / 'build/inputs-before.json',
        W / 'build/inputs-after.json', W / 'build/stdout.log', W / 'build/stderr.log', W / 'build/samples.jsonl']
    inputs.update({str(path.relative_to(R)): pin(path) for path in history})
    actions = doc(W / 'actions.json')
    assert set(actions) == {'build'}
    prior = W / 'build-audit-prior-actions.json'
    with prior.open('xb') as stream:
        stream.write((W / 'actions.json').read_bytes())
    original_actions = doc(W / 'build/inputs-before.json')['out/native-cpu-closed-dp-component/actions.json']
    assert pin(prior) == {key: original_actions[key] for key in ('bytes', 'sha256')}
    inputs[str(prior.relative_to(R))] = pin(prior)
    write(W / 'build-audit-plan.json', dict(inputs=inputs, original_build_plan=ORIGINAL_PLAN,
        original_failed_owner=FAILED_OWNER, original_compile_argv=plan['compile_argv'],
        original_merge_argv=plan['merge_argv'], source_hash_override=override,
        retained_products=RETAINED_PRODUCTS, compilation_repeated=False, relocation_link_repeated=False,
        accepted_artifact_owner='out/native-cpu-closed-dp-component/build-audit/owner.json'))
    actions['build-audit'] = dict(revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        kind='diagnostic', expected_exit=0, software_provider=False,
        limits=dict(seconds=300, rss_kib=4194304, reserve_kib=2097152),
        argv=['/usr/bin/python3', '-B', str(W / 'build.py'), 'audit'],
        frozen_inputs=sorted(inputs) + ['out/native-cpu-closed-dp-component/build-audit-plan.json'],
        scope='Read-only nm/readelf audit of exact retained adapter/12-part link after decimal-only symbol-size reader failed. Original failed owner, plan/argv, logs and all products preserved; only current reader hash replaces its historical input. No compiler/link/provider/workload rerun or numerical qualification.')
    (W / 'actions.json').write_text(json.dumps(actions, indent=2) + '\n')
    print(json.dumps(dict(prepared='build-audit', compilation_repeated=False, relocation_link_repeated=False)))
elif phase in ('build', 'audit'):
    plan = doc(W / 'build-plan.json')
    if phase == 'audit':
        retained_failure()
        audit = doc(W / 'build-audit-plan.json')
        assert audit['original_build_plan'] == ORIGINAL_PLAN and audit['original_failed_owner'] == FAILED_OWNER
        assert audit['original_compile_argv'] == plan['compile_argv'] and audit['original_merge_argv'] == plan['merge_argv']
        assert audit['retained_products'] == RETAINED_PRODUCTS
        assert audit['source_hash_override'] == dict(path=READER, before=plan['inputs'][READER], after=pin(R / READER))
        assert not audit['compilation_repeated'] and not audit['relocation_link_repeated']
        assert audit['accepted_artifact_owner'] == 'out/native-cpu-closed-dp-component/build-audit/owner.json'
        for name, expected in audit['inputs'].items():
            assert pin(R / name) == expected, name
    else:
        for name, expected in plan['inputs'].items():
            assert pin(R / name) == expected, name
        B.mkdir()
        subprocess.run(plan['compile_argv'], check=True)
        subprocess.run(plan['merge_argv'], check=True)
    target = B / 'rows_adapter.o'
    undefined = subprocess.check_output([str(AUDIT_TOOLS['nm']), '-u', str(target)], text=True)
    assert not any(prefix in undefined for prefix in ('sirius_shared_value_', 'sirius_local_value_', 'sirius_scope_', 'sirius_closed_')), undefined
    listing = subprocess.check_output([str(AUDIT_TOOLS['readelf']), '-Ws', str(target)], text=True)
    expected = {prefix + suffix for prefix in ('sirius_shared_value_', 'sirius_local_value_') for suffix in SUFFIXES}
    cores, closed = [], []
    for line in listing.splitlines():
        fields = line.split()
        if len(fields) != 8:
            continue
        if fields[-1] in expected or fields[-1].startswith('sirius_closed_dp_packet'):
            assert fields[3:6] == ['FUNC', 'GLOBAL', 'HIDDEN'] and fields[6] != 'UND' and symbol_size(fields[2]) > 0
            (cores if fields[-1] in expected else closed).append(fields[-1])
    assert set(cores) == expected and len(cores) == 18
    assert closed == ['sirius_closed_dp_packet']
    assert phase == 'audit', 'This retained candidate is accepted only through its separate read-only audit owner.'
    write(W / 'build-result.json', dict(passed=True, output=dict(path=str(target.relative_to(R)), **pin(target)),
        adapter=pin(B / 'adapter-only.o'), core_definitions=sorted(cores), closed_raw_exports=closed,
        all_raw_scope_and_core_references_bound=True, remaining_external_symbols=undefined,
        reused=plan['reused'], reused_au_value_objects=plan['reused_au_value_objects'], new_closed_objects=plan['new_closed_objects'],
        accepted_artifact_owner=audit['accepted_artifact_owner'], original_failed_owner=FAILED_OWNER,
        original_build_plan=ORIGINAL_PLAN, audit_plan=pin(W / 'build-audit-plan.json'),
        compilation_repeated=False, relocation_link_repeated=False,
        compiled_interface_passed=True, numerical_qualification_transferred=False,
        scope='Exact12-part merged object/static symbol binding only. Original genuine-job C8/all32 is the next numerical gate. No ABBA/frame/scientific/release acceptance.'))
    print(json.dumps(dict(passed=True, core_definitions=18, merged_bytes=target.stat().st_size)))
else:
    raise ValueError(phase)
