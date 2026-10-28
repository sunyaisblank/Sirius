"""Seal context/native products, preserved failures and bounded reader audits."""
from pathlib import Path
import hashlib
import json
import subprocess
import zipfile

R = Path.cwd()
W = R / 'out/native-cpu-canonical-closed-remaining'
V = R / 'out/native-cpu-strict-ir-linking'
A = R / 'attestations/qualification/dc584b05/native-cpu-closed-remaining-context-native-gates'
FAILURES = {'cleanup-endpoint', 'cleanup-audit-endpoint'}
ACTIONS = [
    'context-plugin-transport', 'context-transport', 'cleanup-transport', 'native-transport',
    'context-plugin-endpoint', 'context-endpoint', 'cleanup-endpoint',
    'cleanup-audit-plugin-endpoint', 'cleanup-audit-endpoint', 'cleanup-audit-error-capture-endpoint',
    'endpoint-residual-serialization-diagnostic', 'cleanup-canonical-audit-plugin-endpoint',
    'cleanup-canonical-audit-endpoint', 'native-audited-endpoint']


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def doc(path):
    return json.loads(path.read_text())


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not A.exists()
sources = [Path(__file__).resolve()]
sources += [path for path in W.iterdir() if path.is_file() and path.suffix in ('.py', '.json')]
sources += [V / name for name in ['derive_closed_remaining_context.py', 'derive_closed_dp_context.py',
    'audit_closed_endpoint_residual.cpp', 'audit_closed_endpoint_residual_canonical.cpp']]
owners, products, IR = {}, [], {}
current_actions = doc(W / 'actions.json')
for action in ACTIONS:
    folder = W / action
    row = doc(folder / 'owner.json')
    assert row['status'] == 'terminal' and row['passed'] is (action not in FAILURES)
    assert row['child_exit'] == (1 if action in FAILURES or action == 'cleanup-audit-error-capture-endpoint' else 0)
    assert row['source_unchanged'] and row['whole_inputs_unchanged']
    assert not row['remaining_owned_processes'] and not row['cleanup_errors']
    assert row['observed_births'] and all(item['absent'] for item in row['observed_births'])
    assert (folder / 'inputs-before.json').read_bytes() == (folder / 'inputs-after.json').read_bytes()
    sources += [path for path in folder.iterdir() if path.is_file()]
    owners[action] = dict(path=str((folder / 'owner.json').relative_to(R)), identity=pin(folder / 'owner.json'),
        passed=row['passed'], child_exit=row['child_exit'], expected_exit=current_actions[action]['expected_exit'],
        source=row['source'], wall_seconds=row['wall_seconds'], peak_sampled_rss_kib=row['peak_sampled_rss_kib'])
assert len(owners) == 14 and sum(row['passed'] for row in owners.values()) == 12
assert not (W / 'native-endpoint').exists()
for stage in ('transport', 'endpoint'):
    B = R / ('bin/linux-gcc/native-cpu-canonical-closed-' + stage)
    sources += [V / ('closed_' + stage + '_context.cpp'), V / ('closed_' + stage + '_context_contract.json')]
    IR[stage] = [B / (stage + '.' + phase + '.ll') for phase in ('context', 'closed')]
    products += [B / name for name in ('context.so', stage + '.context.bc', stage + '.closed.bc', stage + '.closed.o')]
    cleanup, native = doc(W / (stage + '-cleanup-result.json')), doc(W / (stage + '-native-result.json'))
    assert cleanup['owned_context_allocas'] == 0 and native['emitted_fused_instructions'] == []
    assert native['unwind_table']['full_function_FDE']
    for receipt in (cleanup, native):
        assert pin(R / receipt['output']['path']) == {key:receipt['output'][key] for key in ('bytes', 'sha256')}
products += [R / 'bin/linux-gcc/native-cpu-canonical-closed-endpoint' / name for name in
    ('endpoint.residual-audit.so', 'endpoint.residual-canonical-audit.so')]
proof = doc(W / 'endpoint-canonical-residual-proof.json')
assert proof['passed'] and proof['read_only'] and proof['full_module_fingerprint_unchanged']
assert proof['text_reference_predecessor_comments_normalized'] and proof['residual_bytes'] == 320
assert proof['value_call_total'] == 3954 and proof['staging_bytes'] == 122528
assert doc(W / 'endpoint-cleanup-result.json')['accepted_artifact_owner'] == 'out/native-cpu-canonical-closed-remaining/cleanup-canonical-audit-endpoint/owner.json'
tail = doc(W / 'transport-native-call-proof.json')
assert tail['passed'] and tail['all_16_original_call_paths_reach_the_same_exact_target']
assert tail['IR_call_count'] == 2711 and tail['native_relocation_count'] == 2696
retirement = doc(W / 'endpoint-residual-serialization-diagnostic/serialized-text-retirement-result.json')
assert retirement['payload_bytes_removed'] == 238947628 and retirement['observed_available_bytes_change'] == 238960640
assert all(not (W / 'endpoint-residual-serialization-diagnostic' / name).exists() for name in ('from-bc.ll', 'from-text.ll'))
predecessors = {}
for name in ['attestations/qualification/0de5e321/native-cpu-closed-remaining-strict-value-gates/MANIFEST.json',
    'attestations/qualification/620a0c6e/native-cpu-closed-dp-unwind-original-c8/MANIFEST.json',
    'attestations/qualification/620a0c6e/native-cpu-closed-dp-unwind-original-c8/ALLOCATION_CORRECTION.json']:
    predecessors[name] = pin(R / name)
A.mkdir(parents=True)
bundles = {}
for label, paths in [('source-and-evidence', sources), ('transport-IR', IR['transport']),
                     ('endpoint-IR', IR['endpoint']), ('products', products)]:
    entries = {str(path.relative_to(R)):pin(path) for path in sorted(set(paths))}
    target = A / (label + '.zip')
    with zipfile.ZipFile(target, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name in entries:
            archive.write(R / name, name)
    with zipfile.ZipFile(target) as archive:
        assert len(archive.namelist()) == len(entries) and set(archive.namelist()) == set(entries)
        for name, expected in entries.items():
            with archive.open(name) as stream:
                assert hashlib.file_digest(stream, 'sha256').hexdigest() == expected['sha256']
            assert archive.getinfo(name).file_size == expected['bytes'] and pin(R / name) == expected
    assert target.stat().st_size < 100 * 1024 * 1024
    bundles[label] = dict(path=str(target.relative_to(R)), **pin(target), members=entries)
manifest = dict(status='REMAINING_CONTEXT_AND_NATIVE_COMPILER_GATES_PASS_NUMERICAL_GATES_PENDING',
    custody_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    owners=owners, preserved_failed_owners=sorted(FAILURES), original_unrun_native_action_preserved='native-endpoint',
    accepted_native_owners=dict(transport='native-transport', endpoint='native-audited-endpoint'),
    native_core_relocations=dict(transport=2696, endpoint=3954),
    source_operation_join=tail, bounded_endpoint_residual_proof=proof,
    measured_disposable_cleanup=retirement, predecessors=predecessors, bundles=bundles,
    scope='T92192B/E117792B main Context allocations eliminated. Preserve exact transaction staging T150368B/E122528B and all original fields/masks/FP/loops/whole-state rollback; T retains only its isolated4000B coefficient cache. E retains precisely320B of initialized bounded original selector/projection locals, accepted by a separate read-only audit after the overly strict zero-residual expectation failed. Original cleanup, initial comparator failure and diagnostic costs remain explicit; no optimizer cleanup was repeated. Only nonsemantic predecessor-comment order is normalized for the independently parsed text association; raw bitcode/full read-only fingerprint checks remain exact. T16-to1 static division sharing retains one operation in each original25 cache-loop iteration. Both native objects have no fused instructions, sole hidden C entries and full-range FDE/call-return-CFA inventories; complete GP/prologue/final-linked-code joins and genuine positive-job C8 remain required. No timing benefit, full3-stage numerical/provider adoption, ABBA/Decimal/all8/frame/scientific/hardware/release/main acceptance. The old closed-DP C8 evidence applies only to its own identity; all eight ultimate-goal subgoals remain open.')
with (A / 'MANIFEST.json').open('x') as stream:
    json.dump(manifest, stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(sealed=str(A.relative_to(R)), owners=14, passed=12, preserved_failures=2,
    bundles={name:row['bytes'] for name,row in bundles.items()})))
