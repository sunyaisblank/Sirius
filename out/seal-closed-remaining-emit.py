"""Preserve first ordinary closed Transport/Endpoint IR and exact source proofs."""
from pathlib import Path
import hashlib
import json
import subprocess
import zipfile

R = Path.cwd()
W = R / 'out/native-cpu-canonical-closed-remaining'
V = R / 'out/native-cpu-strict-ir-linking'
A = R / 'attestations/qualification/d970f13e/native-cpu-closed-remaining-first-emit'


def pin(path):
    with path.open('rb') as stream:
        return dict(bytes=path.stat().st_size, sha256=hashlib.file_digest(stream, 'sha256').hexdigest())


def doc(path):
    return json.loads(path.read_text())


assert not subprocess.check_output(['git', 'status', '--porcelain']).strip()
assert not A.exists()
sources = [Path(__file__).resolve()]
sources += [path for path in W.rglob('*') if path.is_file() and '__pycache__' not in path.parts]
sources += [V / name for name in ('canonical_private_remaining.py', 'canonical_private_remaining_contract.json')]
owners, inventories = {}, {}
for action in ('derive', 'emit-transport', 'emit-endpoint'):
    path = W / action / 'owner.json'
    row = doc(path)
    assert row['status'] == 'terminal' and row['passed'] and row['child_exit'] == 0
    assert row['source_unchanged'] and row['whole_inputs_unchanged']
    assert not row['remaining_owned_processes'] and not row['cleanup_errors']
    assert row['observed_births'] and all(item['absent'] for item in row['observed_births'])
    assert path.with_name('inputs-before.json').read_bytes() == path.with_name('inputs-after.json').read_bytes()
    owners[action] = dict(path=str(path.relative_to(R)), identity=pin(path), source=row['source'],
        wall_seconds=row['wall_seconds'], peak_sampled_rss_kib=row['peak_sampled_rss_kib'])
for stage in ('transport', 'endpoint'):
    inventory = doc(W / (stage + '-inventory.json'))
    assert inventory['fp_instruction_flags'] == {}
    for key in ('output', 'source'):
        path = R / inventory[key]['path']
        assert pin(path) == {name:inventory[key][name] for name in ('bytes', 'sha256')}
        sources.append(path)
    inventories[stage] = inventory
predecessor = 'attestations/qualification/620a0c6e/native-cpu-closed-dp-unwind-original-c8/MANIFEST.json'
assert pin(R / predecessor)['sha256'] == '3521b941e2c5aded85b6949846de6efc64cdd98b63b7ac9044a191c51dd6184b'
A.mkdir(parents=True)
entries = {str(path.relative_to(R)):pin(path) for path in sorted(set(sources))}
target = A / 'source-and-ordinary-ir.zip'
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
manifest = dict(status='FIRST_CLOSED_T_E_ORDINARY_IR_PASS_NO_STRICT_OR_RUNTIME_VERDICT',
    custody_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
    owners=owners, predecessor={predecessor:pin(R / predecessor)},
    bundle=dict(path=str(target.relative_to(R)), **pin(target), members=entries), inventories=inventories,
    scope='One source derivation and first ordinary ISPC emission per remaining stage, all three owners PASS. Exact7292 remaining canonical tuples/four table profiles, mixed family/phase order, root validation before publication and Endpoint11241 live guard remain bound. Actual IR function/operation/math inventories are retained. This is ordinary compiler evidence only; no strict conversion, whole-context normalization, native object/FDE, numerical/provider/runtime/ABBA/Decimal/all8/frame/scientific/hardware/release acceptance. The independently sealed closed DP C8 pass applies only to its own provider. All eight ultimate-goal subgoals remain open.')
with (A / 'MANIFEST.json').open('x') as stream:
    json.dump(manifest, stream, indent=2)
    stream.write('\n')
print(json.dumps(dict(sealed=str(A.relative_to(R)), passed_owners=3, archive_bytes=target.stat().st_size)))
